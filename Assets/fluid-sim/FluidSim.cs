using System;
using System.Runtime.CompilerServices;
using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using Unity.Profiling;
using UnityEngine;
using Random = UnityEngine.Random;

public class FluidSim : MonoBehaviour
{
    public int m_ParticleCount;

    public int width;
    public int height;
    public float Gravity = 9.8f;
    public float Mass = 1;
    public float SmoothingLength = 10;
    public float m_PressureMultiplier = 1;
    public bool UseAdaptativeStepTime = false;
    [Range(0, 0.3f)]
    public float SimulationStep = 0.0001f;
    public float BoundaryPushStrength = 1;
    public float CollisionDamping = 0.5f;
    public bool AutoStep = true;
    public float ViscosityFactor = 0.5f;
    public float m_InteractionStrength;
    
    NativeArray<Vector2> m_Position;
    NativeArray<Vector2> m_PredictedPosition;
    NativeArray<float> m_Density;
    NativeArray<float> m_Pressure;
    NativeArray<Vector2> m_Velocity;
    internal GridSpatialLookup m_LookupHelper;

    float m_SquaredSmoothingLength;
    float m_KernelTerm;
    float m_KernelDerivativeTerm;
    public float m_TargetDensity;
    float m_LastStepMaxVelocity;
    Vector2 m_MousePosition;
    float m_MouseRadius;
    InteractionDirection m_InteractionDirection;

    static readonly ProfilerMarker s_UpdatePerfMarker = new ProfilerMarker("FluidSim.Update");
    static readonly ProfilerMarker s_StepPerfMarker = new ProfilerMarker("FluidSim.Step");
    static readonly ProfilerMarker s_PressurePerfMarker = new ProfilerMarker("FluidSim.PressureCalc");
    static readonly ProfilerMarker s_DensityPerfMarker = new ProfilerMarker("FluidSim.DensityeCalc");

    public ReadOnlySpan<Vector2> GetPositions() { return m_Position.AsReadOnlySpan(); }
    public ReadOnlySpan<float> GetDensities() { return m_Density.AsReadOnlySpan(); }
    public ReadOnlySpan<float> GetPressures() { return m_Pressure.AsReadOnlySpan(); }
    public ReadOnlySpan<Vector2> GetVelocities() { return m_Velocity.AsReadOnlySpan(); }
    
    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        InitParticles();
    }

    void OnDisable()
    {
        CleanupParticles();
        CleanupSpatialAcceleration();
    }

    void CleanupParticles()
    {
        if (m_Position.IsCreated) 
            m_Position.Dispose();
        if (m_Density.IsCreated)
            m_Density.Dispose();
        if (m_Pressure.IsCreated)
            m_Pressure.Dispose();
        if (m_Velocity.IsCreated)
            m_Velocity.Dispose();
        if(m_PredictedPosition.IsCreated)
            m_PredictedPosition.Dispose();
    }

    void CleanupSpatialAcceleration()
    {
        if (m_LookupHelper.IsValid)
            m_LookupHelper.Dispose();
    }

    void InitializeSpatialAcceleration()
    {
        CleanupSpatialAcceleration();
        m_LookupHelper = new GridSpatialLookup(SmoothingLength, m_ParticleCount);
    }

    void InitParticles()
    {
        CleanupParticles();
        m_Position = new NativeArray<Vector2>(m_ParticleCount, Allocator.Persistent);
        m_PredictedPosition = new NativeArray<Vector2>(m_ParticleCount, Allocator.Persistent);
        m_Density = new NativeArray<float>(m_ParticleCount, Allocator.Persistent);
        m_Pressure = new NativeArray<float>(m_ParticleCount, Allocator.Persistent);
        m_Velocity = new NativeArray<Vector2>(m_ParticleCount, Allocator.Persistent);
        for ( var i = 0; i < m_ParticleCount; i++ )
        {
            var position = new Vector2(Random.Range(0f, width), Random.Range(0f, height));
            m_Position[i] = position;
            m_PredictedPosition[i] = position;
            m_Density[i] = 0;
            m_Pressure[i] = 0;
            m_Velocity[i] = Vector2.zero;
        }
        
    }
    
    void DoReInitializationIfNecessary()
    {
        var particleCountDifferent = m_ParticleCount != m_Position.Length;
        var smoothingRadiusDifferent = !m_LookupHelper.IsValid || !Mathf.Approximately(SmoothingLength, m_LookupHelper.CellSize);
        if (particleCountDifferent)
        {
            InitParticles();
        }
        if (particleCountDifferent || smoothingRadiusDifferent)
        {
            InitializeSpatialAcceleration();
        }
    }

    // Update is called once per frame
    void Update()
    {
        DoReInitializationIfNecessary();
        
        if (AutoStep)
            DoUpdate();
    }

    public void DoUpdate()
    {
        using var markerScope = s_UpdatePerfMarker.Auto();
        
        CachePrecomputedValues();
        var stepTime = SimulationStep;
        if (UseAdaptativeStepTime)
        {
            stepTime = Mathf.Min(SimulationStep, CalcCFLTimeStep(m_LastStepMaxVelocity));
        }
        
        JobHandle dependencyJobHandle = default;

        var predictPositionJob = new PredictPositionJob
        {
            positions = m_Position,
            velocities = m_Velocity,
            predictedPositions = m_PredictedPosition,
            deltaTime = stepTime,
        };
        var predictPositionJobHandle = predictPositionJob.ScheduleParallelByRef(m_PredictedPosition.Length,
            64, dependencyJobHandle);
        predictPositionJobHandle.Complete();
        
        m_LookupHelper.UpdateParticles(m_PredictedPosition);
        
        CalculateParticlesDensity();
        CalculateParticlePressure();
        SimulateStep(stepTime);
        m_MouseRadius = 0; //clear mouse interaction
    }

    public enum InteractionDirection
    {
        Attract,
        Repel
    }
    public void Interact(Vector2 position, float radius, InteractionDirection direction)
    {
        m_MousePosition = position;
        m_MouseRadius = radius;
        m_InteractionDirection = direction;
    }
    
    [BurstCompile]
    struct PredictPositionJob : IJobFor
    {
        [ReadOnly]
        public NativeArray<Vector2> velocities;
        [ReadOnly]
        public NativeArray<Vector2> positions;
        public NativeArray<Vector2> predictedPositions;
        public float deltaTime;
        public void Execute(int index)
        {
            predictedPositions[index] = positions[index] + velocities[index] * deltaTime;
        }
    }
    
    [BurstCompile]
    struct CalculateAccelerationFromExternalForcesJob : IJobFor
    {
        [ReadOnly]
        public NativeArray<Vector2> positions;
        public NativeArray<Vector2> velocities;
        public float deltaTime;
        public float gravity;
        
        public Vector2 forceCenter;
        public float forceRadius;
        public float forceStrength;
        
        public void Execute(int index)
        {
            var position = positions[index];
            var velocity = velocities[index];
            var interactionForce = Vector2.zero;
            
            if (forceRadius > 0)
                interactionForce = InteractionForce(position, velocity, forceCenter, forceRadius, forceStrength);
            
            velocities[index] += (Vector2.down * gravity + interactionForce) * deltaTime;
        }
        
        static Vector2 InteractionForce(Vector2 particlePosition, Vector2 particleVelocity, Vector2 forceCenter, float radius, float strength)
        {
            var diff = forceCenter - particlePosition;
            var sqrDistance = diff.sqrMagnitude;
            
            if (sqrDistance > radius * radius) return Vector2.zero;
            
            var distance = Mathf.Sqrt(sqrDistance);
            var dirToInputPoint = distance <= float.Epsilon ? Vector2.zero : diff/distance;
            var centerT = 1 - distance / radius;
            return  (dirToInputPoint * strength - particleVelocity) * centerT;
        }
    }
    
    [BurstCompile]
    struct CalculatePositionFromVelocityJobFor : IJobFor
    {
        public NativeArray<Vector2> positions;
        public NativeArray<Vector2> velocities;

        public float deltaTime;
        public float width;
        public float height;
        public float smoothingLength;
        public float boundaryPushStrength;
        public float collisionDamping;
        public void Execute(int index)
        {
            var position = positions[index];
            var velocity = velocities[index];
            position += velocity * deltaTime;
            
            if (position.x < smoothingLength)
            {
                var strength = (smoothingLength - position.x) / smoothingLength;
                velocity.x += boundaryPushStrength * strength * deltaTime;
                if (position.x < 0) {
                    position.x = -position.x;
                    velocity.x *= -collisionDamping;
                }
            }
            if (position.y < smoothingLength)
            {
                var strength = (smoothingLength - position.y) / smoothingLength;
                velocity.y += boundaryPushStrength * strength * deltaTime;
                if (position.y < 0)
                {
                    position.y = - position.y;
                    velocity.y *= -collisionDamping;
                }
            }

            var distToRightWall = Mathf.Abs(width - position.x);
            if (distToRightWall < smoothingLength)
            {
                var strength = (smoothingLength - distToRightWall) / smoothingLength;
                velocity.x -= boundaryPushStrength * strength * deltaTime;
                if (position.x > width)
                {
                    position.x = width - (position.x-width);
                    velocity.x *= -collisionDamping;
                }
            }
            var distToBottomWall = Mathf.Abs(height - position.y);
            if (distToBottomWall < smoothingLength)
            {
                var strength = (smoothingLength - distToBottomWall) / smoothingLength;
                velocity.y -= boundaryPushStrength * strength * deltaTime;
                if (position.y > height)
                {
                    position.y = height;
                    velocity.y *= -collisionDamping;
                }
            }
            positions[index] = position;
            velocities[index] = velocity;
        }
    }
    
    [BurstCompile]
    struct SetVelocityFromPressureJobFor : IJobFor
    {
        public NativeArray<Vector2> velocities;
        [ReadOnly]
        public NativeArray<Vector2> positions;
        [ReadOnly]
        public NativeArray<float> densities;
        [ReadOnly]
        public NativeArray<float> pressures;
        [ReadOnly]
        public GridSpatialLookup lookupHelper;
        public float mass;
        public float squaredSmoothingLength;
        public float smoothingLength;
        public float kernelDerivativeTerm;
        public float deltaTime;
        public void Execute(int index)
        {
            var density = densities[index]; 
            var pressureForce = CalculatePressureGradient(index, lookupHelper, mass, squaredSmoothingLength, smoothingLength, kernelDerivativeTerm,
                positions, pressures, densities);
            var pressureAcceleration = pressureForce / density;
            var velocity = pressureAcceleration * deltaTime;
            velocities[index] += velocity;
        }
    }
    
    [BurstCompile]
    struct CalculateViscosityForceJob : IJobFor
    {
        public NativeArray<Vector2> velocities;
        [ReadOnly]
        public NativeArray<Vector2> positions;
        [ReadOnly]
        public NativeArray<float> densities;
        [ReadOnly]
        public GridSpatialLookup lookupHelper;
        public float mass;
        public float smoothingLength;
        public float squaredSmoothingLength;

        public float viscosityFactor;
        public float deltaTime;
        public void Execute(int index)
        {
            var density = densities[index]; 
            var viscosityForce = CalculateViscosityForce(index, lookupHelper, mass, smoothingLength, squaredSmoothingLength, viscosityFactor,
                positions, velocities, densities);
            var viscosityAcceleration = viscosityForce / density;
            var velocity = viscosityAcceleration * deltaTime;
            velocities[index] += velocity;
        }
    }

    void SimulateStep(float deltaTime)
    {
        using var markerScope = s_StepPerfMarker.Auto();
        
        var strengthSign = m_InteractionDirection == InteractionDirection.Attract ? 1 : -1;
        JobHandle dependencyJobHandle = default;
        var gravityJob = new CalculateAccelerationFromExternalForcesJob()
        {
            positions = m_PredictedPosition,
            velocities = m_Velocity,
            gravity = Gravity,
            deltaTime = deltaTime,
            forceCenter = m_MousePosition,
            forceRadius = m_MouseRadius,
            forceStrength = m_InteractionStrength * strengthSign,
        };
        var velocityFromGravityJobHandle = gravityJob.ScheduleParallelByRef(m_Velocity.Length,
            64, dependencyJobHandle);
        
        var velocityFromPressureJob = new SetVelocityFromPressureJobFor()
        {
            positions = m_PredictedPosition,
            densities = m_Density,
            pressures = m_Pressure,
            velocities = m_Velocity,
            deltaTime = deltaTime,
            lookupHelper = m_LookupHelper,
            mass = Mass,
            squaredSmoothingLength = m_SquaredSmoothingLength,
            smoothingLength = SmoothingLength,
            kernelDerivativeTerm = m_KernelDerivativeTerm,
        };
        var velocityFromPressureJobHandle = velocityFromPressureJob.ScheduleParallelByRef(m_Position.Length,
            64, velocityFromGravityJobHandle);

        var viscosityJob = new CalculateViscosityForceJob()
        {
            positions = m_PredictedPosition,
            densities = m_Density,
            velocities = m_Velocity,
            deltaTime = deltaTime,
            lookupHelper = m_LookupHelper,
            mass = Mass,
            smoothingLength = SmoothingLength,
            squaredSmoothingLength = m_SquaredSmoothingLength,
            viscosityFactor = ViscosityFactor,
        };
        var viscosityJobHandle = viscosityJob.ScheduleParallelByRef(m_Position.Length,
            64, velocityFromPressureJobHandle);
        
        var job = new CalculatePositionFromVelocityJobFor()
        {
            positions = m_Position,
            deltaTime = deltaTime,
            width = width,
            height = height,
            velocities = m_Velocity,
            smoothingLength = SmoothingLength,
            boundaryPushStrength = BoundaryPushStrength,
            collisionDamping = CollisionDamping,
        };
        var velocityJobHandle = job.ScheduleParallelByRef(m_Position.Length,
            64, viscosityJobHandle);
        velocityJobHandle.Complete();

        //TODO: this is only needed next frame so this can be backgrounded.
        m_LastStepMaxVelocity = 0;
        for (int i = 0; i < m_ParticleCount; i++)
        {
            m_LastStepMaxVelocity = Mathf.Max(m_LastStepMaxVelocity, m_Velocity[i].magnitude);
        }
    }

    public float CalcTargetDensity()
    {
        return m_ParticleCount / (float)(width * height);
    }

    float CalcCFLTimeStep(float maxVelocity)
    {
        const float courantNumber = 0.3f;
        const float speedOfSound = 10.0f;
        return courantNumber * SmoothingLength / (speedOfSound + maxVelocity);
    }
    
    void CachePrecomputedValues()
    {
        m_SquaredSmoothingLength = SmoothingLength*SmoothingLength;
        m_KernelTerm = SmoothingKernels.CalcSmoothingKernel2Factor(SmoothingLength);
        m_KernelDerivativeTerm = SmoothingKernels.CalcSmoothingKernelDerivativeNormalization(SmoothingLength);
        //m_TargetDensity = CalcTargetDensity();
    }
    
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    static float CalculateDensity(Vector2 pos, GridSpatialLookup lookup, float mass, float squaredSmoothingLength, float smoothingLength, float precalculatedKernelFactor, ReadOnlySpan<Vector2> positions)
    {
        var particleIndices = new NativeList<int>(positions.Length, Allocator.Temp);
        lookup.GetParticlesAround(pos, particleIndices);
        
        var density = 0f;
        //for (var i = 0; i < positions.Length; i++)
        foreach(var i in particleIndices)
        {
            var position = positions[i];
            var sqrDst = Vector2.SqrMagnitude(position - pos);
            if (sqrDst > squaredSmoothingLength) continue;
            //var influence = SmoothingKernel(sqrDst, squaredSmoothingLength, kernelTerm);
            var influence = SmoothingKernels.SmoothingKernel2(Mathf.Sqrt(sqrDst), smoothingLength, precalculatedKernelFactor);
            density += mass * influence;
        }
        particleIndices.Dispose();
        return density;
    }

    [BurstCompile]
    struct CalculateParticlesDensityJobFor : IJobFor
    {
        [ReadOnly]
        public NativeArray<Vector2> positions;
        [ReadOnly]
        public GridSpatialLookup lookupHelper;
        public NativeArray<float> density;

        public float mass;
        public float squaredSmoothingLength;
        public float smoothingLength;
        public float precalculatedKernelFactor; 
        
        public void Execute(int index)
        {
            density[index] = CalculateDensity(positions[index], lookupHelper, mass, squaredSmoothingLength, smoothingLength, precalculatedKernelFactor, positions.AsReadOnlySpan());
        }
    }

    void CalculateParticlesDensity()
    {
        using var markerScope = s_DensityPerfMarker.Auto();
        
        JobHandle dependencyJobHandle = default;
        var job = new CalculateParticlesDensityJobFor()
        {
            positions = m_PredictedPosition,
            lookupHelper = m_LookupHelper,
            mass = Mass,
            squaredSmoothingLength = m_SquaredSmoothingLength,
            smoothingLength = SmoothingLength,
            precalculatedKernelFactor = m_KernelTerm,
            density = m_Density,
        };
        var jobHandle = job.ScheduleParallelByRef(m_Position.Length,
            64, dependencyJobHandle);
        
        jobHandle.Complete();

        // Parallel.For(0, m_ParticleCount, index =>
        // {
        //     m_Density[index] = CalculateDensity(m_Position[index]);
        // });
    }

    void CalculateParticlePressure()
    {
        using var markerScope = s_PressurePerfMarker.Auto();
        
        // const float stiffness = 300;
        // const float adiabaticComponent = 2;
        // for (var index = 0; index < m_ParticleCount; index++)
        // {
        //     m_Pressure[index] = stiffness * (Mathf.Pow(m_Density[index] / m_TargetDensity, adiabaticComponent) - 1);
        // }
        

        for (var index = 0; index < m_ParticleCount; index++)
        {
            var densityError = m_Density[index] - m_TargetDensity;
            var pressure = densityError * m_PressureMultiplier;
            m_Pressure[index] = pressure;
        }
    }
    
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    static Vector2 CalculatePressureGradient(int i, GridSpatialLookup lookup, float mass,float squaredSmoothingLength, float smoothingLength, float kernelDerivativeTerm,
        ReadOnlySpan<Vector2> position, ReadOnlySpan<float> pressure, ReadOnlySpan<float> density)
    {
        var particleIndices = new NativeList<int>(position.Length, Allocator.Temp);
        lookup.GetParticlesAround(position[i], particleIndices);
        var pressureGradient = Vector2.zero;
        //for (var j = 0; j < position.Length; j++)
        foreach(var j in particleIndices)
        {
            if (i == j) continue;
            
            var dif = position[j] - position[i];
            var dir = dif.normalized;
            var sqrDst = Vector2.SqrMagnitude(dif);
            if (sqrDst > squaredSmoothingLength) continue; 
            var distance = Mathf.Sqrt(sqrDst);
            //var influence = SmoothingKernelDerivative(distance, sqrDst, squaredSmoothingLength, kernelDerivativeTerm);
            var influence = SmoothingKernels.SmoothingKernel2Derivative(distance, smoothingLength);
            var averagedPressure = (pressure[j] + pressure[i]) / 2;
            pressureGradient += dir * (averagedPressure * mass) / density[j] * influence;
        }
        particleIndices.Dispose();
        return pressureGradient;
    }
    
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    static Vector2 CalculateViscosityForce(int i, GridSpatialLookup lookup, float mass, float smoothingLength, float squaredSmoothingLength, float viscosityFactor,
        ReadOnlySpan<Vector2> position, ReadOnlySpan<Vector2> velocity, ReadOnlySpan<float> density)
    {
        var particleIndices = new NativeList<int>(position.Length, Allocator.Temp);
        lookup.GetParticlesAround(position[i], particleIndices);
        var viscosity = Vector2.zero;
        //for (var j = 0; j < position.Length; j++)
        foreach(var j in particleIndices)
        {
            if (i == j) continue;
            
            var dif = position[j] - position[i];
            var sqrDst = Vector2.SqrMagnitude(dif);
            if (sqrDst > squaredSmoothingLength) continue;
            
            var velDif = velocity[j] - velocity[i];
            var distance = Mathf.Sqrt(sqrDst);
            //var influence = SmoothingKernelDerivative(distance, sqrDst, squaredSmoothingLength, kernelDerivativeTerm);
            var influence = SmoothingKernels.SmoothingKernel2Derivative(distance, smoothingLength);
            viscosity +=  velDif * viscosityFactor * mass / density[j] * influence;
        }
        particleIndices.Dispose();
        return viscosity;
    }
    

}
