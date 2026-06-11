using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using Unity.Profiling;
using UnityEngine;
using UnityEngine.Rendering;
using Random = UnityEngine.Random;

public class FluidSim : MonoBehaviour, IFluidSim
{
    public int m_ParticleCount;

    public int width;
    public int height;
    public float Gravity = 9.8f;
    // The fluid's rest density. Particle mass and target density derive from
    // this and the particle spacing, so the knobs stay independent of
    // resolution and domain size.
    public float RestDensity = 1f;
    // Kernel support radius in units of particle spacing. ~2 gives ~12
    // neighbors in 2D; raise for smoother fields at higher cost.
    [Range(1.2f, 4f)]
    public float SmoothingRadiusInSpacings = 2f;
    public float m_PressureMultiplier = 1;
    public bool UseAdaptativeStepTime = false;
    [Range(0, 0.3f)]
    public float SimulationStep = 0.0001f;
    public float BoundaryPushStrength = 1;
    public float CollisionDamping = 0.5f;
    public bool AutoStep = true;
    public int MaxStepsPerFrame = 8;
    public float ViscosityFactor = 0.5f;
    public float m_InteractionStrength;
    public ParticlePlacementMode PlacementMode = ParticlePlacementMode.GridWithJitter;
    // Fraction of the box (from the bottom) the fluid occupies at rest. Below 1
    // there is a free surface, so the fluid can slosh and splash.
    [Range(0.1f, 1f)]
    public float FillFraction = 0.6f;

    
    NativeArray<Vector2> m_Position;
    NativeArray<Vector2> m_PredictedPosition;
    NativeArray<float> m_Density;
    NativeArray<float> m_Pressure;
    NativeArray<Vector2> m_Velocity;
    NativeArray<Vector2> m_ViscosityVelocityDelta;
    internal GridSpatialLookup m_LookupHelper;
    
    NativeArray<InputSimulationSurfacePoints> m_ExternalPoints;
    NativeArray<OutputSimulationSurfacePoints> m_ExternalPointsResults;
    NativeArray<Vector2> m_ExternalPointPositions;
    NativeArrayToComputeAdapter<Vector2> m_ExternalPointsBuffer;

    
    

    float m_SquaredSmoothingLength;
    float m_KernelTerm;
    float m_KernelDerivativeTerm;
    float m_LastStepMaxVelocity;
    float m_TimeAccumulator;
    ParticlePlacementMode m_LastPlacementMode;
    float m_LastFillFraction;
    Vector2 m_MousePosition;
    float m_MouseRadius;
    InteractionDirection m_InteractionDirection;
    NativeArrayToComputeAdapter<Vector2> m_PositionComputeBuffer;
    NativeArrayToComputeAdapter<float> m_DensityComputeBuffer;
    NativeArrayToComputeAdapter<float> m_PressureComputeBuffer;
    NativeArrayToComputeAdapter<Vector2> m_VelocitiesComputeBuffer;

    static readonly ProfilerMarker s_UpdatePerfMarker = new ProfilerMarker("FluidSim.Update");
    static readonly ProfilerMarker s_StepPerfMarker = new ProfilerMarker("FluidSim.Step");
    static readonly ProfilerMarker s_PressurePerfMarker = new ProfilerMarker("FluidSim.PressureCalc");
    static readonly ProfilerMarker s_DensityPerfMarker = new ProfilerMarker("FluidSim.DensityeCalc");
    
    public event Action PreSimulation;
    public event Action PostSimulation;
    public GridSpatialLookup LookupHelper => m_LookupHelper;
    public ComputeBuffer GetPositionComputeBuffer() { return m_PositionComputeBuffer.Buffer; }
    public ComputeBuffer GetDensities() { return m_DensityComputeBuffer.Buffer; }
    public ComputeBuffer GetPressures() { return m_PressureComputeBuffer.Buffer; }
    public ComputeBuffer GetVelocities() { return m_VelocitiesComputeBuffer.Buffer; }
    public int ParticleCount => m_ParticleCount;
    public int Height => height;
    public int Width => width;
    public float SmoothingRadius => SmoothingLength;
    float IFluidSim.TargetDensity => RestDensity;
    Transform IFluidSim.Transform => transform;
    public float LastStepDeltaTime { get; private set; }

    // Derived quantities: the spacing comes from how many particles fill the
    // spawn region, and mass follows so that the fill sits exactly at RestDensity.
    public float ParticleSpacing => Mathf.Sqrt(width * height * FillFraction / Mathf.Max(1, m_ParticleCount));
    public float Mass => RestDensity * ParticleSpacing * ParticleSpacing;
    public float SmoothingLength => SmoothingRadiusInSpacings * ParticleSpacing;

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
        if (m_ViscosityVelocityDelta.IsCreated)
            m_ViscosityVelocityDelta.Dispose();
        if(m_PredictedPosition.IsCreated)
            m_PredictedPosition.Dispose();
        
        m_PositionComputeBuffer.Dispose();
        m_DensityComputeBuffer.Dispose();
        m_PressureComputeBuffer.Dispose();
        m_VelocitiesComputeBuffer.Dispose();
        m_ExternalPointsBuffer.Dispose();
        if (m_ExternalPoints.IsCreated)
            m_ExternalPoints.Dispose();
        if (m_ExternalPointsResults.IsCreated)
            m_ExternalPointsResults.Dispose();
        if (m_ExternalPointPositions.IsCreated)
            m_ExternalPointPositions.Dispose();
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
        m_ViscosityVelocityDelta = new NativeArray<Vector2>(m_ParticleCount, Allocator.Persistent);
        m_LastPlacementMode = PlacementMode;
        m_LastFillFraction = FillFraction;
        for ( var i = 0; i < m_ParticleCount; i++ )
        {
            var position = ParticlePlacement.GetPosition(PlacementMode, i, m_ParticleCount, width, height * FillFraction);
            m_Position[i] = position;
            m_PredictedPosition[i] = position;
            m_Density[i] = 0;
            m_Pressure[i] = 0;
            m_Velocity[i] = Vector2.zero;
        }
        
        m_PositionComputeBuffer = new NativeArrayToComputeAdapter<Vector2>(m_Position);
        m_DensityComputeBuffer = new NativeArrayToComputeAdapter<float>(m_Density);
        m_PressureComputeBuffer = new NativeArrayToComputeAdapter<float>(m_Pressure);
        m_VelocitiesComputeBuffer = new NativeArrayToComputeAdapter<Vector2>(m_Velocity);

        // Zero-length until a rigid body registers surface points, so every
        // job can always be scheduled whether or not a body exists.
        m_ExternalPoints = new NativeArray<InputSimulationSurfacePoints>(0, Allocator.Persistent);
        m_ExternalPointsResults = new NativeArray<OutputSimulationSurfacePoints>(0, Allocator.Persistent);
    }
    
    void DoReInitializationIfNecessary()
    {
        var particleCountDifferent = m_ParticleCount != m_Position.Length;
        var smoothingRadiusDifferent = !m_LookupHelper.IsValid || !Mathf.Approximately(SmoothingLength, m_LookupHelper.CellSize);
        if (particleCountDifferent || PlacementMode != m_LastPlacementMode || !Mathf.Approximately(FillFraction, m_LastFillFraction))
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
            AdvanceTime(Time.deltaTime);
    }

    // Advance the simulation by frameTime of wall-clock time, taking fixed-size
    // substeps so simulation speed is independent of framerate.
    public void AdvanceTime(float frameTime)
    {
        using var markerScope = s_UpdatePerfMarker.Auto();

        if (SimulationStep <= 0)
            return;

        CachePrecomputedValues();
        m_TimeAccumulator += frameTime;
        var steps = 0;
        while (m_TimeAccumulator >= SimulationStep && steps < MaxStepsPerFrame)
        {
            var stepTime = CurrentStepTime();
            StepSimulation(stepTime);
            m_TimeAccumulator -= stepTime;
            steps++;
        }
        // Running behind realtime: drop the surplus so the sim slows down gracefully
        // instead of accumulating an ever-growing debt of steps.
        if (steps == MaxStepsPerFrame)
            m_TimeAccumulator = 0;

        m_MouseRadius = 0; //clear mouse interaction
        UpdateComputeBuffers();
    }

    // Advance exactly one substep, ignoring wall-clock time (editor step button).
    public void DoUpdate()
    {
        using var markerScope = s_UpdatePerfMarker.Auto();

        CachePrecomputedValues();
        StepSimulation(CurrentStepTime());
        m_MouseRadius = 0; //clear mouse interaction
        UpdateComputeBuffers();
    }

    float CurrentStepTime()
    {
        var stepTime = SimulationStep;
        if (UseAdaptativeStepTime)
        {
            stepTime = Mathf.Min(SimulationStep, CalcCFLTimeStep(m_LastStepMaxVelocity));
        }
        return stepTime;
    }

    void StepSimulation(float stepTime)
    {
        LastStepDeltaTime = stepTime;

        // Fired per substep: the rigid-body coupling registers fresh surface
        // data before the step and retrieves the exchanged impulses after it.
        // Each step overwrites the impulse outputs, so retrieving once per
        // frame would drop the other substeps' momentum; and a frame with zero
        // substeps correctly exchanges nothing.
        PreSimulation?.Invoke();

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

        PostSimulation?.Invoke();
    }

    void UpdateComputeBuffers()
    {
        m_PositionComputeBuffer.Update(m_Position);
        m_DensityComputeBuffer .Update(m_Density);
        m_PressureComputeBuffer.Update(m_Pressure);
        m_VelocitiesComputeBuffer.Update(m_Velocity);

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
    
    //[BurstCompile]
    struct CalculateAccelerationFromExternalForcesJob : IJobFor
    {
        [ReadOnly]
        public NativeArray<Vector2> positions;
        public NativeArray<Vector2> velocities;
        [ReadOnly]
        public NativeArray<InputSimulationSurfacePoints> rigidBodyPoints;
        [ReadOnly]
        public NativeArray<OutputSimulationSurfacePoints> outputSurfacePoints;
        [ReadOnly]
        public NativeArray<float> densities;
        public float squaredSmoothingLength;
        public float smoothingLength;
        public float precalculatedKernelFactor;
        public float mass;
        public float contactDistance;

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

            var reactionVelocityChange = ReactionVelocityChangeFromRigidBody(position, densities[index],
                mass, squaredSmoothingLength, smoothingLength, precalculatedKernelFactor,
                rigidBodyPoints, outputSurfacePoints);

            // The reaction is an impulse-derived velocity change, not a force:
            // it must not be scaled by deltaTime again.
            var newVelocity = velocity + (Vector2.down * gravity + interactionForce) * deltaTime + reactionVelocityChange;
            newVelocity = CancelVelocityIntoSurface(position, newVelocity, contactDistance, rigidBodyPoints);
            velocities[index] = newVelocity;
        }

        // Last line of defense against tunneling: the boundary-density barrier
        // is soft, so a fast particle can still cross the hull within one
        // step. Within contact distance of a sample, remove the velocity
        // component moving into the surface, keeping tangential motion (free
        // slip) — the body-surface analog of the wall bounce.
        static Vector2 CancelVelocityIntoSurface(Vector2 particlePosition, Vector2 velocity, float contactDistance,
            NativeArray<InputSimulationSurfacePoints> rigidBodyPoints)
        {
            var sqrContact = contactDistance * contactDistance;
            for (var i = 0; i < rigidBodyPoints.Length; i++)
            {
                var p = rigidBodyPoints[i];
                if (Vector2.SqrMagnitude(particlePosition - p.SimSpacePoint) > sqrContact) continue;
                var relativeNormalSpeed = Vector2.Dot(velocity - p.velocity, p.normal);
                if (relativeNormalSpeed < 0)
                    velocity -= relativeNormalSpeed * p.normal;
            }
            return velocity;
        }

        // Newton's third law. Each surface point received an impulse interpolated
        // from nearby particles with weights w = (m/rho) * W, normalized by the
        // weight sum stored alongside it. Give every particle back the equal and
        // opposite share of that impulse: the shares w/weightSum sum to 1 over
        // the contributing particles, so the momentum handed to the fluid is
        // exactly the negative of what the body received. Conservation holds by
        // bookkeeping; no tuning constant.
        static Vector2 ReactionVelocityChangeFromRigidBody(Vector2 particlePosition, float particleDensity,
            float mass, float squaredSmoothingLength, float smoothingLength, float kernelFactor,
            NativeArray<InputSimulationSurfacePoints> rigidBodyPoints,
            NativeArray<OutputSimulationSurfacePoints> surfacePointImpulses)
        {
            var velocityChange = Vector2.zero;
            for (var i = 0; i < rigidBodyPoints.Length; i++)
            {
                var output = surfacePointImpulses[i];
                if (output.weightSum <= 1e-6f) continue;

                var dif = particlePosition - rigidBodyPoints[i].SimSpacePoint;
                var sqrDst = dif.sqrMagnitude;
                if (sqrDst > squaredSmoothingLength) continue;

                var w = mass / particleDensity * SmoothingKernels.SmoothingKernel2(Mathf.Sqrt(sqrDst), smoothingLength, kernelFactor);
                var share = w / output.weightSum;
                velocityChange += -output.impulse * share / mass;
            }
            return velocityChange;
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
            var distToTopWall = Mathf.Abs(height - position.y);
            if (distToTopWall < smoothingLength)
            {
                var strength = (smoothingLength - distToTopWall) / smoothingLength;
                velocity.y -= boundaryPushStrength * strength * deltaTime;
                if (position.y > height)
                {
                    position.y = height - (position.y - height);
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
    
    // Viscosity is split into two passes: this one reads the velocity field
    // read-only and writes per-particle deltas, and ApplyVelocityDeltasJob adds
    // them afterwards. Writing velocities in place here would race with the
    // neighbor reads happening concurrently on other threads.
    [BurstCompile]
    struct CalculateViscosityForceJob : IJobFor
    {
        [ReadOnly]
        public NativeArray<Vector2> velocities;
        [ReadOnly]
        public NativeArray<Vector2> positions;
        [ReadOnly]
        public NativeArray<float> densities;
        [ReadOnly]
        public GridSpatialLookup lookupHelper;
        [WriteOnly]
        public NativeArray<Vector2> velocityDeltas;
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
            velocityDeltas[index] = viscosityAcceleration * deltaTime;
        }
    }

    [BurstCompile]
    struct ApplyVelocityDeltasJob : IJobFor
    {
        public NativeArray<Vector2> velocities;
        [ReadOnly]
        public NativeArray<Vector2> velocityDeltas;
        public void Execute(int index)
        {
            velocities[index] += velocityDeltas[index];
        }
    }

    //[BurstCompile]
    struct CalculateForcesToExternalObjects : IJobFor
    {
        [ReadOnly]
        public NativeArray<Vector2> positions;
        [ReadOnly]
        public GridSpatialLookup lookupHelper;
        [ReadOnly]
        public NativeArray<Vector2> velocities;
        [ReadOnly]
        public NativeArray<float> pressures;
        [ReadOnly]
        public NativeArray<float> densities;
        [ReadOnly]
        public NativeArray<InputSimulationSurfacePoints> rigidBodyPoints;
        public NativeArray<OutputSimulationSurfacePoints> outputSurfacePoints;
        public float mass;
        public float squaredSmoothingLength;
        public float smoothingLength;
        public float precalculatedKernelFactor;
        public float deltaTime;
        public float viscosityCoefficient;

        public void Execute(int index)
        {
            var surfacePoint = rigidBodyPoints[index];

            var fluidForce = CalculateFluidForceAtSurfacePoint(
                surfacePoint.SimSpacePoint, surfacePoint.normal, surfacePoint.velocity, lookupHelper,
                mass, squaredSmoothingLength, smoothingLength, precalculatedKernelFactor, viscosityCoefficient,
                positions, pressures, densities, velocities, out var weightSum);
            // Integrating the force over this step's dt up front makes the
            // output an impulse: the consumer no longer needs to know (or
            // guess) which timestep the simulation used.
            var impulse = fluidForce * surfacePoint.areaWeight * deltaTime;
            outputSurfacePoints[index] = new OutputSimulationSurfacePoints
            {
                impulse = impulse,
                weightSum = weightSum
            };
        }
    }

    void SimulateStep(float deltaTime)
    {
        using var markerScope = s_StepPerfMarker.Auto();
        
        JobHandle dependencyJobHandle = default;

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
            64, dependencyJobHandle);

        var viscosityJob = new CalculateViscosityForceJob()
        {
            positions = m_PredictedPosition,
            densities = m_Density,
            velocities = m_Velocity,
            velocityDeltas = m_ViscosityVelocityDelta,
            deltaTime = deltaTime,
            lookupHelper = m_LookupHelper,
            mass = Mass,
            smoothingLength = SmoothingLength,
            squaredSmoothingLength = m_SquaredSmoothingLength,
            viscosityFactor = ViscosityFactor,
        };
        var viscosityJobHandle = viscosityJob.ScheduleParallelByRef(m_Position.Length,
            64, velocityFromPressureJobHandle);

        var applyViscosityJob = new ApplyVelocityDeltasJob()
        {
            velocities = m_Velocity,
            velocityDeltas = m_ViscosityVelocityDelta,
        };
        var applyViscosityJobHandle = applyViscosityJob.ScheduleParallelByRef(m_Position.Length,
            64, viscosityJobHandle);

        // The surface-force job reads velocities for drag, so it runs after the
        // viscosity deltas are applied; the gravity/reaction job consumes its
        // outputs and runs after it.
        var externalObjectForcesJob = new CalculateForcesToExternalObjects
        {
            positions = m_PredictedPosition,
            lookupHelper = m_LookupHelper,
            velocities = m_Velocity,
            pressures = m_Pressure,
            densities = m_Density,
            rigidBodyPoints = m_ExternalPoints,
            outputSurfacePoints = m_ExternalPointsResults,
            mass = Mass,
            squaredSmoothingLength = m_SquaredSmoothingLength,
            smoothingLength = SmoothingLength,
            precalculatedKernelFactor = m_KernelTerm,
            viscosityCoefficient = ViscosityFactor,
            deltaTime = deltaTime,

        };
        var externalObjectForcesJobHandle = externalObjectForcesJob.ScheduleParallelByRef(m_ExternalPoints.Length,
            64, applyViscosityJobHandle);
        externalObjectForcesJobHandle.Complete();
        
        var strengthSign = m_InteractionDirection == InteractionDirection.Attract ? 1 : -1;
        var externalForcesJob = new CalculateAccelerationFromExternalForcesJob
        {
            positions = m_PredictedPosition,
            velocities = m_Velocity,
            rigidBodyPoints = m_ExternalPoints,
            outputSurfacePoints = m_ExternalPointsResults,
            densities = m_Density,
            mass = Mass,
            precalculatedKernelFactor = m_KernelTerm,
            squaredSmoothingLength = m_SquaredSmoothingLength,
            smoothingLength = SmoothingLength,
            contactDistance = SmoothingLength * 0.25f,
            gravity = Gravity,
            deltaTime = deltaTime,
            forceCenter = m_MousePosition,
            forceRadius = m_MouseRadius,
            forceStrength = m_InteractionStrength * strengthSign,

        };
        var externalForcesJobHandle = externalForcesJob.ScheduleParallelByRef(m_Velocity.Length,
            64, externalObjectForcesJobHandle);

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
            64, externalForcesJobHandle);
        velocityJobHandle.Complete();

        //TODO: this is only needed next frame so this can be backgrounded.
        m_LastStepMaxVelocity = 0;
        for (int i = 0; i < m_ParticleCount; i++)
        {
            m_LastStepMaxVelocity = Mathf.Max(m_LastStepMaxVelocity, m_Velocity[i].magnitude);
        }
    }

    float CalcCFLTimeStep(float maxVelocity)
    {
        const float courantNumber = 0.3f;
        // For the linear EOS p = k * (density - rest), the speed of sound is sqrt(k).
        var speedOfSound = Mathf.Sqrt(Mathf.Max(1e-6f, m_PressureMultiplier));
        return courantNumber * SmoothingLength / (speedOfSound + maxVelocity);
    }

    void CachePrecomputedValues()
    {
        m_SquaredSmoothingLength = SmoothingLength*SmoothingLength;
        m_KernelTerm = SmoothingKernels.CalcSmoothingKernel2Factor(SmoothingLength);
        m_KernelDerivativeTerm = SmoothingKernels.CalcSmoothingKernelDerivativeNormalization(SmoothingLength);
    }
    
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    static float CalculateDensity(Vector2 pos, GridSpatialLookup lookup, float mass, float squaredSmoothingLength, float smoothingLength, float precalculatedKernelFactor, ReadOnlySpan<Vector2> positions,
        ReadOnlySpan<InputSimulationSurfacePoints> boundaryPoints)
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

        // The rigid body's surface samples contribute their pseudo-mass exactly
        // like fluid neighbors. A particle approaching the hull reads rising
        // density instead of a void, so the standard pressure repulsion keeps
        // fluid out of the body, and pressure sampled at the hull becomes
        // hydrostatically correct (real buoyancy magnitudes).
        for (var b = 0; b < boundaryPoints.Length; b++)
        {
            var sqrDst = Vector2.SqrMagnitude(boundaryPoints[b].SimSpacePoint - pos);
            if (sqrDst > squaredSmoothingLength) continue;
            density += boundaryPoints[b].pseudoMass * SmoothingKernels.SmoothingKernel2(Mathf.Sqrt(sqrDst), smoothingLength, precalculatedKernelFactor);
        }
        return density;
    }

    [BurstCompile]
    struct CalculateParticlesDensityJobFor : IJobFor
    {
        [ReadOnly]
        public NativeArray<Vector2> positions;
        [ReadOnly]
        public GridSpatialLookup lookupHelper;
        [ReadOnly]
        public NativeArray<InputSimulationSurfacePoints> rigidBodyPoints;
        public NativeArray<float> density;

        public float mass;
        public float squaredSmoothingLength;
        public float smoothingLength;
        public float precalculatedKernelFactor;

        public void Execute(int index)
        {
            density[index] = CalculateDensity(positions[index], lookupHelper, mass, squaredSmoothingLength, smoothingLength, precalculatedKernelFactor, positions.AsReadOnlySpan(),
                rigidBodyPoints.AsReadOnlySpan());
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
            rigidBodyPoints = m_ExternalPoints,
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
            // Clamp to non-negative: below-target density (free surfaces, wall-truncated
            // kernels) would otherwise produce attractive pressure, causing particle
            // clumping and sticking to boundaries. For the rigid body this matters
            // doubly — the boundary-density barrier only works if pressure rises
            // with proximity rather than oscillating around zero.
            var densityError = m_Density[index] - RestDensity;
            var pressure = Mathf.Max(0, densityError * m_PressureMultiplier);
            m_Pressure[index] = pressure;
        }
    }
    
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    static Vector2 CalculatePressureGradient(int i, GridSpatialLookup lookup, float mass,float squaredSmoothingLength, float smoothingLength, float kernelDerivativeTerm,
        ReadOnlySpan<Vector2> positions, ReadOnlySpan<float> pressure, ReadOnlySpan<float> density)
    {
        var particleIndices = new NativeList<int>(positions.Length, Allocator.Temp);
        lookup.GetParticlesAround(positions[i], particleIndices);
        var pressureGradient = Vector2.zero;
        //for (var j = 0; j < position.Length; j++)
        foreach(var j in particleIndices)
        {
            if (i == j) continue;
            
            var dif = positions[j] - positions[i];
            var sqrDst = Vector2.SqrMagnitude(dif);
            if (sqrDst > squaredSmoothingLength) continue;
            
            var dir = dif.normalized;
            var distance = Mathf.Sqrt(sqrDst);
            //var influence = SmoothingKernelDerivative(distance, sqrDst, squaredSmoothingLength, kernelDerivativeTerm);
            //TODO: this should be the laplacian of the Kernel.
            var influence = SmoothingKernels.SmoothingKernel2Derivative(distance, smoothingLength);
            var averagedPressure = (pressure[j] + pressure[i]) / 2;
            pressureGradient += dir * (averagedPressure * mass) / density[j] * influence;
        }
        particleIndices.Dispose();
        return pressureGradient;
    }
    
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    // Samples the fluid at a rigid-body surface point and returns the force per
    // unit of surface length the fluid exerts there. Pressure and fluid
    // velocity are SPH-interpolated with Shepard normalization (divide by the
    // weight sum), so the estimate is unbiased even where particle coverage is
    // partial. Pressure pushes along the inward surface normal (F = -p * n);
    // drag opposes the surface's motion relative to the fluid.
    static Vector2 CalculateFluidForceAtSurfacePoint(Vector2 position, Vector2 normal, Vector2 surfaceVelocity,
        GridSpatialLookup lookup, float mass,
        float squaredSmoothingLength, float smoothingLength, float kernelFactor, float dragCoefficient,
        ReadOnlySpan<Vector2> positions, ReadOnlySpan<float> pressure, ReadOnlySpan<float> density, ReadOnlySpan<Vector2> velocities,
        out float weightSum)
    {
        var particleIndices = new NativeList<int>(positions.Length, Allocator.Temp);
        lookup.GetParticlesAround(position, particleIndices);
        weightSum = 0f;
        var pressureSum = 0f;
        var velocitySum = Vector2.zero;
        foreach(var j in particleIndices)
        {
            var sqrDst = Vector2.SqrMagnitude(positions[j] - position);
            if (sqrDst > squaredSmoothingLength) continue;

            var w = mass / density[j] * SmoothingKernels.SmoothingKernel2(Mathf.Sqrt(sqrDst), smoothingLength, kernelFactor);
            weightSum += w;
            // Clamped to non-negative: negative pressure (density deficiency)
            // would read as suction pulling the body into the fluid.
            pressureSum += Mathf.Max(0f, pressure[j]) * w;
            velocitySum += velocities[j] * w;
        }
        particleIndices.Dispose();

        // No fluid within kernel range: no force (and no division by zero).
        if (weightSum <= 1e-6f) return Vector2.zero;

        var interpolatedPressure = pressureSum / weightSum;
        var interpolatedVelocity = velocitySum / weightSum;
        var pressureForce = -normal * interpolatedPressure;
        var dragForce = (interpolatedVelocity - surfaceVelocity) * dragCoefficient;
        return pressureForce + dragForce;
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
            
            // Same conventions as main: velDif = v_j - v_i and a positive
            // kernel weight (negated spiky derivative), so the force pulls
            // velocities together. NOTE: this branch originally had BOTH signs
            // flipped (v_i - v_j against the raw negative derivative), which
            // cancelled out to correct damping — negating only one of them, as
            // an earlier "backport" here did, turns viscosity into energy
            // injection. Sign fixes are relative to local conventions.
            var velDif = velocity[j] - velocity[i];
            var distance = Mathf.Sqrt(sqrDst);
            // Viscosity must pull velocities together (force ∝ +(v_j - v_i)), so it needs a
            // positive kernel weight; the spiky derivative is negative inside the radius.
            var influence = -SmoothingKernels.SmoothingKernel2Derivative(distance, smoothingLength);
            viscosity +=  velDif * mass / density[j] * influence;
        }
        particleIndices.Dispose();
        return viscosity * viscosityFactor;
    }
    
    //TODO: this only handles a single rigidbody
    public void SetRigidBodySurfaceResults(IList<InputSimulationSurfacePoints> points)
    {
        if (!m_ExternalPoints.IsCreated || m_ExternalPoints.Length < points.Count)
        {
            if (m_ExternalPoints.IsCreated)
                m_ExternalPoints.Dispose();
            if (m_ExternalPointsResults.IsCreated)
                m_ExternalPointsResults.Dispose();
            if (m_ExternalPointPositions.IsCreated)
                m_ExternalPointPositions.Dispose();
            // Release the previous adapter's ComputeBuffer before replacing it.
            m_ExternalPointsBuffer.Dispose();

            m_ExternalPoints = new NativeArray<InputSimulationSurfacePoints>(points.Count, Allocator.Persistent);
            m_ExternalPointsResults = new NativeArray<OutputSimulationSurfacePoints>(points.Count, Allocator.Persistent);
            m_ExternalPointPositions = new NativeArray<Vector2>(points.Count, Allocator.Persistent);
            m_ExternalPointsBuffer = new NativeArrayToComputeAdapter<Vector2>(m_ExternalPointPositions);
        }
        for (var i = 0; i < points.Count; i++)
        {
            var point = points[i];
            m_ExternalPoints[i] = point;
            m_ExternalPointPositions[i] = point.SimSpacePoint;
        }
        m_ExternalPointsBuffer.Update(m_ExternalPointPositions);
    }
    public void RetrieveRigidBodySurfaceResults(IList<OutputSimulationSurfacePoints> points)
    {
        // The native array only ever grows, so it can be larger than the
        // caller's list; copy only what the caller asked for.
        var count = Mathf.Min(m_ExternalPointsResults.Length, points.Count);
        for (var i = 0; i < count; i++)
        {
            points[i] = m_ExternalPointsResults[i];
        }
    }
    public ComputeBuffer InputExternalPoints => m_ExternalPointsBuffer.Buffer;
}
