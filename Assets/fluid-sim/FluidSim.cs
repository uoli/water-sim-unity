using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;
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
    public float Mass = 1;
    public float SmoothingLength = 10;
    [Range(0, 1)]
    public float SimulationStep = 0.0001f;
    public bool AutoStep = true;
    
    NativeArray<Vector2> m_Position;
    NativeArray<float> m_Density;
    NativeArray<float> m_Pressure;
    NativeArray<Vector2> m_Velocity;

    float m_SquaredSmoothingLength;
    float m_KernelTerm;
    float m_KernelDerivativeTerm;
    float m_TargetDensity;
    
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
    }

    void InitParticles()
    {
        CleanupParticles();
        m_Position = new NativeArray<Vector2>(m_ParticleCount, Allocator.Persistent);
        m_Density = new NativeArray<float>(m_ParticleCount, Allocator.Persistent);
        m_Pressure = new NativeArray<float>(m_ParticleCount, Allocator.Persistent);
        m_Velocity = new NativeArray<Vector2>(m_ParticleCount, Allocator.Persistent);
        for ( var i = 0; i < m_ParticleCount; i++ )
        {
            var position = new Vector2(Random.Range(0f, width), Random.Range(0f, height));
            m_Position[i] = position;
            m_Density[i] = 0;
            m_Pressure[i] = 0;
            m_Velocity[i] = Vector2.zero;
        }
    }

    // Update is called once per frame
    void Update()
    {
        if (m_ParticleCount != m_Position.Length)
        {
            InitParticles();
        }
        
        if (AutoStep)
            DoUpdate();
    }

    public void DoUpdate()
    {
        using var markerScope = s_UpdatePerfMarker.Auto();
        
        CachePrecomputedValues();
        CalculateParticlesDensity();
        CalculateParticlePressure();
        SimulateStep(SimulationStep);
    }
    
    struct CalculatePositionFromVelocityJobFor : IJobFor
    {
        public NativeArray<Vector2> positions;
        public float deltaTime;
        public float width;
        public float height;
        
        [ReadOnly]
        public NativeArray<Vector2> velocities;
        public void Execute(int index)
        {
            var position = positions[index];
            position += velocities[index] * deltaTime;
            if (position.x < 0)
                position.x = 0;
            if (position.y < 0)
                position.y = 0;
            if (position.x > width)
                position.x = width;
            if (position.y > height)
                position.y = height;
            positions[index] = position;
        }
    }
    
    struct SetVelocityFromPressureJobFor : IJobFor
    {
        [ReadOnly]
        public NativeArray<Vector2> positions;
        [ReadOnly]
        public NativeArray<float> densities;
        [ReadOnly]
        public NativeArray<float> pressures;
        
        public NativeArray<Vector2> velocities;
        public float mass;
        public float squaredSmoothingLength;
        public float smoothingLength;
        public float kernelDerivativeTerm;
        public float deltaTime;
        public void Execute(int index)
        {
            var position = positions[index];
            var density = densities[index]; 
            var pressureForce = CalculatePressureGradient(position, mass, squaredSmoothingLength, smoothingLength, kernelDerivativeTerm,
                positions, pressures, densities);
            var pressureAcceleration = pressureForce / density;
            var velocity = pressureAcceleration * deltaTime;
            velocities[index] = velocity;
        }
    }

    void SimulateStep(float deltaTime)
    {
        using var markerScope = s_StepPerfMarker.Auto();
        
        JobHandle dependencyJobHandle = default;
        var velocityFromPressureJob = new SetVelocityFromPressureJobFor()
        {
            positions = m_Position,
            densities = m_Density,
            pressures = m_Pressure,
            velocities = m_Velocity,
            deltaTime = deltaTime,
            mass = Mass,
            squaredSmoothingLength = m_SquaredSmoothingLength,
            smoothingLength = SmoothingLength,
            kernelDerivativeTerm = m_KernelDerivativeTerm,
        };
        var velocityFromPressureJobHandle = velocityFromPressureJob.ScheduleParallelByRef(m_Position.Length,
            64, dependencyJobHandle);
        
        var job = new CalculatePositionFromVelocityJobFor()
        {
            positions = m_Position,
            deltaTime = deltaTime,
            width = width,
            height = height,
            velocities = m_Velocity
        };
        var velocityJobHandle = job.ScheduleParallelByRef(m_Position.Length,
            64, velocityFromPressureJobHandle);
        velocityJobHandle.Complete();
    }

    public float CalcTargetDensity()
    {
        return m_ParticleCount / (float)(width * height);
    }
    
    void CachePrecomputedValues()
    {
        m_SquaredSmoothingLength = SmoothingLength*SmoothingLength;
        m_KernelTerm = 4 / (Mathf.PI * Mathf.Pow(SmoothingLength, 8));
        m_KernelDerivativeTerm = -24 / (Mathf.PI * Mathf.Pow(SmoothingLength, 8));
        m_TargetDensity = CalcTargetDensity();
    }

    float CalculateDensity(Vector2 pos)
    {
        return CalculateDensity(pos, Mass, m_SquaredSmoothingLength, SmoothingLength, m_KernelTerm, m_Position.AsReadOnlySpan());
    }
    
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    static float CalculateDensity(Vector2 pos, float mass, float squaredSmoothingLength, float smoothingLength, float kernelTerm, ReadOnlySpan<Vector2> positions)
    {
        var density = 0f;
        for (var i = 0; i < positions.Length; i++)
        {
            var position = positions[i];
            var sqrDst = Vector2.SqrMagnitude(position - pos);
            //var influence = SmoothingKernel(sqrDst, squaredSmoothingLength, kernelTerm);
            var influence = SmoothingKernel2(Mathf.Sqrt(sqrDst), smoothingLength);
            density += mass * influence;
        }

        return density;
    }

    struct CalculateParticlesDensityJobFor : IJobFor
    {
        [ReadOnly]
        public NativeArray<Vector2> positions;
        public float mass;
        public float squaredSmoothingLength;
        public float smoothingLength;
        public float kernelTerm;
        
        [WriteOnly]
        public NativeArray<float> density;
        public void Execute(int index)
        {
            density[index] = CalculateDensity(positions[index], mass, squaredSmoothingLength, smoothingLength, kernelTerm, positions.AsReadOnlySpan());
        }
    }

    void CalculateParticlesDensity()
    {
        using var markerScope = s_DensityPerfMarker.Auto();
        
        JobHandle dependencyJobHandle = default;
        var job = new CalculateParticlesDensityJobFor()
        {
            positions = m_Position,
            mass = Mass,
            squaredSmoothingLength = m_SquaredSmoothingLength,
            smoothingLength = SmoothingLength,
            kernelTerm = m_KernelTerm,
            density = m_Density,
        };
        var velocityJobHandle = job.ScheduleParallelByRef(m_Position.Length,
            64, dependencyJobHandle);
        
        velocityJobHandle.Complete();

        // Parallel.For(0, m_ParticleCount, index =>
        // {
        //     m_Density[index] = CalculateDensity(m_Position[index]);
        // });
    }

    void CalculateParticlePressure()
    {
        using var markerScope = s_PressurePerfMarker.Auto();
        
        const float stiffness = 3000;
        const float adiabaticComponent = 7;
        //Parallel.For(0, m_ParticleCount, index =>
        for (var index = 0; index < m_ParticleCount; index++)
        {
            m_Pressure[index] = stiffness * (Mathf.Pow(m_Density[index] / m_TargetDensity, adiabaticComponent) - 1);
        }

        //});
    }

    public float CalculatePressure(Vector2 pos)
    {
        var pressure = 0f;
        for (var i = 0; i < m_ParticleCount; i++)
        {
            
            var sqrDst = Vector2.SqrMagnitude(m_Position[i] - pos);
            var influence = SmoothingKernel2(Mathf.Sqrt(sqrDst), SmoothingLength);
            //var influence = SmoothingKernel(sqrDst);
            //var density = CalculateDensity(pos);
            pressure += m_Pressure[i] * Mass / m_Density[i] * influence;
        }

        return pressure;
    }

    public Vector2 CalculatePressureGradient(Vector2 pos)
    {
        return CalculatePressureGradient(pos, Mass, m_SquaredSmoothingLength, SmoothingLength, m_KernelDerivativeTerm, m_Position, m_Pressure, m_Density);
    }
    
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    static public Vector2 CalculatePressureGradient(Vector2 pos, float mass,float squaredSmoothingLength, float smoothingLength, float kernelDerivativeTerm,
        ReadOnlySpan<Vector2> position, ReadOnlySpan<float> pressure, ReadOnlySpan<float> density)
    {
        var pressureGradient = Vector2.zero;
        for (var i = 0; i < position.Length; i++)
        {
            var dif = position[i] - pos;
            var dir = dif.normalized;
            var sqrDst = Vector2.SqrMagnitude(dif);
            var distance = Mathf.Sqrt(sqrDst);
            //var influence = SmoothingKernelDerivative(distance, sqrDst, squaredSmoothingLength, kernelDerivativeTerm);
            var influence = SmoothingKernel2Derivative(distance, smoothingLength);
            pressureGradient += dir * (pressure[i] * mass) / density[i] * influence;
        }

        return pressureGradient;
    }

    static float SmoothingKernel2(float dst, float radius)
    {
        if (dst > radius) return 0;

        float volume = (Mathf.PI * Mathf.Pow(radius, 4)) / 6;
        return (radius - dst) * (radius - dst) / volume;
    }
    
    static float SmoothingKernel2Derivative(float dst, float radius)
    {
        if (dst > radius) return 0;

        float scale = 12 / (Mathf.Pow(radius, 4) * Mathf.PI);
        return (dst - radius ) * scale;
    }

    float SmoothingKernelDerivative(float dst, float sqrdDistance)
    {
        return SmoothingKernelDerivative(dst, sqrdDistance, m_SquaredSmoothingLength, m_KernelDerivativeTerm);
    }
    
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    static float SmoothingKernelDerivative(float dst, float sqrdDistance, float squaredSmoothingLength, float kernelDerivativeTerm)
    {
        if (sqrdDistance > squaredSmoothingLength) return 0;
        var diff = squaredSmoothingLength - sqrdDistance;
        
        return diff*diff *dst * kernelDerivativeTerm;
    }

    float SmoothingKernel(float sqrdDistance)
    {
        return SmoothingKernel(sqrdDistance, m_SquaredSmoothingLength, m_KernelTerm);
    }
    
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    static float SmoothingKernel(float sqrdDistance, float squaredSmoothingLength, float kernelTerm)
    {
        if (sqrdDistance > squaredSmoothingLength) return 0;
        //var term = 4 / (Mathf.PI * Mathf.Pow(smoothingFactor, 8));
        var diff = squaredSmoothingLength - sqrdDistance;
        //return 315f / (64f * Mathf.PI * Mathf.Pow(smoothingFactor,9)) * Mathf.Pow(smoothingFactor*smoothingFactor - distance*distance, 3);
        return diff*diff*diff * kernelTerm;
    }
}
