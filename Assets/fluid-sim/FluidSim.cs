using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using Unity.Collections;
using Unity.Profiling;
using UnityEngine;
using Random = UnityEngine.Random;

[Serializable]
public struct Particle //probably an Struct of arrays could be beneficial
{
    public Vector2 position;
    public float density;
    public float pressure;
    public Vector2 velocity;
}

public class FluidSim : MonoBehaviour
{
    public int m_ParticleCount;

    public int width;
    public int height;
    public float Mass = 1;
    public float SmoothingLength = 10;
    [Range(0, 1)]
    public float SimulationStep = 0.0001f;
    
    Vector2[] m_Position;
    NativeArray<float> m_Density;
    NativeArray<float> m_Pressure;
    Vector2[] m_Velocity;

    float m_SquaredSmoothingLength;
    float m_KernelTerm;
    float m_KernelDerivativeTerm;
    float m_TargetDensity;
    
    static readonly ProfilerMarker s_UpdatePerfMarker = new ProfilerMarker("FluidSim.Update");
    static readonly ProfilerMarker s_StepPerfMarker = new ProfilerMarker("FluidSim.Step");
    static readonly ProfilerMarker s_PressurePerfMarker = new ProfilerMarker("FluidSim.PressureCalc");
    static readonly ProfilerMarker s_DensityPerfMarker = new ProfilerMarker("FluidSim.DensityeCalc");

    public IReadOnlyList<Vector2> GetPositions() { return m_Position; }
    public NativeArray<float>.ReadOnly GetDensities() { return m_Density.AsReadOnly(); }
    public NativeArray<float>.ReadOnly GetPressures() { return m_Pressure.AsReadOnly(); }
    public IReadOnlyList<Vector2> GetVelocities() { return m_Velocity; }
    
    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        InitParticles();
    }

    void InitParticles()
    {
        m_Position = new Vector2[m_ParticleCount];
        m_Density = new NativeArray<float>(m_ParticleCount, Allocator.Persistent);
        m_Pressure = new NativeArray<float>(m_ParticleCount, Allocator.Persistent);
        m_Velocity = new Vector2[m_ParticleCount];
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
        using var markerScope = s_UpdatePerfMarker.Auto();

        if (m_ParticleCount != m_Position.Length)
        {
            InitParticles();
        }

        CachePrecomputedValues();
        CalculateParticlesDensity();
        CalculateParticlePressure();
        SimulateStep(SimulationStep);
    }

    void SimulateStep(float deltaTime)
    {
        using var markerScope = s_StepPerfMarker.Auto();

        Parallel.For(0, m_ParticleCount, index =>
        {
            var position = m_Position[index];
            var density = m_Density[index];
            var pressureForce = CalculatePressureGradient(position);
            var pressureAcceleration = pressureForce / density;
            var velocity = pressureAcceleration * deltaTime;
            m_Velocity[index] += velocity;
        });

        Parallel.For(0, m_ParticleCount, index =>
        {
            var position = m_Position[index];
            position += m_Velocity[index] * deltaTime;
            if (position.x < 0)
                position.x = 0;
            if (position.y < 0)
                position.y = 0;
            if (position.x > width)
                position.x = width;
            if (position.y > height)
                position.y = height;
            m_Position[index] = position;
        });

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
        var density = 0f;
        for (var i = 0; i < m_ParticleCount; i++)
        {
            var position = m_Position[i];
            var sqrDst = Vector2.SqrMagnitude(position - pos);
            var influence = SmoothingKernel(sqrDst, SmoothingLength);
            density += Mass * influence;
        }

        return density;
    }

    void CalculateParticlesDensity()
    {
        using var markerScope = s_DensityPerfMarker.Auto();

        Parallel.For(0, m_ParticleCount, index =>
        //for (var index = 0; index < m_ParticleCount; index++)
        {
            m_Density[index] = CalculateDensity(m_Position[index]);
        });
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
            var influence = SmoothingKernel(sqrDst, SmoothingLength);
            //var density = CalculateDensity(pos);
            pressure += m_Pressure[i] * Mass / m_Density[i] * influence;
        }

        return pressure;
    }

    public Vector2 CalculatePressureGradient(Vector2 pos)
    {
        var pressureGradient = Vector2.zero;
        for (var i = 0; i < m_ParticleCount; i++)
        {
            var dif = m_Position[i] - pos;
            var dir = dif.normalized;
            var sqrDst = Vector2.SqrMagnitude(dif);
            var distance = Mathf.Sqrt(sqrDst);
            var influence = SmoothingKernelDerivative(distance, sqrDst);
            pressureGradient += dir * (m_Pressure[i] * Mass) / m_Density[i] * influence;
        }

        return pressureGradient;
    }

    float SmoothingKernelDerivative(float dst, float sqrdDistance)
    {
        if (sqrdDistance > m_SquaredSmoothingLength) return 0;
        var diff = m_SquaredSmoothingLength - sqrdDistance;
        
        return diff*diff *dst * m_KernelDerivativeTerm;
    }

    float SmoothingKernel(float sqrdDistance, float smoothingFactor)
    {
        if (sqrdDistance > m_SquaredSmoothingLength) return 0;
        //var term = 4 / (Mathf.PI * Mathf.Pow(smoothingFactor, 8));
        var diff = m_SquaredSmoothingLength - sqrdDistance;
        //return 315f / (64f * Mathf.PI * Mathf.Pow(smoothingFactor,9)) * Mathf.Pow(smoothingFactor*smoothingFactor - distance*distance, 3);
        return diff*diff*diff * m_KernelTerm;
    }
}
