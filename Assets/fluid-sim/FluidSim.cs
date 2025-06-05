using System.Collections.Generic;
using Unity.Profiling;
using UnityEngine;

public struct Particle //probably an Struct of arrays could be beneficial
{
    public Vector2 position;
    public float density;
    public float pressure;
}

public class FluidSim : MonoBehaviour
{
    public int m_ParticleCount;

    public int width;
    public int height;
    public float Mass = 1;
    public float SmoothingLength = 10;
    
    Particle[] m_Particles;

    float m_SquaredSmoothingLength;
    float m_KernelTerm;
    float m_TargetDensity;
    
    static readonly ProfilerMarker s_UpdatePerfMarker = new ProfilerMarker("FluidSim.Update");
    public IReadOnlyList<Particle> GetParticles() { return m_Particles; }
    
    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        InitParticles();
    }

    void InitParticles()
    {
        m_Particles = new Particle[m_ParticleCount];
        for ( var i = 0; i < m_ParticleCount; i++ )
        {
            var position = new Vector2(Random.Range(0f, width), Random.Range(0f, height));
            m_Particles[i] = new Particle() { position = position };
        }
    }

    // Update is called once per frame
    void Update()
    {
        using var markerScope = s_UpdatePerfMarker.Auto();

        if (m_ParticleCount != m_Particles.Length)
        {
            InitParticles();
        }

        CachePrecomputedValues();
        CalculateParticlesDensity();
        CalculateParticlePressure();
    }

    public float CalcTargetDensity()
    {
        return m_ParticleCount / (float)(width * height);
    }
    
    void CachePrecomputedValues()
    {
        m_SquaredSmoothingLength = SmoothingLength*SmoothingLength;
        m_KernelTerm = 4 / (Mathf.PI * Mathf.Pow(SmoothingLength, 8));
        m_TargetDensity = CalcTargetDensity();
    }
    
    float CalculateDensity(Vector2 pos)
    {
        var density = 0f;
        for (var i = 0; i < m_ParticleCount; i++)
        {
            var particle = m_Particles[i];
            
            var sqrDst = Vector2.SqrMagnitude(particle.position - pos);
            var influence = SmoothingKernel(sqrDst, SmoothingLength);
            density += Mass * influence;
        }

        return density;
    }

    void CalculateParticlesDensity()
    {
        for (var index = 0; index < m_Particles.Length; index++)
        {
            var particle = m_Particles[index];
            m_Particles[index].density = CalculateDensity(particle.position);
        }
    }

    void CalculateParticlePressure()
    {
        const float stiffness = 3000;
        const float adiabaticComponent = 7;
        for (var index = 0; index < m_Particles.Length; index++)
        {
            var particle = m_Particles[index];
            m_Particles[index].pressure = stiffness * (Mathf.Pow(particle.density / m_TargetDensity, adiabaticComponent) - 1);
        }
    }

    public float CalculatePressure(Vector2 pos)
    {
        var pressure = 0f;
        for (var i = 0; i < m_ParticleCount; i++)
        {
            var particle = m_Particles[i];
            
            var sqrDst = Vector2.SqrMagnitude(particle.position - pos);
            var influence = SmoothingKernel(sqrDst, SmoothingLength);
            var density = CalculateDensity(pos);
            pressure += particle.pressure * Mass / density * influence;
        }

        return pressure;
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
