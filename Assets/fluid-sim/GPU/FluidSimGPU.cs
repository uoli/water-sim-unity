using System;
using Unity.Collections;
using UnityEngine;
using UnityEngine.Serialization;
using Random = UnityEngine.Random;

public class FluidSimGPU : MonoBehaviour, IFluidSim
{
    public ComputeShader SimComputeShader;
    public int ParticleCount = 100;
    public float Width = 100;
    public float Height = 100;
    public float Mass = 1;
    public float SmoothingRadius = 1.4f;
    public float TargetDensity = 0.01f;
    public float PressureMultiplier = 50;
    public float DeltaTime = 0.001f;
    public float BoundaryPushStrength = 0;
    public float CollisionDamping = 0.5f;

    NativeArray<Vector2> m_Position;
    NativeArray<Vector2> m_PredictedPosition;
    NativeArray<float> m_Density;
    NativeArray<float> m_Pressure;
    NativeArray<Vector2> m_Velocity;
    internal GridSpatialLookup m_LookupHelper;
    
    float[] m_PointPositionData;
    float[] m_PointDensitiesData;
    float[] m_PointPressureData;
    float[] m_PointVelocityData;
    ComputeBuffer m_PointBuffer;
    ComputeBuffer m_PointDensitiesBuffer;
    ComputeBuffer m_PointPressureBuffer;
    ComputeBuffer m_PointVelocityBuffer;
    
    float[] m_DebugData;
    ComputeBuffer m_DebugBuffer;

    

    int IFluidSim.ParticleCount => ParticleCount;
    float IFluidSim.Mass => Mass;
    float IFluidSim.Height => Height;
    float IFluidSim.Width => Width;
    float IFluidSim.SmoothingRadius => SmoothingRadius;
    float IFluidSim.TargetDensity => TargetDensity;
    public GridSpatialLookup LookupHelper => throw new NotImplementedException();
    public ReadOnlySpan<Vector2> GetPositions() { return m_Position; }
    public ReadOnlySpan<float> GetDensities() { return m_Density; }
    public ReadOnlySpan<float> GetPressures() { return m_Pressure; }
    public ReadOnlySpan<Vector2> GetVelocities() { return m_Velocity; }

    void OnEnable()
    {
        InitParticleData();
    }

    void OnDisable()
    {
        CleanupParticles();
        CleanupComputeBuffers();
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

    void InitParticleData()
    {
        CleanupParticles();
        CleanupComputeBuffers();
        m_Position = new NativeArray<Vector2>(ParticleCount, Allocator.Persistent);
        m_PredictedPosition = new NativeArray<Vector2>(ParticleCount, Allocator.Persistent);
        m_Density = new NativeArray<float>(ParticleCount, Allocator.Persistent);
        m_Pressure = new NativeArray<float>(ParticleCount, Allocator.Persistent);
        m_Velocity = new NativeArray<Vector2>(ParticleCount, Allocator.Persistent);
        m_PointPositionData = new float[ParticleCount * 2];
        m_PointDensitiesData = new float[ParticleCount];
        m_PointPressureData = new float[ParticleCount];
        m_DebugData = new float[ParticleCount];
        m_PointVelocityData = new float[ParticleCount * 2];
        for ( var i = 0; i < ParticleCount; i++ )
        {
            var position = new Vector2(Random.Range(0f, Width), Random.Range(0f, Height));
            m_Position[i] = position;
            m_PredictedPosition[i] = position;
            m_Density[i] = 0;
            m_Pressure[i] = 0;
            m_Velocity[i] = Vector2.zero;
      
            m_PointPositionData[i*2] = position.x;
            m_PointPositionData[i*2 + 1] = position.y;
            m_PointDensitiesData[i] = 0;
            m_PointPressureData[i] = 0;
            m_PointVelocityData[i*2] = 0;
            m_PointVelocityData[i*2 + 1] = 0;

            m_DebugData[i] = 0;
        }
        m_PointBuffer = new ComputeBuffer(m_PointPositionData.Length, sizeof(float));
        m_PointDensitiesBuffer = new ComputeBuffer(m_PointDensitiesData.Length, sizeof(float));
        m_PointPressureBuffer = new ComputeBuffer(m_PointPressureData.Length, sizeof(float));
        m_PointVelocityBuffer = new ComputeBuffer(m_PointVelocityData.Length, sizeof(float));
        m_PointBuffer.SetData(m_PointPositionData);
        m_PointDensitiesBuffer.SetData(m_PointDensitiesData);
        m_PointPressureBuffer.SetData(m_PointDensitiesData);
        m_PointVelocityBuffer.SetData(m_PointVelocityData);
        
        m_DebugBuffer = new ComputeBuffer(m_DebugData.Length, sizeof(float));

    }
    void CleanupComputeBuffers()
    {
        if (m_PointBuffer != null)
            m_PointBuffer.Release();
        if (m_PointDensitiesBuffer != null)
            m_PointDensitiesBuffer.Release();
        if (m_PointPressureBuffer != null)
            m_PointPressureBuffer.Release();
        if (m_PointVelocityBuffer != null)
            m_PointVelocityBuffer.Release();
        if (m_DebugBuffer != null)
            m_DebugBuffer.Release();
    }


    // Update is called once per frame
    void Update()
    {
        if (ParticleCount != m_PointDensitiesBuffer.count)
        {
            InitParticleData();
        }
        
        SimComputeShader.SetFloat("ParticleCount", ParticleCount);
        SimComputeShader.SetFloat("Mass", Mass);
        SimComputeShader.SetFloat("Width", Width);
        SimComputeShader.SetFloat("Height", Height);
        SimComputeShader.SetFloat("TargetDensity", TargetDensity);
        SimComputeShader.SetFloat("PressureMultiplier", PressureMultiplier);
        SimComputeShader.SetFloat("SmoothingRadius", SmoothingRadius);
        SimComputeShader.SetFloat("SquaredSmoothingRadius", SmoothingRadius * SmoothingRadius);
        SimComputeShader.SetFloat("DeltaTime", DeltaTime);
        SimComputeShader.SetFloat("BoundaryPushStrength", BoundaryPushStrength);
        SimComputeShader.SetFloat("CollisionDamping", CollisionDamping);

        //Calculate Pressure and Density
        var kernelIndex0 = SimComputeShader.FindKernel("ComputeDensityAndPressure");
        SimComputeShader.GetKernelThreadGroupSizes(kernelIndex0, out var xGroupSize,out var yGroupSize,out var zGroupSize);
        SimComputeShader.SetBuffer(kernelIndex0, "Positions", m_PointBuffer);
        SimComputeShader.SetBuffer(kernelIndex0, "Densities", m_PointDensitiesBuffer);
        SimComputeShader.SetBuffer(kernelIndex0, "Pressure", m_PointPressureBuffer);
        SimComputeShader.SetBuffer(kernelIndex0, "Velocity", m_PointVelocityBuffer);
        SimComputeShader.SetBuffer(kernelIndex0, "DebugBuff", m_DebugBuffer);
        var threadGroupsX = Mathf.CeilToInt((float)ParticleCount / (int)xGroupSize);
        SimComputeShader.Dispatch(kernelIndex0, threadGroupsX,1,1);
        
        var kernelIndex1 = SimComputeShader.FindKernel("ComputePressureForce");
        SimComputeShader.GetKernelThreadGroupSizes(kernelIndex1, out xGroupSize,out yGroupSize,out zGroupSize);
        SimComputeShader.SetBuffer(kernelIndex1, "Positions", m_PointBuffer);
        SimComputeShader.SetBuffer(kernelIndex1, "Densities", m_PointDensitiesBuffer);
        SimComputeShader.SetBuffer(kernelIndex1, "Pressure", m_PointPressureBuffer);
        SimComputeShader.SetBuffer(kernelIndex1, "Velocity", m_PointVelocityBuffer);
        SimComputeShader.SetBuffer(kernelIndex1, "DebugBuff", m_DebugBuffer);
        threadGroupsX = Mathf.CeilToInt((float)ParticleCount / (int)xGroupSize);
        SimComputeShader.Dispatch(kernelIndex1, threadGroupsX,1,1);

        var kernelIndex2 = SimComputeShader.FindKernel("ComputePositionFromVelocityAndHandleCollision");
        SimComputeShader.GetKernelThreadGroupSizes(kernelIndex2, out xGroupSize,out yGroupSize,out zGroupSize);
        SimComputeShader.SetBuffer(kernelIndex2, "Positions", m_PointBuffer);
        SimComputeShader.SetBuffer(kernelIndex2, "Densities", m_PointDensitiesBuffer);
        SimComputeShader.SetBuffer(kernelIndex2, "Pressure", m_PointPressureBuffer);
        SimComputeShader.SetBuffer(kernelIndex2, "Velocity", m_PointVelocityBuffer);
        SimComputeShader.SetBuffer(kernelIndex2, "DebugBuff", m_DebugBuffer);
        threadGroupsX = Mathf.CeilToInt((float)ParticleCount / (int)xGroupSize);
        SimComputeShader.Dispatch(kernelIndex2, threadGroupsX,1,1);

        
        m_PointBuffer.GetData(m_PointPositionData);
        m_PointDensitiesBuffer.GetData(m_PointDensitiesData);
        m_PointPressureBuffer.GetData(m_PointPressureData);
        m_PointVelocityBuffer.GetData(m_PointVelocityData);
        m_DebugBuffer.GetData(m_DebugData);

        for (var i = 0; i < m_Position.Length; i++)
        {
            m_Position[i] = new Vector2(m_PointPositionData[i*2], m_PointPositionData[i*2 + 1]);
            m_Density[i] = m_PointDensitiesData[i];
            m_Pressure[i] = m_PointPressureData[i];
            m_Velocity[i] = new Vector2(m_PointVelocityData[i*2], m_PointVelocityData[i*2 + 1]);
        }
    }
    
    
    public void Interact(Vector2 mouseInSimulationSpace, float scalingFactor, InteractionDirection interactionDirection)
    {
        //throw new NotImplementedException();
    }
}
