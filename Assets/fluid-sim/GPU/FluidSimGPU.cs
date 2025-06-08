using System;
using System.Runtime.InteropServices;
using Unity.Collections;
using UnityEditor;
using UnityEngine;
using UnityEngine.Serialization;
using Random = UnityEngine.Random;

public class FluidSimGPU : MonoBehaviour, IFluidSim
{
    public ComputeShader SimComputeShader;
    public ComputeShader SortComputeShader;
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
    
    float[] m_PointPositionData;
    float[] m_PointDensitiesData;
    float[] m_PointPressureData;
    float[] m_PointVelocityData;
    ComputeBuffer m_PointBuffer;
    ComputeBuffer m_PredictedPositionBuffer;
    ComputeBuffer m_PointDensitiesBuffer;
    ComputeBuffer m_PointPressureBuffer;
    ComputeBuffer m_PointVelocityBuffer;
    
    [Serializable]
    [StructLayout(LayoutKind.Sequential, Pack = 1)] 
    struct SpatialEntry
    {
        internal uint index;
        internal uint cellKey;
    } 
    SpatialEntry[] m_SpatialEntry;
    int[] m_StartIndices;
    ComputeBuffer m_SpatialEntryBuffer;
    ComputeBuffer m_StartIndicesBuffer;

    
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
        m_PointVelocityData = new float[ParticleCount * 2];
        m_SpatialEntry = new SpatialEntry[ParticleCount];
        m_StartIndices= new int[ParticleCount];
        m_DebugData = new float[ParticleCount];

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
        m_PredictedPositionBuffer = new ComputeBuffer(ParticleCount, sizeof(float)*2);
        m_PointDensitiesBuffer = new ComputeBuffer(m_PointDensitiesData.Length, sizeof(float));
        m_PointPressureBuffer = new ComputeBuffer(m_PointPressureData.Length, sizeof(float));
        m_PointVelocityBuffer = new ComputeBuffer(m_PointVelocityData.Length, sizeof(float));
        m_PointBuffer.SetData(m_PointPositionData);
        m_PointDensitiesBuffer.SetData(m_PointDensitiesData);
        m_PointPressureBuffer.SetData(m_PointDensitiesData);
        m_PointVelocityBuffer.SetData(m_PointVelocityData);
        
        m_DebugBuffer = new ComputeBuffer(m_DebugData.Length, sizeof(float));
        m_DebugBuffer.SetData(m_DebugData);
        
        var size = Marshal.SizeOf(typeof(SpatialEntry));
        m_SpatialEntryBuffer = new ComputeBuffer(ParticleCount, size);
        m_SpatialEntryBuffer.SetData(m_SpatialEntry);
        
        m_StartIndicesBuffer = new ComputeBuffer(ParticleCount, sizeof(int));
        m_StartIndicesBuffer.SetData(m_StartIndices);


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
        
        //Predicted Positions
        var kernelIndex0 = SimComputeShader.FindKernel("ComputePredictedPositions");
        SimComputeShader.GetKernelThreadGroupSizes(kernelIndex0, out var xGroupSize,out var yGroupSize,out var zGroupSize);
        SimComputeShader.SetBuffer(kernelIndex0, "Positions", m_PointBuffer);
        SimComputeShader.SetBuffer(kernelIndex0, "PredictedPositions", m_PredictedPositionBuffer);
        SimComputeShader.SetBuffer(kernelIndex0, "Velocity", m_PointVelocityBuffer);
        SimComputeShader.SetBuffer(kernelIndex0, "DebugBuff", m_DebugBuffer);
        var threadGroupsX = Mathf.CeilToInt((float)ParticleCount / (int)xGroupSize);
        SimComputeShader.Dispatch(kernelIndex0, threadGroupsX,1,1);
        
        //Calculate Spatial Acceleration tables
        var kernelSpatial0 = SimComputeShader.FindKernel("UpdateSpatialLookup_HashEntries");
        SimComputeShader.GetKernelThreadGroupSizes(kernelSpatial0, out xGroupSize,out yGroupSize,out zGroupSize);
        SimComputeShader.SetBuffer(kernelSpatial0, "Positions", m_PointBuffer);
        SimComputeShader.SetBuffer(kernelSpatial0, "SpatialEntry", m_SpatialEntryBuffer);
        SimComputeShader.SetBuffer(kernelSpatial0, "StartIndices", m_StartIndicesBuffer);
        threadGroupsX = Mathf.CeilToInt((float)ParticleCount / (int)xGroupSize);
        SimComputeShader.Dispatch(kernelSpatial0, threadGroupsX,1,1);
        
        //Sort Spatial Acceleration Tables
        ComputeSort(ref m_SpatialEntryBuffer);
        
        //Calculate Spatial Acceleration start index
        var kernelSpatial2 = SimComputeShader.FindKernel("UpdateSpatialLookup_UpdateStartIndices");
        SimComputeShader.GetKernelThreadGroupSizes(kernelSpatial2, out xGroupSize,out yGroupSize,out zGroupSize);
        SimComputeShader.SetBuffer(kernelSpatial2, "SpatialEntry", m_SpatialEntryBuffer);
        SimComputeShader.SetBuffer(kernelSpatial2, "StartIndices", m_StartIndicesBuffer);
        SimComputeShader.SetBuffer(kernelSpatial2, "DebugBuff", m_DebugBuffer);
        threadGroupsX = Mathf.CeilToInt((float)ParticleCount / (int)xGroupSize);
        SimComputeShader.Dispatch(kernelSpatial2, threadGroupsX,1,1);
        
        //Calculate Pressure, Density and Predicted Positions
        var kernelIndex1 = SimComputeShader.FindKernel("ComputeDensityAndPressure");
        SimComputeShader.GetKernelThreadGroupSizes(kernelIndex1, out xGroupSize,out yGroupSize,out zGroupSize);
        SimComputeShader.SetBuffer(kernelIndex1, "SpatialEntry", m_SpatialEntryBuffer);
        SimComputeShader.SetBuffer(kernelIndex1, "StartIndices", m_StartIndicesBuffer);
        SimComputeShader.SetBuffer(kernelIndex1, "Positions", m_PointBuffer);
        SimComputeShader.SetBuffer(kernelIndex1, "PredictedPositions", m_PredictedPositionBuffer);
        SimComputeShader.SetBuffer(kernelIndex1, "Densities", m_PointDensitiesBuffer);
        SimComputeShader.SetBuffer(kernelIndex1, "Pressure", m_PointPressureBuffer);
        SimComputeShader.SetBuffer(kernelIndex1, "Velocity", m_PointVelocityBuffer);
        SimComputeShader.SetBuffer(kernelIndex1, "DebugBuff", m_DebugBuffer);
        threadGroupsX = Mathf.CeilToInt((float)ParticleCount / (int)xGroupSize);
        SimComputeShader.Dispatch(kernelIndex1, threadGroupsX,1,1);
        
        //Calculate Pressure Force
        var kernelIndex2 = SimComputeShader.FindKernel("ComputePressureForce");
        SimComputeShader.GetKernelThreadGroupSizes(kernelIndex2, out xGroupSize,out yGroupSize,out zGroupSize);
        SimComputeShader.SetBuffer(kernelIndex2, "SpatialEntry", m_SpatialEntryBuffer);
        SimComputeShader.SetBuffer(kernelIndex2, "StartIndices", m_StartIndicesBuffer);
        SimComputeShader.SetBuffer(kernelIndex2, "Positions", m_PointBuffer);
        SimComputeShader.SetBuffer(kernelIndex2, "Densities", m_PointDensitiesBuffer);
        SimComputeShader.SetBuffer(kernelIndex2, "Pressure", m_PointPressureBuffer);
        SimComputeShader.SetBuffer(kernelIndex2, "Velocity", m_PointVelocityBuffer);
        SimComputeShader.SetBuffer(kernelIndex2, "DebugBuff", m_DebugBuffer);
        threadGroupsX = Mathf.CeilToInt((float)ParticleCount / (int)xGroupSize);
        SimComputeShader.Dispatch(kernelIndex2, threadGroupsX,1,1);

        var kernelIndex3 = SimComputeShader.FindKernel("ComputePositionFromVelocityAndHandleCollision");
        SimComputeShader.GetKernelThreadGroupSizes(kernelIndex3, out xGroupSize,out yGroupSize,out zGroupSize);
        SimComputeShader.SetBuffer(kernelIndex3, "Positions", m_PointBuffer);
        SimComputeShader.SetBuffer(kernelIndex3, "Densities", m_PointDensitiesBuffer);
        SimComputeShader.SetBuffer(kernelIndex3, "Pressure", m_PointPressureBuffer);
        SimComputeShader.SetBuffer(kernelIndex3, "Velocity", m_PointVelocityBuffer);
        SimComputeShader.SetBuffer(kernelIndex3, "DebugBuff", m_DebugBuffer);
        threadGroupsX = Mathf.CeilToInt((float)ParticleCount / (int)xGroupSize);
        SimComputeShader.Dispatch(kernelIndex3, threadGroupsX,1,1);

        
        m_PointBuffer.GetData(m_PointPositionData);
        m_PointDensitiesBuffer.GetData(m_PointDensitiesData);
        m_PointPressureBuffer.GetData(m_PointPressureData);
        m_PointVelocityBuffer.GetData(m_PointVelocityData);
        m_DebugBuffer.GetData(m_DebugData);
        m_StartIndicesBuffer.GetData(m_StartIndices);
        m_SpatialEntryBuffer.GetData(m_SpatialEntry);

        for (var i = 0; i < m_Position.Length; i++)
        {
            m_Position[i] = new Vector2(m_PointPositionData[i*2], m_PointPositionData[i*2 + 1]);
            m_Density[i] = m_PointDensitiesData[i];
            m_Pressure[i] = m_PointPressureData[i];
            m_Velocity[i] = new Vector2(m_PointVelocityData[i*2], m_PointVelocityData[i*2 + 1]);
            if (m_DebugData[i] != 0)
            {
                EditorApplication.isPaused = true;
                Debug.Log($"Item {i} has DebugData {m_DebugData[i]}");
            }
        }
    }
    
    void ComputeSort(ref ComputeBuffer targetBuffer)
    {
        var originalCount = targetBuffer.count;
        var stride = targetBuffer.stride;
        var paddedCount = Mathf.NextPowerOfTwo(originalCount);
        
        //used for ping pong
        var bufferA = new ComputeBuffer(paddedCount, stride);
        var bufferB = new ComputeBuffer(paddedCount, stride);
        
        // Copy input to bufferA
        var temp = new SpatialEntry[paddedCount];
        targetBuffer.GetData(temp, 0, 0, originalCount);
        
        // Fill unused slots with float.MaxValue or another sentinel value
        for (var i = originalCount; i < paddedCount; i++)
        {
            temp[i] = new SpatialEntry { index = uint.MaxValue, cellKey = uint.MaxValue };
        }
        bufferA.SetData(temp);
        
        // Sort using ping-ponging between A and B
        var input = bufferA;
        var output = bufferB;
        
        var kernel = SortComputeShader.FindKernel("BitonicMerge");
        SortComputeShader.GetKernelThreadGroupSizes(kernel, out var xGroupSize, out _, out _);
        var groups = Mathf.CeilToInt((float)paddedCount / xGroupSize);
        
        for (int k = 2; k <= paddedCount; k <<= 1)
        {
            for (int j = k >> 1; j > 0; j >>= 1)
            {
                SortComputeShader.SetInt("k", k);
                SortComputeShader.SetInt("j", j);
                SortComputeShader.SetBuffer(kernel, "InputBuffer", input);
                SortComputeShader.SetBuffer(kernel, "OutputBuffer", output);
                SortComputeShader.Dispatch(kernel, groups, 1, 1);

                // Swap input/output for next pass
                (input, output) = (output, input);
            }
        }
        
        // After sorting, `input` holds the sorted result.
        // Copy the first N elements back to the target buffer.
        targetBuffer.Dispose();
        targetBuffer = new ComputeBuffer(originalCount, stride);
        var sortedData = new SpatialEntry[originalCount];
        input.GetData(sortedData, 0, 0, originalCount);
        targetBuffer.SetData(sortedData);

        // Cleanup
        bufferA.Dispose();
        bufferB.Dispose();
    }
    static void Swap(ref ComputeBuffer a, ref ComputeBuffer b)
    {
        (a, b) = (b, a);
    }

    public void Interact(Vector2 mouseInSimulationSpace, float scalingFactor, InteractionDirection interactionDirection)
    {
        //throw new NotImplementedException();
    }
}
