using System;
using System.Runtime.InteropServices;
using Unity.Collections;
using UnityEditor;
using UnityEngine;
using Random = UnityEngine.Random;

[Serializable]
[StructLayout(LayoutKind.Sequential, Pack = 1)] 
struct SpatialEntry
{
    internal uint index;
    internal uint cellKey;
} 

public class FluidSimGPU : MonoBehaviour, IFluidSim
{
    public ComputeShader SimComputeShader;
    public ComputeShader SortComputeShader;
    public int ParticleCount = 100;
    public int Width = 100;
    public int Height = 100;
    public float Mass = 1;
    public float SmoothingRadius = 1.4f;
    public float TargetDensity = 0.01f;
    public float PressureMultiplier = 50;
    public float DeltaTime = 0.001f;
    public float BoundaryPushStrength = 0;
    public float CollisionDamping = 0.5f;
    public float ExternalForceStrength = 10f;
    public float Gravity = 9.8f;
    
    Vector2[] m_PointPositionData;
    float[] m_PointDensitiesData;
    float[] m_PointPressureData;
    Vector2[] m_PointVelocityData;
    Vector2[] m_PredictedPositionData;
    ComputeBuffer m_PointBuffer;
    ComputeBuffer m_PredictedPositionBuffer;
    ComputeBuffer m_PointDensitiesBuffer;
    ComputeBuffer m_PointPressureBuffer;
    ComputeBuffer m_PointVelocityBuffer;
    
    SpatialEntry[] m_SpatialEntry;
    int[] m_StartIndices;
    ComputeBuffer m_SpatialEntryBuffer;
    ComputeBuffer m_StartIndicesBuffer;
    
    float[] m_DebugData;
    ComputeBuffer m_DebugBuffer;
    Vector2 m_ExternalForceCenter;
    float m_ExternalForceRadius;
    InteractionDirection m_InteractionDirection;
    IFluidSim m_FluidSimImplementation;

    int IFluidSim.ParticleCount => ParticleCount;
    float IFluidSim.Mass => Mass;
    int IFluidSim.Height => Height;
    int IFluidSim.Width => Width;
    float IFluidSim.SmoothingRadius => SmoothingRadius;
    float IFluidSim.TargetDensity => TargetDensity;
    public bool HasDataInCompute => true;
    public GridSpatialLookup LookupHelper => throw new NotImplementedException();
    public ComputeBuffer GetDensities() { return m_PointDensitiesBuffer; }
    public ComputeBuffer GetPressures() { return m_PointPressureBuffer; }
    public ComputeBuffer GetVelocities() { return m_PointVelocityBuffer; }
    
    public ComputeBuffer GetPositionComputeBuffer() { return m_PointBuffer; }

    void OnEnable()
    {
        InitParticleData();
    }

    void OnDisable()
    {
        CleanupComputeBuffers();
    }

    void InitParticleData()
    {
        CleanupComputeBuffers();
        m_PointPositionData = new Vector2[ParticleCount];
        m_PointDensitiesData = new float[ParticleCount];
        m_PointPressureData = new float[ParticleCount];
        m_PointVelocityData = new Vector2[ParticleCount];
        m_SpatialEntry = new SpatialEntry[ParticleCount];
        m_StartIndices= new int[ParticleCount];
        m_DebugData = new float[ParticleCount];
        m_PredictedPositionData = new Vector2[ParticleCount];

        for ( var i = 0; i < ParticleCount; i++ )
        {
            var position = new Vector2(Random.Range(0f, Width), Random.Range(0f, Height));
      
            m_PointPositionData[i] = position;
            m_PointDensitiesData[i] = 0;
            m_PointPressureData[i] = 0;
            m_PointVelocityData[i] = Vector2.zero;

            m_DebugData[i] = 0;
        }
        m_PointBuffer = new ComputeBuffer(ParticleCount, Marshal.SizeOf(typeof(Vector2)));
        m_PredictedPositionBuffer = new ComputeBuffer(ParticleCount, Marshal.SizeOf(typeof(Vector2)));
        m_PointDensitiesBuffer = new ComputeBuffer(ParticleCount, sizeof(float));
        m_PointPressureBuffer = new ComputeBuffer(ParticleCount, sizeof(float));
        m_PointVelocityBuffer = new ComputeBuffer(ParticleCount, Marshal.SizeOf(typeof(Vector2)));
        m_PointBuffer.SetData(m_PointPositionData);
        m_PointDensitiesBuffer.SetData(m_PointDensitiesData);
        m_PointPressureBuffer.SetData(m_PointDensitiesData);
        m_PointVelocityBuffer.SetData(m_PointVelocityData);
        m_PredictedPositionBuffer.SetData(m_PredictedPositionData);
        
        var size = Marshal.SizeOf(typeof(SpatialEntry));
        m_SpatialEntryBuffer = new ComputeBuffer(ParticleCount, size);
        m_StartIndicesBuffer = new ComputeBuffer(ParticleCount, sizeof(int));
        m_SpatialEntryBuffer.SetData(m_SpatialEntry);
        m_StartIndicesBuffer.SetData(m_StartIndices);
        
        m_DebugBuffer = new ComputeBuffer(m_DebugData.Length, sizeof(float));
        m_DebugBuffer.SetData(m_DebugData);
    }
    void CleanupComputeBuffers()
    {
        if (m_PointBuffer != null)
            m_PointBuffer.Release();
        if (m_PredictedPositionBuffer != null)
            m_PredictedPositionBuffer.Release();
        if (m_PointDensitiesBuffer != null)
            m_PointDensitiesBuffer.Release();
        if (m_PointPressureBuffer != null)
            m_PointPressureBuffer.Release();
        if (m_PointVelocityBuffer != null)
            m_PointVelocityBuffer.Release();

        if (m_SpatialEntryBuffer != null)
            m_SpatialEntryBuffer.Release();
        if (m_StartIndicesBuffer != null)
            m_StartIndicesBuffer.Release();
        
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

        var externalForceStrength = ExternalForceStrength;
        if (m_InteractionDirection == InteractionDirection.Repel)
        {
            externalForceStrength *= -1;
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
        SimComputeShader.SetFloat("ForceRadius", m_ExternalForceRadius); 
        SimComputeShader.SetFloat("ForceCenterX", m_ExternalForceCenter.x);
        SimComputeShader.SetFloat("ForceCenterY", m_ExternalForceCenter.y);
        SimComputeShader.SetFloat("ForceStrength", externalForceStrength); 
        SimComputeShader.SetFloat("Gravity", Gravity); 
        
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
        SimComputeShader.SetBuffer(kernelSpatial0, "PredictedPositions", m_PredictedPositionBuffer);
        SimComputeShader.SetBuffer(kernelSpatial0, "SpatialEntry", m_SpatialEntryBuffer);
        SimComputeShader.SetBuffer(kernelSpatial0, "StartIndices", m_StartIndicesBuffer);
        threadGroupsX = Mathf.CeilToInt((float)ParticleCount / (int)xGroupSize);
        SimComputeShader.Dispatch(kernelSpatial0, threadGroupsX,1,1);
        
        //Sort Spatial Acceleration Tables
        //TODO: This is the bottleneck right now
        //It needs to copy things into and out of the GPU
        //Also it dispatches like 50 jobs
        ComputeSort(ref m_SpatialEntryBuffer);
        
        //Calculate Spatial Acceleration start index
        var kernelSpatial2 = SimComputeShader.FindKernel("UpdateSpatialLookup_UpdateStartIndices");
        SimComputeShader.GetKernelThreadGroupSizes(kernelSpatial2, out xGroupSize,out yGroupSize,out zGroupSize);
        SimComputeShader.SetBuffer(kernelSpatial2, "SpatialEntry", m_SpatialEntryBuffer);
        SimComputeShader.SetBuffer(kernelSpatial2, "StartIndices", m_StartIndicesBuffer);
        SimComputeShader.SetBuffer(kernelSpatial2, "DebugBuff", m_DebugBuffer);
        threadGroupsX = Mathf.CeilToInt((float)ParticleCount / (int)xGroupSize);
        SimComputeShader.Dispatch(kernelSpatial2, threadGroupsX,1,1);
        
        //Calculate Pressure and Density 
        var kernelIndex1 = SimComputeShader.FindKernel("ComputeDensityAndPressure");
        SimComputeShader.GetKernelThreadGroupSizes(kernelIndex1, out xGroupSize,out yGroupSize,out zGroupSize);
        SimComputeShader.SetBuffer(kernelIndex1, "SpatialEntry", m_SpatialEntryBuffer);
        SimComputeShader.SetBuffer(kernelIndex1, "StartIndices", m_StartIndicesBuffer);
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
        SimComputeShader.SetBuffer(kernelIndex2, "PredictedPositions", m_PredictedPositionBuffer);
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
        
        /*
        m_PointBuffer.GetData(m_PointPositionData);
        m_PointDensitiesBuffer.GetData(m_PointDensitiesData);
        m_PointPressureBuffer.GetData(m_PointPressureData);
        m_PointVelocityBuffer.GetData(m_PointVelocityData);
        m_DebugBuffer.GetData(m_DebugData);
        m_StartIndicesBuffer.GetData(m_StartIndices);
        m_SpatialEntryBuffer.GetData(m_SpatialEntry);
        m_PredictedPositionBuffer.GetData(m_PredictedPositionData);
*/
        m_ExternalForceRadius = 0; //clear External Force interaction

        
        for (var i = 0; i < ParticleCount; i++)
        {
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

    public void Interact(Vector2 mouseInSimulationSpace, float radius, InteractionDirection interactionDirection)
    {
        m_ExternalForceCenter = mouseInSimulationSpace;
        m_ExternalForceRadius = radius;
        m_InteractionDirection = interactionDirection;
    }
}
