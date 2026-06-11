using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using Unity.Profiling;
using UnityEngine;
using Random = UnityEngine.Random;

[Serializable]
[StructLayout(LayoutKind.Sequential, Pack = 1)]
struct SpatialEntry
{
    internal uint index;
    internal uint cellKey;
    // Must match the int2 cell field of the HLSL SpatialEntryData/SpatialEntry structs
    internal int cellX;
    internal int cellY;
}

public class FluidSimGPU : MonoBehaviour, IFluidSim
{
    public ComputeShader SimComputeShader;
    public int ParticleCount = 100;
    public int Width = 100;
    public int Height = 100;
    // The fluid's rest density. Particle mass and target density derive from
    // this and the particle spacing, so the knobs stay independent of
    // resolution and domain size.
    public float RestDensity = 1f;
    // Kernel support radius in units of particle spacing. ~2 gives ~12
    // neighbors in 2D; raise for smoother fields at higher cost.
    [Range(1.2f, 4f)]
    public float SmoothingRadiusInSpacings = 2f;
    public float PressureMultiplier = 50;
    public int MaxStepsPerFrame = 8;
    // Derive the substep from the CFL condition instead of using DeltaTime.
    public bool AutoDeltaTime = true;
    public float DeltaTime = 0.001f;
    public float BoundaryPushStrength = 0;
    public float CollisionDamping = 0.5f;
    public float ExternalForceStrength = 10f;
    public float Gravity = 9.8f;
    public float ViscosityFactor = 0.5f;
    public ParticlePlacementMode PlacementMode = ParticlePlacementMode.GridWithJitter;
    // Fraction of the box (from the bottom) the fluid occupies at rest. Below 1
    // there is a free surface, so the fluid can slosh and splash.
    [Range(0.1f, 1f)]
    public float FillFraction = 0.6f;
    
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
    
    ComputeBuffer m_SpatialEntryBuffer;
    ComputeBuffer m_UnsortedSpatialEntryBuffer;
    ComputeBuffer m_CellCountsBuffer;
    ComputeBuffer m_CellStartsBuffer;
    ComputeBuffer m_ScatterOffsetsBuffer;
    ComputeBuffer m_BlockSumsBuffer;
    ComputeBuffer m_VelocityDeltaBuffer;

    const int k_ScanBlockSize = 512;

    int m_KernelPredict;
    int m_KernelSpatialClear;
    int m_KernelSpatialHashAndCount;
    int m_KernelScanBlocks;
    int m_KernelScanBlockSums;
    int m_KernelAddBlockOffsets;
    int m_KernelScatter;
    int m_KernelDensityPressure;
    int m_KernelForces;
    int m_KernelIntegrate;
    int m_GroupsParticles;
    int m_GroupsTable;
    int m_ScanBlockCount;
    int m_CellTableSize;
    
    float[] m_DebugData;
    ComputeBuffer m_DebugBuffer;
    Vector2 m_ExternalForceCenter;
    float m_ExternalForceRadius;
    InteractionDirection m_InteractionDirection;
    IFluidSim m_FluidSimImplementation;
    float m_TimeAccumulator;

    static readonly ProfilerMarker s_UpdatePerfMarker = new ProfilerMarker("FluidSimGPU.Update");
    static readonly ProfilerMarker s_SortPerfMarker = new ProfilerMarker("FluidSimGPU.Sort");
    static readonly ProfilerMarker s_SimKernelsPerfMarker = new ProfilerMarker("FluidSimGPU.SimKernels");
    ParticlePlacementMode m_LastPlacementMode;
    float m_LastFillFraction;

    int IFluidSim.ParticleCount => ParticleCount;
    int IFluidSim.Height => Height;
    int IFluidSim.Width => Width;

    // Derived quantities: the spacing comes from how many particles fill the
    // spawn region, and mass follows so that the fill sits exactly at RestDensity.
    public float ParticleSpacing => Mathf.Sqrt(Width * Height * FillFraction / Mathf.Max(1, ParticleCount));
    public float Mass => RestDensity * ParticleSpacing * ParticleSpacing;
    public float SmoothingRadius => SmoothingRadiusInSpacings * ParticleSpacing;
    public float TargetDensity => RestDensity;
    public float EffectiveDeltaTime => AutoDeltaTime ? CalcCFLDeltaTime() : DeltaTime;

    // CFL bound: a substep must not let information travel more than a fraction
    // of the kernel radius. The velocity scale is the speed of sound of the
    // linear EOS (sqrt of the pressure multiplier) plus the worst-case
    // free-fall speed over the box height; the GPU path never reads particle
    // data back, so the actual max velocity is not available.
    float CalcCFLDeltaTime()
    {
        const float courantNumber = 0.3f;
        var speedOfSound = Mathf.Sqrt(Mathf.Max(1e-6f, PressureMultiplier));
        var freeFallSpeed = Mathf.Sqrt(2f * Mathf.Abs(Gravity) * Height);
        return courantNumber * SmoothingRadius / (speedOfSound + freeFallSpeed);
    }
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
        m_DebugData = new float[ParticleCount];
        m_PredictedPositionData = new Vector2[ParticleCount];

        m_LastPlacementMode = PlacementMode;
        m_LastFillFraction = FillFraction;
        for ( var i = 0; i < ParticleCount; i++ )
        {
            var position = ParticlePlacement.GetPosition(PlacementMode, i, ParticleCount, Width, Height * FillFraction);

            m_PointPositionData[i] = position;
            m_PointDensitiesData[i] = 0;
            m_PointPressureData[i] = 0;
            m_PointVelocityData[i] = Vector2.zero;

            m_DebugData[i] = 0;
        }
        m_PointBuffer = new ComputeBuffer(ParticleCount, Marshal.SizeOf(typeof(Vector2)));
        m_PointBuffer.name = "PositionsBuffer";
        m_PredictedPositionBuffer = new ComputeBuffer(ParticleCount, Marshal.SizeOf(typeof(Vector2)));
        m_PredictedPositionBuffer.name = "PredictedPositionsBuffer";
        m_PointDensitiesBuffer = new ComputeBuffer(ParticleCount, sizeof(float));
        m_PointDensitiesBuffer.name = "DensitiesBuffer";
        m_PointPressureBuffer = new ComputeBuffer(ParticleCount, sizeof(float));
        m_PointDensitiesBuffer.name = "PressureBuffer";
        m_PointVelocityBuffer = new ComputeBuffer(ParticleCount, Marshal.SizeOf(typeof(Vector2)));
        m_PointVelocityBuffer.name = "VelocityBuffer";
        m_VelocityDeltaBuffer = new ComputeBuffer(ParticleCount, Marshal.SizeOf(typeof(Vector2)));
        m_VelocityDeltaBuffer.name = "VelocityDeltaBuffer";
        m_PointBuffer.SetData(m_PointPositionData);
        m_PointDensitiesBuffer.SetData(m_PointDensitiesData);
        m_PointPressureBuffer.SetData(m_PointDensitiesData);
        m_PointVelocityBuffer.SetData(m_PointVelocityData);
        m_PredictedPositionBuffer.SetData(m_PredictedPositionData);
        
        var size = Marshal.SizeOf(typeof(SpatialEntry));
        m_SpatialEntryBuffer = new ComputeBuffer(ParticleCount, size);
        m_SpatialEntryBuffer.name = "SpatialEntriesBuffer";
        m_UnsortedSpatialEntryBuffer = new ComputeBuffer(ParticleCount, size);
        m_UnsortedSpatialEntryBuffer.name = "UnsortedSpatialEntriesBuffer";

        // Counting-sort cell table. Sized at ~2x the particle count (rounded to
        // a power of two so the hash modulo is a bit-mask): a load factor of
        // ~0.5 keeps hash buckets mostly unshared, so neighbor queries waste
        // less time rejecting entries from far-away colliding cells.
        m_CellTableSize = Mathf.Min(k_ScanBlockSize * k_ScanBlockSize, Mathf.NextPowerOfTwo(ParticleCount * 2));
        m_CellCountsBuffer = new ComputeBuffer(m_CellTableSize, sizeof(uint));
        m_CellCountsBuffer.name = "CellCountsBuffer";
        m_CellStartsBuffer = new ComputeBuffer(m_CellTableSize, sizeof(uint));
        m_CellStartsBuffer.name = "CellStartsBuffer";
        m_ScatterOffsetsBuffer = new ComputeBuffer(m_CellTableSize, sizeof(uint));
        m_ScatterOffsetsBuffer.name = "ScatterOffsetsBuffer";
        m_BlockSumsBuffer = new ComputeBuffer(k_ScanBlockSize, sizeof(uint));
        m_BlockSumsBuffer.name = "ScanBlockSumsBuffer";
        
        m_DebugBuffer = new ComputeBuffer(m_DebugData.Length, sizeof(float));
        m_DebugBuffer.name = "DebugBuffer";
        m_DebugBuffer.SetData(m_DebugData);

        SetupKernels();
    }

    // Kernel ids and buffer bindings persist on the compute shader, so they are
    // set up once here instead of per dispatch; the substep loop only issues
    // Dispatch calls.
    void SetupKernels()
    {
        m_KernelPredict = SimComputeShader.FindKernel("ComputePredictedPositions");
        m_KernelSpatialClear = SimComputeShader.FindKernel("SpatialClear");
        m_KernelSpatialHashAndCount = SimComputeShader.FindKernel("SpatialHashAndCount");
        m_KernelScanBlocks = SimComputeShader.FindKernel("SpatialScanBlocks");
        m_KernelScanBlockSums = SimComputeShader.FindKernel("SpatialScanBlockSums");
        m_KernelAddBlockOffsets = SimComputeShader.FindKernel("SpatialAddBlockOffsets");
        m_KernelScatter = SimComputeShader.FindKernel("SpatialScatter");
        m_KernelDensityPressure = SimComputeShader.FindKernel("ComputeDensityAndPressure");
        m_KernelForces = SimComputeShader.FindKernel("ComputeForces");
        m_KernelIntegrate = SimComputeShader.FindKernel("ComputePositionFromVelocityAndHandleCollision");

        m_GroupsParticles = Mathf.CeilToInt(ParticleCount / 64f);
        m_GroupsTable = Mathf.CeilToInt(m_CellTableSize / 64f);
        m_ScanBlockCount = Mathf.CeilToInt(m_CellTableSize / (float)k_ScanBlockSize);
        if (ParticleCount * 2 > k_ScanBlockSize * k_ScanBlockSize)
            Debug.LogError($"ParticleCount {ParticleCount} exceeds the cell table limit of {k_ScanBlockSize * k_ScanBlockSize / 2}; expect degraded neighbor lookups.");

        SimComputeShader.SetInt("NumBlocks", m_ScanBlockCount);
        SimComputeShader.SetInt("CellTableSize", m_CellTableSize);

        SimComputeShader.SetBuffer(m_KernelPredict, "Positions", m_PointBuffer);
        SimComputeShader.SetBuffer(m_KernelPredict, "PredictedPositions", m_PredictedPositionBuffer);
        SimComputeShader.SetBuffer(m_KernelPredict, "Velocity", m_PointVelocityBuffer);

        SimComputeShader.SetBuffer(m_KernelSpatialClear, "CellCounts", m_CellCountsBuffer);
        SimComputeShader.SetBuffer(m_KernelSpatialClear, "ScatterOffsets", m_ScatterOffsetsBuffer);

        SimComputeShader.SetBuffer(m_KernelSpatialHashAndCount, "PredictedPositions", m_PredictedPositionBuffer);
        SimComputeShader.SetBuffer(m_KernelSpatialHashAndCount, "UnsortedSpatialEntry", m_UnsortedSpatialEntryBuffer);
        SimComputeShader.SetBuffer(m_KernelSpatialHashAndCount, "CellCounts", m_CellCountsBuffer);

        SimComputeShader.SetBuffer(m_KernelScanBlocks, "CellCounts", m_CellCountsBuffer);
        SimComputeShader.SetBuffer(m_KernelScanBlocks, "CellStarts", m_CellStartsBuffer);
        SimComputeShader.SetBuffer(m_KernelScanBlocks, "BlockSums", m_BlockSumsBuffer);

        SimComputeShader.SetBuffer(m_KernelScanBlockSums, "BlockSums", m_BlockSumsBuffer);

        SimComputeShader.SetBuffer(m_KernelAddBlockOffsets, "CellStarts", m_CellStartsBuffer);
        SimComputeShader.SetBuffer(m_KernelAddBlockOffsets, "BlockSums", m_BlockSumsBuffer);

        SimComputeShader.SetBuffer(m_KernelScatter, "UnsortedSpatialEntry", m_UnsortedSpatialEntryBuffer);
        SimComputeShader.SetBuffer(m_KernelScatter, "SpatialEntry", m_SpatialEntryBuffer);
        SimComputeShader.SetBuffer(m_KernelScatter, "CellStarts", m_CellStartsBuffer);
        SimComputeShader.SetBuffer(m_KernelScatter, "ScatterOffsets", m_ScatterOffsetsBuffer);

        SimComputeShader.SetBuffer(m_KernelDensityPressure, "SpatialEntry", m_SpatialEntryBuffer);
        SimComputeShader.SetBuffer(m_KernelDensityPressure, "CellStarts", m_CellStartsBuffer);
        SimComputeShader.SetBuffer(m_KernelDensityPressure, "CellCounts", m_CellCountsBuffer);
        SimComputeShader.SetBuffer(m_KernelDensityPressure, "PredictedPositions", m_PredictedPositionBuffer);
        SimComputeShader.SetBuffer(m_KernelDensityPressure, "Densities", m_PointDensitiesBuffer);
        SimComputeShader.SetBuffer(m_KernelDensityPressure, "Pressure", m_PointPressureBuffer);

        SimComputeShader.SetBuffer(m_KernelForces, "SpatialEntry", m_SpatialEntryBuffer);
        SimComputeShader.SetBuffer(m_KernelForces, "CellStarts", m_CellStartsBuffer);
        SimComputeShader.SetBuffer(m_KernelForces, "CellCounts", m_CellCountsBuffer);
        SimComputeShader.SetBuffer(m_KernelForces, "PredictedPositions", m_PredictedPositionBuffer);
        SimComputeShader.SetBuffer(m_KernelForces, "Densities", m_PointDensitiesBuffer);
        SimComputeShader.SetBuffer(m_KernelForces, "Pressure", m_PointPressureBuffer);
        SimComputeShader.SetBuffer(m_KernelForces, "Velocity", m_PointVelocityBuffer);
        SimComputeShader.SetBuffer(m_KernelForces, "VelocityDelta", m_VelocityDeltaBuffer);

        SimComputeShader.SetBuffer(m_KernelIntegrate, "Positions", m_PointBuffer);
        SimComputeShader.SetBuffer(m_KernelIntegrate, "Velocity", m_PointVelocityBuffer);
        SimComputeShader.SetBuffer(m_KernelIntegrate, "VelocityDelta", m_VelocityDeltaBuffer);
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
        if (m_VelocityDeltaBuffer != null)
            m_VelocityDeltaBuffer.Release();

        if (m_SpatialEntryBuffer != null)
            m_SpatialEntryBuffer.Release();
        if (m_UnsortedSpatialEntryBuffer != null)
            m_UnsortedSpatialEntryBuffer.Release();
        if (m_CellCountsBuffer != null)
            m_CellCountsBuffer.Release();
        if (m_CellStartsBuffer != null)
            m_CellStartsBuffer.Release();
        if (m_ScatterOffsetsBuffer != null)
            m_ScatterOffsetsBuffer.Release();
        if (m_BlockSumsBuffer != null)
            m_BlockSumsBuffer.Release();
        
        if (m_DebugBuffer != null)
            m_DebugBuffer.Release();
    }


    // Update is called once per frame
    void Update()
    {
        using var markerScope = s_UpdatePerfMarker.Auto();

        if (ParticleCount != m_PointDensitiesBuffer.count || PlacementMode != m_LastPlacementMode || !Mathf.Approximately(FillFraction, m_LastFillFraction))
        {
            InitParticleData();
        }

        var externalForceStrength = ExternalForceStrength;
        if (m_InteractionDirection == InteractionDirection.Repel)
        {
            externalForceStrength *= -1;
        }
        
        SimComputeShader.SetInt("ParticleCount", ParticleCount);
        SimComputeShader.SetFloat("Mass", Mass);
        SimComputeShader.SetFloat("Width", Width);
        SimComputeShader.SetFloat("Height", Height);
        SimComputeShader.SetFloat("TargetDensity", TargetDensity);
        SimComputeShader.SetFloat("PressureMultiplier", PressureMultiplier);
        SimComputeShader.SetFloat("SmoothingRadius", SmoothingRadius);
        SimComputeShader.SetFloat("SquaredSmoothingRadius", SmoothingRadius * SmoothingRadius);
        var deltaTime = EffectiveDeltaTime;
        SimComputeShader.SetFloat("DeltaTime", deltaTime);
        SimComputeShader.SetFloat("BoundaryPushStrength", BoundaryPushStrength);
        SimComputeShader.SetFloat("CollisionDamping", CollisionDamping); 
        SimComputeShader.SetFloat("ForceRadius", m_ExternalForceRadius); 
        SimComputeShader.SetFloat("ForceCenterX", m_ExternalForceCenter.x);
        SimComputeShader.SetFloat("ForceCenterY", m_ExternalForceCenter.y);
        SimComputeShader.SetFloat("ForceStrength", externalForceStrength); 
        SimComputeShader.SetFloat("Gravity", Gravity);
        SimComputeShader.SetFloat("ViscosityFactor", ViscosityFactor);

        // Fixed-timestep accumulator: take deltaTime-sized substeps to cover the
        // elapsed wall-clock time, so simulation speed is independent of framerate.
        m_TimeAccumulator += Time.deltaTime;
        var steps = 0;
        while (deltaTime > 0 && m_TimeAccumulator >= deltaTime && steps < MaxStepsPerFrame)
        {
            steps++;
            m_TimeAccumulator -= deltaTime;

            s_SimKernelsPerfMarker.Begin();
            SimComputeShader.Dispatch(m_KernelPredict, m_GroupsParticles, 1, 1);
            s_SimKernelsPerfMarker.End();

            //Build the spatial lookup with a counting sort keyed by cell
            s_SortPerfMarker.Begin();
            SimComputeShader.Dispatch(m_KernelSpatialClear, m_GroupsTable, 1, 1);
            SimComputeShader.Dispatch(m_KernelSpatialHashAndCount, m_GroupsParticles, 1, 1);
            SimComputeShader.Dispatch(m_KernelScanBlocks, m_ScanBlockCount, 1, 1);
            SimComputeShader.Dispatch(m_KernelScanBlockSums, 1, 1, 1);
            SimComputeShader.Dispatch(m_KernelAddBlockOffsets, m_GroupsTable, 1, 1);
            SimComputeShader.Dispatch(m_KernelScatter, m_GroupsParticles, 1, 1);
            s_SortPerfMarker.End();

            s_SimKernelsPerfMarker.Begin();
            SimComputeShader.Dispatch(m_KernelDensityPressure, m_GroupsParticles, 1, 1);
            SimComputeShader.Dispatch(m_KernelForces, m_GroupsParticles, 1, 1);
            SimComputeShader.Dispatch(m_KernelIntegrate, m_GroupsParticles, 1, 1);
            s_SimKernelsPerfMarker.End();
        }
        // Running behind realtime: drop the surplus so the sim slows down gracefully
        // instead of accumulating an ever-growing debt of steps.
        if (steps == MaxStepsPerFrame)
            m_TimeAccumulator = 0;
        /*
        m_PointBuffer.GetData(m_PointPositionData);
        m_PointDensitiesBuffer.GetData(m_PointDensitiesData);
        m_PointPressureBuffer.GetData(m_PointPressureData);
        m_PointVelocityBuffer.GetData(m_PointVelocityData);
        m_DebugBuffer.GetData(m_DebugData);
        m_PredictedPositionBuffer.GetData(m_PredictedPositionData);
*/
        m_ExternalForceRadius = 0; //clear External Force interaction

        
        // for (var i = 0; i < ParticleCount; i++)
        // {
        //     if (m_DebugData[i] != 0)
        //     {
        //         EditorApplication.isPaused = true;
        //         Debug.Log($"Item {i} has DebugData {m_DebugData[i]}");
        //     }
        // }
    }
    

    public void Interact(Vector2 mouseInSimulationSpace, float radius, InteractionDirection interactionDirection)
    {
        m_ExternalForceCenter = mouseInSimulationSpace;
        m_ExternalForceRadius = radius;
        m_InteractionDirection = interactionDirection;
    }

    // ---- Rigid-body coupling: not implemented on the GPU path yet. ----
    // The interface members exist so FluidRigibodyInteraction can target either
    // sim; on this path registration is accepted and ignored.
#pragma warning disable 0067 // events never invoked on this path
    public event Action PreSimulation;
    public event Action PostSimulation;
#pragma warning restore 0067
    public void SetRigidBodySurfaceResults(IList<InputSimulationSurfacePoints> points) { }
    public void RetrieveRigidBodySurfaceResults(IList<OutputSimulationSurfacePoints> points) { }
    public ComputeBuffer InputExternalPoints => null;
    Transform IFluidSim.Transform => transform;
}
