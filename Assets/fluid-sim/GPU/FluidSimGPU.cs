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
    ComputeBuffer m_SurfacePointsBuffer;
    ComputeBuffer m_SurfaceResultsBuffer;
    ComputeBuffer m_SurfaceAccumBuffer;
    InputSimulationSurfacePoints[] m_SurfacePointData;
    Vector2[] m_RetrievedImpulses;
    readonly Queue<(UnityEngine.Rendering.AsyncGPUReadbackRequest request, float simTime)> m_PendingSurfaceReadbacks = new();
    int m_SurfacePointCount;
    float m_LastCoveredSimTime;

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
    int m_KernelSurfaceForces;
    int m_KernelClearSurfaceAccum;
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
        if (m_SurfacePointsBuffer == null)
            CreateSurfaceBuffers(1); // dummy size so kernels always have valid bindings
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
        m_KernelSurfaceForces = SimComputeShader.FindKernel("ComputeSurfaceForces");
        m_KernelClearSurfaceAccum = SimComputeShader.FindKernel("ClearSurfaceAccumulation");

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

        BindSurfaceBuffers();
    }

    void CreateSurfaceBuffers(int count)
    {
        m_SurfacePointsBuffer?.Release();
        m_SurfaceResultsBuffer?.Release();
        m_SurfaceAccumBuffer?.Release();
        m_SurfacePointsBuffer = new ComputeBuffer(count, Marshal.SizeOf(typeof(InputSimulationSurfacePoints)));
        m_SurfacePointsBuffer.name = "SurfacePointsBuffer";
        m_SurfaceResultsBuffer = new ComputeBuffer(count, Marshal.SizeOf(typeof(OutputSimulationSurfacePoints)));
        m_SurfaceResultsBuffer.name = "SurfaceResultsBuffer";
        m_SurfaceAccumBuffer = new ComputeBuffer(count, Marshal.SizeOf(typeof(Vector2)));
        m_SurfaceAccumBuffer.name = "SurfaceImpulseAccumBuffer";
        m_SurfaceAccumBuffer.SetData(new Vector2[count]);
    }

    void BindSurfaceBuffers()
    {
        SimComputeShader.SetBuffer(m_KernelDensityPressure, "SurfacePoints", m_SurfacePointsBuffer);

        SimComputeShader.SetBuffer(m_KernelSurfaceForces, "SurfacePoints", m_SurfacePointsBuffer);
        SimComputeShader.SetBuffer(m_KernelSurfaceForces, "SurfaceResults", m_SurfaceResultsBuffer);
        SimComputeShader.SetBuffer(m_KernelSurfaceForces, "SurfaceImpulseAccum", m_SurfaceAccumBuffer);
        SimComputeShader.SetBuffer(m_KernelSurfaceForces, "SpatialEntry", m_SpatialEntryBuffer);
        SimComputeShader.SetBuffer(m_KernelSurfaceForces, "CellStarts", m_CellStartsBuffer);
        SimComputeShader.SetBuffer(m_KernelSurfaceForces, "CellCounts", m_CellCountsBuffer);
        SimComputeShader.SetBuffer(m_KernelSurfaceForces, "PredictedPositions", m_PredictedPositionBuffer);
        SimComputeShader.SetBuffer(m_KernelSurfaceForces, "Densities", m_PointDensitiesBuffer);
        SimComputeShader.SetBuffer(m_KernelSurfaceForces, "Pressure", m_PointPressureBuffer);
        SimComputeShader.SetBuffer(m_KernelSurfaceForces, "Velocity", m_PointVelocityBuffer);

        SimComputeShader.SetBuffer(m_KernelForces, "SurfacePoints", m_SurfacePointsBuffer);
        SimComputeShader.SetBuffer(m_KernelForces, "SurfaceResults", m_SurfaceResultsBuffer);

        SimComputeShader.SetBuffer(m_KernelIntegrate, "SurfacePoints", m_SurfacePointsBuffer);

        SimComputeShader.SetBuffer(m_KernelClearSurfaceAccum, "SurfaceImpulseAccum", m_SurfaceAccumBuffer);
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
        if (m_SurfacePointsBuffer != null)
            m_SurfacePointsBuffer.Release();
        if (m_SurfaceResultsBuffer != null)
            m_SurfaceResultsBuffer.Release();
        if (m_SurfaceAccumBuffer != null)
            m_SurfaceAccumBuffer.Release();
        m_SurfacePointsBuffer = null;
        m_SurfaceResultsBuffer = null;
        m_SurfaceAccumBuffer = null;
        m_PendingSurfaceReadbacks.Clear();

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

        // Registers this frame's surface points (body pose is fixed-step
        // driven, so once per frame is exact within the frame).
        PreSimulation?.Invoke();

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
        SimComputeShader.SetInt("SurfacePointCount", m_SurfacePointCount);
        SimComputeShader.SetFloat("ContactDistance", SmoothingRadius * 0.25f);

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
            if (m_SurfacePointCount > 0)
                SimComputeShader.Dispatch(m_KernelSurfaceForces, Mathf.CeilToInt(m_SurfacePointCount / 64f), 1, 1);
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
        // Read the frame's accumulated surface impulses back asynchronously and
        // reset the accumulator for the next frame; deliver any completed
        // readbacks to the coupling.
        if (m_SurfacePointCount > 0 && steps > 0)
        {
            m_PendingSurfaceReadbacks.Enqueue(
                (UnityEngine.Rendering.AsyncGPUReadback.Request(m_SurfaceAccumBuffer), steps * deltaTime));
            SimComputeShader.Dispatch(m_KernelClearSurfaceAccum, Mathf.CeilToInt(m_SurfacePointCount / 64f), 1, 1);
        }
        ProcessSurfaceReadbacks();

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

    // ---- Rigid-body coupling ----
    // PreSimulation fires once per frame (the body's pose is constant within a
    // frame's substeps); the registered surface points are uploaded once and
    // used by every substep. Impulses accumulate on the GPU across substeps
    // and are read back asynchronously: PostSimulation fires when a readback
    // completes, with LastStepDeltaTime reporting the sim time those impulses
    // cover. The frame of latency is absorbed by the coupling's rate
    // estimator, with no pipeline stall.
    public event Action PreSimulation;
    public event Action PostSimulation;

    public void SetRigidBodySurfaceResults(IList<InputSimulationSurfacePoints> points)
    {
        m_SurfacePointCount = points.Count;
        if (m_SurfacePointCount == 0) return;
        if (m_SurfacePointData == null || m_SurfacePointData.Length < m_SurfacePointCount)
        {
            m_SurfacePointData = new InputSimulationSurfacePoints[m_SurfacePointCount];
            m_RetrievedImpulses = new Vector2[m_SurfacePointCount];
            CreateSurfaceBuffers(m_SurfacePointCount);
            BindSurfaceBuffers();
        }
        for (var i = 0; i < m_SurfacePointCount; i++)
            m_SurfacePointData[i] = points[i];
        m_SurfacePointsBuffer.SetData(m_SurfacePointData, 0, 0, m_SurfacePointCount);
    }

    public void RetrieveRigidBodySurfaceResults(IList<OutputSimulationSurfacePoints> points)
    {
        if (m_RetrievedImpulses == null) return;
        var count = Mathf.Min(m_RetrievedImpulses.Length, points.Count);
        for (var i = 0; i < count; i++)
        {
            points[i] = new OutputSimulationSurfacePoints { impulse = m_RetrievedImpulses[i], weightSum = 1f };
        }
    }

    void ProcessSurfaceReadbacks()
    {
        while (m_PendingSurfaceReadbacks.Count > 0 && m_PendingSurfaceReadbacks.Peek().request.done)
        {
            var (request, simTime) = m_PendingSurfaceReadbacks.Dequeue();
            if (request.hasError) continue;
            var data = request.GetData<Vector2>();
            var count = Mathf.Min(data.Length, m_RetrievedImpulses?.Length ?? 0);
            for (var i = 0; i < count; i++)
                m_RetrievedImpulses[i] = data[i];
            m_LastCoveredSimTime = simTime;
            PostSimulation?.Invoke();
        }
    }

    public ComputeBuffer InputExternalPoints => m_SurfacePointCount > 0 ? m_SurfacePointsBuffer : null;
    Transform IFluidSim.Transform => transform;
    public float LastStepDeltaTime => m_LastCoveredSimTime;
}
