using System;
using Unity.Collections;
using Unity.Profiling;
using UnityEditor;
using UnityEngine;
using UnityEngine.Rendering;

public enum FieldVisualizationMode
{
    Density = 0,
    Pressure = 1,
    PressureGradient = 2,
    NoBackground = 3,
}

public enum ParticleVisualizationType
{
    Solid = 0,
    Pressure = 1,
    Velocity = 2,
}

public enum ParticleVisualizationMode {
    None = 0,
    DrawInstanced = 1,
    FieldFragment = 2,
    UIOverlay = 3,
}

[Serializable]
public struct VisualizationRanges
{
    public float minPressure;
    public float maxDensity;
    public float maxVelocity;
    public float minVelocity;
    public float maxPressure;
}

public partial class FluidSimViz : MonoBehaviour
{
    public float m_ScalingFactor = 1.0f;
    public Texture2D m_CircleTexture;
    public Shader m_FluidMaterialDebugShader;
    public Shader m_PointInstancedShader;
    public ComputeShader m_StatsComputeShader;
    public MeshRenderer m_FluidMeshRenderer;
    public bool m_UseDynamicRanges = false;
    public float m_CircleSize = 10;
    public bool m_ShowVelocities = false;
    public bool m_ShowGrid = false;
    public ParticleVisualizationMode m_ShowParticles;
    public float m_VelocityScale = 1.0f;
    public float m_DensityVizFactor = 1;
    public FieldVisualizationMode m_FieldVisualizationMode = FieldVisualizationMode.Density;
    public ParticleVisualizationType m_ParticleVisualizationMode = ParticleVisualizationType.Pressure;
    public Color m_NegativePressureColor;
    public Color m_NeutralPressureColor;
    public Color m_PositivePressureColor;
    public Color m_FluidParticleColor;
    public Color m_ExternalPointsColor;
    public float m_MouseRadius = 1.0f;
    public VisualizationRanges m_Ranges = new VisualizationRanges();

    Material m_FluidMaterialDebugViz;

    IFluidSim m_FluidSim;
    RectTransform m_RectTransform;
    
    Vector2[] m_PointPositionData;
    float[] m_PointDensitiesData;
    float[] m_PointPressureData;
    Vector2[] m_PointVelocityData;
    ComputeBuffer m_PointBuffer;
    ComputeBuffer m_PointDensitiesBuffer;
    ComputeBuffer m_PointPressureBuffer;
    ComputeBuffer m_PointVelocityBuffer;

    // Matches the FluidStats struct in FluidStats.compute.
    struct FluidStatsGPU
    {
        public float maxDensity;
        public float minPressure;
        public float maxPressure;
        public float minVelocitySq;
        public float maxVelocitySq;
        public float velocitySqSum;
    }

    const int k_StatsGroupSize = 256;
    const int k_StatsGroupCount = 64;

    float m_KineticEnergy = 0;
    VisualizationRanges m_CachedDynamicRanges;
    bool m_StatsRequestInFlight;
    AsyncGPUReadbackRequest m_StatsRequest;
    ComputeBuffer m_StatsPartialsBuffer;
    ComputeBuffer m_StatsResultBuffer;
    GUIContent energyContent = new GUIContent();
    Vector3[] m_FourCornersWorldSpace = new Vector3[4];
    Vector3[] m_FourCornersScreenSpace = new Vector3[4];
    
    Vector2 m_MousePos = Vector2.zero;
    bool m_MousePressed = false;
    InteractionDirection m_InteractionDirection;
    
    Material m_PointInstancedMaterial;
    Mesh m_PointMesh;
    DrawSimPointsInstanced m_DrawSimPointsInstanced;
    DrawSimPointsInstanced m_DrawSimPointsInstanced_external;


    static readonly ProfilerMarker s_UpdatePerfMarker = new ProfilerMarker("FluidSimViz.Update");
    static readonly ProfilerMarker s_OnGuiPerfMarker = new ProfilerMarker("FluidSimViz.OnGUi");
    
    void Start()
    {
        m_RectTransform = GetComponent<RectTransform>();
    }

    void OnEnable()
    {
        m_FluidSim = GetComponent<IFluidSim>() as IFluidSim;
        m_PointInstancedMaterial = new Material(m_PointInstancedShader);
        m_FluidMaterialDebugViz = new Material(m_FluidMaterialDebugShader);
        m_FluidMeshRenderer.material = m_FluidMaterialDebugViz;
        var statsStride = System.Runtime.InteropServices.Marshal.SizeOf(typeof(FluidStatsGPU));
        m_StatsPartialsBuffer = new ComputeBuffer(k_StatsGroupCount, statsStride);
        m_StatsPartialsBuffer.name = "FluidStatsPartialsBuffer";
        m_StatsResultBuffer = new ComputeBuffer(1, statsStride);
        m_StatsResultBuffer.name = "FluidStatsResultBuffer";
        m_DrawSimPointsInstanced.InitializePointInstancedRender(m_FluidSim.ParticleCount, m_PointInstancedMaterial);
        m_DrawSimPointsInstanced_external.InitializePointInstancedRender(0, m_PointInstancedMaterial);
    }

    void OnDisable()
    {
        m_StatsPartialsBuffer?.Release();
        m_StatsResultBuffer?.Release();
        m_StatsRequestInFlight = false;
        m_DrawSimPointsInstanced.CleanupPointInstancedRender();
        m_DrawSimPointsInstanced_external.CleanupPointInstancedRender();
        m_FluidMeshRenderer.material = null;
        Destroy(m_FluidMaterialDebugViz);
        Destroy(m_PointInstancedMaterial);
    }

    void OnDestroy() {}

    
    // Update is called once per frame
    void Update()
    {
        using var markerScope = s_UpdatePerfMarker.Auto();
        m_RectTransform.GetWorldCorners(m_FourCornersWorldSpace);
        var cam = Camera.main;
        m_FourCornersScreenSpace[0] = cam.WorldToScreenPoint(m_FourCornersWorldSpace[0]);
        m_FourCornersScreenSpace[1] = cam.WorldToScreenPoint(m_FourCornersWorldSpace[1]);
        m_FourCornersScreenSpace[2] = cam.WorldToScreenPoint(m_FourCornersWorldSpace[2]);
        m_FourCornersScreenSpace[3] = cam.WorldToScreenPoint(m_FourCornersWorldSpace[3]);

        if (m_PointPositionData == null || m_PointPositionData.Length != m_FluidSim.ParticleCount)
        {
            m_PointBuffer = m_FluidSim.GetPositionComputeBuffer();
            m_PointDensitiesBuffer = m_FluidSim.GetDensities();;
            m_PointPressureBuffer = m_FluidSim.GetPressures();
            m_PointVelocityBuffer = m_FluidSim.GetVelocities();
            
            m_PointPositionData = new Vector2[m_FluidSim.ParticleCount];
            m_PointDensitiesData = new float[m_FluidSim.ParticleCount];
            m_PointPressureData = new float[m_FluidSim.ParticleCount];
            m_PointVelocityData = new Vector2[m_FluidSim.ParticleCount];
        }
        var ranges = m_Ranges;
        if (m_UseDynamicRanges)
        {
            ranges = CalculateStats();
        }
        
        var mouseInSimulationSpace = CalcMouseInSimulationSpace();
        var mouseRadiusRatio  =  m_FluidSim.Width / (m_FourCornersScreenSpace[2].x - m_FourCornersScreenSpace[0].x);
        var adjustedRadius = m_MouseRadius * mouseRadiusRatio;

        m_FluidMaterialDebugViz.SetBuffer("_particle_positions", m_PointBuffer);
        m_FluidMaterialDebugViz.SetBuffer("_particle_densities", m_PointDensitiesBuffer);
        m_FluidMaterialDebugViz.SetBuffer("_particle_pressures", m_PointPressureBuffer);
        m_FluidMaterialDebugViz.SetBuffer("_particle_velocities", m_PointVelocityBuffer);
        m_FluidMaterialDebugViz.SetInt("_PointCount", m_FluidSim.ParticleCount);
        m_FluidMaterialDebugViz.SetFloat("_mass", m_FluidSim.Mass);
        m_FluidMaterialDebugViz.SetFloat("_scaling_factor", m_ScalingFactor);
        m_FluidMaterialDebugViz.SetFloat("_sim_width", m_FluidSim.Width);
        m_FluidMaterialDebugViz.SetFloat("_sim_height", m_FluidSim.Height);
        m_FluidMaterialDebugViz.SetFloat("_mousex", mouseInSimulationSpace.x);
        m_FluidMaterialDebugViz.SetFloat("_mousey", mouseInSimulationSpace.y);
        m_FluidMaterialDebugViz.SetInt("_mousepressed", m_MousePressed ? 1 :0);
        m_FluidMaterialDebugViz.SetFloat("_mouseradius", adjustedRadius);
        m_FluidMaterialDebugViz.SetFloat("_smoothingLength", m_FluidSim.SmoothingRadius);
        m_FluidMaterialDebugViz.SetFloat("_DensityVizFactor", m_DensityVizFactor);
        m_FluidMaterialDebugViz.SetFloat("_circleSize", m_ShowParticles == ParticleVisualizationMode.FieldFragment ? m_CircleSize: 0);
        m_FluidMaterialDebugViz.SetInt("_fieldVisMode", (int)m_FieldVisualizationMode);
        m_FluidMaterialDebugViz.SetInt("_particleVisMode", (int)m_ParticleVisualizationMode);
        m_FluidMaterialDebugViz.SetColor("_negativePressureColor", m_NegativePressureColor);
        m_FluidMaterialDebugViz.SetColor("_neutralPressureColor", m_NeutralPressureColor);
        m_FluidMaterialDebugViz.SetColor("_positivePressureColor", m_PositivePressureColor);
        m_FluidMaterialDebugViz.SetFloat("_target_density", m_FluidSim.TargetDensity);
        m_FluidMaterialDebugViz.SetFloat("_max_pressure", ranges.maxPressure);
        m_FluidMaterialDebugViz.SetFloat("_min_pressure", ranges.minPressure);
        m_FluidMaterialDebugViz.SetFloat("_max_density", ranges.maxDensity);
        m_FluidMaterialDebugViz.SetFloat("_max_velocity", ranges.maxVelocity);
        m_FluidMaterialDebugViz.SetFloat("_min_velocity", ranges.minVelocity);

        if (m_ShowParticles == ParticleVisualizationMode.DrawInstanced)
        {
            m_DrawSimPointsInstanced.RenderPointsGPU(ranges.maxVelocity, m_FluidSim.ParticleCount,
                m_PointBuffer, m_PointVelocityBuffer, m_PointDensitiesBuffer, m_PointPressureBuffer,
                (int)m_ParticleVisualizationMode, m_FourCornersWorldSpace[0], m_RectTransform.rect.width, m_RectTransform.rect.height,
                m_FluidSim.Width, m_FluidSim.Height, m_CircleSize, m_ScalingFactor, m_FluidParticleColor);
            // Null until a rigid body registers surface points (and always null
            // on the GPU sim, which has no coupling yet).
            var externalPoints = m_FluidSim.InputExternalPoints;
            if (externalPoints != null)
            {
                m_DrawSimPointsInstanced_external.RenderPointsGPU(ranges.maxVelocity, externalPoints.count,
                    externalPoints, m_PointVelocityBuffer, m_PointDensitiesBuffer, m_PointPressureBuffer,
                    (int)m_ParticleVisualizationMode, m_FourCornersWorldSpace[0], m_RectTransform.rect.width, m_RectTransform.rect.height,
                    m_FluidSim.Width, m_FluidSim.Height, m_CircleSize, m_ScalingFactor, m_ExternalPointsColor);
            }
        }

        if (m_MousePressed)
        {
            m_FluidSim.Interact(mouseInSimulationSpace, adjustedRadius, m_InteractionDirection);
        }
    }

    // The min/max/sum reduction runs on the GPU (FluidStats.compute), so only
    // a single small struct is ever read back, asynchronously. The ranges only
    // normalize the visualization, so the few frames of readback latency are
    // invisible; a synchronous readback here would drain the whole queued GPU
    // pipeline (all sim substeps of the frame) and dominate frame time.
    VisualizationRanges CalculateStats()
    {
        if (m_StatsComputeShader == null)
            return m_CachedDynamicRanges;

        if (m_StatsRequestInFlight && m_StatsRequest.done)
        {
            m_StatsRequestInFlight = false;
            if (!m_StatsRequest.hasError)
            {
                var stats = m_StatsRequest.GetData<FluidStatsGPU>()[0];
                m_KineticEnergy = 0.5f * m_FluidSim.Mass * stats.velocitySqSum;
                m_CachedDynamicRanges = new VisualizationRanges
                {
                    maxDensity = stats.maxDensity,
                    minPressure = stats.minPressure,
                    maxPressure = stats.maxPressure,
                    minVelocity = Mathf.Sqrt(Mathf.Max(0f, stats.minVelocitySq)),
                    maxVelocity = Mathf.Sqrt(Mathf.Max(1f, stats.maxVelocitySq)),
                };
            }
        }

        if (!m_StatsRequestInFlight)
        {
            var reduceKernel = m_StatsComputeShader.FindKernel("ReduceStats");
            var finalKernel = m_StatsComputeShader.FindKernel("ReduceStatsFinal");
            m_StatsComputeShader.SetInt("ParticleCount", m_FluidSim.ParticleCount);
            m_StatsComputeShader.SetInt("GroupCount", k_StatsGroupCount);
            m_StatsComputeShader.SetBuffer(reduceKernel, "Densities", m_PointDensitiesBuffer);
            m_StatsComputeShader.SetBuffer(reduceKernel, "Pressures", m_PointPressureBuffer);
            m_StatsComputeShader.SetBuffer(reduceKernel, "Velocities", m_PointVelocityBuffer);
            m_StatsComputeShader.SetBuffer(reduceKernel, "Partials", m_StatsPartialsBuffer);
            m_StatsComputeShader.Dispatch(reduceKernel, k_StatsGroupCount, 1, 1);
            m_StatsComputeShader.SetBuffer(finalKernel, "Partials", m_StatsPartialsBuffer);
            m_StatsComputeShader.SetBuffer(finalKernel, "Result", m_StatsResultBuffer);
            m_StatsComputeShader.Dispatch(finalKernel, 1, 1, 1);

            m_StatsRequest = AsyncGPUReadback.Request(m_StatsResultBuffer);
            m_StatsRequestInFlight = true;
        }

        return m_CachedDynamicRanges;
    }
    
    Vector2 CalcMouseInSimulationSpace()
    {

        var topLeftCameraSpace = m_FourCornersScreenSpace[0];
        var bottomRightCameraSpace = m_FourCornersScreenSpace[2];
        var tx = Mathf.InverseLerp(topLeftCameraSpace.x, bottomRightCameraSpace.x, m_MousePos.x);
        //The top and bottom are flipped here, because screen space is flipped in relation to simulation space
        var mousey = Camera.main.pixelHeight - m_MousePos.y;
        var ty = Mathf.InverseLerp( topLeftCameraSpace.y, bottomRightCameraSpace.y, mousey);
        var mouseInSimulationSpace = new Vector2(m_FluidSim.Width * tx, m_FluidSim.Height * ty);
        return mouseInSimulationSpace;
    }

    void OnGUI()
    {
        using var markerScope = s_OnGuiPerfMarker.Auto();
        
        if (Event.current.type == EventType.MouseDown || Event.current.type == EventType.MouseDrag)
        {
            m_MousePos = Event.current.mousePosition;
            m_MousePressed = true;
            m_InteractionDirection = InteractionDirection.Attract;
            if (Event.current.button == 1)
            {
                m_InteractionDirection = InteractionDirection.Repel;
            }
        }

        if (Event.current.type == EventType.MouseUp)
        {
            m_MousePressed = false;
        }
        // OnGUI does not receive MouseUp when the release happens outside the
        // game view, which left the interaction force latched on permanently.
        if (m_MousePressed && !Input.GetMouseButton(0) && !Input.GetMouseButton(1))
        {
            m_MousePressed = false;
        }
        
        if (Event.current.type != EventType.Repaint) return;

        // if (m_FieldVisualizationMode != FieldVisualizationMode.NoBackground)
        // {
        //     var fourCorners = new Vector3[4];
        //     m_RectTransform.GetWorldCorners(fourCorners);
        //     var screenSpace = Camera.main.WorldToScreenPoint(fourCorners[0]);
        //     DrawMaterial(new Rect(screenSpace.x, screenSpace.y, m_FluidSim.Width * m_ScalingFactor, m_FluidSim.Height * m_ScalingFactor), m_FluidMaterialDebugViz);
        // }
        energyContent.text = $"Kinetic Energy: {m_KineticEnergy}";
        GUI.skin.label.Draw(new Rect(800,0,1000,20), energyContent, 0 );

        var showUIParticles = m_ShowParticles == ParticleVisualizationMode.UIOverlay;
#if UNITY_EDITOR
        if (!m_ShowGrid && !m_ShowVelocities && !showUIParticles) return;
        
        var thickness = 1;
        // var topleft = new Vector2(0, 0);
        // var topright = new Vector2(0, m_FluidSim.height * m_ScalingFactor);
        // var bottomleft = new Vector2(m_FluidSim.width*m_ScalingFactor, 0);
        // var bottomright = new Vector2(m_FluidSim.width*m_ScalingFactor, m_FluidSim.height*m_ScalingFactor);
        
        var boundaries = new Rect(0, 0, m_FluidSim.Width * m_ScalingFactor, m_FluidSim.Height * m_ScalingFactor);
        DrawRect(boundaries, thickness, Color.gray);
        // DrawLine(topleft,topright, thickness, Color.gray);
        // DrawLine(topright, bottomright, thickness, Color.gray);
        // DrawLine(bottomright, bottomleft, thickness, Color.gray);
        // DrawLine(bottomleft, topleft, thickness, Color.gray);
        
        var cellSize = m_FluidSim.LookupHelper.CellSize;
        var cellsHorizontalCount = m_FluidSim.Width / cellSize;
        var cellsVerticalCount = m_FluidSim.Height / cellSize;
        for (int i = 0; i < cellsHorizontalCount; i++)
        {
            var x = cellSize * i * m_ScalingFactor;
            var top = new Vector2(x, 0);
            var bottom = new Vector2(x, boundaries.yMax);
            DrawLine(top, bottom, thickness, Color.gray);
        }
        for (int i = 0; i < cellsVerticalCount; i++)
        {
            var y = cellSize * i * m_ScalingFactor;
            var top = new Vector2(0, y);
            var bottom = new Vector2(boundaries.xMax, y);
            DrawLine(top, bottom, thickness, Color.gray);
        }

        var selectedGrid = new Vector2Int(-1, -1);
        selectedGrid.x = Mathf.FloorToInt(m_MousePos.x / m_ScalingFactor / cellSize);
        selectedGrid.y = Mathf.FloorToInt(m_MousePos.y / m_ScalingFactor / cellSize);
        var cellRect = new Rect(
            selectedGrid.x * cellSize * m_ScalingFactor,
            selectedGrid.y * cellSize * m_ScalingFactor,
            cellSize * m_ScalingFactor,
            cellSize * m_ScalingFactor
            );
        
        DrawRect(cellRect, thickness*2, Color.red);
        
        var particleIndexes = new NativeList<int>(m_FluidSim.ParticleCount, Allocator.Temp);
        m_FluidSim.LookupHelper.GetParticlesAround(
            new Vector2(
                selectedGrid.x * cellSize + cellSize * 0.5f,
                (cellsVerticalCount -1 - selectedGrid.y ) * cellSize + cellSize * 0.5f),
            particleIndexes
            );

        if (!m_ShowVelocities && !showUIParticles) return;

        m_FluidSim.GetPositionComputeBuffer().GetData(m_PointPositionData);
        m_FluidSim.GetVelocities().GetData(m_PointVelocityData);

        for (var i = 0; i < m_FluidSim.ParticleCount; i++)
        {
            var circleColor = Color.blue;
            if (particleIndexes.Contains(i))
                circleColor = Color.red;
            var position = m_PointPositionData[i] * m_ScalingFactor;
            position.y = m_FluidSim.Height* m_ScalingFactor - position.y;
            var velocity = m_PointVelocityData[i];
            velocity.y *= -1;
            if (showUIParticles)
                DrawCircle(position, m_CircleSize , circleColor);
            if (m_ShowVelocities)
                DrawArrow(position, velocity * m_VelocityScale, Color.green);
        }

        particleIndexes.Dispose();
        #endif
    }
#if UNITY_EDITOR
    static void DrawRect(Rect rect, float thickness, Color color)
    {
        var topLeft = new Vector2(rect.x, rect.y);
        var topRight = new Vector2(rect.xMax, rect.y);
        var bottomLeft = new Vector2(rect.x, rect.yMax);
        var bottomRight = new Vector2(rect.xMax, rect.yMax);
        DrawLine(topLeft, topRight, thickness, color);
        DrawLine(topRight, bottomRight, thickness, color);
        DrawLine(bottomRight, bottomLeft, thickness, color);
        DrawLine(bottomLeft, topLeft, thickness, color);
    }

    static void DrawMaterial(Rect rect, Material mat)
    {
        Graphics.DrawTexture(rect, Texture2D.grayTexture, mat);
    }

    void DrawCircle(Vector2 center, float radius, Color color)
    {
        var originalColor = GUI.color;
        GUI.color = color;
        var rect = new Rect(center.x - radius , center.y - radius, radius+radius, radius+radius);
        GUI.DrawTexture(rect, m_CircleTexture);
        GUI.color = originalColor;
    }

    static void DrawLine(Vector2 start, Vector2 end, float thickness, Color color)
    {
        var prevColor = Handles.color;
        Handles.color = color;
        Handles.DrawLine(start, end, thickness);
        Handles.color = prevColor;
    }

    static void DrawArrow(Vector2 position, Vector2 direction, Color color)
    {
        var prevColor = Handles.color;
        Handles.color = color;
        // var newDir = Quaternion.AngleAxis(45.0f, Vector3.forward) * - direction.normalized;
        // var newDir2 = new Vector2(newDir.x, newDir.y);
        // Handles.DrawLine(position+direction, position+direction+newDir2, 1.0f );
        Handles.DrawLine(position, position+direction, 2.0f);
        Handles.color = prevColor;
    }
#endif
 
}
