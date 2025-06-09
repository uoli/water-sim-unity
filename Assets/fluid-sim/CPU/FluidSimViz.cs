using System;
using System.Runtime.InteropServices;
using Unity.Collections;
using Unity.Profiling;
using UnityEditor;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Serialization;
using Random = UnityEngine.Random;

public enum FieldVisualizationMode
{
    Density = 0,
    Pressure = 1,
    PressureGradient = 2,
    NoBackground = 3,
}

public enum ParticleVisualizationMode
{
    Solid = 0,
    Pressure = 1,
    Velocity = 2,
}

public class FluidSimViz : MonoBehaviour
{
    public float m_ScalingFactor = 1.0f;
    public Texture2D m_CircleTexture;
    public Material m_FluidMaterialDebugViz;
    public Material m_PointInstancedMaterial;

    public float m_CircleSize = 10;
    public bool m_ShowVelocities = false;
    public bool m_ShowGrid = false;
    public bool m_ShowParticleOverlay = false;
    [FormerlySerializedAs("m_ShowGpuParticles")]
    public bool m_ShowParticles;
    public float m_VelocityScale = 1.0f;
    public float m_DensityVizFactor = 1;
    public FieldVisualizationMode m_FieldVisualizationMode = FieldVisualizationMode.Density;
    public ParticleVisualizationMode m_ParticleVisualizationMode = ParticleVisualizationMode.Pressure;
    public Color m_NegativePressureColor;
    public Color m_NeutralPressureColor;
    public Color m_PositivePressureColor;
    public float m_MouseRadius = 1.0f;
    
    Mesh m_PointMesh;
    ComputeBuffer m_PointRenderArgsBuffer;
    uint[] m_PointRenderArgsData;
    
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

    float m_KineticEnergy = 0;
    GUIContent energyContent = new GUIContent();
    
    Vector2 m_MousePos = Vector2.zero;
    bool m_MousePressed = false;
    InteractionDirection m_InteractionDirection;

    static readonly ProfilerMarker s_UpdatePerfMarker = new ProfilerMarker("FluidSimViz.Update");
    static readonly ProfilerMarker s_OnGuiPerfMarker = new ProfilerMarker("FluidSimViz.OnGUi");
    
    void Start()
    {
        m_RectTransform = GetComponent<RectTransform>();
    }

    void OnEnable()
    {
        m_FluidSim = GetComponent<IFluidSim>() as IFluidSim;
        InitializePointInstancedRender();
    }

    void OnDisable()
    {
        CleanupComputeBuffers();
        CleanupPointInstancedRender();
    }

    void OnDestroy() {}

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
    }
    
    // Update is called once per frame
    void Update()
    {
        using var markerScope = s_UpdatePerfMarker.Auto();

        if (m_PointBuffer == null || m_PointDensitiesBuffer.count != m_FluidSim.ParticleCount)
        {
            CleanupComputeBuffers();
            
            m_PointPositionData = new Vector2[m_FluidSim.ParticleCount];
            m_PointBuffer = new ComputeBuffer(m_PointPositionData.Length, Marshal.SizeOf(typeof(Vector2)));
            m_PointDensitiesData = new float[m_FluidSim.ParticleCount];
            m_PointDensitiesBuffer = new ComputeBuffer(m_PointDensitiesData.Length, sizeof(float));
            m_PointPressureData = new float[m_FluidSim.ParticleCount];
            m_PointPressureBuffer = new ComputeBuffer(m_PointPressureData.Length, sizeof(float));
            m_PointVelocityData = new Vector2[m_FluidSim.ParticleCount];
            m_PointVelocityBuffer = new ComputeBuffer(m_PointVelocityData.Length, Marshal.SizeOf(typeof(Vector2)));
        }

        var positions = m_FluidSim.GetPositions();
        for (var index = 0; index < m_FluidSim.ParticleCount; index++)
        {
            var position = positions[index];
            m_PointPositionData[index] = position;
        }

        var densities = m_FluidSim.GetDensities();
        var maxDensity = 0f;
        for (var index = 0; index < m_FluidSim.ParticleCount; index++)
        {
            m_PointDensitiesData[index] = densities[index];
            if (densities[index] > maxDensity)
                maxDensity = densities[index];
        }
        
        var pressure = m_FluidSim.GetPressures();
        var maxPressure = pressure[0];
        var minPressure = pressure[0];
        for (var index = 0; index < m_FluidSim.ParticleCount; index++)
        {
            m_PointPressureData[index] = pressure[index];
            if(m_PointPressureData[index] > maxPressure)
                maxPressure = m_PointPressureData[index];
            if (m_PointPressureData[index] < minPressure)
                minPressure = m_PointPressureData[index];
        }
        
        var velocities = m_FluidSim.GetVelocities();
        m_KineticEnergy = 0;
        var minVelocity = 0f;
        var maxVelocity = 0f;
        for (var index = 0; index < m_FluidSim.ParticleCount; index++)
        {
            var velocity = velocities[index];
            m_PointVelocityData[index] = velocity;
            m_KineticEnergy += 0.5f * m_FluidSim.Mass * Vector2.Dot(velocity, velocity);
            minVelocity = Mathf.Min(velocity.sqrMagnitude, minVelocity);
            maxVelocity = Mathf.Max(velocity.sqrMagnitude, maxVelocity);
        }
        minVelocity = Mathf.Sqrt(minVelocity);
        maxVelocity = Mathf.Sqrt(maxVelocity);
        
        m_PointBuffer.SetData(m_PointPositionData);
        m_PointDensitiesBuffer.SetData(m_PointDensitiesData);
        m_PointPressureBuffer.SetData(m_PointPressureData);
        m_PointVelocityBuffer.SetData(m_PointVelocityData);
        m_FluidMaterialDebugViz.SetBuffer("_particle_positions", m_PointBuffer);
        m_FluidMaterialDebugViz.SetBuffer("_particle_densities", m_PointDensitiesBuffer);
        m_FluidMaterialDebugViz.SetBuffer("_particle_pressures", m_PointPressureBuffer);
        m_FluidMaterialDebugViz.SetBuffer("_particle_velocities", m_PointVelocityBuffer);
        m_FluidMaterialDebugViz.SetInt("_PointCount", m_FluidSim.ParticleCount);

        m_FluidMaterialDebugViz.SetFloat("_scaling_factor", m_ScalingFactor);

        var mouseInSimulationSpace = new Vector2(m_MousePos.x / m_ScalingFactor,
        m_FluidSim.Height - m_MousePos.y / m_ScalingFactor);

        m_FluidMaterialDebugViz.SetFloat("_sim_width", m_FluidSim.Width);
        m_FluidMaterialDebugViz.SetFloat("_sim_height", m_FluidSim.Height);
        m_FluidMaterialDebugViz.SetFloat("_mousex", mouseInSimulationSpace.x);
        m_FluidMaterialDebugViz.SetFloat("_mousey", mouseInSimulationSpace.y);
        m_FluidMaterialDebugViz.SetInt("_mousepressed", m_MousePressed ? 1 :0);
        m_FluidMaterialDebugViz.SetFloat("_mouseradius", m_MouseRadius);
        m_FluidMaterialDebugViz.SetFloat("_smoothingLength", m_FluidSim.SmoothingRadius);
        m_FluidMaterialDebugViz.SetFloat("_DensityVizFactor", m_DensityVizFactor);
        m_FluidMaterialDebugViz.SetFloat("_circleSize", m_CircleSize);
        m_FluidMaterialDebugViz.SetInt("_fieldVisMode", (int)m_FieldVisualizationMode);
        m_FluidMaterialDebugViz.SetInt("_particleVisMode", (int)m_ParticleVisualizationMode);
        m_FluidMaterialDebugViz.SetColor("_negativePressureColor", m_NegativePressureColor);
        m_FluidMaterialDebugViz.SetColor("_neutralPressureColor", m_NeutralPressureColor);
        m_FluidMaterialDebugViz.SetColor("_positivePressureColor", m_PositivePressureColor);
        m_FluidMaterialDebugViz.SetFloat("_max_pressure", maxPressure);
        m_FluidMaterialDebugViz.SetFloat("_min_pressure", minPressure);
        m_FluidMaterialDebugViz.SetFloat("_target_density", m_FluidSim.TargetDensity);
        m_FluidMaterialDebugViz.SetFloat("_max_density", maxDensity);
        m_FluidMaterialDebugViz.SetFloat("_max_velocity", maxVelocity);
        m_FluidMaterialDebugViz.SetFloat("_min_velocity", minVelocity);

        if (m_ShowParticles)
            RenderPointsGPU(maxVelocity);
        
        if (m_MousePressed)
            m_FluidSim.Interact(mouseInSimulationSpace, m_MouseRadius / m_ScalingFactor, m_InteractionDirection);
    }

    void InitializePointInstancedRender()
    {
        m_PointInstancedMaterial.enableInstancing = true;
        m_PointMesh = new Mesh();
        m_PointMesh.SetVertices(new Vector3[]{ new(-0.5f,0.5f,0), new(0.5f,0.5f,0), new(0.5f,-0.5f,0), new(-0.5f,-0.5f,0)});
        m_PointMesh.SetIndices(new []{ 0,1,2, 2,3,0}, MeshTopology.Triangles, 0);
        m_PointMesh.Optimize();
        m_PointRenderArgsData = new uint[]
        {
            m_PointMesh.GetIndexCount(0),
            (uint)m_FluidSim.ParticleCount,
            m_PointMesh.GetIndexStart(0),
            m_PointMesh.GetBaseVertex(0),
            0
        };
        m_PointRenderArgsBuffer = new ComputeBuffer(1, m_PointRenderArgsData.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
    }

    void CleanupPointInstancedRender()
    {
        m_PointRenderArgsBuffer.Dispose();
        Destroy(m_PointMesh);
    }

    void RenderPointsGPU(float maxVelocity)
    {
        
        m_PointRenderArgsData[1] = (uint)m_FluidSim.ParticleCount;
        m_PointRenderArgsBuffer.SetData(m_PointRenderArgsData);

        var fourCorners = new Vector3[4];
        m_RectTransform.GetWorldCorners(fourCorners);
        
        m_PointInstancedMaterial.SetBuffer("positions", m_PointBuffer); 
        m_PointInstancedMaterial.SetBuffer("_particle_densities", m_PointDensitiesBuffer);
        m_PointInstancedMaterial.SetBuffer("_particle_pressures", m_PointPressureBuffer);
        m_PointInstancedMaterial.SetBuffer("_particle_velocities", m_PointVelocityBuffer);
        m_PointInstancedMaterial.SetFloat("_max_velocity", maxVelocity);
        m_PointInstancedMaterial.SetInt("_particleVisMode", (int)m_ParticleVisualizationMode);
            
        m_PointInstancedMaterial.SetVector("world_origin", fourCorners[0]);
        m_PointInstancedMaterial.SetFloat("display_area_width", m_RectTransform.rect.width);
        m_PointInstancedMaterial.SetFloat("display_area_height", m_RectTransform.rect.height);

        m_PointInstancedMaterial.SetInt("point_count", m_FluidSim.ParticleCount);
        m_PointInstancedMaterial.SetFloat("sim_width", m_FluidSim.Width);
        m_PointInstancedMaterial.SetFloat("sim_height", m_FluidSim.Height);
        m_PointInstancedMaterial.SetFloat("circle_size", m_CircleSize);
        m_PointInstancedMaterial.SetFloat("scaling_factor", m_ScalingFactor);
        
        var bounds = new Bounds(Vector3.zero, 10000*Vector3.one); // use tighter bounds
        Graphics.DrawMeshInstancedIndirect(m_PointMesh, 0, m_PointInstancedMaterial, bounds, m_PointRenderArgsBuffer);
        
        // CommandBuffer cmd = new CommandBuffer();
        // cmd.SetRenderTarget(m_PointRenderTexture);
        // //cmd.ClearRenderTarget(true, true, Color.red);
        // Graphics.ExecuteCommandBuffer(cmd);
        // cmd.Release();
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
        
        if (!m_ShowGrid) return;
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
        
        if (!m_ShowVelocities) return;

        var positions = m_FluidSim.GetPositions();
        var velocities = m_FluidSim.GetVelocities();

        for (var i = 0; i < m_FluidSim.ParticleCount; i++)
        {
            var circleColor = Color.blue;
            if (particleIndexes.Contains(i))
                circleColor = Color.red;
            var position = positions[i] * m_ScalingFactor;
            position.y = m_FluidSim.Height* m_ScalingFactor - position.y;
            var velocity = velocities[i];
            velocity.y *= -1;
            if (m_ShowParticleOverlay)
                DrawCircle(position, m_CircleSize , circleColor);
            DrawArrow(position, velocity * m_VelocityScale, Color.green);
        }

        particleIndexes.Dispose();
    }

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

 
}
