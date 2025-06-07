using System;
using Unity.Collections;
using Unity.Profiling;
using UnityEditor;
using UnityEngine;
using Vector2 = UnityEngine.Vector2;
using Vector3 = UnityEngine.Vector3;

public enum VisualizationMode
{
    Density = 0,
    Pressure = 1,
    PressureGradient = 2,
    NoBackground = 3,
}

public class FluidSimViz : MonoBehaviour
{
    static readonly int k_ParticlePositions = Shader.PropertyToID("ParticlePositions");
    public float m_ScalingFactor = 1.0f;
    public Texture2D m_CircleTexture;
    public Material m_FluidMaterialDebugViz;
    public float m_CircleSize = 10;
    public bool m_ShowVelocities = false;
    public float m_velocityScale = 1.0f;
    public bool m_ShowGrid = false;
    public float m_DensityVizFactor = 1;
    public VisualizationMode m_VisualizationMode = VisualizationMode.Density;
    public Color m_NegativePressureColor;
    public Color m_NeutralPressureColor;
    public Color m_PositivePressureColor;
    public float m_MouseRadius = 1.0f;
    
    FluidSim m_FluidSim;
    
    Texture2D m_Texture;
    Mesh m_Mesh;
    
    float[] m_PointPositionData;
    float[] m_PointDensitiesData;
    float[] m_PointPressureData;
    float[] m_PointVelocityData;
    ComputeBuffer m_PointBuffer;
    ComputeBuffer m_PointDensitiesBuffer;
    ComputeBuffer m_PointPressureBuffer;
    ComputeBuffer m_PointVelocityBuffer;

    float m_KineticEnergy = 0;
    
    Vector2 m_MousePos = Vector2.zero;
    GUIContent energyContent = new GUIContent();
    bool m_MousePressed = false;
    FluidSim.InteractionDirection m_InteractionDirection;

    static readonly ProfilerMarker s_UpdatePerfMarker = new ProfilerMarker("FluidSimViz.Update");
    static readonly ProfilerMarker s_OnGuiPerfMarker = new ProfilerMarker("FluidSimViz.OnGUi");
    
    void Start()
    {
        m_FluidSim = GetComponent<FluidSim>();
        InitTexture(); 
        m_Mesh = new Mesh();
        m_Mesh.SetVertices(new Vector3[]{new (0, 0), new (0, -1), new (1, -1), new (1, 0)});
        m_Mesh.SetIndices(new int[]{0, 2, 1, 0, 3, 2}, MeshTopology.Triangles, 0);
    }

    void OnDisable()
    {
        CleanupComputeBuffers();  
    }

    void OnDestroy()
    {
        CleanupTexture();
        Destroy(m_Mesh);
    }

    void CleanupTexture()
    {
        if (m_Texture != null) 
            Destroy(m_Texture);
    }

    void InitTexture()
    {
        CleanupTexture();
        m_Texture = new Texture2D(m_FluidSim.width, m_FluidSim.height);
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
    }
    
    // Update is called once per frame
    void Update()
    {
        using var markerScope = s_UpdatePerfMarker.Auto();
        if (m_FluidSim.width != m_Texture.width || m_FluidSim.height != m_Texture.height)
        {
            InitTexture();
        }
        
        //var particles = m_FluidSim.GetParticles();
        // var colors = new Color[m_Texture.width * m_Texture.height];
        // for (var x = 0; x < m_Texture.width; x++)
        // {
        //     for (var y = 0; y < m_Texture.height; y++)
        //     {
        //         //var color = new Color(x / (float)m_Texture.width, y / (float)m_Texture.height, 0f,1f);
        //         var density = m_FluidSim.CalculateDensity(new Vector2(x,m_Texture.height-y));
        //         var color = new Color(density, 0f, 0f,1f);
        //         colors[y * m_Texture.width + x] = color;
        //         //Debug.Log(colors[y * m_Texture.width + x]);
        //     }
        // }
        // m_Texture.SetPixels(colors);
        // m_Texture.Apply();

        if (m_PointBuffer == null || m_PointDensitiesBuffer.count != m_FluidSim.m_ParticleCount)
        {
            CleanupComputeBuffers();

            m_PointPositionData = new float[m_FluidSim.m_ParticleCount * 2];
            m_PointBuffer = new ComputeBuffer(m_PointPositionData.Length, sizeof(float));
            m_PointDensitiesData = new float[m_FluidSim.m_ParticleCount];
            m_PointDensitiesBuffer = new ComputeBuffer(m_PointDensitiesData.Length, sizeof(float));
            m_PointPressureData = new float[m_FluidSim.m_ParticleCount];
            m_PointPressureBuffer = new ComputeBuffer(m_PointPressureData.Length, sizeof(float));
            m_PointVelocityData = new float[m_FluidSim.m_ParticleCount*2];
            m_PointVelocityBuffer = new ComputeBuffer(m_PointVelocityData.Length, sizeof(float));
        }

        var positions = m_FluidSim.GetPositions();
        for (var index = 0; index < m_FluidSim.m_ParticleCount; index++)
        {
            var position = positions[index];
            m_PointPositionData[index * 2] = position.x;
            m_PointPositionData[index * 2 + 1] = position.y;
        }
        var densities = m_FluidSim.GetDensities();
        var maxDensity = 0f;
        for (var index = 0; index < m_FluidSim.m_ParticleCount; index++)
        {
            m_PointDensitiesData[index] = densities[index];
            if (densities[index] > maxDensity)
                maxDensity = densities[index];
        }
        
        var pressure = m_FluidSim.GetPressures();
        var maxPressure = pressure[0];
        var minPressure = pressure[0];
        for (var index = 0; index < m_FluidSim.m_ParticleCount; index++)
        {
            m_PointPressureData[index] = pressure[index];
            if(m_PointPressureData[index] > maxPressure)
                maxPressure = m_PointPressureData[index];
            if (m_PointPressureData[index] < minPressure)
                minPressure = m_PointPressureData[index];
        }
        
        var velocities = m_FluidSim.GetVelocities();
        m_KineticEnergy = 0;
        for (var index = 0; index < m_FluidSim.m_ParticleCount; index++)
        {
            var velocity = velocities[index];
            m_PointVelocityData[index * 2] = velocity.x;
            m_PointVelocityData[index * 2 + 1] = velocity.y;
            m_KineticEnergy += 0.5f * m_FluidSim.Mass * Vector2.Dot(velocity, velocity);
        }
        
        m_PointBuffer.SetData(m_PointPositionData);
        m_PointDensitiesBuffer.SetData(m_PointDensitiesData);
        m_PointPressureBuffer.SetData(m_PointPressureData);
        m_FluidMaterialDebugViz.SetBuffer("_particle_positions", m_PointBuffer);
        m_FluidMaterialDebugViz.SetBuffer("_particle_densities", m_PointDensitiesBuffer);
        m_FluidMaterialDebugViz.SetBuffer("_particle_pressures", m_PointPressureBuffer);
        m_FluidMaterialDebugViz.SetInt("_PointCount", m_FluidSim.m_ParticleCount);

        m_FluidMaterialDebugViz.SetFloat("_scaling_factor", m_ScalingFactor);

        var mouseInSimulationSpace = new Vector2(m_MousePos.x / m_ScalingFactor,
        m_FluidSim.height - m_MousePos.y / m_ScalingFactor);

        m_FluidMaterialDebugViz.SetFloat("_sizex", m_FluidSim.width);
        m_FluidMaterialDebugViz.SetFloat("_sizey", m_FluidSim.height);
        m_FluidMaterialDebugViz.SetFloat("_mousex", mouseInSimulationSpace.x);
        m_FluidMaterialDebugViz.SetFloat("_mousey", mouseInSimulationSpace.y);
        m_FluidMaterialDebugViz.SetInt("_mousepressed", m_MousePressed ? 1 :0);
        m_FluidMaterialDebugViz.SetFloat("_mouseradius", m_MouseRadius);
        m_FluidMaterialDebugViz.SetFloat("_smoothingLength", m_FluidSim.SmoothingLength);
        m_FluidMaterialDebugViz.SetFloat("_DensityVizFactor", m_DensityVizFactor);
        m_FluidMaterialDebugViz.SetFloat("_circleSize", m_CircleSize);
        m_FluidMaterialDebugViz.SetInt("_visMode", (int)m_VisualizationMode);
        m_FluidMaterialDebugViz.SetColor("_negativePressureColor", m_NegativePressureColor);
        m_FluidMaterialDebugViz.SetColor("_neutralPressureColor", m_NeutralPressureColor);
        m_FluidMaterialDebugViz.SetColor("_positivePressureColor", m_PositivePressureColor);
        m_FluidMaterialDebugViz.SetFloat("_max_pressure", maxPressure);
        m_FluidMaterialDebugViz.SetFloat("_min_pressure", minPressure);
        m_FluidMaterialDebugViz.SetFloat("_target_density", m_FluidSim.CalcTargetDensity());
        m_FluidMaterialDebugViz.SetFloat("_max_density", maxDensity);
        
        if (m_MousePressed)
            m_FluidSim.Interact(mouseInSimulationSpace, m_MouseRadius / m_ScalingFactor, m_InteractionDirection);
        


        //This code is shit, all I wanted is to draw something screen space, I tried using CommandBuffers, but somethings I needed were not compatible with SRP
        // var worldPoint = Camera.main.ScreenToWorldPoint(new Vector3(0, Camera.main.pixelHeight-1, 1));
        // var worldPointE = Camera.main.ScreenToWorldPoint(new Vector3(Camera.main.pixelWidth-1, 0, 1));
        // var meshWidth = m_FluidSim.width / (float)Camera.main.pixelWidth * (float)(worldPointE.x - worldPoint.x);
        // var meshHeight = m_FluidSim.height / (float)Camera.main.pixelHeight* (float)(worldPoint.y - worldPointE.y);
        // var matrix = Matrix4x4.TRS(worldPoint, Camera.main.transform.rotation, new Vector3(meshWidth, meshHeight, 1));
        // Graphics.DrawMesh(m_Mesh, matrix, m_FluidMaterialDebugViz,0);
    }
    

    void OnGUI()
    {
        using var markerScope = s_OnGuiPerfMarker.Auto();
        
        if (Event.current.type == EventType.MouseDown || Event.current.type == EventType.MouseDrag)
        {
            m_MousePos = Event.current.mousePosition;
            m_MousePressed = true;
            m_InteractionDirection = FluidSim.InteractionDirection.Attract;
            if (Event.current.button == 1)
            {
                m_InteractionDirection = FluidSim.InteractionDirection.Repel;
            }
        }

        if (Event.current.type == EventType.MouseUp)
        {
            m_MousePressed = false;
        }
        
        
        
        if (Event.current.type != EventType.Repaint) return;
        
        DrawMaterial(new Rect(0,0, m_FluidSim.width * m_ScalingFactor, m_FluidSim.height * m_ScalingFactor), m_FluidMaterialDebugViz);
        energyContent.text = $"Kinetic Energy: {m_KineticEnergy}";
        GUI.skin.label.Draw(new Rect(800,0,1000,20), energyContent, 0 );
        

        if (!m_ShowGrid) return;
        var thickness = 1;
        // var topleft = new Vector2(0, 0);
        // var topright = new Vector2(0, m_FluidSim.height * m_ScalingFactor);
        // var bottomleft = new Vector2(m_FluidSim.width*m_ScalingFactor, 0);
        // var bottomright = new Vector2(m_FluidSim.width*m_ScalingFactor, m_FluidSim.height*m_ScalingFactor);
        
        //Graphics.DrawTexture(new Rect(0,0, m_FluidSim.width, m_FluidSim.height), m_Texture);
        var boundaries = new Rect(0, 0, m_FluidSim.width * m_ScalingFactor, m_FluidSim.height * m_ScalingFactor);
        DrawRect(boundaries, thickness, Color.gray);
        // DrawLine(topleft,topright, thickness, Color.gray);
        // DrawLine(topright, bottomright, thickness, Color.gray);
        // DrawLine(bottomright, bottomleft, thickness, Color.gray);
        // DrawLine(bottomleft, topleft, thickness, Color.gray);
        
        var cellSize = m_FluidSim.m_LookupHelper.CellSize;
        var cellsHorizontalCount = m_FluidSim.width / cellSize;
        var cellsVerticalCount = m_FluidSim.height / cellSize;
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
            var bottom = new Vector2(boundaries.yMax, y);
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
        
        var particleIndexes = new NativeList<int>(m_FluidSim.m_ParticleCount, Allocator.Temp);
        m_FluidSim.m_LookupHelper.GetParticlesAround(
            new Vector2(
                selectedGrid.x * cellSize + cellSize * 0.5f,
                (cellsHorizontalCount -1 - selectedGrid.y ) * cellSize + cellSize * 0.5f),
            particleIndexes
            );

        
        if (!m_ShowVelocities) return;

        var positions = m_FluidSim.GetPositions();
        var velocities = m_FluidSim.GetVelocities();

        for (var i = 0; i < m_FluidSim.m_ParticleCount; i++)
        {
            var circleColor = Color.blue;
            if (particleIndexes.Contains(i))
                circleColor = Color.red;
            var position = positions[i] * m_ScalingFactor;
            position.y = m_FluidSim.height* m_ScalingFactor - position.y;
            var velocity = velocities[i];
            velocity.y *= -1;

            DrawCircle(position, m_CircleSize , circleColor);
            DrawArrow(position, velocity * m_velocityScale, Color.green);
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
        //var rect = new Rect(0,0, thickness, thickness);
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
