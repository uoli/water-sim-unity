using System;
using System.Collections.Generic;
using Unity.Profiling;
using UnityEditor;
using UnityEngine;

public enum VisualizationMode
{
    Density = 0,
    Pressure = 1,
}

public class FluidSimViz : MonoBehaviour
{
    static readonly int k_ParticlePositions = Shader.PropertyToID("ParticlePositions");
    public Texture2D m_CircleTexture;
    public Material m_CircleMaterial;
    public Material m_FluidMaterialDebugViz;
    public int m_CircleSize = 10;
    public bool m_ShowOverlay = false;
    public float m_DensityVizFactor = 1;
    public VisualizationMode m_VisualizationMode = VisualizationMode.Density;
    public Color m_NegativePressureColor;
    public Color m_NeutralPressureColor;
    public Color m_PositivePressureColor;
    
    FluidSim m_FluidSim;
    
    Texture2D m_Texture;
    Mesh m_Mesh;
    
    float[] m_PointPositionData;
    float[] m_PointDensitiesData;
    float[] m_PointPressureData;
    ComputeBuffer m_PointBuffer;
    ComputeBuffer m_PointDensitiesBuffer;
    ComputeBuffer m_PointPressureBuffer;
    
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
        }

        var positions = m_FluidSim.GetPositions();
        for (var index = 0; index < m_FluidSim.m_ParticleCount; index++)
        {
            var position = positions[index];
            m_PointPositionData[index * 2] = position.x;
            m_PointPositionData[index * 2 + 1] = position.y;
        }
        var densities = m_FluidSim.GetDensities();
        for (var index = 0; index < m_FluidSim.m_ParticleCount; index++)
        {
            m_PointDensitiesData[index] = densities[index];
        }
        var pressure = m_FluidSim.GetPressures();
        for (var index = 0; index < m_FluidSim.m_ParticleCount; index++)
        {
            m_PointPressureData[index] = pressure[index];
        }
        
        m_PointBuffer.SetData(m_PointPositionData);
        m_PointDensitiesBuffer.SetData(m_PointDensitiesData);
        m_PointPressureBuffer.SetData(m_PointPressureData);
        m_FluidMaterialDebugViz.SetBuffer("_particle_positions", m_PointBuffer);
        m_FluidMaterialDebugViz.SetBuffer("_particle_densities", m_PointDensitiesBuffer);
        m_FluidMaterialDebugViz.SetBuffer("_particle_pressures", m_PointPressureBuffer);
        m_FluidMaterialDebugViz.SetInt("_PointCount", m_FluidSim.m_ParticleCount);
        
        m_FluidMaterialDebugViz.SetFloat("_sizex", m_FluidSim.width);
        m_FluidMaterialDebugViz.SetFloat("_sizey", m_FluidSim.height);
        m_FluidMaterialDebugViz.SetFloat("_mousex", m_MousePos.x);
        m_FluidMaterialDebugViz.SetFloat("_mousey", m_MousePos.y);
        m_FluidMaterialDebugViz.SetFloat("_smoothingLength", m_FluidSim.SmoothingLength);
        m_FluidMaterialDebugViz.SetFloat("_DensityVizFactor", m_DensityVizFactor);
        m_FluidMaterialDebugViz.SetFloat("_circleSize", m_CircleSize);
        m_FluidMaterialDebugViz.SetInt("_visMode", (int)m_VisualizationMode);
        m_FluidMaterialDebugViz.SetColor("_negativePressureColor", m_NegativePressureColor);
        m_FluidMaterialDebugViz.SetColor("_neutralPressureColor", m_NeutralPressureColor);
        m_FluidMaterialDebugViz.SetColor("_positivePressureColor", m_PositivePressureColor);
        


        //This code is shit, all I wanted is to draw something screen space, I tried using CommandBuffers, but somethings I needed were not compatible with SRP
        // var worldPoint = Camera.main.ScreenToWorldPoint(new Vector3(0, Camera.main.pixelHeight-1, 1));
        // var worldPointE = Camera.main.ScreenToWorldPoint(new Vector3(Camera.main.pixelWidth-1, 0, 1));
        // var meshWidth = m_FluidSim.width / (float)Camera.main.pixelWidth * (float)(worldPointE.x - worldPoint.x);
        // var meshHeight = m_FluidSim.height / (float)Camera.main.pixelHeight* (float)(worldPoint.y - worldPointE.y);
        // var matrix = Matrix4x4.TRS(worldPoint, Camera.main.transform.rotation, new Vector3(meshWidth, meshHeight, 1));
        // Graphics.DrawMesh(m_Mesh, matrix, m_FluidMaterialDebugViz,0);
    }
    

    Vector2 m_MousePos = Vector2.zero;
    void OnGUI()
    {
        using var markerScope = s_OnGuiPerfMarker.Auto();
        
        if (Event.current.type == EventType.MouseDown)
        {
            m_MousePos = Event.current.mousePosition;
        }
        
        if (Event.current.type != EventType.Repaint) return;
        
        
        
        DrawMaterial(new Rect(0,0, m_FluidSim.width, m_FluidSim.height), m_FluidMaterialDebugViz);

        if (!m_ShowOverlay) return;
        var thickness = 1;
        var topleft = new Vector2(0, 0);
        var topright = new Vector2(0, m_FluidSim.height);
        var bottomleft = new Vector2(m_FluidSim.width, 0);
        var bottomright = new Vector2(m_FluidSim.width, m_FluidSim.height);
        
        //Graphics.DrawTexture(new Rect(0,0, m_FluidSim.width, m_FluidSim.height), m_Texture);

        DrawLine(topleft,topright, thickness, Color.gray);
        DrawLine(topright, bottomright, thickness, Color.gray);
        DrawLine(bottomright, bottomleft, thickness, Color.gray);
        DrawLine(bottomleft, topleft, thickness, Color.gray);
        
        
        var positions = m_FluidSim.GetPositions();
        for (var i = 0; i < m_FluidSim.m_ParticleCount; i++)
        {
            var position = positions[i];
            DrawCircle(position, m_CircleSize, Color.blue);
        }
        
    }

    static void DrawMaterial(Rect rect, Material mat)
    {
        Graphics.DrawTexture(rect, Texture2D.grayTexture, mat);
    }

    void DrawCircle(Vector2 center, float radius, Color color)
    {
        var rect = new Rect(center.x - radius , center.y - radius, radius+radius, radius+radius);
        Graphics.DrawTexture(rect, m_CircleTexture);
    }

    static void DrawLine(Vector2 start, Vector2 end, int thickness, Color color)
    {
        //var rect = new Rect(0,0, thickness, thickness);
        Handles.DrawLine(start, end, thickness);
    }

 
}
