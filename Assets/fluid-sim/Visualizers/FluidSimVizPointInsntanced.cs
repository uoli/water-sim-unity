using UnityEngine;
struct DrawSimPointsInstanced
{
    Material m_PointInstancedMaterial;
    Mesh m_PointMesh;
    ComputeBuffer m_PointRenderArgsBuffer;
    uint[] m_PointRenderArgsData;
    public void InitializePointInstancedRender(int particleCount, Material material)
    {
        if (m_PointInstancedMaterial)
            MonoBehaviour.Destroy(m_PointInstancedMaterial);
        m_PointInstancedMaterial = new Material(material);
        m_PointInstancedMaterial.enableInstancing = true;
        m_PointMesh = new Mesh();
        m_PointMesh.SetVertices(new Vector3[] { new(-0.5f, 0.5f, 0), new(0.5f, 0.5f, 0), new(0.5f, -0.5f, 0), new(-0.5f, -0.5f, 0) });
        m_PointMesh.SetIndices(new[] { 0, 1, 2, 2, 3, 0 }, MeshTopology.Triangles, 0);
        m_PointMesh.Optimize();
        m_PointRenderArgsData = new uint[]
        {
            m_PointMesh.GetIndexCount(0),
            (uint)particleCount,
            m_PointMesh.GetIndexStart(0),
            m_PointMesh.GetBaseVertex(0),
            0
        };
        m_PointRenderArgsBuffer = new ComputeBuffer(1, m_PointRenderArgsData.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        m_PointRenderArgsBuffer.name = "m_PointRenderArgsBuffer";
    }
    public void CleanupPointInstancedRender()
    {
        m_PointRenderArgsBuffer.Dispose();
        MonoBehaviour.Destroy(m_PointMesh);
    }
    public void RenderPointsGPU(float maxVelocity, int particleCount,
        ComputeBuffer positions, ComputeBuffer velocities, ComputeBuffer densities,
        ComputeBuffer pressure, int particleVisualizationMode,
        Vector3 world_origin, float display_area_width, float display_area_height,
        float sim_width, float sim_height, float circle_size, float scaling_factor,
        Color circle_color
        )
    {

        m_PointRenderArgsData[1] = (uint)particleCount;
        m_PointRenderArgsBuffer.SetData(m_PointRenderArgsData);

        m_PointInstancedMaterial.SetColor("_Color", circle_color);
        m_PointInstancedMaterial.SetBuffer("positions", positions);
        m_PointInstancedMaterial.SetBuffer("_particle_densities", densities);
        m_PointInstancedMaterial.SetBuffer("_particle_pressures", pressure);
        m_PointInstancedMaterial.SetBuffer("_particle_velocities", velocities);
        m_PointInstancedMaterial.SetFloat("_max_velocity", maxVelocity);
        m_PointInstancedMaterial.SetInt("_particleVisMode", (int)particleVisualizationMode);

        m_PointInstancedMaterial.SetVector("world_origin", world_origin);
        m_PointInstancedMaterial.SetFloat("display_area_width", display_area_width);
        m_PointInstancedMaterial.SetFloat("display_area_height", display_area_height);

        m_PointInstancedMaterial.SetInt("point_count", particleCount);
        m_PointInstancedMaterial.SetFloat("sim_width", sim_width);
        m_PointInstancedMaterial.SetFloat("sim_height", sim_height);
        m_PointInstancedMaterial.SetFloat("circle_size", circle_size);
        m_PointInstancedMaterial.SetFloat("scaling_factor", scaling_factor);

        var bounds = new Bounds(Vector3.zero, 10000 * Vector3.one); // use tighter bounds
        Graphics.DrawMeshInstancedIndirect(m_PointMesh, 0, m_PointInstancedMaterial, bounds, m_PointRenderArgsBuffer);

        // CommandBuffer cmd = new CommandBuffer();
        // cmd.SetRenderTarget(m_PointRenderTexture);
        // //cmd.ClearRenderTarget(true, true, Color.red);
        // Graphics.ExecuteCommandBuffer(cmd);
        // cmd.Release();
    }
}
