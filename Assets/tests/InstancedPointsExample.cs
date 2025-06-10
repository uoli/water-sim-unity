using System;
using System.Runtime.InteropServices;
using UnityEngine;
using UnityEngine.Rendering;
using Random = UnityEngine.Random;

public class InstancedPointsExample : MonoBehaviour
{
    public Material instancedMaterial;
    public int instanceCount = 1000;
    public float pointSize = 0.1f;

    Mesh quadMesh;
    Matrix4x4[] matrices;
    Vector4[] positions;
    Vector4[] colors;
    MaterialPropertyBlock mpb;
    uint[] m_PointRenderArgsData;
    ComputeBuffer m_PointRenderArgsBuffer;
    ComputeBuffer positionBuffer;
    RenderTexture m_PointRenderTexture;

    void Start()
    {
        
        m_PointRenderTexture = new RenderTexture(100, 100, 0, RenderTextureFormat.ARGB32);
        m_PointRenderTexture.enableRandomWrite = true;
        m_PointRenderTexture.Create();
        
        // Create a simple quad mesh centered at (0,0,0)
        quadMesh = new Mesh();
        quadMesh.SetVertices(new Vector3[] {
            new Vector3(-0.5f, 0.5f, 0),
            new Vector3(0.5f, 0.5f, 0),
            new Vector3(0.5f, -0.5f, 0),
            new Vector3(-0.5f, -0.5f, 0)
        });
        quadMesh.SetIndices(new int[] { 0,1,2, 2,3,0 }, MeshTopology.Triangles, 0);
        quadMesh.RecalculateBounds();

        // Setup arrays
        matrices = new Matrix4x4[instanceCount];
        positions = new Vector4[instanceCount];
        colors = new Vector4[instanceCount];
        positionBuffer = new ComputeBuffer(instanceCount,  Marshal.SizeOf(typeof(Vector4)));
        positionBuffer.name = "DrawingPositionBuffer";

        // Fill with data
        for (int i = 0; i < instanceCount; i++)
        {
            float x = Random.Range(-5f, 5f);
            float y = Random.Range(-5f, 5f);
            positions[i] = new Vector4(x, y, 0, 0);

            colors[i] = new Vector4(Random.value, Random.value, Random.value, 1);

            // Use identity matrices because position is handled in shader
            matrices[i] = Matrix4x4.identity;
        }
        positionBuffer.SetData(positions);
        // Enable instancing on the material
        instancedMaterial.enableInstancing = true;

        mpb = new MaterialPropertyBlock();
        mpb.SetVectorArray("_Position", positions);
        mpb.SetVectorArray("_Color", colors);
        
        m_PointRenderArgsData = new uint[]
        {
            quadMesh.GetIndexCount(0),
            (uint)instanceCount,
            quadMesh.GetIndexStart(0),
            quadMesh.GetBaseVertex(0),
            0
        };
        m_PointRenderArgsBuffer = new ComputeBuffer(1,m_PointRenderArgsData.Length* Marshal.SizeOf(typeof(uint)), ComputeBufferType.IndirectArguments);
        
        instancedMaterial.SetBuffer("positions", positionBuffer);
    }

    void OnGUI()
    {
        GUI.DrawTexture(new Rect(100,100,100,100), m_PointRenderTexture);
    }

    void Update()
    {
        m_PointRenderArgsBuffer.SetData(m_PointRenderArgsData);
        var bounds = new Bounds(Vector3.zero, 10000*Vector3.one); // use tighter bounds
        Graphics.DrawMeshInstancedIndirect(quadMesh, 0, instancedMaterial, bounds, m_PointRenderArgsBuffer, 0, mpb );
        
        // CommandBuffer cmd = new CommandBuffer();
        // cmd.SetRenderTarget(m_PointRenderTexture);
        // cmd.ClearRenderTarget(false, true, Color.blue);
        // cmd.DrawMeshInstancedIndirect(quadMesh, 0, instancedMaterial, 0, m_PointRenderArgsBuffer );
        // Graphics.ExecuteCommandBuffer(cmd);
        // cmd.Release();

    }
}
