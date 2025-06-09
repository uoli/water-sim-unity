using System;
using System.Runtime.InteropServices;
using Unity.Collections;
using UnityEngine;

public struct NativeArrayToComputeAdapter<T> : IDisposable where T : struct
{
    ComputeBuffer m_ComputeBuffer;
    public NativeArrayToComputeAdapter(NativeArray<T> array )
    {
        m_ComputeBuffer = new ComputeBuffer(array.Length, Marshal.SizeOf(typeof(T)));
        m_ComputeBuffer.SetData(array);
    }

    public ComputeBuffer Buffer => m_ComputeBuffer;

    public void Dispose()
    {
        m_ComputeBuffer?.Dispose();
    }
}

