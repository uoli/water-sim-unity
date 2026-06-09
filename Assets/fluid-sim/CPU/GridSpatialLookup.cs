using System;
using Unity.Collections;
using Unity.Jobs;
using UnityEngine;

public struct GridSpatialLookup : IDisposable
{
    float m_CellSize;
    int m_ParticleCount;
    bool m_Valid;
    
    NativeArray<Entry> m_SpatialLookup;
    NativeArray<int> m_StartIndices;

    NativeArray<Vector2Int> m_LookupKernel;
    

    struct Entry : IComparable<Entry>
    {
        internal int index;
        internal int cellKey;
        // The actual cell, kept because distinct cells can hash to the same key.
        internal Vector2Int cell;

        public int CompareTo(Entry other)
        {
            return cellKey.CompareTo(other.cellKey);
        }
    }
    
    public bool IsValid => m_Valid;
    public float CellSize => m_CellSize;
    
    public GridSpatialLookup(float cellSize, int particleCount)
    {
        m_CellSize = cellSize;
        m_ParticleCount = particleCount;
        m_Valid = true;
        
        m_SpatialLookup = new NativeArray<Entry>(particleCount, Allocator.Persistent);
        m_StartIndices = new NativeArray<int>(particleCount, Allocator.Persistent);
        
        var lookupKernel = new Vector2Int[]{
            new (-1,-1), new (-1,0), new (-1,1),
            new (0,-1),  new (0,0),  new (0,1),
            new (1,-1),  new (1,0),  new (1,1),
        };
        m_LookupKernel = new NativeArray<Vector2Int>(lookupKernel, Allocator.Persistent);
    }

    Vector2Int PositionToCell(Vector2 position)
    {
        return new Vector2Int(
            Mathf.FloorToInt(position.x / m_CellSize),
            Mathf.FloorToInt(position.y / m_CellSize)
        );
    }

    int CellPositionToKey(Vector2Int cell)
    {
        var hash = (cell.x * 73856093) ^ (cell.y * 19349663);
        var cellKey = ((hash % m_ParticleCount) + m_ParticleCount) % m_ParticleCount;
        return cellKey;
    }

    int PositionToCellKey(Vector2 position)
    {
        var cell = PositionToCell(position);
        var cellKey = CellPositionToKey(cell);
        return cellKey;
    }

    public void UpdateParticles(ReadOnlySpan<Vector2> positions)
    {
        for (var i = 0; i < positions.Length; i++)
        {
            var position = positions[i];
            var cell = PositionToCell(position);
            var cellKey = CellPositionToKey(cell);
            m_SpatialLookup[i] = new Entry() { index = i, cellKey = cellKey, cell = cell };
        }
        
        //Sort array by the cellKey
        JobHandle dependencyJobHandle = default;
        var sortJob = m_SpatialLookup.SortJob();
        var sortJobHandle = sortJob.Schedule(dependencyJobHandle);
        sortJobHandle.Complete();

        //Keep track where each cell starts on the m_SpatialLookup array
        var entryIndex = 0;
        for (var i = 0; i < m_StartIndices.Length; i++)
        {
            if (entryIndex < m_StartIndices.Length && m_SpatialLookup[entryIndex].cellKey == i)
            {
                m_StartIndices[i] = entryIndex;
                do
                {
                    entryIndex++;
                }
                while(entryIndex < m_StartIndices.Length && m_SpatialLookup[entryIndex].cellKey == i );
            }
            else
            {
                m_StartIndices[i] = -1;
            }

        }
    }

    public void GetParticlesAround(Vector2 position, NativeList<int> particleIndexes)
    {
        var centralCell = PositionToCell(position);
        Span<int> visitedKeys = stackalloc int[9];
        var visitedCount = 0;
        foreach (var cellOffset in m_LookupKernel)
        {
            var cell = centralCell + cellOffset;
            var cellKey = CellPositionToKey(cell);

            // Two of the nine neighbor cells can hash to the same key; iterating
            // that bucket twice would double-count its particles.
            var alreadyVisited = false;
            for (var v = 0; v < visitedCount; v++)
            {
                if (visitedKeys[v] == cellKey)
                    alreadyVisited = true;
            }
            if (alreadyVisited) continue;
            visitedKeys[visitedCount++] = cellKey;

            //do the lookup
            var lookupStart = m_StartIndices[cellKey];
            if (lookupStart == -1) continue;
            for (var i = lookupStart; i < m_SpatialLookup.Length; i++)
            {
                var entry = m_SpatialLookup[i];
                if (entry.cellKey != cellKey) break;
                // The bucket can also hold entries from far-away cells that hash
                // to the same key; only accept the 3x3 neighborhood. Chebyshev
                // distance <= 1 also picks up the particles of a colliding
                // neighbor cell on this single visit.
                var cellDif = entry.cell - centralCell;
                if (Mathf.Abs(cellDif.x) > 1 || Mathf.Abs(cellDif.y) > 1) continue;
                particleIndexes.Add(entry.index);
            }
        }
    }

    public void Dispose()
    {
        m_SpatialLookup.Dispose();
        m_StartIndices.Dispose();
        m_LookupKernel.Dispose();
        m_Valid = false;
    }
}
