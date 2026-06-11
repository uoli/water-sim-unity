using System;
using System.Collections.Generic;
using UnityEngine;

[Serializable]
public struct SurfacePoint
{
    public Vector2 position;
    public Vector2 normal;
    // Arc length of the perimeter segment this sample represents, in local
    // units. Zero for corner samples (they reinforce the density barrier and
    // the contact backstop at normal discontinuities but carry no force
    // weight) and for legacy hand-placed points.
    public float edgeLength;
}

public class SurfaceSampler2D : MonoBehaviour
{
    public List<SurfacePoint> SurfacePoints;

    // Distance between samples along the perimeter, in the body's local units.
    // Must stay at or below the fluid's particle spacing (in the same units)
    // or particles can slip between samples.
    public float SampleSpacing = 0.1f;
    // Regenerate from the attached collider on enable, keeping the samples in
    // sync with the collider. Disable to preserve a hand-edited list.
    public bool GenerateOnEnable = true;

    // Corner samples are skipped when adjacent edge normals are this aligned
    // (finely tessellated curves don't need them; real corners do).
    const float k_CornerNormalDotThreshold = 0.99f;

    void OnEnable()
    {
        if (GenerateOnEnable)
            GenerateFromCollider();
    }

    // ---- Layer 1: shape -> closed vertex loops ----
    // Every supported shape reduces to closed polygon outlines; the walker
    // below never knows what produced them. Future custom shapes (sprite
    // outlines via Sprite.GetPhysicsShape, mesh boundary edges) are new
    // adapters here, nothing else changes.
    public void GenerateFromCollider()
    {
        var loops = new List<List<Vector2>>();
        var collider = GetComponent<Collider2D>();
        switch (collider)
        {
            case BoxCollider2D box:
            {
                var half = box.size * 0.5f;
                loops.Add(new List<Vector2>
                {
                    box.offset + new Vector2(-half.x, -half.y),
                    box.offset + new Vector2(half.x, -half.y),
                    box.offset + new Vector2(half.x, half.y),
                    box.offset + new Vector2(-half.x, half.y),
                });
                break;
            }
            case CircleCollider2D circle:
            {
                // Tessellate finely relative to the sample spacing; the walker
                // resamples to even arc length anyway.
                var vertexCount = Mathf.Max(12, Mathf.CeilToInt(2f * Mathf.PI * circle.radius / (SampleSpacing * 0.5f)));
                var loop = new List<Vector2>(vertexCount);
                for (var i = 0; i < vertexCount; i++)
                {
                    var angle = i * 2f * Mathf.PI / vertexCount;
                    loop.Add(circle.offset + circle.radius * new Vector2(Mathf.Cos(angle), Mathf.Sin(angle)));
                }
                loops.Add(loop);
                break;
            }
            case PolygonCollider2D polygon:
            {
                // Arbitrary custom shapes, including multiple paths (holes).
                for (var p = 0; p < polygon.pathCount; p++)
                {
                    var path = polygon.GetPath(p);
                    var loop = new List<Vector2>(path.Length);
                    foreach (var v in path)
                        loop.Add(v + polygon.offset);
                    loops.Add(loop);
                }
                break;
            }
            default:
                Debug.LogWarning($"SurfaceSampler2D: no supported Collider2D found on {name} " +
                                 "(Box, Circle and Polygon are supported); keeping existing samples.", this);
                return;
        }

        SurfacePoints ??= new List<SurfacePoint>();
        SurfacePoints.Clear();
        foreach (var loop in loops)
            SampleLoop(loop, SampleSpacing, SurfacePoints);
    }

    // ---- Layer 2: loop -> samples (shape-agnostic) ----
    static void SampleLoop(List<Vector2> loop, float spacing, List<SurfacePoint> output)
    {
        if (loop.Count < 3 || spacing <= 0f) return;

        // Winding decides which perpendicular points outward: for a
        // counter-clockwise loop (positive signed area) it is (dy, -dx).
        // Holes wind opposite to outer boundaries, so their normals
        // automatically point into the hole, where the fluid is.
        var sign = SignedArea(loop) >= 0f ? 1f : -1f;

        for (var i = 0; i < loop.Count; i++)
        {
            var a = loop[i];
            var b = loop[(i + 1) % loop.Count];
            var edge = b - a;
            var length = edge.magnitude;
            if (length < 1e-6f) continue;
            var dir = edge / length;
            var normal = sign * new Vector2(dir.y, -dir.x);

            // Corner sample at the shared vertex: the normal is discontinuous
            // there and fluid can wedge into the crease between the two edge
            // normals. Bisector normal, zero length (no force weight).
            var prev = loop[(i - 1 + loop.Count) % loop.Count];
            var prevEdge = a - prev;
            if (prevEdge.sqrMagnitude > 1e-12f)
            {
                var prevDir = prevEdge.normalized;
                var prevNormal = sign * new Vector2(prevDir.y, -prevDir.x);
                var bisector = prevNormal + normal;
                if (Vector2.Dot(prevNormal, normal) < k_CornerNormalDotThreshold && bisector.sqrMagnitude > 1e-6f)
                    output.Add(new SurfacePoint { position = a, normal = bisector.normalized, edgeLength = 0f });
            }

            // Evenly spaced samples along the edge (midpoint rule): each
            // carries exactly the arc length it represents, so the boundary
            // integral weights are measured, not hand-tuned.
            var count = Mathf.Max(1, Mathf.RoundToInt(length / spacing));
            var ds = length / count;
            for (var s = 0; s < count; s++)
            {
                var t = (s + 0.5f) * ds;
                output.Add(new SurfacePoint { position = a + dir * t, normal = normal, edgeLength = ds });
            }
        }
    }

    static float SignedArea(List<Vector2> loop)
    {
        var area = 0f;
        for (var i = 0; i < loop.Count; i++)
        {
            var a = loop[i];
            var b = loop[(i + 1) % loop.Count];
            area += a.x * b.y - b.x * a.y;
        }
        return area * 0.5f;
    }

    void OnDrawGizmosSelected()
    {
        if (SurfacePoints == null) return;
        var normalLength = SampleSpacing * 1.5f;
        foreach (var p in SurfacePoints)
        {
            var world = transform.TransformPoint(p.position);
            var worldNormal = transform.TransformDirection(p.normal);
            Gizmos.color = p.edgeLength > 0f ? Color.cyan : Color.yellow;
            Gizmos.DrawSphere(world, SampleSpacing * 0.15f);
            Gizmos.DrawLine(world, world + worldNormal * normalLength);
        }
    }
}
