using System;
using System.Collections.Generic;
using UnityEngine;

[Serializable]
public struct SurfacePoint
{
    public Vector2 position;
    public Vector2 normal;
}

public class SurfaceSampler2D : MonoBehaviour
{
    public List<SurfacePoint> SurfacePoints;
    
    //TODO: do some clever decimation here to extract evenly distributed points along the surface
    // For now its an editor tool that allows me to place those manually

}
