
using System;
using Unity.Collections;
using UnityEngine;

public enum InteractionDirection
{
    Attract,
    Repel
}

public interface IFluidSim
{
    int ParticleCount { get; }
    float Mass { get; }
    int Height { get; }
    int Width { get; }
    float SmoothingRadius { get; }
    float TargetDensity { get; }
    
    bool HasDataInCompute { get; }
    
    GridSpatialLookup LookupHelper { get; }
    ComputeBuffer GetPositionComputeBuffer();
    ComputeBuffer GetDensities();
    ComputeBuffer GetPressures();
    ComputeBuffer GetVelocities();

    void Interact(Vector2 mouseInSimulationSpace, float scalingFactor, InteractionDirection interactionDirection);
}
