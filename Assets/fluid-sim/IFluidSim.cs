using System;
using System.Collections.Generic;
using UnityEngine;

public enum InteractionDirection
{
    Attract,
    Repel
}

public struct InputSimulationSurfacePoints
{
    public Vector2 SimSpacePoint;
    public Vector2 normal;
    public Vector2 velocity;
    public float areaWeight;
}
public struct OutputSimulationSurfacePoints
{
    // Momentum exchanged with the fluid during the sim step (the surface force
    // already integrated over the step's dt) — apply with ForceMode2D.Impulse.
    public Vector2 impulse;
}


public interface IFluidSim
{
    int ParticleCount { get; }
    float Mass { get; }
    int Height { get; }
    int Width { get; }
    float SmoothingRadius { get; }
    float TargetDensity { get; }
    
    GridSpatialLookup LookupHelper { get; }
    ComputeBuffer GetPositionComputeBuffer();
    ComputeBuffer GetDensities();
    ComputeBuffer GetPressures();
    ComputeBuffer GetVelocities();

    void Interact(Vector2 mouseInSimulationSpace, float scalingFactor, InteractionDirection interactionDirection);
    
    event Action PreSimulation;
    event Action PostSimulation;
    void SetRigidBodySurfaceResults(IList<InputSimulationSurfacePoints> points);
    void RetrieveRigidBodySurfaceResults(IList<OutputSimulationSurfacePoints> points);
    ComputeBuffer InputExternalPoints { get; }
    
    Transform Transform { get; }
}
