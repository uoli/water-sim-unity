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
    // Akinci-style boundary pseudo-mass: the patch of body volume this sample
    // represents, filled at the fluid's rest density. Contributes to fluid
    // density so the body excludes volume.
    public float pseudoMass;
}
public struct OutputSimulationSurfacePoints
{
    // Momentum exchanged with the fluid during the sim step (the surface force
    // already integrated over the step's dt) — apply with ForceMode2D.Impulse.
    public Vector2 impulse;
    // Shepard weight sum of the interpolation that produced the impulse. Doubles
    // as a smooth submersion measure (~0 in air, ~1 fully in fluid): the
    // analytic coupling uses it as per-point coverage, and the back-reaction
    // uses it to return each particle its exact share.
    public float weightSum;
    // Shepard-interpolated fluid velocity at the point (sim units), for drag.
    public Vector2 fluidVelocity;
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
    // Duration of the most recent simulation substep, so impulse consumers can
    // convert accumulated momentum into a force rate over actual sim time.
    float LastStepDeltaTime { get; }
    // Sim-space gravity magnitude, needed by analytic (hydrostatic) coupling.
    float Gravity { get; }
}
