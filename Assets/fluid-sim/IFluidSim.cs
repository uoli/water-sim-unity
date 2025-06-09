
using System;
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
    GridSpatialLookup LookupHelper { get; }
    ReadOnlySpan<Vector2> GetPositions();
    ReadOnlySpan<float> GetDensities();
    ReadOnlySpan<float> GetPressures();
    ReadOnlySpan<Vector2> GetVelocities();

    void Interact(Vector2 mouseInSimulationSpace, float scalingFactor, InteractionDirection interactionDirection);
}
