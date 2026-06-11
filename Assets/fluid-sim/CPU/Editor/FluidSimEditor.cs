using UnityEditor;
using UnityEngine;

[CustomEditor(typeof(FluidSim))]
public class FluidSimEditor : Editor
{
    public override void OnInspectorGUI()
    {
        var sim = target as FluidSim;

        base.OnInspectorGUI();

        if (GUILayout.Button("Step Simulation"))
        {
            sim.DoUpdate();
        }

        EditorGUILayout.Space();
        EditorGUILayout.LabelField("Derived values", EditorStyles.boldLabel);
        EditorGUILayout.LabelField($"Particle spacing: {sim.ParticleSpacing:F3}");
        EditorGUILayout.LabelField($"Smoothing length: {sim.SmoothingLength:F3}");
        EditorGUILayout.LabelField($"Particle mass: {sim.Mass:F4}");

        // For the linear EOS the speed of sound is sqrt(PressureMultiplier).
        // Keeping it ~10x the free-fall speed over the box height keeps
        // compression around 1%; the CFL step then bounds the usable substep.
        var freeFallSpeed = Mathf.Sqrt(2f * sim.Gravity * sim.height);
        var suggestedStiffness = 100f * freeFallSpeed * freeFallSpeed;
        var speedOfSound = Mathf.Sqrt(Mathf.Max(1e-6f, sim.m_PressureMultiplier));
        var cflStep = 0.3f * sim.SmoothingLength / speedOfSound;
        EditorGUILayout.Space();
        EditorGUILayout.LabelField("Tuning hints", EditorStyles.boldLabel);
        EditorGUILayout.LabelField($"Suggested PressureMultiplier: ~{suggestedStiffness:F0}");
        EditorGUILayout.LabelField($"Max stable SimulationStep at current stiffness: ~{cflStep:E2}");
    }
}
