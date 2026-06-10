using UnityEditor;
using UnityEngine;

[CustomEditor(typeof(FluidSimGPU))]
public class FluidSimGPUEditor : Editor
{
    public override void OnInspectorGUI()
    {
        var sim = target as FluidSimGPU;

        base.OnInspectorGUI();

        EditorGUILayout.Space();
        EditorGUILayout.LabelField("Derived values", EditorStyles.boldLabel);
        EditorGUILayout.LabelField($"Particle spacing: {sim.ParticleSpacing:F3}");
        EditorGUILayout.LabelField($"Smoothing radius: {sim.SmoothingRadius:F3}");
        EditorGUILayout.LabelField($"Particle mass: {sim.Mass:F4}");
        EditorGUILayout.LabelField($"Effective DeltaTime: {sim.EffectiveDeltaTime:E2}");

        var stepsForRealtime = 1f / 60f / Mathf.Max(1e-9f, sim.EffectiveDeltaTime);
        EditorGUILayout.LabelField($"Substeps/frame for realtime at 60fps: {stepsForRealtime:F1} (cap: {sim.MaxStepsPerFrame})");
        if (stepsForRealtime > sim.MaxStepsPerFrame)
        {
            EditorGUILayout.HelpBox(
                "MaxStepsPerFrame is below what realtime needs at this substep size; " +
                "the simulation will run in slow motion. Raise MaxStepsPerFrame or " +
                "lower PressureMultiplier/Gravity.",
                MessageType.Warning);
        }
    }
}
