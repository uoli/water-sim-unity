using UnityEditor;
using UnityEngine;

[CustomEditor(typeof(FluidSim))]
public class FluidSimEditor : Editor
{
    public override void OnInspectorGUI()
    {
        var sim = target as FluidSim;

        base.OnInspectorGUI();
        if (GUILayout.Button("Calculate ideal smooth length"))
        {
            sim.SmoothingLength = 1.3f * Mathf.Sqrt(sim.width * sim.height/ (float) sim.m_ParticleCount);
        }

        if (GUILayout.Button("Step Simulation"))
        {
            sim.DoUpdate();
        }
        EditorGUILayout.LabelField($"TargetDensity:{sim.TargetDensity()}");
    }
}

