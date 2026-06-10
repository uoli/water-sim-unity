using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Assertions;

public class FluidRigibodyInteraction : MonoBehaviour
{
    public Rigidbody2D Rigidbody;
    // Sim-space surface length the sample points represent in total (2D
    // "area" is a perimeter); each point gets an equal share.
    public float Area;
    // Conversion between the simulation's arbitrary mass unit and kilograms.
    // Lengths and times are converted exactly; mass is the one free unit.
    public float ForceScale = 1f;

    SurfaceSampler2D m_SurfaceSampler;
    IFluidSim m_FluidSim;
    RectTransform m_SimRect;
    OutputSimulationSurfacePoints[] m_SimResults;
    InputSimulationSurfacePoints[] m_SimInput;

    void OnEnable()
    {
        m_SurfaceSampler = Rigidbody.GetComponent<SurfaceSampler2D>();
        m_FluidSim = GetComponent<IFluidSim>();
        m_SimRect = m_FluidSim.Transform.GetComponent<RectTransform>();
        m_FluidSim.PreSimulation += PreSimulation;
        m_FluidSim.PostSimulation += PostSimulation;

        m_SimInput = new InputSimulationSurfacePoints[m_SurfaceSampler.SurfacePoints.Count];
        m_SimResults = new OutputSimulationSurfacePoints[m_SurfaceSampler.SurfacePoints.Count];
    }

    void OnDisable()
    {
        m_FluidSim.PreSimulation -= PreSimulation;
        m_FluidSim.PostSimulation -= PostSimulation;
    }

    // ---- World <-> simulation space ----
    // The simulation is Width x Height sim units, displayed on the quad whose
    // RectTransform is cached in m_SimRect. Everything crossing the boundary
    // (positions, velocities, forces) must go through these conversions; the
    // sim-units-per-world-unit factor is Width / rect.width, which is far
    // from 1.

    Vector2 WorldToSimPoint(Vector2 worldPoint)
    {
        Vector2 local = m_FluidSim.Transform.InverseTransformPoint(worldPoint);
        var rect = m_SimRect.rect;
        // Deliberately unclamped (unlike Rect.PointToNormalized) so points
        // outside the tank still map continuously.
        return new Vector2(
            (local.x - rect.xMin) / rect.width * m_FluidSim.Width,
            (local.y - rect.yMin) / rect.height * m_FluidSim.Height);
    }

    Vector2 WorldToSimVector(Vector2 worldVector)
    {
        Vector2 local = m_FluidSim.Transform.InverseTransformVector(worldVector);
        var rect = m_SimRect.rect;
        return new Vector2(
            local.x * m_FluidSim.Width / rect.width,
            local.y * m_FluidSim.Height / rect.height);
    }

    Vector2 SimToWorldVector(Vector2 simVector)
    {
        var rect = m_SimRect.rect;
        var local = new Vector2(
            simVector.x * rect.width / m_FluidSim.Width,
            simVector.y * rect.height / m_FluidSim.Height);
        return m_FluidSim.Transform.TransformVector(local);
    }

    Vector2 CalcSurfacePointVelocity(Vector2 worldPoint)
    {
        // Lever arm in world space about the actual center of mass. The sampled
        // point is stored in the body's local space and must not be mixed with
        // world-space positions.
        Vector2 offset = worldPoint - Rigidbody.worldCenterOfMass;
        // Rigidbody2D.angularVelocity is in degrees per second.
        var omega = Rigidbody.angularVelocity * Mathf.Deg2Rad;
        return Rigidbody.linearVelocity + omega * Vector2.Perpendicular(offset);
    }

    void PreSimulation()
    {
        var points = m_SurfaceSampler.SurfacePoints;
        var areaWeight = Area / points.Count;
        for (var i = 0; i < points.Count; i++)
        {
            var p = points[i];
            var worldPoint = (Vector2)Rigidbody.transform.TransformPoint(p.position);

            m_SimInput[i] = new InputSimulationSurfacePoints
            {
                SimSpacePoint = WorldToSimPoint(worldPoint),
                // Re-normalized because a non-square mapping skews directions.
                normal = WorldToSimVector(Rigidbody.transform.TransformDirection(p.normal)).normalized,
                velocity = WorldToSimVector(CalcSurfacePointVelocity(worldPoint)),
                areaWeight = areaWeight
            };
        }
        m_FluidSim.SetRigidBodySurfaceResults(m_SimInput);
    }
    
    void PostSimulation()
    {

        m_FluidSim.RetrieveRigidBodySurfaceResults(m_SimResults);

        for (var index = 0; index < m_SimResults.Length; index++)
        {
            var simResult = m_SimResults[index];
            // The sampled point is in the body's local space; the force must be
            // applied at its world position or the resulting torque is wrong.
            var worldPoint = (Vector2)Rigidbody.transform.TransformPoint(m_SurfaceSampler.SurfacePoints[index].position);
            // The sim reports the momentum exchanged during its step; Impulse
            // mode transfers it exactly, independent of how render frames align
            // with physics steps. No further dt belongs here.
            var worldImpulse = SimToWorldVector(simResult.impulse) * ForceScale;
            Rigidbody.AddForceAtPosition(worldImpulse, worldPoint, ForceMode2D.Impulse);
        }
    }

}
