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
    float[] m_PseudoMasses;
    Vector4 m_PseudoMassCacheKey = new Vector4(float.NaN, 0, 0, 0);

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

        UpdatePseudoMassesIfNeeded();
        for (var i = 0; i < m_SimInput.Length; i++)
            m_SimInput[i].pseudoMass = m_PseudoMasses[i];

        m_FluidSim.SetRigidBodySurfaceResults(m_SimInput);
    }

    // Psi depends only on sim-space distances between samples and the kernel
    // radius. Under rigid motion those distances are invariant as long as the
    // world->sim scale and the body's scale stay fixed, so the computation is
    // cached against exactly those invariants. The one case rigid motion does
    // change them is an anisotropic mapping (different sim-units-per-world in
    // x and y), where rotating the body skews distances: the key includes the
    // rotation then, recomputing whenever the body turns.
    void UpdatePseudoMassesIfNeeded()
    {
        var rect = m_SimRect.rect;
        var scaleX = m_FluidSim.Width / rect.width;
        var scaleY = m_FluidSim.Height / rect.height;
        var anisotropic = !Mathf.Approximately(scaleX, scaleY);
        var key = new Vector4(
            m_FluidSim.SmoothingRadius,
            scaleX * Rigidbody.transform.lossyScale.x,
            scaleY * Rigidbody.transform.lossyScale.y,
            anisotropic ? Rigidbody.rotation : 0f);
        if (key == m_PseudoMassCacheKey) return;
        m_PseudoMassCacheKey = key;

        // Self-calibrating boundary pseudo-mass (Akinci et al. 2012):
        // psi_b = rho0 / sum_k W(|x_b - x_k|), summed over the boundary
        // samples themselves (including b). The denominator measures local
        // sample crowding; its reciprocal is the patch of body volume the
        // sample represents, so psi is that patch filled at rest density.
        // Sampling the surface more densely makes each psi proportionally
        // smaller — the wall's effect does not depend on the sampling rate.
        if (m_PseudoMasses == null || m_PseudoMasses.Length != m_SimInput.Length)
            m_PseudoMasses = new float[m_SimInput.Length];
        var radius = m_FluidSim.SmoothingRadius;
        var sqrRadius = radius * radius;
        var kernelFactor = SmoothingKernels.CalcSmoothingKernel2Factor(radius);
        for (var i = 0; i < m_SimInput.Length; i++)
        {
            var selfPosition = m_SimInput[i].SimSpacePoint;
            var crowding = 0f;
            for (var k = 0; k < m_SimInput.Length; k++)
            {
                var sqrDst = (m_SimInput[k].SimSpacePoint - selfPosition).sqrMagnitude;
                if (sqrDst > sqrRadius) continue;
                crowding += SmoothingKernels.SmoothingKernel2(Mathf.Sqrt(sqrDst), radius, kernelFactor);
            }
            // crowding >= W(0) because the sample counts itself; never zero.
            m_PseudoMasses[i] = m_FluidSim.TargetDensity / crowding;
        }
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
