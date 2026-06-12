using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Assertions;

public enum FluidCouplingMode
{
    // Body force = analytic hydrostatic buoyancy + drag, computed from smooth
    // measured quantities (per-point submersion and local flow velocity).
    // Noiseless by construction and directly tunable; momentum with the fluid
    // agrees on average rather than per pair.
    AnalyticBuoyancy,
    // Body force = the sim's sampled surface impulses. Exact pairwise
    // Newton's third law, but inherits WCSPH pressure noise.
    SampledImpulses,
}

public class FluidRigibodyInteraction : MonoBehaviour
{
    public Rigidbody2D Rigidbody;
    public FluidCouplingMode CouplingMode = FluidCouplingMode.AnalyticBuoyancy;
    // Drag strength against relative flow, used by the analytic mode.
    public float DragCoefficient = 1f;
    // Sim-space surface length the sample points represent in total (2D
    // "area" is a perimeter); each point gets an equal share.
    public float Area;
    // Conversion between the simulation's arbitrary mass unit and kilograms.
    // Lengths and times are converted exactly; mass is the one free unit.
    public float ForceScale = 1f;
    // Time constant (seconds) of the exponential smoothing applied to the
    // impulses before they reach the body. The high-frequency content of the
    // sampled fluid force is discretization noise — particles entering and
    // leaving kernels, WCSPH pressure fluctuations — not fluid dynamics; a
    // real hull's inertia would never respond to it. 0 disables smoothing.
    public float ImpulseSmoothingTime = 0.08f;

    SurfaceSampler2D m_SurfaceSampler;
    IFluidSim m_FluidSim;
    RectTransform m_SimRect;
    OutputSimulationSurfacePoints[] m_SimResults;
    InputSimulationSurfacePoints[] m_SimInput;
    Vector2[] m_AccumulatedImpulses;
    Vector2[] m_SmoothedForces;
    float[] m_LatestCoverage;
    Vector2[] m_LatestFluidVelocity;
    float[] m_SmoothedCoverage;
    Vector2[] m_SmoothedFluidVelocity;
    bool m_HasFreshSample;
    float m_AccumulatedSimTime;
    float[] m_PseudoMasses;
    Vector4 m_PseudoMassCacheKey = new Vector4(float.NaN, 0, 0, 0);

    void OnEnable()
    {
        m_SurfaceSampler = Rigidbody.GetComponent<SurfaceSampler2D>();
        m_FluidSim = GetComponent<IFluidSim>();
        m_SimRect = m_FluidSim.Transform.GetComponent<RectTransform>();
        m_FluidSim.PreSimulation += PreSimulation;
        m_FluidSim.PostSimulation += PostSimulation;
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

    Vector2 SimToWorldPoint(Vector2 simPoint)
    {
        var rect = m_SimRect.rect;
        var local = new Vector2(
            simPoint.x / m_FluidSim.Width * rect.width + rect.xMin,
            simPoint.y / m_FluidSim.Height * rect.height + rect.yMin);
        return m_FluidSim.Transform.TransformPoint(local);
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

        // Sized lazily: the sampler may regenerate its points (OnEnable order
        // between components is undefined, and the editor button can run at
        // any time), so the count is only trusted at use.
        if (m_SimInput == null || m_SimInput.Length != points.Count)
        {
            m_SimInput = new InputSimulationSurfacePoints[points.Count];
            m_SimResults = new OutputSimulationSurfacePoints[points.Count];
        }

        // Generated samples carry the exact arc length they represent; convert
        // it to sim units along the edge tangent. Legacy hand-placed lists
        // (edgeLength never set) fall back to splitting the manual Area evenly.
        var hasMeasuredLengths = false;
        for (var i = 0; i < points.Count; i++)
        {
            if (points[i].edgeLength > 0f)
            {
                hasMeasuredLengths = true;
                break;
            }
        }
        var uniformWeight = Area / points.Count;

        for (var i = 0; i < points.Count; i++)
        {
            var p = points[i];
            var worldPoint = (Vector2)Rigidbody.transform.TransformPoint(p.position);

            float areaWeight;
            if (hasMeasuredLengths)
            {
                var localTangent = Vector2.Perpendicular(p.normal) * p.edgeLength;
                var worldTangent = (Vector2)Rigidbody.transform.TransformVector(localTangent);
                areaWeight = WorldToSimVector(worldTangent).magnitude;
            }
            else
            {
                areaWeight = uniformWeight;
            }

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
        if (m_PseudoMasses != null && m_PseudoMasses.Length == m_SimInput.Length && key == m_PseudoMassCacheKey) return;
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

        if (m_AccumulatedImpulses == null || m_AccumulatedImpulses.Length != m_SimResults.Length)
        {
            m_AccumulatedImpulses = new Vector2[m_SimResults.Length];
            m_SmoothedForces = new Vector2[m_SimResults.Length];
            m_LatestCoverage = new float[m_SimResults.Length];
            m_LatestFluidVelocity = new Vector2[m_SimResults.Length];
            m_SmoothedCoverage = new float[m_SimResults.Length];
            m_SmoothedFluidVelocity = new Vector2[m_SimResults.Length];
        }

        // Only accumulate here. The sim's substeps fire irregularly relative
        // to the physics clock (some render frames run zero substeps), so
        // applying impulses at this point delivers buoyancy in packets with
        // gaps — the body falls between packets and gets kicked back ("jumpy").
        // FixedUpdate converts the accumulated momentum into a continuous
        // force where the physics engine actually integrates.
        for (var index = 0; index < m_SimResults.Length; index++)
        {
            var impulse = m_SimResults[index].impulse;
            // Never let non-finite values into the accumulator: if the fluid
            // blows up (NaN/Infinity), drop the sample and say why.
            if (!float.IsFinite(impulse.x) || !float.IsFinite(impulse.y))
            {
                Debug.LogWarning($"Non-finite fluid impulse at surface point {index}; the fluid simulation has likely gone unstable. Skipping.", this);
                continue;
            }
            m_AccumulatedImpulses[index] += impulse;

            // Measured states for the analytic mode: per-point submersion
            // (the Shepard weight sum, a smooth 0..1 ramp) and local flow.
            var coverage = m_SimResults[index].weightSum;
            var fluidVelocity = m_SimResults[index].fluidVelocity;
            if (float.IsFinite(coverage) && float.IsFinite(fluidVelocity.x) && float.IsFinite(fluidVelocity.y))
            {
                m_LatestCoverage[index] = Mathf.Clamp01(coverage);
                m_LatestFluidVelocity[index] = fluidVelocity;
            }
        }
        m_HasFreshSample = true;
        m_AccumulatedSimTime += m_FluidSim.LastStepDeltaTime;
    }

    void FixedUpdate()
    {
        if (m_AccumulatedImpulses == null) return;
        if (CouplingMode == FluidCouplingMode.AnalyticBuoyancy)
            ApplyAnalyticForces();
        else
            ApplySampledImpulses();
    }

    // Body force from smooth measurements instead of sampled impulses: the
    // fluid only MEASURES per-point submersion and local flow velocity (both
    // Shepard interpolations, smooth by construction); the forces are the
    // textbook formulas. Buoyancy via the divergence theorem over the hull:
    //   A_submerged = sum cov_i * x_i * n.x_i * ds_i
    // — the waterline chord that closes the submerged region's boundary has a
    // vertical normal (n.x = 0) and contributes nothing, so the sum is exact
    // Archimedes in the sampling limit. The centroid uses the matching flux
    // forms; drag opposes relative flow per point. The fluid still receives
    // the sim's sampled reaction, so momentum agrees on average, not per pair.
    void ApplyAnalyticForces()
    {
        if (m_SimInput == null || m_SmoothedCoverage == null) return;
        var count = Mathf.Min(m_SimInput.Length, m_SmoothedCoverage.Length);

        var alpha = ImpulseSmoothingTime > 0f
            ? 1f - Mathf.Exp(-Time.fixedDeltaTime / ImpulseSmoothingTime)
            : 1f;
        if (m_HasFreshSample)
        {
            for (var i = 0; i < count; i++)
            {
                m_SmoothedCoverage[i] = Mathf.Lerp(m_SmoothedCoverage[i], m_LatestCoverage[i], alpha);
                m_SmoothedFluidVelocity[i] = Vector2.Lerp(m_SmoothedFluidVelocity[i], m_LatestFluidVelocity[i], alpha);
            }
            m_HasFreshSample = false;
        }

        // Body-centered coordinates keep the flux sums numerically tame.
        var center = WorldToSimPoint(Rigidbody.worldCenterOfMass);
        var area = 0f;
        var momentX = 0f;
        var momentY = 0f;
        for (var i = 0; i < count; i++)
        {
            var sp = m_SimInput[i];
            var x = sp.SimSpacePoint.x - center.x;
            var y = sp.SimSpacePoint.y - center.y;
            var flux = m_SmoothedCoverage[i] * sp.normal.x * sp.areaWeight;
            area += x * flux;
            momentX += 0.5f * x * x * flux;
            momentY += x * y * flux;
        }

        if (area > 1e-4f)
        {
            var buoyancySim = new Vector2(0f, m_FluidSim.TargetDensity * m_FluidSim.Gravity * area);
            var centroidSim = center + new Vector2(momentX / area, momentY / area);
            var worldForce = SimToWorldVector(buoyancySim) * ForceScale;
            Rigidbody.AddForceAtPosition(worldForce, SimToWorldPoint(centroidSim), ForceMode2D.Force);
        }

        var points = m_SurfaceSampler.SurfacePoints;
        var pointCount = Mathf.Min(count, points.Count);
        for (var i = 0; i < pointCount; i++)
        {
            var coverage = m_SmoothedCoverage[i];
            if (coverage <= 1e-4f) continue;
            var sp = m_SimInput[i];
            var dragSim = (m_SmoothedFluidVelocity[i] - sp.velocity) * (DragCoefficient * coverage * sp.areaWeight);
            var worldPoint = (Vector2)Rigidbody.transform.TransformPoint(points[i].position);
            Rigidbody.AddForceAtPosition(SimToWorldVector(dragSim) * ForceScale, worldPoint, ForceMode2D.Force);
        }
    }

    void ApplySampledImpulses()
    {
        var points = m_SurfaceSampler.SurfacePoints;
        var count = Mathf.Min(m_AccumulatedImpulses.Length, points.Count);

        // A physics step with no sim substeps since the last one carries no
        // information about the fluid force — it does NOT mean the force is
        // zero. Updating the average with a zero rate on those steps sags the
        // support by alpha each time (a ~10Hz droop train at typical settings,
        // felt as jumpiness); instead, hold the previous estimate and keep
        // applying it.
        var hasNewData = m_AccumulatedSimTime > 0f;

        // 1 - exp(-dt/tau): exponential moving average with a framerate-
        // independent time constant, evaluated at the fixed-step cadence.
        // Smoothed per point so the torque distribution over the hull is
        // preserved.
        var alpha = ImpulseSmoothingTime > 0f
            ? 1f - Mathf.Exp(-Time.fixedDeltaTime / ImpulseSmoothingTime)
            : 1f;

        for (var index = 0; index < count; index++)
        {
            if (hasNewData)
            {
                // Accumulated momentum divided by the sim time it was exchanged
                // over: a force rate that is independent of how many substeps
                // happened to land between physics steps.
                var forceRate = m_AccumulatedImpulses[index] / m_AccumulatedSimTime;
                m_AccumulatedImpulses[index] = Vector2.zero;
                m_SmoothedForces[index] = Vector2.Lerp(m_SmoothedForces[index], forceRate, alpha);
            }
            if (m_SmoothedForces[index] == Vector2.zero) continue;

            // The sampled point is in the body's local space; the force must be
            // applied at its world position or the resulting torque is wrong.
            var worldPoint = (Vector2)Rigidbody.transform.TransformPoint(points[index].position);
            var worldForce = SimToWorldVector(m_SmoothedForces[index]) * ForceScale;
            Rigidbody.AddForceAtPosition(worldForce, worldPoint, ForceMode2D.Force);
        }

        if (hasNewData)
            m_AccumulatedSimTime = 0f;
    }

}
