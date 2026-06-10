using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Assertions;

public class FluidRigibodyInteraction : MonoBehaviour
{
    public Rigidbody2D Rigidbody;
    public float Area;
    
    SurfaceSampler2D m_SurfaceSampler;
    IFluidSim m_FluidSim;
    OutputSimulationSurfacePoints[] m_SimResults;
    InputSimulationSurfacePoints[] m_SimInput;

    void OnEnable()
    {
        m_SurfaceSampler = Rigidbody.GetComponent<SurfaceSampler2D>();
        m_FluidSim = GetComponent<IFluidSim>();
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

    Vector2 TransformPoint(Vector2 rigidBodyPoint)
    {
        var worldPoint = Rigidbody.transform.TransformPoint(rigidBodyPoint);
        var simLocalPoint = m_FluidSim.Transform.InverseTransformPoint(worldPoint);
        const float simWorldSize = 4f; //TODO: get this data instead of hard coding it
        const float cubeWorldSize = 1f; //TODO: get this data instead of hard coding it
        var rigidBodyToSimRatio = cubeWorldSize / simWorldSize; 
        simLocalPoint.x = m_FluidSim.Width / 2f + simLocalPoint.x * m_FluidSim.Width  * rigidBodyToSimRatio;
        simLocalPoint.y = m_FluidSim.Width / 2f + simLocalPoint.y * m_FluidSim.Height * rigidBodyToSimRatio;
        //Assert.AreEqual(simLocalPoint.z, 0);
        return simLocalPoint;
    }
    Vector2 TransformDirection(Vector2 rigidBodyDirection)
    {
        var worldDirection = Rigidbody.transform.TransformDirection(rigidBodyDirection);
        var simLocalDirection = m_FluidSim.Transform.InverseTransformDirection(worldDirection);
        //Assert.AreEqual(simLocalDirection.z, 0);
        return simLocalDirection;
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
                SimSpacePoint = TransformPoint(p.position),
                normal = TransformDirection(p.normal),
                velocity = CalcSurfacePointVelocity(worldPoint),
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
            Rigidbody.AddForceAtPosition(simResult.force * Time.deltaTime, worldPoint);
        }
    }

}
