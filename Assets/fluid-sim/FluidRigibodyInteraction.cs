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
        Assert.AreEqual(simLocalPoint.z, 0);
        return simLocalPoint;
    }
    Vector2 TransformDirection(Vector2 rigidBodyDirection)
    {
        var worldDirection = Rigidbody.transform.TransformDirection(rigidBodyDirection);
        var simLocalDirection = m_FluidSim.Transform.InverseTransformDirection(worldDirection);
        Assert.AreEqual(simLocalDirection.z, 0);
        return simLocalDirection;
    }

    Vector2 CalcSurfacePointVelocity(Vector2 rigidBodyPoint)
    {
        Vector2 localOffset = rigidBodyPoint - Rigidbody.position;
        Vector2 surfacePointVelocity = Rigidbody.linearVelocity + 
            Rigidbody.angularVelocity * Vector2.Perpendicular(localOffset);
        return surfacePointVelocity;
    }

    void PreSimulation()
    {
        var points = m_SurfaceSampler.SurfacePoints;
        var areaWeight = Area / points.Count;
        for (var i = 0; i < points.Count; i++)
        {
            var p = points[i];
            
            m_SimInput[i] = new InputSimulationSurfacePoints
            {
                SimSpacePoint = TransformPoint(p.position),
                normal = TransformDirection(p.normal),
                velocity = CalcSurfacePointVelocity(p.position),
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
            Rigidbody.AddForceAtPosition(simResult.force,m_SurfaceSampler.SurfacePoints[index].position); //TODO: position in its own array would be beneficial here
        }
    }

}
