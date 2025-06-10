using System;
using UnityEditor;
using UnityEngine;

[CustomEditor(typeof(SurfaceSampler2D))]
public class SurfaceSampler2DEditor : Editor
{
    protected virtual void OnSceneGUI()
    {
        SurfaceSampler2D t = (SurfaceSampler2D)target;
        
        using var hc = new HandlesColorContext(new Color(1,1,0,1f));
        float size = HandleUtility.GetHandleSize(t.transform.position) * 0.2f;
        var snap = 0.1f;
        var handleDir1 = t.transform.TransformDirection(Vector3.forward);
        var slideDir1 = t.transform.TransformDirection(Vector3.left);
        var slideDir2 = t.transform.TransformDirection(Vector3.up);


        for (var i = 0; i < t.SurfacePoints.Count; i++)
        {
            var surfacePoint = t.SurfacePoints[i];
            var center = t.transform.TransformPoint(surfacePoint.position);

            EditorGUI.BeginChangeCheck();
            var newPosition = Handles.Slider2D(center, handleDir1, slideDir1, slideDir2, size, Handles.CircleHandleCap, snap);
            if (EditorGUI.EndChangeCheck())
            {
                Undo.RecordObject(t, "Change Surface Sample position");
                surfacePoint.position = t.transform.InverseTransformPoint(newPosition);
                t.SurfacePoints[i] = surfacePoint;
            }
            
            EditorGUI.BeginChangeCheck();
            
            var newNormal = Normal2DRotationHandle(center, handleDir1, surfacePoint.normal, size*2);
            if (EditorGUI.EndChangeCheck())
            {
                Undo.RecordObject(t, "Change Surface Sample position");
                surfacePoint.normal = newNormal.normalized;
                t.SurfacePoints[i] = surfacePoint;
            }
            
            var line2a = center - slideDir2 * size;
            var line2b = center + slideDir2 * size;
            var line1a = center - slideDir1 * size;
            var line1b = center + slideDir1 * size;
            Handles.DrawLine(line1a, line1b);
            Handles.DrawLine(line2a, line2b);
        }
    }
    
    Vector2 Normal2DRotationHandle(Vector3 position, Vector3 planeNormal, Vector2 normal2D, float size)
    {
        var controlID = GUIUtility.GetControlID(FocusType.Passive);
        var evt = Event.current;
        
        switch (evt.GetTypeForControl(controlID))
        {
            case EventType.MouseDown:
                if (controlID == HandleUtility.nearestControl && evt.button == 0 && GUIUtility.hotControl == 0)
                {
                    GUIUtility.hotControl = controlID;
                    evt.Use();
                }
                break; 
            case EventType.MouseUp:
                if (GUIUtility.hotControl == controlID && evt.button == 0)
                {
                    GUIUtility.hotControl = 0;
                    evt.Use();
                }
                break;
            case EventType.MouseDrag:
                if (GUIUtility.hotControl == controlID)
                {
                    Plane plane = new Plane(planeNormal, position);
                    Ray ray = HandleUtility.GUIPointToWorldRay(Event.current.mousePosition);
                    plane.Raycast(ray, out float enter);
                    var planePoint =  ray.GetPoint(enter);
                    var diff = planePoint - position;
                    GUI.changed = true;
                    normal2D = diff;
                    evt.Use();
                }
                break;
            case EventType.Layout:
            case EventType.MouseMove:
                // Set the nearest control ID to "controlID" if the mouse cursor is
                // near or directly above the solid disc handle.
                
                //Transform normal to worldSpace (there must be a better way)
                Vector3 tangent1 = Vector3.Cross(planeNormal, Vector3.up);
                if (tangent1.magnitude < 0.1f) // Handle case where normal is up/down
                    tangent1 = Vector3.Cross(planeNormal, Vector3.right);
                tangent1.Normalize();

                Vector3 tangent2 = Vector3.Cross(planeNormal, tangent1).normalized;
                Vector3 normalWorldPos = position - 
                    normal2D.x * size * tangent1 - 
                    normal2D.y * size * tangent2;
                
                var distanceToHandle = HandleUtility.DistanceToLine(position, normalWorldPos);
                HandleUtility.AddControl(controlID, distanceToHandle);
                break;
            case EventType.Repaint:
                var quat = Quaternion.LookRotation(normal2D, planeNormal);

                // Display an orange solid disc where the object is.
                //Handles.color = new Color(0.5f, 0.8f, 0.4f, 1);
                //Handles.DrawSolidDisc(pos, tr.up, t.value);
                Handles.ArrowHandleCap(controlID, position, quat, size, EventType.Repaint);

                break;
        }

        return normal2D;
    }
}

    




struct HandlesColorContext : IDisposable
{
    Color m_PrevColor;
    public HandlesColorContext(Color color)
    {
        m_PrevColor = Handles.color;
        Handles.color = color;
    }

    public void Dispose()
    {
        Handles.color = m_PrevColor;
    }
}