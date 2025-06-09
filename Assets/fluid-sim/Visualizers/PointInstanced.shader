Shader "Unlit/PointsInstanced"
{
    Properties
    {
        _Color("Color", Color) = (1,1,1,1)
    }
    SubShader
    {
        Tags { "RenderType"="Transparent" "Queue"="Transparent" }
        Blend SrcAlpha OneMinusSrcAlpha
        ZWrite Off
        ZTest Always

        Cull Off
        LOD 100

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma multi_compile_instancing
            #include "UnityCG.cginc"
            #include "vizCommon.cginc"

            UNITY_INSTANCING_BUFFER_START(Props)
                UNITY_DEFINE_INSTANCED_PROP(float4, _Color)
            UNITY_INSTANCING_BUFFER_END(Props)

            StructuredBuffer<float2> positions;
            StructuredBuffer<float> _particle_densities;
            StructuredBuffer<float> _particle_pressures;
            StructuredBuffer<float2> _particle_velocities;
            float circle_size;
            float scaling_factor;
            float sim_height;
            float sim_width;
            float4 world_origin;
            float display_area_width;
            float display_area_height;
            int _particleVisMode;
            float _max_velocity;


            
            struct appdata
            {
                float4 vertex : POSITION;
                UNITY_VERTEX_INPUT_INSTANCE_ID
            };

            struct v2f
            {
                float4 vertex : SV_POSITION;
                float2 local_position : TEXCOORD0;
                fixed4 color : COLOR;
                UNITY_VERTEX_OUTPUT_STEREO
            };

            v2f vert(appdata v, uint instanceID : SV_InstanceID)
            {
                UNITY_SETUP_INSTANCE_ID(v);

                float2 point_position = positions[instanceID];
                
                float4 pos = float4(point_position * (display_area_width /sim_width  ),0,0) * scaling_factor;
                float4 worldPos = world_origin + v.vertex * circle_size * display_area_width /sim_width + pos;
                float4 fvertex = UnityObjectToClipPos(worldPos);
                
                
                // float2 screenSize = _ScreenParams.xy;
                // float2 pixelPos = positions[instanceID] * scaling_factor; // e.g. (100, 200)
                //
                // float2 centerClip = (pixelPos / screenSize) * 2.0 - 1.0;  // map pixel pos to [-1,1]
                // centerClip.y *= -1.0; // flip y axis because clip y goes up, screen y down
                //
                // float2 offset = v.vertex.xy * circle_size; 
                // // scale quad vertices by circle_size (in pixels) then convert to clip space units
                //
                // float2 vertexClip = centerClip + offset;
                // float4 fvertex = float4(vertexClip, 0.5, 1);

                fixed4 color = _Color;
                
                v2f o;
                UNITY_TRANSFER_INSTANCE_ID(v, o);
                o.vertex = fvertex;
                o.local_position = v.vertex;
                
                color = GetPointColor(_particleVisMode, _particle_pressures[instanceID],  _particle_velocities[instanceID],  _max_velocity);
                o.color = color;
                
                return o;
            }

            fixed4 frag(v2f i) : SV_Target
            {
                float distance_to_center = 0.5 - distance(i.local_position,0);
                float s = smoothstep(0.0,0.1,distance_to_center);
                return i.color * s;
            }
            ENDCG
        }
    }
}
