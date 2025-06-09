Shader "Unlit/SimpleInstancedPoints"
{
    Properties
    {
        _Color("Color", Color) = (1,1,1,1)
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100

        Pass
        {
            ZTest Always
            ZWrite Off
    
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma multi_compile_instancing
            #include "UnityCG.cginc"

            UNITY_INSTANCING_BUFFER_START(Props)
                UNITY_DEFINE_INSTANCED_PROP(float4, _Color)
            UNITY_INSTANCING_BUFFER_END(Props)

            StructuredBuffer<float4> positions;

            struct appdata
            {
                float4 vertex : POSITION;
                UNITY_VERTEX_INPUT_INSTANCE_ID
            };

            struct v2f
            {
                float4 vertex : SV_POSITION;
                fixed4 color : COLOR;
                UNITY_VERTEX_OUTPUT_STEREO
            };

            v2f vert(appdata v, uint instanceID : SV_InstanceID)
            {
                UNITY_SETUP_INSTANCE_ID(v);
                v2f o;
                UNITY_INITIALIZE_VERTEX_OUTPUT_STEREO(o);

                float4 color = UNITY_ACCESS_INSTANCED_PROP(Props, _Color);
                o.color = color;
                
                //float4 pos = UNITY_ACCESS_INSTANCED_PROP(Props, _Position);
                float4 pos = positions[instanceID];
                float4 worldPos = v.vertex + pos;
                o.vertex = UnityObjectToClipPos(worldPos);

                // float2 screenSize = _ScreenParams.xy;
                // float2 pixelPos = positions[instanceID]; // e.g. (100, 200)
                //
                // float2 centerClip = (pixelPos / screenSize) * 2.0 - 1.0;  // map pixel pos to [-1,1]
                // centerClip.y *= -1.0; // flip y axis because clip y goes up, screen y down
                //
                // float2 offset = v.vertex.xy; 
                // // scale quad vertices by circle_size (in pixels) then convert to clip space units
                //
                // float2 vertexClip = centerClip + offset;
                // float4 fvertex = float4(vertexClip, 0.5, 1);
                // o.vertex = fvertex;

                return o;
            }

            fixed4 frag(v2f i) : SV_Target
            {
                return i.color;
            }
            ENDCG
        }
    }
}
