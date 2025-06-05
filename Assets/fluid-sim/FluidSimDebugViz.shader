Shader "Unlit/FluidSimDebugViz"
{
    Properties {}
    
    CGINCLUDE
    #pragma vertex vert
    #pragma fragment frag
    #pragma target 3.0

    #include "UnityCG.cginc"

    struct appdata
    {
        float4 vertex : POSITION;
        fixed4 color : COLOR;
        float2 uv : TEXCOORD0;
        UNITY_VERTEX_INPUT_INSTANCE_ID
    };

    struct v2f
    {
        float2 uv : TEXCOORD0;
        float4 vertex : SV_POSITION;
        float4 scrPos: UNITY_VPOS_TYPE ;
    };

    StructuredBuffer<float> _particle_positions;
    StructuredBuffer<float> _particle_densities;
    StructuredBuffer<float> _particle_pressures;
    int _PointCount;
    float _sizex;
    float _sizey;
    float _mousex;
    float _mousey;
    float _smoothingLength;
    float _DensityVizFactor;
    float _circleSize;
    int _visMode;
    float4 _negativePressureColor;
    float4 _neutralPressureColor;
    float4 _positivePressureColor;

    v2f vert (appdata v)
    {
        v2f o;
        UNITY_SETUP_INSTANCE_ID(v);
        o.vertex = UnityObjectToClipPos(v.vertex);
        o.scrPos = ComputeScreenPos(o.vertex);
        //o.uv = v.uv;
        o.uv = v.uv;
        return o;
    }

    static const float PI = 3.1415926535897932384626433;

    static const float smoothing_kernel_const = 4 / (PI * pow(_smoothingLength, 8));
    static const float squared_smooth_length = _smoothingLength * _smoothingLength;

    float SmoothingKernel(float squared_distance)
    {
        if (squared_distance > squared_smooth_length) return 0;
        //poly 6 kernel
        //float term = 315 / (64 * PI * pow(smoothingLength,9));
        //float term = 4 / (PI * pow(smoothingLength, 8));
        //return 315f / (64f * Mathf.PI * Mathf.Pow(smoothingFactor,9)) * Mathf.Pow(smoothingFactor*smoothingFactor - distance*distance, 3);
        float dif = squared_smooth_length - squared_distance; 
        return dif * dif * dif * smoothing_kernel_const;
    }

    float squaredDistance(float2 a, float2 b)
    {
        float2 diff = a-b;
        return dot(diff,diff);
    }

    float CalculateDensity(float2 pos)
    {
        float density = 0;
        float mass = 1;
        for (int i=0; i < _PointCount; i++)
        {
            float2 particlePos = float2(_particle_positions[i*2], _particle_positions[i*2+1]);
            float squaredDist = squaredDistance(particlePos,pos);
            //float dst = distance(particlePos, pos);
            float influence = SmoothingKernel(squaredDist);
            density += mass * influence;
        }

        return density;
    }

    float CalculatePressure(float2 pos)
    {
        float pressure = 0;
        float mass = 1;
        for (int i=0; i < _PointCount; i++)
        {
            float2 particlePos = float2(_particle_positions[i*2], _particle_positions[i*2+1]);
            float squaredDist = squaredDistance(particlePos,pos);
            //float dst = distance(particlePos, pos);
            float influence = SmoothingKernel(squaredDist);
            float density = _particle_densities[i];
            float particle_pressure = _particle_pressures[i];
            pressure += particle_pressure * mass / density * influence;
        }

        return pressure;
    }

    fixed4 frag (v2f i) : SV_Target
    {
        float2 size = float2(_sizex, _sizey);
        float2 mouse = float2(_mousex, _mousey);
        float2 normalizedScreenCoord = (i.scrPos.xy / i.scrPos.w); // [0-1] in viewport
        float2 screenPixel = normalizedScreenCoord *  _ScreenParams.xy ;
        float2 localPos = float2(screenPixel.x, _ScreenParams.y - screenPixel.y)  ;
        //fixed4 col = float4(1, frac(localPos / _ScreenParams.xy) , 1);

        fixed4 col;
        if (_visMode == 0) {
            float prop = CalculateDensity(localPos);
            col = float4(1,prop*_DensityVizFactor,0,1);
        } else if (_visMode == 1) {
            float prop = CalculatePressure(localPos);
            if (prop < -10)
                col = _negativePressureColor;
            else if (prop > 10)
                col = _positivePressureColor;
            else
                col = _neutralPressureColor;
        }

        

        float2 pos = mouse;
        if (distance(pos, localPos) <100)
        {
            col.b = 1;
        }
        

        for (int i = 0; i < _PointCount; i++)
        {
            float2 pos = float2(_particle_positions[i*2], _particle_positions[i*2+1]);
            
            if (distance(pos, localPos) <_circleSize)
            {
                float pressure = _particle_pressures[i];
                col = float4(0, pressure, 1, 1);
            }
        }

       
        return col;
    }

    

    
    
    ENDCG

    SubShader
    {
        Tags { "RenderType"="Opaque" }
        //LOD 100

        Pass
        {
            CGPROGRAM
            ENDCG
        }
    }
}
