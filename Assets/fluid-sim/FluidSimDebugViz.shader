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
    StructuredBuffer<float> _particle_velocities;
    int _PointCount;
    float _scaling_factor;
    float _sizex;
    float _sizey;
    float _mousex;
    float _mousey;
    int _mousepressed;
    float _mouseradius;
    float _smoothingLength;
    float _DensityVizFactor;
    float _circleSize;
    int _fieldVisMode;
    int _particleVisMode;
    float _min_pressure;
    float _max_pressure;
    float _target_density;
    float _max_density;
    float _max_velocity;
    float _min_velocity;
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
    static const float mass = 1;

    static const float smoothing_kernel_const = 4 / (PI * pow(_smoothingLength, 8));
    static const float squared_smooth_length = _smoothingLength * _smoothingLength;
    static const float kernel_derivative_term = -24 / (PI * pow(_smoothingLength, 8));

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
    
    static float SmoothingKernelDerivative(float dst, float sqrdDistance)
    {
        if (sqrdDistance > squared_smooth_length) return 0;
        float diff = squared_smooth_length - sqrdDistance;
        
        return diff*diff *dst * kernel_derivative_term;
    }

    static float SmoothingKernel2(float dst, float radius)
    {
        if (dst > radius) return 0;

        float volume = (PI * pow(radius, 4)) / 6;
        return (radius - dst) * (radius - dst) / volume;
    }
    
    static float SmoothingKernel2Derivative(float dst, float radius)
    {
        if (dst > radius) return 0;

        float scale = 12 / (pow(radius, 4) * PI);
        return (dst - radius ) * scale;
    }

    float squaredDistance(float2 a, float2 b)
    {
        float2 diff = a-b;
        return dot(diff,diff);
    }

    float CalculateDensity(float2 pos)
    {
        float density = 0;
        for (int i=0; i < _PointCount; i++)
        {
            float2 particlePos = float2(_particle_positions[i*2], _particle_positions[i*2+1]);
            float squaredDist = squaredDistance(particlePos,pos);
            //float dst = distance(particlePos, pos);
            //float influence = SmoothingKernel(squaredDist);
            float influence = SmoothingKernel2(sqrt(squaredDist), _smoothingLength);
            density += mass * influence;
        }

        return density;
    }

    float CalculatePressure(float2 pos)
    {
        float pressure = 0;
        for (int i=0; i < _PointCount; i++)
        {
            float2 particlePos = float2(_particle_positions[i*2], _particle_positions[i*2+1]);
            float squaredDist = squaredDistance(particlePos,pos);
            //float dst = distance(particlePos, pos);
            //float influence = SmoothingKernel(squaredDist);
            float influence = SmoothingKernel2(sqrt(squaredDist), _smoothingLength);
            float density = _particle_densities[i];
            float particle_pressure = _particle_pressures[i];
            pressure += particle_pressure * mass / density * influence;
        }

        return pressure;
    }

    float2 CalculatePressureGradient(float2 pos)
    {
        float2 pressureGradient = 0;
        for (int i = 0; i < _PointCount; i++)
        {
            float2 particlePos = float2(_particle_positions[i*2], _particle_positions[i*2+1]);
            float sqrDst = squaredDistance(particlePos,pos );
            if (sqrDst < 0.001 || sqrDst > squared_smooth_length) continue;
            float2 dif = particlePos - pos;
            float2 dir = normalize(dif);
            float pressure = _particle_pressures[i];
            float density = _particle_densities[i];
            float distance = sqrt(sqrDst);
            //float influence = SmoothingKernelDerivative(distance, sqrDst);
            float influence = SmoothingKernel2Derivative(distance, _smoothingLength);
            pressureGradient += dir * (pressure * mass) / density * influence;
        }

        return pressureGradient;
    }

    float inverseLerp(float min, float max, float val)
    {
        return (val - min) /  (max - min);
    }

    fixed4 frag (v2f i) : SV_Target
    {
        //float2 size = float2(_sizex, _sizey);
        float2 mouse = float2(_mousex, _mousey);
        float2 normalizedScreenCoord = (i.scrPos.xy / i.scrPos.w); // [0-1] in viewport
        float2 screenPixel = normalizedScreenCoord *  _ScreenParams.xy ;
        float2 localPos = float2(screenPixel.x,  screenPixel.y - _ScreenParams.y + _sizey * _scaling_factor ) / _scaling_factor;
        //fixed4 col = float4(1, frac(localPos / _ScreenParams.xy) , 1);
        _min_pressure = min(_min_pressure, -1);

        fixed4 col;
        if (_fieldVisMode == 0) {
            float prop = CalculateDensity(localPos);
            if (prop < _target_density)
            {
                float t = inverseLerp(0, _target_density, prop);
                col = lerp(_negativePressureColor,_neutralPressureColor, t);
            }
            else if (prop > _target_density)
            {
                float t = inverseLerp(_target_density, _max_density, prop);
                col = lerp(_neutralPressureColor,_positivePressureColor, t);
            }
            else
            {
                col = _neutralPressureColor;
            }
        } else if (_fieldVisMode == 1) {
            float prop = CalculatePressure(localPos);
            
            float t = inverseLerp(_min_pressure, _max_pressure, prop);
            col = lerp(_negativePressureColor,_positivePressureColor, t);
           
            //col = float4(1,prop*_DensityVizFactor,0,1);
        } else if (_fieldVisMode == 2) {
            float2 prop = CalculatePressureGradient(localPos);
            col = fixed4(prop * _DensityVizFactor, 0, 1);
        } else if (_fieldVisMode == 3)
        {
            col = fixed4(0,0,0,1);
        }

        

        if (_mousepressed)
        {
            float2 pos = mouse;
            if (distance(pos, localPos) < _mouseradius / _scaling_factor)
            {
                col.b = 1;
            }
        }
        

        for (int i = 0; i < _PointCount; i++)
        {
            float2 pos = float2(_particle_positions[i*2], _particle_positions[i*2+1]);
            
            if (distance(pos, localPos) <_circleSize / _scaling_factor)
            {
                // Solid = 0,
                // Pressure = 1,
                // Velocity = 2,
                switch (_particleVisMode)
                {
                    case 0:
                    default:
                    {
                        col = float4(0, 1, 1, 1);
                        break;
                    }
                    case 1:
                    {
                        float pressure = _particle_pressures[i];
                        col = float4(0, pressure, 1, 1);
                        break;
                    }
                    case 2:
                    {
                        float2 vel = float2(
                            _particle_velocities[i*2],
                            _particle_velocities[i*2+1]);
                        float velMag = length(vel);
                        //velMag *= .5;
                        //velMag = velMag + _min_velocity;
                        velMag = velMag / _max_velocity;
                        //vel = 1 + vel * 0.5;
                        col = float4(velMag, 0, 1-velMag, 1);
                        break;
                    }
                        
                }

                
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
