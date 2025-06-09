fixed4 GetPointColor(int _particleVisMode, float pressure, float2 velocity, float _max_velocity)
{
    fixed4 color;
    switch (_particleVisMode)
    {
        case 0:
            default:
        {
            color = fixed4(0, 1, 1, 1);
            break;
        }
        case 1:
        {
            //float pressure = _particle_pressures[instanceID];
            color = fixed4(0, pressure, 1, 1);
            break;
        }
        case 2:
        {
            //float2 velocity = _particle_velocities[instanceID];
            float velMag = length(velocity);
            //velMag = velMag / _max_velocity;
            velMag = log(velMag + 1) / log(_max_velocity + 1); // maps [0, maxValue] -> [0,1]
            //velMag = pow(velMag / _max_velocity, 2);
            color = fixed4(velMag, 0, 1-velMag, 1);
            break;
        }
                        
    }
    return color;
}