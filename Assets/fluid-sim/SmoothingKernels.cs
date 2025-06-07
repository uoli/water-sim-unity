using System.Runtime.CompilerServices;
using UnityEngine;
namespace fluid_sim
{
    public static class SmoothingKernels
    {
        
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static float SmoothingKernel2(float dst, float radius)
        {
            if (dst >= radius) return 0;

            var scalingFactor = 6f /(Mathf.PI * Mathf.Pow(radius, 4));
            var dif = radius - dst;
            return dif * dif * scalingFactor; 
        }
    
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static float SmoothingKernel2Derivative(float dst, float radius)
        {
            if (dst >= radius) return 0;

            var scalingFactor = 12f / (Mathf.Pow(radius, 4) * Mathf.PI);
            var dif = radius - dst;

            return -dif * scalingFactor;
        }

        public static float CalcSmoothingKernelNormalization(float smoothingRadius)
        {
            return 4 / (Mathf.PI * Mathf.Pow(smoothingRadius, 8));
        }

        public static float CalcSmoothingKernelDerivativeNormalization(float smoothingRadius)
        {
            return -24 / (Mathf.PI * Mathf.Pow(smoothingRadius, 8));
        } 
    
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static float SmoothingKernel(float sqrdDistance, float squaredSmoothingLength, float kernelTerm)
        {
            if (sqrdDistance > squaredSmoothingLength) return 0;
            //var term = 4 / (Mathf.PI * Mathf.Pow(smoothingFactor, 8));
            var diff = squaredSmoothingLength - sqrdDistance;
            //return 315f / (64f * Mathf.PI * Mathf.Pow(smoothingFactor,9)) * Mathf.Pow(smoothingFactor*smoothingFactor - distance*distance, 3);
            return diff * diff * diff * kernelTerm;
        }
        
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static float SmoothingKernelDerivative(float dst, float sqrdDistance, float squaredSmoothingLength, float kernelDerivativeTerm)
        {
            if (sqrdDistance > squaredSmoothingLength) return 0;
            var diff = squaredSmoothingLength - sqrdDistance;
        
            return diff * diff *dst * kernelDerivativeTerm;
        }

    }
}
