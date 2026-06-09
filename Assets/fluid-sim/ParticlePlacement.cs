using UnityEngine;

public enum ParticlePlacementMode
{
    // Regular lattice at the spacing that fills the box, plus a small random
    // offset. Initial densities start near uniform, so the sim settles gently
    // instead of exploding from random overlaps; the jitter breaks the lattice
    // symmetry to avoid crystalline artifacts.
    GridWithJitter,
    // Uniform random scatter. Produces Poisson clustering and large initial
    // density spikes; mostly useful for stress-testing stability.
    Random,
}

public static class ParticlePlacement
{
    const float k_JitterFraction = 0.05f;

    public static Vector2 GetPosition(ParticlePlacementMode mode, int index, int count, float width, float height)
    {
        switch (mode)
        {
            case ParticlePlacementMode.GridWithJitter:
                return GridWithJitter(index, count, width, height);
            case ParticlePlacementMode.Random:
            default:
                return new Vector2(UnityEngine.Random.Range(0f, width), UnityEngine.Random.Range(0f, height));
        }
    }

    static Vector2 GridWithJitter(int index, int count, float width, float height)
    {
        var spacing = Mathf.Sqrt(width * height / count);
        var columns = Mathf.Max(1, Mathf.CeilToInt(width / spacing));
        var jitter = spacing * k_JitterFraction;
        var x = (index % columns + 0.5f) * spacing + UnityEngine.Random.Range(-jitter, jitter);
        var y = (index / columns + 0.5f) * spacing + UnityEngine.Random.Range(-jitter, jitter);
        // The last row can land slightly past the top when the grid doesn't
        // divide the box exactly.
        return new Vector2(Mathf.Clamp(x, 0, width), Mathf.Clamp(y, 0, height));
    }
}
