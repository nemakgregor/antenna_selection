# Fading Distributions and the TWDP Stress Profile

## Direct Answer

Fading is random variation in wireless channel gain. It appears because the
received signal is usually a sum of many propagated copies of the transmitted
signal. Those copies can have different attenuation, delay, phase, reflection
paths, movement-induced Doppler shifts, and shadowing. Sometimes they add
constructively; sometimes they cancel and create a deep fade.

We use fading distributions because the exact channel state is not stable or
fully known before transmission. For antenna-selection experiments, a
distribution is a practical way to generate many plausible channel matrices and
ask whether a heuristic works statistically, not only on one lucky matrix.

## How This Connects to Our Profiles

The current experiment profiles are synthetic channel-matrix generators. They
are not measured channel data, but each one stresses a different plausible
shape of row power and channel gain.

| profile | interpretation |
|---|---|
| `gaussian` | Normalized complex Gaussian baseline. This is the original repo generator and is kept unchanged. |
| `rayleigh` | No dominant line-of-sight component; many scattered components with random phases. |
| `rician` | One dominant line-of-sight or strong reflected component plus diffuse scatter. |
| `nakagami` | Flexible fading-severity model; useful because it can emulate several fading severities with one shape parameter. |
| `lognormal` | Shadowing or large-scale row-power variation stress profile. |
| `thin_tail` | Synthetic bounded-tail stress profile for checking whether the threshold rule depends on very heavy right tails. |
| `twdp` | Two-Wave with Diffuse Power fading: two strong specular paths plus diffuse scatter. This is the multimodal/deep-fade stress case. |

The `twdp` profile is experimental. It should not replace the existing
profiles; it expands the threshold study to a physically motivated channel
shape that can be more severe than Rayleigh or Rician.

## Figure 6 in the Paper

The local paper is `publiucations/Power allocation algorithm for capacity maximization in 5G MIMO systems.pdf`.
Its Fig. 6 caption is "Comparisons of MIMO channels with PDFs." The visible
y-axis says it is the PDF of `lambda` elements in the SVD decomposition of
channel matrix `H`, and the legend overlays several `(nt, nr)` MIMO
configurations.

So Fig. 6 should be read carefully:

- The multi-peak visual shape is not necessarily one fading-envelope
  distribution.
- It can come from overlaying several MIMO configurations on one plot.
- It can also come from singular-value/eigenmode statistics of the channel
  matrix, not directly from the scalar fading amplitude distribution.
- For our experiments, the figure is good motivation to test non-simple channel
  shapes, but it is not enough to infer a measured multimodal fading law.

This is why `twdp` is a better next profile than an arbitrary mixture of the
existing generators: it is a known fading model with an explicit physical
interpretation.

## TWDP Parameters

In this report, `twdp_K` is different from the antenna-selection active count
`K`.

The TWDP channel coefficient is generated as:

```text
h = V1 * exp(j phi1) + V2 * exp(j phi2) + diffuse
```

where `phi1` and `phi2` are uniform random phases and `diffuse` is complex
Gaussian scatter.

The two main TWDP parameters are:

```text
twdp_K = (V1^2 + V2^2) / (2 sigma^2)
twdp_Delta = 2 V1 V2 / (V1^2 + V2^2)
```

- `twdp_K` is the specular-to-diffuse power ratio. Larger values mean the two
  strong paths dominate the diffuse scatter more.
- `twdp_Delta` measures the balance between the two specular paths. `0` reduces
  toward one dominant specular path; `1` means two equal-strength specular
  paths.

The implemented experimental defaults are:

```text
twdp_K = 8.0
twdp_Delta = 0.9
Omega = 1.0
```

This means the profile has two strong, nearly balanced paths with enough
diffuse scatter to remain a fading model. It is intentionally a stress case:
the two strong paths can cancel when phases oppose each other, producing deeper
fades and a less simple distribution shape than Rayleigh/Rician/Nakagami.

The implementation derives:

```text
V1^2 + V2^2 = Omega * twdp_K / (twdp_K + 1)
2 sigma^2 = Omega / (twdp_K + 1)
```

Then `V1` and `V2` are chosen deterministically from `twdp_Delta`, so the
profile has no hidden random parameter drift. The generated matrix is finally
passed through the same normalization as the other profiles.

## Experiment Commands

Smoke run:

```bash
venv/bin/python -m experiments.algorithm_comparison \
  --threshold-full-sweep \
  --N 100 \
  --L 2 \
  --K-values 50 25 \
  --samples 2 \
  --generator-seeds 42 \
  --data-profiles twdp \
  --sigma 1 \
  --out-dir results/smoke_threshold_full_sweep_twdp
```

Full run with the new profile included:

```bash
venv/bin/python -m experiments.algorithm_comparison \
  --threshold-full-sweep \
  --N 1000 \
  --L 2 \
  --K-values 500 250 \
  --samples 1000 \
  --generator-seeds 42 \
  --data-profiles gaussian rayleigh rician nakagami lognormal thin_tail twdp \
  --sigma 1 \
  --out-dir results/threshold_full_sweep_L2_N1000_K250_500_with_twdp_s1000
```

Rule CDF/report regeneration:

```bash
venv/bin/python -m experiments.algorithm_comparison \
  --threshold-rule-cdf \
  --N 1000 \
  --L 2 \
  --K-values 500 250 \
  --samples 1000 \
  --generator-seeds 42 \
  --data-profiles gaussian rayleigh rician nakagami lognormal thin_tail twdp \
  --sigma 1 \
  --out-dir results/threshold_full_sweep_L2_N1000_K250_500_with_twdp_s1000
```

## Sources

- Local paper: `publiucations/Power allocation algorithm for capacity maximization in 5G MIMO systems.pdf`
- Fading overview: https://en.wikipedia.org/wiki/Fading
- Rayleigh fading: https://en.wikipedia.org/wiki/Rayleigh_fading
- Rician fading: https://en.wikipedia.org/wiki/Rician_fading
- Nakagami distribution: https://en.wikipedia.org/wiki/Nakagami_distribution
- Two-Wave with Diffuse Power fading: https://en.wikipedia.org/wiki/Two-wave_with_diffuse_power_fading
- Wishart distribution and MIMO eigenvalue context: https://en.wikipedia.org/wiki/Wishart_distribution
