# Channel Distribution Profiles and Fading Connection

This report documents every synthetic channel profile currently used by the
experiments, the exact generator parameters in `utils/data.py`, and how those
profiles connect to fading models.

## Direct Answer

Fading is random variation of a wireless channel coefficient caused by
multipath propagation, phase addition/cancellation, movement, shadowing, and
path loss. In a MIMO channel matrix `V`, each entry is a complex channel
coefficient between one antenna and one stream/layer. A fading distribution is
therefore a statistical model for those complex coefficients or for their
amplitudes/powers.

In this repo, the distribution profile is used to generate a raw complex matrix
and then the matrix is normalized:

```python
column_norms = np.linalg.norm(V, axis=0)
V /= column_norms
antenna_max = np.max(np.linalg.norm(V, axis=1))
V /= antenna_max
```

So the final solver input `V` is not an exact textbook distribution sample.
Instead, the selected profile controls the relative channel/row-power shape
before normalization. This is enough for threshold experiments because the
heuristic depends mainly on sorted row powers, tail shape, and relative
structure.

## Profiles and Parameters

| profile | fading/distribution meaning | exact parameters used |
|---|---|---|
| `gaussian` | Normalized complex Gaussian baseline. Useful as a rich-scattering reference; complex Gaussian coefficients also produce Rayleigh-like magnitudes before normalization. Links: [normal distribution](https://en.wikipedia.org/wiki/Normal_distribution), [Rayleigh fading](https://en.wikipedia.org/wiki/Rayleigh_fading). | `real ~ N(0,1)`, `imag ~ N(0,1)`, `V = real + j*imag`, then `normalize_channel_matrix(V)`. |
| `rayleigh` | No dominant line-of-sight path; many scattered components with random phase. Link: [Rayleigh fading](https://en.wikipedia.org/wiki/Rayleigh_fading). | `amplitude ~ Rayleigh(scale=1.0)`, `phase ~ Uniform(0, 2*pi)`, `V = amplitude * exp(j*phase)`, then normalization. |
| `rician` | Dominant LOS/specular component plus diffuse scatter. Link: [Rician fading](https://en.wikipedia.org/wiki/Rician_fading). | `k_factor = 4.0`; `scatter = N(0,1) + j*N(0,1)`; one LOS phase per stream, shape `(1,L)`, shared across antennas; `V = sqrt(4/5)*los + sqrt(1/5)*scatter`, then normalization. |
| `nakagami` | Flexible fading-severity amplitude model. Link: [Nakagami distribution](https://en.wikipedia.org/wiki/Nakagami_distribution). | `m_shape = 2.0`, `omega = 1.0`; `amplitude_sq ~ Gamma(shape=2.0, scale=0.5)`; `amplitude = sqrt(amplitude_sq)`; `phase ~ Uniform(0,2*pi)`, then normalization. |
| `lognormal` | Shadowing / large-scale row-power variation stress profile. Links: [log-normal distribution](https://en.wikipedia.org/wiki/Log-normal_distribution), [shadow fading](https://en.wikipedia.org/wiki/Shadow_fading). | Per-row `shadow ~ LogNormal(mean=0.0, sigma=0.6)`, shape `(N,1)`; `scatter = N(0,1)+j*N(0,1)`; `V = shadow * scatter`, then normalization. |
| `thin_tail` | Synthetic bounded-tail stress case. It is not a standard named fading law; it checks whether the threshold rule still works when row-power spread is deliberately compressed. Link for the row scale law: [beta distribution](https://en.wikipedia.org/wiki/Beta_distribution). | Per-row `row_scale = 0.2 + 0.8 * Beta(a=2.0, b=8.0)`, shape `(N,1)`; `scatter = N(0,1)+j*N(0,1)`; `V = row_scale * scatter`, then normalization. |
| `twdp` | Two-Wave with Diffuse Power: two specular paths plus diffuse scatter. It is a multimodal/deep-fade stress model. Link: [Two-wave with diffuse power fading](https://en.wikipedia.org/wiki/Two-wave_with_diffuse_power_fading). | `twdp_K = 8.0`, `twdp_Delta = 0.9`, `Omega = 1.0`; `h = V1*exp(j*phi1) + V2*exp(j*phi2) + diffuse`; phases are independent uniform; diffuse is complex Gaussian; then normalization. |

## Derived TWDP Constants

The implemented TWDP defaults are:

```text
twdp_K = 8.0
twdp_Delta = 0.9
Omega = 1.0
```

The code derives:

```text
specular_power = Omega * twdp_K / (twdp_K + 1) = 0.8888888889
diffuse average power = Omega / (twdp_K + 1) = 0.1111111111
diffuse_sigma = sqrt(Omega / (2*(twdp_K + 1))) = 0.2357022604
V1 = 0.7988574882
V2 = 0.5007150912
V1^2 = 0.6381732864
V2^2 = 0.2507156025
```

Here `twdp_K` is not the antenna-selection active count `K`. It is the
specular-to-diffuse power ratio:

```text
twdp_K = (V1^2 + V2^2) / (2*sigma^2)
```

`twdp_Delta` describes balance between the two specular waves:

```text
twdp_Delta = 2*V1*V2 / (V1^2 + V2^2)
```

With `twdp_Delta = 0.9`, the two specular paths are strong and fairly balanced,
which makes destructive phase cancellation more likely than in a single-LOS
Rician model.

## Fading to `V` Matrix

For a MIMO system, a channel matrix entry can be written abstractly as:

```text
V[n, l] = channel coefficient from antenna n to stream/layer l
```

The coefficient is complex:

```text
V[n, l] = amplitude[n, l] * exp(j * phase[n, l])
```

Different fading assumptions define how `amplitude`, `phase`, and any
line-of-sight/specular components are sampled:

- Rayleigh: random amplitude, random phase, no dominant specular path.
- Rician: deterministic/specular component plus random scatter.
- Nakagami: flexible amplitude severity with random phase.
- Lognormal: slower row-level shadowing multiplied by small-scale scatter.
- TWDP: two specular paths plus diffuse scatter, allowing stronger cancellation
  and multimodal-looking amplitude/power behavior.

For antenna selection, the important derived quantity is often row power:

```text
p_n = sum_l |V[n, l]|^2
```

The threshold-window heuristic sorts antennas by `p_n`. Therefore, changing the
fading profile changes the distribution of row powers, the tail thickness, and
where a good threshold `T` tends to sit.

## Why These Profiles Are Useful

The profiles are not intended as measured channel data. They are a controlled
set of stress tests:

- `gaussian` and `rayleigh` test common diffuse-scattering baselines.
- `rician` tests a strong LOS/specular case; in the scaling experiments it was
  the profile that shifted best thresholds farthest right.
- `nakagami` tests a smoother fading-severity family.
- `lognormal` tests row-level shadowing and heavy row-power imbalance.
- `thin_tail` tests whether the rule survives a compressed right tail.
- `twdp` tests a physically motivated two-specular-path case with possible deep
  fades.

The normalization step makes comparisons fairer across profiles by preventing
simple total-power scale differences from dominating the objective. What
remains is mostly the relative geometry and row-power order, which is exactly
what threshold selection uses.

## Source Links

- Fading overview: https://en.wikipedia.org/wiki/Fading
- Rayleigh fading: https://en.wikipedia.org/wiki/Rayleigh_fading
- Rician fading: https://en.wikipedia.org/wiki/Rician_fading
- Nakagami distribution: https://en.wikipedia.org/wiki/Nakagami_distribution
- Log-normal distribution: https://en.wikipedia.org/wiki/Log-normal_distribution
- Shadow fading: https://en.wikipedia.org/wiki/Shadow_fading
- Beta distribution: https://en.wikipedia.org/wiki/Beta_distribution
- Normal distribution: https://en.wikipedia.org/wiki/Normal_distribution
- Two-wave with diffuse power fading: https://en.wikipedia.org/wiki/Two-wave_with_diffuse_power_fading
