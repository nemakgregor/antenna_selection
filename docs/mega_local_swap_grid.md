# Mega Local-Swap Grid

This run extends the unified local-swap comparison across multiple dimensions:

- `N=1000`
- `L=2 4 6 8 10`
- `off-pcts=25 50`
- `samples=100`
- `generator-seeds=42`
- profiles: `gaussian rayleigh rician nakagami lognormal thin_tail`
- sigmas: `1 10 100 1000 10000`
- methods: fixed threshold, cyclic best-T, `CapWindow-Gen`, `CapWindowFull-Gen`, `CapSubmod-Gen`, `Frame-Gen`, and `StrongWeak`, each with `0` and `1` local swap.

The sigma grid was chosen from `results/noise_scale_calibration_N1000_profiles6_s8/`.
For `N=1000`, the median RMS off-diagonal scale of `V_eff V_eff*` is about
`70..343` across `L=2..10`, so `1` is low noise, `100` is same-order noise,
`1000` is high noise, and `10000` is an extra high-noise guard point.

Canonical command:

```bash
venv/bin/python -m experiments.algorithm_comparison \
  --unified-local-swap-comparison \
  --N 1000 \
  --L-values 2 4 6 8 10 \
  --sigmas 1 10 100 1000 10000 \
  --off-pcts 25 50 \
  --samples 100 \
  --generator-seeds 42 \
  --data-profiles gaussian rayleigh rician nakagami lognormal thin_tail \
  --workers 1 \
  --compact-runs \
  --checkpoint-every 1000 \
  --resume \
  --out-dir results/mega_unified_local_swap_N1000_L2_4_6_8_10_profiles6_sigma1_10_100_1000_10000_s100
```

All generated CSV, JSON, report, and plot files should stay under `results/`.
