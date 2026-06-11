# CapWindow-Gen

`CapWindow-Gen` is a fast cap-aware window scan for the raw general objective
`U_G = det(V_eq V_eq* + sigma I)`.

The solver is intentionally much simpler than `FrameOnly-Gen`:

1. Sort rows by antenna power.
2. Build a compact set of contiguous power windows of size `K`.
3. Always include the H3 strong/weak split window.
4. Evaluate each window with the exact `log det(sigma I + z(S)^2 G_S^2)` score.
5. Return the best-scoring window.

There is no local swap search and no external H3/threshold call. The expensive
part is one prefix sum of row Gram matrices, so the intended runtime class is
close to H3-style sorting rather than Frame local search.

## Reproducible Checks

Focused CDF run against the main dB competitors:

```powershell
venv/bin/python -m experiments.algorithm_comparison --N 1000 --L 4 --samples 20 --generator-seeds 10 42 --off-pcts 25 50 --algorithms H3 FrameOnly-Gen CapWindow-Gen --checkpoint-every 20 --out-dir results/cdf_cap_window_gen_smoke
```

Focused grid run across representative layer counts:

```powershell
venv/bin/python -m deprecated.experiments.grid_benchmark --N-values 1000 --L-values 1 2 4 8 10 --off-pcts 25 50 --samples 5 --save-runs --out-dir results/cap_window_gen_grid_smoke
```

Required unit test after algorithm or benchmark changes:

```powershell
venv/bin/python -m unittest motor_challenge_1205.py
```
