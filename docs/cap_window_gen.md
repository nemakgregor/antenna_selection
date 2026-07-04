# CapWindow-Gen and CapWindowFull-Gen

`CapWindow-Gen` and `CapWindowFull-Gen` are cap-aware power-window scans for the
raw general objective `U_G = det(V_eq V_eq* + sigma I)`.

Both solvers are intentionally simpler than `FrameOnly-Gen`:

1. Sort rows by antenna power.
2. Build contiguous power windows of size `K`.
3. Include the H3 strong/weak split window.
4. Evaluate each window with the exact `log det(sigma I + z(S)^2 G_S^2)` score.
5. Return the best-scoring window.

`CapWindow-Gen` evaluates a compact candidate grid. `CapWindowFull-Gen`
evaluates every contiguous window and is the current focused `L=2` requested-gen
variant. There is no local swap search and no external H3/threshold call.

## Reproducible Checks

Focused requested-gen CDF run:

```powershell
venv/bin/python -m experiments.algorithm_comparison --solver-set requested-gen --N 1000 --L 2 --samples 20 --generator-seeds 10 42 --off-pcts 25 50 --algorithms H3 H3ThresholdT123-Gen CapWindow-Gen CapWindowFull-Gen CapSubmod-Gen --checkpoint-every 20 --out-dir results/cdf_requested_gen_smoke
```

Broader CDF run using the default comparison set:

```powershell
venv/bin/python -m experiments.algorithm_comparison --N 1000 --L 4 --samples 20 --generator-seeds 10 42 --off-pcts 25 50 --algorithms H3 FrameOnly-Gen CapWindow-Gen CapSubmod-Gen --checkpoint-every 20 --out-dir results/cdf_cap_window_gen_smoke
```

Required unit test after algorithm or benchmark changes:

```powershell
venv/bin/python -m unittest motor_challenge_1205.py
```
