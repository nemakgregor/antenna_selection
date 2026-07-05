# Sigma Sweep

> Historical K semantics note: this report uses active-K semantics. Here `K` is the number of selected/kept antennas, not the number turned off. A `25% active` or `K=0.25N` case means `75% off`, not the real `25% off` task. For real off-percent experiments, `25% off => K_active=0.75N` and `50% off => K_active=0.50N`.

Use `--K-pcts` when `K` is an active antenna percentage:

```bash
venv/bin/python -m experiments.sigma_variation \
  --N 1000 --L 4 \
  --K-pcts 25 50 \
  --samples 30 \
  --sigmas 0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000 3000 10000 30000 100000 \
  --seed 42 \
  --checkpoint-every 1 \
  --summary-every 1 \
  --plot-every 0 \
  --out-dir results/sigma_sweep_N1000_L4_Kpct25_50_30samples
```

For direct comparison with `experiments.algorithm_comparison`, use
`--off-pcts 25 50` instead. That matches the CDF benchmark convention where
`25` means 25 percent disabled, so `K = round(N * 0.75)`.

Power-of-ten sigma repeat for 10 samples with 25 and 50 percent disabled:

```bash
venv/bin/python -m experiments.sigma_variation \
  --N 1000 --L 4 \
  --off-pcts 25 50 \
  --samples 10 \
  --sigmas 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1 1e1 1e2 1e3 \
  --seed 42 \
  --checkpoint-every 1 \
  --summary-every 1 \
  --plot-every 10 \
  --out-dir results/sigma_sweep_N1000_L4_offpct25_50_10samples_pow10
```

The script rewrites `sigma_sweep_runs.csv` and `sigma_sweep_progress.json` after
every completed `(sample, K, sigma)` case. Summary CSV files and the report are
controlled by `--summary-every`; plots are always written at the end and can also
be refreshed periodically with `--plot-every`. To rebuild report/plots from an
existing run without recomputing algorithms, use:

```bash
venv/bin/python -m experiments.sigma_variation \
  --N 1000 --L 4 \
  --K-pcts 25 50 \
  --samples 30 \
  --seed 42 \
  --out-dir results/sigma_sweep_N1000_L4_Kpct25_50_30samples \
  --refresh-from-runs results/sigma_sweep_N1000_L4_Kpct25_50_30samples/sigma_sweep_runs.csv
```

`U_G = det(V_eq V_eq* + sigma I)`. Equivalently, if `lambda_i` are eigenvalues of
`V_eq V_eq*`, then `U_G = product_i(lambda_i + sigma)`. `U_BF` and `U_I` do not
contain `sigma`; if their curves move across `sigma`, the algorithm selected a
different antenna set. With `L = 4`,
`U_G = sigma^4 + sigma^3 sum(lambda_i) + sigma^2 e2(lambda) + sigma e3(lambda)
+ product(lambda_i)`, so the first algorithm-dependent term at large `sigma`
tracks `U_BF`, not the interference objective.

All generated CSV, JSON, report, and plot files belong under `results/`, which
is ignored by Git.
