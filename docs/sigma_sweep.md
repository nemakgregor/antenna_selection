# Sigma Sweep

Use `--K-pcts` when `K` is an active antenna percentage:

```bash
venv/bin/python sigma_sweep_analysis.py \
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

For direct comparison with `grid_benchmark`, use `--off-pcts 25 50` instead.
That matches the benchmark convention where `25` means 25 percent disabled, so
`K = round(N * 0.75)`.

The script rewrites `sigma_sweep_runs.csv` and `sigma_sweep_progress.json` after
every completed `(sample, K, sigma)` case. Summary CSV files and the report are
controlled by `--summary-every`; plots are always written at the end and can also
be refreshed periodically with `--plot-every`. To rebuild report/plots from an
existing run without recomputing algorithms, use:

```bash
venv/bin/python sigma_sweep_analysis.py \
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
