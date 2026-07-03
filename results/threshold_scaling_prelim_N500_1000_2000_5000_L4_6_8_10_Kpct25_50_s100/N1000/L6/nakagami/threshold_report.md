# Threshold Full Sweep: nakagami

- N: 1000
- L: 6
- K values: 250, 500
- Samples: 100
- Generator seeds: 42
- Sigma: 1.0

The experiment sweeps every integer `T` from `0` to `K` and evaluates raw `U_G`.

## Answer

- `K=250`: best fixed `T=50`; 99% mean-`U_G` diapason `41..57`; best tested `T` median `46.0` (p05..p95 `28.9..79.3`).
- `K=500`: best fixed `T=71`; 99% mean-`U_G` diapason `57..88`; best tested `T` median `72.5` (p05..p95 `41.0..100.1`).

## Best Fixed Thresholds And Formula Checks

| K | best fixed T | 99% diapason | best tested T median | best tested T std | best formula | formula T | formula fraction |
|---:|---:|---|---:|---:|---|---:|---:|
| 250 | 50 | 41..57 | 46.000 | 16.431 | T_0p05N | 50 | 0.9527 |
| 500 | 71 | 57..88 | 72.500 | 19.664 | T_0p075N | 75 | 0.9527 |

## Plots

![U_G vs T](threshold_u_g_by_T.png)

![Best T histogram](threshold_best_T_hist.png)

![Best T boxplot](threshold_best_T_boxplot.png)

![Raw U_G CDF](threshold_raw_u_g_cdf.png)

![Fraction of best tested U_G CDF](threshold_fraction_best_cdf.png)

## Artifacts

- `threshold_runs.csv.gz`
- `best_thresholds.csv`
- `threshold_summary.csv`
- `threshold_best_t_stats.csv`
- `threshold_formula_comparison.csv`
