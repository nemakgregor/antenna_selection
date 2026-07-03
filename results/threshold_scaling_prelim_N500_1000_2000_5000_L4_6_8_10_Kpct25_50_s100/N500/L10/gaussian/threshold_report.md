# Threshold Full Sweep: gaussian

- N: 500
- L: 10
- K values: 125, 250
- Samples: 100
- Generator seeds: 42
- Sigma: 1.0

The experiment sweeps every integer `T` from `0` to `K` and evaluates raw `U_G`.

## Answer

- `K=125`: best fixed `T=22`; 99% mean-`U_G` diapason `21..22`; best tested `T` median `21.5` (p05..p95 `9.0..46.0`).
- `K=250`: best fixed `T=34`; 99% mean-`U_G` diapason `33..38`; best tested `T` median `32.0` (p05..p95 `18.0..54.1`).

## Best Fixed Thresholds And Formula Checks

| K | best fixed T | 99% diapason | best tested T median | best tested T std | best formula | formula T | formula fraction |
|---:|---:|---|---:|---:|---|---:|---:|
| 125 | 22 | 21..22 | 21.500 | 11.456 | T_0p05N | 25 | 0.8617 |
| 250 | 34 | 33..38 | 32.000 | 12.208 | T_0p075N | 38 | 0.8823 |

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
