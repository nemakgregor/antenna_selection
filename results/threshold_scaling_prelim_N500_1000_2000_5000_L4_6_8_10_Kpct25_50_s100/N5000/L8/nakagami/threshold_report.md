# Threshold Full Sweep: nakagami

- N: 5000
- L: 8
- K values: 1250, 2500
- Samples: 100
- Generator seeds: 42
- Sigma: 1.0

The experiment sweeps every integer `T` from `0` to `K` and evaluates raw `U_G`.

## Answer

- `K=1250`: best fixed `T=256`; 99% mean-`U_G` diapason `206..312`; best tested `T` median `263.0` (p05..p95 `185.5..371.1`).
- `K=2500`: best fixed `T=358`; 99% mean-`U_G` diapason `300..438`; best tested `T` median `360.0` (p05..p95 `267.0..446.3`).

## Best Fixed Thresholds And Formula Checks

| K | best fixed T | 99% diapason | best tested T median | best tested T std | best formula | formula T | formula fraction |
|---:|---:|---|---:|---:|---|---:|---:|
| 1250 | 256 | 206..312 | 263.000 | 55.006 | T_0p05N | 250 | 0.9755 |
| 2500 | 358 | 300..438 | 360.000 | 59.403 | T_0p075N | 375 | 0.9792 |

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
