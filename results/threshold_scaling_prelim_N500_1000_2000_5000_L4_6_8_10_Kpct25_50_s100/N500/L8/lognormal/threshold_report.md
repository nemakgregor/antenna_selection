# Threshold Full Sweep: lognormal

- N: 500
- L: 8
- K values: 125, 250
- Samples: 100
- Generator seeds: 42
- Sigma: 1.0

The experiment sweeps every integer `T` from `0` to `K` and evaluates raw `U_G`.

## Answer

- `K=125`: best fixed `T=15`; 99% mean-`U_G` diapason `14..16`; best tested `T` median `19.0` (p05..p95 `7.0..33.0`).
- `K=250`: best fixed `T=22`; 99% mean-`U_G` diapason `22..22`; best tested `T` median `27.0` (p05..p95 `8.9..50.1`).

## Best Fixed Thresholds And Formula Checks

| K | best fixed T | 99% diapason | best tested T median | best tested T std | best formula | formula T | formula fraction |
|---:|---:|---|---:|---:|---|---:|---:|
| 125 | 15 | 14..16 | 19.000 | 8.435 | T_0p05NL_over_Lp2 | 20 | 0.6572 |
| 250 | 22 | 22..22 | 27.000 | 11.662 | T_0p075NL_over_Lp2 | 30 | 0.6972 |

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
