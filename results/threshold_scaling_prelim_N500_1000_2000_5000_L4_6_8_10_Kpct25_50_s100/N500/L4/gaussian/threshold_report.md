# Threshold Full Sweep: gaussian

- N: 500
- L: 4
- K values: 125, 250
- Samples: 100
- Generator seeds: 42
- Sigma: 1.0

The experiment sweeps every integer `T` from `0` to `K` and evaluates raw `U_G`.

## Answer

- `K=125`: best fixed `T=21`; 99% mean-`U_G` diapason `18..25`; best tested `T` median `20.5` (p05..p95 `8.0..37.0`).
- `K=250`: best fixed `T=30`; 99% mean-`U_G` diapason `24..37`; best tested `T` median `28.5` (p05..p95 `14.0..51.0`).

## Best Fixed Thresholds And Formula Checks

| K | best fixed T | 99% diapason | best tested T median | best tested T std | best formula | formula T | formula fraction |
|---:|---:|---|---:|---:|---|---:|---:|
| 125 | 21 | 18..25 | 20.500 | 8.936 | T_0p05N | 25 | 0.9079 |
| 250 | 30 | 24..37 | 28.500 | 11.733 | T_0p10NL_over_Lp2 | 33 | 0.9225 |

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
