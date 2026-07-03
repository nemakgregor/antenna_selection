# Threshold Full Sweep: nakagami

- N: 500
- L: 10
- K values: 125, 250
- Samples: 100
- Generator seeds: 42
- Sigma: 1.0

The experiment sweeps every integer `T` from `0` to `K` and evaluates raw `U_G`.

## Answer

- `K=125`: best fixed `T=24`; 99% mean-`U_G` diapason `21..27`; best tested `T` median `27.0` (p05..p95 `9.9..46.0`).
- `K=250`: best fixed `T=38`; 99% mean-`U_G` diapason `33..41`; best tested `T` median `35.0` (p05..p95 `21.0..58.0`).

## Best Fixed Thresholds And Formula Checks

| K | best fixed T | 99% diapason | best tested T median | best tested T std | best formula | formula T | formula fraction |
|---:|---:|---|---:|---:|---|---:|---:|
| 125 | 24 | 21..27 | 27.000 | 10.905 | T_0p05N | 25 | 0.8992 |
| 250 | 38 | 33..41 | 35.000 | 11.709 | T_0p075N | 38 | 0.9165 |

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
