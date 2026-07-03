# Threshold Full Sweep: twdp

- N: 500
- L: 6
- K values: 125, 250
- Samples: 100
- Generator seeds: 42
- Sigma: 1.0

The experiment sweeps every integer `T` from `0` to `K` and evaluates raw `U_G`.

## Answer

- `K=125`: best fixed `T=20`; 99% mean-`U_G` diapason `17..27`; best tested `T` median `21.0` (p05..p95 `7.0..35.0`).
- `K=250`: best fixed `T=27`; 99% mean-`U_G` diapason `24..32`; best tested `T` median `27.0` (p05..p95 `12.9..54.0`).

## Best Fixed Thresholds And Formula Checks

| K | best fixed T | 99% diapason | best tested T median | best tested T std | best formula | formula T | formula fraction |
|---:|---:|---|---:|---:|---|---:|---:|
| 125 | 20 | 17..27 | 21.000 | 8.871 | T_0p05N | 25 | 0.9292 |
| 250 | 27 | 24..32 | 27.000 | 11.583 | T_0p075NL_over_Lp2 | 28 | 0.9388 |

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
