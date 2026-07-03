# Threshold Full Sweep: nakagami

- N: 1000
- L: 2
- K values: 250, 500
- Samples: 100
- Generator seeds: 42
- Sigma: 1.0

The experiment sweeps every integer `T` from `0` to `K` and evaluates raw `U_G`.

## Answer

- `K=250`: best fixed `T=45`; 99% mean-`U_G` diapason `35..61`; best tested `T` median `42.5` (p05..p95 `20.9..70.0`).
- `K=500`: best fixed `T=59`; 99% mean-`U_G` diapason `45..79`; best tested `T` median `57.0` (p05..p95 `31.9..95.0`).

## Best Fixed Thresholds And Formula Checks

| K | best fixed T | 99% diapason | best tested T median | best tested T std | best formula | formula T | formula fraction |
|---:|---:|---|---:|---:|---|---:|---:|
| 250 | 45 | 35..61 | 42.500 | 16.281 | T_0p05N | 50 | 0.9695 |
| 500 | 59 | 45..79 | 57.000 | 19.792 | T_0p125NL_over_Lp2 | 62 | 0.9743 |

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
