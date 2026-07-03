# Threshold Full Sweep: nakagami

- N: 2000
- L: 2
- K values: 500, 1000
- Samples: 100
- Generator seeds: 42
- Sigma: 1.0

The experiment sweeps every integer `T` from `0` to `K` and evaluates raw `U_G`.

## Answer

- `K=500`: best fixed `T=92`; 99% mean-`U_G` diapason `65..121`; best tested `T` median `89.5` (p05..p95 `53.9..138.1`).
- `K=1000`: best fixed `T=129`; 99% mean-`U_G` diapason `89..159`; best tested `T` median `116.5` (p05..p95 `80.7..171.1`).

## Best Fixed Thresholds And Formula Checks

| K | best fixed T | 99% diapason | best tested T median | best tested T std | best formula | formula T | formula fraction |
|---:|---:|---|---:|---:|---|---:|---:|
| 500 | 92 | 65..121 | 89.500 | 25.677 | T_0p05N | 100 | 0.9812 |
| 1000 | 129 | 89..159 | 116.500 | 28.030 | T_0p125NL_over_Lp2 | 125 | 0.9825 |

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
