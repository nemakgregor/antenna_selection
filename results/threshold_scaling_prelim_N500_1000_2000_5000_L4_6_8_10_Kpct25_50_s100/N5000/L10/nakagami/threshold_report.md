# Threshold Full Sweep: nakagami

- N: 5000
- L: 10
- K values: 1250, 2500
- Samples: 100
- Generator seeds: 42
- Sigma: 1.0

The experiment sweeps every integer `T` from `0` to `K` and evaluates raw `U_G`.

## Answer

- `K=1250`: best fixed `T=255`; 99% mean-`U_G` diapason `218..296`; best tested `T` median `251.5` (p05..p95 `177.9..329.3`).
- `K=2500`: best fixed `T=391`; 99% mean-`U_G` diapason `308..449`; best tested `T` median `374.5` (p05..p95 `274.9..477.7`).

## Best Fixed Thresholds And Formula Checks

| K | best fixed T | 99% diapason | best tested T median | best tested T std | best formula | formula T | formula fraction |
|---:|---:|---|---:|---:|---|---:|---:|
| 1250 | 255 | 218..296 | 251.500 | 48.754 | T_0p05N | 250 | 0.9773 |
| 2500 | 391 | 308..449 | 374.500 | 60.622 | T_0p075N | 375 | 0.9760 |

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
