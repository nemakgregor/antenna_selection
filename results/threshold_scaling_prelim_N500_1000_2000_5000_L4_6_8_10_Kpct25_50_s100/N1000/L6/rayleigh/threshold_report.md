# Threshold Full Sweep: rayleigh

- N: 1000
- L: 6
- K values: 250, 500
- Samples: 100
- Generator seeds: 42
- Sigma: 1.0

The experiment sweeps every integer `T` from `0` to `K` and evaluates raw `U_G`.

## Answer

- `K=250`: best fixed `T=48`; 99% mean-`U_G` diapason `36..56`; best tested `T` median `47.5` (p05..p95 `25.9..79.0`).
- `K=500`: best fixed `T=69`; 99% mean-`U_G` diapason `55..76`; best tested `T` median `62.5` (p05..p95 `38.0..102.3`).

## Best Fixed Thresholds And Formula Checks

| K | best fixed T | 99% diapason | best tested T median | best tested T std | best formula | formula T | formula fraction |
|---:|---:|---|---:|---:|---|---:|---:|
| 250 | 48 | 36..56 | 47.500 | 17.041 | T_0p05N | 50 | 0.9254 |
| 500 | 69 | 55..76 | 62.500 | 19.975 | T_0p075N | 75 | 0.9315 |

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
