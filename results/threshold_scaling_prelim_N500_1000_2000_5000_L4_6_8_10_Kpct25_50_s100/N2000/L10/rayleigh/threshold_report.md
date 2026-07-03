# Threshold Full Sweep: rayleigh

- N: 2000
- L: 10
- K values: 500, 1000
- Samples: 100
- Generator seeds: 42
- Sigma: 1.0

The experiment sweeps every integer `T` from `0` to `K` and evaluates raw `U_G`.

## Answer

- `K=500`: best fixed `T=103`; 99% mean-`U_G` diapason `84..118`; best tested `T` median `95.5` (p05..p95 `53.9..152.1`).
- `K=1000`: best fixed `T=144`; 99% mean-`U_G` diapason `130..168`; best tested `T` median `139.5` (p05..p95 `86.9..189.1`).

## Best Fixed Thresholds And Formula Checks

| K | best fixed T | 99% diapason | best tested T median | best tested T std | best formula | formula T | formula fraction |
|---:|---:|---|---:|---:|---|---:|---:|
| 500 | 103 | 84..118 | 95.500 | 30.033 | T_0p05N | 100 | 0.9314 |
| 1000 | 144 | 130..168 | 139.500 | 32.453 | T_0p075N | 150 | 0.9433 |

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
