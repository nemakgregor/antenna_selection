# Threshold Full Sweep: rician

- N: 2000
- L: 10
- K values: 500, 1000
- Samples: 100
- Generator seeds: 42
- Sigma: 1.0

The experiment sweeps every integer `T` from `0` to `K` and evaluates raw `U_G`.

## Answer

- `K=500`: best fixed `T=307`; 99% mean-`U_G` diapason `271..329`; best tested `T` median `309.5` (p05..p95 `149.8..490.1`).
- `K=1000`: best fixed `T=321`; 99% mean-`U_G` diapason `293..378`; best tested `T` median `358.5` (p05..p95 `206.8..570.9`).

## Best Fixed Thresholds And Formula Checks

| K | best fixed T | 99% diapason | best tested T median | best tested T std | best formula | formula T | formula fraction |
|---:|---:|---|---:|---:|---|---:|---:|
| 500 | 307 | 271..329 | 309.500 | 111.466 | T_0p15NL_over_Lp2 | 250 | 0.8315 |
| 1000 | 321 | 293..378 | 358.500 | 106.958 | T_0p15NL_over_Lp2 | 250 | 0.8903 |

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
