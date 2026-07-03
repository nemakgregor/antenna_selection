# Threshold Full Sweep: rayleigh

- N: 5000
- L: 10
- K values: 1250, 2500
- Samples: 100
- Generator seeds: 42
- Sigma: 1.0

The experiment sweeps every integer `T` from `0` to `K` and evaluates raw `U_G`.

## Answer

- `K=1250`: best fixed `T=250`; 99% mean-`U_G` diapason `208..293`; best tested `T` median `242.5` (p05..p95 `181.8..358.1`).
- `K=2500`: best fixed `T=340`; 99% mean-`U_G` diapason `308..394`; best tested `T` median `336.0` (p05..p95 `255.8..448.4`).

## Best Fixed Thresholds And Formula Checks

| K | best fixed T | 99% diapason | best tested T median | best tested T std | best formula | formula T | formula fraction |
|---:|---:|---|---:|---:|---|---:|---:|
| 1250 | 250 | 208..293 | 242.500 | 55.756 | T_0p05N | 250 | 0.9623 |
| 2500 | 340 | 308..394 | 336.000 | 57.886 | T_0p075N | 375 | 0.9639 |

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
