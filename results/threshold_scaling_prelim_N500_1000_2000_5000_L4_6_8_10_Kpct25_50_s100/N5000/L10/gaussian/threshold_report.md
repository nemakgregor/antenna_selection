# Threshold Full Sweep: gaussian

- N: 5000
- L: 10
- K values: 1250, 2500
- Samples: 100
- Generator seeds: 42
- Sigma: 1.0

The experiment sweeps every integer `T` from `0` to `K` and evaluates raw `U_G`.

## Answer

- `K=1250`: best fixed `T=247`; 99% mean-`U_G` diapason `205..283`; best tested `T` median `239.0` (p05..p95 `174.8..327.8`).
- `K=2500`: best fixed `T=338`; 99% mean-`U_G` diapason `296..405`; best tested `T` median `332.5` (p05..p95 `247.0..428.6`).

## Best Fixed Thresholds And Formula Checks

| K | best fixed T | 99% diapason | best tested T median | best tested T std | best formula | formula T | formula fraction |
|---:|---:|---|---:|---:|---|---:|---:|
| 1250 | 247 | 205..283 | 239.000 | 52.752 | T_0p05N | 250 | 0.9638 |
| 2500 | 338 | 296..405 | 332.500 | 58.339 | T_0p075N | 375 | 0.9670 |

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
