# Threshold Full Sweep: rayleigh

- N: 1000
- L: 10
- K values: 250, 500
- Samples: 100
- Generator seeds: 42
- Sigma: 1.0

The experiment sweeps every integer `T` from `0` to `K` and evaluates raw `U_G`.

## Answer

- `K=250`: best fixed `T=53`; 99% mean-`U_G` diapason `50..58`; best tested `T` median `46.0` (p05..p95 `23.9..74.0`).
- `K=500`: best fixed `T=78`; 99% mean-`U_G` diapason `65..88`; best tested `T` median `71.5` (p05..p95 `46.0..110.0`).

## Best Fixed Thresholds And Formula Checks

| K | best fixed T | 99% diapason | best tested T median | best tested T std | best formula | formula T | formula fraction |
|---:|---:|---|---:|---:|---|---:|---:|
| 250 | 53 | 50..58 | 46.000 | 16.631 | T_0p05N | 50 | 0.9039 |
| 500 | 78 | 65..88 | 71.500 | 21.047 | T_0p075N | 75 | 0.9189 |

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
