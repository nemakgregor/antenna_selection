# Threshold Exploration

- N: 1000
- L: 2
- Samples per seed/profile: 200
- Generator seeds: 10, 42
- Data profiles: gaussian, rayleigh, rician, nakagami, lognormal
- Off percentages: 25, 50
- Sigma: 1.0
- P: 1.0
- Threshold target objective: gen

Primary metric is `U_G`. `u_g_vs_best_T_mean` is the mean fraction of the per-case best threshold objective.

## Best Mean Thresholds

| profile | off % | K | best mean T | source | mean U_G / best_T | winner rate | mean gap % |
|---|---:|---:|---:|---|---:|---:|---:|
| gaussian | 25.000 | 750 | 31 | gap_top3 | 1.00000 | 1.000 | 0.000 |
| gaussian | 50.000 | 500 | 31 | gap_top3 | 1.00000 | 1.000 | 0.000 |
| lognormal | 25.000 | 750 | 65 | mean_plus_1std | 0.98755 | 0.538 | 1.245 |
| lognormal | 50.000 | 500 | 65 | mean_plus_1std | 0.97658 | 0.385 | 2.342 |
| nakagami | 25.000 | 750 | 65 | gap_top2 | 1.00000 | 1.000 | 0.000 |
| nakagami | 50.000 | 500 | 70 | mean_plus_1.5std | 0.99982 | 0.667 | 0.018 |
| rayleigh | 25.000 | 750 | 64 | mean_plus_1.5std | 1.00000 | 1.000 | 0.000 |
| rayleigh | 50.000 | 500 | 34 | gap_top3 | 1.00000 | 1.000 | 0.000 |
| rician | 25.000 | 750 | 65 | mean_plus_1.5std | 1.00000 | 1.000 | 0.000 |
| rician | 50.000 | 500 | 65 | mean_plus_1.5std | 1.00000 | 1.000 | 0.000 |

## Distribution-Metric Signal

The table lists the strongest absolute Pearson correlation between a recorded distribution metric and the per-case best `T`.

| profile | off % | K | strongest metric | corr | interpretation |
|---|---:|---:|---|---:|---|
| gaussian | 25.000 | 750 | tail_mass_p95 | 0.298 | weak/no single-metric signal |
| gaussian | 50.000 | 500 | tail_mass_p95 | 0.246 | weak/no single-metric signal |
| lognormal | 25.000 | 750 | tail_mass_p95 | 0.207 | weak/no single-metric signal |
| lognormal | 50.000 | 500 | tail_mass_p95 | 0.177 | weak/no single-metric signal |
| nakagami | 25.000 | 750 | tail_mass_p95 | 0.307 | weak/no single-metric signal |
| nakagami | 50.000 | 500 | tail_mass_p95 | 0.309 | weak/no single-metric signal |
| rayleigh | 25.000 | 750 | tail_mass_p95 | 0.248 | weak/no single-metric signal |
| rayleigh | 50.000 | 500 | tail_mass_p95 | 0.247 | weak/no single-metric signal |
| rician | 25.000 | 750 | row_power_p95_p50 | 0.190 | weak/no single-metric signal |
| rician | 50.000 | 500 | tail_mass_p95 | 0.168 | weak/no single-metric signal |

## Artifacts

- `threshold_runs.csv`: one row per `(profile, seed, sample, K, T)`.
- `threshold_summary.csv`: threshold rankings by mean `U_G` and normalized gap.
- `threshold_metric_correlations.csv`: simple distribution-metric correlations with best `T`.
- `threshold_cdf_u_g.png`: CDF of raw `U_G` for top thresholds.
- `threshold_cdf_u_g_db.png`: CDF of `10 lg(U_G)` for top thresholds.
- `threshold_cdf_u_g_vs_best.png`: CDF of `U_G(T) / max_T U_G(T)`.