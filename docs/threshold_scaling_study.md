# Preliminary Threshold Scaling Study

> Historical K semantics note: this report uses active-K semantics. Here `K` is the number of selected/kept antennas, not the number turned off. A `25% active` or `K=0.25N` case means `75% off`, not the real `25% off` task. For real off-percent experiments, `25% off => K_active=0.75N` and `50% off => K_active=0.50N`.

This experiment studies whether the threshold-window parameter `T` scales better
with `N`, active `K`, stream count `L`, or row-power distribution metrics.
Production solvers are unchanged.

## Smoke Run

```bash
venv/bin/python -m experiments.algorithm_comparison \
  --threshold-scaling-study \
  --N-values 100 200 \
  --L-values 4 6 \
  --K-pcts 25 50 \
  --samples 2 \
  --generator-seeds 42 \
  --data-profiles gaussian twdp \
  --sigma 1 \
  --out-dir results/smoke_threshold_scaling
```

## Preliminary Run

```bash
venv/bin/python -m experiments.algorithm_comparison \
  --threshold-scaling-study \
  --N-values 500 1000 2000 5000 \
  --L-values 4 6 8 10 \
  --K-pcts 25 50 \
  --samples 100 \
  --generator-seeds 42 \
  --data-profiles gaussian rayleigh rician nakagami lognormal thin_tail twdp \
  --sigma 1 \
  --out-dir results/threshold_scaling_prelim_N500_1000_2000_5000_L4_6_8_10_Kpct25_50_s100
```

## Outputs

Each `(N, L, profile)` shard writes compressed raw threshold rows under:

```text
results/.../N{N}/L{L}/{profile}/threshold_runs.csv.gz
```

Each shard also writes `threshold_strong_weak_runs.csv`, which evaluates the
strong/weak H3 rule against that sample's best tested threshold result.

Root outputs include:

- `all_best_thresholds.csv`
- `all_threshold_best_t_stats.csv`
- `all_distribution_comparison.csv`
- `all_formula_comparison.csv`
- `all_scaling_formula_summary.csv`
- `all_scaling_metric_correlations.csv`
- `all_formula_selected_runs.csv`
- `all_strong_weak_runs.csv`
- `all_strong_weak_summary.csv`
- `all_scaling_report.md`
- combined CDF and formula/scale plots

## Refresh Combined Outputs

Use this when shards already exist and only the combined report/plots need to be
rebuilt, or when adding derived comparisons such as strong/weak H3:

```bash
venv/bin/python -m experiments.algorithm_comparison \
  --threshold-scaling-study \
  --plot-only \
  --samples 100 \
  --generator-seeds 42 \
  --K-pcts 25 50 \
  --sigma 1 \
  --out-dir results/threshold_scaling_prelim_N500_1000_2000_5000_L4_6_8_10_Kpct25_50_s100
```
