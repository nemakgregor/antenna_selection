# Threshold Exploration

This experiment studies whether the threshold value `T` in the threshold
heuristic can be selected from matrix distribution metrics before trying a
hand-picked threshold list.

The current production threshold solver, `solve_h3`, is unchanged. The
experiment uses copied helper logic in `algorithms/h3_threshold_explore.py`.

## Smoke Run

```bash
venv/bin/python -m experiments.algorithm_comparison \
  --threshold-explore \
  --N 300 \
  --L 2 \
  --samples 20 \
  --off-pcts 25 50 \
  --out-dir results/threshold_exploration_smoke
```

## Full L2 Study

```bash
venv/bin/python -m experiments.algorithm_comparison \
  --threshold-explore \
  --N 1000 \
  --L 2 \
  --samples 200 \
  --off-pcts 25 50 \
  --data-profiles gaussian rayleigh rician nakagami lognormal \
  --out-dir results/threshold_exploration_L2
```

## Optional L4 Check

```bash
venv/bin/python -m experiments.algorithm_comparison \
  --threshold-explore \
  --N 1000 \
  --L 4 \
  --samples 100 \
  --off-pcts 25 50 \
  --data-profiles gaussian rayleigh rician nakagami lognormal thin_tail \
  --out-dir results/threshold_exploration_L4
```

## Full Dense L2 Sweep

This mode studies every integer threshold `T=0..K` for shifted power-window
selection. It is intended to answer where the best tested threshold lies and
whether a simple threshold formula is justified.

Smoke run:

```bash
venv/bin/python -m experiments.algorithm_comparison \
  --threshold-full-sweep \
  --N 100 \
  --L 2 \
  --K-values 50 25 \
  --samples 2 \
  --generator-seeds 42 \
  --data-profiles gaussian rayleigh \
  --out-dir results/smoke_threshold_full_sweep
```

Full run:

```bash
venv/bin/python -m experiments.algorithm_comparison \
  --threshold-full-sweep \
  --N 1000 \
  --L 2 \
  --K-values 500 250 \
  --samples 100 \
  --generator-seeds 42 \
  --data-profiles gaussian rayleigh rician nakagami lognormal thin_tail \
  --sigma 1 \
  --out-dir results/threshold_full_sweep_L2_N1000_K250_500
```

Larger 1000-sample run:

```bash
venv/bin/python -m experiments.algorithm_comparison \
  --threshold-full-sweep \
  --N 1000 \
  --L 2 \
  --K-values 500 250 \
  --samples 1000 \
  --generator-seeds 42 \
  --data-profiles gaussian rayleigh rician nakagami lognormal thin_tail \
  --sigma 1 \
  --out-dir results/threshold_full_sweep_L2_N1000_K250_500_s1000
```

CDF comparison report for formula, best fixed `T`, and strong/weak H3:

```bash
venv/bin/python -m experiments.algorithm_comparison \
  --threshold-rule-cdf \
  --N 1000 \
  --L 2 \
  --K-values 500 250 \
  --samples 1000 \
  --generator-seeds 42 \
  --data-profiles gaussian rayleigh rician nakagami lognormal thin_tail \
  --sigma 1 \
  --out-dir results/threshold_full_sweep_L2_N1000_K250_500_s1000
```

Reports avoid the word "oracle" and use "best tested threshold" instead.

## Outputs

- `threshold_runs.csv`: one row per profile/sample/K/T with objectives,
  selected candidate kind, runtime, and row-power distribution metrics.
- `threshold_summary.csv`: threshold ranking by mean `U_G`, winner rate, and
  normalized gap to the per-case best `T`.
- `threshold_metric_correlations.csv`: Pearson correlations between recorded
  distribution metrics and the per-case best `T`.
- `threshold_exploration_report.md`: compact summary of best thresholds and
  strongest observed metric signals.
- `threshold_cdf_u_g.png`: empirical CDF of raw `U_G`.
- `threshold_cdf_u_g_db.png`: empirical CDF of `10 lg(U_G)`.
- `threshold_cdf_u_g_vs_best.png`: empirical CDF of `U_G(T) / max_T U_G(T)`.

## Interpretation

The primary objective is `U_G`, because it is the capacity-related objective.
If the same distribution metric has a strong correlation with best `T` and the
best mean threshold has a small normalized gap across profiles, that is evidence
for an a priori threshold rule. If correlations are weak or best `T` varies
widely, the result supports keeping a small threshold grid rather than choosing
one threshold from distribution metrics.
