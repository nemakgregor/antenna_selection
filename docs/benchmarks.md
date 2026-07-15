# Benchmarks

Generated CSV, JSON, archive, report, and plot files should be written under
`results/`. The directory is ignored by Git and is not part of the source
handoff.

## Regression Tests

```bash
python -m unittest discover
```

## Default Comparison

```bash
python -m experiments.algorithm_comparison \
  --N 1000 --L 2 \
  --samples 100 \
  --generator-seeds 10 42 \
  --off-pcts 25 50 \
  --checkpoint-every 25
```

## Focused General-Objective Comparison

```bash
python -m experiments.algorithm_comparison \
  --solver-set requested-gen \
  --N 1000 --L 2 \
  --samples 100 \
  --generator-seeds 10 42 \
  --off-pcts 25 50 \
  --checkpoint-every 25
```

`requested-gen` includes `R2Delta-Gen`, `CapWindow*`, `CapSubmod*`, H3
threshold variants, and their active one-swap refinements.

## Unified Local-Swap Grid

```bash
python -m experiments.algorithm_comparison \
  --unified-local-swap-comparison \
  --N 1000 \
  --L-values 2 4 6 8 10 \
  --sigmas 1 10 100 1000 10000 \
  --off-pcts 25 50 \
  --samples 100 \
  --generator-seeds 42 \
  --data-profiles gaussian rayleigh rician nakagami lognormal thin_tail \
  --compact-runs \
  --checkpoint-every 1000 \
  --resume
```

## Sigma Sweep

```bash
python -m experiments.sigma_variation \
  --N 1000 --L 4 \
  --off-pcts 25 50 \
  --samples 10 \
  --sigmas 1e-6 1e-4 1e-2 1 100 10000 \
  --checkpoint-every 1
```

## Exact Small Cases

This command requires `gurobipy`:

```bash
python -m experiments.gurobi_exact --N 10 --L 2 --samples 5
```
