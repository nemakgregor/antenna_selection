# Small Gurobi Benchmark

Use `venv/bin/python -m experiments.gurobi_exact` for exact small-case checks.
The experiment enumerates feasible subsets, solves the exact Gurobi selection
model for the requested objective, verifies the selected subset with a direct
scan, and compares it with the registered heuristic set from
`utils.solver_sets.SMALL_GUROBI_HEURISTICS`.

## Current Heuristic Set

The exact benchmark compares against:

- `H1`, `H2`, `BackwardTrueGreedy`, `CoutinoSchur-Gen`
- `MISO-EE`, `Pareto-H2`
- `H3-threshold-BF`, `H3-threshold-Int`, `H3-threshold-Gen`
- `Frame-BF`, `Frame-Int`, `Frame-Gen`
- `CapWindow-Gen`
- `CapSubmod-Gen`, `CapSubmodPort-Gen`
- `ThreshDOpt-Gen`, `ThreshWLogdet-Gen`, `ThreshDOptSwap-Gen`
- `H3-Fast`

Multi-objective runs expose the exact solver as:

- `Gurobi-BF`
- `Gurobi-Int`
- `Gurobi-Gen`

## Smoke Run

```bash
venv/bin/python -m experiments.gurobi_exact \
  --N 18 \
  --L 4 \
  --active-frac 0.5 \
  --objective general \
  --seed 42 \
  --out-dir results/gurobi_small_N18_L4_general
```

## Multi-Objective Run

```bash
venv/bin/python -m experiments.gurobi_exact \
  --N 18 \
  --L 4 \
  --active-frac 0.5 \
  --samples 10 \
  --objectives bf int gen \
  --seed 42 \
  --out-dir results/gurobi_small_N18_L4_3obj_10samples
```

Generated CSV, report, and plot files are written under `results/`, which is
ignored by Git.

## Larger Exact Run

```bash
venv/bin/python -m experiments.gurobi_exact \
  --N 22 \
  --L 5 \
  --active-frac 0.5 \
  --samples 10 \
  --objectives bf int gen \
  --seed 42 \
  --max-candidates 3000000 \
  --out-dir results/gurobi_small_N22_L5_3obj_10samples
```

This run is intentionally heavier and should be treated as a benchmark job, not
as a routine regression check. Use `motor_challenge_1205.py` for the fast
post-change regression.
