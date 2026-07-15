# Antenna Selection

This repository contains antenna-selection algorithms and reproducible
benchmarks for the objectives from the original problem statement. The main
large-scale comparison target is the general objective
`U_G = det(V_eq V_eq* + sigma I)`.

## Repository Layout

- `algorithms/`: solver implementations. Public solver functions keep the
  common `solve_*(V, K, sigma, P)` API shape.
- `utils/`: data generation, objective evaluation, solver registries, reporting
  helpers, and atomic result writers.
- `experiments/`: runnable benchmark entry points. Run them as Python modules.
- `visualization/`: plotting code used by the experiment modules.
- `docs/`: maintained method notes and reproducible benchmark commands.
- `results/`: generated benchmark artifacts. This directory is intentionally
  ignored by Git and is not required for a clean source handoff.

## Setup

Create a virtual environment and install the runtime dependencies:

```bash
python -m venv venv
python -m pip install -r requirements.txt
```

On Windows, use `venv\Scripts\python.exe` if your virtual environment uses the
standard Windows layout.

`gurobipy` is optional and is needed only for exact Gurobi benchmarks in
`experiments.gurobi_exact`.

## Active Algorithms

The active implementation surface is `algorithms/__init__.py`. The benchmark
registries live in `utils/solver_sets.py`.

Current solver families include:

- Baselines: `H1`, `H2`, strong/weak `H3`, `H3-Fast`,
  `TrueBackwardGreedy`.
- Threshold methods: `H3` objective variants and the `H3ThresholdT123-Gen`
  portfolio used in focused `U_G` comparisons.
- Frame methods: `Frame-*` and `FrameOnly-*` portfolios.
- Cap-aware methods: `CapWindow-Gen`, `CapWindowFull-Gen`, `CapSubmod-Gen`,
  and `CapSubmodPort-Gen`.
- Local-search refinements: one-swap `U_G` repair variants used by the focused
  requested-gen and unified local-swap comparisons.
- Reference heuristics: `CoutinoSchur-Gen`, `MISO-EE`, `Pareto-H2`, and
  thresholded log-det variants.

## Main Commands

Show command-line options:

```bash
python -m experiments.algorithm_comparison --help
python -m experiments.sigma_variation --help
python -m experiments.gurobi_exact --help
```

Default CDF/runtime comparison:

```bash
python -m experiments.algorithm_comparison \
  --N 1000 --L 2 \
  --samples 100 \
  --generator-seeds 10 42 \
  --off-pcts 25 50 \
  --checkpoint-every 25
```

Focused `U_G` comparison for the requested general-objective solvers:

```bash
python -m experiments.algorithm_comparison \
  --solver-set requested-gen \
  --N 1000 --L 2 \
  --samples 100 \
  --generator-seeds 10 42 \
  --off-pcts 25 50 \
  --checkpoint-every 25
```

Unified local-swap grid:

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

Sigma sweep:

```bash
python -m experiments.sigma_variation \
  --N 1000 --L 4 \
  --off-pcts 25 50 \
  --samples 10 \
  --sigmas 1e-6 1e-4 1e-2 1 100 10000 \
  --checkpoint-every 1
```

Exact small-case benchmark, when Gurobi is installed:

```bash
python -m experiments.gurobi_exact --N 10 --L 2 --samples 5
```

## Validation

Run the regression suite after changing algorithms or benchmark scripts:

```bash
python -m unittest discover
```

The tests cover objective calculation, solver constraints, benchmark registries,
focused comparison modes, archive creation, and unified local-swap output.
