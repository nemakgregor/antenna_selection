# Antenna Selection

This project compares antenna-selection algorithms for the objectives defined in
the original problem statement, with the main experimental focus on the general
objective `U_G`.

## Active Layout

- `algorithms/`: solver implementations. Public solver functions keep the
  `solve_*(V, K, sigma, P)` API shape.
- `experiments/`: active runnable experiments only.
  - `algorithm_comparison.py`: compares all registered algorithms, CDF-style,
    with primary plots for `U_G` and runtime.
  - `sigma_variation.py`: varies `sigma` and writes objective/winner plots.
  - `gurobi_exact.py`: exact small-case Gurobi benchmark and multi-objective
    Gurobi comparisons.
- `visualization/`: plotting functions used by active experiments.
- `utils/`: data generation, solver evaluation, solver sets, IO helpers, and
  matplotlib backend setup.
- `docs/`: reproducible commands and notes.
- `results/`: generated benchmark artifacts. This directory is ignored by Git.

## Commands

Use the project virtual environment:

```bash
venv/bin/python -m experiments.algorithm_comparison --help
venv/bin/python -m experiments.sigma_variation --help
venv/bin/python -m experiments.gurobi_exact --help
```

Required regression check after algorithm or benchmark changes:

```bash
venv/bin/python -m unittest motor_challenge_1205.py
```
