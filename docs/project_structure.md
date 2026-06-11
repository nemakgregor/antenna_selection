# Project Structure

The repository separates solver code, reusable helpers, active experiments,
plotting code, and archived experiments.

## Directories

- `algorithms/`: antenna-selection solvers. Public solver functions keep the common `solve_*(V, K, sigma, P)` API shape.
- `utils/`: shared helpers. Put data generation in `utils/data.py`, generic solver evaluation in `utils/evaluation.py`, benchmark solver sets in `utils/solver_sets.py`, file-writing helpers in `utils/io.py`, and plotting setup in `utils/plotting.py`.
- `visualization/`: all plotting functions used by active experiments.
- `experiments/`: active runnable experiment modules only:
  - `experiments.algorithm_comparison`: compare all algorithms with CDF-style `U_G` and runtime plots.
  - `experiments.sigma_variation`: vary `sigma` and plot objective/winner curves.
  - `experiments.gurobi_exact`: small exact Gurobi comparisons.
- `deprecated/experiments/`: old experiments kept for reference, not part of the active workflow.
- `results/`: generated benchmark artifacts.
- `docs/`: workflow notes and reproducible commands.

Run experiments through their package modules:

```bash
venv/bin/python -m experiments.algorithm_comparison --help
venv/bin/python -m experiments.sigma_variation --help
venv/bin/python -m experiments.gurobi_exact --help
```

After changing algorithms or benchmark scripts, run:

```bash
venv/bin/python -m unittest motor_challenge_1205.py
```
