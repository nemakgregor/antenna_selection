# Project Instructions

- Use the local virtual environment: `venv/bin/python`.
- Keep generated benchmark artifacts under `results/`.
- Keep workflow notes and reproducible commands under `docs/`.
- Solver functions in `algorithms/` should keep the common API shape: `solve_*(V, K, sigma, P)`.
- Keep shared helpers in `utils/`; data generation belongs in `utils/data.py`, generic run/evaluation helpers belong in `utils/evaluation.py`, and benchmark solver sets belong in `utils/solver_sets.py`.
- Keep only active experiment implementations in `experiments/`: `algorithm_comparison`, `sigma_variation`, and `gurobi_exact`. Invoke them with `venv/bin/python -m experiments.<module>`. Do not add experiment entrypoints in the repository root.
- For exact small benchmarks, use `venv/bin/python -m experiments.gurobi_exact`; multi-objective runs compare `Gurobi-BF`, `Gurobi-Int`, `Gurobi-Gen` against H1, H2, Coutino, MISO-EE, Pareto-H2, H3-Fast, threshold modes, and the active Frame variants.
- After changing algorithms or benchmark scripts, run `venv/bin/python -m unittest motor_challenge_1205.py`.
