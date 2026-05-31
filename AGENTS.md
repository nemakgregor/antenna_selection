# Project Instructions

- Use the local virtual environment: `venv/bin/python`.
- Keep generated benchmark artifacts under `results/`.
- Keep workflow notes and reproducible commands under `docs/`.
- Solver functions in `algorithms/` should keep the common API shape: `solve_*(V, K, sigma, P)`.
- For exact small benchmarks, use `small_gurobi_optimum.py`; multi-objective runs compare `Gurobi-BF`, `Gurobi-Int`, `Gurobi-Gen` against H1, H2, Coutino, MISO-EE, Pareto-H2, H3-Fast, threshold modes, and the active Frame variants.
- After changing algorithms or benchmark scripts, run `venv/bin/python -m unittest motor_challenge_1205.py`.
