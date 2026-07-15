# Project Structure

The repository is organized as a source handoff: solver code, reusable helpers,
benchmark entry points, plotting code, and reproducible notes are kept separate.

## Directories

- `algorithms/`: antenna-selection solvers. Public solver functions keep the
  common `solve_*(V, K, sigma, P)` API shape.
- `utils/`: shared helpers. Data generation belongs in `utils/data.py`, generic
  solver evaluation in `utils/evaluation.py`, benchmark solver registries in
  `utils/solver_sets.py`, file-writing helpers in `utils/io.py`, and plotting
  setup in `utils/plotting.py`.
- `experiments/`: active runnable experiment modules only:
  - `experiments.algorithm_comparison`: CDF/runtime comparisons, focused
    requested-gen comparisons, threshold-local studies, and unified local-swap
    grid runs.
  - `experiments.sigma_variation`: sigma sweeps with checkpointing and
    report/plot refresh.
  - `experiments.gurobi_exact`: exact small-case Gurobi comparisons.
- `visualization/`: plotting functions used by active experiments.
- `docs/`: maintained method notes and reproducible commands.
- `results/`: generated benchmark artifacts. This directory is ignored by Git
  and can be deleted without affecting source code or tests.

Do not add one-off experiment entry points in the repository root. New reusable
benchmark modes should be added to the existing `experiments/` modules, and
generated CSV/JSON/plots should go under `results/`.

## Commands

Run experiments through their package modules:

```bash
python -m experiments.algorithm_comparison --help
python -m experiments.sigma_variation --help
python -m experiments.gurobi_exact --help
```

After changing algorithms or benchmark scripts, run:

```bash
python -m unittest discover
```
