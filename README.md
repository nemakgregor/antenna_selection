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

## Tested Algorithms

The final comparison set for the customer handoff is the tested
general-objective family around `U_G`. Public solvers are exported from
`algorithms/__init__.py`; benchmark groups are defined in `utils/solver_sets.py`.

| Benchmark name | Implementation | How it is used |
|---|---|---|
| `H3` | `algorithms/h3_strong_weak.py::solve_h3_strong_weak` | Strong/weak power split baseline. |
| `H3Threshold-T0.05N-Gen` | `algorithms/h3_threshold.py::solve_h3` with `t_tests=(round(0.05 * N),)` | Fixed H3 threshold with `T = 0.05N`; wired in `experiments/algorithm_comparison.py` for unified local-search runs. |
| `CapWindow-Gen`, `CapWindowFull-Gen` | `algorithms/cap_window.py` | Cap-aware power-window scans for `U_G`. |
| `CapSubmod-Gen` | `algorithms/cap_submodular.py::solve_cap_submodular_gen` | Cap-aware submodular/log-det candidate for `U_G`. |
| `*-1SwapLS-Gen` and unified `+1swap` runs | `algorithms/local_search.py::refine_general_1swap`, `algorithms/ug_swap_local.py` | Local search applied to the tested seeds. |

`H3ThresholdT123-Gen` in `utils/solver_sets.py` is a broader threshold
portfolio (`0.05N`, `0.15K`, and `0.125NL/(L+2)`). Use the unified local-search
command below when the benchmark must be exactly `T = 0.05N`.

Other implemented heuristics are kept for future research and reference only.
They did not pass the final quality criteria for the customer delivery set, or
were not revalidated under the final benchmark gate. This includes `H1`, `H2`,
`TrueBackwardGreedy`, `CoutinoSchur-Gen`, `MISO-EE`, `Pareto-H2`, `Frame-*`,
`FrameOnly-*`, `H3-Fast`, thresholded log-det variants, `CapSubmodPort-Gen`,
and `R2Delta-Gen`.

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

Run only the main tested registered solvers:

```bash
python -m experiments.algorithm_comparison \
  --solver-set requested-gen \
  --algorithms H3 CapWindow-Gen CapWindowFull-Gen CapSubmod-Gen \
  --N 1000 --L 2 \
  --samples 100 \
  --generator-seeds 10 42 \
  --off-pcts 25 50 \
  --checkpoint-every 25
```

Run the tested local-search variants:

```bash
python -m experiments.algorithm_comparison \
  --solver-set requested-gen \
  --algorithms H3-1SwapLS-Gen CapWindow-1SwapLS-Gen \
    CapWindowFull-1SwapLS-Gen CapSubmod-1SwapLS-Gen \
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

This unified mode is the benchmark path for `H3Threshold-T0.05N-Gen`,
`CapWindowFull-Gen`, `CapSubmod-Gen`, `strong_weak`, and their local-search
`+0swap`/`+1swap` variants.

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
