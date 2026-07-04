# Small Gurobi Benchmark

> Historical K semantics note: this report uses active-K semantics. Here `K` is the number of selected/kept antennas, not the number turned off. A `25% active` or `K=0.25N` case means `75% off`, not the real `25% off` task. For real off-percent experiments, `25% off => K_active=0.75N` and `50% off => K_active=0.50N`.

This repository uses `venv/bin/python -m experiments.gurobi_exact` for exact small-case comparisons.
The script enumerates feasible subsets, asks Gurobi to select the exact optimum
for the requested objective, verifies the answer with a direct scan, then compares
the optimum with all heuristic algorithms.

## Algorithms Compared

- `H1`
- `H2`
- `Coutino`
- `MISO-EE`
- `Pareto-H2`
- `H3-threshold-BF`: `solve_h3(..., target_obj="bf")`
- `H3-threshold-Int`: `solve_h3(..., target_obj="int")`
- `H3-threshold-Gen`: `solve_h3(..., target_obj="gen")`
- `H3-Fast`

Multi-objective runs expose the exact solver as three separate algorithms:

- `Gurobi-BF`
- `Gurobi-Int`
- `Gurobi-Gen`

## Reproduce Current Run

```bash
venv/bin/python -m experiments.gurobi_exact \
  --N 18 \
  --L 4 \
  --active-frac 0.5 \
  --objective general \
  --seed 42 \
  --out-dir results/gurobi_small_N18_L4_general
```

Generated files:

- `results/gurobi_small_N18_L4_general/gurobi_small_report.md`
- `results/gurobi_small_N18_L4_general/gurobi_small_comparison.csv`
- `results/gurobi_small_N18_L4_general/gurobi_small_optimum.csv`
- `results/gurobi_small_N18_L4_general/gurobi_small_comparison.png`

## Current Result Summary

Parameters: `N=18`, `L=4`, `K<=9`, `objective=general`, `seed=42`.

- Enumerated subsets: `154394`
- Gurobi optimum subset: `3 7 9 11 12 13 15 16 17`
- Gurobi optimum metrics: `BF=15.1666`, `Interference=19.9859`, `U_G=307.0654`
- Timing: enumeration `5.602428s`, Gurobi model build/solve `1.364855s`, direct scan verification `0.056316s`

For this run, `Coutino` and `MISO-EE` matched the Gurobi `U_G` optimum. The
best `h3_threshold` variants by `U_G` were `H3-threshold-BF` and
`H3-threshold-Gen`, each with a `1.41%` objective gap to the exact optimum.

## Multi-Objective Averaged Run

A single Gurobi optimum is exact only for the objective it optimizes. To compare
against the three `h3_threshold.py` target modes, run three separate exact
Gurobi objectives on each random sample:

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

Generated files:

- `results/gurobi_small_N18_L4_3obj_10samples/gurobi_multi_objective_report.md`
- `results/gurobi_small_N18_L4_3obj_10samples/gurobi_multi_objective_summary.csv`
- `results/gurobi_small_N18_L4_3obj_10samples/gurobi_multi_objective_runs.csv`
- `results/gurobi_small_N18_L4_3obj_10samples/gurobi_multi_objective_optima.csv`
- `results/gurobi_small_N18_L4_3obj_10samples/gurobi_multi_objective_wins.csv`
- `results/gurobi_small_N18_L4_3obj_10samples/gurobi_multi_objective_summary.png`

For `N=18`, `L=4`, `K<=9`, seeds `42..51`, the averaged exact target means are:

- BF objective: `Gurobi-BF` has `u_bf=13.2846`.
- Interference objective: `Gurobi-Int` has `u_i=0.1439`.
- General objective: `Gurobi-Gen` has `u_g=126.1325`.

The summary plot includes all three Gurobi variants in each objective panel, so
it shows the cross-objective tradeoff between exact BF, exact interference, and
exact general-objective solutions.

`winner_rate` means the algorithm is in the tied winner set for a sample.
`split win share` divides one sample's win across tied winners.

## Larger One-Hour Run

The next heavier exact case increases both antenna count and layer count while
staying below the default laptop-scale exact-search ceiling:

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

Parameters: `N=22`, `L=5`, `K<=11`, seeds `42..51`.

- Enumerated subsets per seed: `2440759`.
- BF objective: `Gurobi-BF` has mean `u_bf=16.4390` and `winner_rate=1.00`.
- Interference objective: `Gurobi-Int` has mean `u_i=0.3277` and `winner_rate=1.00`.
- General objective: `Gurobi-Gen` has mean `u_g=376.1847` and `winner_rate=1.00`.
- Mean timing per Gurobi variant: enumeration `95.07s`, model build/solve about `21.6..21.9s`, direct scan about `0.9..1.0s`.

Generated files are under `results/gurobi_small_N22_L5_3obj_10samples/`.

## Sigma Sweep Without Gurobi

For larger non-exact experiments, use `venv/bin/python -m experiments.sigma_variation`. It compares all
non-Gurobi algorithms. The current percent-based run uses `--K-pcts`, so `25`
means `K=round(0.25*N)`, not absolute `K=25`.

```bash
venv/bin/python -m experiments.sigma_variation \
  --N 1000 \
  --L 4 \
  --K-pcts 25 50 \
  --samples 30 \
  --sigmas 0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000 3000 10000 30000 100000 \
  --seed 42 \
  --checkpoint-every 1 \
  --summary-every 1 \
  --plot-every 0 \
  --out-dir results/sigma_sweep_N1000_L4_Kpct25_50_30samples
```

Generated files:

- `results/sigma_sweep_N1000_L4_Kpct25_50_30samples/sigma_sweep_report.md`
- `results/sigma_sweep_N1000_L4_Kpct25_50_30samples/sigma_sweep_summary.csv`
- `results/sigma_sweep_N1000_L4_Kpct25_50_30samples/sigma_sweep_winners.csv`
- `results/sigma_sweep_N1000_L4_Kpct25_50_30samples/sigma_sweep_mean_leaders.csv`
- `results/sigma_sweep_N1000_L4_Kpct25_50_30samples/sigma_sweep_selection_stability.csv`
- `results/sigma_sweep_N1000_L4_Kpct25_50_30samples/sigma_sweep_K250.png`
- `results/sigma_sweep_N1000_L4_Kpct25_50_30samples/sigma_winners_K250.png`
- `results/sigma_sweep_N1000_L4_Kpct25_50_30samples/sigma_sweep_K500.png`
- `results/sigma_sweep_N1000_L4_Kpct25_50_30samples/sigma_winners_K500.png`

Observed mean `U_G` leader segments:

- `K<=250` (`25% active`): `H3-threshold-Gen` from `sigma=0.001` through
  `100000`.
- `K<=500` (`50% active`): `H3-threshold-Gen` from `0.001..3000`, then
  `H3-threshold-BF` from `10000..100000`.
