# Algorithms

All public solvers live in `algorithms/` and return a binary vector `x`, where
`x[n] = 1` means antenna `n` is active. The common call shape is:

```python
solve_*(V, K, sigma=1.0, P=1.0, ...)
```

`V` is a complex channel matrix with shape `(N, L)`, `K` is the exact number of
active antennas for the main benchmark solvers, `sigma` is the noise term, and
`P` is the power parameter.

The objective helper is `algorithms.calculate_objectives(V, x, sigma, P)`. It
returns:

```text
U_BF = Tr(M)
U_I  = ||offdiag(M)||_F^2
U_G  = det(M + sigma I)
```

where `M = V_eq V_eq*`.

## Solver Families

- `solve_h1`: row-power baseline.
- `solve_h2`: correlation-aware baseline.
- `solve_h3` and `solve_h3_strong_weak`: threshold-window family.
- `solve_h3_fast`: fast randomized/beam-search heuristic.
- `solve_true_backward_greedy`: backward greedy deletion using the true general
  objective.
- `solve_coutino_schur_greedy`: Schur-complement greedy heuristic.
- `solve_frame_portfolio`: frame-based multi-start portfolio for `bf`, `int`,
  and `gen` targets.
- `solve_cap_window_gen` and `solve_cap_window_full_gen`: cap-aware contiguous
  power-window scans for `U_G`.
- `solve_r2_delta_gen`: analytical R2+Delta rule for `U_G`. It starts from a
  two-moment power-window estimate and then accepts only exact improving pair
  swaps guided by the log-det derivative.
- `solve_cap_submodular_gen` and `solve_cap_submodular_portfolio_gen`: cap-aware
  submodular/log-det candidates for `U_G`.
- `solve_thresholded_logdet_greedy`: thresholded log-det heuristic.
- `solve_miso_energy_greedy` and `solve_pareto_interference_greedy`: reference
  heuristics for energy and interference comparisons.

Local-search helpers are explicit functions:

- `refine_general_1swap`
- `refine_selection_by_ug_swaps`
- `refine_selection_by_ug_swaps_steps`

## Benchmark Registries

Benchmark solver sets are defined in `utils/solver_sets.py`. Use those
registries for experiments instead of constructing ad hoc lists in scripts.
