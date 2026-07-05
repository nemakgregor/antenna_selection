# Threshold Heuristic Report

> Historical K semantics note: this report uses active-K semantics. Here `K` is the number of selected/kept antennas, not the number turned off. A `25% active` or `K=0.25N` case means `75% off`, not the real `25% off` task. For real off-percent experiments, `25% off => K_active=0.75N` and `50% off => K_active=0.50N`.

## Repository And Environment

The repository follows the active layout described in `README.md` and
`docs/project_structure.md`:

- solver implementations live in `algorithms/`;
- active experiment entrypoints live in `experiments/`;
- shared generation/evaluation helpers live in `utils/`;
- generated benchmark artifacts live in `results/`;
- reproducible notes live in `docs/`.

No `requirements.txt` or `pyproject.toml` was present. The active Python
dependency surface is:

- `numpy`: solver math, data generation, objective evaluation;
- `pandas`: experiment summaries and report tables;
- `matplotlib`: plot generation;
- optional `gurobipy`: needed only when actually solving exact Gurobi runs in
  `experiments.gurobi_exact`.

Local setup used:

```bash
python3 -m venv venv
source venv/bin/activate
python -m pip install numpy pandas matplotlib
```

In this workspace the venv uses Python 3.12.3. Commands were run through
`venv/bin/python`, which is equivalent to using the activated environment in
non-interactive shells.

Regression check:

```bash
venv/bin/python -m unittest motor_challenge_1205.py
```

Result: 14 tests passed.

## Problem Context

The PDF problem statement asks to strike out rows of a complex matrix `V` so the
effective selected matrix remains close to an orthonormal basis under a
per-antenna power cap. The implementation evaluates the same three objective
families:

- `U_BF`: maximize beamforming gain;
- `U_I`: minimize off-diagonal interference;
- `U_G`: maximize `det(V_eq V_eq* + sigma I)`.

The shared implementation is `calculate_objectives` in `algorithms/common.py`.
It rescales selected rows by the active maximum row power before calculating all
three metrics, so candidate scoring honors the per-antenna power constraint.

## Where The Threshold Heuristic Lives

The threshold heuristic is `solve_h3` in `algorithms/h3_threshold.py`. It is
registered three times through `utils/solver_sets.py`:

- `H3-threshold-BF` / `S-threshold-BF`: `target_obj="bf"`;
- `H3-threshold-Int` / `S-threshold-Int`: `target_obj="int"`;
- `H3-threshold-Gen` / `S-threshold-Gen`: `target_obj="gen"`.

Frame portfolio solvers also use threshold modes as external starts, especially
the interference start for `Frame-Int`.

## How It Works

The solver is deterministic and has two phases: generate a small family of
thresholded candidate sets, then score those candidates using the exact project
objective function.

1. Validate `V`, `K`, and `target_obj`.
   `target_obj` must be `bf`, `int`, or `gen`. Edge cases return all-zero or
   all-one selections for `K=0` and `K=N`.

2. Rank rows by row power.
   It computes `p_n = sum(abs(V[n, :]) ** 2)` and sorts rows in descending
   power. This ranking defines the thresholds.

3. Precompute per-row off-diagonal interference tensors.
   For each row it builds `row_interference[n] = conj(V[n]) outer V[n]`, then
   zeros the diagonal. The greedy interference routine works on sums of these
   off-diagonal matrices.

4. Add a pure power candidate.
   The first candidate is simply the top `K` rows by row power.

5. Add a full-pool phase-nulling candidate.
   If `N > K`, it starts from all rows and repeatedly deletes one active row
   until exactly `K` remain. At each deletion it chooses the row whose removal
   leaves the smallest squared norm of the residual off-diagonal sum.

6. Sweep hard-coded thresholds.
   The threshold list is `[1, 2, 5, 10, 25, 50, N // 10]`, filtered to
   `0 < T <= N - K` and deduplicated.

   For each threshold `T`:

   - BF and Gen modes add a shifted power window:
     choose sorted rows `idx_desc[T : T + K]`. This deliberately skips the top
     `T` strongest rows in case a few high-power antennas hurt balance.
   - Int and Gen modes add a tail-pool phase-nulling candidate:
     start from `idx_desc[T:]`, lock the strongest remaining row as an anchor,
     then greedily delete down to `K`.
   - Gen mode also adds a small buffered phase-nulling candidate:
     start from `idx_desc[T : T + K + buffer]`, with `buffer <= 30`, lock the
     strongest row in that window, and delete only the buffer rows.

7. Score candidates.
   Every generated candidate is evaluated with `calculate_objectives`.
   `target_obj="bf"` maximizes `U_BF`; `target_obj="int"` maximizes `-U_I`;
   `target_obj="gen"` maximizes `U_G`.

The result always contains exactly `K` active rows for ordinary `0 < K < N`
cases, even though the formal constraint is `active_count <= K`.

## Mode Behavior

`target_obj="bf"` is the cheapest mode. It mostly searches power-ranked windows:
top `K`, full-pool phase-nulling, and shifted top-`K` windows. It often gives
excellent BF gain and can also be very good for `U_G` at large `sigma`, where
the determinant objective increasingly tracks BF gain.

`target_obj="int"` focuses on off-diagonal cancellation. It can produce much
lower interference than BF mode, but it often gives up too much BF gain for the
general objective. The locked anchor prevents the deletion loop from drifting
into an overly weak or degenerate set.

`target_obj="gen"` is the hybrid. It evaluates power windows, full/tail
phase-nulling, and buffered candidates with the actual determinant objective.
On the existing large sigma sweep, it is the main `U_G` winner at small and
moderate `sigma`.

## Good Sides

- Simple and deterministic: no random state, no optimizer state, and no fragile
  convergence criteria.
- Scores with the real project objective after candidate generation, including
  the per-antenna power rescaling.
- Gives a useful bridge between the two extremes from the PDF: power-heavy
  BF candidates and interference-cancelling candidates.
- Strong empirical `U_G` baseline. In the existing `N=1000, L=4` sigma sweep,
  `H3-threshold-Gen` leads mean `U_G` for `K=250` across the full tested sigma
  grid and for `K=500` until very large `sigma`, where `H3-threshold-BF` takes
  over.
- Useful as a seed generator. The Frame portfolio code uses threshold starts,
  then local refinement closes much of the remaining exact-optimality gap.
- Fast enough for benchmark-scale use at `N=1000`, especially BF mode.

## Pitfalls

- The threshold grid is hard-coded. If the best set requires skipping, say, 15
  or 80 top-power rows, the solver only finds it indirectly through a nearby
  candidate or not at all.
- It has no local swap refinement. A near-miss threshold window cannot repair a
  few bad inclusions/exclusions unless another generated candidate already
  captures them.
- It always returns exactly `K` active rows. That is usually sensible for BF and
  `U_G`, but pure interference minimization under the formal `<= K` constraint
  could prefer fewer active antennas.
- Candidate generation is not the same as objective optimization. The
  phase-nulling deletion proxy minimizes off-diagonal residuals before final
  scoring; it does not directly optimize determinant, BF, or the rescaled
  objective at every deletion.
- `target_obj="int"` can be too conservative for `U_G`. Existing exact small
  benchmarks show it can have low interference but large general-objective gap.
- Runtime grows quadratically in the size of phase-nulling pools. The precompute
  step stores `O(N L^2)` row tensors, and each greedy deletion scans active rows
  with `L x L` residuals. Existing `N=1000, L=4` sigma-sweep averages show:
  `H3-threshold-BF` about 0.04-0.05 s per solve, while `H3-threshold-Int` and
  `H3-threshold-Gen` are about 0.39-0.54 s per solve.

## Evidence From Existing Results

From `results/sigma_sweep_N1000_L4_Kpct25_50_30samples/sigma_sweep_report.md`:

- For `K=250`, `H3-threshold-Gen` is the mean `U_G` leader from
  `sigma=0.001` through `sigma=100000`.
- For `K=500`, `H3-threshold-Gen` is the mean `U_G` leader from
  `sigma=0.001` through `sigma=3000`; `H3-threshold-BF` leads from
  `sigma=10000` through `sigma=100000`.
- `H3-threshold-BF` has the best mean BF gain among the compared heuristics in
  both `K` cases.
- `H3-Fast`, not threshold mode, has the lowest mean interference in that
  large sweep.

From exact small reports:

- In the `N=18, L=4, K=9` general-objective run, the best threshold variants
  were about 1.41 percent below exact optimum.
- In the `N=22, L=5, K=11`, 3-sample exact comparison with Frame variants,
  `H3-threshold-BF` and `H3-threshold-Gen` were about 2.17 percent below exact
  `U_G` on average, while Frame refinement matched exact `U_G`.

## Practical Guidance

- Use `H3-threshold-Gen` as the standalone default when the goal is `U_G` at
  normal or low/moderate `sigma`.
- Check `H3-threshold-BF` at large `sigma`, because `U_G` becomes increasingly
  BF-driven.
- Do not use `H3-threshold-Int` as a general-objective default. Use it when the
  explicit target is interference or as a start for a refinement method.
- For stronger final answers, feed threshold candidates into Frame/refinement
  rather than treating threshold selection as the last step.
- If improving this heuristic, the most direct upgrade is adaptive thresholds
  plus a bounded swap polish over the best few candidates.
