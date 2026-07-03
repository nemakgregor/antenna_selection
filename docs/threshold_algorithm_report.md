# Threshold-Window Algorithm Report

Checked on 2026-07-02.

This report describes the threshold-window antenna-selection algorithm used in
the current experiments, including pseudocode, objective equations, asymptotic
cost, and best/worst input regimes.

## Direct Answer

The threshold approach sorts antennas by row power and selects a contiguous
window of exactly `K` antennas after skipping `T` strongest rows:

```text
active antennas = sorted_by_row_power_desc[T : T + K]
```

The best single tested formula in the current scaling run is:

```text
T = round(0.15 K)
```

The best 3-threshold portfolio found from the existing scaling result is:

```text
T1 = round(0.05 N)
T2 = round(0.15 K)
T3 = round(0.125 N L / (L + 2))
```

For the portfolio, evaluate all three windows and return the one with largest
`U_G`. This is still cheap because it uses only three candidates.

Observed on `results/threshold_scaling_prelim_N500_1000_2000_5000_L4_6_8_10_Kpct25_50_s100`:

| rule | mean fraction of best tested `U_G` |
|---|---:|
| best single: `T = 0.15K` | 0.8927 |
| best 2-threshold portfolio: `0.05N + 0.15K` | 0.9234 |
| best 3-threshold portfolio: `0.05N + 0.125NL/(L+2) + 0.15K` | 0.9366 |

## Objective Equations

Let `V` be the `N x L` complex channel matrix and `S` be the active antenna set.
The row power is:

```text
p_n = sum_l |V[n,l]|^2
```

The selected-set Gram matrix is:

```text
G_S = V_S^H V_S
```

The code scales by the maximum selected row power:

```text
z^2 = P / max_{n in S} p_n
```

The equivalent matrix in `calculate_objectives` is:

```text
V_eq = z G_S
V_eq V_eq^H = z^2 G_S G_S^H
```

The three project objectives are:

```text
U_BF = trace(V_eq V_eq^H)

U_I = sum_ij |(V_eq V_eq^H)_ij|^2
      - sum_i |(V_eq V_eq^H)_ii|^2

U_G = det(V_eq V_eq^H + sigma I_L)
```

For a Hermitian Gram matrix `G_S`, this is computed in the faster helper as:

```text
G2 = G_S G_S^H
U_BF = z^2 trace(G2)
U_I  = z^4 offdiag_energy(G2)
U_G  = det(z^2 G2 + sigma I_L)
```

Numerical check against `calculate_objectives` on random threshold windows:

```text
overall max absolute error = 2.98e-08
```

The determinant objective is capacity-related because MIMO capacity formulas
commonly use a log-determinant form such as:

```text
log2 det(I + rho H Q H^H)
```

Our `U_G` is not itself the final bit-rate formula because it omits the log and
uses the project-specific scaling `z`, but maximizing `U_G` still favors large
and balanced eigenmodes of the selected equivalent channel.

## Pseudocode: Single Formula Threshold

```text
Input:
    V      complex N x L channel matrix
    K      active antenna count
    sigma  noise parameter
    P      power parameter
    formula, for example T = round(0.15 K)

Output:
    x      binary length-N selection vector with exactly K active antennas
    U_G    selected objective value

Algorithm:
    N, L <- shape(V)

    for n in 0..N-1:
        p[n] <- sum_l |V[n,l]|^2

    order <- argsort(p, descending=True)

    T <- round(formula(N, L, K))
    T <- clip(T, 0, min(K, N-K))

    active_idx <- order[T : T + K]

    x <- zeros(N)
    x[active_idx] <- 1

    G <- V[active_idx]^H V[active_idx]
    z2 <- P / max_{n in active_idx} p[n]
    U_G <- det(z2 * G * G^H + sigma * I_L)

    return x, U_G
```

## Pseudocode: 3-Threshold Portfolio

```text
Input:
    V, K, sigma, P

Candidate thresholds:
    T1 = round(0.05 N)
    T2 = round(0.15 K)
    T3 = round(0.125 N L / (L + 2))

Algorithm:
    compute row powers p[n]
    order <- argsort(p, descending=True)

    best_score <- -infinity
    best_x <- None

    for T in [T1, T2, T3]:
        T <- clip(T, 0, min(K, N-K))
        active_idx <- order[T : T + K]
        score <- U_G(active_idx)

        if score > best_score:
            best_score <- score
            best_x <- selection vector for active_idx

    return best_x, best_score
```

This is not a dense search. It evaluates only three shifted windows.

## Pseudocode: Dense Experimental Sweep

The dense sweep is for analysis only:

```text
Input:
    V, K, sigma, P

compute p[n]
order <- argsort(p, descending=True)
V_sorted <- V[order]

prefix[0] <- zero L x L matrix
for i in 0..N-1:
    prefix[i+1] <- prefix[i] + V_sorted[i]^H V_sorted[i]

for T in 0..K:
    G_T <- prefix[T+K] - prefix[T]
    z2 <- P / p[order[T]]
    U_G(T) <- det(z2 * G_T * G_T^H + sigma I_L)

return all U_G(T)
```

This is how the experiments find the "best tested threshold" for comparison.

## Asymptotic Complexity

Let:

- `N`: number of antennas / rows.
- `L`: number of streams / columns.
- `K`: active antennas selected.
- `q`: number of candidate thresholds evaluated.

For a direct formula/portfolio implementation:

```text
row powers:             O(N L)
sorting row powers:     O(N log N)
per candidate Gram:     O(K L^2)
per candidate objective O(L^3)
```

Total time:

```text
O(N L + N log N + q (K L^2 + L^3))
```

For the single formula, `q = 1`. For the recommended portfolio, `q = 3`.
When `L` is small and fixed, the practical scaling is close to:

```text
O(N log N + q K)
```

Auxiliary space for the direct version:

```text
O(N + L^2)
```

If returning the binary vector `x`, output space is also `O(N)`.

For the dense experimental prefix-Gram sweep:

```text
row powers and sort:    O(N L + N log N)
prefix Gram build:      O(N L^2)
per threshold scoring:  O(L^3)
```

Total dense-sweep time:

```text
O(N L + N log N + N L^2 + q L^3)
```

Dense-sweep auxiliary space:

```text
O(N L^2 + N)
```

The current dense sweep uses `q = K + 1` because it evaluates every `T = 0..K`.
This is for experiments, not for the final heuristic.

For the exact brute-force baseline added for small `N`, every exact-`K` subset
is evaluated:

```text
number of subsets = C(N, K)
```

Each subset needs one selected Gram sum and one objective evaluation:

```text
per subset: O(K L^2 + L^3)
```

So the exact baseline costs:

```text
O(C(N,K) (K L^2 + L^3))
```

and uses:

```text
O(N L^2 + L^2 + N)
```

auxiliary space in the implementation, because row Gram contributions are
precomputed once and reused across subsets. This is intentionally only a
small-`N` validation tool.

## When the Threshold Approach Is Optimal

For exact-`K` selection, the threshold algorithm is exactly optimal if an
optimal set is a contiguous block in descending row-power order:

```text
S* = { order[T], order[T+1], ..., order[T+K-1] }
```

for one of the tested threshold values.

The single formula is optimal when the optimal block starts at the formula
threshold, for example:

```text
T* = round(0.15 K)
```

The 3-threshold portfolio is optimal when the optimal block starts at one of:

```text
round(0.05N)
round(0.15K)
round(0.125NL/(L+2))
```

## Approach View: When Sort + Shift Can Be Best Or Worst

The threshold approach restricts the search space from all exact-`K` subsets:

```text
C(N, K)
```

to the shifted row-power windows:

```text
{ order[T : T+K] }
```

For the dense experiment in this repo, `T = 0..K`, so there are only `K+1`
tested windows. A formula or portfolio tests even fewer windows.

This means the approach can achieve the exact global optimum only when at
least one globally optimal subset is representable as one of those contiguous
windows. In symbols, for some tested `T`:

```text
S_opt = { order[T], order[T+1], ..., order[T+K-1] }
```

If the exact optimum is not contiguous in row-power rank, no choice of a single
threshold `T` can reproduce it. In that case, the best possible threshold result
is the best representable window, not the true subset optimum.

The new exact study measures this directly by saving:

```text
U_G(best tested threshold) / U_G(exact brute-force optimum)
```

and a boolean flag:

```text
exact_is_threshold_window
```

So the experiment separates two questions:

1. Is the threshold-window search space itself enough?
2. If yes, does the proposed formula choose the right window?

That is the key distinction between evaluating an approach and evaluating a
formula.

## Constructed Exact Best Case

Here is a small artificial case where the threshold formula is exactly optimal.
Use:

```text
N = 4
L = 2
K = 2
sigma = 1
P = 1
```

Let:

```text
V =
[
  [1.0, 0.0],
  [0.0, 1.0],
  [0.1, 0.0],
  [0.0, 0.1],
]
```

Row powers are:

```text
p = [1.0, 1.0, 0.01, 0.01]
```

The descending row-power order is:

```text
order = [0, 1, 2, 3]
```

For the single formula:

```text
T = round(0.15K) = round(0.3) = 0
```

For the 3-threshold portfolio:

```text
round(0.05N) = round(0.2) = 0
round(0.15K) = round(0.3) = 0
round(0.125NL/(L+2)) = round(0.25) = 0
```

So both the single formula and portfolio select:

```text
S_threshold = {0, 1}
```

Exhaustive exact-`K` search gives:

| subset | `U_G` |
|---|---:|
| `{0, 1}` | 4.0000 |
| `{1, 3}` | 2.0201 |
| `{0, 2}` | 2.0201 |
| `{1, 2}` | 2.0002 |
| `{0, 3}` | 2.0002 |
| `{2, 3}` | 1.0201 |

Therefore:

```text
S_threshold = S_opt
U_G(threshold) / U_G(optimum) = 1.0
```

This is the clean best case: the strongest two rows are orthogonal, have the
largest powers, and form the capacity-relevant balanced channel.

## Best-Case Data Regimes

The approach tends to be optimal or near-optimal when:

1. Row power ranking is informative.
   High row power generally means the row contributes useful channel energy.

2. The best subset is approximately a contiguous power window.
   The algorithm cannot choose arbitrary separated rows, so it works best when
   there is one good band in the sorted power order.

3. The top few rows are too costly or too correlated.
   Because `z^2 = P / max_row_power`, a single extreme row can reduce the scale
   applied to all selected rows. Skipping a small number of strongest rows can
   improve `U_G`.

4. Rows inside the selected window have good angular diversity.
   `U_G` is a determinant, so it benefits from balanced eigenvalues of
   `z^2 G_S G_S^H + sigma I`. A window with high energy spread across all `L`
   dimensions is much better than a high-power but nearly rank-one window.

5. The distribution resembles the non-Rician profiles from the scaling study.
   In the current experiments, `gaussian`, `rayleigh`, `nakagami`, `thin_tail`,
   `lognormal`, and `twdp` usually had best-tested threshold scale around:

```text
T/N ~= 0.048..0.059
```

The Rician profile shifted much farther right:

```text
T/N ~= 0.139
```

This is why the 3-threshold portfolio is safer than one fixed rule.

## Worst-Case Data Regimes

The threshold approach can be poor when row-power order is misleading.

## Constructed Exact Bad Case

Here is a deterministic artificial case where the 3-threshold portfolio is far
from the exact optimum. It is intentionally adversarial, not a realistic fading
sample.

Use:

```text
N = 20
L = 10
K = 10
sigma = 1
P = 1
```

Define standard basis vectors `e_1, ..., e_10`. Rows `0..9` are almost unit
power and all point in the same direction `e_1`:

```text
V[i] = sqrt(1.000 - 0.001 i) * e_1,  for i = 0..9
```

Rows `10..18` have slightly smaller row power but cover the missing orthogonal
directions:

```text
V[10+j] = sqrt(0.980 - 0.001 j) * e_{j+2},  for j = 0..8
```

The last row is a duplicate low-power `e_2` row:

```text
V[19] = sqrt(0.971) * e_2
```

The row powers are strictly descending:

```text
p =
[
  1.000, 0.999, 0.998, 0.997, 0.996,
  0.995, 0.994, 0.993, 0.992, 0.991,
  0.980, 0.979, 0.978, 0.977, 0.976,
  0.975, 0.974, 0.973, 0.972, 0.971
]
```

So the row-power order is exactly:

```text
order = [0, 1, 2, ..., 19]
```

For the recommended portfolio:

```text
round(0.05N) = round(1.0) = 1
round(0.15K) = round(1.5) = 2
round(0.125NL/(L+2)) = round(2.0833...) = 2
```

Thus it tests only:

```text
T in {1, 2}
```

The threshold windows are:

```text
T = 1: S = {1,2,3,4,5,6,7,8,9,10}
T = 2: S = {2,3,4,5,6,7,8,9,10,11}
```

Exhaustive exact-`K` search over all `C(20,10) = 184756` subsets gives:

```text
S_opt = {7,8,9,10,11,12,13,14,15,19}
U_G(S_opt) = 1391.167888237965
```

The best portfolio window is:

```text
T = 2
S_threshold = {2,3,4,5,6,7,8,9,10,11}
U_G(S_threshold) = 247.83352564028513
```

Therefore:

```text
U_G(S_threshold) / U_G(S_opt) = 0.17814781935068083
```

Why it fails:

- The sorted top rows are high-power but almost completely collinear.
- The threshold windows with `T=1` or `T=2` still contain too many collinear
  `e_1` rows.
- The optimum deliberately uses fewer top collinear rows and more slightly
  weaker orthogonal rows, producing a much more balanced determinant.

This family can be made worse as `K` and `L` grow: put the first `K` rows
almost collinear in one direction, and put the next `K` rows slightly lower in
power but spread across many orthogonal dimensions. A fixed small threshold
portfolio will keep too many collinear high-power rows, while the optimum can
trade some row power for many more active eigenmodes.

### Worst Case 1: Good Rows Are Not Contiguous

Suppose the best exact-`K` set alternates through the sorted row-power list:

```text
good rows:  order[0], order[2], order[4], ...
bad rows:   order[1], order[3], order[5], ...
```

No single contiguous window can select only the good rows. A 3-window portfolio
still cannot represent a highly interleaved optimum.

### Worst Case 2: Selected Window Is Nearly Rank-One

For `L > 1`, a selected window can have high row powers but poor spatial
diversity. If the selected rows are almost parallel:

```text
rank(G_S) ~= 1
```

then most eigenmodes are weak. The determinant objective becomes small because
it multiplies eigenmode terms. A lower-power set with more orthogonal rows can
produce much larger `U_G`.

### Worst Case 3: Extreme Power Outlier Inside the Window

If a very strong row lies inside the chosen window, it sets:

```text
z^2 = P / max_row_power
```

This can shrink the contribution of all selected rows. If that outlier is also
correlated with other rows, the window can have both bad scaling and poor
eigenvalue balance.

### Worst Case 4: Distribution Shift

A fixed formula can miss when the best threshold scale changes strongly with
the channel profile. In the scaling experiment, Rician shifted much farther
right than the other profiles. This does not make the threshold idea invalid,
but it means one formula is less robust than a small portfolio.

### Worst Case 5: Flat Row Powers and Unstable Ordering

If many rows have almost equal row power, small noise can permute their order.
The selected window may become unstable even though the threshold `T` is fixed.
In this case the row-power rank is a weak signal, and local refinement or a
diversity-aware candidate generator can help.

## What To Use Next

For the next experiment, use the 3-threshold portfolio:

```text
T candidates:
    round(0.05 N)
    round(0.15 K)
    round(0.125 N L / (L + 2))
```

Then score all candidates by exact `U_G` and return the best.

This keeps the algorithm simple and has a clear empirical gain over one
threshold:

```text
0.8927 -> 0.9366 mean fraction of best tested U_G
```

If we later want stronger quality, the natural next step is not a larger
threshold list; it is local refinement after the best threshold-window seed.

## Exact-Study Command

The small Gaussian exact study uses direct brute force, so it is limited to
small `N`. The accepted max point is `N=24`; calibration measured
`N=24, K=12` at about `87.041s` with the initial Python loop and about
`3.586s` after batched exact enumeration.

```bash
venv/bin/python -m experiments.algorithm_comparison \
  --threshold-exact-study \
  --N-values 8 12 16 20 24 \
  --L 2 \
  --K-pcts 25 50 \
  --samples 100 \
  --generator-seeds 42 \
  --data-profiles gaussian \
  --sigma 1 \
  --exact-time-limit 120 \
  --out-dir results/threshold_exact_gaussian_L2_N8_12_16_20_24_Kpct25_50_s100
```

Expected artifacts:

```text
threshold_runs.csv.gz
exact_runs.csv
exact_formula_runs.csv
exact_summary.csv
exact_formula_summary.csv
threshold_exact_report.md
```

## Checked Links

External references:

- Fading overview: https://en.wikipedia.org/wiki/Fading
- Rayleigh fading: https://en.wikipedia.org/wiki/Rayleigh_fading
- Rician fading: https://en.wikipedia.org/wiki/Rician_fading
- TWDP fading: https://en.wikipedia.org/wiki/Two-wave_with_diffuse_power_fading
- MIMO capacity/log-det context: https://en.wikipedia.org/wiki/MIMO
- Channel capacity background: https://en.wikipedia.org/wiki/Channel_capacity
- Sorting complexity background: https://en.wikipedia.org/wiki/Sorting_algorithm
- Comparison sort lower bound: https://en.wikipedia.org/wiki/Comparison_sort
- Determinant/LU complexity background: https://en.wikipedia.org/wiki/Leibniz_formula_for_determinants
- Space complexity definition: https://en.wikipedia.org/wiki/Space_complexity

Local code references:

- Objective equations: `algorithms/common.py`
- Threshold window and prefix sweep: `algorithms/h3_threshold_explore.py`
- Brute-force exact helper: `utils/brute_force.py`
- Scaling, exact-study orchestration, and reports: `experiments/algorithm_comparison.py`
- Scaling results: `results/threshold_scaling_prelim_N500_1000_2000_5000_L4_6_8_10_Kpct25_50_s100/`
