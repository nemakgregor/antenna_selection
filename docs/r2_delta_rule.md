# R2+Delta Analytical Rule

`R2Delta-Gen` is implemented in `algorithms/r2_delta.py` as
`solve_r2_delta_gen(V, K, sigma=1.0, P=1.0)`. It maximizes the general
objective

```text
U_G(S) = det(sigma I + z^2 G_S^2)
G_S = sum_{n in S} v_n^H v_n
z^2 = P / max_{n in S} ||v_n||^2
```

The method is deterministic and returns exactly `K` active antennas.

## Step 1: R2 Power Window

Sort antennas by row power `p_n = ||v_n||^2` in descending order. For every
contiguous window of length `K`, compute

```text
T = sum p_n
Q = ||G_window||_F^2
cap = max p_n in the window
c = P / cap
```

R2 approximates the spectrum of `G_window` by one spike `a` and a flat
background `b`:

```text
a + (L - 1)b = T
a^2 + (L - 1)b^2 = Q
a = (T + sqrt((L - 1)(LQ - T^2))) / L
b = (T - a) / (L - 1)
```

For `L = 1`, the score is exact:

```text
score = log(sigma + c T^2)
```

For `L > 1`, the R2 window score is:

```text
score = log(sigma + c a^2) + (L - 1) log(sigma + c b^2)
```

The best-scoring window is the initial active set.

## Step 2: Delta Pair Swaps

For the current active set, compute the Gram matrix `G`, the cap, and
`c = P / cap`. Eigendecompose `G`:

```text
G u_i = lambda_i u_i
h_i = lambda_i / (sigma + c lambda_i^2)
```

Each antenna receives a derivative usefulness score:

```text
s(n) = sum_i h_i |v_n u_i|^2
```

At each round:

- remove candidates are the active antennas with the smallest `s(n)`;
- add candidates are inactive antennas with the largest `s(n)` and
  `p_n <= current cap`;
- every shortlisted remove/add pair is evaluated with the exact `U_G`;
- only the best strictly improving exact swap is accepted.

This keeps the local search monotone in the true objective. The default
parameters are `shortlist_size=8` and `max_rounds=250`.

## Validation Status

The historical report supplied for this project described strong synthetic
benchmark results for R2+Delta. The current repository keeps the algorithm and
its formulas, but large historical benchmark numbers should be treated as
historical until rerun from this checkout. The maintained regression suite
checks the scalar analytical example, exact cardinality, registry integration,
and that Delta swaps never reduce the exact objective relative to the R2
window on the small test instance.
