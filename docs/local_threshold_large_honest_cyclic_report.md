# Honest All-Inactive One-Swap Large-N Study

This report studies the threshold-window approach without brute-force exact enumeration.

## K Semantics

- `K` means active/kept antennas selected by the solver.
- `K_off = N - K` is the number of disabled antennas.
- For `N=1000`, `K=750` means `25% off`; `K=500` means `50% off`.

## Setup

- Profiles: gaussian
- N values: 1000
- L values: 2
- K_active values: 500, 750
- Seeds compared: best tested cyclic window, `T=round(0.05N)`, and strong/weak.
- Local search: greedy remove-one/add-one refinement by raw `U_G`, with swaps 0, 1.

## Local Swap Scheme And Cost

- A `1-swap` removes one currently active antenna and adds one currently inactive antenna, preserving exact `K`.
- One local-search pass evaluates every `(remove, add)` pair from the current active set and the configured add-candidate pool, applies only the single pair with the largest positive `U_G` improvement, and stops if no pair improves `U_G`.
- Honest mode uses all inactive antennas as add candidates, so `A = N - K`.
- One honest pass costs `O(K * (N-K) * L^3)` time after row Gram matrices are built.
- Best cyclic seed construction scans `N` cyclic windows before local refinement; `T=0.05N` and strong/weak construct one sorted-window seed.
- Extra working space is `O(N * L^2 + K + (N-K) + L^2)` for row Gram matrices, active/add sets, and the current Gram matrix.

## Direct Answer

- Best observed tested method by mean `U_G`: `cyclic best + 1 swap` with mean fraction of best observed `1.0000`.
- Cyclic best + 1 swap changes mean `U_G` by `0.013`% of the cyclic seed baseline.
- `T=0.05N + 1 swap` reaches mean fraction of cyclic seed `0.9675`.

## Honest Versus Previous Bounded Local Search

- `cyclic best + 1 honest swap` changes mean fraction versus bounded by `0.000000` on average across K values.
- `T=0.05N + 1 honest swap` changes mean fraction versus bounded by `0.000000` on average across K values.
- `strong/weak + 1 honest swap` changes mean fraction versus bounded by `0.000000` on average across K values.

| K_active | seed | bounded frac cyclic seed | honest frac cyclic seed | delta | bounded runtime s | honest runtime s | runtime ratio | bounded pairs | honest pairs |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 500 | T=0.05N + 1 swap | 0.9669 | 0.9669 | 0.000000 | 0.6392 | 6.2274 | 9.7422 | 25000.0000 | 250000.0000 |
| 500 | cyclic best + 1 swap | 1.0002 | 1.0002 | 0.000000 | 0.8269 | 6.4401 | 7.7882 | 25000.0000 | 250000.0000 |
| 500 | strong/weak + 1 swap | 0.5648 | 0.5648 | 0.000000 | 0.6412 | 6.2515 | 9.7501 | 25000.0000 | 250000.0000 |
| 750 | T=0.05N + 1 swap | 0.9681 | 0.9681 | 0.000000 | 1.4616 | 4.7056 | 3.2196 | 57000.0000 | 187500.0000 |
| 750 | cyclic best + 1 swap | 1.0001 | 1.0001 | 0.000000 | 1.7117 | 4.9393 | 2.8857 | 57000.0000 | 187500.0000 |
| 750 | strong/weak + 1 swap | 0.8391 | 0.8391 | 0.000000 | 1.4532 | 4.7030 | 3.2363 | 57000.0000 | 187500.0000 |

## Best Cyclic T

`T` is the cyclic start position in descending row-power order. `T=0` starts at the strongest row; larger `T` shifts the selected cyclic window.

| K_active | off % | cases | T mean | T std | T p05 | T median | T p95 | T/N mean | T/N median |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 500 | 50 | 100 | 50.9800 | 16.6963 | 28.8000 | 48.0000 | 80.1000 | 0.0510 | 0.0480 |
| 750 | 25 | 100 | 50.3700 | 15.9054 | 29.0000 | 48.0000 | 79.0500 | 0.0504 | 0.0480 |

## Method Summary

| K_active | seed | swaps | cases | mean U_G | p05 U_G | mean frac cyclic seed | p05 frac cyclic seed | mean frac best observed | mean total runtime s | mean swaps applied |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 500 | T=0.05N | 0 | 100 | 8.979e+07 | 5.712e+07 | 0.9611 | 0.9168 | 0.9609 | 7.095e-04 | 0.000 |
| 500 | T=0.05N + 1 swap | 1 | 100 | 9.033e+07 | 5.714e+07 | 0.9669 | 0.9242 | 0.9667 | 6.2274 | 0.990 |
| 500 | cyclic best | 0 | 100 | 9.353e+07 | 5.945e+07 | 1.0000 | 1.0000 | 0.9998 | 0.1812 | 0.000 |
| 500 | cyclic best + 1 swap | 1 | 100 | 9.355e+07 | 5.946e+07 | 1.0002 | 1.0000 | 1.0000 | 6.4401 | 0.990 |
| 500 | strong/weak | 0 | 100 | 5.242e+07 | 3.186e+07 | 0.5620 | 0.5069 | 0.5619 | 6.657e-04 | 0.000 |
| 500 | strong/weak + 1 swap | 1 | 100 | 5.268e+07 | 3.207e+07 | 0.5648 | 0.5093 | 0.5647 | 6.2515 | 1.000 |
| 750 | T=0.05N | 0 | 100 | 2.012e+08 | 1.249e+08 | 0.9624 | 0.9208 | 0.9623 | 7.786e-04 | 0.000 |
| 750 | T=0.05N + 1 swap | 1 | 100 | 2.024e+08 | 1.249e+08 | 0.9681 | 0.9250 | 0.9680 | 4.7056 | 0.970 |
| 750 | cyclic best | 0 | 100 | 2.093e+08 | 1.300e+08 | 1.0000 | 1.0000 | 0.9999 | 0.2513 | 0.000 |
| 750 | cyclic best + 1 swap | 1 | 100 | 2.093e+08 | 1.300e+08 | 1.0001 | 1.0000 | 1.0000 | 4.9393 | 0.970 |
| 750 | strong/weak | 0 | 100 | 1.745e+08 | 1.075e+08 | 0.8351 | 0.7725 | 0.8350 | 7.234e-04 | 0.000 |
| 750 | strong/weak + 1 swap | 1 | 100 | 1.754e+08 | 1.079e+08 | 0.8391 | 0.7759 | 0.8390 | 4.7030 | 0.990 |

## Plots

![Raw U_G CDF](large_raw_u_g_cdf_by_K.png)

![Fraction cyclic seed CDF](large_fraction_cyclic_seed_cdf_by_K.png)

![Fraction best observed CDF](large_fraction_best_observed_cdf_by_K.png)

![Mean fraction cyclic seed](large_mean_fraction_cyclic_seed_by_K.png)

![Best cyclic T boxplot](large_best_cyclic_T_boxplot.png)

![Best cyclic T over N](large_best_cyclic_T_over_N.png)

![Runtime](large_runtime_by_method.png)

## Artifacts

- Detailed CSVs are packed in `csv_data.tar.gz` after the run.
- Main report: `local_threshold_large_honest_report.md`.
