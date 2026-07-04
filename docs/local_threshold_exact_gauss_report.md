# Local Threshold Exact Gaussian Study

> Historical K semantics note: this report uses active-K semantics. Here `K` is the number of selected/kept antennas, not the number turned off. A `25% active` or `K=0.25N` case means `75% off`, not the real `25% off` task. For real off-percent experiments, `25% off => K_active=0.75N` and `50% off => K_active=0.50N`.

This experiment starts from a pure row-power threshold window and then applies greedy one-swap or two-swap local search by `U_G`.
It uses saved exact Gaussian cases and reconstructs matrices from `(profile, seed, sample)`; brute force is not rerun here.

## Setup

- Exact source: `results/threshold_exact_gaussian_L2_N8_12_16_20_24_Kpct25_to_50_s100`
- Profiles: gaussian
- N values: 8, 12, 16, 20
- L values: 2
- Requested active K percentages: 25, 30, 35, 40, 45, 50
- Seed rules: `best_tested_T`, `T=0`, `T=0.025N`, `T=0.05N`.
- Local search: all active rows are removable; add candidates are inactive rows near the window boundary.

## Direct Answer

- `best_tested_T + 0 swaps`: mean fraction exact `0.9904`, p05 `0.9354`, exact recovery `74.2%`, mean swaps applied `0.000`.
- `best_tested_T + 1 swaps`: mean fraction exact `0.9993`, p05 `1.0000`, exact recovery `97.5%`, mean swaps applied `0.248`.
- `best_tested_T + 2 swaps`: mean fraction exact `0.9996`, p05 `1.0000`, exact recovery `98.8%`, mean swaps applied `0.262`.
- Two swaps improve the headline threshold-window approach by `0.915` percentage points in mean exact fraction.
- `T=0.05N + 2 swaps` reaches mean fraction exact `0.9993` with exact recovery `98.3%`.
- After two swaps from `best_tested_T`, remaining misses have mean overlap `0.614` with exact and mean swap distance `2.310` rows.

## Global Rule Summary

| seed rule | swaps allowed | cases | mean fraction exact | p05 | p50 | p95 | exact rate | near-99 rate | mean runtime s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| T_0 | 0 | 2400 | 0.9593 | 0.7690 | 1.0000 | 1.0000 | 58.8% | 63.7% | 1.302e-04 |
| T_0 | 1 | 2400 | 0.9954 | 0.9771 | 1.0000 | 1.0000 | 92.1% | 93.8% | 0.0012 |
| T_0 | 2 | 2400 | 0.9975 | 1.0000 | 1.0000 | 1.0000 | 96.6% | 97.2% | 0.0016 |
| T_0p025N | 0 | 2400 | 0.9593 | 0.7690 | 1.0000 | 1.0000 | 58.8% | 63.7% | 1.308e-04 |
| T_0p025N | 1 | 2400 | 0.9954 | 0.9771 | 1.0000 | 1.0000 | 92.1% | 93.8% | 0.0012 |
| T_0p025N | 2 | 2400 | 0.9975 | 1.0000 | 1.0000 | 1.0000 | 96.6% | 97.2% | 0.0016 |
| T_0p05N | 0 | 2400 | 0.8686 | 0.6369 | 0.9016 | 1.0000 | 30.1% | 33.0% | 1.297e-04 |
| T_0p05N | 1 | 2400 | 0.9972 | 0.9839 | 1.0000 | 1.0000 | 91.9% | 93.8% | 0.0012 |
| T_0p05N | 2 | 2400 | 0.9993 | 1.0000 | 1.0000 | 1.0000 | 98.3% | 98.7% | 0.0021 |
| best_tested_T | 0 | 2400 | 0.9904 | 0.9354 | 1.0000 | 1.0000 | 74.2% | 80.3% | 1.515e-04 |
| best_tested_T | 1 | 2400 | 0.9993 | 1.0000 | 1.0000 | 1.0000 | 97.5% | 98.4% | 0.0012 |
| best_tested_T | 2 | 2400 | 0.9996 | 1.0000 | 1.0000 | 1.0000 | 98.8% | 99.1% | 0.0015 |

## Best Tested T Local Search By K%

| requested K% | swaps | mean fraction exact | p05 | exact rate | mean swaps applied |
|---:|---:|---:|---:|---:|---:|
| 25 | 0 | 0.9911 | 0.9155 | 79.0% | 0.000 |
| 25 | 1 | 0.9996 | 1.0000 | 98.5% | 0.200 |
| 25 | 2 | 0.9996 | 1.0000 | 98.8% | 0.203 |
| 30 | 0 | 0.9881 | 0.8997 | 74.8% | 0.000 |
| 30 | 1 | 0.9995 | 1.0000 | 98.5% | 0.247 |
| 30 | 2 | 0.9997 | 1.0000 | 99.5% | 0.258 |
| 35 | 0 | 0.9901 | 0.9121 | 74.0% | 0.000 |
| 35 | 1 | 0.9990 | 0.9997 | 97.8% | 0.247 |
| 35 | 2 | 0.9993 | 1.0000 | 98.5% | 0.258 |
| 40 | 0 | 0.9914 | 0.9289 | 74.5% | 0.000 |
| 40 | 1 | 0.9990 | 0.9934 | 96.8% | 0.240 |
| 40 | 2 | 0.9994 | 1.0000 | 98.2% | 0.258 |
| 45 | 0 | 0.9913 | 0.9126 | 71.8% | 0.000 |
| 45 | 1 | 0.9994 | 0.9936 | 96.8% | 0.273 |
| 45 | 2 | 0.9997 | 1.0000 | 98.8% | 0.295 |
| 50 | 0 | 0.9906 | 0.9231 | 71.2% | 0.000 |
| 50 | 1 | 0.9994 | 0.9973 | 96.8% | 0.280 |
| 50 | 2 | 0.9997 | 1.0000 | 99.0% | 0.302 |

## Worst Remaining Best-T + 2-Swap Cases

| N | K | requested K% | sample | T | fraction exact | overlap | swap distance | exact subset | local subset |
|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| 20 | 7 | 35 | 61 | 4 | 0.9115 | 0.429 | 4 | 0 3 5 8 9 12 17 | 3 6 9 13 15 17 19 |
| 16 | 6 | 40 | 37 | 3 | 0.9189 | 0.833 | 1 | 0 2 3 4 10 12 | 0 3 4 10 11 12 |
| 16 | 6 | 35 | 37 | 3 | 0.9189 | 0.833 | 1 | 0 2 3 4 10 12 | 0 3 4 10 11 12 |
| 8 | 2 | 25 | 35 | 1 | 0.9214 | 0.000 | 2 | 0 6 | 2 7 |
| 8 | 2 | 30 | 35 | 1 | 0.9214 | 0.000 | 2 | 0 6 | 2 7 |
| 20 | 8 | 40 | 61 | 4 | 0.9380 | 0.500 | 4 | 0 3 5 6 8 9 12 17 | 1 3 6 9 13 15 17 19 |
| 20 | 10 | 50 | 18 | 3 | 0.9405 | 0.700 | 3 | 1 3 4 5 6 7 8 14 15 19 | 1 3 4 6 9 13 14 15 16 19 |
| 20 | 8 | 40 | 41 | 3 | 0.9566 | 0.750 | 2 | 1 2 7 8 9 12 18 19 | 2 6 7 8 9 12 16 18 |
| 20 | 7 | 35 | 71 | 2 | 0.9609 | 0.714 | 2 | 2 3 7 12 14 16 19 | 0 2 7 8 14 16 19 |
| 20 | 10 | 50 | 61 | 4 | 0.9632 | 0.600 | 4 | 0 3 5 6 8 9 12 15 17 19 | 1 3 4 6 9 13 15 17 18 19 |

## Interpretation

- If `best_tested_T + 1/2 swaps` closes most remaining exact gap, the threshold-window approach is structurally strong: sorting by row power finds almost the right neighborhood, and local swaps repair non-contiguous exact subsets.
- Remaining misses mean the exact subset is outside the fixed boundary neighborhood, or needs a coordinated multi-row replacement that greedy one-at-a-time swaps cannot see.
- `T=0.05N + local search` measures whether the practical formula can benefit from the same repair step, while `best_tested_T + local search` measures the approach itself independent of formula choice.

## Plots

![Raw U_G CDF by K percentage](local_threshold_raw_u_g_cdf_by_k_pct.png)

![Fraction exact CDF](local_threshold_fraction_exact_cdf.png)

![Mean fraction by active percentage](local_threshold_mean_fraction_by_active_pct.png)

![Exact recovery rate](local_threshold_exact_rate_by_active_pct.png)

![Seed threshold dependence](local_threshold_seed_dependence.png)

![Failure diagnostics](local_threshold_failure_diagnostics.png)

![Runtime](local_threshold_runtime_by_N_K.png)

## Artifacts

- `local_threshold_runs.csv`
- `local_threshold_summary.csv`
- `local_threshold_failure_cases.csv`
- `local_threshold_diagnostics.csv`
