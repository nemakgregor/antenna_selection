# Real-Off Cyclic Threshold Local Search Study

This report uses the real task condition: turn off a requested percentage of antennas while keeping the rest active.

## K Semantics

- `K_off = round(N * off_pct / 100)` is the number of disabled antennas.
- `K_active = N - K_off` is the number of selected/kept antennas.
- The solver variable `K` equals `K_active`; it does not mean the number of antennas turned off.
- Therefore `25% off` means `K_active = 0.75N`, and `50% off` means `K_active = 0.50N`.

## Setup

- Exact source: `results/threshold_exact_gaussian_L2_N8_12_16_20_24_Kpct25_to_50_s100`
- Profiles: gaussian
- N values: 8, 12
- L values: 2
- Off percentages: 25, 50
- Seeds compared: best tested cyclic window, `T=round(0.05N)`, and strong/weak.
- Local search: greedy remove-one/add-one refinement by raw `U_G`, with 0, 1, or 2 swaps.

## Direct Answer

- `cyclic best`: mean fraction exact `0.9953`, p05 `0.9755`, exact recovery `87.5%`, mean swaps applied `0.000`.
- `cyclic best + 1 swap`: mean fraction exact `1.0000`, p05 `1.0000`, exact recovery `100.0%`, mean swaps applied `0.125`.
- `cyclic best + 2 swaps`: mean fraction exact `1.0000`, p05 `1.0000`, exact recovery `100.0%`, mean swaps applied `0.125`.
- `T=0.05N`: mean fraction exact `0.8965`, p05 `0.6320`, exact recovery `50.0%`, mean swaps applied `0.000`.
- `T=0.05N + 1 swap`: mean fraction exact `1.0000`, p05 `1.0000`, exact recovery `100.0%`, mean swaps applied `0.500`.
- `T=0.05N + 2 swaps`: mean fraction exact `1.0000`, p05 `1.0000`, exact recovery `100.0%`, mean swaps applied `0.500`.
- `strong/weak`: mean fraction exact `0.6052`, p05 `0.3491`, exact recovery `0.0%`, mean swaps applied `0.000`.
- `strong/weak + 1 swap`: mean fraction exact `0.8242`, p05 `0.5330`, exact recovery `25.0%`, mean swaps applied `1.000`.
- `strong/weak + 2 swaps`: mean fraction exact `0.9558`, p05 `0.7703`, exact recovery `87.5%`, mean swaps applied `1.750`.
- Two swaps improve cyclic-threshold mean exact fraction by `0.471` percentage points.
- After two swaps from cyclic best, every analyzed case recovered exact `U_G`.

## Global Summary

| seed | swaps | cases | mean fraction exact | p05 | p50 | p95 | exact rate | near-99 rate | runtime s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| T=0.05N | 0 | 8 | 0.8965 | 0.6320 | 0.9812 | 1.0000 | 50.0% | 50.0% | 1.289e-04 |
| T=0.05N + 1 swap | 1 | 8 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 100.0% | 100.0% | 7.800e-04 |
| T=0.05N + 2 swaps | 2 | 8 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 100.0% | 100.0% | 0.0011 |
| cyclic best | 0 | 8 | 0.9953 | 0.9755 | 1.0000 | 1.0000 | 87.5% | 87.5% | 1.248e-04 |
| cyclic best + 1 swap | 1 | 8 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 100.0% | 100.0% | 7.610e-04 |
| cyclic best + 2 swaps | 2 | 8 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 100.0% | 100.0% | 7.980e-04 |
| strong/weak | 0 | 8 | 0.6052 | 0.3491 | 0.6355 | 0.8039 | 0.0% | 0.0% | 1.312e-04 |
| strong/weak + 1 swap | 1 | 8 | 0.8242 | 0.5330 | 0.8885 | 1.0000 | 25.0% | 25.0% | 7.562e-04 |
| strong/weak + 2 swaps | 2 | 8 | 0.9558 | 0.7703 | 1.0000 | 1.0000 | 87.5% | 87.5% | 0.0012 |

## By Real Off Percentage

| off % | K semantics | seed | swaps | mean fraction exact | p05 | exact rate |
|---:|---|---|---:|---:|---:|---:|
| 25 | K_active mean 7.5, K_off mean 2.5 | T=0.05N | 0 | 0.8907 | 0.6398 | 50.0% |
| 25 | K_active mean 7.5, K_off mean 2.5 | T=0.05N + 1 swap | 1 | 1.0000 | 1.0000 | 100.0% |
| 25 | K_active mean 7.5, K_off mean 2.5 | T=0.05N + 2 swaps | 2 | 1.0000 | 1.0000 | 100.0% |
| 25 | K_active mean 7.5, K_off mean 2.5 | cyclic best | 0 | 1.0000 | 1.0000 | 100.0% |
| 25 | K_active mean 7.5, K_off mean 2.5 | cyclic best + 1 swap | 1 | 1.0000 | 1.0000 | 100.0% |
| 25 | K_active mean 7.5, K_off mean 2.5 | cyclic best + 2 swaps | 2 | 1.0000 | 1.0000 | 100.0% |
| 25 | K_active mean 7.5, K_off mean 2.5 | strong/weak | 0 | 0.6743 | 0.4145 | 0.0% |
| 25 | K_active mean 7.5, K_off mean 2.5 | strong/weak + 1 swap | 1 | 0.8907 | 0.6398 | 50.0% |
| 25 | K_active mean 7.5, K_off mean 2.5 | strong/weak + 2 swaps | 2 | 1.0000 | 1.0000 | 100.0% |
| 50 | K_active mean 5.0, K_off mean 5.0 | T=0.05N | 0 | 0.9022 | 0.6642 | 50.0% |
| 50 | K_active mean 5.0, K_off mean 5.0 | T=0.05N + 1 swap | 1 | 1.0000 | 1.0000 | 100.0% |
| 50 | K_active mean 5.0, K_off mean 5.0 | T=0.05N + 2 swaps | 2 | 1.0000 | 1.0000 | 100.0% |
| 50 | K_active mean 5.0, K_off mean 5.0 | cyclic best | 0 | 0.9906 | 0.9642 | 75.0% |
| 50 | K_active mean 5.0, K_off mean 5.0 | cyclic best + 1 swap | 1 | 1.0000 | 1.0000 | 100.0% |
| 50 | K_active mean 5.0, K_off mean 5.0 | cyclic best + 2 swaps | 2 | 1.0000 | 1.0000 | 100.0% |
| 50 | K_active mean 5.0, K_off mean 5.0 | strong/weak | 0 | 0.5362 | 0.3358 | 0.0% |
| 50 | K_active mean 5.0, K_off mean 5.0 | strong/weak + 1 swap | 1 | 0.7577 | 0.4977 | 0.0% |
| 50 | K_active mean 5.0, K_off mean 5.0 | strong/weak + 2 swaps | 2 | 0.9116 | 0.6642 | 75.0% |

## Notes

- Exact enumeration completed for `100.0%` of cases.
- Exact rows loaded from previous artifacts for `50.0%` of cases; missing real-off cases were recomputed.
- Strong/weak under real-off semantics disables `K_off` antennas from the weakest and strongest row-power tails, then keeps the middle `K_active` antennas.
- Historical reports that say `25% active` are not the real `25% off` task. `25% active` means `75% off`.

## Plots

![Raw U_G CDF](real_off_raw_u_g_cdf_by_off_pct.png)

![Fraction exact CDF](real_off_fraction_exact_cdf.png)

![Mean fraction by off percentage](real_off_mean_fraction_by_off_pct.png)

![Exact recovery by off percentage](real_off_exact_recovery_by_off_pct.png)

![Best cyclic start](real_off_best_cyclic_start_hist.png)

![Failure diagnostics](real_off_failure_diagnostics.png)

![Runtime](real_off_runtime_by_method.png)

## Artifacts

- Detailed CSVs are packed in `csv_data.tar.gz` after the run.
- Main report: `local_threshold_real_off_report.md`.
