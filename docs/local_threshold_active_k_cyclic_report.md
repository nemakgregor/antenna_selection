# Active-K Cyclic Threshold Local Search Study

This exploratory report sweeps requested active antenna fractions directly.

## K Semantics

- `K_active = round(N * requested_active_pct / 100)` is the number of selected/kept antennas.
- `K_off = N - K_active` is the number of disabled antennas.
- The solver variable `K` equals `K_active`; it does not mean the number of antennas turned off.

## Setup

- Exact source: `results/threshold_exact_gaussian_L2_N8_12_16_20_24_Kpct25_to_50_s100`
- Profiles: gaussian
- N values: 8, 12, 16, 20, 24
- L values: 2
- Requested active K percentages: 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75
- Seeds compared: best tested cyclic window, `T=round(0.05N)`, and strong/weak.
- Local search: greedy remove-one/add-one refinement by raw `U_G`, with 0, 1, or 2 swaps.

## Direct Answer

- `cyclic best`: mean fraction exact `0.9915`, p05 `0.9487`, exact recovery `71.7%`, mean swaps applied `0.000`.
- `cyclic best + 1 swap`: mean fraction exact `0.9994`, p05 `1.0000`, exact recovery `97.5%`, mean swaps applied `0.273`.
- `cyclic best + 2 swaps`: mean fraction exact `0.9996`, p05 `1.0000`, exact recovery `98.9%`, mean swaps applied `0.287`.
- `T=0.05N`: mean fraction exact `0.8652`, p05 `0.6273`, exact recovery `29.2%`, mean swaps applied `0.000`.
- `T=0.05N + 1 swap`: mean fraction exact `0.9960`, p05 `0.9810`, exact recovery `90.7%`, mean swaps applied `0.704`.
- `T=0.05N + 2 swaps`: mean fraction exact `0.9985`, p05 `1.0000`, exact recovery `97.7%`, mean swaps applied `0.783`.
- `strong/weak`: mean fraction exact `0.5344`, p05 `0.2760`, exact recovery `0.6%`, mean swaps applied `0.000`.
- `strong/weak + 1 swap`: mean fraction exact `0.6616`, p05 `0.3412`, exact recovery `7.8%`, mean swaps applied `0.933`.
- `strong/weak + 2 swaps`: mean fraction exact `0.7746`, p05 `0.4029`, exact recovery `28.9%`, mean swaps applied `1.696`.

## Best Cyclic T Usually

`T` is the cyclic start position in descending row-power order. `T=0` starts at the strongest row; larger `T` shifts the active window toward weaker rows.

| requested active % | actual active % mean | T p05 | T median | T p95 | T/N median | T/K median |
|---:|---:|---:|---:|---:|---:|---:|
| 25 | 25.00 | 0.0 | 0.0 | 2.0 | 0.000 | 0.000 |
| 30 | 29.75 | 0.0 | 0.0 | 2.0 | 0.000 | 0.000 |
| 35 | 35.33 | 0.0 | 0.0 | 2.0 | 0.000 | 0.000 |
| 40 | 39.67 | 0.0 | 0.0 | 2.0 | 0.000 | 0.000 |
| 45 | 45.25 | 0.0 | 0.0 | 2.0 | 0.000 | 0.000 |
| 50 | 50.00 | 0.0 | 0.0 | 2.0 | 0.000 | 0.000 |
| 55 | 54.75 | 0.0 | 0.0 | 2.0 | 0.000 | 0.000 |
| 60 | 60.33 | 0.0 | 0.0 | 2.0 | 0.000 | 0.000 |
| 65 | 64.67 | 0.0 | 0.0 | 2.0 | 0.000 | 0.000 |
| 70 | 70.25 | 0.0 | 0.0 | 2.0 | 0.000 | 0.000 |
| 75 | 75.00 | 0.0 | 0.0 | 2.0 | 0.000 | 0.000 |

## By Requested Active K Percentage

| requested active % | actual active % mean | K_active mean | K_off mean | seed | swaps | mean fraction exact | p05 | exact rate |
|---:|---:|---:|---:|---|---:|---:|---:|---:|
| 25 | 25.00 | 4.0 | 12.0 | T=0.05N | 0 | 0.8746 | 0.6535 | 22.6% |
| 25 | 25.00 | 4.0 | 12.0 | T=0.05N + 1 swap | 1 | 0.9975 | 0.9435 | 92.0% |
| 25 | 25.00 | 4.0 | 12.0 | T=0.05N + 2 swaps | 2 | 0.9998 | 1.0000 | 98.8% |
| 25 | 25.00 | 4.0 | 12.0 | cyclic best | 0 | 0.9888 | 0.9121 | 74.4% |
| 25 | 25.00 | 4.0 | 12.0 | cyclic best + 1 swap | 1 | 0.9992 | 0.9892 | 97.4% |
| 25 | 25.00 | 4.0 | 12.0 | cyclic best + 2 swaps | 2 | 0.9993 | 1.0000 | 98.4% |
| 25 | 25.00 | 4.0 | 12.0 | strong/weak | 0 | 0.4793 | 0.2013 | 0.0% |
| 25 | 25.00 | 4.0 | 12.0 | strong/weak + 1 swap | 1 | 0.5768 | 0.2401 | 1.6% |
| 25 | 25.00 | 4.0 | 12.0 | strong/weak + 2 swaps | 2 | 0.6692 | 0.2777 | 15.6% |
| 30 | 29.75 | 4.8 | 11.2 | T=0.05N | 0 | 0.8658 | 0.6269 | 25.2% |
| 30 | 29.75 | 4.8 | 11.2 | T=0.05N + 1 swap | 1 | 0.9971 | 0.9564 | 91.8% |
| 30 | 29.75 | 4.8 | 11.2 | T=0.05N + 2 swaps | 2 | 0.9995 | 0.9998 | 97.6% |
| 30 | 29.75 | 4.8 | 11.2 | cyclic best | 0 | 0.9862 | 0.8978 | 70.0% |
| 30 | 29.75 | 4.8 | 11.2 | cyclic best + 1 swap | 1 | 0.9992 | 0.9996 | 97.6% |
| 30 | 29.75 | 4.8 | 11.2 | cyclic best + 2 swaps | 2 | 0.9994 | 1.0000 | 99.0% |
| 30 | 29.75 | 4.8 | 11.2 | strong/weak | 0 | 0.4766 | 0.1923 | 0.0% |
| 30 | 29.75 | 4.8 | 11.2 | strong/weak + 1 swap | 1 | 0.5768 | 0.2340 | 1.8% |
| 30 | 29.75 | 4.8 | 11.2 | strong/weak + 2 swaps | 2 | 0.6635 | 0.2696 | 16.2% |
| 35 | 35.33 | 5.6 | 10.4 | T=0.05N | 0 | 0.8671 | 0.6266 | 28.4% |
| 35 | 35.33 | 5.6 | 10.4 | T=0.05N + 1 swap | 1 | 0.9953 | 0.9325 | 89.2% |
| 35 | 35.33 | 5.6 | 10.4 | T=0.05N + 2 swaps | 2 | 0.9992 | 1.0000 | 98.2% |
| 35 | 35.33 | 5.6 | 10.4 | cyclic best | 0 | 0.9878 | 0.9121 | 68.8% |
| 35 | 35.33 | 5.6 | 10.4 | cyclic best + 1 swap | 1 | 0.9989 | 0.9997 | 97.4% |
| 35 | 35.33 | 5.6 | 10.4 | cyclic best + 2 swaps | 2 | 0.9992 | 1.0000 | 98.2% |
| 35 | 35.33 | 5.6 | 10.4 | strong/weak | 0 | 0.4697 | 0.1959 | 0.0% |
| 35 | 35.33 | 5.6 | 10.4 | strong/weak + 1 swap | 1 | 0.5673 | 0.2398 | 0.4% |
| 35 | 35.33 | 5.6 | 10.4 | strong/weak + 2 swaps | 2 | 0.6681 | 0.2938 | 2.8% |
| 40 | 39.67 | 6.4 | 9.6 | T=0.05N | 0 | 0.8642 | 0.6002 | 28.6% |
| 40 | 39.67 | 6.4 | 9.6 | T=0.05N + 1 swap | 1 | 0.9948 | 0.9433 | 88.0% |
| 40 | 39.67 | 6.4 | 9.6 | T=0.05N + 2 swaps | 2 | 0.9982 | 0.9844 | 96.8% |
| 40 | 39.67 | 6.4 | 9.6 | cyclic best | 0 | 0.9906 | 0.9234 | 71.2% |
| 40 | 39.67 | 6.4 | 9.6 | cyclic best + 1 swap | 1 | 0.9987 | 0.9847 | 95.4% |
| 40 | 39.67 | 6.4 | 9.6 | cyclic best + 2 swaps | 2 | 0.9991 | 0.9938 | 97.4% |
| 40 | 39.67 | 6.4 | 9.6 | strong/weak | 0 | 0.4791 | 0.2110 | 0.2% |
| 40 | 39.67 | 6.4 | 9.6 | strong/weak + 1 swap | 1 | 0.5811 | 0.2808 | 0.8% |
| 40 | 39.67 | 6.4 | 9.6 | strong/weak + 2 swaps | 2 | 0.6873 | 0.3268 | 3.6% |
| 45 | 45.25 | 7.2 | 8.8 | T=0.05N | 0 | 0.8658 | 0.6002 | 31.2% |
| 45 | 45.25 | 7.2 | 8.8 | T=0.05N + 1 swap | 1 | 0.9951 | 0.9350 | 88.0% |
| 45 | 45.25 | 7.2 | 8.8 | T=0.05N + 2 swaps | 2 | 0.9988 | 0.9796 | 96.8% |
| 45 | 45.25 | 7.2 | 8.8 | cyclic best | 0 | 0.9905 | 0.9126 | 68.6% |
| 45 | 45.25 | 7.2 | 8.8 | cyclic best + 1 swap | 1 | 0.9992 | 0.9936 | 96.4% |
| 45 | 45.25 | 7.2 | 8.8 | cyclic best + 2 swaps | 2 | 0.9996 | 1.0000 | 98.6% |
| 45 | 45.25 | 7.2 | 8.8 | strong/weak | 0 | 0.4818 | 0.1943 | 0.2% |
| 45 | 45.25 | 7.2 | 8.8 | strong/weak + 1 swap | 1 | 0.5878 | 0.2638 | 0.8% |
| 45 | 45.25 | 7.2 | 8.8 | strong/weak + 2 swaps | 2 | 0.7029 | 0.3268 | 20.2% |
| 50 | 50.00 | 8.0 | 8.0 | T=0.05N | 0 | 0.8646 | 0.5915 | 30.6% |
| 50 | 50.00 | 8.0 | 8.0 | T=0.05N + 1 swap | 1 | 0.9951 | 0.8969 | 88.4% |
| 50 | 50.00 | 8.0 | 8.0 | T=0.05N + 2 swaps | 2 | 0.9986 | 0.9990 | 97.6% |
| 50 | 50.00 | 8.0 | 8.0 | cyclic best | 0 | 0.9900 | 0.9231 | 68.0% |
| 50 | 50.00 | 8.0 | 8.0 | cyclic best + 1 swap | 1 | 0.9995 | 0.9973 | 96.6% |
| 50 | 50.00 | 8.0 | 8.0 | cyclic best + 2 swaps | 2 | 0.9997 | 1.0000 | 98.8% |
| 50 | 50.00 | 8.0 | 8.0 | strong/weak | 0 | 0.5291 | 0.2563 | 0.2% |
| 50 | 50.00 | 8.0 | 8.0 | strong/weak + 1 swap | 1 | 0.6543 | 0.3321 | 2.0% |
| 50 | 50.00 | 8.0 | 8.0 | strong/weak + 2 swaps | 2 | 0.7788 | 0.3995 | 26.4% |
| 55 | 54.75 | 8.8 | 7.2 | T=0.05N | 0 | 0.8621 | 0.5635 | 29.0% |
| 55 | 54.75 | 8.8 | 7.2 | T=0.05N + 1 swap | 1 | 0.9956 | 0.9273 | 89.4% |
| 55 | 54.75 | 8.8 | 7.2 | T=0.05N + 2 swaps | 2 | 0.9981 | 0.9967 | 97.6% |
| 55 | 54.75 | 8.8 | 7.2 | cyclic best | 0 | 0.9913 | 0.9448 | 67.4% |
| 55 | 54.75 | 8.8 | 7.2 | cyclic best + 1 swap | 1 | 0.9997 | 0.9998 | 98.0% |
| 55 | 54.75 | 8.8 | 7.2 | cyclic best + 2 swaps | 2 | 0.9998 | 1.0000 | 99.4% |
| 55 | 54.75 | 8.8 | 7.2 | strong/weak | 0 | 0.5142 | 0.2538 | 0.2% |
| 55 | 54.75 | 8.8 | 7.2 | strong/weak + 1 swap | 1 | 0.6432 | 0.3316 | 1.8% |
| 55 | 54.75 | 8.8 | 7.2 | strong/weak + 2 swaps | 2 | 0.7737 | 0.4017 | 25.4% |
| 60 | 60.33 | 9.6 | 6.4 | T=0.05N | 0 | 0.8636 | 0.5635 | 31.6% |
| 60 | 60.33 | 9.6 | 6.4 | T=0.05N + 1 swap | 1 | 0.9962 | 0.9464 | 91.0% |
| 60 | 60.33 | 9.6 | 6.4 | T=0.05N + 2 swaps | 2 | 0.9979 | 0.9938 | 97.6% |
| 60 | 60.33 | 9.6 | 6.4 | cyclic best | 0 | 0.9936 | 0.9448 | 74.6% |
| 60 | 60.33 | 9.6 | 6.4 | cyclic best + 1 swap | 1 | 0.9996 | 1.0000 | 97.8% |
| 60 | 60.33 | 9.6 | 6.4 | cyclic best + 2 swaps | 2 | 0.9997 | 1.0000 | 99.2% |
| 60 | 60.33 | 9.6 | 6.4 | strong/weak | 0 | 0.5547 | 0.2729 | 0.2% |
| 60 | 60.33 | 9.6 | 6.4 | strong/weak + 1 swap | 1 | 0.6894 | 0.4110 | 3.8% |
| 60 | 60.33 | 9.6 | 6.4 | strong/weak + 2 swaps | 2 | 0.8330 | 0.4966 | 33.0% |
| 65 | 64.67 | 10.4 | 5.6 | T=0.05N | 0 | 0.8640 | 0.5540 | 33.0% |
| 65 | 64.67 | 10.4 | 5.6 | T=0.05N + 1 swap | 1 | 0.9965 | 0.9440 | 93.6% |
| 65 | 64.67 | 10.4 | 5.6 | T=0.05N + 2 swaps | 2 | 0.9978 | 0.9990 | 97.8% |
| 65 | 64.67 | 10.4 | 5.6 | cyclic best | 0 | 0.9952 | 0.9612 | 74.8% |
| 65 | 64.67 | 10.4 | 5.6 | cyclic best + 1 swap | 1 | 0.9997 | 1.0000 | 98.4% |
| 65 | 64.67 | 10.4 | 5.6 | cyclic best + 2 swaps | 2 | 0.9999 | 1.0000 | 99.4% |
| 65 | 64.67 | 10.4 | 5.6 | strong/weak | 0 | 0.5925 | 0.2967 | 0.8% |
| 65 | 64.67 | 10.4 | 5.6 | strong/weak + 1 swap | 1 | 0.7386 | 0.4066 | 7.6% |
| 65 | 64.67 | 10.4 | 5.6 | strong/weak + 2 swaps | 2 | 0.8837 | 0.5467 | 48.2% |
| 70 | 70.25 | 11.2 | 4.8 | T=0.05N | 0 | 0.8635 | 0.5540 | 31.4% |
| 70 | 70.25 | 11.2 | 4.8 | T=0.05N + 1 swap | 1 | 0.9964 | 0.9616 | 92.8% |
| 70 | 70.25 | 11.2 | 4.8 | T=0.05N + 2 swaps | 2 | 0.9977 | 1.0000 | 98.4% |
| 70 | 70.25 | 11.2 | 4.8 | cyclic best | 0 | 0.9960 | 0.9645 | 75.0% |
| 70 | 70.25 | 11.2 | 4.8 | cyclic best + 1 swap | 1 | 0.9998 | 1.0000 | 98.8% |
| 70 | 70.25 | 11.2 | 4.8 | cyclic best + 2 swaps | 2 | 0.9999 | 1.0000 | 99.4% |
| 70 | 70.25 | 11.2 | 4.8 | strong/weak | 0 | 0.6344 | 0.3242 | 1.6% |
| 70 | 70.25 | 11.2 | 4.8 | strong/weak + 1 swap | 1 | 0.8092 | 0.4578 | 29.4% |
| 70 | 70.25 | 11.2 | 4.8 | strong/weak + 2 swaps | 2 | 0.9070 | 0.5668 | 53.8% |
| 75 | 75.00 | 12.0 | 4.0 | T=0.05N | 0 | 0.8618 | 0.5384 | 29.6% |
| 75 | 75.00 | 12.0 | 4.0 | T=0.05N + 1 swap | 1 | 0.9967 | 0.9635 | 93.6% |
| 75 | 75.00 | 12.0 | 4.0 | T=0.05N + 2 swaps | 2 | 0.9977 | 0.9997 | 98.0% |
| 75 | 75.00 | 12.0 | 4.0 | cyclic best | 0 | 0.9968 | 0.9772 | 76.4% |
| 75 | 75.00 | 12.0 | 4.0 | cyclic best + 1 swap | 1 | 0.9999 | 1.0000 | 98.8% |
| 75 | 75.00 | 12.0 | 4.0 | cyclic best + 2 swaps | 2 | 1.0000 | 1.0000 | 99.6% |
| 75 | 75.00 | 12.0 | 4.0 | strong/weak | 0 | 0.6673 | 0.3738 | 3.4% |
| 75 | 75.00 | 12.0 | 4.0 | strong/weak + 1 swap | 1 | 0.8531 | 0.5299 | 36.0% |
| 75 | 75.00 | 12.0 | 4.0 | strong/weak + 2 swaps | 2 | 0.9533 | 0.6530 | 72.6% |

## Notes

- Exact enumeration completed for `100.0%` of cases.
- Exact rows loaded from previous artifacts for `56.4%` of cases; missing active-K cases were recomputed.
- Strong/weak is one non-cyclic sorted row-power window: it disables weakest and strongest tails and keeps the middle rows.
- The cyclic best seed includes the strong/weak window as one possible cyclic or non-cyclic window, but chooses by tested `U_G`.

## Plots

![Raw U_G CDF](active_k_raw_u_g_cdf_by_requested_active_pct.png)

![Fraction exact CDF](active_k_fraction_exact_cdf.png)

![Mean fraction by requested active percentage](active_k_mean_fraction_by_requested_active_pct.png)

![Exact recovery by requested active percentage](active_k_exact_recovery_by_requested_active_pct.png)

![Best cyclic T boxplot](active_k_best_cyclic_T_boxplot.png)

![Best cyclic T over N](active_k_best_cyclic_T_over_N.png)

![Failure diagnostics](real_off_failure_diagnostics.png)

![Runtime](active_k_runtime_by_method.png)

## Artifacts

- Detailed CSVs are packed in `csv_data.tar.gz` after the run.
- Main report: `local_threshold_active_k_report.md`.
