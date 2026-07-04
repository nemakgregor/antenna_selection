# Sigma Sweep

> Historical K semantics note: this report uses active-K semantics. Here `K` is the number of selected/kept antennas, not the number turned off. A `25% active` or `K=0.25N` case means `75% off`, not the real `25% off` task. For real off-percent experiments, `25% off => K_active=0.75N` and `50% off => K_active=0.50N`.

- N: 1000
- L: 4
- K cases: K<=250 (25% active), K<=500 (50% active)
- Samples: 30
- Seed range: 42..71
- Sigma values: 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000
- Algorithms: H1, H2, Coutino, MISO-EE, Pareto-H2, H3-threshold-BF, H3-threshold-Int, H3-threshold-Gen, H3-Fast

Plots show raw mean objective values only: `U_BF`, `U_I`, and `U_G`.
`U_G = det(V_eq V_eq* + sigma I)` is the general objective.
If `lambda_i` are eigenvalues of `V_eq V_eq*`, then `U_G = product_i(lambda_i + sigma)`.
For `L=4`, `U_G = sigma^4 + sigma^3 sum(lambda_i) + sigma^2 e2(lambda) + sigma e3(lambda) + product(lambda_i)`.
The first algorithm-dependent term at large sigma is `sigma^3 sum(lambda_i)`, i.e. it tracks `U_BF = trace(V_eq V_eq*)`.
`U_BF` and `U_I` formulas do not contain `sigma`; if their curves move across sigma, the selected antenna set changed.

## K<=250 (25% active)

Mean `U_G` leader segments: 0.001..100000: H3-threshold-Gen.

| sigma | U_G mean leader | U_G mean | BF mean leader | BF mean | Interference mean leader | Interference mean |
|---:|:---|---:|:---|---:|:---|---:|
| 0.001 | H3-threshold-Gen | 2.8820e+12 | H3-threshold-BF | 5.2322e+03 | H3-Fast | 1.2177e+01 |
| 0.003 | H3-threshold-Gen | 2.8820e+12 | H3-threshold-BF | 5.2322e+03 | H3-Fast | 1.2177e+01 |
| 0.01 | H3-threshold-Gen | 2.8821e+12 | H3-threshold-BF | 5.2322e+03 | H3-Fast | 1.2177e+01 |
| 0.03 | H3-threshold-Gen | 2.8823e+12 | H3-threshold-BF | 5.2322e+03 | H3-Fast | 1.2177e+01 |
| 0.1 | H3-threshold-Gen | 2.8829e+12 | H3-threshold-BF | 5.2322e+03 | H3-Fast | 1.2177e+01 |
| 0.3 | H3-threshold-Gen | 2.8847e+12 | H3-threshold-BF | 5.2322e+03 | H3-Fast | 1.2177e+01 |
| 1 | H3-threshold-Gen | 2.8909e+12 | H3-threshold-BF | 5.2322e+03 | H3-Fast | 1.2177e+01 |
| 3 | H3-threshold-Gen | 2.9087e+12 | H3-threshold-BF | 5.2322e+03 | H3-Fast | 1.2177e+01 |
| 10 | H3-threshold-Gen | 2.9718e+12 | H3-threshold-BF | 5.2322e+03 | H3-Fast | 1.2177e+01 |
| 30 | H3-threshold-Gen | 3.1576e+12 | H3-threshold-BF | 5.2322e+03 | H3-Fast | 1.2177e+01 |
| 100 | H3-threshold-Gen | 3.8772e+12 | H3-threshold-BF | 5.2322e+03 | H3-Fast | 1.2177e+01 |
| 300 | H3-threshold-Gen | 6.6152e+12 | H3-threshold-BF | 5.2322e+03 | H3-Fast | 1.2177e+01 |
| 1000 | H3-threshold-Gen | 2.8225e+13 | H3-threshold-BF | 5.2322e+03 | H3-Fast | 1.2177e+01 |
| 3000 | H3-threshold-Gen | 3.4390e+14 | H3-threshold-BF | 5.2322e+03 | H3-Fast | 1.2177e+01 |
| 10000 | H3-threshold-Gen | 1.6347e+16 | H3-threshold-BF | 5.2322e+03 | H3-Fast | 1.2177e+01 |
| 30000 | H3-threshold-Gen | 9.6075e+17 | H3-threshold-BF | 5.2322e+03 | H3-Fast | 1.2177e+01 |
| 100000 | H3-threshold-Gen | 1.0534e+20 | H3-threshold-BF | 5.2322e+03 | H3-Fast | 1.2177e+01 |

### Raw Endpoint Comparison

| algorithm | BF at min sigma | BF at max sigma | Interference at min sigma | Interference at max sigma | U_G at min sigma | U_G at max sigma | unique selected sets/sample |
|:---|---:|---:|---:|---:|---:|---:|---:|
| H1 | 3.8802e+03 | 3.8802e+03 | 1.5219e+05 | 1.5219e+05 | 9.9360e+11 | 1.0394e+20 | 1.00 (1..1) |
| H2 | 3.4866e+03 | 3.4866e+03 | 2.0861e+02 | 2.0861e+02 | 6.8442e+11 | 1.0353e+20 | 1.00 (1..1) |
| Coutino | 4.6571e+03 | 4.7327e+03 | 1.0280e+05 | 6.4242e+05 | 1.9318e+12 | 1.0482e+20 | 7.87 (6..10) |
| MISO-EE | 3.9352e+03 | 4.7327e+03 | 8.0232e+04 | 6.4242e+05 | 1.0790e+12 | 1.0482e+20 | 8.33 (7..11) |
| Pareto-H2 | 3.6209e+03 | 3.6209e+03 | 1.5827e+02 | 1.5827e+02 | 7.6521e+11 | 1.0367e+20 | 1.00 (1..1) |
| H3-threshold-BF | 5.2322e+03 | 5.2322e+03 | 2.3114e+05 | 2.3114e+05 | 2.8726e+12 | 1.0534e+20 | 1.00 (1..1) |
| H3-threshold-Int | 3.6938e+03 | 3.6938e+03 | 1.7984e+02 | 1.7984e+02 | 8.7129e+11 | 1.0375e+20 | 1.00 (1..1) |
| H3-threshold-Gen | 5.2303e+03 | 5.2322e+03 | 2.1612e+05 | 2.2602e+05 | 2.8820e+12 | 1.0534e+20 | 1.10 (1..2) |
| H3-Fast | 1.5948e+03 | 1.5948e+03 | 1.2177e+01 | 1.2177e+01 | 2.8650e+10 | 1.0160e+20 | 1.00 (1..1) |

### Selection Stability

The count below is the number of distinct selected antenna sets across the sigma grid for the same random sample.

| algorithm | mean unique sets/sample | min | max |
|:---|---:|---:|---:|
| Coutino | 7.87 | 6 | 10 |
| H1 | 1.00 | 1 | 1 |
| H2 | 1.00 | 1 | 1 |
| H3-Fast | 1.00 | 1 | 1 |
| H3-threshold-BF | 1.00 | 1 | 1 |
| H3-threshold-Gen | 1.10 | 1 | 2 |
| H3-threshold-Int | 1.00 | 1 | 1 |
| MISO-EE | 8.33 | 7 | 11 |
| Pareto-H2 | 1.00 | 1 | 1 |

## K<=500 (50% active)

Mean `U_G` leader segments: 0.001..3000: H3-threshold-Gen, 10000..100000: H3-threshold-BF.

| sigma | U_G mean leader | U_G mean | BF mean leader | BF mean | Interference mean leader | Interference mean |
|---:|:---|---:|:---|---:|:---|---:|
| 0.001 | H3-threshold-Gen | 1.9658e+14 | H3-threshold-BF | 1.4887e+04 | H3-Fast | 4.7202e+01 |
| 0.003 | H3-threshold-Gen | 1.9658e+14 | H3-threshold-BF | 1.4887e+04 | H3-Fast | 4.7202e+01 |
| 0.01 | H3-threshold-Gen | 1.9658e+14 | H3-threshold-BF | 1.4887e+04 | H3-Fast | 4.7202e+01 |
| 0.03 | H3-threshold-Gen | 1.9658e+14 | H3-threshold-BF | 1.4887e+04 | H3-Fast | 4.7202e+01 |
| 0.1 | H3-threshold-Gen | 1.9660e+14 | H3-threshold-BF | 1.4887e+04 | H3-Fast | 4.7202e+01 |
| 0.3 | H3-threshold-Gen | 1.9664e+14 | H3-threshold-BF | 1.4887e+04 | H3-Fast | 4.7202e+01 |
| 1 | H3-threshold-Gen | 1.9678e+14 | H3-threshold-BF | 1.4887e+04 | H3-Fast | 4.7202e+01 |
| 3 | H3-threshold-Gen | 1.9720e+14 | H3-threshold-BF | 1.4887e+04 | H3-Fast | 4.7202e+01 |
| 10 | H3-threshold-Gen | 1.9867e+14 | H3-threshold-BF | 1.4887e+04 | H3-Fast | 4.7202e+01 |
| 30 | H3-threshold-Gen | 2.0291e+14 | H3-threshold-BF | 1.4887e+04 | H3-Fast | 4.7202e+01 |
| 100 | H3-threshold-Gen | 2.1829e+14 | H3-threshold-BF | 1.4887e+04 | H3-Fast | 4.7202e+01 |
| 300 | H3-threshold-Gen | 2.6710e+14 | H3-threshold-BF | 1.4887e+04 | H3-Fast | 4.7202e+01 |
| 1000 | H3-threshold-Gen | 5.0457e+14 | H3-threshold-BF | 1.4887e+04 | H3-Fast | 4.7202e+01 |
| 3000 | H3-threshold-Gen | 2.0564e+15 | H3-threshold-BF | 1.4887e+04 | H3-Fast | 4.7202e+01 |
| 10000 | H3-threshold-BF | 3.5513e+16 | H3-threshold-BF | 1.4887e+04 | H3-Fast | 4.7202e+01 |
| 30000 | H3-threshold-BF | 1.2935e+18 | H3-threshold-BF | 1.4887e+04 | H3-Fast | 4.7202e+01 |
| 100000 | H3-threshold-BF | 1.1574e+20 | H3-threshold-BF | 1.4887e+04 | H3-Fast | 4.7202e+01 |

### Raw Endpoint Comparison

| algorithm | BF at min sigma | BF at max sigma | Interference at min sigma | Interference at max sigma | U_G at min sigma | U_G at max sigma | unique selected sets/sample |
|:---|---:|---:|---:|---:|---:|---:|---:|
| H1 | 1.0477e+04 | 1.0477e+04 | 5.8887e+05 | 5.8887e+05 | 5.5777e+13 | 1.1091e+20 | 1.00 (1..1) |
| H2 | 9.6023e+03 | 9.6023e+03 | 2.5122e+02 | 2.5122e+02 | 3.9996e+13 | 1.0997e+20 | 1.00 (1..1) |
| Coutino | 1.2764e+04 | 1.2814e+04 | 5.5114e+05 | 1.2663e+06 | 1.1139e+14 | 1.1345e+20 | 7.23 (5..9) |
| MISO-EE | 1.0713e+04 | 1.2814e+04 | 4.1463e+05 | 1.2663e+06 | 6.0644e+13 | 1.1345e+20 | 7.77 (7..11) |
| Pareto-H2 | 1.0216e+04 | 1.0216e+04 | 2.0232e+02 | 2.0232e+02 | 4.8496e+13 | 1.1062e+20 | 1.00 (1..1) |
| H3-threshold-BF | 1.4887e+04 | 1.4887e+04 | 1.0699e+06 | 1.0699e+06 | 1.9646e+14 | 1.1574e+20 | 1.00 (1..1) |
| H3-threshold-Int | 9.9650e+03 | 9.9650e+03 | 2.2173e+02 | 2.2173e+02 | 4.7516e+13 | 1.1036e+20 | 1.00 (1..1) |
| H3-threshold-Gen | 1.4885e+04 | 1.4887e+04 | 1.0354e+06 | 1.0699e+06 | 1.9658e+14 | 1.1574e+20 | 1.07 (1..2) |
| H3-Fast | 5.7734e+03 | 5.7734e+03 | 4.7202e+01 | 4.7202e+01 | 5.1633e+12 | 1.0590e+20 | 1.00 (1..1) |

### Selection Stability

The count below is the number of distinct selected antenna sets across the sigma grid for the same random sample.

| algorithm | mean unique sets/sample | min | max |
|:---|---:|---:|---:|
| Coutino | 7.23 | 5 | 9 |
| H1 | 1.00 | 1 | 1 |
| H2 | 1.00 | 1 | 1 |
| H3-Fast | 1.00 | 1 | 1 |
| H3-threshold-BF | 1.00 | 1 | 1 |
| H3-threshold-Gen | 1.07 | 1 | 2 |
| H3-threshold-Int | 1.00 | 1 | 1 |
| MISO-EE | 7.77 | 7 | 11 |
| Pareto-H2 | 1.00 | 1 | 1 |

Generated plots: `sigma_sweep_K*.png` for metric curves and `sigma_winners_K*.png` for per-sample `U_G` winner shares.