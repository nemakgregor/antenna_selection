# Sigma Sweep

- N: 1000
- L: 4
- K cases: K<=500 (50% off, 50% active), K<=750 (25% off, 75% active)
- Samples: 10
- Seed range: 42..51
- Sigma values: 1e-10, 1e-09, 1e-08, 1e-07, 1e-06, 1e-05, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000
- Algorithms: H1, H2, Coutino, MISO-EE, Pareto-H2, H3-threshold-BF, H3-threshold-Int, H3-threshold-Gen, H3-Fast

Plots show raw mean objective values only: `U_BF`, `U_I`, and `U_G`.
`U_G = det(V_eq V_eq* + sigma I)` is the general objective.
If `lambda_i` are eigenvalues of `V_eq V_eq*`, then `U_G = product_i(lambda_i + sigma)`.
For `L=4`, `U_G = sigma^4 + sigma^3 sum(lambda_i) + sigma^2 e2(lambda) + sigma e3(lambda) + product(lambda_i)`.
The first algorithm-dependent term at large sigma is `sigma^3 sum(lambda_i)`, i.e. it tracks `U_BF = trace(V_eq V_eq*)`.
`U_BF` and `U_I` formulas do not contain `sigma`; if their curves move across sigma, the selected antenna set changed.

## K<=500 (50% off, 50% active)

Mean `U_G` leader segments: 1e-10..1000: H3-threshold-BF.

| sigma | U_G mean leader | U_G mean | BF mean leader | BF mean | Interference mean leader | Interference mean |
|---:|:---|---:|:---|---:|:---|---:|
| 1e-10 | H3-threshold-BF | 2.1951e+14 | H3-threshold-BF | 1.5164e+04 | H3-Fast | 4.4145e+01 |
| 1e-09 | H3-threshold-BF | 2.1951e+14 | H3-threshold-BF | 1.5164e+04 | H3-Fast | 4.4145e+01 |
| 1e-08 | H3-threshold-BF | 2.1951e+14 | H3-threshold-BF | 1.5164e+04 | H3-Fast | 4.4145e+01 |
| 1e-07 | H3-threshold-BF | 2.1951e+14 | H3-threshold-BF | 1.5164e+04 | H3-Fast | 4.4145e+01 |
| 1e-06 | H3-threshold-BF | 2.1951e+14 | H3-threshold-BF | 1.5164e+04 | H3-Fast | 4.4145e+01 |
| 1e-05 | H3-threshold-BF | 2.1951e+14 | H3-threshold-BF | 1.5164e+04 | H3-Fast | 4.4145e+01 |
| 0.0001 | H3-threshold-BF | 2.1951e+14 | H3-threshold-BF | 1.5164e+04 | H3-Fast | 4.4145e+01 |
| 0.001 | H3-threshold-BF | 2.1951e+14 | H3-threshold-BF | 1.5164e+04 | H3-Fast | 4.4145e+01 |
| 0.01 | H3-threshold-BF | 2.1951e+14 | H3-threshold-BF | 1.5164e+04 | H3-Fast | 4.4145e+01 |
| 0.1 | H3-threshold-BF | 2.1953e+14 | H3-threshold-BF | 1.5164e+04 | H3-Fast | 4.4145e+01 |
| 1 | H3-threshold-BF | 2.1973e+14 | H3-threshold-BF | 1.5164e+04 | H3-Fast | 4.4145e+01 |
| 10 | H3-threshold-BF | 2.2176e+14 | H3-threshold-BF | 1.5164e+04 | H3-Fast | 4.4145e+01 |
| 100 | H3-threshold-BF | 2.4286e+14 | H3-threshold-BF | 1.5164e+04 | H3-Fast | 4.4145e+01 |
| 1000 | H3-threshold-BF | 5.4749e+14 | H3-threshold-BF | 1.5164e+04 | H3-Fast | 4.4145e+01 |

### Raw Endpoint Comparison

| algorithm | BF at min sigma | BF at max sigma | Interference at min sigma | Interference at max sigma | U_G at min sigma | U_G at max sigma | unique selected sets/sample |
|:---|---:|---:|---:|---:|---:|---:|---:|
| H1 | 1.0876e+04 | 1.0876e+04 | 6.2636e+05 | 6.2636e+05 | 6.9756e+13 | 2.1915e+14 | 1.00 (1..1) |
| H2 | 9.8931e+03 | 9.8931e+03 | 2.6475e+02 | 2.6475e+02 | 4.7387e+13 | 1.6491e+14 | 1.00 (1..1) |
| Coutino | 1.3578e+04 | 1.3594e+04 | 6.1652e+05 | 7.3515e+05 | 1.4231e+14 | 3.8895e+14 | 2.60 (2..3) |
| MISO-EE | 1.1137e+04 | 1.1214e+04 | 4.3577e+05 | 5.3991e+05 | 7.5484e+13 | 2.3688e+14 | 2.70 (2..3) |
| Pareto-H2 | 1.0674e+04 | 1.0674e+04 | 1.8824e+02 | 1.8824e+02 | 6.0813e+13 | 1.9974e+14 | 1.00 (1..1) |
| H3-threshold-BF | 1.5164e+04 | 1.5164e+04 | 9.9552e+05 | 9.9552e+05 | 2.1951e+14 | 5.4749e+14 | 1.00 (1..1) |
| H3-threshold-Int | 1.0187e+04 | 1.0187e+04 | 2.2300e+02 | 2.2300e+02 | 5.5580e+13 | 1.8350e+14 | 1.00 (1..1) |
| H3-threshold-Gen | 1.5164e+04 | 1.5164e+04 | 9.9552e+05 | 9.9552e+05 | 2.1951e+14 | 5.4749e+14 | 1.00 (1..1) |
| H3-Fast | 6.0974e+03 | 6.0974e+03 | 4.4145e+01 | 4.4145e+01 | 6.7184e+12 | 4.4117e+13 | 1.00 (1..1) |

### Selection Stability

The count below is the number of distinct selected antenna sets across the sigma grid for the same random sample.

| algorithm | mean unique sets/sample | min | max |
|:---|---:|---:|---:|
| Coutino | 2.60 | 2 | 3 |
| H1 | 1.00 | 1 | 1 |
| H2 | 1.00 | 1 | 1 |
| H3-Fast | 1.00 | 1 | 1 |
| H3-threshold-BF | 1.00 | 1 | 1 |
| H3-threshold-Gen | 1.00 | 1 | 1 |
| H3-threshold-Int | 1.00 | 1 | 1 |
| MISO-EE | 2.70 | 2 | 3 |
| Pareto-H2 | 1.00 | 1 | 1 |

## K<=750 (25% off, 75% active)

Mean `U_G` leader segments: 1e-10..1000: H3-threshold-BF.

| sigma | U_G mean leader | U_G mean | BF mean leader | BF mean | Interference mean leader | Interference mean |
|---:|:---|---:|:---|---:|:---|---:|
| 1e-10 | H3-threshold-BF | 1.6966e+15 | H3-threshold-BF | 2.5207e+04 | H3-Fast | 8.2938e+01 |
| 1e-09 | H3-threshold-BF | 1.6966e+15 | H3-threshold-BF | 2.5207e+04 | H3-Fast | 8.2938e+01 |
| 1e-08 | H3-threshold-BF | 1.6966e+15 | H3-threshold-BF | 2.5207e+04 | H3-Fast | 8.2938e+01 |
| 1e-07 | H3-threshold-BF | 1.6966e+15 | H3-threshold-BF | 2.5207e+04 | H3-Fast | 8.2938e+01 |
| 1e-06 | H3-threshold-BF | 1.6966e+15 | H3-threshold-BF | 2.5207e+04 | H3-Fast | 8.2938e+01 |
| 1e-05 | H3-threshold-BF | 1.6966e+15 | H3-threshold-BF | 2.5207e+04 | H3-Fast | 8.2938e+01 |
| 0.0001 | H3-threshold-BF | 1.6966e+15 | H3-threshold-BF | 2.5207e+04 | H3-Fast | 8.2938e+01 |
| 0.001 | H3-threshold-BF | 1.6967e+15 | H3-threshold-BF | 2.5207e+04 | H3-Fast | 8.2938e+01 |
| 0.01 | H3-threshold-BF | 1.6967e+15 | H3-threshold-BF | 2.5207e+04 | H3-Fast | 8.2938e+01 |
| 0.1 | H3-threshold-BF | 1.6968e+15 | H3-threshold-BF | 2.5207e+04 | H3-Fast | 8.2938e+01 |
| 1 | H3-threshold-BF | 1.6977e+15 | H3-threshold-BF | 2.5207e+04 | H3-Fast | 8.2938e+01 |
| 10 | H3-threshold-BF | 1.7071e+15 | H3-threshold-BF | 2.5207e+04 | H3-Fast | 8.2938e+01 |
| 100 | H3-threshold-BF | 1.8030e+15 | H3-threshold-BF | 2.5207e+04 | H3-Fast | 8.2938e+01 |
| 1000 | H3-threshold-BF | 3.0030e+15 | H3-threshold-BF | 2.5207e+04 | H3-Fast | 8.2938e+01 |

### Raw Endpoint Comparison

| algorithm | BF at min sigma | BF at max sigma | Interference at min sigma | Interference at max sigma | U_G at min sigma | U_G at max sigma | unique selected sets/sample |
|:---|---:|---:|---:|---:|---:|---:|---:|
| H1 | 1.7857e+04 | 1.7857e+04 | 1.2325e+06 | 1.2325e+06 | 5.1532e+14 | 1.0668e+15 | 1.00 (1..1) |
| H2 | 1.6382e+04 | 1.6382e+04 | 1.9643e+02 | 1.9643e+02 | 3.6061e+14 | 7.9566e+14 | 1.00 (1..1) |
| Coutino | 2.2410e+04 | 2.2415e+04 | 1.5812e+06 | 1.6499e+06 | 1.0610e+15 | 2.0051e+15 | 2.10 (1..3) |
| MISO-EE | 1.8360e+04 | 1.8409e+04 | 1.1342e+06 | 1.2319e+06 | 5.5854e+14 | 1.1531e+15 | 2.50 (2..3) |
| Pareto-H2 | 1.7639e+04 | 1.7639e+04 | 1.3536e+02 | 1.3536e+02 | 4.6067e+14 | 9.7894e+14 | 1.00 (1..1) |
| H3-threshold-BF | 2.5207e+04 | 2.5207e+04 | 2.1372e+06 | 2.1372e+06 | 1.6966e+15 | 3.0030e+15 | 1.00 (1..1) |
| H3-threshold-Int | 1.7268e+04 | 1.7268e+04 | 1.6003e+02 | 1.6003e+02 | 5.1554e+14 | 1.0470e+15 | 1.00 (1..1) |
| H3-threshold-Gen | 2.5207e+04 | 2.5207e+04 | 2.1372e+06 | 2.1372e+06 | 1.6966e+15 | 3.0030e+15 | 1.00 (1..1) |
| H3-Fast | 1.3083e+04 | 1.3083e+04 | 8.2938e+01 | 8.2938e+01 | 1.5226e+14 | 3.9705e+14 | 1.00 (1..1) |

### Selection Stability

The count below is the number of distinct selected antenna sets across the sigma grid for the same random sample.

| algorithm | mean unique sets/sample | min | max |
|:---|---:|---:|---:|
| Coutino | 2.10 | 1 | 3 |
| H1 | 1.00 | 1 | 1 |
| H2 | 1.00 | 1 | 1 |
| H3-Fast | 1.00 | 1 | 1 |
| H3-threshold-BF | 1.00 | 1 | 1 |
| H3-threshold-Gen | 1.00 | 1 | 1 |
| H3-threshold-Int | 1.00 | 1 | 1 |
| MISO-EE | 2.50 | 2 | 3 |
| Pareto-H2 | 1.00 | 1 | 1 |

Generated plots: `sigma_sweep_K*.png` for metric curves and `sigma_winners_K*.png` for per-sample `U_G` winner shares.