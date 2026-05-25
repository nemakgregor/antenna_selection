# Sigma Sweep

- N: 1000
- L: 4
- K max: 500
- Minimum antennas off: 50%
- Seed: 42

Important: `sigma` appears only in `U_G = det(V_eq V_eq* + sigma I)`.
`U_BF` and `U_I` do not contain sigma, so they are expected to be flat unless the selected antenna set changes.
The report below also checks whether the selected antenna set changed across sigma.

| sigma | best U_G | best BF | best interference | best EE proxy | Coutino U_G / best H1H2 | Pareto BF / H2 BF | Pareto Int / H2 Int | MISO-EE active | MISO-EE off |
|---:|:---:|:---:|:---:|:---:|---:|---:|---:|---:|---:|
| 0.001 | Coutino | Coutino | Pareto-H2 | MISO-EE | 1.465 | 1.067 | 0.653 | 470 | 53.0% |
| 0.003 | Coutino | Coutino | Pareto-H2 | MISO-EE | 1.465 | 1.067 | 0.653 | 470 | 53.0% |
| 0.01 | Coutino | Coutino | Pareto-H2 | MISO-EE | 1.465 | 1.067 | 0.653 | 470 | 53.0% |
| 0.03 | Coutino | Coutino | Pareto-H2 | MISO-EE | 1.465 | 1.067 | 0.653 | 470 | 53.0% |
| 0.1 | Coutino | Coutino | Pareto-H2 | MISO-EE | 1.465 | 1.067 | 0.653 | 470 | 53.0% |
| 0.3 | Coutino | Coutino | Pareto-H2 | MISO-EE | 1.465 | 1.067 | 0.653 | 470 | 53.0% |
| 1 | Coutino | Coutino | Pareto-H2 | MISO-EE | 1.465 | 1.067 | 0.653 | 470 | 53.0% |
| 3 | Coutino | Coutino | Pareto-H2 | MISO-EE | 1.464 | 1.067 | 0.653 | 470 | 53.0% |
| 10 | Coutino | Coutino | Pareto-H2 | MISO-EE | 1.463 | 1.067 | 0.653 | 470 | 53.0% |
| 30 | Coutino | Coutino | Pareto-H2 | MISO-EE | 1.461 | 1.067 | 0.653 | 470 | 53.0% |
| 100 | Coutino | Coutino | Pareto-H2 | MISO-EE | 1.451 | 1.067 | 0.653 | 470 | 53.0% |

## U_G Sigma Sensitivity

| algorithm | U_G at min sigma | U_G at max sigma | relative change | max sigma / min eigenvalue |
|:---|---:|---:|---:|---:|
| H1 | 2.1737e+14 | 2.4112e+14 | 10.92% | 0.0307 |
| H2 | 1.3927e+14 | 1.5624e+14 | 12.18% | 0.0319 |
| Coutino | 3.1840e+14 | 3.4983e+14 | 9.87% | 0.0273 |
| MISO-EE | 2.2858e+14 | 2.5320e+14 | 10.77% | 0.0301 |
| Pareto-H2 | 1.7795e+14 | 1.9836e+14 | 11.47% | 0.0315 |

Selected-set check: unchanged for H1, H2, MISO-EE, Pareto-H2; changed for Coutino.

Interpretation: H2 is expected to win pure interference; Coutino usually wins raw `U_G`; MISO-EE is the energy-saving compromise; Pareto-H2 targets H2-like interference while protecting BF.