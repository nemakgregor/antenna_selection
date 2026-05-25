# Small Gurobi Optimum

- N: 18
- L: 4
- K max: 9
- Minimum active antennas: 4
- Exact K only: False
- Objective: general
- BF floor: 0.85
- H2 interference ceiling factor: 2.0
- Enumerated subsets: 154394
- Direct scan verified: True
- Optimum subset: 3 7 9 11 12 13 15 16 17
- Algorithms: Gurobi optimum, H1, H2, Coutino, MISO-EE, Pareto-H2, H3-threshold-BF, H3-threshold-Int, H3-threshold-Gen, H3-Fast

## Timing

- Enumeration: 5.602428s
- Gurobi model build/solve: 1.364855s
- Direct scan verification: 0.056316s
- Gurobi total compared below: 6.967284s

| algorithm | active | BF gain | interference | U_G | time, s | meets BF floor | meets int ceiling | objective gap |
|:---|---:|---:|---:|---:|---:|:---:|:---:|---:|
| Gurobi optimum | 9 | 15.1666 | 19.9859 | 3.0707e+02 | 6.967284 | True | True | 0.00% |
| H1 | 9 | 16.1412 | 36.9397 | 2.5174e+02 | 0.000310 | True | True | 18.02% |
| H2 | 9 | 9.6541 | 3.4007 | 1.1587e+02 | 0.000475 | True | True | 62.27% |
| Coutino | 9 | 15.1666 | 19.9859 | 3.0707e+02 | 0.000702 | True | True | 0.00% |
| MISO-EE | 9 | 15.1666 | 19.9859 | 3.0707e+02 | 0.001376 | True | True | 0.00% |
| Pareto-H2 | 9 | 9.6541 | 3.4007 | 1.1587e+02 | 0.001387 | True | True | 62.27% |
| H3-threshold-BF | 9 | 17.2596 | 43.7568 | 3.0272e+02 | 0.000313 | True | True | 1.41% |
| H3-threshold-Int | 9 | 9.6541 | 3.4007 | 1.1587e+02 | 0.000985 | True | True | 62.27% |
| H3-threshold-Gen | 9 | 17.2596 | 43.7568 | 3.0272e+02 | 0.001728 | True | True | 1.41% |
| H3-Fast | 9 | 8.7826 | 3.5191 | 5.3973e+01 | 0.000735 | True | True | 82.42% |