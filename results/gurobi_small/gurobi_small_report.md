# Small Gurobi Optimum

- N: 10
- L: 2
- K max: 5
- Minimum active antennas: 2
- Exact K only: False
- Objective: general
- BF floor: 0.85
- H2 interference ceiling factor: 2.0
- Enumerated subsets: 627
- Direct scan verified: True
- Optimum subset: 1 3 6 8 9

| algorithm | active | BF gain | interference | U_G | meets BF floor | meets int ceiling | objective gap |
|:---|---:|---:|---:|---:|:---:|:---:|---:|
| Gurobi optimum | 5 | 9.0356 | 2.6310 | 2.6929e+01 | True | True | 0.00% |
| H1 | 5 | 9.6118 | 12.3553 | 2.6491e+01 | True | True | 1.63% |
| H2 | 5 | 7.6727 | 0.1354 | 2.1723e+01 | True | True | 19.33% |
| Coutino | 5 | 9.0356 | 2.6310 | 2.6929e+01 | True | True | 0.00% |
| MISO-EE | 5 | 9.0356 | 2.6310 | 2.6929e+01 | True | True | 0.00% |
| Pareto-H2 | 5 | 7.6727 | 0.1354 | 2.1723e+01 | True | True | 19.33% |