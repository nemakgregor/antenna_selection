# Grid Benchmark

- N values: [1000]
- L values: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
- off percentages: [25, 50]
- samples per case: 10
- algorithms: ['H1', 'H2', 'Coutino', 'MISO-EE', 'Pareto-H2', 'S-threshold-BF', 'S-threshold-Int', 'S-threshold-Gen', 'Frame-BF', 'Frame-Int', 'Frame-Gen', 'N-H3-Fast']

K_active is computed as round(N * (1 - off_pct / 100)).
All algorithm entries call the corresponding implementation directly.

## Winners By Mean

### 25% off, N=1000, K_active=750
- BF gain: Frame-Gen (22661.1)
- Interference: Frame-Int (343.997)
- General objective: Frame-Gen (1.30948e+31)
- Fastest avg runtime: H1 (0.0290422s)
- Energy-efficiency proxy: MISO-EE (0.0944769 log2(U_G)/active)

### 50% off, N=1000, K_active=500
- BF gain: Frame-Gen (13900.6)
- Interference: Frame-Int (152.803)
- General objective: Frame-Gen (2.86472e+28)
- Fastest avg runtime: H1 (0.0192519s)
- Energy-efficiency proxy: MISO-EE (0.128748 log2(U_G)/active)
