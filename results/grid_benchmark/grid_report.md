# Grid Benchmark

- N values: [1000, 5000, 10000]
- L values: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
- off percentages: [25, 50]
- samples per case: 50
- algorithms: ['H1', 'H2', 'Coutino', 'MISO-EE', 'Pareto-H2', 'S-threshold-BF', 'S-threshold-Int', 'S-threshold-Gen', 'N-H3-Fast']

K_active is computed as round(N * (1 - off_pct / 100)).
All algorithm entries call the corresponding implementation directly.

## Winners By Mean

### 25% off, N=1000, K_active=750
- BF gain: S-threshold-BF (22622.6)
- Interference: S-threshold-Int (475.971)
- General objective: S-threshold-Gen (2.63603e+31)
- Fastest avg runtime: H1 (0.0164776s)

### 25% off, N=5000, K_active=3750
- BF gain: S-threshold-BF (478183)
- Interference: Pareto-H2 (6688)
- General objective: S-threshold-BF (7.34121e+44)
- Fastest avg runtime: H1 (0.0776725s)

### 25% off, N=10000, K_active=7500
- BF gain: S-threshold-BF (1.82457e+06)
- Interference: Pareto-H2 (21375.4)
- General objective: S-threshold-BF (5.99801e+50)
- Fastest avg runtime: H1 (0.137354s)

### 50% off, N=1000, K_active=500
- BF gain: S-threshold-BF (13883.8)
- Interference: N-H3-Fast (176.059)
- General objective: S-threshold-Gen (5.66393e+28)
- Fastest avg runtime: H1 (0.00695334s)

### 50% off, N=5000, K_active=2500
- BF gain: S-threshold-BF (291797)
- Interference: N-H3-Fast (3254.05)
- General objective: S-threshold-BF (1.62462e+42)
- Fastest avg runtime: H1 (0.0390313s)

### 50% off, N=10000, K_active=5000
- BF gain: S-threshold-BF (1.11236e+06)
- Interference: N-H3-Fast (12357.1)
- General objective: S-threshold-BF (1.34362e+48)
- Fastest avg runtime: H1 (0.0920151s)
