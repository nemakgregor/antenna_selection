# H1/H2/Coutino/MISO-EE/Pareto-H2 Objective Values

- N: 1000
- L values: 2, 3, 4, 5, 6, 7, 8, 9, 10
- Active limits: 0.75, 0.5
- Seeds: 42
- Sigma: 1.0
- P: 1.0

Objectives: `u_bf` and `u_g` are maximized; `u_i` is minimized.
The constraint is `active_count <= K`, so algorithms may turn off more than the required minimum.
`H1`, `H2`, and `Coutino` scan their feasible deletion paths and keep the best `U_G` set with `active_count <= K`.
`MISO-EE` returns the smallest greedy log-det set that still targets +5% `U_G` over best(H1,H2).
`Pareto-H2` minimizes interference on a BF-protected deletion path; lower `Pareto Int / H2 Int` is better, while `Pareto BF / H2 BF >= 1` means no BF loss versus H2.
Energy-efficiency proxy: `log2(U_G) / active antenna`; larger is better.

Recommendation: for Task 3 (`sigma=1`, maximize `U_G`), use `Coutino`; in the default `N=1000`, `L=2..10` run it wins `U_G` for every 25%+ and 50%+ off case.
Use `MISO-EE` only if the priority is energy efficiency with extra antennas switched off while keeping about +5% `U_G` over best(H1,H2).
Use `Pareto-H2` only for the interference/BF compromise; it is deliberately conservative and falls back to H2 when the protected path is not close enough on interference.

## 25% Antennas Off

| L | K max | Coutino active | MISO-EE active | MISO-EE off | Pareto active | Pareto BF / H2 BF | Pareto Int / H2 Int | H1 Gen | H2 Gen | Coutino Gen | MISO-EE Gen | Coutino Gen vs best H1/H2 | MISO-EE Gen vs best H1/H2 | Int winner | EE proxy winner |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|
| 2 | 750 | 750 | 709 | 29.1% | 750 | 1.000 | 1.000 | 1.1911e+08 | 1.1168e+08 | 1.3706e+08 | 1.2525e+08 | +15.06% | +5.15% | Tie | MISO-EE |
| 3 | 750 | 750 | 587 | 41.3% | 750 | 1.000 | 1.000 | 8.1769e+10 | 1.4265e+11 | 3.3562e+11 | 1.4990e+11 | +135.28% | +5.08% | Tie | MISO-EE |
| 4 | 750 | 750 | 696 | 30.4% | 737 | 1.063 | 0.659 | 1.5786e+15 | 1.0656e+15 | 2.3178e+15 | 1.6629e+15 | +46.83% | +5.34% | Pareto-H2 | MISO-EE |
| 5 | 750 | 750 | 687 | 31.3% | 750 | 1.002 | 0.806 | 3.6341e+17 | 2.5868e+17 | 6.4738e+17 | 3.8328e+17 | +78.14% | +5.47% | Pareto-H2 | MISO-EE |
| 6 | 750 | 750 | 693 | 30.7% | 719 | 1.069 | 1.283 | 4.0675e+20 | 1.5108e+20 | 7.6997e+20 | 4.3104e+20 | +89.30% | +5.97% | H2 | MISO-EE |
| 7 | 750 | 750 | 711 | 28.9% | 750 | 1.000 | 1.000 | 3.5595e+23 | 1.1719e+23 | 6.0156e+23 | 3.7610e+23 | +69.00% | +5.66% | Tie | MISO-EE |
| 8 | 750 | 750 | 726 | 27.4% | 750 | 1.000 | 0.977 | 4.2150e+26 | 1.3549e+26 | 6.2939e+26 | 4.4406e+26 | +49.32% | +5.35% | Pareto-H2 | MISO-EE |
| 9 | 750 | 750 | 643 | 35.7% | 750 | 1.000 | 1.000 | 1.7314e+28 | 3.7844e+27 | 1.3141e+29 | 1.8229e+28 | +658.96% | +5.28% | Tie | MISO-EE |
| 10 | 750 | 750 | 702 | 29.8% | 750 | 1.000 | 1.000 | 2.2444e+30 | 1.0957e+30 | 6.1260e+30 | 2.4000e+30 | +172.95% | +6.94% | Tie | MISO-EE |

## 50% Antennas Off

| L | K max | Coutino active | MISO-EE active | MISO-EE off | Pareto active | Pareto BF / H2 BF | Pareto Int / H2 Int | H1 Gen | H2 Gen | Coutino Gen | MISO-EE Gen | Coutino Gen vs best H1/H2 | MISO-EE Gen vs best H1/H2 | Int winner | EE proxy winner |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|
| 2 | 500 | 500 | 481 | 51.9% | 500 | 1.000 | 1.000 | 5.4233e+07 | 4.6770e+07 | 6.2301e+07 | 5.7055e+07 | +14.88% | +5.20% | Tie | MISO-EE |
| 3 | 500 | 500 | 406 | 59.4% | 500 | 1.000 | 1.000 | 2.0562e+10 | 3.4461e+10 | 8.2741e+10 | 3.6447e+10 | +140.10% | +5.76% | Tie | MISO-EE |
| 4 | 500 | 500 | 470 | 53.0% | 492 | 1.067 | 0.653 | 2.1760e+14 | 1.3944e+14 | 3.1871e+14 | 2.2882e+14 | +46.46% | +5.16% | Pareto-H2 | MISO-EE |
| 5 | 500 | 500 | 464 | 53.6% | 500 | 1.000 | 1.000 | 2.6964e+16 | 1.9614e+16 | 4.7818e+16 | 2.8487e+16 | +77.34% | +5.65% | Tie | MISO-EE |
| 6 | 500 | 500 | 468 | 53.2% | 483 | 1.068 | 0.814 | 1.6164e+19 | 6.4184e+18 | 3.0455e+19 | 1.7119e+19 | +88.41% | +5.91% | Pareto-H2 | MISO-EE |
| 7 | 500 | 500 | 478 | 52.2% | 500 | 1.000 | 1.000 | 7.2045e+21 | 2.6702e+21 | 1.2121e+22 | 7.6190e+21 | +68.24% | +5.75% | Tie | MISO-EE |
| 8 | 500 | 500 | 485 | 51.5% | 500 | 1.000 | 1.000 | 4.0742e+24 | 1.7769e+24 | 6.2095e+24 | 4.3147e+24 | +52.41% | +5.90% | Tie | MISO-EE |
| 9 | 500 | 500 | 435 | 56.5% | 500 | 1.000 | 1.000 | 7.8571e+25 | 2.2183e+25 | 5.8986e+26 | 8.2925e+25 | +650.73% | +5.54% | Tie | MISO-EE |
| 10 | 500 | 500 | 471 | 52.9% | 500 | 1.000 | 1.000 | 5.2487e+27 | 2.9980e+27 | 1.4338e+28 | 5.5778e+27 | +173.17% | +6.27% | Tie | MISO-EE |
