# Frame vs FrameOnly Comparison

Full grid: N=1000, L=1..10, off_pct in {25,50}, samples=10.

- BF: FrameOnly is identical to Frame on quality in this grid; mean time ratio 0.985x.
- GEN: FrameOnly is identical to Frame on quality in this grid; mean time ratio 0.997x.
- INT: FrameOnly-Int mean time ratio 0.090x, median time ratio 0.068x.
- INT: direct mean runtime 0.1045s vs 1.6417s for Frame-Int.
- INT: excluding L=1 zero-interference ties, FrameOnly-Int U_I / Frame-Int U_I median 4.179x, mean 10.036x, min 1.559x, max 76.445x. Lower is better, so the pure variant is faster but weaker on interference.

Mean winners from grid_summary:
- BF: Frame-Gen 12/20, Frame-BF 8/20
- Interference: Frame-Int 18/20, Coutino 2/20
- General objective: Frame-Gen 17/20, Frame-BF 3/20
- Energy proxy: MISO-EE 20/20