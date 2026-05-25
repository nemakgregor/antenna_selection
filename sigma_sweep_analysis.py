import argparse
import hashlib
import os
import time
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib-cache").resolve()))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(".cache").resolve()))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from algorithms import (
    calculate_objectives,
    check_constraints,
    solve_coutino_greedy,
    solve_h1,
    solve_h2,
    solve_miso_energy_greedy,
    solve_pareto_interference_greedy,
)
from motor_challenge_1205 import generate_V


HEURISTICS = (
    ("H1", lambda V, K, sigma, P: solve_h1(V, K, sigma=sigma, P=P)),
    ("H2", lambda V, K, sigma, P: solve_h2(V, K, sigma=sigma, P=P)),
    ("Coutino", lambda V, K, sigma, P: solve_coutino_greedy(V, K, sigma=sigma, P=P)),
    (
        "MISO-EE",
        lambda V, K, sigma, P: solve_miso_energy_greedy(
            V, K, sigma=sigma, P=P, target_margin=0.05
        ),
    ),
    (
        "Pareto-H2",
        lambda V, K, sigma, P: solve_pareto_interference_greedy(
            V, K, sigma=sigma, P=P
        ),
    ),
)


def effective_channel_eigvals(V, x, P):
    active_idx = np.flatnonzero(x)
    if len(active_idx) == 0:
        return np.zeros(V.shape[1], dtype=float)
    row_power = np.sum(np.abs(V[active_idx, :]) ** 2, axis=1).real
    max_power = np.max(row_power) if len(row_power) else 0.0
    if max_power <= 0:
        return np.zeros(V.shape[1], dtype=float)
    z2 = P / max_power
    gram = V[active_idx, :].conj().T @ V[active_idx, :]
    matrix = z2 * (gram @ gram.conj().T)
    matrix = 0.5 * (matrix + matrix.conj().T)
    return np.linalg.eigvalsh(matrix).real


def parse_args():
    parser = argparse.ArgumentParser(
        description="Sweep sigma for a fixed antenna-selection problem."
    )
    parser.add_argument("--N", type=int, default=1000)
    parser.add_argument("--L", type=int, default=4)
    parser.add_argument("--active-frac", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--P", type=float, default=1.0)
    parser.add_argument(
        "--sigmas",
        type=float,
        nargs="+",
        default=[0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0],
    )
    parser.add_argument(
        "--out-dir", type=Path, default=Path("results/sigma_sweep")
    )
    return parser.parse_args()


def run_case(V, K, sigma, P):
    rows = []
    for heuristic, solver in HEURISTICS:
        with np.errstate(all="ignore"):
            started_at = time.perf_counter()
            x = solver(V, K, sigma, P)
            elapsed_seconds = time.perf_counter() - started_at
            valid, active_count = check_constraints(x, K)
            u_bf, u_i, u_g = calculate_objectives(V, x, sigma=sigma, P=P)
            eigvals = effective_channel_eigvals(V, x, P)
        if not valid or not np.isfinite([u_bf, u_i, u_g]).all():
            raise RuntimeError(f"Invalid result for {heuristic} at sigma={sigma}")
        rows.append(
            {
                "sigma": sigma,
                "heuristic": heuristic,
                "selection_hash": hashlib.sha1(x.astype(np.int8).tobytes()).hexdigest()[:12],
                "active_count": active_count,
                "turned_off_fraction": 1.0 - active_count / V.shape[0],
                "u_bf": u_bf,
                "u_i": u_i,
                "u_g": u_g,
                "log10_u_g": np.log10(max(u_g, np.finfo(float).tiny)),
                "log2_u_g_per_active": np.log2(max(u_g, np.finfo(float).tiny))
                / active_count,
                "min_channel_eig": float(np.min(eigvals)),
                "max_channel_eig": float(np.max(eigvals)),
                "sigma_over_min_channel_eig": sigma
                / max(float(np.min(eigvals)), np.finfo(float).eps),
                "elapsed_seconds": elapsed_seconds,
            }
        )
    return rows


def add_relative_metrics(runs):
    runs = runs.copy()
    min_sigma = runs["sigma"].min()
    baseline_u_g = (
        runs[runs["sigma"] == min_sigma]
        .set_index("heuristic")["u_g"]
        .to_dict()
    )
    best_h2_bf = runs[runs["heuristic"] == "H2"].set_index("sigma")["u_bf"]
    h2_interference = runs[runs["heuristic"] == "H2"].set_index("sigma")["u_i"]
    best_h12_u_g = (
        runs[runs["heuristic"].isin(["H1", "H2"])]
        .groupby("sigma")["u_g"]
        .max()
    )
    runs["u_g_vs_min_sigma"] = runs.apply(
        lambda row: row["u_g"] / baseline_u_g[row["heuristic"]],
        axis=1,
    )
    runs["u_g_vs_best_h12"] = runs.apply(
        lambda row: row["u_g"] / best_h12_u_g.loc[row["sigma"]],
        axis=1,
    )
    runs["bf_vs_h2"] = runs.apply(
        lambda row: row["u_bf"] / best_h2_bf.loc[row["sigma"]],
        axis=1,
    )
    runs["interference_score_vs_h2"] = runs.apply(
        lambda row: h2_interference.loc[row["sigma"]] / max(row["u_i"], np.finfo(float).eps),
        axis=1,
    )
    return runs


def build_winners(runs):
    rows = []
    for sigma, chunk in runs.groupby("sigma"):
        by_h = chunk.set_index("heuristic")
        bf_winner = by_h["u_bf"].idxmax()
        int_winner = by_h["u_i"].idxmin()
        gen_winner = by_h["u_g"].idxmax()
        ee_winner = by_h["log2_u_g_per_active"].idxmax()
        rows.append(
            {
                "sigma": sigma,
                "bf_winner": bf_winner,
                "interference_winner": int_winner,
                "general_winner": gen_winner,
                "energy_proxy_winner": ee_winner,
                "best_u_g": by_h["u_g"].max(),
                "best_interference": by_h["u_i"].min(),
            }
        )
    return pd.DataFrame(rows).sort_values("sigma")


def plot_sweep(runs, winners, out_dir, args):
    colors = {
        "H1": "#1f77b4",
        "H2": "#ff7f0e",
        "Coutino": "#2ca02c",
        "MISO-EE": "#4f6d7a",
        "Pareto-H2": "#9467bd",
    }
    markers = {
        "H1": "o",
        "H2": "s",
        "Coutino": "^",
        "MISO-EE": "D",
        "Pareto-H2": "P",
    }
    fig, axes = plt.subplots(2, 3, figsize=(17, 9), constrained_layout=True)
    sigma_values = sorted(runs["sigma"].unique())
    sigma_to_x = {sigma: idx for idx, sigma in enumerate(sigma_values)}

    panels = [
        ("u_g_vs_min_sigma", "U_G / U_G at min sigma, higher is better"),
        ("u_g_vs_best_h12", "U_G / best(H1,H2), higher is better"),
        ("bf_vs_h2", "BF / H2 BF, higher is better; sigma-independent"),
        (
            "interference_score_vs_h2",
            "H2 interference / algorithm interference, higher is better; sigma-independent",
        ),
        ("log2_u_g_per_active", "Energy proxy log2(U_G)/active, higher is better"),
        (
            "sigma_over_min_channel_eig",
            "sigma / min eig(V_eq V_eq*), larger means sigma matters more",
        ),
    ]

    for ax, (column, title) in zip(axes.flat, panels):
        for heuristic, _ in HEURISTICS:
            data = runs[runs["heuristic"] == heuristic].sort_values("sigma")
            ax.plot(
                data["sigma"].map(sigma_to_x),
                data[column],
                marker=markers[heuristic],
                linewidth=2,
                color=colors[heuristic],
                label=heuristic,
            )
        ax.set_title(title)
        ax.set_xlabel("sigma")
        ax.set_xticks(range(len(sigma_values)), [f"{sigma:g}" for sigma in sigma_values])
        ax.tick_params(axis="x", labelrotation=30)
        ax.grid(True, alpha=0.25)
    axes[0, 0].legend(loc="best")
    fig.suptitle(
        f"Sigma sweep, N={args.N}, L={args.L}, active <= {int(round(args.N * args.active_frac))} "
        f"({round((1.0 - args.active_frac) * 100)}%+ antennas off), seed={args.seed}"
    )
    fig.savefig(out_dir / "sigma_sweep.png", dpi=180)
    plt.close(fig)


def write_report(runs, winners, out_dir, args):
    lines = [
        "# Sigma Sweep",
        "",
        f"- N: {args.N}",
        f"- L: {args.L}",
        f"- K max: {int(round(args.N * args.active_frac))}",
        f"- Minimum antennas off: {round((1.0 - args.active_frac) * 100)}%",
        f"- Seed: {args.seed}",
        "",
        "Important: `sigma` appears only in `U_G = det(V_eq V_eq* + sigma I)`.",
        "`U_BF` and `U_I` do not contain sigma, so they are expected to be flat unless the selected antenna set changes.",
        "The report below also checks whether the selected antenna set changed across sigma.",
        "",
        "| sigma | best U_G | best BF | best interference | best EE proxy | Coutino U_G / best H1H2 | Pareto BF / H2 BF | Pareto Int / H2 Int | MISO-EE active | MISO-EE off |",
        "|---:|:---:|:---:|:---:|:---:|---:|---:|---:|---:|---:|",
    ]
    for sigma in sorted(runs["sigma"].unique()):
        chunk = runs[runs["sigma"] == sigma].set_index("heuristic")
        win = winners[winners["sigma"] == sigma].iloc[0]
        best_h12_ug = max(chunk.loc["H1", "u_g"], chunk.loc["H2", "u_g"])
        coutino_ratio = chunk.loc["Coutino", "u_g"] / best_h12_ug
        pareto_bf_ratio = chunk.loc["Pareto-H2", "u_bf"] / chunk.loc["H2", "u_bf"]
        pareto_int_ratio = chunk.loc["Pareto-H2", "u_i"] / chunk.loc["H2", "u_i"]
        lines.append(
            "| "
            + " | ".join(
                [
                    f"{sigma:g}",
                    win["general_winner"],
                    win["bf_winner"],
                    win["interference_winner"],
                    win["energy_proxy_winner"],
                    f"{coutino_ratio:.3f}",
                    f"{pareto_bf_ratio:.3f}",
                    f"{pareto_int_ratio:.3f}",
                    str(int(chunk.loc["MISO-EE", "active_count"])),
                    f"{100.0 * chunk.loc['MISO-EE', 'turned_off_fraction']:.1f}%",
                ]
            )
            + " |"
        )
    lines.append("")
    lines.append("## U_G Sigma Sensitivity")
    lines.append("")
    lines.append("| algorithm | U_G at min sigma | U_G at max sigma | relative change | max sigma / min eigenvalue |")
    lines.append("|:---|---:|---:|---:|---:|")
    min_sigma = runs["sigma"].min()
    max_sigma = runs["sigma"].max()
    for heuristic in [name for name, _ in HEURISTICS]:
        low = runs[(runs["heuristic"] == heuristic) & (runs["sigma"] == min_sigma)].iloc[0]
        high = runs[(runs["heuristic"] == heuristic) & (runs["sigma"] == max_sigma)].iloc[0]
        relative_change = high["u_g"] / low["u_g"] - 1.0
        lines.append(
            "| "
            + " | ".join(
                [
                    heuristic,
                    f"{low['u_g']:.4e}",
                    f"{high['u_g']:.4e}",
                    f"{100.0 * relative_change:.2f}%",
                    f"{high['sigma_over_min_channel_eig']:.4f}",
                ]
            )
            + " |"
        )
    lines.append("")
    selection_counts = runs.groupby("heuristic")["selection_hash"].nunique()
    unchanged = selection_counts[selection_counts == 1].index.tolist()
    changed = selection_counts[selection_counts > 1].index.tolist()
    lines.append(
        "Selected-set check: "
        + (
            f"unchanged for {', '.join(unchanged)}"
            if unchanged
            else "no unchanged algorithms"
        )
        + (
            f"; changed for {', '.join(changed)}."
            if changed
            else "; no algorithm changed selection across sigma."
        )
    )
    lines.append("")
    lines.append("Interpretation: H2 is expected to win pure interference; Coutino usually wins raw `U_G`; MISO-EE is the energy-saving compromise; Pareto-H2 targets H2-like interference while protecting BF.")
    (out_dir / "sigma_sweep_report.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    np.random.seed(args.seed)
    V = generate_V(args.N, args.L)
    K = int(round(args.N * args.active_frac))

    rows = []
    for sigma in args.sigmas:
        print(f"sigma={sigma:g}", flush=True)
        rows.extend(run_case(V, K, sigma, args.P))

    runs = pd.DataFrame(rows)
    runs = add_relative_metrics(runs)
    winners = build_winners(runs)
    runs.to_csv(args.out_dir / "sigma_sweep_runs.csv", index=False)
    winners.to_csv(args.out_dir / "sigma_sweep_winners.csv", index=False)
    plot_sweep(runs, winners, args.out_dir, args)
    write_report(runs, winners, args.out_dir, args)

    print("\nSaved:")
    for path in [
        "sigma_sweep_runs.csv",
        "sigma_sweep_winners.csv",
        "sigma_sweep_report.md",
        "sigma_sweep.png",
    ]:
        print(f"  {args.out_dir / path}")


if __name__ == "__main__":
    main()
