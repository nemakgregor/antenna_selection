import argparse
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
    solve_h3,
    solve_h3_fast,
    solve_miso_energy_greedy,
    solve_pareto_interference_greedy,
)
from motor_challenge_1205 import generate_V


METRICS = (
    ("u_bf", "BF gain", "max"),
    ("u_i", "Interference", "min"),
    ("u_g", "General objective", "max"),
)

ALGORITHMS = (
    ("H1", lambda V, K, sigma, P: solve_h1(V, K, sigma=sigma, P=P)),
    ("H2", lambda V, K, sigma, P: solve_h2(V, K, sigma=sigma, P=P)),
    (
        "Coutino",
        lambda V, K, sigma, P: solve_coutino_greedy(V, K, sigma=sigma, P=P),
    ),
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
    (
        "S-threshold-BF",
        lambda V, K, sigma, P: solve_h3(
            V, K, target_obj="bf", sigma=sigma, P=P
        ),
    ),
    (
        "S-threshold-Int",
        lambda V, K, sigma, P: solve_h3(
            V, K, target_obj="int", sigma=sigma, P=P
        ),
    ),
    (
        "S-threshold-Gen",
        lambda V, K, sigma, P: solve_h3(
            V, K, target_obj="gen", sigma=sigma, P=P
        ),
    ),
    (
        "N-H3-Fast",
        lambda V, K, sigma, P: solve_h3_fast(V, K),
    ),
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Exact algorithm grid benchmark with compact plots."
    )
    parser.add_argument("--N-values", type=int, nargs="+", default=[1000, 5000, 10000])
    parser.add_argument("--L-values", type=int, nargs="+", default=list(range(1, 11)))
    parser.add_argument("--off-pcts", type=int, nargs="+", default=[25, 50])
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--sigma", type=float, default=1.0)
    parser.add_argument("--P", type=float, default=1.0)
    parser.add_argument("--save-runs", action="store_true")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("results/grid_benchmark"),
    )
    return parser.parse_args()


def run_solver(name, solver, V, K, sigma, P):
    started_at = time.perf_counter()
    with np.errstate(all="ignore"):
        x = solver(V, K, sigma, P)
        elapsed_seconds = time.perf_counter() - started_at
        valid, active_count = check_constraints(x, K)
        u_bf, u_i, u_g = calculate_objectives(V, x, sigma=sigma, P=P)

    if not valid or not np.isfinite([u_bf, u_i, u_g]).all():
        raise RuntimeError(f"Invalid result for {name}.")

    return {
        "heuristic": name,
        "active_count": int(active_count),
        "u_bf": float(u_bf),
        "u_i": float(u_i),
        "u_g": float(u_g),
        "elapsed_seconds": float(elapsed_seconds),
    }


def run_grid(args):
    rows = []
    total_cases = len(args.off_pcts) * len(args.N_values) * len(args.L_values)
    case_no = 0

    for off_pct in args.off_pcts:
        for N in args.N_values:
            if not (0 <= off_pct < 100):
                raise ValueError(f"off_pct={off_pct} must satisfy 0 <= off_pct < 100.")
            K = int(round(N * (1.0 - off_pct / 100.0)))

            for L in args.L_values:
                case_no += 1
                print(
                    f"[{case_no:>3}/{total_cases}] N={N}, L={L}, "
                    f"off={off_pct}%, K_active={K}",
                    flush=True,
                )

                for sample in range(args.samples):
                    seed = args.seed + sample
                    np.random.seed(seed)
                    V = generate_V(N, L)

                    for name, solver in ALGORITHMS:
                        result = run_solver(name, solver, V, K, args.sigma, args.P)
                        rows.append(
                            {
                                "N": N,
                                "L": L,
                                "K": K,
                                "off_pct": off_pct,
                                "sample": sample,
                                "seed": seed,
                                "sigma": args.sigma,
                                "P": args.P,
                                **result,
                            }
                        )

    return pd.DataFrame(rows)


def build_summary(runs):
    group_cols = ["N", "L", "K", "off_pct", "heuristic"]
    return (
        runs.groupby(group_cols, as_index=False)
        .agg(
            active_count_mean=("active_count", "mean"),
            u_bf_mean=("u_bf", "mean"),
            u_bf_std=("u_bf", "std"),
            u_i_mean=("u_i", "mean"),
            u_i_std=("u_i", "std"),
            u_g_mean=("u_g", "mean"),
            u_g_std=("u_g", "std"),
            elapsed_seconds_mean=("elapsed_seconds", "mean"),
            elapsed_seconds_std=("elapsed_seconds", "std"),
            samples=("sample", "count"),
        )
        .fillna(0.0)
    )


def build_wins(runs):
    rows = []
    for keys, chunk in runs.groupby(["N", "L", "K", "off_pct", "sample"]):
        N, L, K, off_pct, sample = keys
        by_h = chunk.set_index("heuristic")
        for metric, _, direction in METRICS:
            values = by_h[metric]
            best_value = values.max() if direction == "max" else values.min()
            winners = values[np.isclose(values, best_value)].index.tolist()
            share = 1.0 / len(winners)
            for winner in winners:
                rows.append(
                    {
                        "N": N,
                        "L": L,
                        "K": K,
                        "off_pct": off_pct,
                        "sample": sample,
                        "metric": metric,
                        "heuristic": winner,
                        "win_share": share,
                    }
                )

    wins = pd.DataFrame(rows)
    return (
        wins.groupby(["N", "L", "K", "off_pct", "metric", "heuristic"], as_index=False)
        .agg(win_share=("win_share", "sum"))
        .assign(win_fraction=lambda df: df["win_share"] / runs["sample"].nunique())
    )


def plot_dashboard(summary, wins, out_dir, N, off_pct):
    heuristics = [name for name, _ in ALGORITHMS]
    colors = dict(zip(heuristics, plt.cm.tab10.colors[: len(heuristics)]))
    data = summary[(summary["N"] == N) & (summary["off_pct"] == off_pct)]
    win_data = wins[(wins["N"] == N) & (wins["off_pct"] == off_pct)]
    L_values = sorted(data["L"].unique())
    K_active = int(data["K"].iloc[0])

    fig, axes = plt.subplots(2, 4, figsize=(22, 10), constrained_layout=True)

    for ax, (metric, title, direction) in zip(axes[0, :3], METRICS):
        mean_col = f"{metric}_mean"
        for heuristic in heuristics:
            curve = data[data["heuristic"] == heuristic].sort_values("L")
            y_values = curve[mean_col]
            if metric == "u_i":
                y_values = y_values.clip(lower=1e-12)
            ax.plot(
                curve["L"],
                y_values,
                marker="o",
                linewidth=2,
                color=colors[heuristic],
                label=heuristic,
            )
        ax.set_title(f"{title} mean ({'higher' if direction == 'max' else 'lower'} is better)")
        ax.set_xlabel("L")
        ax.set_xticks(L_values)
        ax.grid(True, alpha=0.25)
        if metric == "u_i":
            ax.set_yscale("log")
        elif metric == "u_g":
            ax.set_yscale("log")

    ax_time = axes[0, 3]
    for heuristic in heuristics:
        curve = data[data["heuristic"] == heuristic].sort_values("L")
        ax_time.plot(
            curve["L"],
            curve["elapsed_seconds_mean"],
            marker="o",
            linewidth=2,
            color=colors[heuristic],
            label=heuristic,
        )
    ax_time.set_title("Mean runtime per example")
    ax_time.set_xlabel("L")
    ax_time.set_ylabel("seconds")
    ax_time.set_xticks(L_values)
    ax_time.set_yscale("log")
    ax_time.grid(True, alpha=0.25)

    for ax, (metric, title, _) in zip(axes[1, :3], METRICS):
        bottom = np.zeros(len(L_values), dtype=float)
        for heuristic in heuristics:
            shares = []
            for L in L_values:
                row = win_data[
                    (win_data["L"] == L)
                    & (win_data["metric"] == metric)
                    & (win_data["heuristic"] == heuristic)
                ]
                shares.append(float(row["win_fraction"].iloc[0]) if len(row) else 0.0)
            ax.bar(
                L_values,
                shares,
                bottom=bottom,
                color=colors[heuristic],
                label=heuristic,
            )
            bottom += np.asarray(shares)
        ax.set_title(f"{title} win share")
        ax.set_xlabel("L")
        ax.set_ylim(0.0, 1.0)
        ax.set_xticks(L_values)
        ax.grid(True, axis="y", alpha=0.25)

    axes[1, 3].axis("off")
    axes[1, 3].legend(
        handles=[
            plt.Line2D([0], [0], color=colors[name], linewidth=6, label=name)
            for name in heuristics
        ],
        loc="center",
    )

    fig.suptitle(
        f"Grid benchmark, N={N}, {off_pct}% off, K_active={K_active}, "
        f"samples={int(data['samples'].iloc[0])}"
    )
    fig.savefig(out_dir / f"dashboard_{off_pct}pct_off_N{N}.png", dpi=180)
    plt.close(fig)


def write_report(summary, wins, out_dir, args):
    lines = [
        "# Grid Benchmark",
        "",
        f"- N values: {args.N_values}",
        f"- L values: {args.L_values}",
        f"- off percentages: {args.off_pcts}",
        f"- samples per case: {args.samples}",
        f"- algorithms: {[name for name, _ in ALGORITHMS]}",
        "",
        "K_active is computed as round(N * (1 - off_pct / 100)).",
        "All algorithm entries call the corresponding implementation directly.",
        "",
        "## Winners By Mean",
        "",
    ]

    for off_pct in args.off_pcts:
        for N in args.N_values:
            chunk = summary[(summary["off_pct"] == off_pct) & (summary["N"] == N)]
            K_active = int(chunk["K"].iloc[0])
            lines.append(f"### {off_pct}% off, N={N}, K_active={K_active}")
            for metric, title, direction in METRICS:
                mean_col = f"{metric}_mean"
                by_h = chunk.groupby("heuristic", as_index=False)[mean_col].mean()
                row = (
                    by_h.loc[by_h[mean_col].idxmax()]
                    if direction == "max"
                    else by_h.loc[by_h[mean_col].idxmin()]
                )
                lines.append(f"- {title}: {row['heuristic']} ({row[mean_col]:.6g})")
            time_row = (
                chunk.groupby("heuristic", as_index=False)["elapsed_seconds_mean"]
                .mean()
                .sort_values("elapsed_seconds_mean")
                .iloc[0]
            )
            lines.append(
                f"- Fastest avg runtime: {time_row['heuristic']} "
                f"({time_row['elapsed_seconds_mean']:.6g}s)"
            )
            lines.append("")

    (out_dir / "grid_report.md").write_text("\n".join(lines))


def main():
    args = parse_args()
    if args.samples <= 0:
        raise ValueError("--samples must be positive.")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    runs = run_grid(args)
    summary = build_summary(runs)
    wins = build_wins(runs)

    summary.to_csv(args.out_dir / "grid_summary.csv", index=False)
    wins.to_csv(args.out_dir / "grid_wins.csv", index=False)
    if args.save_runs:
        runs.to_csv(args.out_dir / "grid_runs.csv", index=False)

    for off_pct in args.off_pcts:
        for N in args.N_values:
            plot_dashboard(summary, wins, args.out_dir, N, off_pct)

    write_report(summary, wins, args.out_dir, args)
    print(f"wrote results to {args.out_dir}")


if __name__ == "__main__":
    main()
