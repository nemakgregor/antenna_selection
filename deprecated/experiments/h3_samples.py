import argparse
from pathlib import Path

from utils.plotting import use_agg_backend

use_agg_backend()
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from utils.solver_sets import MOTOR_SOLVERS
from utils.data import generate_V
from utils.evaluation import evaluate_algorithms


METRICS = (
    ("u_bf", "BF gain", "max"),
    ("u_i", "Interference", "min"),
    ("u_g", "General objective", "max"),
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run 100-sample H3 comparison and build plots."
    )
    parser.add_argument("--N", type=int, default=1000)
    parser.add_argument("--L", type=int, default=2)
    parser.add_argument("--K", type=int, default=500)
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--sigma", type=float, default=1.0)
    parser.add_argument("--P", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("results/h3_samples"),
    )
    return parser.parse_args()


def run_samples(args):
    rows = []
    progress_step = max(1, args.samples // 10)

    for sample_idx in range(args.samples):
        np.random.seed(args.seed + sample_idx)
        V = generate_V(args.N, args.L)
        sample_results = evaluate_algorithms(
            V,
            args.K,
            sigma=args.sigma,
            P=args.P,
            random_state=args.seed + sample_idx,
            solvers=MOTOR_SOLVERS,
        )

        for heuristic, result in sample_results.items():
            rows.append(
                {
                    "sample": sample_idx,
                    "seed": args.seed + sample_idx,
                    "N": args.N,
                    "L": args.L,
                    "K": args.K,
                    "sigma": args.sigma,
                    "P": args.P,
                    "heuristic": heuristic,
                    **result,
                }
            )

        if (sample_idx + 1) % progress_step == 0 or sample_idx + 1 == args.samples:
            print(f"completed {sample_idx + 1}/{args.samples}", flush=True)

    return pd.DataFrame(rows)


def build_summary(runs):
    summary = (
        runs.groupby("heuristic", as_index=False)
        .agg(
            valid_rate=("valid", "mean"),
            active_count_mean=("active_count", "mean"),
            active_count_std=("active_count", "std"),
            u_bf_mean=("u_bf", "mean"),
            u_bf_std=("u_bf", "std"),
            u_i_mean=("u_i", "mean"),
            u_i_std=("u_i", "std"),
            u_g_mean=("u_g", "mean"),
            u_g_std=("u_g", "std"),
        )
        .fillna(0.0)
    )

    win_rows = []
    for metric, _, direction in METRICS:
        for sample, chunk in runs.groupby("sample"):
            values = chunk.set_index("heuristic")[metric]
            winner = values.idxmax() if direction == "max" else values.idxmin()
            win_rows.append({"sample": sample, "metric": metric, "winner": winner})

    wins = pd.DataFrame(win_rows)
    win_counts = (
        wins.groupby(["metric", "winner"])
        .size()
        .rename("wins")
        .reset_index()
        .rename(columns={"winner": "heuristic"})
    )
    summary = summary.merge(
        win_counts.pivot(index="heuristic", columns="metric", values="wins")
        .fillna(0)
        .reset_index()
        .rename(
            columns={
                "u_bf": "bf_wins",
                "u_i": "interference_wins",
                "u_g": "general_wins",
            }
        ),
        on="heuristic",
        how="left",
    ).fillna({"bf_wins": 0, "interference_wins": 0, "general_wins": 0})
    return summary


def plot_metric_distributions(runs, out_dir):
    heuristics = runs["heuristic"].drop_duplicates().tolist()
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), constrained_layout=True)

    for ax, (metric, title, direction) in zip(axes, METRICS):
        data = [runs[runs["heuristic"] == heuristic][metric].to_numpy() for heuristic in heuristics]
        ax.boxplot(data, tick_labels=heuristics, showfliers=False)
        ax.set_title(f"{title} distribution ({'higher' if direction == 'max' else 'lower'} is better)")
        ax.tick_params(axis="x", labelrotation=35)
        ax.grid(True, axis="y", alpha=0.25)
        if metric in {"u_i", "u_g"}:
            ax.set_yscale("log")

    fig.savefig(out_dir / "objective_distributions.png", dpi=180)
    plt.close(fig)


def plot_metric_means(summary, out_dir):
    heuristics = summary["heuristic"].tolist()
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), constrained_layout=True)

    for ax, (metric, title, direction) in zip(axes, METRICS):
        mean_col = f"{metric}_mean"
        std_col = f"{metric}_std"
        ax.bar(heuristics, summary[mean_col], yerr=summary[std_col], capsize=3)
        ax.set_title(f"{title} mean +/- std ({'higher' if direction == 'max' else 'lower'} is better)")
        ax.tick_params(axis="x", labelrotation=35)
        ax.grid(True, axis="y", alpha=0.25)
        if metric in {"u_i", "u_g"}:
            ax.set_yscale("log")

    fig.savefig(out_dir / "objective_means.png", dpi=180)
    plt.close(fig)


def plot_win_counts(summary, out_dir):
    heuristics = summary["heuristic"].tolist()
    x = np.arange(len(heuristics))
    width = 0.27

    fig, ax = plt.subplots(figsize=(14, 5.5), constrained_layout=True)
    ax.bar(x - width, summary["bf_wins"], width=width, label="BF wins")
    ax.bar(x, summary["interference_wins"], width=width, label="Interference wins")
    ax.bar(x + width, summary["general_wins"], width=width, label="General wins")
    ax.set_xticks(x, heuristics, rotation=35)
    ax.set_ylabel("Wins over samples")
    ax.set_title("Per-sample winners by objective")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend()

    fig.savefig(out_dir / "win_counts.png", dpi=180)
    plt.close(fig)


def write_report(summary, out_dir, args):
    lines = [
        "# H3 Sample Comparison",
        "",
        f"- N: {args.N}",
        f"- L: {args.L}",
        f"- K: {args.K}",
        f"- samples: {args.samples}",
        f"- sigma: {args.sigma}",
        f"- P: {args.P}",
        f"- seed range: {args.seed}..{args.seed + args.samples - 1}",
        "",
        "## Winners By Mean",
        "",
    ]

    for metric, title, direction in METRICS:
        mean_col = f"{metric}_mean"
        row = (
            summary.loc[summary[mean_col].idxmax()]
            if direction == "max"
            else summary.loc[summary[mean_col].idxmin()]
        )
        lines.append(f"- {title}: {row['heuristic']} ({row[mean_col]:.6g})")

    lines.extend(["", "## Summary", ""])
    columns = summary.columns.tolist()
    lines.append("| " + " | ".join(columns) + " |")
    lines.append("| " + " | ".join(["---"] * len(columns)) + " |")
    for _, row in summary.iterrows():
        values = []
        for column in columns:
            value = row[column]
            if isinstance(value, float):
                values.append(f"{value:.6g}")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    lines.append("")
    (out_dir / "h3_samples_report.md").write_text("\n".join(lines))


def main():
    args = parse_args()
    if not (0 <= args.K <= args.N):
        raise ValueError("K must satisfy 0 <= K <= N.")
    if args.samples <= 0:
        raise ValueError("samples must be positive.")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    runs = run_samples(args)
    summary = build_summary(runs)

    runs.to_csv(args.out_dir / "h3_samples_runs.csv", index=False)
    summary.to_csv(args.out_dir / "h3_samples_summary.csv", index=False)
    plot_metric_distributions(runs, args.out_dir)
    plot_metric_means(summary, args.out_dir)
    plot_win_counts(summary, args.out_dir)
    write_report(summary, args.out_dir, args)

    print(f"wrote results to {args.out_dir}")


if __name__ == "__main__":
    main()
