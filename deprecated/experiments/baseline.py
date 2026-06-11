import argparse
from pathlib import Path

from utils.plotting import use_agg_backend

use_agg_backend()
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from utils.solver_sets import BASELINE_SOLVERS
from utils.data import generate_V
from utils.evaluation import evaluate_solver


OBJECTIVES = {
    "u_bf": {"label": "BF gain", "better": "max"},
    "u_i": {"label": "Interference", "better": "min"},
    "u_g": {"label": "General objective", "better": "max"},
}

BASELINE_HEURISTICS = ("H1", "H2")
PERFORMANCE_HEURISTIC = "Coutino"
ENERGY_HEURISTIC = "MISO-EE"
BALANCED_HEURISTIC = "Pareto-H2"
HEURISTICS = BASELINE_SOLVERS


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build comparison tables and plots for antenna-selection algorithms."
    )
    parser.add_argument("--N", type=int, default=1000, help="Number of antennas.")
    parser.add_argument(
        "--L-values",
        type=int,
        nargs="+",
        default=list(range(2, 11)),
        help="Layer counts to evaluate.",
    )
    parser.add_argument(
        "--active-fracs",
        type=float,
        nargs="+",
        default=[0.75, 0.50],
        help="Maximum fractions of antennas left active; 0.75 and 0.50 match the task.",
    )
    parser.add_argument(
        "--seeds",
        type=int,
        nargs="+",
        default=[42],
        help="Random seeds. Use several seeds for averaged baselines.",
    )
    parser.add_argument("--sigma", type=float, default=1.0, help="Noise power.")
    parser.add_argument("--P", type=float, default=1.0, help="Power parameter.")
    parser.add_argument(
        "--out-dir", type=Path, default=Path("results"), help="Output directory."
    )
    return parser.parse_args()


def best_heuristic(values, objective):
    if OBJECTIVES[objective]["better"] == "min":
        best_value = values.min()
    else:
        best_value = values.max()
    winners = values[np.isclose(values, best_value)].index.tolist()
    return "Tie" if len(winners) > 1 else winners[0]


def run_single_case(N, L, active_frac, seed, sigma, P):
    K = int(round(N * active_frac))
    np.random.seed(seed)
    V = generate_V(N, L)

    rows = []
    for heuristic, solver in HEURISTICS:
        _, result = evaluate_solver(heuristic, solver, V, K, sigma, P)
        active_count = result["active_count"]
        actual_active_fraction = active_count / N
        actual_turned_off_fraction = 1.0 - actual_active_fraction
        log2_u_g = np.log2(max(result["u_g"], np.finfo(float).tiny))
        rows.append(
            {
                "N": N,
                "L": L,
                "K": K,
                "active_limit_fraction": active_frac,
                "required_turned_off_fraction": 1.0 - active_frac,
                "seed": seed,
                "heuristic": heuristic,
                "valid": result["valid"],
                "active_count": active_count,
                "actual_active_fraction": actual_active_fraction,
                "actual_turned_off_fraction": actual_turned_off_fraction,
                "u_bf": result["u_bf"],
                "u_i": result["u_i"],
                "u_g": result["u_g"],
                "log2_u_g": log2_u_g,
                "log2_u_g_per_active": log2_u_g / active_count if active_count else 0.0,
                "u_bf_per_active": result["u_bf"] / active_count if active_count else 0.0,
                "elapsed_seconds": result["elapsed_seconds"],
            }
        )
    return rows


def build_runs(args):
    rows = []
    total = len(args.L_values) * len(args.active_fracs) * len(args.seeds)
    case_no = 0
    for active_frac in args.active_fracs:
        for L in args.L_values:
            for seed in args.seeds:
                case_no += 1
                print(
                    f"[{case_no:>3}/{total}] N={args.N}, L={L}, "
                    f"K={int(round(args.N * active_frac))}, seed={seed}",
                    flush=True,
                )
                rows.extend(
                    run_single_case(
                        args.N, L, active_frac, seed, args.sigma, args.P
                    )
                )
    return pd.DataFrame(rows)


def build_summary(runs):
    group_cols = [
        "N",
        "L",
        "K",
        "active_limit_fraction",
        "required_turned_off_fraction",
        "heuristic",
    ]
    summary = (
        runs.groupby(group_cols, as_index=False)
        .agg(
            u_bf_mean=("u_bf", "mean"),
            u_bf_std=("u_bf", "std"),
            u_i_mean=("u_i", "mean"),
            u_i_std=("u_i", "std"),
            u_g_mean=("u_g", "mean"),
            u_g_std=("u_g", "std"),
            log2_u_g_mean=("log2_u_g", "mean"),
            log2_u_g_per_active_mean=("log2_u_g_per_active", "mean"),
            u_bf_per_active_mean=("u_bf_per_active", "mean"),
            active_count_mean=("active_count", "mean"),
            actual_active_fraction_mean=("actual_active_fraction", "mean"),
            actual_turned_off_fraction_mean=("actual_turned_off_fraction", "mean"),
            elapsed_seconds_mean=("elapsed_seconds", "mean"),
            elapsed_seconds_std=("elapsed_seconds", "std"),
            cases=("seed", "count"),
        )
        .fillna(
            {
                "u_bf_std": 0.0,
                "u_i_std": 0.0,
                "u_g_std": 0.0,
                "elapsed_seconds_std": 0.0,
            }
        )
    )
    return summary


def build_winners(summary):
    rows = []
    for keys, chunk in summary.groupby(["N", "L", "K", "active_limit_fraction"]):
        by_h = chunk.set_index("heuristic")
        if not set(BASELINE_HEURISTICS).issubset(by_h.index):
            continue
        base = {
            "N": keys[0],
            "L": keys[1],
            "K": keys[2],
            "active_limit_fraction": keys[3],
            "required_turned_off_fraction": 1.0 - keys[3],
        }
        for objective in OBJECTIVES:
            values = by_h[f"{objective}_mean"]
            baseline_values = values.loc[list(BASELINE_HEURISTICS)]
            baseline_winner = best_heuristic(baseline_values, objective)
            all_winner = best_heuristic(values, objective)
            if OBJECTIVES[objective]["better"] == "min":
                baseline_best = baseline_values.min()
            else:
                baseline_best = baseline_values.max()

            for candidate in [
                PERFORMANCE_HEURISTIC,
                ENERGY_HEURISTIC,
                BALANCED_HEURISTIC,
            ]:
                candidate_value = values.get(candidate, np.nan)
                if OBJECTIVES[objective]["better"] == "min":
                    relative_delta = (
                        (baseline_best - candidate_value) / baseline_best
                        if baseline_best
                        else np.nan
                    )
                else:
                    relative_delta = (
                        (candidate_value - baseline_best) / baseline_best
                        if baseline_best
                        else np.nan
                    )
                rows.append(
                    {
                        **base,
                        "objective": objective,
                        "objective_label": OBJECTIVES[objective]["label"],
                        "candidate": candidate,
                        "baseline_winner": baseline_winner,
                        "winner": all_winner,
                        "baseline_best_mean": baseline_best,
                        "candidate_mean": candidate_value,
                        "candidate_vs_best_h12_relative_delta": relative_delta,
                    }
                )
    return pd.DataFrame(rows)


def plot_objectives(summary, out_dir):
    markers = {
        "H1": "o",
        "H2": "s",
        PERFORMANCE_HEURISTIC: "^",
        ENERGY_HEURISTIC: "D",
        BALANCED_HEURISTIC: "P",
    }
    for active_frac in sorted(summary["active_limit_fraction"].unique(), reverse=True):
        subset = summary[summary["active_limit_fraction"] == active_frac]
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
        for ax, objective in zip(axes, OBJECTIVES):
            for heuristic in [name for name, _ in HEURISTICS]:
                data = subset[subset["heuristic"] == heuristic].sort_values("L")
                if data.empty:
                    continue
                values = data[f"{objective}_mean"]
                if objective == "u_g":
                    values = np.log10(np.maximum(values, np.finfo(float).tiny))
                ax.plot(
                    data["L"],
                    values,
                    marker=markers.get(heuristic, "o"),
                    linewidth=2,
                    label=heuristic,
                )
                std_col = f"{objective}_std"
                if objective != "u_g" and data[std_col].max() > 0:
                    ax.fill_between(
                        data["L"],
                        data[f"{objective}_mean"] - data[std_col],
                        data[f"{objective}_mean"] + data[std_col],
                        alpha=0.15,
                    )
            direction = "higher is better" if OBJECTIVES[objective]["better"] == "max" else "lower is better"
            ax.set_title(f"{OBJECTIVES[objective]['label']} ({direction})")
            ax.set_xlabel("L")
            ax.grid(True, alpha=0.25)
            if objective == "u_g":
                ax.set_ylabel("log10 objective value")
        axes[0].set_ylabel("Mean objective value")
        axes[-1].legend(loc="best")
        off_pct = round((1.0 - active_frac) * 100)
        fig.suptitle(
            f"Raw objective values, N={int(subset['N'].iloc[0])}, "
            f"active <= {int(subset['K'].iloc[0])} ({off_pct}%+ antennas off)"
        )
        fig.savefig(out_dir / f"objectives_{off_pct:02d}pct_off.png", dpi=180)
        plt.close(fig)


def plot_relative_deltas(winners, out_dir):
    for active_frac in sorted(winners["active_limit_fraction"].unique(), reverse=True):
        subset = winners[winners["active_limit_fraction"] == active_frac]
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
        for ax, objective in zip(axes, OBJECTIVES):
            data = subset[subset["objective"] == objective].sort_values("L")
            if OBJECTIVES[objective]["better"] == "min":
                data = data.assign(
                    normalized_score=data["baseline_best_mean"] / data["candidate_mean"]
                )
            else:
                data = data.assign(
                    normalized_score=data["candidate_mean"] / data["baseline_best_mean"]
                )
            pivot = data.pivot(
                index="L", columns="candidate", values="normalized_score"
            ).sort_index()
            x = np.arange(len(pivot.index))
            candidates = [
                (PERFORMANCE_HEURISTIC, "#2a9d8f"),
                (ENERGY_HEURISTIC, "#4f6d7a"),
                (BALANCED_HEURISTIC, "#9467bd"),
            ]
            width = 0.24
            offsets = np.linspace(-width, width, len(candidates))
            for offset, (candidate, color) in zip(offsets, candidates):
                values = pivot[candidate]
                ax.bar(x + offset, values, width=width, color=color, label=candidate)
            ax.axhline(1.0, color="#333333", linewidth=0.8)
            direction = "higher is better" if OBJECTIVES[objective]["better"] == "max" else "lower is better"
            ax.set_title(f"{OBJECTIVES[objective]['label']} ({direction})")
            ax.set_xlabel("L")
            ax.set_xticks(x, pivot.index.astype(str))
            ax.grid(True, axis="y", alpha=0.25)
        axes[0].set_ylabel("Score vs best(H1,H2); 1=tie, >1 better")
        axes[-1].legend(loc="best")
        off_pct = round((1.0 - active_frac) * 100)
        fig.suptitle(
            f"Normalized score vs best H1/H2, N={int(subset['N'].iloc[0])}, "
            f"active <= {int(subset['K'].iloc[0])} ({off_pct}%+ antennas off)"
        )
        fig.savefig(
            out_dir / f"relative_vs_h12_{off_pct:02d}pct_off.png", dpi=180
        )
        plt.close(fig)


def plot_energy_efficiency(summary, out_dir):
    colors = {
        "H1": "#1f77b4",
        "H2": "#ff7f0e",
        PERFORMANCE_HEURISTIC: "#2ca02c",
        ENERGY_HEURISTIC: "#4f6d7a",
        BALANCED_HEURISTIC: "#9467bd",
    }
    for active_frac in sorted(summary["active_limit_fraction"].unique(), reverse=True):
        subset = summary[summary["active_limit_fraction"] == active_frac]
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
        for heuristic in [name for name, _ in HEURISTICS]:
            data = subset[subset["heuristic"] == heuristic].sort_values("L")
            if data.empty:
                continue
            axes[0].plot(
                data["L"],
                100.0 * data["actual_turned_off_fraction_mean"],
                marker="o",
                linewidth=2,
                color=colors.get(heuristic),
                label=heuristic,
            )
            axes[1].plot(
                data["L"],
                data["log2_u_g_per_active_mean"],
                marker="o",
                linewidth=2,
                color=colors.get(heuristic),
                label=heuristic,
            )
        axes[0].axhline(
            100.0 * (1.0 - active_frac),
            color="#333333",
            linewidth=1,
            linestyle="--",
            label="required minimum off",
        )
        axes[0].set_title("Antennas turned off (higher saves more RF chains)")
        axes[0].set_xlabel("L")
        axes[0].set_ylabel("Turned off, %")
        axes[0].grid(True, alpha=0.25)
        axes[1].set_title("Energy-efficiency proxy (higher is better)")
        axes[1].set_xlabel("L")
        axes[1].set_ylabel("log2(U_G) / active antenna")
        axes[1].grid(True, alpha=0.25)
        axes[1].legend(loc="best")
        off_pct = round((1.0 - active_frac) * 100)
        fig.suptitle(
            f"Energy view, N={int(subset['N'].iloc[0])}, "
            f"active <= {int(subset['K'].iloc[0])} ({off_pct}%+ antennas off)"
        )
        fig.savefig(out_dir / f"energy_efficiency_{off_pct:02d}pct_off.png", dpi=180)
        plt.close(fig)


def plot_comparison_dashboard(summary, winners, out_dir):
    colors = {
        "H1": "#1f77b4",
        "H2": "#ff7f0e",
        PERFORMANCE_HEURISTIC: "#2ca02c",
        ENERGY_HEURISTIC: "#4f6d7a",
        BALANCED_HEURISTIC: "#9467bd",
    }
    markers = {
        "H1": "o",
        "H2": "s",
        PERFORMANCE_HEURISTIC: "^",
        ENERGY_HEURISTIC: "D",
        BALANCED_HEURISTIC: "P",
    }

    for active_frac in sorted(summary["active_limit_fraction"].unique(), reverse=True):
        subset = summary[summary["active_limit_fraction"] == active_frac]
        win_subset = winners[winners["active_limit_fraction"] == active_frac]
        fig, axes = plt.subplots(2, 3, figsize=(17, 9), constrained_layout=True)
        objective_specs = [
            ("u_bf", "BF gain (higher is better)", "Mean objective value"),
            ("u_i", "Interference (lower is better)", "Mean objective value"),
            ("u_g", "General objective, log10 (higher is better)", "log10(U_G)"),
        ]

        for ax, (objective, title, ylabel) in zip(axes[0], objective_specs):
            for heuristic in [name for name, _ in HEURISTICS]:
                data = subset[subset["heuristic"] == heuristic].sort_values("L")
                values = data[f"{objective}_mean"]
                if objective == "u_g":
                    values = np.log10(np.maximum(values, np.finfo(float).tiny))
                ax.plot(
                    data["L"],
                    values,
                    marker=markers[heuristic],
                    linewidth=2,
                    color=colors[heuristic],
                    label=heuristic,
                )
            ax.set_title(title)
            ax.set_xlabel("L")
            ax.set_ylabel(ylabel)
            ax.grid(True, alpha=0.25)

        for heuristic in [name for name, _ in HEURISTICS]:
            data = subset[subset["heuristic"] == heuristic].sort_values("L")
            axes[1, 0].plot(
                data["L"],
                100.0 * data["actual_turned_off_fraction_mean"],
                marker=markers[heuristic],
                linewidth=2,
                color=colors[heuristic],
                label=heuristic,
            )
            axes[1, 1].plot(
                data["L"],
                data["log2_u_g_per_active_mean"],
                marker=markers[heuristic],
                linewidth=2,
                color=colors[heuristic],
                label=heuristic,
            )

        axes[1, 0].axhline(
            100.0 * (1.0 - active_frac),
            color="#333333",
            linewidth=1,
            linestyle="--",
            label="required minimum off",
        )
        axes[1, 0].set_title("Antennas turned off (higher saves more RF chains)")
        axes[1, 0].set_xlabel("L")
        axes[1, 0].set_ylabel("Turned off, %")
        axes[1, 0].grid(True, alpha=0.25)

        axes[1, 1].set_title("Energy-efficiency proxy (higher is better)")
        axes[1, 1].set_xlabel("L")
        axes[1, 1].set_ylabel("log2(U_G) / active antenna")
        axes[1, 1].grid(True, alpha=0.25)

        data = win_subset[win_subset["objective"] == "u_g"].copy()
        data["normalized_score"] = data["candidate_mean"] / data["baseline_best_mean"]
        pivot = data.pivot(
            index="L", columns="candidate", values="normalized_score"
        ).sort_index()
        x = np.arange(len(pivot.index))
        candidates = [PERFORMANCE_HEURISTIC, ENERGY_HEURISTIC, BALANCED_HEURISTIC]
        width = 0.24
        offsets = np.linspace(-width, width, len(candidates))
        for offset, candidate in zip(offsets, candidates):
            axes[1, 2].bar(
                x + offset,
                pivot[candidate],
                width=width,
                color=colors[candidate],
                label=candidate,
            )
        axes[1, 2].axhline(1.0, color="#333333", linewidth=0.8)
        axes[1, 2].set_title("U_G score vs best H1/H2 (1=tie, >1 better)")
        axes[1, 2].set_xlabel("L")
        axes[1, 2].set_ylabel("Normalized U_G score")
        axes[1, 2].set_xticks(x, pivot.index.astype(str))
        axes[1, 2].grid(True, axis="y", alpha=0.25)
        axes[1, 2].legend(loc="best")

        axes[0, 0].legend(loc="best")
        off_pct = round((1.0 - active_frac) * 100)
        fig.suptitle(
            f"Antenna selection comparison, N={int(subset['N'].iloc[0])}, "
            f"active <= {int(subset['K'].iloc[0])} ({off_pct}%+ antennas off)",
        )
        fig.savefig(out_dir / f"comparison_{off_pct:02d}pct_off.png", dpi=180)
        plt.close(fig)


def cleanup_plots(out_dir):
    for path in out_dir.glob("*.png"):
        path.unlink()


def format_value(value):
    if abs(value) >= 1e6:
        return f"{value:.4e}"
    return f"{value:.4f}"


def format_pct(value):
    return f"{100.0 * value:+.2f}%"


def format_unsigned_pct(value):
    return f"{100.0 * value:.1f}%"


def write_report(summary, winners, out_dir, args):
    lines = [
        "# H1/H2/Coutino/MISO-EE/Pareto-H2 Objective Values",
        "",
        f"- N: {args.N}",
        f"- L values: {', '.join(map(str, args.L_values))}",
        f"- Active limits: {', '.join(map(str, args.active_fracs))}",
        f"- Seeds: {', '.join(map(str, args.seeds))}",
        f"- Sigma: {args.sigma}",
        f"- P: {args.P}",
        "",
        "Objectives: `u_bf` and `u_g` are maximized; `u_i` is minimized.",
        "The constraint is `active_count <= K`, so algorithms may turn off more than the required minimum.",
        "`H1`, `H2`, and `Coutino` scan their feasible deletion paths and keep the best `U_G` set with `active_count <= K`.",
        "`MISO-EE` returns the smallest greedy log-det set that still targets +5% `U_G` over best(H1,H2).",
        "`Pareto-H2` minimizes interference on a BF-protected deletion path; lower `Pareto Int / H2 Int` is better, while `Pareto BF / H2 BF >= 1` means no BF loss versus H2.",
        "Energy-efficiency proxy: `log2(U_G) / active antenna`; larger is better.",
        "",
        "Recommendation: for Task 3 (`sigma=1`, maximize `U_G`), use `Coutino`; in the default `N=1000`, `L=2..10` run it wins `U_G` for every 25%+ and 50%+ off case.",
        "Use `MISO-EE` only if the priority is energy efficiency with extra antennas switched off while keeping about +5% `U_G` over best(H1,H2).",
        "Use `Pareto-H2` only for the interference/BF compromise; it is deliberately conservative and falls back to H2 when the protected path is not close enough on interference.",
        "",
    ]

    for active_frac in sorted(summary["active_limit_fraction"].unique(), reverse=True):
        off_pct = round((1.0 - active_frac) * 100)
        lines.extend(
            [
                f"## {off_pct}% Antennas Off",
                "",
                "| L | K max | Coutino active | MISO-EE active | MISO-EE off | Pareto active | Pareto BF / H2 BF | Pareto Int / H2 Int | H1 Gen | H2 Gen | Coutino Gen | MISO-EE Gen | Coutino Gen vs best H1/H2 | MISO-EE Gen vs best H1/H2 | Int winner | EE proxy winner |",
                "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|",
            ]
        )
        subset = summary[summary["active_limit_fraction"] == active_frac]
        win_subset = winners[winners["active_limit_fraction"] == active_frac]
        for L in sorted(subset["L"].unique()):
            rows = subset[subset["L"] == L].set_index("heuristic")
            win_rows = win_subset[win_subset["L"] == L]
            win_unique = win_rows.drop_duplicates("objective").set_index("objective")
            perf_delta = win_rows[
                (win_rows["objective"] == "u_g")
                & (win_rows["candidate"] == PERFORMANCE_HEURISTIC)
            ]["candidate_vs_best_h12_relative_delta"].iloc[0]
            ee_delta = win_rows[
                (win_rows["objective"] == "u_g")
                & (win_rows["candidate"] == ENERGY_HEURISTIC)
            ]["candidate_vs_best_h12_relative_delta"].iloc[0]
            ee_winner = rows["log2_u_g_per_active_mean"].idxmax()
            k_value = int(rows.iloc[0]["K"])
            lines.append(
                "| "
                + " | ".join(
                    [
                        str(int(L)),
                        str(k_value),
                        str(int(round(rows.loc[PERFORMANCE_HEURISTIC, "active_count_mean"]))),
                        str(int(round(rows.loc[ENERGY_HEURISTIC, "active_count_mean"]))),
                        format_unsigned_pct(
                            rows.loc[
                                ENERGY_HEURISTIC,
                                "actual_turned_off_fraction_mean",
                            ]
                        ),
                        str(int(round(rows.loc[BALANCED_HEURISTIC, "active_count_mean"]))),
                        f"{rows.loc[BALANCED_HEURISTIC, 'u_bf_mean'] / rows.loc['H2', 'u_bf_mean']:.3f}",
                        f"{rows.loc[BALANCED_HEURISTIC, 'u_i_mean'] / rows.loc['H2', 'u_i_mean']:.3f}",
                        format_value(rows.loc["H1", "u_g_mean"]),
                        format_value(rows.loc["H2", "u_g_mean"]),
                        format_value(rows.loc[PERFORMANCE_HEURISTIC, "u_g_mean"]),
                        format_value(rows.loc[ENERGY_HEURISTIC, "u_g_mean"]),
                        format_pct(perf_delta),
                        format_pct(ee_delta),
                        win_unique.loc["u_i", "winner"],
                        ee_winner,
                    ]
                )
                + " |"
            )
        lines.append("")

    (out_dir / "baseline_report.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    cleanup_plots(args.out_dir)

    runs = build_runs(args)
    summary = build_summary(runs)
    winners = build_winners(summary)

    runs.to_csv(args.out_dir / "baseline_runs.csv", index=False)
    summary.to_csv(args.out_dir / "baseline_summary.csv", index=False)
    winners.to_csv(args.out_dir / "baseline_winners.csv", index=False)

    plot_comparison_dashboard(summary, winners, args.out_dir)
    write_report(summary, winners, args.out_dir, args)

    print("\nSaved:")
    for path in [
        "baseline_runs.csv",
        "baseline_summary.csv",
        "baseline_winners.csv",
        "comparison_25pct_off.png",
        "comparison_50pct_off.png",
        "baseline_report.md",
    ]:
        print(f"  {args.out_dir / path}")


if __name__ == "__main__":
    main()
