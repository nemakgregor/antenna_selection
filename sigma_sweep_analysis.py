import argparse
import hashlib
import json
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
    solve_cap_window_gen,
    solve_frame_portfolio,
    solve_h1,
    solve_h2,
    solve_h3,
    solve_h3_fast,
    solve_miso_energy_greedy,
    solve_pareto_interference_greedy,
)
from motor_challenge_1205 import generate_V


FRAME_FAST_KWARGS = {
    "max_refined_starts": 3,
    "max_passes": 2,
    "remove_limit": 60,
    "add_limit": 60,
    "lambdas": (),
}

HEURISTICS = (
    ("H1", lambda V, K, sigma, P, random_state: solve_h1(V, K, sigma=sigma, P=P)),
    ("H2", lambda V, K, sigma, P, random_state: solve_h2(V, K, sigma=sigma, P=P)),
    (
        "Coutino",
        lambda V, K, sigma, P, random_state: solve_coutino_greedy(
            V, K, sigma=sigma, P=P
        ),
    ),
    (
        "MISO-EE",
        lambda V, K, sigma, P, random_state: solve_miso_energy_greedy(
            V, K, sigma=sigma, P=P, target_margin=0.05
        ),
    ),
    (
        "Pareto-H2",
        lambda V, K, sigma, P, random_state: solve_pareto_interference_greedy(
            V, K, sigma=sigma, P=P
        ),
    ),
    (
        "H3-threshold-BF",
        lambda V, K, sigma, P, random_state: solve_h3(
            V, K, target_obj="bf", sigma=sigma, P=P
        ),
    ),
    (
        "H3-threshold-Int",
        lambda V, K, sigma, P, random_state: solve_h3(
            V, K, target_obj="int", sigma=sigma, P=P
        ),
    ),
    (
        "H3-threshold-Gen",
        lambda V, K, sigma, P, random_state: solve_h3(
            V, K, target_obj="gen", sigma=sigma, P=P
        ),
    ),
    (
        "Frame-BF",
        lambda V, K, sigma, P, random_state: solve_frame_portfolio(
            V,
            K,
            target_obj="bf",
            sigma=sigma,
            P=P,
            random_state=random_state,
            **FRAME_FAST_KWARGS,
        ),
    ),
    (
        "Frame-Int",
        lambda V, K, sigma, P, random_state: solve_frame_portfolio(
            V,
            K,
            target_obj="int",
            sigma=sigma,
            P=P,
            random_state=random_state,
            **FRAME_FAST_KWARGS,
        ),
    ),
    (
        "Frame-Gen",
        lambda V, K, sigma, P, random_state: solve_frame_portfolio(
            V,
            K,
            target_obj="gen",
            sigma=sigma,
            P=P,
            random_state=random_state,
            **FRAME_FAST_KWARGS,
        ),
    ),
    (
        "CapWindow-Gen",
        lambda V, K, sigma, P, random_state: solve_cap_window_gen(
            V,
            K,
            sigma=sigma,
            P=P,
            random_state=random_state,
        ),
    ),
    (
        "H3-Fast",
        lambda V, K, sigma, P, random_state: solve_h3_fast(
            V, K, random_state=random_state
        ),
    ),
)

METRICS = (
    ("u_bf", "BF gain", "max"),
    ("u_i", "Interference", "min"),
    ("u_g", "General objective", "max"),
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
    parser.add_argument(
        "--K-values",
        type=int,
        nargs="+",
        default=None,
        help="Explicit active antenna limits. Overrides --active-frac.",
    )
    parser.add_argument(
        "--K-pcts",
        type=float,
        nargs="+",
        default=None,
        help="Active antenna percentages, e.g. 25 means K=round(0.25*N).",
    )
    parser.add_argument(
        "--off-pcts",
        type=float,
        nargs="+",
        default=None,
        help="Disabled antenna percentages, matching grid_benchmark off_pct.",
    )
    parser.add_argument("--samples", type=int, default=1)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--P", type=float, default=1.0)
    parser.add_argument(
        "--sigmas",
        type=float,
        nargs="+",
        default=[
            0.001,
            0.003,
            0.01,
            0.03,
            0.1,
            0.3,
            1.0,
            3.0,
            10.0,
            30.0,
            100.0,
            300.0,
            1000.0,
            3000.0,
            10000.0,
            30000.0,
            100000.0,
        ],
    )
    parser.add_argument(
        "--out-dir", type=Path, default=Path("results/sigma_sweep")
    )
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=1,
        help="Rewrite sigma_sweep_runs.csv and progress JSON every N completed cases.",
    )
    parser.add_argument(
        "--summary-every",
        type=int,
        default=10,
        help="Refresh summary/report CSV files every N completed cases. Use 0 for final only.",
    )
    parser.add_argument(
        "--plot-every",
        type=int,
        default=0,
        help="Refresh plots every N completed cases. Use 0 for final only.",
    )
    parser.add_argument(
        "--refresh-from-runs",
        type=Path,
        default=None,
        help="Read an existing sigma_sweep_runs.csv and rebuild reports/plots only.",
    )
    return parser.parse_args()


def parse_k_cases(args):
    explicit_modes = [
        args.K_values is not None,
        args.K_pcts is not None,
        args.off_pcts is not None,
    ]
    if sum(explicit_modes) > 1:
        raise ValueError("Use only one of --K-values, --K-pcts, or --off-pcts.")

    cases = []
    if args.K_values is not None:
        for K in args.K_values:
            active_pct = 100.0 * K / args.N
            cases.append(
                {
                    "K": int(K),
                    "K_mode": "absolute",
                    "K_value": int(K),
                    "active_pct": active_pct,
                    "off_pct": 100.0 - active_pct,
                    "K_label": f"K<={int(K)} ({active_pct:.3g}% active, absolute K)",
                }
            )
    elif args.K_pcts is not None:
        for pct in args.K_pcts:
            K = int(round(args.N * pct / 100.0))
            cases.append(
                {
                    "K": K,
                    "K_mode": "active_pct",
                    "K_value": float(pct),
                    "active_pct": float(pct),
                    "off_pct": 100.0 - float(pct),
                    "K_label": f"K<={K} ({pct:g}% active)",
                }
            )
    elif args.off_pcts is not None:
        for pct in args.off_pcts:
            K = int(round(args.N * (1.0 - pct / 100.0)))
            cases.append(
                {
                    "K": K,
                    "K_mode": "off_pct",
                    "K_value": float(pct),
                    "active_pct": 100.0 - float(pct),
                    "off_pct": float(pct),
                    "K_label": f"K<={K} ({pct:g}% off, {100.0 - pct:g}% active)",
                }
            )
    else:
        K = int(round(args.N * args.active_frac))
        active_pct = 100.0 * args.active_frac
        cases.append(
            {
                "K": K,
                "K_mode": "active_frac",
                "K_value": float(args.active_frac),
                "active_pct": active_pct,
                "off_pct": 100.0 - active_pct,
                "K_label": f"K<={K} ({active_pct:.3g}% active)",
            }
        )

    seen = {}
    for case in cases:
        K = case["K"]
        if K in seen:
            raise ValueError(
                f"Duplicate K={K} from {seen[K]!r} and {case['K_label']!r}."
            )
        seen[K] = case["K_label"]
    return cases


def geometric_mean(values):
    values = np.asarray(values, dtype=float)
    values = np.maximum(values, np.finfo(float).eps)
    return float(np.exp(np.mean(np.log(values))))


def run_case(V, K, sigma, P, random_state):
    rows = []
    for heuristic, solver in HEURISTICS:
        with np.errstate(all="ignore"):
            started_at = time.perf_counter()
            x = solver(V, K, sigma, P, random_state)
            elapsed_seconds = time.perf_counter() - started_at
            valid, active_count = check_constraints(x, K)
            u_bf, u_i, u_g = calculate_objectives(V, x, sigma=sigma, P=P)
            eigvals = effective_channel_eigvals(V, x, P)
        if not valid or not np.isfinite([u_bf, u_i, u_g]).all():
            raise RuntimeError(f"Invalid result for {heuristic} at sigma={sigma}")
        min_eig = float(np.min(eigvals))
        max_eig = float(np.max(eigvals))
        geom_eig = geometric_mean(eigvals)
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
                "min_channel_eig": min_eig,
                "geom_channel_eig": geom_eig,
                "max_channel_eig": max_eig,
                "condition_channel_eig": max_eig / max(min_eig, np.finfo(float).eps),
                "sigma_over_min_channel_eig": sigma
                / max(min_eig, np.finfo(float).eps),
                "sigma_over_geom_channel_eig": sigma
                / max(geom_eig, np.finfo(float).eps),
                "sigma_over_max_channel_eig": sigma
                / max(max_eig, np.finfo(float).eps),
                "elapsed_seconds": elapsed_seconds,
            }
        )
    return rows


def add_relative_metrics(runs):
    runs = runs.copy()
    case_cols = ["K", "sample", "sigma"]
    runs["best_u_g_case"] = runs.groupby(case_cols)["u_g"].transform("max")
    runs["best_u_bf_case"] = runs.groupby(case_cols)["u_bf"].transform("max")
    runs["best_u_i_case"] = runs.groupby(case_cols)["u_i"].transform("min")
    runs["u_g_vs_best_case"] = runs["u_g"] / runs["best_u_g_case"]
    runs["u_bf_vs_best_case"] = runs["u_bf"] / runs["best_u_bf_case"]
    runs["u_i_vs_best_case"] = runs["best_u_i_case"] / np.maximum(
        runs["u_i"], np.finfo(float).eps
    )

    best_h12_u_g = (
        runs[runs["heuristic"].isin(["H1", "H2"])]
        .groupby(case_cols)["u_g"]
        .max()
        .rename("best_h12_u_g")
        .reset_index()
    )
    runs = runs.merge(best_h12_u_g, on=case_cols, how="left")
    runs["u_g_vs_best_h12"] = runs["u_g"] / runs["best_h12_u_g"]

    h2_metrics = (
        runs[runs["heuristic"] == "H2"][case_cols + ["u_bf", "u_i"]]
        .rename(columns={"u_bf": "h2_u_bf", "u_i": "h2_u_i"})
    )
    runs = runs.merge(h2_metrics, on=case_cols, how="left")
    runs["bf_vs_h2"] = runs["u_bf"] / runs["h2_u_bf"]
    runs["interference_score_vs_h2"] = runs["h2_u_i"] / np.maximum(
        runs["u_i"], np.finfo(float).eps
    )
    return runs


def build_winners(runs):
    rows = []
    for keys, chunk in runs.groupby(["K", "sample", "sigma"]):
        K, sample, sigma = keys
        by_h = chunk.set_index("heuristic")
        for metric, label, direction in METRICS:
            values = by_h[metric]
            best_value = values.min() if direction == "min" else values.max()
            winners = values[np.isclose(values, best_value)].index.tolist()
            share = 1.0 / len(winners)
            for winner in winners:
                rows.append(
                    {
                        "K": K,
                        "sample": sample,
                        "sigma": sigma,
                        "metric": metric,
                        "metric_label": label,
                        "heuristic": winner,
                        "win_share": share,
                        "winner_hit": 1.0,
                    }
                )

    wins = pd.DataFrame(rows)
    return (
        wins.groupby(["K", "sigma", "metric", "metric_label", "heuristic"], as_index=False)
        .agg(win_share=("win_share", "sum"), winner_hits=("winner_hit", "sum"))
        .merge(
            runs.groupby(["K", "sigma"], as_index=False)["sample"]
            .nunique()
            .rename(columns={"sample": "samples"}),
            on=["K", "sigma"],
            how="left",
        )
        .assign(
            win_fraction=lambda df: df["win_share"] / df["samples"],
            winner_rate=lambda df: df["winner_hits"] / df["samples"],
        )
        .sort_values(["K", "metric", "sigma", "heuristic"])
    )


def build_summary(runs):
    group_cols = ["K", "sigma", "heuristic"]
    return (
        runs.groupby(group_cols, as_index=False)
        .agg(
            K_label=("K_label", "first"),
            K_mode=("K_mode", "first"),
            K_value=("K_value", "first"),
            active_pct=("active_pct", "first"),
            off_pct=("off_pct", "first"),
            active_count_mean=("active_count", "mean"),
            active_count_std=("active_count", "std"),
            u_bf_mean=("u_bf", "mean"),
            u_bf_std=("u_bf", "std"),
            u_i_mean=("u_i", "mean"),
            u_i_std=("u_i", "std"),
            u_g_mean=("u_g", "mean"),
            u_g_std=("u_g", "std"),
            log10_u_g_mean=("log10_u_g", "mean"),
            log10_u_g_std=("log10_u_g", "std"),
            u_g_vs_best_case_mean=("u_g_vs_best_case", "mean"),
            u_bf_vs_best_case_mean=("u_bf_vs_best_case", "mean"),
            u_i_vs_best_case_mean=("u_i_vs_best_case", "mean"),
            u_g_vs_best_h12_mean=("u_g_vs_best_h12", "mean"),
            log2_u_g_per_active_mean=("log2_u_g_per_active", "mean"),
            min_channel_eig_mean=("min_channel_eig", "mean"),
            geom_channel_eig_mean=("geom_channel_eig", "mean"),
            max_channel_eig_mean=("max_channel_eig", "mean"),
            condition_channel_eig_mean=("condition_channel_eig", "mean"),
            sigma_over_geom_channel_eig_mean=("sigma_over_geom_channel_eig", "mean"),
            elapsed_seconds_mean=("elapsed_seconds", "mean"),
            elapsed_seconds_std=("elapsed_seconds", "std"),
            samples=("sample", "count"),
            selection_hashes=("selection_hash", "nunique"),
        )
        .fillna(0.0)
    )


def build_selection_stability(runs):
    by_sample = (
        runs.groupby(["K", "heuristic", "sample"], as_index=False)
        .agg(
            K_label=("K_label", "first"),
            unique_selected_sets=("selection_hash", "nunique"),
        )
    )
    return (
        by_sample.groupby(["K", "heuristic"], as_index=False)
        .agg(
            K_label=("K_label", "first"),
            unique_selected_sets_mean=("unique_selected_sets", "mean"),
            unique_selected_sets_min=("unique_selected_sets", "min"),
            unique_selected_sets_max=("unique_selected_sets", "max"),
        )
        .sort_values(["K", "heuristic"])
    )


def build_mean_leaders(summary):
    rows = []
    for keys, chunk in summary.groupby(["K", "sigma"]):
        K, sigma = keys
        by_h = chunk.set_index("heuristic")
        for metric, label, direction in METRICS:
            mean_col = f"{metric}_mean"
            row = (
                by_h.loc[by_h[mean_col].idxmin()]
                if direction == "min"
                else by_h.loc[by_h[mean_col].idxmax()]
            )
            rows.append(
                {
                    "K": K,
                    "sigma": sigma,
                    "metric": metric,
                    "metric_label": label,
                    "direction": direction,
                    "leader": row.name,
                    "leader_value": row[mean_col],
                    "leader_geom_eig": row["geom_channel_eig_mean"],
                    "leader_sigma_over_geom_eig": row[
                        "sigma_over_geom_channel_eig_mean"
                    ],
                }
            )
    return pd.DataFrame(rows).sort_values(["K", "metric", "sigma"])


def atomic_write_csv(frame, path):
    tmp_path = path.with_name(path.name + ".tmp")
    frame.to_csv(tmp_path, index=False)
    os.replace(tmp_path, path)


def atomic_write_json(payload, path):
    tmp_path = path.with_name(path.name + ".tmp")
    tmp_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(tmp_path, path)


def write_runs_checkpoint(rows, out_dir, progress):
    runs = pd.DataFrame(rows)
    runs = add_relative_metrics(runs)
    atomic_write_csv(runs, out_dir / "sigma_sweep_runs.csv")
    atomic_write_json(progress, out_dir / "sigma_sweep_progress.json")
    return runs


def write_derived_outputs(runs, out_dir, args, include_plots):
    winners = build_winners(runs)
    summary = build_summary(runs)
    leaders = build_mean_leaders(summary)
    stability = build_selection_stability(runs)

    atomic_write_csv(summary, out_dir / "sigma_sweep_summary.csv")
    atomic_write_csv(winners, out_dir / "sigma_sweep_winners.csv")
    atomic_write_csv(leaders, out_dir / "sigma_sweep_mean_leaders.csv")
    atomic_write_csv(stability, out_dir / "sigma_sweep_selection_stability.csv")
    write_report(summary, winners, leaders, stability, out_dir, args)
    if include_plots:
        plot_sweep(summary, winners, leaders, out_dir, args)

    return summary, winners, leaders, stability


def plot_sweep(summary, winners, leaders, out_dir, args):
    colors = {
        "H1": "#1f77b4",
        "H2": "#ff7f0e",
        "Coutino": "#2ca02c",
        "MISO-EE": "#4f6d7a",
        "Pareto-H2": "#9467bd",
        "H3-threshold-BF": "#8c564b",
        "H3-threshold-Int": "#e377c2",
        "H3-threshold-Gen": "#7f7f7f",
        "Frame-BF": "#17a398",
        "Frame-Int": "#b56576",
        "Frame-Gen": "#3366cc",
        "H3-Fast": "#bcbd22",
    }
    markers = {
        "H1": "o",
        "H2": "s",
        "Coutino": "^",
        "MISO-EE": "D",
        "Pareto-H2": "P",
        "H3-threshold-BF": "v",
        "H3-threshold-Int": "X",
        "H3-threshold-Gen": "*",
        "Frame-BF": ">",
        "Frame-Int": "<",
        "Frame-Gen": "8",
        "H3-Fast": "h",
    }
    for K in sorted(summary["K"].unique()):
        data_k = summary[summary["K"] == K]
        k_label = data_k["K_label"].iloc[0]
        wins_k = winners[(winners["K"] == K) & (winners["metric"] == "u_g")]
        sigma_values = sorted(data_k["sigma"].unique())
        sigma_to_x = {sigma: idx for idx, sigma in enumerate(sigma_values)}
        fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), constrained_layout=True)

        panels = [
            ("u_bf_mean", "BF gain, higher is better"),
            ("u_i_mean", "Interference, lower is better"),
            ("u_g_mean", "U_G = det(V_eq V_eq* + sigma I), higher is better"),
        ]

        for ax, (column, title) in zip(axes.flat, panels):
            for heuristic, _ in HEURISTICS:
                curve = data_k[data_k["heuristic"] == heuristic].sort_values("sigma")
                ax.plot(
                    curve["sigma"].map(sigma_to_x),
                    curve[column],
                    marker=markers[heuristic],
                    linewidth=1.8,
                    markersize=4,
                    color=colors[heuristic],
                    label=heuristic,
                )
            ax.set_title(title)
            ax.set_xlabel("sigma")
            ax.set_xticks(
                range(len(sigma_values)), [f"{sigma:g}" for sigma in sigma_values]
            )
            ax.tick_params(axis="x", labelrotation=35)
            ax.grid(True, alpha=0.25)
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))

        axes[0].legend(loc="best", fontsize=8)
        fig.suptitle(
            f"Raw objective values, N={args.N}, L={args.L}, {k_label}, samples={args.samples}"
        )
        fig.savefig(out_dir / f"sigma_sweep_K{K}.png", dpi=180)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(14, 4.8), constrained_layout=True)
        bottom = np.zeros(len(sigma_values), dtype=float)
        for heuristic, _ in HEURISTICS:
            shares = []
            for sigma in sigma_values:
                row = wins_k[
                    (wins_k["sigma"] == sigma) & (wins_k["heuristic"] == heuristic)
                ]
                shares.append(float(row["win_fraction"].iloc[0]) if len(row) else 0.0)
            ax.bar(
                range(len(sigma_values)),
                shares,
                bottom=bottom,
                color=colors[heuristic],
                label=heuristic,
            )
            bottom += np.asarray(shares)
        ax.set_title(f"Per-sample U_G winner share, N={args.N}, L={args.L}, {k_label}")
        ax.set_xlabel("sigma")
        ax.set_ylabel("split win share")
        ax.set_xticks(range(len(sigma_values)), [f"{sigma:g}" for sigma in sigma_values])
        ax.tick_params(axis="x", labelrotation=35)
        ax.set_ylim(0.0, 1.0)
        ax.grid(True, axis="y", alpha=0.25)
        ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8)
        fig.savefig(out_dir / f"sigma_winners_K{K}.png", dpi=180)
        plt.close(fig)


def leader_segments(leader_rows):
    segments = []
    current = None
    start_sigma = None
    end_sigma = None
    for _, row in leader_rows.sort_values("sigma").iterrows():
        leader = row["leader"]
        sigma = row["sigma"]
        if current is None:
            current = leader
            start_sigma = sigma
            end_sigma = sigma
            continue
        if leader == current:
            end_sigma = sigma
            continue
        segments.append((start_sigma, end_sigma, current))
        current = leader
        start_sigma = sigma
        end_sigma = sigma
    if current is not None:
        segments.append((start_sigma, end_sigma, current))
    return segments


def format_sigma(value):
    return f"{value:g}"


def write_report(summary, winners, leaders, stability, out_dir, args):
    k_cases = (
        summary[["K", "K_label"]]
        .drop_duplicates()
        .sort_values("K")["K_label"]
        .tolist()
    )
    lines = [
        "# Sigma Sweep",
        "",
        f"- N: {args.N}",
        f"- L: {args.L}",
        f"- K cases: {', '.join(k_cases)}",
        f"- Samples: {args.samples}",
        f"- Seed range: {args.seed}..{args.seed + args.samples - 1}",
        f"- Sigma values: {', '.join(format_sigma(s) for s in sorted(summary['sigma'].unique()))}",
        f"- Algorithms: {', '.join(name for name, _ in HEURISTICS)}",
        "",
        "Plots show raw mean objective values only: `U_BF`, `U_I`, and `U_G`.",
        "`U_G = det(V_eq V_eq* + sigma I)` is the general objective.",
        "If `lambda_i` are eigenvalues of `V_eq V_eq*`, then `U_G = product_i(lambda_i + sigma)`.",
        "For `L=4`, `U_G = sigma^4 + sigma^3 sum(lambda_i) + sigma^2 e2(lambda) + sigma e3(lambda) + product(lambda_i)`.",
        "The first algorithm-dependent term at large sigma is `sigma^3 sum(lambda_i)`, i.e. it tracks `U_BF = trace(V_eq V_eq*)`.",
        "`U_BF` and `U_I` formulas do not contain `sigma`; if their curves move across sigma, the selected antenna set changed.",
        "",
    ]

    for K in sorted(summary["K"].unique()):
        k_label = summary[summary["K"] == K]["K_label"].iloc[0]
        lines.extend([f"## {k_label}", ""])
        gen_leaders = leaders[(leaders["K"] == K) & (leaders["metric"] == "u_g")]
        segments = leader_segments(gen_leaders)
        segment_text = ", ".join(
            (
                f"{format_sigma(start)}: {leader}"
                if start == end
                else f"{format_sigma(start)}..{format_sigma(end)}: {leader}"
            )
            for start, end, leader in segments
        )
        lines.append(f"Mean `U_G` leader segments: {segment_text}.")
        lines.append("")
        lines.append("| sigma | U_G mean leader | U_G mean | BF mean leader | BF mean | Interference mean leader | Interference mean |")
        lines.append("|---:|:---|---:|:---|---:|:---|---:|")

        for sigma in sorted(summary[summary["K"] == K]["sigma"].unique()):
            leader_rows = leaders[(leaders["K"] == K) & (leaders["sigma"] == sigma)]
            gen = leader_rows[leader_rows["metric"] == "u_g"].iloc[0]
            bf = leader_rows[leader_rows["metric"] == "u_bf"].iloc[0]
            ui = leader_rows[leader_rows["metric"] == "u_i"].iloc[0]
            lines.append(
                "| "
                + " | ".join(
                    [
                        format_sigma(sigma),
                        gen["leader"],
                        f"{gen['leader_value']:.4e}",
                        bf["leader"],
                        f"{bf['leader_value']:.4e}",
                        ui["leader"],
                        f"{ui['leader_value']:.4e}",
                    ]
                )
                + " |"
            )

        lines.extend(["", "### Raw Endpoint Comparison", ""])
        min_sigma = summary[summary["K"] == K]["sigma"].min()
        max_sigma = summary[summary["K"] == K]["sigma"].max()
        lines.append("| algorithm | BF at min sigma | BF at max sigma | Interference at min sigma | Interference at max sigma | U_G at min sigma | U_G at max sigma | unique selected sets/sample |")
        lines.append("|:---|---:|---:|---:|---:|---:|---:|---:|")
        for heuristic, _ in HEURISTICS:
            low = summary[
                (summary["K"] == K)
                & (summary["sigma"] == min_sigma)
                & (summary["heuristic"] == heuristic)
            ].iloc[0]
            high = summary[
                (summary["K"] == K)
                & (summary["sigma"] == max_sigma)
                & (summary["heuristic"] == heuristic)
            ].iloc[0]
            stability_row = stability[
                (stability["K"] == K) & (stability["heuristic"] == heuristic)
            ].iloc[0]
            set_range = (
                f"{stability_row['unique_selected_sets_mean']:.2f} "
                f"({int(stability_row['unique_selected_sets_min'])}.."
                f"{int(stability_row['unique_selected_sets_max'])})"
            )
            lines.append(
                "| "
                + " | ".join(
                    [
                        heuristic,
                        f"{low['u_bf_mean']:.4e}",
                        f"{high['u_bf_mean']:.4e}",
                        f"{low['u_i_mean']:.4e}",
                        f"{high['u_i_mean']:.4e}",
                        f"{low['u_g_mean']:.4e}",
                        f"{high['u_g_mean']:.4e}",
                        set_range,
                    ]
                )
                + " |"
            )
        lines.append("")

        lines.extend(["### Selection Stability", ""])
        lines.append("The count below is the number of distinct selected antenna sets across the sigma grid for the same random sample.")
        lines.append("")
        lines.append("| algorithm | mean unique sets/sample | min | max |")
        lines.append("|:---|---:|---:|---:|")
        for _, row in stability[stability["K"] == K].iterrows():
            lines.append(
                "| "
                + " | ".join(
                    [
                        row["heuristic"],
                        f"{row['unique_selected_sets_mean']:.2f}",
                        str(int(row["unique_selected_sets_min"])),
                        str(int(row["unique_selected_sets_max"])),
                    ]
                )
                + " |"
            )
        lines.append("")

    lines.append("Generated plots: `sigma_sweep_K*.png` for metric curves and `sigma_winners_K*.png` for per-sample `U_G` winner shares.")
    (out_dir / "sigma_sweep_report.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    args = parse_args()
    if args.samples <= 0:
        raise ValueError("--samples must be positive.")
    if args.checkpoint_every <= 0:
        raise ValueError("--checkpoint-every must be positive.")
    if args.summary_every < 0:
        raise ValueError("--summary-every must be non-negative.")
    if args.plot_every < 0:
        raise ValueError("--plot-every must be non-negative.")
    K_cases = parse_k_cases(args)
    for case in K_cases:
        K = case["K"]
        if not (0 <= K <= args.N):
            raise ValueError(f"K must satisfy 0 <= K <= N, got K={K}, N={args.N}.")
        if K > 0 and K < args.L:
            raise ValueError(f"K should be at least L for this sweep, got K={K}, L={args.L}.")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.refresh_from_runs is not None:
        runs = pd.read_csv(args.refresh_from_runs)
        write_derived_outputs(runs, args.out_dir, args, include_plots=True)
        print("\nRefreshed from existing runs:")
        print(f"  {args.refresh_from_runs}")
        return

    rows = []
    total_cases = args.samples * len(K_cases) * len(args.sigmas)
    case_no = 0
    for sample in range(args.samples):
        seed = args.seed + sample
        np.random.seed(seed)
        V = generate_V(args.N, args.L)
        for case in K_cases:
            K = case["K"]
            for sigma in args.sigmas:
                case_no += 1
                print(
                    f"[{case_no}/{total_cases}] sample={sample}, seed={seed}, "
                    f"{case['K_label']}, sigma={sigma:g}",
                    flush=True,
                )
                case_rows = run_case(
                    V, K, sigma, args.P, random_state=seed + 1000003 * K
                )
                for row in case_rows:
                    row.update(
                        {
                            "sample": sample,
                            "seed": seed,
                            "N": args.N,
                            "L": args.L,
                            "K": K,
                            "K_mode": case["K_mode"],
                            "K_value": case["K_value"],
                            "active_pct": case["active_pct"],
                            "off_pct": case["off_pct"],
                            "K_label": case["K_label"],
                            "P": args.P,
                        }
                    )
                rows.extend(case_rows)

                progress = {
                    "completed_cases": case_no,
                    "total_cases": total_cases,
                    "completed_fraction": case_no / total_cases,
                    "last_sample": sample,
                    "last_seed": seed,
                    "last_K": K,
                    "last_K_label": case["K_label"],
                    "last_sigma": sigma,
                    "status": "running",
                }
                checkpoint_runs = None
                if case_no % args.checkpoint_every == 0:
                    checkpoint_runs = write_runs_checkpoint(
                        rows, args.out_dir, progress
                    )
                if args.summary_every and case_no % args.summary_every == 0:
                    runs_for_summary = (
                        checkpoint_runs
                        if checkpoint_runs is not None
                        else add_relative_metrics(pd.DataFrame(rows))
                    )
                    write_derived_outputs(
                        runs_for_summary,
                        args.out_dir,
                        args,
                        include_plots=(
                            bool(args.plot_every)
                            and case_no % args.plot_every == 0
                        ),
                    )

    final_progress = {
        "completed_cases": total_cases,
        "total_cases": total_cases,
        "completed_fraction": 1.0,
        "status": "complete",
    }
    runs = write_runs_checkpoint(rows, args.out_dir, final_progress)
    write_derived_outputs(runs, args.out_dir, args, include_plots=True)

    print("\nSaved:")
    paths = [
        "sigma_sweep_runs.csv",
        "sigma_sweep_progress.json",
        "sigma_sweep_summary.csv",
        "sigma_sweep_winners.csv",
        "sigma_sweep_mean_leaders.csv",
        "sigma_sweep_selection_stability.csv",
        "sigma_sweep_report.md",
    ]
    for K in sorted(case["K"] for case in K_cases):
        paths.extend([f"sigma_sweep_K{K}.png", f"sigma_winners_K{K}.png"])
    for path in paths:
        print(f"  {args.out_dir / path}")


if __name__ == "__main__":
    main()
