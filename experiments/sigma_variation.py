import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

from utils.solver_sets import SIGMA_SWEEP_SOLVERS as HEURISTICS
from utils.data import generate_V
from utils.evaluation import evaluate_solver
from utils.io import atomic_write_csv, atomic_write_json
from utils.reporting import (
    format_leader_segments,
    format_number,
    format_sigma,
    leader_segments,
    split_win_shares,
    write_markdown,
)
from visualization.sigma_variation import plot_sigma_sweep


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
        help="Disabled antenna percentages, matching algorithm_comparison --off-pcts.",
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
        x, result = evaluate_solver(heuristic, solver, V, K, sigma, P, random_state)
        active_count = result["active_count"]
        eigvals = effective_channel_eigvals(V, x, P)
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
                "u_bf": result["u_bf"],
                "u_i": result["u_i"],
                "u_g": result["u_g"],
                "log10_u_g": np.log10(max(result["u_g"], np.finfo(float).tiny)),
                "log2_u_g_per_active": np.log2(max(result["u_g"], np.finfo(float).tiny))
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
                "elapsed_seconds": result["elapsed_seconds"],
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
    return split_win_shares(
        runs,
        case_cols=["K", "sample", "sigma"],
        winner_group_cols=["K", "sigma"],
        sample_group_cols=["K", "sigma"],
        algorithm_col="heuristic",
        metric_resolver=lambda case, chunk: METRICS,
        sort_cols=["K", "metric", "sigma", "heuristic"],
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
    write_report(summary, leaders, out_dir, args)
    if include_plots:
        plot_sigma_sweep(summary, winners, out_dir, args, HEURISTICS)

    return summary, winners, leaders, stability


def write_report(summary, leaders, out_dir, args):
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
        "Full data is in the CSV files. This report only summarizes `U_G` leaders.",
        "",
        "## U_G Leader Segments",
        "",
    ]

    for K in sorted(summary["K"].unique()):
        k_label = summary[summary["K"] == K]["K_label"].iloc[0]
        gen_leaders = leaders[(leaders["K"] == K) & (leaders["metric"] == "u_g")]
        segment_text = format_leader_segments(
            leader_segments(gen_leaders),
            format_x=format_sigma,
        )
        end_row = (
            summary[(summary["K"] == K) & (summary["sigma"] == summary[summary["K"] == K]["sigma"].max())]
            .sort_values("u_g_mean", ascending=False)
            .iloc[0]
        )
        lines.append(f"- {k_label}: {segment_text}.")
        lines.append(
            f"  Final sigma leader: {end_row['heuristic']} "
            f"with mean `U_G={format_number(end_row['u_g_mean'])}`."
        )

    lines.extend(
        [
            "",
            "## Artifacts",
            "",
            "- `sigma_sweep_runs.csv`: per-run raw results.",
            "- `sigma_sweep_summary.csv`: means/stds by `K`, `sigma`, algorithm.",
            "- `sigma_sweep_winners.csv`: split winner shares.",
            "- `sigma_sweep_mean_leaders.csv`: mean leaders by metric.",
            "- `sigma_sweep_selection_stability.csv`: selected-set stability.",
            "- `sigma_sweep_K*.png`, `sigma_winners_K*.png`: plots.",
        ]
    )
    write_markdown(out_dir / "sigma_sweep_report.md", lines)


def validate_args(args, K_cases):
    if args.samples <= 0:
        raise ValueError("--samples must be positive.")
    if args.checkpoint_every <= 0:
        raise ValueError("--checkpoint-every must be positive.")
    if args.summary_every < 0:
        raise ValueError("--summary-every must be non-negative.")
    if args.plot_every < 0:
        raise ValueError("--plot-every must be non-negative.")
    for case in K_cases:
        K = case["K"]
        if not (0 <= K <= args.N):
            raise ValueError(f"K must satisfy 0 <= K <= N, got K={K}, N={args.N}.")
        if K > 0 and K < args.L:
            raise ValueError(f"K should be at least L for this sweep, got K={K}, L={args.L}.")


def refresh_outputs(args):
    runs = pd.read_csv(args.refresh_from_runs)
    write_derived_outputs(runs, args.out_dir, args, include_plots=True)
    print("\nRefreshed from existing runs:")
    print(f"  {args.refresh_from_runs}")


def progress_payload(case_no, total_cases, sample, seed, case, sigma, status="running"):
    payload = {
        "completed_cases": case_no,
        "total_cases": total_cases,
        "completed_fraction": case_no / total_cases,
        "status": status,
    }
    if status == "running":
        payload.update(
            {
                "last_sample": sample,
                "last_seed": seed,
                "last_K": case["K"],
                "last_K_label": case["K_label"],
                "last_sigma": sigma,
            }
        )
    return payload


def add_case_metadata(case_rows, sample, seed, args, case):
    for row in case_rows:
        row.update(
            {
                "sample": sample,
                "seed": seed,
                "N": args.N,
                "L": args.L,
                "K": case["K"],
                "K_mode": case["K_mode"],
                "K_value": case["K_value"],
                "active_pct": case["active_pct"],
                "off_pct": case["off_pct"],
                "K_label": case["K_label"],
                "P": args.P,
            }
        )
    return case_rows


def maybe_write_intermediate_outputs(rows, args, progress, case_no):
    checkpoint_runs = None
    if case_no % args.checkpoint_every == 0:
        checkpoint_runs = write_runs_checkpoint(rows, args.out_dir, progress)
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
            include_plots=bool(args.plot_every) and case_no % args.plot_every == 0,
        )


def run_sweep(args, K_cases):
    rows = []
    total_cases = args.samples * len(K_cases) * len(args.sigmas)
    case_no = 0

    for sample in range(args.samples):
        seed = args.seed + sample
        np.random.seed(seed)
        V = generate_V(args.N, args.L)
        for case in K_cases:
            for sigma in args.sigmas:
                case_no += 1
                print(
                    f"[{case_no}/{total_cases}] sample={sample}, seed={seed}, "
                    f"{case['K_label']}, sigma={sigma:g}",
                    flush=True,
                )
                case_rows = run_case(
                    V,
                    case["K"],
                    sigma,
                    args.P,
                    random_state=seed + 1000003 * case["K"],
                )
                rows.extend(add_case_metadata(case_rows, sample, seed, args, case))
                progress = progress_payload(case_no, total_cases, sample, seed, case, sigma)
                maybe_write_intermediate_outputs(rows, args, progress, case_no)

    final_progress = {
        "completed_cases": total_cases,
        "total_cases": total_cases,
        "completed_fraction": 1.0,
        "status": "complete",
    }
    runs = write_runs_checkpoint(rows, args.out_dir, final_progress)
    write_derived_outputs(runs, args.out_dir, args, include_plots=True)


def print_saved_paths(args, K_cases):
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


def main():
    args = parse_args()
    K_cases = parse_k_cases(args)
    validate_args(args, K_cases)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.refresh_from_runs is not None:
        refresh_outputs(args)
        return

    run_sweep(args, K_cases)
    print_saved_paths(args, K_cases)


if __name__ == "__main__":
    main()
