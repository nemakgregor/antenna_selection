import argparse
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from algorithms.common import objective_from_gram
from algorithms.h3_threshold_explore import (
    build_threshold_grid,
    dense_thresholds,
    evaluate_threshold_T,
    evaluate_power_window_thresholds,
    row_power_distribution_metrics,
)
try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None

from utils.solver_sets import CDF_SOLVERS
from utils.data import DATA_PROFILES, generate_v_from_rng, generate_v_profile_from_rng
from utils.evaluation import evaluate_solver
from utils.brute_force import (
    brute_force_exact_u_g,
    contiguous_threshold_window_T,
    subset_to_string,
    threshold_window_subset_string,
)
from utils.io import atomic_write_csv
from utils.local_threshold_analysis import run_local_threshold_exact_analysis
from utils.reporting import format_number_slug
from utils.threshold_exact_analysis import write_threshold_exact_k_analysis
from visualization.algorithm_comparison import (
    write_algorithm_comparison_plots,
    write_threshold_exploration_plots,
)


OUR_ALGORITHMS = (
    "FrameOnly-Gen",
    "CapWindow-Gen",
    "CapSubmod-Gen",
    "CapSubmodPort-Gen",
    "ThreshDOpt-Gen",
    "ThreshWLogdet-Gen",
    "ThreshDOptSwap-Gen",
    "S-threshold-Gen",
    "BackwardTrueGreedy",
)
FOCUSED_H3_CAP_WINDOW = (
    "H3",
    "FrameOnly-Gen",
    "CapWindow-Gen",
    "CapSubmod-Gen",
    "CapSubmodPort-Gen",
    "ThreshDOpt-Gen",
    "ThreshWLogdet-Gen",
    "ThreshDOptSwap-Gen",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build empirical CDF plots for U_G in dB and solver runtime on the "
            "antenna-selection task."
        )
    )
    parser.add_argument("--N", type=int, default=1000)
    parser.add_argument("--L", type=int, default=2)
    parser.add_argument(
        "--samples",
        type=int,
        default=1000,
        help="Number of generated matrices per generator seed.",
    )
    parser.add_argument(
        "--generator-seeds",
        type=int,
        nargs="+",
        default=[10, 42],
        help="Independent RandomState seeds used to generate sample streams.",
    )
    parser.add_argument(
        "--off-pcts",
        type=float,
        nargs="+",
        default=[25.0, 50.0],
        help="Percent of antennas switched off. K = round(N * (1 - off_pct/100)).",
    )
    parser.add_argument(
        "--off-counts",
        type=int,
        nargs="+",
        default=None,
        help="Absolute number of antennas switched off. Overrides --off-pcts.",
    )
    parser.add_argument(
        "--K-values",
        type=int,
        nargs="+",
        default=None,
        help="Explicit active antenna limits. Overrides --off-pcts.",
    )
    parser.add_argument("--sigma", type=float, default=1.0)
    parser.add_argument("--P", type=float, default=1.0)
    parser.add_argument(
        "--algorithms",
        nargs="+",
        default=None,
        help="Optional subset of algorithm names. Default: every registered comparison algorithm.",
    )
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=25,
        help="Write cdf_runs.csv every N newly completed algorithm runs. Use 0 for final only.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Skip runs already present in cdf_runs.csv under --out-dir.",
    )
    parser.add_argument(
        "--rerun-algorithms",
        nargs="+",
        default=None,
        help="When resuming, drop existing rows for these algorithms and recompute them.",
    )
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="Read cdf_runs.csv and rebuild summary/plots without running algorithms.",
    )
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument(
        "--threshold-explore",
        action="store_true",
        help="Run threshold-value exploration instead of algorithm CDF comparison.",
    )
    parser.add_argument(
        "--threshold-full-sweep",
        action="store_true",
        help="Run dense T=0..K power-window threshold sweep instead of algorithm CDF comparison.",
    )
    parser.add_argument(
        "--threshold-rule-cdf",
        action="store_true",
        help=(
            "Build CDF/report comparing best formula, per-sample best tested T, "
            "and H3 strong/weak from an existing --threshold-full-sweep result."
        ),
    )
    parser.add_argument(
        "--threshold-scaling-study",
        action="store_true",
        help=(
            "Run preliminary dense threshold sweeps across multiple N/L values "
            "and active-K percentages."
        ),
    )
    parser.add_argument(
        "--threshold-exact-study",
        action="store_true",
        help=(
            "Run small-N exact brute-force U_G comparison against threshold windows."
        ),
    )
    parser.add_argument(
        "--threshold-local-exact-analysis",
        action="store_true",
        help=(
            "Use saved exact Gaussian cases to compare pure threshold windows "
            "against one-swap and two-swap local threshold refinements."
        ),
    )
    parser.add_argument(
        "--exact-source-dir",
        type=Path,
        default=Path(
            "results/threshold_exact_gaussian_L2_N8_12_16_20_24_Kpct25_to_50_s100"
        ),
        help="Exact-study result directory used by --threshold-local-exact-analysis.",
    )
    parser.add_argument(
        "--N-values",
        type=int,
        nargs="+",
        default=[500, 1000, 2000, 5000],
        help="N grid used by --threshold-scaling-study.",
    )
    parser.add_argument(
        "--L-values",
        type=int,
        nargs="+",
        default=[4, 6, 8, 10],
        help="L grid used by --threshold-scaling-study.",
    )
    parser.add_argument(
        "--K-pcts",
        type=float,
        nargs="+",
        default=[25.0, 50.0],
        help=(
            "Active antenna percentages used by --threshold-scaling-study. "
            "For example 25 means K=round(0.25*N)."
        ),
    )
    parser.add_argument(
        "--data-profiles",
        nargs="+",
        choices=DATA_PROFILES,
        default=["gaussian"],
        help="Synthetic channel profiles used by --threshold-explore.",
    )
    parser.add_argument(
        "--exact-time-limit",
        type=float,
        default=120.0,
        help="Per-case brute-force exact enumeration guard in seconds.",
    )
    parser.add_argument(
        "--threshold-target",
        choices=("bf", "int", "gen"),
        default="gen",
        help="Objective used to score each explicit threshold in --threshold-explore.",
    )
    return parser.parse_args()


def default_out_dir(args):
    if args.threshold_local_exact_analysis:
        return Path("results/local_threshold_exact_gauss_L2_N8_12_16_20_Kpct25_to_50_s100")

    if args.threshold_exact_study:
        n_label = "_".join(str(value) for value in args.N_values)
        k_label = "_".join(format_number_slug(value) for value in args.K_pcts)
        seed_label = "_".join(str(value) for value in args.generator_seeds)
        profile_label = "_".join(args.data_profiles)
        return Path(
            f"results/threshold_exact_{profile_label}_L{args.L}_N{n_label}_"
            f"Kpct{k_label}_s{args.samples}_seeds{seed_label}"
        )

    if args.threshold_scaling_study:
        n_label = "_".join(str(value) for value in args.N_values)
        l_label = "_".join(str(value) for value in args.L_values)
        k_label = "_".join(format_number_slug(value) for value in args.K_pcts)
        return Path(
            f"results/threshold_scaling_prelim_N{n_label}_L{l_label}_"
            f"Kpct{k_label}_s{args.samples}"
        )

    if args.threshold_full_sweep or args.threshold_rule_cdf:
        if args.K_values is not None:
            k_label = "K" + "_".join(str(value) for value in args.K_values)
        elif args.off_counts is not None:
            k_label = "offcount" + "_".join(str(value) for value in args.off_counts)
        else:
            k_label = "off" + "_".join(format_number_slug(value) for value in args.off_pcts)
        return Path(f"results/threshold_full_sweep_L{args.L}_N{args.N}_{k_label}")

    if args.threshold_explore:
        if args.off_counts is not None:
            off_label = "count" + "_".join(str(value) for value in args.off_counts)
        elif args.K_values is not None:
            off_label = "K" + "_".join(str(value) for value in args.K_values)
        else:
            off_label = "_".join(format_number_slug(value) for value in args.off_pcts)
        profile_label = "_".join(args.data_profiles)
        seed_label = "_".join(str(value) for value in args.generator_seeds)
        return Path(
            f"results/threshold_exploration_N{args.N}_L{args.L}_off{off_label}_"
            f"profiles{profile_label}_seeds{seed_label}_{args.samples}samples"
        )

    if args.off_counts is not None:
        off_label = "count" + "_".join(str(value) for value in args.off_counts)
    else:
        off_label = "_".join(format_number_slug(value) for value in args.off_pcts)
    seed_label = "_".join(str(value) for value in args.generator_seeds)
    return Path(
        f"results/cdf_N{args.N}_L{args.L}_off{off_label}_"
        f"seeds{seed_label}_{args.samples}samples"
    )


def build_off_cases(args):
    if args.K_values is not None and args.off_counts is not None:
        raise ValueError("Use only one of --K-values or --off-counts.")

    cases = []
    if args.K_values is not None:
        for K in args.K_values:
            if not (0 <= K <= args.N):
                raise ValueError("Each K must satisfy 0 <= K <= N.")
            off_pct = 100.0 * float(args.N - K) / float(args.N)
            cases.append({"K": int(K), "off_pct": off_pct})
        return cases

    if args.off_counts is not None:
        for off_count in args.off_counts:
            if not (0 <= off_count < args.N):
                raise ValueError("Each off_count must satisfy 0 <= off_count < N.")
            K = int(args.N - off_count)
            cases.append(
                {
                    "K": K,
                    "off_pct": 100.0 * float(off_count) / float(args.N),
                }
            )
        return cases

    for off_pct in args.off_pcts:
        if not (0.0 <= off_pct < 100.0):
            raise ValueError("Each off_pct must satisfy 0 <= off_pct < 100.")
        K = int(round(args.N * (1.0 - off_pct / 100.0)))
        cases.append({"K": K, "off_pct": float(off_pct)})
    return cases


def build_active_pct_cases(N, active_pcts):
    cases = []
    for active_pct in active_pcts:
        if not (0.0 < active_pct <= 100.0):
            raise ValueError("Each active K percentage must satisfy 0 < pct <= 100.")
        K = int(round(N * active_pct / 100.0))
        if not (0 <= K <= N):
            raise ValueError("Each active K percentage must produce 0 <= K <= N.")
        if K > N - K:
            raise ValueError(
                "The dense scaling study sweeps T=0..K, which requires K <= N-K. "
                f"Got N={N}, active_pct={active_pct:g}, K={K}."
            )
        cases.append(
            {
                "K": K,
                "off_pct": float(100.0 - active_pct),
                "active_pct": float(active_pct),
            }
        )
    return cases


def load_existing_runs(path):
    if not path.exists():
        return pd.DataFrame(), set()

    runs = pd.read_csv(path)
    keys = set(
        zip(
            runs["generator_seed"].astype(int),
            runs["sample"].astype(int),
            runs["off_pct"].astype(float),
            runs["algorithm"].astype(str),
        )
    )
    return runs, keys


def select_algorithms(all_algorithms, args):
    if args.algorithms is None:
        return all_algorithms

    selected_names = args.algorithms
    available = {name for name, _ in all_algorithms}
    unknown = sorted(set(selected_names) - available)
    if unknown:
        raise ValueError(f"Unknown algorithms: {', '.join(unknown)}")

    selected = set(selected_names)
    return tuple((name, solver) for name, solver in all_algorithms if name in selected)


def solver_random_state(generator_seed, sample, K, algorithm_index):
    return int(
        (
            generator_seed * 1_000_003
            + sample * 10_007
            + K * 389
            + algorithm_index * 97
        )
        % (2**32 - 1)
    )


def run_algorithm(name, solver, V, K, off_pct, args, random_state):
    _, result = evaluate_solver(
        name,
        solver,
        V,
        K,
        args.sigma,
        args.P,
        random_state,
    )
    u_g_safe = max(result["u_g"], np.finfo(float).tiny)
    return {
        "algorithm": name,
        "active_count": result["active_count"],
        "u_bf": result["u_bf"],
        "u_i": result["u_i"],
        "u_g": result["u_g"],
        "u_g_db": float(10.0 * np.log10(u_g_safe)),
        "elapsed_seconds": result["elapsed_seconds"],
    }


def run_benchmark(args, algorithms, runs_path):
    existing, all_completed = (
        load_existing_runs(runs_path) if args.resume else (pd.DataFrame(), set())
    )
    rows = existing.to_dict("records") if not existing.empty else []
    rerun_algorithms = set(args.rerun_algorithms or [])
    if rerun_algorithms:
        rows = [row for row in rows if row["algorithm"] not in rerun_algorithms]
        all_completed = {key for key in all_completed if key[3] not in rerun_algorithms}

    new_since_checkpoint = 0
    total_new = 0

    selected_names = [name for name, _ in algorithms]
    selected_name_set = set(selected_names)
    completed = {
        key for key in all_completed
        if key[3] in selected_name_set
    }
    total_cases = (
        len(args.generator_seeds)
        * args.samples
        * len(args.off_cases)
        * len(selected_names)
    )
    case_no = len(completed)
    progress = (
        tqdm(
            total=total_cases,
            initial=min(case_no, total_cases),
            unit="run",
            dynamic_ncols=True,
        )
        if tqdm is not None
        else None
    )

    try:
        for generator_seed in args.generator_seeds:
            rng = np.random.RandomState(generator_seed)
            for sample in range(args.samples):
                V = generate_v_from_rng(rng, args.N, args.L)
                for off_case in args.off_cases:
                    off_pct = off_case["off_pct"]
                    K = off_case["K"]

                    for algorithm_index, (name, solver) in enumerate(algorithms):
                        key = (int(generator_seed), int(sample), float(off_pct), name)
                        if key in completed:
                            continue

                        case_no += 1
                        message = (
                            f"seed={generator_seed}, sample={sample}, "
                            f"off={off_pct:g}%, K={K}, algorithm={name}"
                        )
                        if progress is not None:
                            progress.set_postfix_str(message)
                        else:
                            print(f"[{case_no}/{total_cases}] {message}", flush=True)
                        random_state = solver_random_state(
                            generator_seed, sample, K, algorithm_index
                        )
                        result = run_algorithm(
                            name, solver, V, K, off_pct, args, random_state
                        )
                        rows.append(
                            {
                                "generator_seed": int(generator_seed),
                                "sample": int(sample),
                                "N": int(args.N),
                                "L": int(args.L),
                                "K": int(K),
                                "off_pct": float(off_pct),
                                "active_pct": float(100.0 - off_pct),
                                "sigma": float(args.sigma),
                                "P": float(args.P),
                                **result,
                            }
                        )
                        total_new += 1
                        new_since_checkpoint += 1

                        if progress is not None:
                            progress.update(1)

                        if (
                            args.checkpoint_every
                            and new_since_checkpoint >= args.checkpoint_every
                        ):
                            atomic_write_csv(pd.DataFrame(rows), runs_path)
                            new_since_checkpoint = 0
    finally:
        if progress is not None:
            progress.close()

    runs = pd.DataFrame(rows)
    atomic_write_csv(runs, runs_path)
    print(f"Completed {total_new} new algorithm runs.", flush=True)
    return runs


def build_summary(runs):
    def q(value):
        return lambda data: data.quantile(value)

    summary = (
        runs.groupby(["off_pct", "K", "algorithm"], as_index=False)
        .agg(
            samples=("u_g", "count"),
            u_g_db_mean=("u_g_db", "mean"),
            u_g_db_p05=("u_g_db", q(0.05)),
            u_g_db_p50=("u_g_db", q(0.50)),
            u_g_db_p95=("u_g_db", q(0.95)),
            elapsed_mean=("elapsed_seconds", "mean"),
            elapsed_p05=("elapsed_seconds", q(0.05)),
            elapsed_p50=("elapsed_seconds", q(0.50)),
            elapsed_p95=("elapsed_seconds", q(0.95)),
        )
        .sort_values(["off_pct", "algorithm"])
    )
    return summary


def summarize_values(values):
    values = pd.Series(values, dtype=float).replace([np.inf, -np.inf], np.nan).dropna()
    if values.empty:
        return {
            "min": np.nan,
            "avg": np.nan,
            "max": np.nan,
            "std": np.nan,
        }
    return {
        "min": float(values.min()),
        "avg": float(values.mean()),
        "max": float(values.max()),
        "std": float(values.std(ddof=1)) if len(values) > 1 else 0.0,
    }


def build_baseline_improvement(runs, baseline_names=("H1", "H2", "H3")):
    case_cols = ["generator_seed", "sample", "off_pct", "K"]
    algorithms = list(dict.fromkeys(runs["algorithm"]))
    missing = [name for name in baseline_names if name not in algorithms]
    if missing:
        raise RuntimeError(
            "Baseline algorithms are required for improvement statistics: "
            + ", ".join(missing)
        )

    metric_table = runs.pivot_table(
        index=case_cols,
        columns="algorithm",
        values=["u_g", "u_g_db"],
        aggfunc="first",
    )

    rows = []
    baseline_specs = list(baseline_names) + ["best_H123"]
    grouped_cases = runs[case_cols].drop_duplicates().groupby(
        ["generator_seed", "off_pct", "K"], sort=True
    )

    for group_key, group_cases in grouped_cases:
        generator_seed, off_pct, K = group_key
        index = pd.MultiIndex.from_frame(group_cases[case_cols])
        group_metrics = metric_table.reindex(index)
        baseline_u_g_values = group_metrics["u_g"].loc[:, list(baseline_names)]
        baseline_u_g_db_values = group_metrics["u_g_db"].loc[:, list(baseline_names)]

        for algorithm in algorithms:
            if algorithm not in group_metrics["u_g"]:
                continue
            algorithm_u_g = group_metrics["u_g"][algorithm]
            algorithm_u_g_db = group_metrics["u_g_db"][algorithm]

            for baseline in baseline_specs:
                if baseline == "best_H123":
                    baseline_u_g = baseline_u_g_values.max(axis=1)
                    baseline_u_g_db = baseline_u_g_db_values.max(axis=1)
                else:
                    baseline_u_g = group_metrics["u_g"][baseline]
                    baseline_u_g_db = group_metrics["u_g_db"][baseline]

                valid = (
                    algorithm_u_g.notna()
                    & algorithm_u_g_db.notna()
                    & baseline_u_g.notna()
                    & baseline_u_g_db.notna()
                    & (baseline_u_g > 0.0)
                    & (baseline_u_g_db != 0.0)
                )
                if not valid.any():
                    continue

                u_g_delta = algorithm_u_g[valid] - baseline_u_g[valid]
                u_g_pct = 100.0 * (algorithm_u_g[valid] / baseline_u_g[valid] - 1.0)
                u_g_db_delta = algorithm_u_g_db[valid] - baseline_u_g_db[valid]
                u_g_db_pct = 100.0 * (
                    algorithm_u_g_db[valid] / baseline_u_g_db[valid] - 1.0
                )

                row = {
                    "generator_seed": int(generator_seed),
                    "off_pct": float(off_pct),
                    "K": int(K),
                    "algorithm": algorithm,
                    "baseline": baseline,
                    "cases": int(valid.sum()),
                }
                for prefix, values in (
                    ("u_g_delta", u_g_delta),
                    ("u_g_pct", u_g_pct),
                    ("u_g_db_delta", u_g_db_delta),
                    ("u_g_db_pct", u_g_db_pct),
                ):
                    stats = summarize_values(values)
                    for stat_name, stat_value in stats.items():
                        row[f"{prefix}_{stat_name}"] = stat_value
                rows.append(row)

    return pd.DataFrame(rows).sort_values(
        ["generator_seed", "off_pct", "baseline", "algorithm"]
    )


def format_float(value, precision=3):
    if pd.isna(value):
        return ""
    if abs(value) >= 1e6 or (0 < abs(value) < 1e-3):
        return f"{value:.{precision}e}"
    return f"{value:.{precision}f}"


def write_improvement_report(improvement, out_path):
    focus = improvement[improvement["baseline"] == "best_H123"].copy()
    lines = [
        "# CDF objective improvement vs best(H1,H2,H3)",
        "",
        "Positive values mean the algorithm has a larger `U_G` than the baseline on the same generated case.",
        "`u_g_pct_*` is percent improvement in absolute `U_G`; `u_g_db_delta_*` is the dB difference.",
        "",
    ]

    for (generator_seed, off_pct, K), chunk in focus.groupby(
        ["generator_seed", "off_pct", "K"], sort=True
    ):
        chunk = chunk.sort_values("u_g_pct_avg", ascending=False)
        lines.extend(
            [
                f"## seed={int(generator_seed)}, {off_pct:g}% off, K={int(K)}",
                "",
                "| algorithm | cases | U_G % avg | U_G % min | U_G % max | U_G % std | dB avg | dB min | dB max | dB std | abs U_G avg |",
                "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for _, row in chunk.iterrows():
            lines.append(
                "| "
                + " | ".join(
                    [
                        str(row["algorithm"]),
                        str(int(row["cases"])),
                        format_float(row["u_g_pct_avg"]),
                        format_float(row["u_g_pct_min"]),
                        format_float(row["u_g_pct_max"]),
                        format_float(row["u_g_pct_std"]),
                        format_float(row["u_g_db_delta_avg"]),
                        format_float(row["u_g_db_delta_min"]),
                        format_float(row["u_g_db_delta_max"]),
                        format_float(row["u_g_db_delta_std"]),
                        format_float(row["u_g_delta_avg"]),
                    ]
                )
                + " |"
            )
        lines.append("")

    out_path.write_text("\n".join(lines), encoding="utf-8")


def write_our_vs_h123_report(
    improvement,
    out_dir,
    our_algorithms=OUR_ALGORITHMS,
    baselines=("H1", "H2", "H3"),
):
    report = improvement[
        improvement["algorithm"].isin(our_algorithms)
        & improvement["baseline"].isin(baselines)
    ].copy()
    if report.empty:
        return

    report = report.sort_values(
        ["generator_seed", "off_pct", "algorithm", "baseline"]
    )
    atomic_write_csv(report, out_dir / "cdf_our_vs_h123.csv")

    lines = [
        "# Our algorithms vs H1/H2/H3",
        "",
        "Each row compares one of our algorithms with one baseline heuristic on the same generated cases.",
        "Positive values mean our algorithm has a larger target objective `U_G`.",
        "`U_G %` is percent improvement in absolute `U_G`; `dB %` is percent improvement in `10 lg(U_G)`; `dB delta` is `10 lg(U_G_alg) - 10 lg(U_G_baseline)`; `abs U_G delta` is `U_G_alg - U_G_baseline`.",
        "",
    ]

    for (generator_seed, off_pct, K), chunk in report.groupby(
        ["generator_seed", "off_pct", "K"], sort=True
    ):
        lines.extend(
            [
                f"## seed={int(generator_seed)}, {off_pct:g}% off, K={int(K)}",
                "",
                "| our algorithm | vs heuristic | cases | U_G % avg | U_G % min | U_G % max | U_G % std | dB % avg | dB % min | dB % max | dB % std | dB delta avg | dB delta min | dB delta max | dB delta std | abs U_G delta avg | abs U_G delta min | abs U_G delta max | abs U_G delta std |",
                "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for _, row in chunk.iterrows():
            lines.append(
                "| "
                + " | ".join(
                    [
                        str(row["algorithm"]),
                        str(row["baseline"]),
                        str(int(row["cases"])),
                        format_float(row["u_g_pct_avg"]),
                        format_float(row["u_g_pct_min"]),
                        format_float(row["u_g_pct_max"]),
                        format_float(row["u_g_pct_std"]),
                        format_float(row["u_g_db_pct_avg"]),
                        format_float(row["u_g_db_pct_min"]),
                        format_float(row["u_g_db_pct_max"]),
                        format_float(row["u_g_db_pct_std"]),
                        format_float(row["u_g_db_delta_avg"]),
                        format_float(row["u_g_db_delta_min"]),
                        format_float(row["u_g_db_delta_max"]),
                        format_float(row["u_g_db_delta_std"]),
                        format_float(row["u_g_delta_avg"]),
                        format_float(row["u_g_delta_min"]),
                        format_float(row["u_g_delta_max"]),
                        format_float(row["u_g_delta_std"]),
                    ]
                )
                + " |"
            )
        lines.append("")

    (out_dir / "cdf_our_vs_h123.md").write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def write_outputs(runs, algorithms, out_dir):
    selected = {name for name, _ in algorithms}
    runs = runs[runs["algorithm"].isin(selected)].copy()
    if runs.empty:
        raise RuntimeError("No runs are available for the selected algorithms.")

    summary = build_summary(runs)
    atomic_write_csv(summary, out_dir / "cdf_summary.csv")
    baseline_names = ("H1", "H2", "H3")
    if set(baseline_names).issubset(set(runs["algorithm"])):
        improvement = build_baseline_improvement(runs, baseline_names=baseline_names)
        atomic_write_csv(improvement, out_dir / "cdf_baseline_improvement.csv")
        write_improvement_report(improvement, out_dir / "cdf_baseline_improvement.md")
        write_our_vs_h123_report(improvement, out_dir)
    write_algorithm_comparison_plots(
        runs,
        algorithms,
        out_dir,
        focused_names=FOCUSED_H3_CAP_WINDOW,
    )


THRESHOLD_CASE_COLS = [
    "data_profile",
    "generator_seed",
    "sample",
    "off_pct",
    "K",
    "sigma",
    "P",
]

THRESHOLD_METRIC_COLS = [
    "row_power_cv",
    "row_power_skew3",
    "row_power_p95_p50",
    "row_power_p99_p50",
    "row_power_max_p95",
    "log_power_gap_max",
    "log_power_gap_rank",
    "log_power_gap_rank_pct",
    "tail_mass_p80",
    "tail_mass_p90",
    "tail_mass_p95",
    "tail_mass_p99",
]

FULL_SWEEP_CASE_COLS = [
    "data_profile",
    "generator_seed",
    "sample",
    "N",
    "L",
    "K",
    "sigma",
    "P",
]

FORMULA_SPECS = (
    ("T_0p05N", lambda N, L, K: 0.05 * N),
    ("T_0p075N", lambda N, L, K: 0.075 * N),
    ("T_0p10N", lambda N, L, K: 0.10 * N),
    ("T_0p10K", lambda N, L, K: 0.10 * K),
    ("T_0p20K", lambda N, L, K: 0.20 * K),
    ("legacy_T25", lambda N, L, K: 25),
    ("legacy_T50", lambda N, L, K: 50),
    ("legacy_T100", lambda N, L, K: 100),
)

SCALING_FORMULA_SPECS = (
    ("T_0p025N", lambda N, L, K: 0.025 * N),
    *FORMULA_SPECS,
    ("T_0p05K", lambda N, L, K: 0.05 * K),
    ("T_0p15K", lambda N, L, K: 0.15 * K),
    ("T_0p05NL_over_Lp2", lambda N, L, K: 0.05 * N * L / (L + 2.0)),
    ("T_0p075NL_over_Lp2", lambda N, L, K: 0.075 * N * L / (L + 2.0)),
    ("T_0p10NL_over_Lp2", lambda N, L, K: 0.10 * N * L / (L + 2.0)),
    ("T_0p125NL_over_Lp2", lambda N, L, K: 0.125 * N * L / (L + 2.0)),
    ("T_0p15NL_over_Lp2", lambda N, L, K: 0.15 * N * L / (L + 2.0)),
)

STRONG_WEAK_RULE = "strong_weak"


def threshold_formula_T(label, N, L, K, formula_specs=SCALING_FORMULA_SPECS):
    for formula_label, resolver in formula_specs:
        if formula_label == label:
            return int(np.clip(int(round(resolver(N, L, K))), 0, K))
    raise ValueError(f"Unknown threshold formula: {label}")


def formula_family(label):
    if label == STRONG_WEAK_RULE:
        return "heuristic"
    if "NL_over_Lp2" in label:
        return "N_L"
    if label.startswith("T_") and label.endswith("K"):
        return "K"
    if label.startswith("T_") and label.endswith("N"):
        return "N"
    if label.startswith("legacy_"):
        return "legacy"
    return "other"


def run_threshold_exploration(args, runs_path):
    rows = []
    new_since_checkpoint = 0
    total_cases = (
        len(args.data_profiles)
        * len(args.generator_seeds)
        * args.samples
        * len(args.off_cases)
    )
    case_no = 0
    progress = (
        tqdm(total=total_cases, unit="case", dynamic_ncols=True)
        if tqdm is not None
        else None
    )

    try:
        for profile in args.data_profiles:
            for generator_seed in args.generator_seeds:
                rng = np.random.RandomState(generator_seed)
                for sample in range(args.samples):
                    V = generate_v_profile_from_rng(rng, args.N, args.L, profile)
                    distribution_metrics = row_power_distribution_metrics(V)

                    for off_case in args.off_cases:
                        case_no += 1
                        K = off_case["K"]
                        off_pct = off_case["off_pct"]
                        threshold_grid = build_threshold_grid(V, K)
                        if progress is not None:
                            progress.set_postfix_str(
                                f"profile={profile}, seed={generator_seed}, "
                                f"sample={sample}, off={off_pct:g}%"
                            )
                        else:
                            print(
                                f"[{case_no}/{total_cases}] profile={profile}, "
                                f"seed={generator_seed}, sample={sample}, "
                                f"off={off_pct:g}%, K={K}, thresholds={len(threshold_grid)}",
                                flush=True,
                            )

                        for threshold in threshold_grid:
                            T = int(threshold["T"])
                            result = evaluate_threshold_T(
                                V,
                                K,
                                T,
                                target_obj=args.threshold_target,
                                sigma=args.sigma,
                                P=args.P,
                            )
                            x = result["x"]
                            u_g_safe = max(result["u_g"], np.finfo(float).tiny)
                            rows.append(
                                {
                                    "data_profile": profile,
                                    "generator_seed": int(generator_seed),
                                    "sample": int(sample),
                                    "N": int(args.N),
                                    "L": int(args.L),
                                    "K": int(K),
                                    "off_pct": float(off_pct),
                                    "active_pct": float(100.0 - off_pct),
                                    "sigma": float(args.sigma),
                                    "P": float(args.P),
                                    "target_obj": args.threshold_target,
                                    "T": T,
                                    "threshold_source": threshold["threshold_source"],
                                    "candidate_kind": result["candidate_kind"],
                                    "candidate_count": result["candidate_count"],
                                    "active_count": int(np.sum(x)),
                                    "u_bf": result["u_bf"],
                                    "u_i": result["u_i"],
                                    "u_g": result["u_g"],
                                    "u_g_db": float(10.0 * np.log10(u_g_safe)),
                                    "score": result["score"],
                                    "elapsed_seconds": result["elapsed_seconds"],
                                    **distribution_metrics,
                                }
                            )
                            new_since_checkpoint += 1

                            if (
                                args.checkpoint_every
                                and new_since_checkpoint >= args.checkpoint_every
                            ):
                                checkpoint = add_threshold_relative_metrics(
                                    pd.DataFrame(rows)
                                )
                                atomic_write_csv(checkpoint, runs_path)
                                new_since_checkpoint = 0

                        if progress is not None:
                            progress.update(1)
    finally:
        if progress is not None:
            progress.close()

    runs = add_threshold_relative_metrics(pd.DataFrame(rows))
    atomic_write_csv(runs, runs_path)
    return runs


def add_threshold_relative_metrics(runs):
    if runs.empty:
        return runs
    runs = runs.copy()
    runs["best_u_g_for_case"] = runs.groupby(THRESHOLD_CASE_COLS)["u_g"].transform(
        "max"
    )
    runs["u_g_vs_best_T"] = runs["u_g"] / np.maximum(
        runs["best_u_g_for_case"], np.finfo(float).eps
    )
    runs["u_g_gap_to_best_pct"] = 100.0 * (1.0 - runs["u_g_vs_best_T"])
    runs["is_best_T"] = np.isclose(runs["u_g"], runs["best_u_g_for_case"])
    return runs


def build_threshold_summary(runs):
    def q(value):
        return lambda data: data.quantile(value)

    return (
        runs.groupby(["data_profile", "off_pct", "K", "T"], as_index=False)
        .agg(
            threshold_source=("threshold_source", _join_threshold_sources),
            samples=("u_g", "count"),
            winner_rate=("is_best_T", "mean"),
            u_g_mean=("u_g", "mean"),
            u_g_p05=("u_g", q(0.05)),
            u_g_p50=("u_g", q(0.50)),
            u_g_p95=("u_g", q(0.95)),
            u_g_db_mean=("u_g_db", "mean"),
            u_g_vs_best_T_mean=("u_g_vs_best_T", "mean"),
            u_g_vs_best_T_p05=("u_g_vs_best_T", q(0.05)),
            u_g_gap_to_best_pct_mean=("u_g_gap_to_best_pct", "mean"),
            elapsed_seconds_mean=("elapsed_seconds", "mean"),
            candidate_kind_mode=("candidate_kind", _mode_value),
        )
        .sort_values(
            [
                "data_profile",
                "off_pct",
                "u_g_vs_best_T_mean",
                "winner_rate",
                "T",
            ],
            ascending=[True, True, False, False, True],
        )
    )


def build_threshold_metric_correlations(runs):
    if runs.empty:
        return pd.DataFrame()

    case_best = (
        runs.sort_values(
            [*THRESHOLD_CASE_COLS, "u_g", "T"],
            ascending=[True] * len(THRESHOLD_CASE_COLS) + [False, True],
        )
        .drop_duplicates(THRESHOLD_CASE_COLS)
        .copy()
    )
    rows = []
    for (profile, off_pct, K), chunk in case_best.groupby(
        ["data_profile", "off_pct", "K"], sort=True
    ):
        best_T = chunk["T"].astype(float)
        for metric in THRESHOLD_METRIC_COLS:
            values = chunk[metric].astype(float)
            if len(chunk) < 2 or best_T.std(ddof=0) == 0 or values.std(ddof=0) == 0:
                corr = np.nan
            else:
                corr = float(np.corrcoef(values, best_T)[0, 1])
            rows.append(
                {
                    "data_profile": profile,
                    "off_pct": float(off_pct),
                    "K": int(K),
                    "metric": metric,
                    "pearson_corr_with_best_T": corr,
                    "cases": int(len(chunk)),
                    "best_T_mean": float(best_T.mean()),
                    "best_T_std": float(best_T.std(ddof=1)) if len(chunk) > 1 else 0.0,
                    "best_T_min": int(best_T.min()),
                    "best_T_max": int(best_T.max()),
                }
            )
    return pd.DataFrame(rows).sort_values(
        ["data_profile", "off_pct", "K", "metric"]
    )


def write_threshold_exploration_outputs(runs, out_dir, args):
    summary = build_threshold_summary(runs)
    correlations = build_threshold_metric_correlations(runs)
    atomic_write_csv(summary, out_dir / "threshold_summary.csv")
    atomic_write_csv(correlations, out_dir / "threshold_metric_correlations.csv")
    write_threshold_exploration_report(runs, summary, correlations, out_dir, args)
    write_threshold_exploration_plots(runs, out_dir)
    return summary, correlations


def write_threshold_exploration_report(runs, summary, correlations, out_dir, args):
    off_values = ", ".join(f"{case['off_pct']:g}" for case in args.off_cases)
    lines = [
        "# Threshold Exploration",
        "",
        f"- N: {args.N}",
        f"- L: {args.L}",
        f"- Samples per seed/profile: {args.samples}",
        f"- Generator seeds: {', '.join(str(seed) for seed in args.generator_seeds)}",
        f"- Data profiles: {', '.join(args.data_profiles)}",
        f"- Off percentages: {off_values}",
        f"- Sigma: {args.sigma}",
        f"- P: {args.P}",
        f"- Threshold target objective: {args.threshold_target}",
        "",
        "Primary metric is `U_G`. `u_g_vs_best_T_mean` is the mean fraction of the per-case best threshold objective.",
        "",
        "## Best Mean Thresholds",
        "",
        "| profile | off % | K | best mean T | source | mean U_G / best_T | winner rate | mean gap % |",
        "|---|---:|---:|---:|---|---:|---:|---:|",
    ]

    top_rows = (
        summary.sort_values(
            ["data_profile", "off_pct", "K", "u_g_vs_best_T_mean", "winner_rate"],
            ascending=[True, True, True, False, False],
        )
        .groupby(["data_profile", "off_pct", "K"], as_index=False)
        .head(1)
    )
    for _, row in top_rows.iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["data_profile"]),
                    format_float(row["off_pct"]),
                    str(int(row["K"])),
                    str(int(row["T"])),
                    str(row["threshold_source"]),
                    format_float(row["u_g_vs_best_T_mean"], precision=5),
                    format_float(row["winner_rate"], precision=3),
                    format_float(row["u_g_gap_to_best_pct_mean"], precision=3),
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Distribution-Metric Signal",
            "",
            "The table lists the strongest absolute Pearson correlation between a recorded distribution metric and the per-case best `T`.",
            "",
            "| profile | off % | K | strongest metric | corr | interpretation |",
            "|---|---:|---:|---|---:|---|",
        ]
    )

    if correlations.empty:
        lines.append("| n/a |  |  | n/a |  | no correlation data |")
    else:
        corr_rows = correlations.dropna(subset=["pearson_corr_with_best_T"]).copy()
        if corr_rows.empty:
            lines.append("| n/a |  |  | n/a |  | best T was constant or sample count was too small |")
        else:
            corr_rows["abs_corr"] = corr_rows["pearson_corr_with_best_T"].abs()
            best_corr = (
                corr_rows.sort_values(
                    ["data_profile", "off_pct", "K", "abs_corr"],
                    ascending=[True, True, True, False],
                )
                .groupby(["data_profile", "off_pct", "K"], as_index=False)
                .head(1)
            )
            for _, row in best_corr.iterrows():
                corr = row["pearson_corr_with_best_T"]
                strength = (
                    "strong candidate signal"
                    if abs(corr) >= 0.6
                    else "weak/no single-metric signal"
                )
                lines.append(
                    "| "
                    + " | ".join(
                        [
                            str(row["data_profile"]),
                            format_float(row["off_pct"]),
                            str(int(row["K"])),
                            str(row["metric"]),
                            format_float(corr, precision=3),
                            strength,
                        ]
                    )
                    + " |"
                )

    lines.extend(
        [
            "",
            "## Artifacts",
            "",
            "- `threshold_runs.csv`: one row per `(profile, seed, sample, K, T)`.",
            "- `threshold_summary.csv`: threshold rankings by mean `U_G` and normalized gap.",
            "- `threshold_metric_correlations.csv`: simple distribution-metric correlations with best `T`.",
            "- `threshold_cdf_u_g.png`: CDF of raw `U_G` for top thresholds.",
            "- `threshold_cdf_u_g_db.png`: CDF of `10 lg(U_G)` for top thresholds.",
            "- `threshold_cdf_u_g_vs_best.png`: CDF of `U_G(T) / max_T U_G(T)`.",
        ]
    )

    (out_dir / "threshold_exploration_report.md").write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def _join_threshold_sources(values):
    parts = []
    for value in values:
        for part in str(value).split("+"):
            if part and part not in parts:
                parts.append(part)
    return "+".join(parts)


def _mode_value(values):
    modes = values.mode()
    if modes.empty:
        return ""
    return str(modes.iloc[0])


def run_threshold_full_sweep(args):
    all_runs = []
    per_profile_outputs = []
    total_cases = (
        len(args.data_profiles)
        * len(args.generator_seeds)
        * args.samples
        * len(args.off_cases)
    )
    case_no = 0
    progress = (
        tqdm(total=total_cases, unit="case", dynamic_ncols=True)
        if tqdm is not None
        else None
    )

    try:
        for profile in args.data_profiles:
            profile_rows = []
            profile_dir = args.out_dir / profile
            for generator_seed in args.generator_seeds:
                rng = np.random.RandomState(generator_seed)
                for sample in range(args.samples):
                    V = generate_v_profile_from_rng(rng, args.N, args.L, profile)
                    distribution_metrics = row_power_distribution_metrics(V)

                    for off_case in args.off_cases:
                        case_no += 1
                        K = int(off_case["K"])
                        off_pct = float(off_case["off_pct"])
                        thresholds = dense_thresholds(K)
                        if progress is not None:
                            progress.set_postfix_str(
                                f"profile={profile}, seed={generator_seed}, "
                                f"sample={sample}, K={K}, thresholds=0..{K}"
                            )
                        else:
                            print(
                                f"[{case_no}/{total_cases}] profile={profile}, "
                                f"seed={generator_seed}, sample={sample}, K={K}, "
                                f"thresholds=0..{K}",
                                flush=True,
                            )

                        threshold_rows = evaluate_power_window_thresholds(
                            V,
                            K,
                            thresholds,
                            sigma=args.sigma,
                            P=args.P,
                        )
                        for row in threshold_rows:
                            profile_rows.append(
                                {
                                    "data_profile": profile,
                                    "generator_seed": int(generator_seed),
                                    "sample": int(sample),
                                    "N": int(args.N),
                                    "L": int(args.L),
                                    "K": K,
                                    "off_pct": off_pct,
                                    "active_pct": float(100.0 - off_pct),
                                    "sigma": float(args.sigma),
                                    "P": float(args.P),
                                    "target_obj": "gen",
                                    "threshold_source": "dense_0_to_K",
                                    **row,
                                    **distribution_metrics,
                                }
                            )

                        if progress is not None:
                            progress.update(1)

            profile_runs = add_full_sweep_relative_metrics(pd.DataFrame(profile_rows))
            atomic_write_csv(profile_runs, profile_dir / "threshold_runs.csv")
            outputs = write_full_sweep_profile_outputs(profile_runs, profile_dir, args)
            per_profile_outputs.append(outputs)
            all_runs.append(profile_runs)
    finally:
        if progress is not None:
            progress.close()

    combined_runs = pd.concat(all_runs, ignore_index=True) if all_runs else pd.DataFrame()
    return write_full_sweep_combined_outputs(combined_runs, per_profile_outputs, args)


def add_full_sweep_relative_metrics(runs):
    if runs.empty:
        return runs
    runs = runs.copy()
    runs["best_tested_u_g"] = runs.groupby(FULL_SWEEP_CASE_COLS)["u_g"].transform(
        "max"
    )
    runs["fraction_best_tested_u_g"] = runs["u_g"] / np.maximum(
        runs["best_tested_u_g"], np.finfo(float).eps
    )
    runs["gap_to_best_tested_pct"] = 100.0 * (1.0 - runs["fraction_best_tested_u_g"])
    runs["is_best_tested_T"] = np.isclose(runs["u_g"], runs["best_tested_u_g"])
    return runs


def write_full_sweep_profile_outputs(runs, profile_dir, args, formula_specs=FORMULA_SPECS):
    profile_dir.mkdir(parents=True, exist_ok=True)
    best_thresholds = build_full_sweep_best_thresholds(runs)
    summary = build_full_sweep_threshold_summary(runs)
    stats = build_full_sweep_best_t_stats(best_thresholds)
    formula = build_full_sweep_formula_comparison(runs, formula_specs=formula_specs)
    distribution = build_full_sweep_distribution_comparison(summary, stats, formula)

    atomic_write_csv(best_thresholds, profile_dir / "best_thresholds.csv")
    atomic_write_csv(summary, profile_dir / "threshold_summary.csv")
    atomic_write_csv(stats, profile_dir / "threshold_best_t_stats.csv")
    atomic_write_csv(formula, profile_dir / "threshold_formula_comparison.csv")
    write_full_sweep_profile_plots(runs, best_thresholds, summary, formula, profile_dir)
    write_full_sweep_profile_report(
        runs,
        best_thresholds,
        summary,
        stats,
        formula,
        distribution,
        profile_dir,
        args,
    )
    return {
        "profile": runs["data_profile"].iloc[0] if not runs.empty else "",
        "best_thresholds": best_thresholds,
        "summary": summary,
        "stats": stats,
        "formula": formula,
        "distribution": distribution,
    }


def build_full_sweep_best_thresholds(runs):
    best = (
        runs.sort_values(
            [*FULL_SWEEP_CASE_COLS, "u_g", "T"],
            ascending=[True] * len(FULL_SWEEP_CASE_COLS) + [False, True],
        )
        .drop_duplicates(FULL_SWEEP_CASE_COLS)
        .copy()
    )
    best = best.rename(
        columns={
            "T": "best_tested_T",
            "u_bf": "best_tested_u_bf",
            "u_i": "best_tested_u_i",
            "u_g": "best_tested_u_g_value",
            "u_g_db": "best_tested_u_g_db",
        }
    )

    lookup = runs.set_index([*FULL_SWEEP_CASE_COLS, "T"])["u_g"]
    prev_values = []
    next_values = []
    for _, row in best.iterrows():
        case_key = tuple(row[col] for col in FULL_SWEEP_CASE_COLS)
        T = int(row["best_tested_T"])
        prev_values.append(lookup.get((*case_key, T - 1), np.nan))
        next_values.append(lookup.get((*case_key, T + 1), np.nan))
    best["prev_T"] = best["best_tested_T"] - 1
    best["prev_T_u_g"] = prev_values
    best["next_T"] = best["best_tested_T"] + 1
    best["next_T_u_g"] = next_values
    best["best_tested_T_over_N"] = best["best_tested_T"] / best["N"].astype(float)
    best["best_tested_T_over_K"] = best["best_tested_T"] / np.maximum(
        best["K"].astype(float),
        np.finfo(float).eps,
    )

    keep_cols = [
        *FULL_SWEEP_CASE_COLS,
        "off_pct",
        "active_pct",
        "best_tested_T",
        "best_tested_T_over_N",
        "best_tested_T_over_K",
        "best_tested_u_bf",
        "best_tested_u_i",
        "best_tested_u_g_value",
        "best_tested_u_g_db",
        "prev_T",
        "prev_T_u_g",
        "next_T",
        "next_T_u_g",
        *THRESHOLD_METRIC_COLS,
    ]
    return best[keep_cols].sort_values(["data_profile", "K", "generator_seed", "sample"])


FULL_SWEEP_CONTEXT_COLS = [
    "data_profile",
    "N",
    "L",
    "K",
    "off_pct",
    "active_pct",
]


def build_full_sweep_threshold_summary(runs):
    def q(value):
        return lambda data: data.quantile(value)

    return (
        runs.groupby([*FULL_SWEEP_CONTEXT_COLS, "T"], as_index=False)
        .agg(
            samples=("u_g", "count"),
            winner_rate=("is_best_tested_T", "mean"),
            u_bf_mean=("u_bf", "mean"),
            u_i_mean=("u_i", "mean"),
            u_g_mean=("u_g", "mean"),
            u_g_p05=("u_g", q(0.05)),
            u_g_p50=("u_g", q(0.50)),
            u_g_p95=("u_g", q(0.95)),
            u_g_db_mean=("u_g_db", "mean"),
            fraction_best_tested_u_g_mean=("fraction_best_tested_u_g", "mean"),
            fraction_best_tested_u_g_p05=("fraction_best_tested_u_g", q(0.05)),
            gap_to_best_tested_pct_mean=("gap_to_best_tested_pct", "mean"),
            elapsed_seconds_mean=("elapsed_seconds", "mean"),
        )
        .sort_values(["data_profile", "N", "L", "K", "T"])
    )


def build_full_sweep_best_t_stats(best_thresholds):
    rows = []
    for key, chunk in best_thresholds.groupby(FULL_SWEEP_CONTEXT_COLS, sort=True):
        profile, N, L, K, off_pct, active_pct = key
        values = chunk["best_tested_T"].astype(float)
        mean = float(values.mean())
        std = float(values.std(ddof=1)) if len(values) > 1 else 0.0
        rows.append(
            {
                "data_profile": profile,
                "N": int(N),
                "L": int(L),
                "K": int(K),
                "off_pct": float(off_pct),
                "active_pct": float(active_pct),
                "cases": int(len(values)),
                "best_T_min": int(values.min()),
                "best_T_p05": float(values.quantile(0.05)),
                "best_T_p25": float(values.quantile(0.25)),
                "best_T_median": float(values.median()),
                "best_T_p75": float(values.quantile(0.75)),
                "best_T_p95": float(values.quantile(0.95)),
                "best_T_max": int(values.max()),
                "best_T_mean": mean,
                "best_T_std": std,
                "best_T_iqr": float(values.quantile(0.75) - values.quantile(0.25)),
                "best_T_cv": std / max(abs(mean), np.finfo(float).eps),
                "best_T_over_N_mean": float((values / float(N)).mean()),
                "best_T_over_K_mean": float((values / float(K)).mean()) if K else 0.0,
            }
        )
    return pd.DataFrame(rows).sort_values(["data_profile", "N", "L", "K"])


def build_full_sweep_formula_comparison(runs, formula_specs=FORMULA_SPECS):
    rows = []
    for key, chunk in runs.groupby(FULL_SWEEP_CONTEXT_COLS, sort=True):
        profile, N, L, K, off_pct, active_pct = key
        N = int(N)
        L = int(L)
        K = int(K)
        for label, _ in formula_specs:
            T = threshold_formula_T(label, N, L, K, formula_specs=formula_specs)
            selected = chunk[chunk["T"] == T]
            if selected.empty:
                continue
            rows.append(
                {
                    "data_profile": profile,
                    "N": N,
                    "L": L,
                    "K": int(K),
                    "off_pct": float(off_pct),
                    "active_pct": float(active_pct),
                    "formula": label,
                    "formula_family": formula_family(label),
                    "T": T,
                    "cases": int(selected["u_g"].count()),
                    "u_g_mean": float(selected["u_g"].mean()),
                    "fraction_best_tested_u_g_mean": float(
                        selected["fraction_best_tested_u_g"].mean()
                    ),
                    "fraction_best_tested_u_g_p05": float(
                        selected["fraction_best_tested_u_g"].quantile(0.05)
                    ),
                    "gap_to_best_tested_pct_mean": float(
                        selected["gap_to_best_tested_pct"].mean()
                    ),
                    "winner_rate": float(selected["is_best_tested_T"].mean()),
                }
            )
    return pd.DataFrame(rows).sort_values(["data_profile", "N", "L", "K", "formula"])


def build_full_sweep_distribution_comparison(summary, stats, formula):
    rows = []
    for key, chunk in summary.groupby(FULL_SWEEP_CONTEXT_COLS, sort=True):
        profile, N, L, K, off_pct, active_pct = key
        best_fixed = chunk.sort_values(["u_g_mean", "T"], ascending=[False, True]).iloc[0]
        band_low, band_high = full_sweep_recommended_band(chunk)
        stat_row = stats[
            (stats["data_profile"] == profile)
            & (stats["N"] == N)
            & (stats["L"] == L)
            & (stats["K"] == K)
        ].iloc[0]
        formula_chunk = formula[
            (formula["data_profile"] == profile)
            & (formula["N"] == N)
            & (formula["L"] == L)
            & (formula["K"] == K)
        ]
        best_formula = (
            formula_chunk.sort_values(
                ["fraction_best_tested_u_g_mean", "winner_rate"],
                ascending=[False, False],
            ).iloc[0]
            if not formula_chunk.empty
            else None
        )
        rows.append(
            {
                "data_profile": profile,
                "N": int(N),
                "L": int(L),
                "K": int(K),
                "off_pct": float(off_pct),
                "active_pct": float(active_pct),
                "best_fixed_T": int(best_fixed["T"]),
                "best_fixed_u_g_mean": float(best_fixed["u_g_mean"]),
                "recommended_band_low": int(band_low),
                "recommended_band_high": int(band_high),
                "recommended_band_width": int(band_high - band_low + 1),
                "best_tested_T_median": float(stat_row["best_T_median"]),
                "best_tested_T_p05": float(stat_row["best_T_p05"]),
                "best_tested_T_p95": float(stat_row["best_T_p95"]),
                "best_tested_T_std": float(stat_row["best_T_std"]),
                "best_tested_T_cv": float(stat_row["best_T_cv"]),
                "best_tested_T_over_N_mean": float(stat_row["best_T_over_N_mean"]),
                "best_tested_T_over_K_mean": float(stat_row["best_T_over_K_mean"]),
                "best_formula": "" if best_formula is None else str(best_formula["formula"]),
                "best_formula_T": np.nan if best_formula is None else int(best_formula["T"]),
                "best_formula_fraction": np.nan
                if best_formula is None
                else float(best_formula["fraction_best_tested_u_g_mean"]),
            }
        )
    return pd.DataFrame(rows).sort_values(["data_profile", "N", "L", "K"])


def full_sweep_recommended_band(chunk, threshold=0.99):
    ordered = chunk.sort_values("T").reset_index(drop=True)
    best_idx = int(ordered["u_g_mean"].idxmax())
    cutoff = float(ordered.loc[best_idx, "u_g_mean"]) * threshold
    mask = ordered["u_g_mean"] >= cutoff
    left = best_idx
    while left > 0 and bool(mask.iloc[left - 1]):
        left -= 1
    right = best_idx
    while right + 1 < len(ordered) and bool(mask.iloc[right + 1]):
        right += 1
    return int(ordered.loc[left, "T"]), int(ordered.loc[right, "T"])


def write_full_sweep_combined_outputs(combined_runs, per_profile_outputs, args):
    if combined_runs.empty:
        raise RuntimeError("No threshold full-sweep runs were produced.")

    all_best = pd.concat(
        [item["best_thresholds"] for item in per_profile_outputs],
        ignore_index=True,
    )
    all_stats = pd.concat([item["stats"] for item in per_profile_outputs], ignore_index=True)
    all_distribution = pd.concat(
        [item["distribution"] for item in per_profile_outputs],
        ignore_index=True,
    )
    all_formula = pd.concat(
        [item["formula"] for item in per_profile_outputs],
        ignore_index=True,
    )

    atomic_write_csv(all_best, args.out_dir / "all_best_thresholds.csv")
    atomic_write_csv(all_stats, args.out_dir / "all_threshold_best_t_stats.csv")
    atomic_write_csv(all_distribution, args.out_dir / "all_distribution_comparison.csv")
    atomic_write_csv(all_formula, args.out_dir / "all_formula_comparison.csv")
    write_full_sweep_combined_plots(all_best, all_distribution, all_formula, args.out_dir)
    write_full_sweep_combined_report(all_distribution, all_stats, all_formula, args.out_dir, args)
    return all_distribution


def write_full_sweep_profile_report(
    runs,
    best_thresholds,
    summary,
    stats,
    formula,
    distribution,
    profile_dir,
    args,
):
    profile = runs["data_profile"].iloc[0]
    lines = [
        f"# Threshold Full Sweep: {profile}",
        "",
        f"- N: {args.N}",
        f"- L: {args.L}",
        f"- K values: {', '.join(str(case['K']) for case in args.off_cases)}",
        f"- Samples: {args.samples}",
        f"- Generator seeds: {', '.join(str(seed) for seed in args.generator_seeds)}",
        f"- Sigma: {args.sigma}",
        "",
        "The experiment sweeps every integer `T` from `0` to `K` and evaluates raw `U_G`.",
        "",
        "## Answer",
        "",
    ]

    for _, row in distribution.iterrows():
        lines.append(
            f"- `K={int(row['K'])}`: best fixed `T={int(row['best_fixed_T'])}`; "
            f"99% mean-`U_G` diapason `{int(row['recommended_band_low'])}..{int(row['recommended_band_high'])}`; "
            f"best tested `T` median `{row['best_tested_T_median']:.1f}` "
            f"(p05..p95 `{row['best_tested_T_p05']:.1f}..{row['best_tested_T_p95']:.1f}`)."
        )

    lines.extend(
        [
            "",
            "## Best Fixed Thresholds And Formula Checks",
            "",
            "| K | best fixed T | 99% diapason | best tested T median | best tested T std | best formula | formula T | formula fraction |",
            "|---:|---:|---|---:|---:|---|---:|---:|",
        ]
    )
    for _, row in distribution.iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    str(int(row["K"])),
                    str(int(row["best_fixed_T"])),
                    f"{int(row['recommended_band_low'])}..{int(row['recommended_band_high'])}",
                    format_float(row["best_tested_T_median"]),
                    format_float(row["best_tested_T_std"]),
                    str(row["best_formula"]),
                    "" if pd.isna(row["best_formula_T"]) else str(int(row["best_formula_T"])),
                    format_float(row["best_formula_fraction"], precision=4),
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Plots",
            "",
            "![U_G vs T](threshold_u_g_by_T.png)",
            "",
            "![Best T histogram](threshold_best_T_hist.png)",
            "",
            "![Best T boxplot](threshold_best_T_boxplot.png)",
            "",
            "![Raw U_G CDF](threshold_raw_u_g_cdf.png)",
            "",
            "![Fraction of best tested U_G CDF](threshold_fraction_best_cdf.png)",
            "",
            "## Artifacts",
            "",
            f"- `{threshold_runs_artifact_name(profile_dir)}`",
            "- `best_thresholds.csv`",
            "- `threshold_summary.csv`",
            "- `threshold_best_t_stats.csv`",
            "- `threshold_formula_comparison.csv`",
        ]
    )
    (profile_dir / "threshold_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def threshold_runs_artifact_name(profile_dir):
    if (profile_dir / "threshold_runs.csv.gz").exists():
        return "threshold_runs.csv.gz"
    return "threshold_runs.csv"


def write_full_sweep_combined_report(all_distribution, all_stats, all_formula, out_dir, args):
    best_global_formula = (
        all_formula.groupby("formula", as_index=False)
        .agg(
            fraction_best_tested_u_g_mean=("fraction_best_tested_u_g_mean", "mean"),
            gap_to_best_tested_pct_mean=("gap_to_best_tested_pct_mean", "mean"),
        )
        .sort_values("fraction_best_tested_u_g_mean", ascending=False)
    )
    best_formula_row = best_global_formula.iloc[0] if not best_global_formula.empty else None

    lines = [
        "# Threshold Full Sweep: All Distributions",
        "",
        f"- N: {args.N}",
        f"- L: {args.L}",
        f"- K values: {', '.join(str(case['K']) for case in args.off_cases)}",
        f"- Samples: {args.samples}",
        f"- Generator seeds: {', '.join(str(seed) for seed in args.generator_seeds)}",
        f"- Profiles: {', '.join(args.data_profiles)}",
        f"- Sigma: {args.sigma}",
        "",
        "## Direct Answer",
        "",
    ]
    if best_formula_row is not None:
        lines.append(
            f"- Best simple formula overall: `{best_formula_row['formula']}` "
            f"with mean fraction of best tested `U_G = {best_formula_row['fraction_best_tested_u_g_mean']:.4f}`."
        )
    lines.extend(
        [
            "- A single global formula is acceptable only if its per-distribution rows stay close in the table below; otherwise use the reported 99% diapason per distribution.",
            "- Scaling is reported both as `T/N` and `T/K`; compare those columns to see which is more stable.",
            "",
            "## Distribution Comparison",
            "",
            "| profile | K | best fixed T | 99% diapason | best tested T median | T/N mean | T/K mean | best formula | formula fraction |",
            "|---|---:|---:|---|---:|---:|---:|---|---:|",
        ]
    )

    for _, row in all_distribution.iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["data_profile"]),
                    str(int(row["K"])),
                    str(int(row["best_fixed_T"])),
                    f"{int(row['recommended_band_low'])}..{int(row['recommended_band_high'])}",
                    format_float(row["best_tested_T_median"]),
                    format_float(row["best_tested_T_over_N_mean"], precision=4),
                    format_float(row["best_tested_T_over_K_mean"], precision=4),
                    str(row["best_formula"]),
                    format_float(row["best_formula_fraction"], precision=4),
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Combined Plots",
            "",
            "![Best tested T boxplot](all_best_T_boxplot.png)",
            "",
            "![Best fixed T by distribution](all_best_fixed_T.png)",
            "",
            "![Formula comparison](all_formula_fraction.png)",
            "",
            "## Artifacts",
            "",
            "- `all_best_thresholds.csv`",
            "- `all_threshold_best_t_stats.csv`",
            "- `all_distribution_comparison.csv`",
            "- `all_formula_comparison.csv`",
        ]
    )
    (out_dir / "all_experiments_report.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def write_full_sweep_profile_plots(runs, best_thresholds, summary, formula, out_dir):
    from utils.plotting import use_agg_backend

    use_agg_backend()
    import matplotlib.pyplot as plt

    Ks = sorted(summary["K"].unique())
    fig, axes = plt.subplots(1, len(Ks), figsize=(6.4 * len(Ks), 4.2), squeeze=False)
    for col, K in enumerate(Ks):
        ax = axes[0, col]
        chunk = summary[summary["K"] == K].sort_values("T")
        ax.plot(chunk["T"], chunk["u_g_mean"], color="#0072B2", linewidth=1.8, label="mean")
        ax.fill_between(
            chunk["T"].to_numpy(),
            chunk["u_g_p05"].to_numpy(),
            chunk["u_g_p95"].to_numpy(),
            color="#0072B2",
            alpha=0.18,
            label="p05..p95",
        )
        ax.set_title(f"K={int(K)}")
        ax.set_xlabel("T")
        ax.set_ylabel("U_G")
        ax.set_yscale("log")
        ax.grid(True, alpha=0.25)
        ax.legend()
    fig.suptitle("Mean raw U_G by threshold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out_dir / "threshold_u_g_by_T.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, len(Ks), figsize=(6.4 * len(Ks), 4.0), squeeze=False)
    for col, K in enumerate(Ks):
        ax = axes[0, col]
        values = best_thresholds[best_thresholds["K"] == K]["best_tested_T"]
        ax.hist(values, bins=30, color="#009E73", alpha=0.85)
        ax.set_title(f"K={int(K)}")
        ax.set_xlabel("best tested T")
        ax.set_ylabel("count")
        ax.grid(True, alpha=0.25)
    fig.suptitle("Best tested T histogram")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out_dir / "threshold_best_T_hist.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    box_data = [best_thresholds[best_thresholds["K"] == K]["best_tested_T"] for K in Ks]
    ax.boxplot(box_data, tick_labels=[f"K={int(K)}" for K in Ks], showmeans=True)
    ax.set_ylabel("best tested T")
    ax.set_title("Best tested T by K")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_dir / "threshold_best_T_boxplot.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, len(Ks), figsize=(6.4 * len(Ks), 4.2), squeeze=False)
    for col, K in enumerate(Ks):
        ax = axes[0, col]
        selected = selected_formula_thresholds(formula[formula["K"] == K])
        for label, T in selected:
            values = runs[(runs["K"] == K) & (runs["T"] == T)]["u_g"]
            x_values, y_values = empirical_cdf_local(values)
            if len(x_values):
                ax.step(y_values, x_values, where="post", label=f"{label} (T={T})")
        ax.set_title(f"K={int(K)}")
        ax.set_xlabel("cumulative fraction")
        ax.set_ylabel("U_G")
        ax.set_yscale("log")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
    fig.suptitle("Raw U_G CDF for formula thresholds")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out_dir / "threshold_raw_u_g_cdf.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, len(Ks), figsize=(6.4 * len(Ks), 4.2), squeeze=False)
    for col, K in enumerate(Ks):
        ax = axes[0, col]
        selected = selected_formula_thresholds(formula[formula["K"] == K])
        for label, T in selected:
            values = runs[(runs["K"] == K) & (runs["T"] == T)][
                "fraction_best_tested_u_g"
            ]
            x_values, y_values = empirical_cdf_local(values)
            if len(x_values):
                ax.step(y_values, x_values, where="post", label=f"{label} (T={T})")
        ax.set_title(f"K={int(K)}")
        ax.set_xlabel("cumulative fraction")
        ax.set_ylabel("fraction of best tested U_G")
        ax.set_ylim(0.0, 1.05)
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
    fig.suptitle("Fraction of best tested U_G CDF for formula thresholds")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out_dir / "threshold_fraction_best_cdf.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_full_sweep_combined_plots(all_best, all_distribution, all_formula, out_dir):
    from utils.plotting import use_agg_backend

    use_agg_backend()
    import matplotlib.pyplot as plt

    profiles = sorted(all_best["data_profile"].unique())
    Ks = sorted(all_best["K"].unique())

    fig, axes = plt.subplots(1, len(Ks), figsize=(7.0 * len(Ks), 4.6), squeeze=False)
    for col, K in enumerate(Ks):
        ax = axes[0, col]
        data = [
            all_best[(all_best["data_profile"] == profile) & (all_best["K"] == K)][
                "best_tested_T"
            ]
            for profile in profiles
        ]
        ax.boxplot(data, tick_labels=profiles, showmeans=True)
        ax.tick_params(axis="x", rotation=30)
        ax.set_title(f"K={int(K)}")
        ax.set_ylabel("best tested T")
        ax.grid(True, axis="y", alpha=0.25)
    fig.suptitle("Best tested T by distribution")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out_dir / "all_best_T_boxplot.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, len(Ks), figsize=(7.0 * len(Ks), 4.2), squeeze=False)
    for col, K in enumerate(Ks):
        ax = axes[0, col]
        chunk = all_distribution[all_distribution["K"] == K]
        ax.bar(chunk["data_profile"], chunk["best_fixed_T"], color="#D55E00")
        ax.tick_params(axis="x", rotation=30)
        ax.set_title(f"K={int(K)}")
        ax.set_ylabel("best fixed T")
        ax.grid(True, axis="y", alpha=0.25)
    fig.suptitle("Best fixed T by distribution")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out_dir / "all_best_fixed_T.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    formula_mean = (
        all_formula.groupby("formula", as_index=False)["fraction_best_tested_u_g_mean"]
        .mean()
        .sort_values("fraction_best_tested_u_g_mean", ascending=False)
    )
    fig, ax = plt.subplots(figsize=(8.4, 4.4))
    ax.bar(formula_mean["formula"], formula_mean["fraction_best_tested_u_g_mean"], color="#009E73")
    ax.tick_params(axis="x", rotation=35)
    ax.set_ylim(0.0, 1.02)
    ax.set_ylabel("mean fraction of best tested U_G")
    ax.set_title("Formula threshold comparison")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_dir / "all_formula_fraction.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def selected_formula_thresholds(formula):
    if formula.empty:
        return []
    selected = []
    best = formula.sort_values("fraction_best_tested_u_g_mean", ascending=False).iloc[0]
    selected.append(("best formula", int(best["T"])))
    for name in ("T_0p05N", "T_0p075N", "T_0p10N", "T_0p10K", "T_0p20K"):
        row = formula[formula["formula"] == name]
        if not row.empty:
            selected.append((name, int(row.iloc[0]["T"])))
    result = []
    seen = set()
    for label, T in selected:
        if T not in seen:
            seen.add(T)
            result.append((label, T))
    return result[:6]


def threshold_scaling_shard_args(args, N, L):
    values = vars(args).copy()
    values["N"] = int(N)
    values["L"] = int(L)
    values["off_cases"] = build_active_pct_cases(int(N), args.K_pcts)
    return SimpleNamespace(**values)


def collect_threshold_scaling_profile_runs(args, N, L, profile, progress=None):
    shard_args = threshold_scaling_shard_args(args, N, L)
    rows = []

    for generator_seed in shard_args.generator_seeds:
        rng = np.random.RandomState(generator_seed)
        for sample in range(shard_args.samples):
            V = generate_v_profile_from_rng(rng, shard_args.N, shard_args.L, profile)
            distribution_metrics = row_power_distribution_metrics(V)

            for off_case in shard_args.off_cases:
                K = int(off_case["K"])
                thresholds = dense_thresholds(K)
                threshold_rows = evaluate_power_window_thresholds(
                    V,
                    K,
                    thresholds,
                    sigma=shard_args.sigma,
                    P=shard_args.P,
                )
                for row in threshold_rows:
                    rows.append(
                        {
                            "data_profile": profile,
                            "generator_seed": int(generator_seed),
                            "sample": int(sample),
                            "N": int(shard_args.N),
                            "L": int(shard_args.L),
                            "K": K,
                            "off_pct": float(off_case["off_pct"]),
                            "active_pct": float(off_case["active_pct"]),
                            "sigma": float(shard_args.sigma),
                            "P": float(shard_args.P),
                            "target_obj": "gen",
                            "threshold_source": "dense_0_to_K",
                            **row,
                            **distribution_metrics,
                        }
                    )

                if progress is not None:
                    progress.set_postfix_str(
                        f"N={shard_args.N}, L={shard_args.L}, profile={profile}, K={K}"
                    )
                    progress.update(1)

    return add_full_sweep_relative_metrics(pd.DataFrame(rows))


def run_threshold_scaling_study(args):
    total_cases = (
        len(args.N_values)
        * len(args.L_values)
        * len(args.data_profiles)
        * len(args.generator_seeds)
        * args.samples
        * len(args.K_pcts)
    )
    progress = (
        tqdm(total=total_cases, unit="case", dynamic_ncols=True)
        if tqdm is not None
        else None
    )

    try:
        for N in args.N_values:
            for L in args.L_values:
                shard_args = threshold_scaling_shard_args(args, N, L)
                for profile in args.data_profiles:
                    profile_dir = args.out_dir / f"N{int(N)}" / f"L{int(L)}" / profile
                    print(
                        f"Running scaling shard N={int(N)}, L={int(L)}, profile={profile}",
                        flush=True,
                    )
                    runs = collect_threshold_scaling_profile_runs(
                        args,
                        int(N),
                        int(L),
                        profile,
                        progress=progress,
                    )
                    atomic_write_csv(runs, profile_dir / "threshold_runs.csv.gz")
                    outputs = write_full_sweep_profile_outputs(
                        runs,
                        profile_dir,
                        shard_args,
                        formula_specs=SCALING_FORMULA_SPECS,
                    )
                    selected = build_scaling_formula_selected_runs(
                        runs,
                        formula_specs=SCALING_FORMULA_SPECS,
                    )
                    atomic_write_csv(selected, profile_dir / "threshold_formula_selected_runs.csv")
                    strong_weak = build_scaling_strong_weak_runs(outputs["best_thresholds"])
                    atomic_write_csv(strong_weak, profile_dir / "threshold_strong_weak_runs.csv")
    finally:
        if progress is not None:
            progress.close()

    per_shard_outputs, formula_selected_parts = load_threshold_scaling_shard_outputs(
        args.out_dir
    )
    return write_threshold_scaling_combined_outputs(
        per_shard_outputs,
        formula_selected_parts,
        args,
    )


def load_threshold_scaling_shard_outputs(out_dir):
    per_shard_outputs = []
    formula_selected_parts = []
    for shard_dir in sorted(out_dir.glob("N*/L*/*")):
        if not shard_dir.is_dir():
            continue
        required = {
            "best_thresholds": shard_dir / "best_thresholds.csv",
            "summary": shard_dir / "threshold_summary.csv",
            "stats": shard_dir / "threshold_best_t_stats.csv",
            "formula": shard_dir / "threshold_formula_comparison.csv",
        }
        if not all(path.exists() for path in required.values()):
            continue

        best_thresholds = pd.read_csv(required["best_thresholds"])
        summary = pd.read_csv(required["summary"])
        stats = pd.read_csv(required["stats"])
        formula = pd.read_csv(required["formula"])
        strong_weak_path = shard_dir / "threshold_strong_weak_runs.csv"
        if strong_weak_path.exists():
            strong_weak = pd.read_csv(strong_weak_path)
        else:
            strong_weak = build_scaling_strong_weak_runs(best_thresholds)
            atomic_write_csv(strong_weak, strong_weak_path)
        strong_weak_formula = build_scaling_strong_weak_formula_comparison(strong_weak)
        if not strong_weak_formula.empty:
            formula = pd.concat([formula, strong_weak_formula], ignore_index=True)
        distribution = build_full_sweep_distribution_comparison(
            summary,
            stats,
            formula,
        )
        per_shard_outputs.append(
            {
                "profile": str(summary["data_profile"].iloc[0])
                if not summary.empty
                else shard_dir.name,
                "best_thresholds": best_thresholds,
                "summary": summary,
                "stats": stats,
                "formula": formula,
                "distribution": distribution,
            }
        )

        selected_path = shard_dir / "threshold_formula_selected_runs.csv"
        if selected_path.exists():
            selected = pd.read_csv(selected_path)
            if not strong_weak.empty:
                selected = pd.concat([selected, strong_weak], ignore_index=True, sort=False)
            formula_selected_parts.append(selected)

    return per_shard_outputs, formula_selected_parts


def build_scaling_formula_selected_runs(runs, formula_specs=SCALING_FORMULA_SPECS):
    keep_cols = [
        *FULL_SWEEP_CASE_COLS,
        "off_pct",
        "active_pct",
        "T",
        "u_bf",
        "u_i",
        "u_g",
        "u_g_db",
        "best_tested_u_g",
        "fraction_best_tested_u_g",
        "gap_to_best_tested_pct",
        "is_best_tested_T",
    ]
    parts = []
    best = (
        runs[boolean_series(runs["is_best_tested_T"])]
        .sort_values([*FULL_SWEEP_CASE_COLS, "T"])
        .drop_duplicates(FULL_SWEEP_CASE_COLS, keep="first")
        .copy()
    )
    if not best.empty:
        best = best[keep_cols].copy()
        best["formula"] = "best_tested_T"
        best["formula_family"] = "best_tested"
        best["formula_T"] = best["T"].astype(int)
        best["formula_label"] = "best tested T per sample"
        best["fraction_best_tested_u_g"] = 1.0
        best["gap_to_best_tested_pct"] = 0.0
        best["is_best_tested_T"] = True
        parts.append(best)

    for key, chunk in runs.groupby(FULL_SWEEP_CONTEXT_COLS, sort=True):
        profile, N, L, K, _off_pct, _active_pct = key
        for label, _ in formula_specs:
            T = threshold_formula_T(label, int(N), int(L), int(K), formula_specs=formula_specs)
            selected = chunk[chunk["T"] == T]
            if selected.empty:
                continue
            selected = selected[keep_cols].copy()
            selected["formula"] = label
            selected["formula_family"] = formula_family(label)
            selected["formula_T"] = int(T)
            selected["formula_label"] = label
            parts.append(selected)

    if not parts:
        return pd.DataFrame(columns=[*keep_cols, "formula", "formula_family", "formula_T", "formula_label"])
    return pd.concat(parts, ignore_index=True).sort_values(
        ["data_profile", "N", "L", "K", "formula", "generator_seed", "sample"]
    )


def build_scaling_strong_weak_runs(best_thresholds):
    if best_thresholds.empty:
        return pd.DataFrame()

    profile = str(best_thresholds["data_profile"].iloc[0])
    N = int(best_thresholds["N"].iloc[0])
    L = int(best_thresholds["L"].iloc[0])
    rows = []

    for generator_seed in sorted(best_thresholds["generator_seed"].unique()):
        seed_cases = best_thresholds[best_thresholds["generator_seed"] == generator_seed]
        sample_max = int(seed_cases["sample"].max())
        rng = np.random.RandomState(int(generator_seed))
        for sample in range(sample_max + 1):
            V = generate_v_profile_from_rng(rng, N, L, profile)
            sample_cases = seed_cases[seed_cases["sample"] == sample]
            if sample_cases.empty:
                continue
            for _, best_row in sample_cases.iterrows():
                K = int(best_row["K"])
                sigma = float(best_row["sigma"])
                P = float(best_row["P"])
                u_bf, u_i, u_g, T = evaluate_strong_weak_window(
                    V,
                    K,
                    sigma=sigma,
                    P=P,
                )
                best_u_g = float(best_row["best_tested_u_g_value"])
                fraction = u_g / best_u_g if best_u_g > 0 else np.nan
                rows.append(
                    {
                        "data_profile": profile,
                        "generator_seed": int(generator_seed),
                        "sample": int(sample),
                        "N": N,
                        "L": L,
                        "K": K,
                        "sigma": sigma,
                        "P": P,
                        "off_pct": float(best_row["off_pct"]),
                        "active_pct": float(best_row["active_pct"]),
                        "T": int(T),
                        "u_bf": float(u_bf),
                        "u_i": float(u_i),
                        "u_g": float(u_g),
                        "u_g_db": 10.0 * np.log10(max(float(u_g), np.finfo(float).tiny)),
                        "best_tested_u_g": best_u_g,
                        "fraction_best_tested_u_g": fraction,
                        "gap_to_best_tested_pct": 100.0 * (1.0 - fraction),
                        "is_best_tested_T": bool(fraction >= 1.0 - 1e-12),
                        "formula": STRONG_WEAK_RULE,
                        "formula_family": formula_family(STRONG_WEAK_RULE),
                        "formula_T": int(T),
                        "formula_label": "strong/weak H3",
                        "rule_note": "outside_dense_0_to_K" if T > K else "",
                    }
                )

    return pd.DataFrame(rows).sort_values(
        ["data_profile", "N", "L", "K", "generator_seed", "sample"]
    )


def build_scaling_strong_weak_formula_comparison(strong_weak_runs):
    if strong_weak_runs.empty:
        return pd.DataFrame()

    rows = []
    for key, chunk in strong_weak_runs.groupby(FULL_SWEEP_CONTEXT_COLS, sort=True):
        profile, N, L, K, off_pct, active_pct = key
        rows.append(
            {
                "data_profile": profile,
                "N": int(N),
                "L": int(L),
                "K": int(K),
                "off_pct": float(off_pct),
                "active_pct": float(active_pct),
                "formula": STRONG_WEAK_RULE,
                "formula_family": formula_family(STRONG_WEAK_RULE),
                "T": int(round(float(chunk["formula_T"].median()))),
                "cases": int(chunk["u_g"].count()),
                "u_g_mean": float(chunk["u_g"].mean()),
                "fraction_best_tested_u_g_mean": float(
                    chunk["fraction_best_tested_u_g"].mean()
                ),
                "fraction_best_tested_u_g_p05": float(
                    chunk["fraction_best_tested_u_g"].quantile(0.05)
                ),
                "gap_to_best_tested_pct_mean": float(
                    chunk["gap_to_best_tested_pct"].mean()
                ),
                "winner_rate": float(chunk["is_best_tested_T"].mean()),
            }
        )
    return pd.DataFrame(rows).sort_values(["data_profile", "N", "L", "K"])


def write_threshold_scaling_combined_outputs(per_shard_outputs, formula_selected_parts, args):
    if not per_shard_outputs:
        raise RuntimeError("No threshold scaling shards were produced.")

    all_best = pd.concat(
        [item["best_thresholds"] for item in per_shard_outputs],
        ignore_index=True,
    )
    all_stats = pd.concat([item["stats"] for item in per_shard_outputs], ignore_index=True)
    all_distribution = pd.concat(
        [item["distribution"] for item in per_shard_outputs],
        ignore_index=True,
    )
    all_formula = pd.concat(
        [item["formula"] for item in per_shard_outputs],
        ignore_index=True,
    )
    formula_selected = pd.concat(formula_selected_parts, ignore_index=True)
    strong_weak_runs = formula_selected[
        formula_selected["formula"] == STRONG_WEAK_RULE
    ].copy()
    strong_weak_summary = build_scaling_strong_weak_summary(strong_weak_runs)
    scaling_formula_summary = build_scaling_formula_summary(all_formula)
    metric_correlations = build_scaling_metric_correlations(all_best)

    atomic_write_csv(all_best, args.out_dir / "all_best_thresholds.csv")
    atomic_write_csv(all_stats, args.out_dir / "all_threshold_best_t_stats.csv")
    atomic_write_csv(all_distribution, args.out_dir / "all_distribution_comparison.csv")
    atomic_write_csv(all_formula, args.out_dir / "all_formula_comparison.csv")
    atomic_write_csv(scaling_formula_summary, args.out_dir / "all_scaling_formula_summary.csv")
    atomic_write_csv(metric_correlations, args.out_dir / "all_scaling_metric_correlations.csv")
    atomic_write_csv(formula_selected, args.out_dir / "all_formula_selected_runs.csv")
    atomic_write_csv(strong_weak_runs, args.out_dir / "all_strong_weak_runs.csv")
    atomic_write_csv(strong_weak_summary, args.out_dir / "all_strong_weak_summary.csv")

    write_threshold_scaling_combined_plots(
        all_best,
        all_formula,
        formula_selected,
        args.out_dir,
    )
    write_threshold_scaling_report(
        all_distribution,
        all_stats,
        all_formula,
        scaling_formula_summary,
        metric_correlations,
        strong_weak_summary,
        args.out_dir,
        args,
    )
    return all_distribution


def build_scaling_strong_weak_summary(strong_weak_runs):
    if strong_weak_runs.empty:
        return pd.DataFrame()

    def q(value):
        return lambda data: data.quantile(value)

    return (
        strong_weak_runs.groupby(FULL_SWEEP_CONTEXT_COLS, as_index=False)
        .agg(
            cases=("u_g", "count"),
            T_min=("formula_T", "min"),
            T_max=("formula_T", "max"),
            T_mean=("formula_T", "mean"),
            outside_dense_rate=(
                "rule_note",
                lambda values: (values.astype(str) == "outside_dense_0_to_K").mean(),
            ),
            u_g_mean=("u_g", "mean"),
            u_g_p05=("u_g", q(0.05)),
            u_g_p50=("u_g", q(0.50)),
            u_g_p95=("u_g", q(0.95)),
            fraction_mean=("fraction_best_tested_u_g", "mean"),
            fraction_p05=("fraction_best_tested_u_g", q(0.05)),
            fraction_p50=("fraction_best_tested_u_g", q(0.50)),
            fraction_p95=("fraction_best_tested_u_g", q(0.95)),
            gap_pct_mean=("gap_to_best_tested_pct", "mean"),
            winner_rate=("is_best_tested_T", "mean"),
        )
        .sort_values(["data_profile", "N", "L", "K"])
    )


def build_scaling_formula_summary(all_formula):
    frames = []
    group_specs = [
        ("overall", ["formula", "formula_family"]),
        ("by_N", ["N", "formula", "formula_family"]),
        ("by_L", ["L", "formula", "formula_family"]),
        ("by_active_pct", ["active_pct", "formula", "formula_family"]),
        ("by_profile", ["data_profile", "formula", "formula_family"]),
        ("by_N_L", ["N", "L", "formula", "formula_family"]),
        ("by_N_L_active_pct", ["N", "L", "active_pct", "formula", "formula_family"]),
        ("by_profile_L", ["data_profile", "L", "formula", "formula_family"]),
    ]
    for scope, columns in group_specs:
        frame = (
            all_formula.groupby(columns, as_index=False)
            .agg(
                contexts=("fraction_best_tested_u_g_mean", "count"),
                fraction_mean=("fraction_best_tested_u_g_mean", "mean"),
                fraction_p05_mean=("fraction_best_tested_u_g_p05", "mean"),
                gap_pct_mean=("gap_to_best_tested_pct_mean", "mean"),
                winner_rate_mean=("winner_rate", "mean"),
            )
            .sort_values("fraction_mean", ascending=False)
        )
        frame.insert(0, "scope", scope)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True, sort=False)


def build_scaling_metric_correlations(all_best):
    targets = ["best_tested_T_over_N", "best_tested_T_over_K"]
    rows = []
    group_specs = [
        ("overall", []),
        ("by_profile", ["data_profile"]),
        ("by_L", ["L"]),
        ("by_active_pct", ["active_pct"]),
    ]
    for scope, group_cols in group_specs:
        grouped = [((), all_best)] if not group_cols else all_best.groupby(group_cols, sort=True)
        for key, chunk in grouped:
            if not isinstance(key, tuple):
                key = (key,)
            context = dict(zip(group_cols, key))
            for target in targets:
                target_values = chunk[target].astype(float)
                for metric in THRESHOLD_METRIC_COLS:
                    metric_values = chunk[metric].astype(float)
                    valid = np.isfinite(target_values) & np.isfinite(metric_values)
                    if int(valid.sum()) < 3:
                        corr = np.nan
                    elif float(metric_values[valid].std()) == 0.0:
                        corr = np.nan
                    else:
                        corr = float(np.corrcoef(target_values[valid], metric_values[valid])[0, 1])
                    rows.append(
                        {
                            "scope": scope,
                            **context,
                            "target": target,
                            "metric": metric,
                            "pearson_corr": corr,
                            "abs_pearson_corr": np.nan if pd.isna(corr) else abs(corr),
                            "cases": int(valid.sum()),
                        }
                    )
    return pd.DataFrame(rows).sort_values(
        ["scope", "target", "abs_pearson_corr"],
        ascending=[True, True, False],
    )


def top_scaling_formula_labels(all_formula, max_count=6):
    ranked = (
        all_formula.groupby("formula", as_index=False)["fraction_best_tested_u_g_mean"]
        .mean()
        .sort_values("fraction_best_tested_u_g_mean", ascending=False)
    )
    labels = ["best_tested_T"]
    for label in ranked["formula"]:
        if label not in labels:
            labels.append(label)
        if len(labels) >= max_count:
            break
    if "T_0p05N" not in labels:
        labels.append("T_0p05N")
    if STRONG_WEAK_RULE not in labels:
        labels.append(STRONG_WEAK_RULE)
    return labels


def write_threshold_scaling_combined_plots(all_best, all_formula, formula_selected, out_dir):
    from utils.plotting import use_agg_backend

    use_agg_backend()
    import matplotlib.pyplot as plt

    labels = top_scaling_formula_labels(all_formula)
    selected = formula_selected[formula_selected["formula"].isin(labels)]

    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    for label in labels:
        chunk = selected[selected["formula"] == label]
        if chunk.empty:
            continue
        values, probs = empirical_cdf_local(chunk["u_g"])
        if len(values):
            ax.step(probs, values, where="post", linewidth=1.5, label=label)
    ax.set_xlabel("cumulative fraction")
    ax.set_ylabel("raw U_G")
    ax.set_yscale("log")
    ax.set_xlim(0.0, 1.0)
    ax.set_title("Scaling study CDF: raw U_G by threshold rule")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "scaling_cdf_raw_u_g.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    for label in labels:
        chunk = selected[selected["formula"] == label]
        if chunk.empty:
            continue
        values, probs = empirical_cdf_local(chunk["fraction_best_tested_u_g"])
        if len(values):
            ax.step(probs, values, where="post", linewidth=1.5, label=label)
    ax.set_xlabel("cumulative fraction")
    ax.set_ylabel("fraction of best tested U_G")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.05)
    ax.set_title("Scaling study CDF: fraction of best tested U_G")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "scaling_cdf_fraction_best_u_g.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    write_threshold_scaling_ratio_plot(
        all_best,
        out_dir / "scaling_best_T_over_N.png",
        "best_tested_T_over_N",
        "best tested T / N",
    )
    write_threshold_scaling_ratio_plot(
        all_best,
        out_dir / "scaling_best_T_over_K.png",
        "best_tested_T_over_K",
        "best tested T / K",
    )
    write_threshold_scaling_formula_plot(
        all_formula,
        labels,
        out_dir / "scaling_formula_fraction_by_N_L.png",
    )
    write_threshold_scaling_profile_formula_plot(
        all_formula,
        labels,
        out_dir / "scaling_formula_fraction_by_profile.png",
    )


def write_threshold_scaling_ratio_plot(all_best, out_path, value_col, ylabel):
    import matplotlib.pyplot as plt

    grouped = (
        all_best.groupby(["N", "L", "data_profile", "active_pct"], as_index=False)[value_col]
        .median()
        .sort_values(["N", "active_pct", "data_profile", "L"])
    )
    Ns = sorted(grouped["N"].unique())
    active_pcts = sorted(grouped["active_pct"].unique())
    profiles = sorted(grouped["data_profile"].unique())
    fig, axes = plt.subplots(
        len(Ns),
        len(active_pcts),
        figsize=(5.4 * len(active_pcts), 2.8 * len(Ns)),
        squeeze=False,
        sharey=True,
    )
    for row_index, N in enumerate(Ns):
        for col_index, active_pct in enumerate(active_pcts):
            ax = axes[row_index, col_index]
            chunk = grouped[(grouped["N"] == N) & (grouped["active_pct"] == active_pct)]
            for profile in profiles:
                profile_chunk = chunk[chunk["data_profile"] == profile]
                if profile_chunk.empty:
                    continue
                ax.plot(
                    profile_chunk["L"],
                    profile_chunk[value_col],
                    marker="o",
                    linewidth=1.2,
                    label=profile,
                )
            ax.set_title(f"N={int(N)}, K={active_pct:g}% active")
            ax.set_xlabel("L")
            ax.set_ylabel(ylabel)
            ax.grid(True, alpha=0.25)
            if row_index == 0 and col_index == len(active_pcts) - 1:
                ax.legend(fontsize=7, bbox_to_anchor=(1.02, 1.0), loc="upper left")
    fig.suptitle(ylabel)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_threshold_scaling_formula_plot(all_formula, labels, out_path):
    import matplotlib.pyplot as plt

    labels = [label for label in labels if label != "best_tested_T"]
    selected = all_formula[all_formula["formula"].isin(labels)]
    active_pcts = sorted(selected["active_pct"].unique())
    Ns = sorted(selected["N"].unique())
    fig, axes = plt.subplots(
        len(Ns),
        len(active_pcts),
        figsize=(5.4 * len(active_pcts), 2.8 * len(Ns)),
        squeeze=False,
        sharey=True,
    )
    for row_index, N in enumerate(Ns):
        for col_index, active_pct in enumerate(active_pcts):
            ax = axes[row_index, col_index]
            chunk = selected[(selected["N"] == N) & (selected["active_pct"] == active_pct)]
            grouped = (
                chunk.groupby(["L", "formula"], as_index=False)[
                    "fraction_best_tested_u_g_mean"
                ]
                .mean()
                .sort_values(["formula", "L"])
            )
            for label in labels:
                formula_chunk = grouped[grouped["formula"] == label]
                if formula_chunk.empty:
                    continue
                ax.plot(
                    formula_chunk["L"],
                    formula_chunk["fraction_best_tested_u_g_mean"],
                    marker="o",
                    linewidth=1.2,
                    label=label,
                )
            ax.set_title(f"N={int(N)}, K={active_pct:g}% active")
            ax.set_xlabel("L")
            ax.set_ylabel("mean fraction of best tested U_G")
            ax.set_ylim(0.0, 1.02)
            ax.grid(True, alpha=0.25)
            if row_index == 0 and col_index == len(active_pcts) - 1:
                ax.legend(fontsize=7, bbox_to_anchor=(1.02, 1.0), loc="upper left")
    fig.suptitle("Formula performance by N and L")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_threshold_scaling_profile_formula_plot(all_formula, labels, out_path):
    import matplotlib.pyplot as plt

    labels = [label for label in labels if label != "best_tested_T"]
    selected = all_formula[all_formula["formula"].isin(labels)]
    profiles = sorted(selected["data_profile"].unique())
    columns = 2
    rows = int(np.ceil(len(profiles) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(7.2 * columns, 3.8 * rows), squeeze=False)
    for index, profile in enumerate(profiles):
        ax = axes[index // columns, index % columns]
        chunk = (
            selected[selected["data_profile"] == profile]
            .groupby("formula", as_index=False)["fraction_best_tested_u_g_mean"]
            .mean()
            .sort_values("fraction_best_tested_u_g_mean", ascending=True)
        )
        positions = np.arange(len(chunk))
        ax.barh(
            positions,
            chunk["fraction_best_tested_u_g_mean"],
            color="#009E73",
            height=0.72,
        )
        ax.set_yticks(positions)
        ax.set_yticklabels(chunk["formula"], fontsize=8)
        ax.set_xlim(0.0, 1.02)
        ax.set_title(profile)
        ax.set_xlabel("mean fraction")
        ax.grid(True, axis="x", alpha=0.25)
        for y_pos, value in zip(positions, chunk["fraction_best_tested_u_g_mean"]):
            ax.text(
                min(float(value) + 0.012, 1.0),
                y_pos,
                f"{value:.3f}",
                va="center",
                fontsize=7,
            )
    for index in range(len(profiles), rows * columns):
        axes[index // columns, index % columns].axis("off")
    fig.suptitle("Rule performance by distribution")
    fig.tight_layout(rect=(0, 0, 1, 0.96), h_pad=2.0, w_pad=3.2)
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_threshold_scaling_report(
    all_distribution,
    all_stats,
    all_formula,
    scaling_formula_summary,
    metric_correlations,
    strong_weak_summary,
    out_dir,
    args,
):
    overall = scaling_formula_summary[scaling_formula_summary["scope"] == "overall"]
    best_formula = overall.sort_values("fraction_mean", ascending=False).iloc[0]
    best_family = (
        all_formula.groupby("formula_family", as_index=False)["fraction_best_tested_u_g_mean"]
        .mean()
        .sort_values("fraction_best_tested_u_g_mean", ascending=False)
        .iloc[0]
    )
    global_formula_contexts = all_formula[all_formula["formula"] == best_formula["formula"]]
    context_best = all_formula.groupby(FULL_SWEEP_CONTEXT_COLS, as_index=False)[
        "fraction_best_tested_u_g_mean"
    ].max()
    global_with_best = global_formula_contexts.merge(
        context_best,
        on=FULL_SWEEP_CONTEXT_COLS,
        suffixes=("", "_context_best"),
    )
    within_99_rate = float(
        (
            global_with_best["fraction_best_tested_u_g_mean"]
            >= 0.99 * global_with_best["fraction_best_tested_u_g_mean_context_best"]
        ).mean()
    )

    profile_shift = (
        all_stats.groupby("data_profile", as_index=False)["best_T_over_N_mean"]
        .mean()
        .sort_values("best_T_over_N_mean")
    )
    global_t_over_n = float(all_stats["best_T_over_N_mean"].mean())
    left_profiles = profile_shift[
        profile_shift["best_T_over_N_mean"] < global_t_over_n * 0.9
    ]["data_profile"].tolist()
    right_profiles = profile_shift[
        profile_shift["best_T_over_N_mean"] > global_t_over_n * 1.1
    ]["data_profile"].tolist()

    metric_overall = metric_correlations[
        (metric_correlations["scope"] == "overall")
        & (metric_correlations["target"] == "best_tested_T_over_N")
    ].dropna(subset=["abs_pearson_corr"])
    best_metric = metric_overall.iloc[0] if not metric_overall.empty else None
    strong_weak_overall = None
    if not strong_weak_summary.empty:
        strong_weak_overall = strong_weak_summary.agg(
            {
                "fraction_mean": "mean",
                "fraction_p05": "mean",
                "gap_pct_mean": "mean",
                "outside_dense_rate": "mean",
            }
        )
    report_N_values = sorted(all_stats["N"].drop_duplicates().astype(int))
    report_L_values = sorted(all_stats["L"].drop_duplicates().astype(int))
    report_active_pcts = sorted(all_stats["active_pct"].drop_duplicates().astype(float))
    report_profiles = list(all_stats["data_profile"].drop_duplicates())

    lines = [
        "# Preliminary Threshold Scaling Study",
        "",
        f"- N values: {', '.join(str(value) for value in report_N_values)}",
        f"- L values: {', '.join(str(value) for value in report_L_values)}",
        f"- Active K percentages: {', '.join(format_float(value, precision=1) for value in report_active_pcts)}",
        f"- Samples: {args.samples}",
        f"- Generator seeds: {', '.join(str(seed) for seed in args.generator_seeds)}",
        f"- Profiles: {', '.join(report_profiles)}",
        f"- Sigma: {args.sigma}",
        "",
        "## Direct Answer",
        "",
        f"- Best global tested rule: `{best_formula['formula']}` with mean fraction of best tested `U_G` `{best_formula['fraction_mean']:.4f}` and mean gap `{best_formula['gap_pct_mean']:.2f}%`.",
        f"- The same rule is within 99% of the best rule for `{within_99_rate:.1%}` of `(N, L, distribution, K%)` contexts.",
        f"- The strongest rule family is `{best_family['formula_family']}` with mean fraction `{best_family['fraction_best_tested_u_g_mean']:.4f}`.",
        f"- Mean best tested threshold scale over all contexts: `T/N={global_t_over_n:.4f}`.",
    ]
    if strong_weak_overall is not None:
        lines.append(
            f"- Strong/weak H3: mean fraction `{strong_weak_overall['fraction_mean']:.4f}`, "
            f"mean p05 fraction `{strong_weak_overall['fraction_p05']:.4f}`, "
            f"mean gap `{strong_weak_overall['gap_pct_mean']:.2f}%`; "
            f"outside dense `T=0..K` in `{strong_weak_overall['outside_dense_rate']:.1%}` of contexts."
        )
    if left_profiles or right_profiles:
        lines.append(
            "- Distribution shift by `T/N`: "
            f"left/lower `{', '.join(left_profiles) if left_profiles else 'none'}`; "
            f"right/higher `{', '.join(right_profiles) if right_profiles else 'none'}`."
        )
    else:
        lines.append("- Distribution shift by `T/N`: no profile moved more than 10% from the global mean.")
    if best_metric is not None:
        strength = "weak"
        if best_metric["abs_pearson_corr"] >= 0.5:
            strength = "strong"
        elif best_metric["abs_pearson_corr"] >= 0.25:
            strength = "moderate"
        lines.append(
            f"- Best single distribution metric predictor is `{best_metric['metric']}` "
            f"with {strength} Pearson correlation `{best_metric['pearson_corr']:.3f}` to `T/N`."
        )
    lines.extend(
        [
            "",
            "## Rule Ranking",
            "",
            "| rule | family | contexts | mean fraction | p05 fraction avg | mean gap % |",
            "|---|---|---:|---:|---:|---:|",
        ]
    )
    for _, row in overall.sort_values("fraction_mean", ascending=False).head(12).iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["formula"]),
                    str(row["formula_family"]),
                    str(int(row["contexts"])),
                    format_float(row["fraction_mean"], precision=4),
                    format_float(row["fraction_p05_mean"], precision=4),
                    format_float(row["gap_pct_mean"], precision=2),
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Best Tested Threshold Scale By Distribution",
            "",
            "| profile | mean T/N | mean T/K |",
            "|---|---:|---:|",
        ]
    )
    profile_scale = (
        all_stats.groupby("data_profile", as_index=False)
        .agg(
            best_T_over_N_mean=("best_T_over_N_mean", "mean"),
            best_T_over_K_mean=("best_T_over_K_mean", "mean"),
        )
        .sort_values("best_T_over_N_mean")
    )
    for _, row in profile_scale.iterrows():
        lines.append(
            f"| {row['data_profile']} | {row['best_T_over_N_mean']:.4f} | {row['best_T_over_K_mean']:.4f} |"
        )

    if not strong_weak_summary.empty:
        lines.extend(
            [
                "",
                "## Strong/Weak H3 By Distribution",
                "",
                "| profile | mean fraction | p05 fraction avg | mean gap % | outside dense rate |",
                "|---|---:|---:|---:|---:|",
            ]
        )
        strong_by_profile = (
            strong_weak_summary.groupby("data_profile", as_index=False)
            .agg(
                fraction_mean=("fraction_mean", "mean"),
                fraction_p05=("fraction_p05", "mean"),
                gap_pct_mean=("gap_pct_mean", "mean"),
                outside_dense_rate=("outside_dense_rate", "mean"),
            )
            .sort_values("fraction_mean", ascending=False)
        )
        for _, row in strong_by_profile.iterrows():
            lines.append(
                "| "
                + " | ".join(
                    [
                        str(row["data_profile"]),
                        format_float(row["fraction_mean"], precision=4),
                        format_float(row["fraction_p05"], precision=4),
                        format_float(row["gap_pct_mean"], precision=2),
                        f"{row['outside_dense_rate']:.1%}",
                    ]
                )
                + " |"
            )

    lines.extend(
        [
            "",
            "## Plots",
            "",
            "![Raw U_G CDF](scaling_cdf_raw_u_g.png)",
            "",
            "![Fraction CDF](scaling_cdf_fraction_best_u_g.png)",
            "",
            "![Best tested T over N](scaling_best_T_over_N.png)",
            "",
            "![Best tested T over K](scaling_best_T_over_K.png)",
            "",
            "![Formula by N and L](scaling_formula_fraction_by_N_L.png)",
            "",
            "![Formula by distribution](scaling_formula_fraction_by_profile.png)",
            "",
            "## Artifacts",
            "",
            "- `all_best_thresholds.csv`",
            "- `all_threshold_best_t_stats.csv`",
            "- `all_distribution_comparison.csv`",
            "- `all_formula_comparison.csv`",
            "- `all_scaling_formula_summary.csv`",
            "- `all_scaling_metric_correlations.csv`",
            "- `all_formula_selected_runs.csv`",
            "- `all_strong_weak_runs.csv`",
            "- `all_strong_weak_summary.csv`",
            "- Per-shard `threshold_runs.csv.gz` and report files under `N{N}/L{L}/{profile}/`.",
        ]
    )
    (out_dir / "all_scaling_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def empirical_cdf_local(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    values = np.sort(values)
    if len(values) == 0:
        return values, values
    probs = np.arange(1, len(values) + 1, dtype=float) / len(values)
    return np.r_[values[0], values], np.r_[0.0, probs]


RULE_CDF_USECOLS = [
    *FULL_SWEEP_CASE_COLS,
    "off_pct",
    "active_pct",
    "T",
    "u_bf",
    "u_i",
    "u_g",
    "u_g_db",
    "best_tested_u_g",
    "fraction_best_tested_u_g",
    "gap_to_best_tested_pct",
    "is_best_tested_T",
]

RULE_CDF_LABELS = {
    "best_tested_T": "best tested T per sample",
    "best_formula": "best formula",
    "global_0p05N": "global 0.05N",
    "strong_weak": "strong/weak H3",
}

RULE_CDF_COLORS = {
    "best_tested_T": "#0072B2",
    "best_formula": "#009E73",
    "global_0p05N": "#E69F00",
    "strong_weak": "#D55E00",
}


def run_threshold_rule_cdf_report(args):
    distribution_path = args.out_dir / "all_distribution_comparison.csv"
    formula_path = args.out_dir / "all_formula_comparison.csv"
    if not distribution_path.exists():
        raise FileNotFoundError(f"Missing full-sweep distribution file: {distribution_path}")
    if not formula_path.exists():
        raise FileNotFoundError(f"Missing full-sweep formula file: {formula_path}")

    distribution = pd.read_csv(distribution_path)
    formulas = pd.read_csv(formula_path)

    comparison_parts = []
    for profile in args.data_profiles:
        runs_path = args.out_dir / profile / "threshold_runs.csv"
        if not runs_path.exists():
            raise FileNotFoundError(f"Missing profile full-sweep runs file: {runs_path}")

        runs = pd.read_csv(runs_path, usecols=RULE_CDF_USECOLS)
        if "fraction_best_tested_u_g" not in runs or runs["fraction_best_tested_u_g"].isna().any():
            runs = add_full_sweep_relative_metrics(runs)

        profile_distribution = distribution[distribution["data_profile"] == profile]
        comparison_parts.append(
            select_threshold_rule_rows(runs, profile_distribution, args)
        )
        comparison_parts.append(
            build_strong_weak_rule_rows(runs, profile, args)
        )

    comparison = pd.concat(comparison_parts, ignore_index=True)
    comparison = comparison.sort_values(["data_profile", "K", "rule", "generator_seed", "sample"])
    summary = build_threshold_rule_summary(comparison)

    atomic_write_csv(comparison, args.out_dir / "threshold_rule_comparison.csv")
    atomic_write_csv(summary, args.out_dir / "threshold_rule_summary.csv")
    write_threshold_rule_cdf_plots(comparison, args.out_dir)
    write_threshold_rule_cdf_report(
        comparison,
        summary,
        distribution,
        formulas,
        args.out_dir,
        args,
    )


def select_threshold_rule_rows(runs, distribution, args):
    selected_parts = []
    for _, row in distribution.iterrows():
        K = int(row["K"])
        N = int(args.N)
        k_runs = runs[runs["K"] == K]
        best_mask = boolean_series(k_runs["is_best_tested_T"])
        best_chunk = (
            k_runs[best_mask]
            .sort_values([*FULL_SWEEP_CASE_COLS, "T"])
            .drop_duplicates(FULL_SWEEP_CASE_COLS, keep="first")
            .copy()
        )
        if not best_chunk.empty:
            best_chunk["rule"] = "best_tested_T"
            best_chunk["rule_label"] = RULE_CDF_LABELS["best_tested_T"]
            best_chunk["rule_T"] = best_chunk["T"].astype(int)
            best_chunk["rule_note"] = "sample_specific_best_tested_T"
            best_chunk["fraction_best_tested_u_g"] = 1.0
            best_chunk["gap_to_best_tested_pct"] = 0.0
            best_chunk["is_best_tested_T"] = True
            selected_parts.append(best_chunk)

        rules = [
            ("best_formula", int(row["best_formula_T"])),
            ("global_0p05N", int(np.clip(round(0.05 * N), 0, K))),
        ]
        for rule, T in rules:
            chunk = k_runs[k_runs["T"] == T].copy()
            if chunk.empty:
                continue
            chunk["rule"] = rule
            chunk["rule_label"] = RULE_CDF_LABELS[rule]
            chunk["rule_T"] = int(T)
            chunk["rule_note"] = ""
            selected_parts.append(chunk)

    if not selected_parts:
        return pd.DataFrame(columns=[*RULE_CDF_USECOLS, "rule", "rule_label", "rule_T", "rule_note"])
    return pd.concat(selected_parts, ignore_index=True)


def boolean_series(series):
    if series.dtype == object:
        return series.astype(str).str.lower().isin(("true", "1", "yes"))
    return series.astype(bool)


def build_strong_weak_rule_rows(reference_runs, profile, args):
    best_lookup = (
        reference_runs.drop_duplicates(FULL_SWEEP_CASE_COLS)
        .set_index(FULL_SWEEP_CASE_COLS)["best_tested_u_g"]
        .to_dict()
    )
    rows = []

    total_cases = len(args.generator_seeds) * args.samples * len(args.off_cases)
    case_no = 0
    for generator_seed in args.generator_seeds:
        rng = np.random.RandomState(generator_seed)
        for sample in range(args.samples):
            V = generate_v_profile_from_rng(rng, args.N, args.L, profile)
            for off_case in args.off_cases:
                case_no += 1
                K = int(off_case["K"])
                off_pct = float(off_case["off_pct"])
                u_bf, u_i, u_g, T = evaluate_strong_weak_window(
                    V,
                    K,
                    sigma=args.sigma,
                    P=args.P,
                )
                case_key = (
                    profile,
                    int(generator_seed),
                    int(sample),
                    int(args.N),
                    int(args.L),
                    K,
                    float(args.sigma),
                    float(args.P),
                )
                best_u_g = float(best_lookup[case_key])
                fraction = u_g / best_u_g if best_u_g > 0 else np.nan
                if case_no % 1000 == 0:
                    print(
                        f"strong/weak profile={profile}: {case_no}/{total_cases} cases",
                        flush=True,
                    )
                rows.append(
                    {
                        "data_profile": profile,
                        "generator_seed": int(generator_seed),
                        "sample": int(sample),
                        "N": int(args.N),
                        "L": int(args.L),
                        "K": K,
                        "sigma": float(args.sigma),
                        "P": float(args.P),
                        "off_pct": off_pct,
                        "active_pct": float(100.0 - off_pct),
                        "T": int(T),
                        "u_bf": float(u_bf),
                        "u_i": float(u_i),
                        "u_g": float(u_g),
                        "u_g_db": 10.0 * np.log10(max(float(u_g), np.finfo(float).tiny)),
                        "best_tested_u_g": best_u_g,
                        "fraction_best_tested_u_g": fraction,
                        "gap_to_best_tested_pct": 100.0 * (1.0 - fraction),
                        "is_best_tested_T": bool(fraction >= 1.0 - 1e-12),
                        "rule": "strong_weak",
                        "rule_label": RULE_CDF_LABELS["strong_weak"],
                        "rule_T": int(T),
                        "rule_note": "outside_dense_0_to_K" if T > K else "",
                    }
                )

    return pd.DataFrame(rows)


def evaluate_strong_weak_window(V, K, sigma=1.0, P=1.0):
    V = np.asarray(V)
    N, L = V.shape
    K = int(np.clip(K, 0, N))
    if K == 0:
        return 0.0, 0.0, 0.0, 0
    if K == N:
        gram = V.conj().T @ V
        p_n = np.sum(np.abs(V) ** 2, axis=1).real
        return (*objective_from_gram(gram, float(np.max(p_n)), L, sigma=sigma, P=P), 0)

    off_count = N - K
    weak_drop = off_count // 2
    strong_drop = off_count - weak_drop
    T = int(strong_drop)

    p_n = np.sum(np.abs(V) ** 2, axis=1).real
    idx_desc = np.argsort(p_n)[::-1]
    active_idx = idx_desc[T : T + K]
    gram = V[active_idx].conj().T @ V[active_idx]
    max_row_power = float(p_n[idx_desc[T]])
    return (*objective_from_gram(gram, max_row_power, L, sigma=sigma, P=P), T)


def collect_threshold_exact_study_runs(args):
    threshold_rows_all = []
    exact_rows = []
    formula_rows = []
    total_cases = (
        len(args.N_values)
        * len(args.data_profiles)
        * len(args.generator_seeds)
        * args.samples
        * len(args.K_pcts)
    )
    progress = (
        tqdm(total=total_cases, unit="case", dynamic_ncols=True)
        if tqdm is not None
        else None
    )

    try:
        for N in args.N_values:
            N = int(N)
            off_cases = build_active_pct_cases(N, args.K_pcts)
            for profile in args.data_profiles:
                for generator_seed in args.generator_seeds:
                    rng = np.random.RandomState(int(generator_seed))
                    for sample in range(args.samples):
                        V = generate_v_profile_from_rng(rng, N, args.L, profile)
                        distribution_metrics = row_power_distribution_metrics(V)

                        for off_case in off_cases:
                            K = int(off_case["K"])
                            message = (
                                f"N={N}, L={args.L}, profile={profile}, "
                                f"seed={generator_seed}, sample={sample}, K={K}"
                            )
                            if progress is not None:
                                progress.set_postfix_str(message)
                            else:
                                print(message, flush=True)

                            exact = brute_force_exact_u_g(
                                V,
                                K,
                                sigma=args.sigma,
                                P=args.P,
                                time_limit_seconds=args.exact_time_limit,
                            )
                            exact_u_g = float(exact["u_g"])
                            exact_completed = bool(exact["completed"])
                            exact_window_T = (
                                contiguous_threshold_window_T(V, exact["subset"])
                                if exact_completed
                                else None
                            )

                            thresholds = dense_thresholds(K)
                            threshold_rows = evaluate_power_window_thresholds(
                                V,
                                K,
                                thresholds,
                                sigma=args.sigma,
                                P=args.P,
                            )
                            threshold_rows = sorted(
                                threshold_rows,
                                key=lambda row: (-float(row["u_g"]), int(row["T"])),
                            )
                            best_threshold = threshold_rows[0]
                            threshold_by_T = {
                                int(row["T"]): row for row in threshold_rows
                            }

                            case_base = {
                                "data_profile": profile,
                                "generator_seed": int(generator_seed),
                                "sample": int(sample),
                                "N": N,
                                "L": int(args.L),
                                "K": K,
                                "off_pct": float(off_case["off_pct"]),
                                "active_pct": float(off_case["active_pct"]),
                                "sigma": float(args.sigma),
                                "P": float(args.P),
                            }
                            exact_base = {
                                "exact_completed": exact_completed,
                                "exact_timed_out": bool(exact["timed_out"]),
                                "exact_candidate_count": int(exact["candidate_count"]),
                                "exact_evaluated_count": int(exact["evaluated_count"]),
                                "exact_elapsed_seconds": float(exact["elapsed_seconds"]),
                                "exact_subset": subset_to_string(exact["subset"]),
                                "exact_window_T": np.nan
                                if exact_window_T is None
                                else int(exact_window_T),
                                "exact_is_threshold_window": bool(
                                    exact_window_T is not None and exact_completed
                                ),
                                "exact_u_bf": float(exact["u_bf"]),
                                "exact_u_i": float(exact["u_i"]),
                                "exact_u_g": exact_u_g,
                                "exact_u_g_db": 10.0
                                * np.log10(max(exact_u_g, np.finfo(float).tiny))
                                if exact_completed
                                else np.nan,
                            }

                            best_u_g = float(best_threshold["u_g"])
                            threshold_fraction_exact = (
                                best_u_g / exact_u_g
                                if exact_completed and exact_u_g > 0
                                else np.nan
                            )
                            exact_rows.append(
                                {
                                    **case_base,
                                    **exact_base,
                                    "best_tested_T": int(best_threshold["T"]),
                                    "best_tested_u_bf": float(best_threshold["u_bf"]),
                                    "best_tested_u_i": float(best_threshold["u_i"]),
                                    "best_tested_u_g": best_u_g,
                                    "best_tested_u_g_db": 10.0
                                    * np.log10(max(best_u_g, np.finfo(float).tiny)),
                                    "best_tested_subset": threshold_window_subset_string(
                                        V,
                                        K,
                                        int(best_threshold["T"]),
                                    ),
                                    "best_tested_fraction_exact_u_g": threshold_fraction_exact,
                                    "best_tested_gap_to_exact_pct": 100.0
                                    * (1.0 - threshold_fraction_exact)
                                    if np.isfinite(threshold_fraction_exact)
                                    else np.nan,
                                }
                            )

                            for row in threshold_rows:
                                row_u_g = float(row["u_g"])
                                threshold_rows_all.append(
                                    {
                                        **case_base,
                                        "target_obj": "gen",
                                        "threshold_source": "dense_0_to_K",
                                        **row,
                                        "exact_completed": exact_completed,
                                        "exact_u_g": exact_u_g
                                        if exact_completed
                                        else np.nan,
                                        "fraction_exact_u_g": row_u_g / exact_u_g
                                        if exact_completed and exact_u_g > 0
                                        else np.nan,
                                        "gap_to_exact_pct": 100.0
                                        * (1.0 - row_u_g / exact_u_g)
                                        if exact_completed and exact_u_g > 0
                                        else np.nan,
                                        **distribution_metrics,
                                    }
                                )

                            formula_rows.extend(
                                exact_formula_rows_for_case(
                                    case_base,
                                    threshold_by_T,
                                    best_threshold,
                                    exact,
                                    V,
                                    K,
                                    args,
                                )
                            )

                            if progress is not None:
                                progress.update(1)
    finally:
        if progress is not None:
            progress.close()

    threshold_runs = pd.DataFrame(threshold_rows_all)
    if not threshold_runs.empty:
        threshold_runs = add_full_sweep_relative_metrics(threshold_runs)
    return threshold_runs, pd.DataFrame(exact_rows), pd.DataFrame(formula_rows)


def exact_formula_rows_for_case(case_base, threshold_by_T, best_threshold, exact, V, K, args):
    rows = []
    exact_completed = bool(exact["completed"])
    exact_u_g = float(exact["u_g"])
    best_rule_specs = [
        (
            "best_tested_T",
            "best_tested",
            "best tested T per sample",
            int(best_threshold["T"]),
            best_threshold,
            "",
        )
    ]
    for label, _ in SCALING_FORMULA_SPECS:
        T = threshold_formula_T(
            label,
            int(case_base["N"]),
            int(case_base["L"]),
            int(K),
            formula_specs=SCALING_FORMULA_SPECS,
        )
        selected = threshold_by_T.get(T)
        if selected is None:
            continue
        best_rule_specs.append((label, formula_family(label), label, T, selected, ""))

    u_bf, u_i, u_g, strong_weak_T = evaluate_strong_weak_window(
        V,
        K,
        sigma=args.sigma,
        P=args.P,
    )
    best_rule_specs.append(
        (
            STRONG_WEAK_RULE,
            formula_family(STRONG_WEAK_RULE),
            "strong/weak H3",
            int(strong_weak_T),
            {
                "T": int(strong_weak_T),
                "u_bf": float(u_bf),
                "u_i": float(u_i),
                "u_g": float(u_g),
            },
            "outside_dense_0_to_K" if int(strong_weak_T) > int(K) else "",
        )
    )

    for formula, family, label, T, selected, note in best_rule_specs:
        selected_u_g = float(selected["u_g"])
        fraction = (
            selected_u_g / exact_u_g
            if exact_completed and exact_u_g > 0.0
            else np.nan
        )
        rows.append(
            {
                **case_base,
                "formula": formula,
                "formula_family": family,
                "formula_label": label,
                "formula_T": int(T),
                "rule_note": note,
                "u_bf": float(selected["u_bf"]),
                "u_i": float(selected["u_i"]),
                "u_g": selected_u_g,
                "u_g_db": 10.0 * np.log10(max(selected_u_g, np.finfo(float).tiny)),
                "exact_completed": exact_completed,
                "exact_u_g": exact_u_g if exact_completed else np.nan,
                "fraction_exact_u_g": fraction,
                "gap_to_exact_pct": 100.0 * (1.0 - fraction)
                if np.isfinite(fraction)
                else np.nan,
                "is_exact_best": bool(fraction >= 1.0 - 1e-12)
                if np.isfinite(fraction)
                else False,
            }
        )
    return rows


def build_threshold_exact_summary(exact_runs):
    if exact_runs.empty:
        return pd.DataFrame()

    def q(value):
        return lambda data: data.quantile(value)

    return (
        exact_runs.groupby(FULL_SWEEP_CONTEXT_COLS, as_index=False)
        .agg(
            cases=("best_tested_u_g", "count"),
            exact_completed_rate=("exact_completed", "mean"),
            exact_timeout_rate=("exact_timed_out", "mean"),
            exact_candidate_count_mean=("exact_candidate_count", "mean"),
            exact_elapsed_mean=("exact_elapsed_seconds", "mean"),
            exact_elapsed_p95=("exact_elapsed_seconds", q(0.95)),
            best_tested_T_min=("best_tested_T", "min"),
            best_tested_T_p05=("best_tested_T", q(0.05)),
            best_tested_T_p50=("best_tested_T", q(0.50)),
            best_tested_T_p95=("best_tested_T", q(0.95)),
            best_tested_T_max=("best_tested_T", "max"),
            best_tested_T_mean=("best_tested_T", "mean"),
            threshold_fraction_exact_mean=("best_tested_fraction_exact_u_g", "mean"),
            threshold_fraction_exact_p05=(
                "best_tested_fraction_exact_u_g",
                q(0.05),
            ),
            threshold_fraction_exact_p50=(
                "best_tested_fraction_exact_u_g",
                q(0.50),
            ),
            threshold_fraction_exact_p95=(
                "best_tested_fraction_exact_u_g",
                q(0.95),
            ),
            threshold_exact_rate=(
                "best_tested_fraction_exact_u_g",
                lambda values: (values >= 1.0 - 1e-9).mean(),
            ),
            threshold_near_99_rate=(
                "best_tested_fraction_exact_u_g",
                lambda values: (values >= 0.99).mean(),
            ),
            exact_is_threshold_window_rate=("exact_is_threshold_window", "mean"),
        )
        .sort_values(["data_profile", "N", "L", "K"])
    )


def build_threshold_exact_formula_summary(formula_runs):
    if formula_runs.empty:
        return pd.DataFrame()

    def q(value):
        return lambda data: data.quantile(value)

    return (
        formula_runs.groupby([*FULL_SWEEP_CONTEXT_COLS, "formula", "formula_family"], as_index=False)
        .agg(
            cases=("u_g", "count"),
            T_min=("formula_T", "min"),
            T_p50=("formula_T", q(0.50)),
            T_max=("formula_T", "max"),
            outside_dense_rate=(
                "rule_note",
                lambda values: (values.astype(str) == "outside_dense_0_to_K").mean(),
            ),
            u_g_mean=("u_g", "mean"),
            fraction_exact_mean=("fraction_exact_u_g", "mean"),
            fraction_exact_p05=("fraction_exact_u_g", q(0.05)),
            fraction_exact_p50=("fraction_exact_u_g", q(0.50)),
            fraction_exact_p95=("fraction_exact_u_g", q(0.95)),
            gap_to_exact_pct_mean=("gap_to_exact_pct", "mean"),
            exact_rate=("is_exact_best", "mean"),
        )
        .sort_values(["data_profile", "N", "L", "K", "formula"])
    )


def run_threshold_exact_study(args):
    threshold_runs, exact_runs, formula_runs = collect_threshold_exact_study_runs(args)
    summary = build_threshold_exact_summary(exact_runs)
    formula_summary = build_threshold_exact_formula_summary(formula_runs)

    atomic_write_csv(threshold_runs, args.out_dir / "threshold_runs.csv.gz")
    atomic_write_csv(exact_runs, args.out_dir / "exact_runs.csv")
    atomic_write_csv(formula_runs, args.out_dir / "exact_formula_runs.csv")
    atomic_write_csv(summary, args.out_dir / "exact_summary.csv")
    atomic_write_csv(formula_summary, args.out_dir / "exact_formula_summary.csv")
    write_threshold_exact_plots(exact_runs, formula_runs, formula_summary, args.out_dir)
    write_threshold_exact_report(exact_runs, summary, formula_summary, args.out_dir, args)
    write_threshold_exact_k_analysis(exact_runs, formula_runs, args.out_dir)
    return summary


def write_threshold_exact_plots(exact_runs, formula_runs, formula_summary, out_dir):
    if exact_runs.empty:
        return

    from utils.plotting import use_agg_backend

    use_agg_backend()
    import matplotlib.pyplot as plt

    selected_formulas = [
        "best_tested_T",
        "T_0p05N",
        "T_0p15K",
        "T_0p125NL_over_Lp2",
        STRONG_WEAK_RULE,
    ]
    selected = formula_runs[formula_runs["formula"].isin(selected_formulas)]

    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    for formula in selected_formulas:
        chunk = selected[selected["formula"] == formula]
        x_values, y_values = empirical_cdf_local(chunk["fraction_exact_u_g"])
        if len(x_values):
            ax.step(y_values, x_values, where="post", label=formula)
    ax.set_xlabel("cumulative fraction")
    ax.set_ylabel("U_G / exact best U_G")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.05)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    ax.set_title("Exact-study CDF: fraction of exact U_G")
    fig.tight_layout()
    fig.savefig(out_dir / "exact_fraction_cdf.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    raw_u_g_formulas = [
        "best_tested_T",
        "T_0p05N",
        STRONG_WEAK_RULE,
    ]
    raw_selected = formula_runs[formula_runs["formula"].isin(raw_u_g_formulas)]

    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    for formula in raw_u_g_formulas:
        chunk = raw_selected[raw_selected["formula"] == formula]
        x_values, y_values = empirical_cdf_local(chunk["u_g"])
        if len(x_values):
            ax.step(y_values, x_values, where="post", label=formula)
    exact_completed = exact_runs[boolean_series(exact_runs["exact_completed"])]
    x_values, y_values = empirical_cdf_local(exact_completed["exact_u_g"])
    if len(x_values):
        ax.step(y_values, x_values, where="post", label="exact best", color="#000000")
    ax.set_xlabel("cumulative fraction")
    ax.set_ylabel("raw U_G")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    ax.set_title("Exact-study CDF: raw U_G")
    fig.tight_layout()
    fig.savefig(out_dir / "exact_raw_u_g_cdf.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    labels = []
    data = []
    for (N, K), chunk in exact_runs.groupby(["N", "K"], sort=True):
        labels.append(f"N={int(N)}\nK={int(K)}")
        data.append(chunk["best_tested_T"].astype(float))
    fig, ax = plt.subplots(figsize=(max(7.4, 0.85 * len(labels)), 4.6))
    ax.boxplot(data, tick_labels=labels, showmeans=True)
    ax.set_ylabel("best tested T")
    ax.set_title("Best tested threshold by exact-study case")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_dir / "exact_best_T_boxplot.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    for K, chunk in exact_runs.groupby("K", sort=True):
        grouped = chunk.groupby("N", as_index=False)["exact_elapsed_seconds"].mean()
        ax.plot(grouped["N"], grouped["exact_elapsed_seconds"], marker="o", label=f"K={int(K)}")
    ax.set_xlabel("N")
    ax.set_ylabel("mean exact elapsed seconds")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.25)
    ax.legend()
    ax.set_title("Brute-force exact runtime")
    fig.tight_layout()
    fig.savefig(out_dir / "exact_runtime_by_N.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    if not formula_summary.empty:
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        summary = formula_summary[formula_summary["formula"].isin(selected_formulas)]
        for formula in selected_formulas:
            chunk = summary[summary["formula"] == formula]
            if chunk.empty:
                continue
            grouped = chunk.groupby("N", as_index=False)["fraction_exact_mean"].mean()
            ax.plot(grouped["N"], grouped["fraction_exact_mean"], marker="o", label=formula)
        ax.set_xlabel("N")
        ax.set_ylabel("mean U_G / exact best U_G")
        ax.set_ylim(0.0, 1.02)
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
        ax.set_title("Formula quality against exact optimum")
        fig.tight_layout()
        fig.savefig(out_dir / "exact_formula_fraction_by_N.png", dpi=180, bbox_inches="tight")
        plt.close(fig)


def write_threshold_exact_report(exact_runs, summary, formula_summary, out_dir, args):
    lines = [
        "# Exact Threshold Approach Study",
        "",
        f"- N values: {', '.join(str(value) for value in args.N_values)}",
        f"- L: {args.L}",
        f"- Active K percentages: {', '.join(format_float(value) for value in args.K_pcts)}",
        f"- Samples: {args.samples}",
        f"- Generator seeds: {', '.join(str(value) for value in args.generator_seeds)}",
        f"- Profiles: {', '.join(args.data_profiles)}",
        f"- Sigma: {args.sigma}",
        f"- Exact time limit: {args.exact_time_limit} seconds",
        "",
        "The exact solver enumerates every subset of size `K` and maximizes raw `U_G`.",
        "The threshold comparison uses the best tested shifted window from `T=0..K`.",
        "",
        "## Direct Answer",
        "",
    ]

    completed_rate = float(exact_runs["exact_completed"].mean()) if not exact_runs.empty else 0.0
    if not summary.empty:
        threshold_mean = float(summary["threshold_fraction_exact_mean"].mean())
        window_rate = float(summary["exact_is_threshold_window_rate"].mean())
        near_99 = float(summary["threshold_near_99_rate"].mean())
        lines.extend(
            [
                f"- Exact enumeration completed for `{completed_rate:.1%}` of cases.",
                f"- Best tested threshold-window mean fraction of exact `U_G`: `{threshold_mean:.4f}`.",
                f"- Fraction of cases where threshold window is within 99% of exact: `{near_99:.1%}` on average by context.",
                f"- Exact optimum was itself a contiguous row-power window in `{window_rate:.1%}` of cases on average by context.",
            ]
        )
    else:
        lines.append("- No exact-study rows were produced.")

    lines.extend(
        [
            "",
            "## Threshold-vs-Exact Summary",
            "",
            "| profile | N | K | exact completed | candidates | exact time mean | best T p50 | threshold/exact mean | threshold/exact p05 | exact-window rate |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for _, row in summary.iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["data_profile"]),
                    str(int(row["N"])),
                    str(int(row["K"])),
                    f"{row['exact_completed_rate']:.1%}",
                    format_float(row["exact_candidate_count_mean"], precision=0),
                    format_float(row["exact_elapsed_mean"]),
                    format_float(row["best_tested_T_p50"]),
                    format_float(row["threshold_fraction_exact_mean"], precision=4),
                    format_float(row["threshold_fraction_exact_p05"], precision=4),
                    f"{row['exact_is_threshold_window_rate']:.1%}",
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Formula And Strong/Weak Comparison",
            "",
            "| formula | mean fraction exact | p05 fraction | exact rate | outside dense rate |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    if not formula_summary.empty:
        global_formula = (
            formula_summary.groupby("formula", as_index=False)
            .agg(
                fraction_exact_mean=("fraction_exact_mean", "mean"),
                fraction_exact_p05=("fraction_exact_p05", "mean"),
                exact_rate=("exact_rate", "mean"),
                outside_dense_rate=("outside_dense_rate", "mean"),
            )
            .sort_values("fraction_exact_mean", ascending=False)
        )
        for _, row in global_formula.iterrows():
            lines.append(
                "| "
                + " | ".join(
                    [
                        str(row["formula"]),
                        format_float(row["fraction_exact_mean"], precision=4),
                        format_float(row["fraction_exact_p05"], precision=4),
                        f"{row['exact_rate']:.1%}",
                        f"{row['outside_dense_rate']:.1%}",
                    ]
                )
                + " |"
            )

    completed = exact_runs[boolean_series(exact_runs["exact_completed"])] if not exact_runs.empty else pd.DataFrame()
    if not completed.empty:
        best_examples = completed[
            completed["best_tested_fraction_exact_u_g"] >= 1.0 - 1e-9
        ].head(5)
        worst_examples = completed.sort_values("best_tested_fraction_exact_u_g").head(5)
        lines.extend(["", "## Exact Best Cases Found", ""])
        if best_examples.empty:
            lines.append("- No exact equality cases were found in this run.")
        else:
            for _, row in best_examples.iterrows():
                lines.append(
                    f"- `N={int(row['N'])}`, `K={int(row['K'])}`, "
                    f"sample `{int(row['sample'])}`: best tested `T={int(row['best_tested_T'])}` "
                    f"matches exact `U_G`."
                )

        lines.extend(["", "## Worst Threshold-Window Cases Found", ""])
        for _, row in worst_examples.iterrows():
            lines.append(
                f"- `N={int(row['N'])}`, `K={int(row['K'])}`, sample `{int(row['sample'])}`: "
                f"best tested `T={int(row['best_tested_T'])}`, "
                f"fraction exact `{row['best_tested_fraction_exact_u_g']:.4f}`, "
                f"exact-window `{bool(row['exact_is_threshold_window'])}`."
            )

    lines.extend(
        [
            "",
            "## Plots",
            "",
            "![Exact fraction CDF](exact_fraction_cdf.png)",
            "",
            "![Raw U_G CDF](exact_raw_u_g_cdf.png)",
            "",
            "![Best T boxplot](exact_best_T_boxplot.png)",
            "",
            "![Exact runtime by N](exact_runtime_by_N.png)",
            "",
            "![Formula fraction by N](exact_formula_fraction_by_N.png)",
            "",
            "## Artifacts",
            "",
            "- `threshold_runs.csv.gz`",
            "- `exact_runs.csv`",
            "- `exact_formula_runs.csv`",
            "- `exact_summary.csv`",
            "- `exact_formula_summary.csv`",
        ]
    )
    (out_dir / "threshold_exact_report.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def build_threshold_rule_summary(comparison):
    def q(value):
        return lambda data: data.quantile(value)

    summary = (
        comparison.groupby(["data_profile", "K", "off_pct", "rule", "rule_label"], as_index=False)
        .agg(
            cases=("u_g", "count"),
            T_min=("rule_T", "min"),
            T_max=("rule_T", "max"),
            T_mean=("rule_T", "mean"),
            T_std=("rule_T", "std"),
            T_p05=("rule_T", q(0.05)),
            T_p50=("rule_T", q(0.50)),
            T_p95=("rule_T", q(0.95)),
            u_g_mean=("u_g", "mean"),
            u_g_p05=("u_g", q(0.05)),
            u_g_p50=("u_g", q(0.50)),
            u_g_p95=("u_g", q(0.95)),
            u_g_db_mean=("u_g_db", "mean"),
            fraction_mean=("fraction_best_tested_u_g", "mean"),
            fraction_p05=("fraction_best_tested_u_g", q(0.05)),
            fraction_p50=("fraction_best_tested_u_g", q(0.50)),
            fraction_p95=("fraction_best_tested_u_g", q(0.95)),
            gap_pct_mean=("gap_to_best_tested_pct", "mean"),
            winner_rate=("is_best_tested_T", "mean"),
            outside_dense_rate=("rule_note", lambda values: (values == "outside_dense_0_to_K").mean()),
        )
    )
    order = {rule: index for index, rule in enumerate(RULE_CDF_LABELS)}
    summary["_rule_order"] = summary["rule"].map(order).fillna(99)
    return summary.sort_values(["data_profile", "K", "_rule_order"]).drop(columns="_rule_order")


def write_threshold_rule_cdf_plots(comparison, out_dir):
    from utils.plotting import use_agg_backend

    use_agg_backend()
    import matplotlib.pyplot as plt

    profiles = list(comparison["data_profile"].drop_duplicates())
    Ks = sorted(comparison["K"].unique())

    write_threshold_rule_cdf_grid(
        comparison,
        profiles,
        Ks,
        value_col="fraction_best_tested_u_g",
        ylabel="fraction of best tested U_G",
        title="CDF: rule performance vs best tested threshold",
        out_path=out_dir / "threshold_rule_cdf_fraction.png",
        yscale="linear",
        ylim=(0.0, max(1.05, min(1.35, comparison["fraction_best_tested_u_g"].max() * 1.02))),
    )
    write_threshold_rule_cdf_grid(
        comparison,
        profiles,
        Ks,
        value_col="u_g",
        ylabel="raw U_G",
        title="CDF: raw U_G by threshold rule",
        out_path=out_dir / "threshold_rule_cdf_u_g.png",
        yscale="log",
        ylim=None,
    )

    for profile in profiles:
        profile_dir = out_dir / profile
        profile_data = comparison[comparison["data_profile"] == profile]
        write_threshold_rule_cdf_grid(
            profile_data,
            [profile],
            Ks,
            value_col="fraction_best_tested_u_g",
            ylabel="fraction of best tested U_G",
            title=f"CDF: {profile} rule performance",
            out_path=profile_dir / "threshold_rule_cdf_fraction.png",
            yscale="linear",
            ylim=(0.0, max(1.05, min(1.35, profile_data["fraction_best_tested_u_g"].max() * 1.02))),
        )
        write_threshold_rule_cdf_grid(
            profile_data,
            [profile],
            Ks,
            value_col="u_g",
            ylabel="raw U_G",
            title=f"CDF: {profile} raw U_G",
            out_path=profile_dir / "threshold_rule_cdf_u_g.png",
            yscale="log",
            ylim=None,
        )


def write_threshold_rule_cdf_grid(
    comparison,
    profiles,
    Ks,
    value_col,
    ylabel,
    title,
    out_path,
    yscale="linear",
    ylim=None,
):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(
        len(profiles),
        len(Ks),
        figsize=(6.2 * len(Ks), 3.2 * len(profiles)),
        squeeze=False,
    )
    for row_index, profile in enumerate(profiles):
        for col_index, K in enumerate(Ks):
            ax = axes[row_index, col_index]
            chunk = comparison[(comparison["data_profile"] == profile) & (comparison["K"] == K)]
            plotted_T = {}
            for rule in RULE_CDF_LABELS:
                rule_chunk = chunk[chunk["rule"] == rule]
                if rule_chunk.empty:
                    continue
                if rule == "best_tested_T":
                    T = None
                else:
                    T = int(rule_chunk["rule_T"].iloc[0])
                if rule == "global_0p05N" and T in plotted_T.values():
                    continue
                values, probs = empirical_cdf_local(rule_chunk[value_col])
                if len(values):
                    if rule == "best_tested_T":
                        label = f"{RULE_CDF_LABELS[rule]} (T varies)"
                    else:
                        label = f"{RULE_CDF_LABELS[rule]} (T={T})"
                    ax.step(
                        probs,
                        values,
                        where="post",
                        linewidth=1.6,
                        color=RULE_CDF_COLORS.get(rule),
                        label=label,
                    )
                    if T is not None:
                        plotted_T[rule] = T
            ax.set_title(f"{profile}, K={int(K)}")
            ax.set_xlabel("cumulative fraction")
            ax.set_ylabel(ylabel)
            ax.set_yscale(yscale)
            ax.set_xlim(0.0, 1.0)
            if ylim is not None:
                ax.set_ylim(*ylim)
            ax.grid(True, alpha=0.25)
            ax.legend(fontsize=7)
    fig.suptitle(title)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_threshold_rule_cdf_report(comparison, summary, distribution, formulas, out_dir, args):
    global_formula = (
        formulas.groupby("formula", as_index=False)
        .agg(fraction=("fraction_best_tested_u_g_mean", "mean"))
        .sort_values("fraction", ascending=False)
    )
    best_global_formula = global_formula.iloc[0] if not global_formula.empty else None
    profile_formula_choices = describe_profile_formula_choices(distribution)

    rule_overall = (
        summary.groupby("rule", as_index=False)
        .agg(
            fraction_mean=("fraction_mean", "mean"),
            fraction_p05=("fraction_p05", "mean"),
            gap_pct_mean=("gap_pct_mean", "mean"),
        )
        .sort_values("fraction_mean", ascending=False)
    )

    lines = [
        "# Threshold Rule CDF Comparison",
        "",
        f"- N: {args.N}",
        f"- L: {args.L}",
        f"- K values: {', '.join(str(case['K']) for case in args.off_cases)}",
        f"- Samples: {args.samples}",
        f"- Generator seeds: {', '.join(str(seed) for seed in args.generator_seeds)}",
        f"- Profiles: {', '.join(args.data_profiles)}",
        f"- Sigma: {args.sigma}",
        "",
        "## Direct Answer",
        "",
    ]
    if best_global_formula is not None:
        lines.append(
            f"- Best simple global formula in this run is `{best_global_formula['formula']}` "
            f"with mean fraction of best tested `U_G` `{best_global_formula['fraction']:.4f}`."
        )
    if profile_formula_choices:
        lines.append(
            "- Profile-aware best simple formula choices in this run: "
            f"{profile_formula_choices}."
        )
    lines.extend(
        [
            "- The blue reference curve is the sample-specific best tested `T` from the dense sweep.",
            "- Strong/weak H3 is a weak match for this `U_G` objective in the dense-threshold study; it often shifts too far into the middle of the power ranking.",
            "",
            "## Overall Rule Ranking",
            "",
            "| rule | mean fraction | p05 fraction avg | mean gap % |",
            "|---|---:|---:|---:|",
        ]
    )
    for _, row in rule_overall.iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    RULE_CDF_LABELS.get(row["rule"], str(row["rule"])),
                    format_float(row["fraction_mean"], precision=4),
                    format_float(row["fraction_p05"], precision=4),
                    format_float(row["gap_pct_mean"], precision=2),
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Distribution Summary",
            "",
            "| profile | K | best tested T median | best tested T p05..p95 | best formula | formula T | formula mean frac | strong/weak T | strong/weak mean frac | strong/weak p05 frac |",
            "|---|---:|---:|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    for _, dist_row in distribution.sort_values(["data_profile", "K"]).iterrows():
        profile = dist_row["data_profile"]
        K = int(dist_row["K"])
        best_tested = summary[
            (summary["data_profile"] == profile)
            & (summary["K"] == K)
            & (summary["rule"] == "best_tested_T")
        ].iloc[0]
        best_formula = summary[
            (summary["data_profile"] == profile)
            & (summary["K"] == K)
            & (summary["rule"] == "best_formula")
        ].iloc[0]
        strong_weak = summary[
            (summary["data_profile"] == profile)
            & (summary["K"] == K)
            & (summary["rule"] == "strong_weak")
        ].iloc[0]
        lines.append(
            "| "
            + " | ".join(
                [
                    str(profile),
                    str(K),
                    format_float(best_tested["T_p50"], precision=1),
                    f"{format_float(best_tested['T_p05'], precision=1)}..{format_float(best_tested['T_p95'], precision=1)}",
                    str(dist_row["best_formula"]),
                    str(int(dist_row["best_formula_T"])),
                    format_float(best_formula["fraction_mean"], precision=4),
                    str(int(strong_weak["T_min"])),
                    format_float(strong_weak["fraction_mean"], precision=4),
                    format_float(strong_weak["fraction_p05"], precision=4),
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## CDF Evidence",
            "",
            "![CDF fraction comparison](threshold_rule_cdf_fraction.png)",
            "",
            "![CDF raw U_G comparison](threshold_rule_cdf_u_g.png)",
            "",
            "## Artifacts",
            "",
            "- `threshold_rule_comparison.csv`",
            "- `threshold_rule_summary.csv`",
            "- `threshold_rule_cdf_fraction.png`",
            "- `threshold_rule_cdf_u_g.png`",
            "- Per-distribution `threshold_rule_cdf_fraction.png` and `threshold_rule_cdf_u_g.png` files.",
        ]
    )

    (out_dir / "threshold_rule_comparison_report.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def describe_profile_formula_choices(distribution):
    choices = []
    for profile, chunk in distribution.sort_values(["data_profile", "K"]).groupby("data_profile"):
        per_k = []
        for _, row in chunk.iterrows():
            if pd.isna(row["best_formula_T"]):
                continue
            per_k.append(
                f"K={int(row['K'])}: {row['best_formula']} (T={int(row['best_formula_T'])})"
            )
        if per_k:
            choices.append(f"{profile}: {', '.join(per_k)}")
    return "; ".join(choices)


def main():
    args = parse_args()
    args.off_cases = build_off_cases(args)
    args.out_dir = args.out_dir or default_out_dir(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.threshold_local_exact_analysis:
        if (
            args.threshold_exact_study
            or args.threshold_scaling_study
            or args.threshold_full_sweep
            or args.threshold_explore
            or args.threshold_rule_cdf
        ):
            raise ValueError(
                "Use --threshold-local-exact-analysis by itself, not with other threshold modes."
            )
        run_local_threshold_exact_analysis(
            exact_dir=args.exact_source_dir,
            out_dir=args.out_dir,
            docs_path=Path("docs/local_threshold_exact_gauss_report.md"),
            n_values=args.N_values,
            k_pcts=args.K_pcts,
            profiles=args.data_profiles,
        )
        print("Wrote:", flush=True)
        for output_name in (
            "local_threshold_runs.csv",
            "local_threshold_summary.csv",
            "local_threshold_failure_cases.csv",
            "local_threshold_diagnostics.csv",
            "local_threshold_exact_gauss_report.md",
            "local_threshold_raw_u_g_cdf_by_k_pct.png",
            "local_threshold_fraction_exact_cdf.png",
            "local_threshold_mean_fraction_by_active_pct.png",
            "local_threshold_exact_rate_by_active_pct.png",
            "local_threshold_seed_dependence.png",
            "local_threshold_failure_diagnostics.png",
            "local_threshold_runtime_by_N_K.png",
        ):
            output_path = args.out_dir / output_name
            if output_path.exists():
                print(f"  {output_path}", flush=True)
        docs_output = Path("docs/local_threshold_exact_gauss_report.md")
        if docs_output.exists():
            print(f"  {docs_output}", flush=True)
        return

    if args.threshold_exact_study:
        if (
            args.threshold_scaling_study
            or args.threshold_full_sweep
            or args.threshold_explore
            or args.threshold_rule_cdf
        ):
            raise ValueError(
                "Use --threshold-exact-study by itself, not with other threshold modes."
            )
        if args.plot_only:
            exact_path = args.out_dir / "exact_runs.csv"
            formula_path = args.out_dir / "exact_formula_runs.csv"
            if not exact_path.exists():
                raise FileNotFoundError(f"No existing exact-study file: {exact_path}")
            if not formula_path.exists():
                raise FileNotFoundError(f"No existing exact formula file: {formula_path}")
            exact_runs = pd.read_csv(exact_path)
            formula_runs = pd.read_csv(formula_path)
            summary = build_threshold_exact_summary(exact_runs)
            formula_summary = build_threshold_exact_formula_summary(formula_runs)
            atomic_write_csv(summary, args.out_dir / "exact_summary.csv")
            atomic_write_csv(formula_summary, args.out_dir / "exact_formula_summary.csv")
            write_threshold_exact_plots(
                exact_runs,
                formula_runs,
                formula_summary,
                args.out_dir,
            )
            write_threshold_exact_report(
                exact_runs,
                summary,
                formula_summary,
                args.out_dir,
                args,
            )
            write_threshold_exact_k_analysis(exact_runs, formula_runs, args.out_dir)
        else:
            run_threshold_exact_study(args)
        print("Wrote:", flush=True)
        for output_name in (
            "threshold_runs.csv.gz",
            "exact_runs.csv",
            "exact_formula_runs.csv",
            "exact_summary.csv",
            "exact_formula_summary.csv",
            "threshold_exact_report.md",
            "exact_fraction_cdf.png",
            "exact_raw_u_g_cdf.png",
            "exact_best_T_boxplot.png",
            "exact_runtime_by_N.png",
            "exact_formula_fraction_by_N.png",
            "threshold_exact_k_pct_analysis.md",
            "exact_k_pct_raw_u_g_cdf.png",
            "exact_k_pct_fraction_exact_cdf.png",
            "exact_k_pct_fraction_by_active_pct.png",
            "exact_k_pct_best_T_dependence.png",
            "exact_k_pct_rule_summary.csv",
            "exact_k_pct_best_t_summary.csv",
        ):
            output_path = args.out_dir / output_name
            if output_path.exists():
                print(f"  {output_path}", flush=True)
        return

    if args.threshold_scaling_study:
        if args.threshold_full_sweep or args.threshold_explore or args.threshold_rule_cdf:
            raise ValueError(
                "Use --threshold-scaling-study by itself, not with other threshold modes."
            )
        if args.plot_only:
            per_shard_outputs, formula_selected_parts = load_threshold_scaling_shard_outputs(
                args.out_dir
            )
            write_threshold_scaling_combined_outputs(
                per_shard_outputs,
                formula_selected_parts,
                args,
            )
        else:
            run_threshold_scaling_study(args)
        print("Wrote:", flush=True)
        for output_name in (
            "all_best_thresholds.csv",
            "all_threshold_best_t_stats.csv",
            "all_distribution_comparison.csv",
            "all_formula_comparison.csv",
            "all_scaling_formula_summary.csv",
            "all_scaling_metric_correlations.csv",
            "all_formula_selected_runs.csv",
            "all_strong_weak_runs.csv",
            "all_strong_weak_summary.csv",
            "all_scaling_report.md",
            "scaling_cdf_raw_u_g.png",
            "scaling_cdf_fraction_best_u_g.png",
            "scaling_best_T_over_N.png",
            "scaling_best_T_over_K.png",
            "scaling_formula_fraction_by_N_L.png",
            "scaling_formula_fraction_by_profile.png",
        ):
            output_path = args.out_dir / output_name
            if output_path.exists():
                print(f"  {output_path}", flush=True)
        return

    if args.threshold_rule_cdf:
        if args.threshold_full_sweep or args.threshold_explore:
            raise ValueError(
                "Use --threshold-rule-cdf by itself on an existing full-sweep result directory."
            )
        run_threshold_rule_cdf_report(args)
        print("Wrote:", flush=True)
        for output_name in (
            "threshold_rule_comparison.csv",
            "threshold_rule_summary.csv",
            "threshold_rule_comparison_report.md",
            "threshold_rule_cdf_fraction.png",
            "threshold_rule_cdf_u_g.png",
        ):
            output_path = args.out_dir / output_name
            if output_path.exists():
                print(f"  {output_path}", flush=True)
        for profile in args.data_profiles:
            for output_name in (
                "threshold_rule_cdf_fraction.png",
                "threshold_rule_cdf_u_g.png",
            ):
                output_path = args.out_dir / profile / output_name
                if output_path.exists():
                    print(f"  {output_path}", flush=True)
        return

    if args.threshold_full_sweep:
        if args.threshold_explore:
            raise ValueError("Use only one of --threshold-full-sweep or --threshold-explore.")
        if args.plot_only:
            per_profile_outputs = []
            all_runs = []
            for profile in args.data_profiles:
                profile_dir = args.out_dir / profile
                runs_path = profile_dir / "threshold_runs.csv"
                if not runs_path.exists():
                    raise FileNotFoundError(f"No existing full-sweep runs file: {runs_path}")
                runs = pd.read_csv(runs_path)
                runs = add_full_sweep_relative_metrics(runs)
                atomic_write_csv(runs, runs_path)
                per_profile_outputs.append(
                    write_full_sweep_profile_outputs(runs, profile_dir, args)
                )
                all_runs.append(runs)
            combined_runs = pd.concat(all_runs, ignore_index=True)
            write_full_sweep_combined_outputs(combined_runs, per_profile_outputs, args)
        else:
            run_threshold_full_sweep(args)
        print("Wrote:", flush=True)
        for output_name in (
            "all_best_thresholds.csv",
            "all_threshold_best_t_stats.csv",
            "all_distribution_comparison.csv",
            "all_formula_comparison.csv",
            "all_experiments_report.md",
            "all_best_T_boxplot.png",
            "all_best_fixed_T.png",
            "all_formula_fraction.png",
        ):
            output_path = args.out_dir / output_name
            if output_path.exists():
                print(f"  {output_path}", flush=True)
        for profile in args.data_profiles:
            profile_dir = args.out_dir / profile
            runs_path = profile_dir / "threshold_runs.csv"
            if runs_path.exists():
                print(f"  {runs_path}", flush=True)
        return

    if args.threshold_explore:
        runs_path = args.out_dir / "threshold_runs.csv"
        if args.plot_only:
            if not runs_path.exists():
                raise FileNotFoundError(f"No existing runs file: {runs_path}")
            runs = pd.read_csv(runs_path)
            if "u_g_vs_best_T" not in runs:
                runs = add_threshold_relative_metrics(runs)
                atomic_write_csv(runs, runs_path)
        else:
            runs = run_threshold_exploration(args, runs_path)
        write_threshold_exploration_outputs(runs, args.out_dir, args)
        print("Wrote:", flush=True)
        for output_name in (
            "threshold_runs.csv",
            "threshold_summary.csv",
            "threshold_metric_correlations.csv",
            "threshold_exploration_report.md",
            "threshold_cdf_u_g.png",
            "threshold_cdf_u_g_db.png",
            "threshold_cdf_u_g_vs_best.png",
        ):
            output_path = args.out_dir / output_name
            if output_path.exists():
                print(f"  {output_path}", flush=True)
        return

    algorithms = select_algorithms(CDF_SOLVERS, args)

    runs_path = args.out_dir / "cdf_runs.csv"
    if args.plot_only:
        if not runs_path.exists():
            raise FileNotFoundError(f"No existing runs file: {runs_path}")
        runs = pd.read_csv(runs_path)
    else:
        runs = run_benchmark(args, algorithms, runs_path)
    write_outputs(runs, algorithms, args.out_dir)
    print("Wrote:", flush=True)
    for output_name in (
        "cdf_runs.csv",
        "cdf_summary.csv",
        "cdf_baseline_improvement.csv",
        "cdf_baseline_improvement.md",
        "cdf_our_vs_h123.csv",
        "cdf_our_vs_h123.md",
        "cdf_u_g_db.png",
        "cdf_runtime_seconds.png",
        "cdf_u_g_db_h3_submodular_gen.png",
        "cdf_runtime_seconds_h3_submodular_gen.png",
    ):
        output_path = args.out_dir / output_name
        if output_path.exists():
            print(f"  {output_path}", flush=True)


if __name__ == "__main__":
    main()
