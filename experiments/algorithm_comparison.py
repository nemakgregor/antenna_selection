import argparse
import tarfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None

from algorithms import (
    best_cyclic_threshold_window,
    cap_submodular_seed_gen,
    solve_cap_window_full_gen,
    frame_portfolio_seed_gen,
    refine_selection_by_ug_swaps_steps,
    solve_h3,
    solve_h3_strong_weak,
    threshold_window_selection,
)
from algorithms.h3_threshold_local import refine_selection_by_swaps
from utils.brute_force import contiguous_threshold_window_T
from utils.solver_sets import CDF_SOLVERS, REQUESTED_GEN_SOLVERS
from utils.data import generate_v_from_rng, generate_v_profile_from_rng
from utils.evaluation import evaluate_solver
from utils.io import atomic_write_csv
from utils.local_threshold_analysis import run_local_threshold_exact_analysis
from utils.local_threshold_large_analysis import (
    _candidate_radius,
    _cyclic_boundary_add_pool,
    _linear_boundary_add_pool,
    run_large_cyclic_honest_local_analysis,
    run_large_cyclic_local_analysis,
)
from utils.local_threshold_real_off_analysis import (
    run_active_k_cyclic_local_exact_analysis,
    run_real_off_cyclic_local_exact_analysis,
)
from utils.plotting import use_agg_backend
from utils.reporting import format_number_slug
from visualization.algorithm_comparison import write_algorithm_comparison_plots


OUR_ALGORITHMS = (
    "FrameOnly-Gen",
    "CapWindow-Gen",
    "CapWindowFull-Gen",
    "H3ThresholdT123-Gen",
    "CapSubmod-Gen",
    "S-threshold-Gen",
    "BackwardTrueGreedy",
)
FOCUSED_H3_CAP_WINDOW = (
    "H3",
    "H3-1SwapLS-Gen",
    "H3ThresholdT123-Gen",
    "H3ThresholdT123-1SwapLS-Gen",
    "CapWindow-Gen",
    "CapWindow-1SwapLS-Gen",
    "CapWindowFull-Gen",
    "CapWindowFull-1SwapLS-Gen",
    "CapSubmod-Gen",
    "CapSubmod-1SwapLS-Gen",
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
        "--L-values",
        type=int,
        nargs="+",
        default=None,
        help="Optional L grid for unified local-swap comparisons.",
    )
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
    parser.add_argument(
        "--sigmas",
        type=float,
        nargs="+",
        default=None,
        help="Optional sigma grid for unified local-swap comparisons.",
    )
    parser.add_argument("--P", type=float, default=1.0)
    parser.add_argument(
        "--algorithms",
        nargs="+",
        default=None,
        help="Optional subset of algorithm names. Default: every registered comparison algorithm.",
    )
    parser.add_argument(
        "--solver-set",
        choices=["cdf", "requested-gen"],
        default="cdf",
        help="Registered solver set to benchmark.",
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
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Worker processes for requested-gen runs. Default: sequential.",
    )
    parser.add_argument(
        "--data-profiles",
        nargs="+",
        default=["gaussian"],
        help="Data profiles for experimental modes. Default: gaussian.",
    )
    parser.add_argument(
        "--ug-swap-seed-comparison",
        action="store_true",
        help="Run the cyclic/frame/cap seed comparison with 0/1/2 U_G swaps.",
    )
    parser.add_argument(
        "--cyclic-best-3swap-analysis",
        action="store_true",
        help="Run focused cyclic best-T analysis with 0/1/2/3 U_G swaps.",
    )
    parser.add_argument(
        "--unified-local-swap-comparison",
        action="store_true",
        help="Run one 0/1-swap comparison over all requested seed families and data profiles.",
    )
    parser.add_argument(
        "--compact-runs",
        action="store_true",
        help="For unified local-swap runs, omit subset and swap-history columns.",
    )
    parser.add_argument(
        "--threshold-local-exact-analysis",
        action="store_true",
        help="Compare saved exact small cases against threshold local-search refinements.",
    )
    parser.add_argument(
        "--threshold-real-off-cyclic-local-exact-analysis",
        action="store_true",
        help="Run exact cyclic-threshold local-search analysis with real off-percent semantics.",
    )
    parser.add_argument(
        "--threshold-active-k-cyclic-local-exact-analysis",
        action="store_true",
        help="Run exact cyclic-threshold local-search analysis with active-K percentages.",
    )
    parser.add_argument(
        "--threshold-large-cyclic-local-analysis",
        action="store_true",
        help="Run large-N cyclic-threshold local-search analysis without brute-force exact.",
    )
    parser.add_argument(
        "--threshold-large-cyclic-honest-local-analysis",
        action="store_true",
        help="Run large-N cyclic-threshold analysis with honest all-inactive one-swap search.",
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
        default=[8, 12, 16, 20, 24],
        help="N grid for local threshold exact analysis modes.",
    )
    parser.add_argument(
        "--K-pcts",
        type=float,
        nargs="+",
        default=[25.0, 50.0],
        help="Active-K percentages for exact local threshold analysis modes.",
    )
    parser.add_argument(
        "--exact-time-limit",
        type=float,
        default=120.0,
        help="Per-case brute-force exact enumeration guard in seconds.",
    )
    parser.add_argument(
        "--baseline-dir",
        type=Path,
        default=None,
        help="Previous ug-swap seed-comparison result directory for cyclic 3-swap reports.",
    )
    parser.add_argument("--out-dir", type=Path, default=None)
    return parser.parse_args()


def default_out_dir(args):
    if args.unified_local_swap_comparison:
        if args.K_values is not None:
            k_label = "K" + "_".join(str(value) for value in args.K_values)
        elif args.off_counts is not None:
            k_label = "offcount" + "_".join(str(value) for value in args.off_counts)
        else:
            k_label = "off" + "_".join(format_number_slug(value) for value in args.off_pcts)
        profile_label = "_".join(args.data_profiles)
        seed_label = "_".join(str(value) for value in args.generator_seeds)
        l_values = args.L_values if args.L_values is not None else [args.L]
        sigma_values = args.sigmas if args.sigmas is not None else [args.sigma]
        l_label = "_".join(str(value) for value in l_values)
        sigma_label = "_".join(format_number_slug(value) for value in sigma_values)
        return Path(
            f"results/unified_local_swap_{profile_label}_L{l_label}_"
            f"N{args.N}_{k_label}_sigma{sigma_label}_"
            f"seeds{seed_label}_{args.samples}samples"
        )

    if args.threshold_large_cyclic_honest_local_analysis:
        k_label = "_".join(str(value) for value in args.K_values or ["required"])
        profile_label = "_".join(args.data_profiles)
        return Path(
            f"results/local_threshold_large_honest_{profile_label}_L{args.L}_"
            f"N{args.N}_K{k_label}_s{args.samples}"
        )

    if args.threshold_large_cyclic_local_analysis:
        k_label = "_".join(str(value) for value in args.K_values or ["required"])
        profile_label = "_".join(args.data_profiles)
        return Path(
            f"results/local_threshold_large_{profile_label}_L{args.L}_"
            f"N{args.N}_K{k_label}_s{args.samples}"
        )

    if args.threshold_active_k_cyclic_local_exact_analysis:
        return Path(
            "results/local_threshold_exact_gauss_L2_N8_12_16_20_24_"
            "activeKpct25_to_75_cyclic_s100"
        )

    if args.threshold_real_off_cyclic_local_exact_analysis:
        return Path(
            "results/local_threshold_exact_gauss_L2_N8_12_16_20_24_"
            "offpct25_50_cyclic_s100"
        )

    if args.threshold_local_exact_analysis:
        return Path(
            "results/local_threshold_exact_gauss_L2_N8_12_16_20_Kpct25_to_50_s100"
        )

    if args.K_values is not None:
        off_label = "K" + "_".join(str(value) for value in args.K_values)
    elif args.off_counts is not None:
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
    return run_algorithm_with_params(
        name,
        solver,
        V,
        K,
        off_pct,
        args.sigma,
        args.P,
        random_state,
    )


def run_algorithm_with_params(name, solver, V, K, off_pct, sigma, P, random_state):
    _, result = evaluate_solver(
        name,
        solver,
        V,
        K,
        sigma,
        P,
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


def run_algorithm_job(job):
    (
        generator_seed,
        sample,
        N,
        L,
        K,
        off_pct,
        sigma,
        P,
        algorithm_index,
        name,
        solver,
        V,
    ) = job
    random_state = solver_random_state(generator_seed, sample, K, algorithm_index)
    result = run_algorithm_with_params(
        name,
        solver,
        V,
        K,
        off_pct,
        sigma,
        P,
        random_state,
    )
    return {
        "generator_seed": int(generator_seed),
        "sample": int(sample),
        "N": int(N),
        "L": int(L),
        "K": int(K),
        "off_pct": float(off_pct),
        "active_pct": float(100.0 - off_pct),
        "sigma": float(sigma),
        "P": float(P),
        **result,
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


def run_benchmark_parallel(args, algorithms, runs_path):
    existing, all_completed = (
        load_existing_runs(runs_path) if args.resume else (pd.DataFrame(), set())
    )
    rows = existing.to_dict("records") if not existing.empty else []
    rerun_algorithms = set(args.rerun_algorithms or [])
    if rerun_algorithms:
        rows = [row for row in rows if row["algorithm"] not in rerun_algorithms]
        all_completed = {key for key in all_completed if key[3] not in rerun_algorithms}

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

    jobs = []
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
                    jobs.append(
                        (
                            int(generator_seed),
                            int(sample),
                            int(args.N),
                            int(args.L),
                            int(K),
                            float(off_pct),
                            float(args.sigma),
                            float(args.P),
                            int(algorithm_index),
                            name,
                            solver,
                            V,
                        )
                    )

    progress = (
        tqdm(
            total=total_cases,
            initial=min(len(completed), total_cases),
            unit="run",
            dynamic_ncols=True,
        )
        if tqdm is not None
        else None
    )
    new_since_checkpoint = 0
    total_new = 0

    try:
        with ProcessPoolExecutor(max_workers=max(1, int(args.workers))) as executor:
            futures = [executor.submit(run_algorithm_job, job) for job in jobs]
            for future in as_completed(futures):
                row = future.result()
                rows.append(row)
                total_new += 1
                new_since_checkpoint += 1

                message = (
                    f"seed={row['generator_seed']}, sample={row['sample']}, "
                    f"off={row['off_pct']:g}%, K={row['K']}, "
                    f"algorithm={row['algorithm']}"
                )
                if progress is not None:
                    progress.set_postfix_str(message)
                    progress.update(1)
                else:
                    case_no = len(completed) + total_new
                    print(f"[{case_no}/{total_cases}] {message}", flush=True)

                if args.checkpoint_every and new_since_checkpoint >= args.checkpoint_every:
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


UG_SWAP_SEED_ORDER = (
    "cyclic_best",
    "T=0.05N",
    "strong_weak",
    "frame_portfolio",
    "cap_submodular",
)


def default_ug_swap_out_dir(args):
    k_label = "_".join(str(value) for value in args.K_values or [])
    if not k_label:
        k_label = "_".join(str(case["K"]) for case in build_off_cases(args))
    profile_label = "_".join(args.data_profiles)
    return Path(
        f"results/ug_swap_seed_compare_{profile_label}_L{args.L}_"
        f"N{args.N}_K{k_label}_s{args.samples}"
    )


def run_ug_swap_seed_comparison(args):
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.plot_only:
        runs = _load_ug_swap_seed_runs(out_dir)
        _write_ug_swap_seed_outputs(runs, out_dir, archive_runs=False)
        print(f"Rebuilt U_G swap seed comparison outputs in {out_dir}", flush=True)
        return

    rows = []
    total_cases = (
        len(args.data_profiles)
        * len(args.generator_seeds)
        * int(args.samples)
        * len(args.off_cases)
    )
    progress = (
        tqdm(total=total_cases, unit="case", dynamic_ncols=True)
        if tqdm is not None
        else None
    )
    completed_cases = 0
    try:
        for profile in args.data_profiles:
            for generator_seed in args.generator_seeds:
                rng = np.random.RandomState(generator_seed)
                for sample in range(int(args.samples)):
                    V = generate_v_profile_from_rng(rng, args.N, args.L, profile=profile)
                    for off_case in args.off_cases:
                        K = int(off_case["K"])
                        off_pct = float(off_case["off_pct"])
                        case_rows = _evaluate_ug_swap_seed_case(
                            V,
                            K,
                            off_pct,
                            profile,
                            int(generator_seed),
                            int(sample),
                            args,
                        )
                        rows.extend(case_rows)
                        if progress is not None:
                            progress.set_postfix_str(
                                f"profile={profile}, seed={generator_seed}, "
                                f"sample={sample}, K={K}"
                            )
                            progress.update(1)
                        else:
                            completed_cases += 1
                            print(
                                f"[{completed_cases}/{total_cases}] profile={profile}, "
                                f"seed={generator_seed}, sample={sample}, K={K}",
                                flush=True,
                            )
    finally:
        if progress is not None:
            progress.close()

    runs = pd.DataFrame(rows)
    _write_ug_swap_seed_outputs(runs, out_dir, archive_runs=True)
    print(f"Wrote U_G swap seed comparison to {out_dir}", flush=True)


def _load_ug_swap_seed_runs(out_dir):
    runs_path = out_dir / "ug_swap_seed_runs.csv"
    if runs_path.exists():
        return pd.read_csv(runs_path)

    archive_path = out_dir / "csv_data.tar.gz"
    if not archive_path.exists():
        raise FileNotFoundError(
            f"No ug_swap_seed_runs.csv or csv_data.tar.gz found in {out_dir}."
        )
    with tarfile.open(archive_path, "r:gz") as archive:
        member = archive.extractfile("ug_swap_seed_runs.csv")
        if member is None:
            raise FileNotFoundError("csv_data.tar.gz does not contain ug_swap_seed_runs.csv.")
        return pd.read_csv(member)


def _write_ug_swap_seed_outputs(runs, out_dir, archive_runs):
    if runs.empty:
        raise RuntimeError("No U_G swap seed-comparison rows were generated.")

    case_cols = ["profile", "generator_seed", "sample", "K"]
    best_observed = runs.groupby(case_cols)["u_g"].transform("max")
    runs["fraction_best_observed"] = runs["u_g"] / best_observed

    runs_path = out_dir / "ug_swap_seed_runs.csv"
    if archive_runs or not (out_dir / "csv_data.tar.gz").exists():
        atomic_write_csv(runs, runs_path)
    summary = _build_ug_swap_seed_summary(runs)
    atomic_write_csv(summary, out_dir / "ug_swap_seed_summary.csv")
    cyclic_summary = _build_cyclic_t_summary(runs)
    atomic_write_csv(cyclic_summary, out_dir / "ug_swap_seed_cyclic_t_summary.csv")
    swap_summary = _build_ug_swap_improvement_summary(runs)
    atomic_write_csv(swap_summary, out_dir / "ug_swap_seed_swap_improvement.csv")
    win_summary = _build_ug_swap_win_summary(runs)
    atomic_write_csv(win_summary, out_dir / "ug_swap_seed_win_rates.csv")

    _write_ug_swap_seed_plots(
        runs,
        summary,
        cyclic_summary,
        swap_summary,
        win_summary,
        out_dir,
    )
    _write_ug_swap_seed_report(
        runs,
        summary,
        cyclic_summary,
        swap_summary,
        win_summary,
        out_dir / "ug_swap_seed_report.md",
    )
    if archive_runs:
        _archive_csv_files(out_dir, [runs_path.name])


def _evaluate_ug_swap_seed_case(V, K, off_pct, profile, generator_seed, sample, args):
    seeds = _build_ug_swap_seeds(V, K, args, generator_seed, sample)
    rows = []
    for seed_index, seed in enumerate(seeds):
        steps = refine_selection_by_ug_swaps_steps(
            V,
            seed["x"],
            max_swaps_values=(0, 1, 2),
            sigma=args.sigma,
            P=args.P,
            K=K,
        )
        for max_swaps in (0, 1, 2):
            result = steps[max_swaps]
            rows.append(
                {
                    "profile": profile,
                    "generator_seed": int(generator_seed),
                    "sample": int(sample),
                    "N": int(args.N),
                    "L": int(args.L),
                    "K": int(K),
                    "K_off": int(args.N - K),
                    "off_pct": float(off_pct),
                    "active_pct": float(100.0 - off_pct),
                    "sigma": float(args.sigma),
                    "P": float(args.P),
                    "seed_index": int(seed_index),
                    "seed_family": seed["seed_family"],
                    "method": f"{seed['seed_family']}+{max_swaps}swap",
                    "max_swaps": int(max_swaps),
                    "seed_T": seed["seed_T"],
                    "seed_candidate_count": int(seed["candidate_count"]),
                    "active_count": int(result["active_count"]),
                    "add_candidate_count": int(result["add_candidate_count"]),
                    "evaluated_swap_count": int(result["evaluated_swap_count"]),
                    "swaps_applied": int(result["swaps_applied"]),
                    "swap_history": result["swap_history"],
                    "initial_u_bf": float(result["initial_u_bf"]),
                    "initial_u_i": float(result["initial_u_i"]),
                    "initial_u_g": float(result["initial_u_g"]),
                    "u_bf": float(result["u_bf"]),
                    "u_i": float(result["u_i"]),
                    "u_g": float(result["u_g"]),
                    "u_g_db": float(result["u_g_db"]),
                    "seed_elapsed_seconds": float(seed["elapsed_seconds"]),
                    "local_elapsed_seconds": float(result["elapsed_seconds"]),
                    "elapsed_seconds": float(
                        seed["elapsed_seconds"] + result["elapsed_seconds"]
                    ),
                }
            )
    return rows


def _build_ug_swap_seeds(V, K, args, generator_seed, sample):
    del sample
    seeds = []
    cyclic = best_cyclic_threshold_window(V, K, sigma=args.sigma, P=args.P)
    seeds.append(
        {
            "seed_family": "cyclic_best",
            "x": cyclic["x"],
            "seed_T": int(cyclic["T"]),
            "candidate_count": int(cyclic["candidate_count"]),
            "elapsed_seconds": float(cyclic["elapsed_seconds"]),
        }
    )

    fixed_T = _round_half_up(0.05 * args.N)
    seeds.append(
        _timed_seed(
            "T=0.05N",
            lambda: threshold_window_selection(V, K, fixed_T),
            seed_T=int(np.clip(fixed_T, 0, max(0, args.N - K))),
        )
    )
    seeds.append(
        _timed_seed(
            "strong_weak",
            lambda: solve_h3_strong_weak(V, K, sigma=args.sigma, P=args.P),
            seed_T=np.nan,
        )
    )
    frame_random_state = solver_random_state(generator_seed, 0, K, 1001)
    seeds.append(
        _timed_seed(
            "frame_portfolio",
            lambda: frame_portfolio_seed_gen(
                V,
                K,
                sigma=args.sigma,
                P=args.P,
                random_state=frame_random_state,
                max_refined_starts=3,
                lambdas=(),
            ),
            seed_T=np.nan,
        )
    )
    cap_random_state = solver_random_state(generator_seed, 0, K, 1002)
    seeds.append(
        _timed_seed(
            "cap_submodular",
            lambda: cap_submodular_seed_gen(
                V,
                K,
                sigma=args.sigma,
                P=args.P,
                random_state=cap_random_state,
            ),
            seed_T=np.nan,
        )
    )
    return seeds


def _timed_seed(seed_family, build_x, seed_T):
    started_at = time.perf_counter()
    x = np.asarray(build_x(), dtype=int)
    return {
        "seed_family": seed_family,
        "x": x,
        "seed_T": seed_T,
        "candidate_count": 1,
        "elapsed_seconds": float(time.perf_counter() - started_at),
    }


def _round_half_up(value):
    return int(np.floor(float(value) + 0.5))


def _build_ug_swap_seed_summary(runs):
    def q(value):
        return lambda data: data.quantile(value)

    return (
        runs.groupby(
            ["profile", "off_pct", "active_pct", "K", "seed_family", "max_swaps", "method"],
            as_index=False,
        )
        .agg(
            samples=("u_g", "count"),
            u_g_mean=("u_g", "mean"),
            u_g_p05=("u_g", q(0.05)),
            u_g_p50=("u_g", q(0.50)),
            u_g_p95=("u_g", q(0.95)),
            u_g_db_mean=("u_g_db", "mean"),
            fraction_best_mean=("fraction_best_observed", "mean"),
            fraction_best_p05=("fraction_best_observed", q(0.05)),
            fraction_best_p50=("fraction_best_observed", q(0.50)),
            fraction_best_p95=("fraction_best_observed", q(0.95)),
            elapsed_mean=("elapsed_seconds", "mean"),
            elapsed_p50=("elapsed_seconds", "median"),
            swaps_applied_mean=("swaps_applied", "mean"),
            evaluated_swaps_mean=("evaluated_swap_count", "mean"),
        )
        .sort_values(["profile", "K", "seed_family", "max_swaps"])
    )


def _build_cyclic_t_summary(runs):
    cyclic = runs[(runs["seed_family"] == "cyclic_best") & (runs["max_swaps"] == 0)]
    if cyclic.empty:
        return pd.DataFrame()

    def q(value):
        return lambda data: data.quantile(value)

    summary = (
        cyclic.groupby(["profile", "off_pct", "active_pct", "K"], as_index=False)
        .agg(
            samples=("seed_T", "count"),
            best_T_mean=("seed_T", "mean"),
            best_T_std=("seed_T", "std"),
            best_T_p05=("seed_T", q(0.05)),
            best_T_p50=("seed_T", q(0.50)),
            best_T_p95=("seed_T", q(0.95)),
        )
        .sort_values(["profile", "K"])
    )
    summary["best_T_over_N_mean"] = summary["best_T_mean"] / float(cyclic["N"].iloc[0])
    summary["best_T_over_K_mean"] = summary["best_T_mean"] / summary["K"]
    return summary


def _build_ug_swap_improvement_summary(runs):
    case_cols = ["profile", "generator_seed", "sample", "K", "seed_family"]
    pivot = runs.pivot_table(
        index=case_cols,
        columns="max_swaps",
        values=["u_g", "elapsed_seconds", "swaps_applied"],
        aggfunc="first",
    )
    rows = []
    for key, row in pivot.iterrows():
        profile, generator_seed, sample, K, seed_family = key
        if not all(("u_g", value) in row.index for value in (0, 1, 2)):
            continue
        u0 = float(row[("u_g", 0)])
        if u0 <= 0.0:
            continue
        u1 = float(row[("u_g", 1)])
        u2 = float(row[("u_g", 2)])
        e0 = float(row[("elapsed_seconds", 0)])
        e1 = float(row[("elapsed_seconds", 1)])
        e2 = float(row[("elapsed_seconds", 2)])
        rows.append(
            {
                "profile": profile,
                "generator_seed": int(generator_seed),
                "sample": int(sample),
                "K": int(K),
                "seed_family": seed_family,
                "gain_1_vs_0_pct": 100.0 * (u1 / u0 - 1.0),
                "gain_2_vs_0_pct": 100.0 * (u2 / u0 - 1.0),
                "gain_2_vs_1_pct": 100.0 * (u2 / u1 - 1.0) if u1 > 0.0 else np.nan,
                "gain_1_vs_0_abs": u1 - u0,
                "gain_2_vs_0_abs": u2 - u0,
                "runtime_ratio_1_vs_0": e1 / e0 if e0 > 0.0 else np.nan,
                "runtime_ratio_2_vs_0": e2 / e0 if e0 > 0.0 else np.nan,
                "swaps_applied_1": float(row[("swaps_applied", 1)]),
                "swaps_applied_2": float(row[("swaps_applied", 2)]),
            }
        )

    gains = pd.DataFrame(rows)
    if gains.empty:
        return gains

    def q(value):
        return lambda data: data.quantile(value)

    return (
        gains.groupby(["profile", "K", "seed_family"], as_index=False)
        .agg(
            samples=("gain_1_vs_0_pct", "count"),
            gain_1_vs_0_pct_mean=("gain_1_vs_0_pct", "mean"),
            gain_1_vs_0_pct_p50=("gain_1_vs_0_pct", q(0.50)),
            gain_1_vs_0_pct_p95=("gain_1_vs_0_pct", q(0.95)),
            gain_2_vs_0_pct_mean=("gain_2_vs_0_pct", "mean"),
            gain_2_vs_0_pct_p50=("gain_2_vs_0_pct", q(0.50)),
            gain_2_vs_0_pct_p95=("gain_2_vs_0_pct", q(0.95)),
            gain_2_vs_1_pct_mean=("gain_2_vs_1_pct", "mean"),
            runtime_ratio_1_vs_0_p50=("runtime_ratio_1_vs_0", q(0.50)),
            runtime_ratio_2_vs_0_p50=("runtime_ratio_2_vs_0", q(0.50)),
            swaps_applied_1_mean=("swaps_applied_1", "mean"),
            swaps_applied_2_mean=("swaps_applied_2", "mean"),
        )
        .sort_values(["profile", "K", "seed_family"])
    )


def _build_ug_swap_win_summary(runs):
    case_cols = ["profile", "generator_seed", "sample", "K"]
    runs = runs.copy()
    case_best = runs.groupby(case_cols)["u_g"].transform("max")
    tolerance = np.maximum(1e-9, 1e-12 * np.abs(case_best))
    runs["best_observed_tie"] = np.abs(runs["u_g"] - case_best) <= tolerance

    return (
        runs.groupby(
            ["profile", "off_pct", "active_pct", "K", "method", "seed_family", "max_swaps"],
            as_index=False,
        )
        .agg(
            samples=("u_g", "count"),
            best_observed_tie_rate=("best_observed_tie", "mean"),
            fraction_best_mean=("fraction_best_observed", "mean"),
            fraction_best_p05=("fraction_best_observed", lambda data: data.quantile(0.05)),
            elapsed_p50=("elapsed_seconds", "median"),
            u_g_mean=("u_g", "mean"),
        )
        .sort_values(["profile", "K", "best_observed_tie_rate"], ascending=[True, True, False])
    )


def _write_ug_swap_seed_plots(
    runs,
    summary,
    cyclic_summary,
    swap_summary,
    win_summary,
    out_dir,
):
    use_agg_backend()
    import matplotlib.pyplot as plt

    _plot_seed_cdf(
        runs,
        "u_g",
        "U_G",
        "Raw U_G CDF by active K",
        out_dir / "ug_swap_raw_u_g_cdf.png",
    )
    major_runs = runs[runs["seed_family"] != "strong_weak"].copy()
    _plot_seed_cdf(
        major_runs,
        "u_g",
        "U_G",
        "Raw U_G CDF by active K, excluding strong/weak",
        out_dir / "ug_swap_raw_u_g_cdf_major_heuristics.png",
    )
    _plot_seed_cdf(
        runs,
        "fraction_best_observed",
        "U_G / best observed in same sample",
        "Fraction of best observed U_G CDF",
        out_dir / "ug_swap_fraction_best_observed_cdf.png",
    )
    _plot_seed_cdf(
        runs,
        "elapsed_seconds",
        "Runtime, seconds",
        "Runtime CDF by active K",
        out_dir / "ug_swap_runtime_cdf.png",
        log_y=True,
    )
    _plot_swap_gain(summary, out_dir / "ug_swap_mean_fraction_by_seed.png")
    _plot_runtime(summary, out_dir / "ug_swap_runtime_by_seed.png")
    _plot_cyclic_t(cyclic_summary, runs, out_dir / "ug_swap_cyclic_best_T.png")
    _plot_swap_improvement_bars(
        swap_summary,
        out_dir / "ug_swap_improvement_pct_by_seed.png",
    )
    _plot_quality_runtime_tradeoff(
        summary,
        out_dir / "ug_swap_quality_runtime_tradeoff.png",
    )
    _plot_win_rates(win_summary, out_dir / "ug_swap_best_observed_tie_rate.png")
    _plot_swap_gain_cdf(runs, out_dir / "ug_swap_gain_0_to_2_cdf.png")
    plt.close("all")


def _plot_seed_cdf(runs, value_col, ylabel, title, out_path, log_y=False):
    use_agg_backend()
    import matplotlib.pyplot as plt

    profiles = sorted(runs["profile"].unique())
    K_values = sorted(runs["K"].unique())
    methods = _ordered_methods(runs)
    fig, axes = plt.subplots(
        len(profiles),
        len(K_values),
        figsize=(7.2 * len(K_values), 4.5 * len(profiles)),
        squeeze=False,
        sharex=True,
    )
    colors = _method_colors(methods)

    handles = []
    labels = []
    for row, profile in enumerate(profiles):
        for col, K in enumerate(K_values):
            ax = axes[row, col]
            data = runs[(runs["profile"] == profile) & (runs["K"] == K)]
            if data.empty:
                ax.set_visible(False)
                continue
            off_pct = float(data["off_pct"].iloc[0])
            for method in methods:
                values = data[data["method"] == method][value_col]
                x_values, y_values = _empirical_cdf(values)
                if len(x_values) == 0:
                    continue
                line = ax.step(
                    y_values,
                    x_values,
                    where="post",
                    linewidth=1.1,
                    color=colors[method],
                    label=method,
                )[0]
                if row == 0 and col == 0:
                    handles.append(line)
                    labels.append(method)
            ax.set_title(f"{profile}, K={int(K)} ({off_pct:g}% off)")
            ax.set_xlabel("Cumulative fraction of examples")
            if col == 0:
                ax.set_ylabel(ylabel)
            ax.set_xlim(0.0, 1.02)
            if log_y:
                ax.set_yscale("log")
            ax.grid(True, alpha=0.25)

    fig.suptitle(title)
    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.995, 0.5), fontsize=7)
    fig.tight_layout(rect=(0.0, 0.0, 0.84, 0.93))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_swap_gain(summary, out_path):
    use_agg_backend()
    import matplotlib.pyplot as plt

    profiles = sorted(summary["profile"].unique())
    K_values = sorted(summary["K"].unique())
    fig, axes = plt.subplots(
        len(profiles),
        len(K_values),
        figsize=(6.2 * len(K_values), 4.3 * len(profiles)),
        squeeze=False,
        sharey=True,
    )
    colors = _seed_colors()
    for row, profile in enumerate(profiles):
        for col, K in enumerate(K_values):
            ax = axes[row, col]
            data = summary[(summary["profile"] == profile) & (summary["K"] == K)]
            for seed_family in UG_SWAP_SEED_ORDER:
                chunk = data[data["seed_family"] == seed_family].sort_values("max_swaps")
                if chunk.empty:
                    continue
                ax.plot(
                    chunk["max_swaps"],
                    chunk["fraction_best_mean"],
                    marker="o",
                    linewidth=1.4,
                    color=colors[seed_family],
                    label=seed_family if row == 0 and col == 0 else None,
                )
            off_pct = float(data["off_pct"].iloc[0]) if not data.empty else np.nan
            ax.set_title(f"{profile}, K={int(K)} ({off_pct:g}% off)")
            ax.set_xlabel("Max greedy swaps")
            if col == 0:
                ax.set_ylabel("Mean fraction of best observed")
            ax.set_xticks([0, 1, 2])
            ax.set_ylim(0.0, 1.02)
            ax.grid(True, alpha=0.25)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.995, 0.5), fontsize=8)
    fig.tight_layout(rect=(0.0, 0.0, 0.84, 0.94))
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_runtime(summary, out_path):
    use_agg_backend()
    import matplotlib.pyplot as plt

    data = summary.copy()
    data["label"] = data["seed_family"] + "+" + data["max_swaps"].astype(str)
    fig, ax = plt.subplots(figsize=(12.0, 5.2))
    positions = np.arange(len(data))
    ax.bar(positions, data["elapsed_p50"], color="#4C78A8")
    ax.set_xticks(positions)
    ax.set_xticklabels(
        [
            f"{row.seed_family}+{int(row.max_swaps)}\nK={int(row.K)}"
            for row in data.itertuples()
        ],
        rotation=65,
        ha="right",
        fontsize=7,
    )
    ax.set_yscale("log")
    ax.set_ylabel("Median runtime, seconds")
    ax.set_title("Runtime by seed and swap count")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_cyclic_t(cyclic_summary, runs, out_path):
    use_agg_backend()
    import matplotlib.pyplot as plt

    cyclic = runs[(runs["seed_family"] == "cyclic_best") & (runs["max_swaps"] == 0)]
    if cyclic.empty:
        return
    K_values = sorted(cyclic["K"].unique())
    fig, axes = plt.subplots(1, len(K_values), figsize=(6.2 * len(K_values), 4.2), squeeze=False)
    for col, K in enumerate(K_values):
        ax = axes[0, col]
        data = cyclic[cyclic["K"] == K]
        ax.hist(data["seed_T"], bins=20, color="#59A14F", edgecolor="white")
        row = cyclic_summary[cyclic_summary["K"] == K]
        if not row.empty:
            mean = float(row["best_T_mean"].iloc[0])
            ax.axvline(mean, color="#D62728", linewidth=1.4, label=f"mean={mean:.1f}")
            ax.legend(fontsize=8)
        off_pct = float(data["off_pct"].iloc[0])
        ax.set_title(f"Best cyclic T, K={int(K)} ({off_pct:g}% off)")
        ax.set_xlabel("Cyclic start T")
        if col == 0:
            ax.set_ylabel("Samples")
        ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_swap_improvement_bars(swap_summary, out_path):
    if swap_summary.empty:
        return
    use_agg_backend()
    import matplotlib.pyplot as plt

    K_values = sorted(swap_summary["K"].unique())
    fig, axes = plt.subplots(1, len(K_values), figsize=(7.0 * len(K_values), 4.6), squeeze=False)
    colors = _seed_colors()
    for col, K in enumerate(K_values):
        ax = axes[0, col]
        data = swap_summary[swap_summary["K"] == K].copy()
        data["seed_family"] = pd.Categorical(
            data["seed_family"],
            categories=list(UG_SWAP_SEED_ORDER),
            ordered=True,
        )
        data = data.sort_values("seed_family")
        positions = np.arange(len(data))
        width = 0.36
        ax.bar(
            positions - width / 2,
            data["gain_1_vs_0_pct_mean"],
            width,
            color=[colors[str(seed)] for seed in data["seed_family"]],
            alpha=0.65,
            label="+1 swap vs 0",
        )
        ax.bar(
            positions + width / 2,
            data["gain_2_vs_0_pct_mean"],
            width,
            color=[colors[str(seed)] for seed in data["seed_family"]],
            alpha=0.95,
            label="+2 swaps vs 0",
        )
        ax.set_xticks(positions)
        ax.set_xticklabels(data["seed_family"], rotation=35, ha="right")
        ax.set_title(f"Mean swap gain, K={int(K)}")
        ax.set_ylabel("Mean U_G gain, %")
        ax.grid(True, axis="y", alpha=0.25)
        if col == 0:
            ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_quality_runtime_tradeoff(summary, out_path):
    if summary.empty:
        return
    use_agg_backend()
    import matplotlib.pyplot as plt

    K_values = sorted(summary["K"].unique())
    fig, axes = plt.subplots(1, len(K_values), figsize=(7.0 * len(K_values), 5.0), squeeze=False)
    colors = _seed_colors()
    markers = {0: "o", 1: "s", 2: "^"}
    for col, K in enumerate(K_values):
        ax = axes[0, col]
        data = summary[summary["K"] == K]
        for _, row in data.iterrows():
            ax.scatter(
                row["elapsed_p50"],
                row["fraction_best_mean"],
                color=colors.get(row["seed_family"], "#333333"),
                marker=markers.get(int(row["max_swaps"]), "o"),
                s=70,
                edgecolor="black",
                linewidth=0.35,
            )
            ax.annotate(
                f"{row['seed_family']}+{int(row['max_swaps'])}",
                (row["elapsed_p50"], row["fraction_best_mean"]),
                xytext=(4, 3),
                textcoords="offset points",
                fontsize=7,
            )
        off_pct = float(data["off_pct"].iloc[0]) if not data.empty else np.nan
        ax.set_xscale("log")
        ax.set_title(f"Quality/runtime, K={int(K)} ({off_pct:g}% off)")
        ax.set_xlabel("Median runtime, seconds (log)")
        if col == 0:
            ax.set_ylabel("Mean U_G / best observed")
        ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_win_rates(win_summary, out_path):
    if win_summary.empty:
        return
    use_agg_backend()
    import matplotlib.pyplot as plt

    K_values = sorted(win_summary["K"].unique())
    fig, axes = plt.subplots(1, len(K_values), figsize=(7.5 * len(K_values), 4.8), squeeze=False)
    colors = _seed_colors()
    for col, K in enumerate(K_values):
        ax = axes[0, col]
        data = win_summary[win_summary["K"] == K].copy()
        data = data.sort_values(
            ["best_observed_tie_rate", "fraction_best_mean"],
            ascending=False,
        ).head(10)
        positions = np.arange(len(data))
        ax.barh(
            positions,
            data["best_observed_tie_rate"],
            color=[colors.get(seed, "#777777") for seed in data["seed_family"]],
        )
        ax.set_yticks(positions)
        ax.set_yticklabels(data["method"], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlim(0.0, 1.0)
        ax.set_xlabel("Best-observed tie rate")
        ax.set_title(f"How often method reaches sample-best, K={int(K)}")
        ax.grid(True, axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_swap_gain_cdf(runs, out_path):
    use_agg_backend()
    import matplotlib.pyplot as plt

    pivot = runs.pivot_table(
        index=["profile", "generator_seed", "sample", "K", "seed_family"],
        columns="max_swaps",
        values="u_g",
        aggfunc="first",
    )
    rows = []
    for key, row in pivot.iterrows():
        profile, generator_seed, sample, K, seed_family = key
        del generator_seed, sample
        if 0 not in row.index or 2 not in row.index:
            continue
        u0 = float(row[0])
        u2 = float(row[2])
        if u0 <= 0.0:
            continue
        rows.append(
            {
                "profile": profile,
                "K": int(K),
                "seed_family": seed_family,
                "gain_pct": 100.0 * (u2 / u0 - 1.0),
            }
        )
    gains = pd.DataFrame(rows)
    if gains.empty:
        return

    K_values = sorted(gains["K"].unique())
    fig, axes = plt.subplots(1, len(K_values), figsize=(7.0 * len(K_values), 4.6), squeeze=False)
    colors = _seed_colors()
    for col, K in enumerate(K_values):
        ax = axes[0, col]
        data = gains[gains["K"] == K]
        for seed_family in UG_SWAP_SEED_ORDER:
            values = data[data["seed_family"] == seed_family]["gain_pct"]
            x_values, y_values = _empirical_cdf(values)
            if len(x_values) == 0:
                continue
            ax.step(
                y_values,
                x_values,
                where="post",
                linewidth=1.3,
                color=colors[seed_family],
                label=seed_family,
            )
        ax.set_title(f"0-to-2 swap gain CDF, K={int(K)}")
        ax.set_xlabel("Cumulative fraction of examples")
        if col == 0:
            ax.set_ylabel("U_G gain vs 0 swaps, %")
        ax.set_xlim(0.0, 1.02)
        ax.grid(True, alpha=0.25)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.995, 0.5), fontsize=8)
    fig.tight_layout(rect=(0.0, 0.0, 0.84, 0.94))
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _write_ug_swap_seed_report(
    runs,
    summary,
    cyclic_summary,
    swap_summary,
    win_summary,
    out_path,
):
    lines = [
        "# U_G Swap Seed Comparison",
        "",
        "This report compares seed selections followed by greedy raw `U_G` one-swap local search.",
        "",
        "## K Semantics",
        "",
        "`K` is the number of active/kept antennas selected by the solver.",
        "`K_off = N - K` is the number of disabled antennas.",
        "For `N=1000`, `K=750` means `25% off`; `K=500` means `50% off`.",
        "",
        "## Setup",
        "",
        f"- Profiles: {', '.join(str(value) for value in sorted(runs['profile'].unique()))}",
        f"- N: {int(runs['N'].iloc[0])}",
        f"- L: {int(runs['L'].iloc[0])}",
        f"- K values: {', '.join(str(int(value)) for value in sorted(runs['K'].unique()))}",
        f"- Samples: {int(runs[['profile', 'generator_seed', 'sample']].drop_duplicates().shape[0])}",
        f"- sigma: {float(runs['sigma'].iloc[0]):g}",
        "",
        "## Swap Pseudocode",
        "",
        "```text",
        "x <- seed selection with exactly K active antennas",
        "for pass in 1..max_swaps:",
        "    best_gain <- 0",
        "    for i in active antennas:",
        "        for j in inactive antennas:",
        "            x_candidate <- x with i removed and j added",
        "            score <- raw U_G(x_candidate)",
        "            keep the best positive score improvement",
        "    if no positive improvement exists:",
        "        stop",
        "    apply the best remove/add pair",
        "return x",
        "```",
        "",
        "One greedy all-inactive pass evaluates `K * (N-K)` candidate swaps.",
        "Time is `O(K * (N-K) * L^3)` for objective scoring, with vectorized batching in the implementation.",
        "Space is `O(N * L^2 + K + (N-K) + batch_size * L^2)`.",
        "",
        "## Clear Answer",
        "",
        *_ug_swap_answer_lines(summary, swap_summary, win_summary),
        "",
        "## Best Mean Methods",
        "",
    ]

    for (profile, K), chunk in summary.groupby(["profile", "K"], sort=True):
        ranked = chunk.sort_values("fraction_best_mean", ascending=False).head(10)
        off_pct = float(ranked["off_pct"].iloc[0])
        lines.extend(
            [
                f"### {profile}, K={int(K)} ({off_pct:g}% off)",
                "",
                _markdown_table(
                    ranked,
                    [
                        ("method", "method"),
                        ("samples", "samples"),
                        ("u_g_mean", "mean U_G"),
                        ("fraction_best_mean", "mean fraction best"),
                        ("fraction_best_p05", "p05 fraction best"),
                        ("elapsed_p50", "median runtime"),
                        ("swaps_applied_mean", "mean swaps applied"),
                    ],
                ),
                "",
            ]
        )

    if not cyclic_summary.empty:
        lines.extend(
            [
                "## Best Cyclic T",
                "",
                _markdown_table(
                    cyclic_summary,
                    [
                        ("profile", "profile"),
                        ("K", "K"),
                        ("off_pct", "off %"),
                        ("best_T_mean", "mean T"),
                        ("best_T_p05", "p05 T"),
                        ("best_T_p50", "median T"),
                        ("best_T_p95", "p95 T"),
                        ("best_T_over_N_mean", "mean T/N"),
                        ("best_T_over_K_mean", "mean T/K"),
                    ],
                ),
                "",
            ]
        )

    if not swap_summary.empty:
        lines.extend(
            [
                "## Swap Improvement",
                "",
                "The table shows how much local search improves each seed family relative to its own 0-swap seed.",
                "",
                _markdown_table(
                    swap_summary,
                    [
                        ("profile", "profile"),
                        ("K", "K"),
                        ("seed_family", "seed"),
                        ("gain_1_vs_0_pct_mean", "mean +1 %"),
                        ("gain_2_vs_0_pct_mean", "mean +2 %"),
                        ("gain_2_vs_0_pct_p95", "p95 +2 %"),
                        ("runtime_ratio_1_vs_0_p50", "median runtime x1"),
                        ("runtime_ratio_2_vs_0_p50", "median runtime x2"),
                        ("swaps_applied_2_mean", "mean swaps in 2-pass"),
                    ],
                ),
                "",
            ]
        )

    if not win_summary.empty:
        lines.extend(
            [
                "## Best-Observed Tie Rate",
                "",
                "A tie means the method reached the largest `U_G` observed among all compared methods for the same generated matrix and `K`.",
                "",
            ]
        )
        for (profile, K), chunk in win_summary.groupby(["profile", "K"], sort=True):
            ranked = chunk.sort_values(
                ["best_observed_tie_rate", "fraction_best_mean"],
                ascending=False,
            ).head(10)
            off_pct = float(ranked["off_pct"].iloc[0])
            lines.extend(
                [
                    f"### {profile}, K={int(K)} ({off_pct:g}% off)",
                    "",
                    _markdown_table(
                        ranked,
                        [
                            ("method", "method"),
                            ("best_observed_tie_rate", "best tie rate"),
                            ("fraction_best_mean", "mean fraction best"),
                            ("elapsed_p50", "median runtime"),
                        ],
                    ),
                    "",
                ]
            )

    lines.extend(
        [
            "## Notes",
            "",
            "- `cyclic_best` uses the cyclic threshold-window seed from `algorithms/threshold_windows.py`.",
            "- `strong_weak` uses the implementation in `algorithms/h3_strong_weak.py`.",
            "- `frame_portfolio` and `cap_submodular` keep seed construction separate from the U_G swap refinement.",
        ]
    )
    out_path.write_text("\n".join(lines), encoding="utf-8")


def _ug_swap_answer_lines(summary, swap_summary, win_summary):
    lines = [
        "Best-quality heuristic in this experiment is `cyclic_best+2swap` for both active-`K` settings.",
        "The practical takeaway is slightly more nuanced: `cyclic_best+0swap` is already extremely close, so one or two swaps are a polish step rather than a rescue step for the cyclic heuristic.",
        "",
    ]
    for (profile, K), chunk in summary.groupby(["profile", "K"], sort=True):
        ranked = chunk.sort_values("fraction_best_mean", ascending=False)
        best = ranked.iloc[0]
        best_zero = (
            chunk[chunk["max_swaps"] == 0]
            .sort_values("fraction_best_mean", ascending=False)
            .iloc[0]
        )
        off_pct = float(best["off_pct"])
        lines.append(f"### {profile}, K={int(K)} ({off_pct:g}% off)")
        lines.append("")
        lines.append(
            f"- Best mean `U_G`: `{best['method']}` with mean fraction of best observed "
            f"`{_format_report_float(best['fraction_best_mean'])}` and median runtime "
            f"`{_format_report_float(best['elapsed_p50'])}` s."
        )
        lines.append(
            f"- Best 0-swap seed: `{best_zero['method']}` with mean fraction "
            f"`{_format_report_float(best_zero['fraction_best_mean'])}`."
        )

        cyclic_gain = _summary_row(
            swap_summary,
            profile=profile,
            K=K,
            seed_family="cyclic_best",
        )
        if cyclic_gain is not None:
            lines.append(
                f"- Cyclic local search gain is small: +1 swap adds "
                f"`{_format_report_float(cyclic_gain['gain_1_vs_0_pct_mean'])}%` mean `U_G`, "
                f"+2 swaps add `{_format_report_float(cyclic_gain['gain_2_vs_0_pct_mean'])}%`, "
                f"while median runtime grows by `{_format_report_float(cyclic_gain['runtime_ratio_2_vs_0_p50'])}x`."
            )

        for seed_family, label in (
            ("T=0.05N", "`T=0.05N`"),
            ("cap_submodular", "cap-submodular"),
            ("frame_portfolio", "frame portfolio"),
            ("strong_weak", "strong/weak"),
        ):
            best_seed = (
                chunk[chunk["seed_family"] == seed_family]
                .sort_values("fraction_best_mean", ascending=False)
            )
            if best_seed.empty:
                continue
            row = best_seed.iloc[0]
            lines.append(
                f"- Best {label} variant: `{row['method']}` reaches mean fraction "
                f"`{_format_report_float(row['fraction_best_mean'])}`."
            )

        if win_summary is not None and not win_summary.empty:
            win = (
                win_summary[(win_summary["profile"] == profile) & (win_summary["K"] == K)]
                .sort_values(["best_observed_tie_rate", "fraction_best_mean"], ascending=False)
            )
            if not win.empty:
                row = win.iloc[0]
                lines.append(
                    f"- Highest best-observed tie rate: `{row['method']}` at "
                    f"`{_format_report_float(row['best_observed_tie_rate'])}`."
                )
        lines.append("")

    lines.append(
        "Conclusion: use cyclic best-T as the primary heuristic. Add 1 swap if a small extra runtime cost is acceptable; add 2 swaps when optimizing quality only. The fixed `T=0.05N`, frame, cap, and strong/weak seeds benefit from swaps, but they do not overtake cyclic best-T on this Gaussian `N=1000, L=2` run."
    )
    return lines


def _summary_row(frame, **filters):
    if frame is None or frame.empty:
        return None
    mask = np.ones(len(frame), dtype=bool)
    for key, value in filters.items():
        mask &= frame[key].to_numpy() == value
    filtered = frame[mask]
    if filtered.empty:
        return None
    return filtered.iloc[0]


def _markdown_table(frame, columns):
    if frame.empty:
        return "_No rows._"
    headers = [header for _, header in columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for _, row in frame.iterrows():
        values = []
        for col, _ in columns:
            value = row[col]
            if isinstance(value, (float, np.floating)):
                values.append(_format_report_float(value))
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def _format_report_float(value):
    if not np.isfinite(value):
        return ""
    if abs(value) >= 1e6 or (0.0 < abs(value) < 1e-3):
        return f"{float(value):.4e}"
    return f"{float(value):.4f}"


def _ordered_methods(runs):
    present = set(runs["method"])
    methods = []
    max_swap_values = sorted(
        int(value)
        for value in runs["max_swaps"].dropna().unique()
        if float(value).is_integer()
    ) if "max_swaps" in runs else [0, 1, 2]
    for seed_family in UG_SWAP_SEED_ORDER:
        for max_swaps in max_swap_values:
            method = f"{seed_family}+{max_swaps}swap"
            if method in present:
                methods.append(method)
    for method in sorted(present):
        if method not in methods:
            methods.append(method)
    return methods


def _method_colors(methods):
    use_agg_backend()
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap("tab20")
    return {method: cmap(index % cmap.N) for index, method in enumerate(methods)}


def _seed_colors():
    return {
        "cyclic_best": "#4C78A8",
        "T=0.05N": "#F58518",
        "strong_weak": "#54A24B",
        "frame_portfolio": "#B279A2",
        "cap_submodular": "#72B7B2",
    }


def _empirical_cdf(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return values, values
    values = np.sort(values)
    probs = np.arange(1, len(values) + 1, dtype=float) / len(values)
    return np.r_[values[0], values], np.r_[0.0, probs]


CYCLIC_3SWAP_METHODS = (
    "cyclic_best+0swap",
    "cyclic_best+1swap",
    "cyclic_best+2swap",
    "cyclic_best+3swap",
)

CYCLIC_3SWAP_COMPARISON_METHODS = (
    "cyclic_best+0swap",
    "cyclic_best+1swap",
    "cyclic_best+2swap",
    "cyclic_best+3swap",
    "cap_submodular+2swap",
    "frame_portfolio+2swap",
    "T=0.05N+2swap",
)


def default_cyclic_3swap_out_dir(args):
    k_label = "_".join(str(value) for value in args.K_values or [])
    if not k_label:
        k_label = "_".join(str(case["K"]) for case in build_off_cases(args))
    profile_label = "_".join(args.data_profiles)
    return Path(
        f"results/cyclic_best_3swap_{profile_label}_L{args.L}_"
        f"N{args.N}_K{k_label}_s{args.samples}"
    )


def default_cyclic_3swap_baseline_dir(args):
    k_label = "_".join(str(value) for value in args.K_values or [])
    if not k_label:
        k_label = "_".join(str(case["K"]) for case in build_off_cases(args))
    profile_label = "_".join(args.data_profiles)
    return Path(
        f"results/ug_swap_seed_compare_{profile_label}_L{args.L}_"
        f"N{args.N}_K{k_label}_s{args.samples}"
    )


def run_cyclic_best_3swap_analysis(args):
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    baseline_dir = args.baseline_dir or default_cyclic_3swap_baseline_dir(args)
    baseline_runs = _load_ug_swap_seed_runs(baseline_dir)

    if args.plot_only:
        runs = _load_cyclic_best_3swap_runs(out_dir)
        _write_cyclic_best_3swap_outputs(
            runs,
            baseline_runs,
            out_dir,
            archive_runs=False,
        )
        print(f"Rebuilt cyclic best 3-swap outputs in {out_dir}", flush=True)
        return

    rows = []
    total_cases = (
        len(args.data_profiles)
        * len(args.generator_seeds)
        * int(args.samples)
        * len(args.off_cases)
    )
    progress = (
        tqdm(total=total_cases, unit="case", dynamic_ncols=True)
        if tqdm is not None
        else None
    )
    completed_cases = 0
    try:
        for profile in args.data_profiles:
            for generator_seed in args.generator_seeds:
                rng = np.random.RandomState(generator_seed)
                for sample in range(int(args.samples)):
                    V = generate_v_profile_from_rng(rng, args.N, args.L, profile=profile)
                    for off_case in args.off_cases:
                        K = int(off_case["K"])
                        off_pct = float(off_case["off_pct"])
                        rows.extend(
                            _evaluate_cyclic_best_3swap_case(
                                V,
                                K,
                                off_pct,
                                profile,
                                int(generator_seed),
                                int(sample),
                                args,
                            )
                        )
                        if progress is not None:
                            progress.set_postfix_str(
                                f"profile={profile}, seed={generator_seed}, "
                                f"sample={sample}, K={K}"
                            )
                            progress.update(1)
                        else:
                            completed_cases += 1
                            print(
                                f"[{completed_cases}/{total_cases}] profile={profile}, "
                                f"seed={generator_seed}, sample={sample}, K={K}",
                                flush=True,
                            )
    finally:
        if progress is not None:
            progress.close()

    runs = pd.DataFrame(rows)
    _write_cyclic_best_3swap_outputs(
        runs,
        baseline_runs,
        out_dir,
        archive_runs=True,
    )
    print(f"Wrote cyclic best 3-swap analysis to {out_dir}", flush=True)


def _evaluate_cyclic_best_3swap_case(V, K, off_pct, profile, generator_seed, sample, args):
    cyclic = best_cyclic_threshold_window(V, K, sigma=args.sigma, P=args.P)
    steps = refine_selection_by_ug_swaps_steps(
        V,
        cyclic["x"],
        max_swaps_values=(0, 1, 2, 3),
        sigma=args.sigma,
        P=args.P,
        K=K,
    )
    rows = []
    for max_swaps in (0, 1, 2, 3):
        result = steps[max_swaps]
        rows.append(
            {
                "profile": profile,
                "generator_seed": int(generator_seed),
                "sample": int(sample),
                "N": int(args.N),
                "L": int(args.L),
                "K": int(K),
                "K_off": int(args.N - K),
                "off_pct": float(off_pct),
                "active_pct": float(100.0 - off_pct),
                "sigma": float(args.sigma),
                "P": float(args.P),
                "seed_family": "cyclic_best",
                "method": f"cyclic_best+{max_swaps}swap",
                "max_swaps": int(max_swaps),
                "seed_T": int(cyclic["T"]),
                "seed_candidate_count": int(cyclic["candidate_count"]),
                "active_count": int(result["active_count"]),
                "add_candidate_count": int(result["add_candidate_count"]),
                "evaluated_swap_count": int(result["evaluated_swap_count"]),
                "swaps_applied": int(result["swaps_applied"]),
                "swap_history": result["swap_history"],
                "initial_u_bf": float(result["initial_u_bf"]),
                "initial_u_i": float(result["initial_u_i"]),
                "initial_u_g": float(result["initial_u_g"]),
                "u_bf": float(result["u_bf"]),
                "u_i": float(result["u_i"]),
                "u_g": float(result["u_g"]),
                "u_g_db": float(result["u_g_db"]),
                "seed_elapsed_seconds": float(cyclic["elapsed_seconds"]),
                "local_elapsed_seconds": float(result["elapsed_seconds"]),
                "elapsed_seconds": float(cyclic["elapsed_seconds"] + result["elapsed_seconds"]),
            }
        )
    return rows


def _load_cyclic_best_3swap_runs(out_dir):
    return _load_named_runs_from_dir(out_dir, "cyclic_best_3swap_runs.csv")


def _load_named_runs_from_dir(out_dir, csv_name):
    csv_path = out_dir / csv_name
    if csv_path.exists():
        return pd.read_csv(csv_path)

    archive_path = out_dir / "csv_data.tar.gz"
    if not archive_path.exists():
        raise FileNotFoundError(f"No {csv_name} or csv_data.tar.gz found in {out_dir}.")
    with tarfile.open(archive_path, "r:gz") as archive:
        member = archive.extractfile(csv_name)
        if member is None:
            raise FileNotFoundError(f"csv_data.tar.gz does not contain {csv_name}.")
        return pd.read_csv(member)


def _write_cyclic_best_3swap_outputs(runs, baseline_runs, out_dir, archive_runs):
    if runs.empty:
        raise RuntimeError("No cyclic best 3-swap rows were generated.")

    combined = _combined_previous_and_cyclic3_runs(baseline_runs, runs)
    case_cols = ["profile", "generator_seed", "sample", "K"]
    combined_best = combined.groupby(case_cols)["u_g"].transform("max")
    combined["fraction_best_combined"] = combined["u_g"] / combined_best
    cyclic_keys = set(
        zip(
            runs["profile"].astype(str),
            runs["generator_seed"].astype(int),
            runs["sample"].astype(int),
            runs["K"].astype(int),
            runs["method"].astype(str),
        )
    )
    cyclic_mask = [
        key in cyclic_keys
        for key in zip(
            combined["profile"].astype(str),
            combined["generator_seed"].astype(int),
            combined["sample"].astype(int),
            combined["K"].astype(int),
            combined["method"].astype(str),
        )
    ]
    runs = combined.loc[cyclic_mask].copy()

    runs_path = out_dir / "cyclic_best_3swap_runs.csv"
    if archive_runs or not (out_dir / "csv_data.tar.gz").exists():
        atomic_write_csv(runs, runs_path)

    summary = _build_cyclic_best_3swap_summary(runs)
    improvement = _build_cyclic_best_3swap_improvement(runs)
    combined_summary = _build_cyclic_best_3swap_combined_summary(combined)
    atomic_write_csv(summary, out_dir / "cyclic_best_3swap_summary.csv")
    atomic_write_csv(improvement, out_dir / "cyclic_best_3swap_improvement.csv")
    atomic_write_csv(combined_summary, out_dir / "combined_previous_vs_cyclic3_summary.csv")

    _write_cyclic_best_3swap_plots(
        runs,
        combined,
        summary,
        improvement,
        combined_summary,
        out_dir,
    )
    _write_cyclic_best_3swap_report(
        runs,
        summary,
        improvement,
        combined_summary,
        out_dir / "cyclic_best_3swap_report.md",
    )
    if archive_runs:
        _archive_csv_files(out_dir, [runs_path.name])


def _combined_previous_and_cyclic3_runs(baseline_runs, cyclic_runs):
    baseline = baseline_runs.copy()
    baseline = baseline[baseline["seed_family"] != "cyclic_best"].copy()
    common_cols = sorted(set(baseline.columns) | set(cyclic_runs.columns))
    baseline = baseline.reindex(columns=common_cols)
    cyclic = cyclic_runs.reindex(columns=common_cols)
    combined = pd.concat([baseline, cyclic], ignore_index=True)
    methods = set(CYCLIC_3SWAP_COMPARISON_METHODS)
    return combined[combined["method"].isin(methods)].copy()


def _build_cyclic_best_3swap_summary(runs):
    def q(value):
        return lambda data: data.quantile(value)

    return (
        runs.groupby(["profile", "off_pct", "active_pct", "K", "method", "max_swaps"], as_index=False)
        .agg(
            samples=("u_g", "count"),
            u_g_mean=("u_g", "mean"),
            u_g_p05=("u_g", q(0.05)),
            u_g_p50=("u_g", q(0.50)),
            u_g_p95=("u_g", q(0.95)),
            fraction_best_combined_mean=("fraction_best_combined", "mean"),
            fraction_best_combined_p05=("fraction_best_combined", q(0.05)),
            fraction_best_combined_p50=("fraction_best_combined", q(0.50)),
            fraction_best_combined_p95=("fraction_best_combined", q(0.95)),
            elapsed_mean=("elapsed_seconds", "mean"),
            elapsed_p50=("elapsed_seconds", "median"),
            swaps_applied_mean=("swaps_applied", "mean"),
            evaluated_swaps_mean=("evaluated_swap_count", "mean"),
            seed_T_mean=("seed_T", "mean"),
            seed_T_p50=("seed_T", "median"),
        )
        .sort_values(["profile", "K", "max_swaps"])
    )


def _build_cyclic_best_3swap_improvement(runs):
    pivot = runs.pivot_table(
        index=["profile", "generator_seed", "sample", "K"],
        columns="max_swaps",
        values=["u_g", "elapsed_seconds", "swaps_applied"],
        aggfunc="first",
    )
    rows = []
    for key, row in pivot.iterrows():
        profile, generator_seed, sample, K = key
        required = [("u_g", value) for value in (0, 1, 2, 3)]
        if not all(col in row.index for col in required):
            continue
        u0, u1, u2, u3 = (float(row[("u_g", value)]) for value in (0, 1, 2, 3))
        e0, e1, e2, e3 = (
            float(row[("elapsed_seconds", value)]) for value in (0, 1, 2, 3)
        )
        if u0 <= 0.0:
            continue
        rows.append(
            {
                "profile": profile,
                "generator_seed": int(generator_seed),
                "sample": int(sample),
                "K": int(K),
                "gain_1_vs_0_pct": 100.0 * (u1 / u0 - 1.0),
                "gain_2_vs_0_pct": 100.0 * (u2 / u0 - 1.0),
                "gain_3_vs_0_pct": 100.0 * (u3 / u0 - 1.0),
                "gain_2_vs_1_pct": 100.0 * (u2 / u1 - 1.0) if u1 > 0.0 else np.nan,
                "gain_3_vs_2_pct": 100.0 * (u3 / u2 - 1.0) if u2 > 0.0 else np.nan,
                "gain_3_vs_2_abs": u3 - u2,
                "runtime_ratio_3_vs_2": e3 / e2 if e2 > 0.0 else np.nan,
                "runtime_ratio_3_vs_0": e3 / e0 if e0 > 0.0 else np.nan,
                "swaps_applied_3": float(row[("swaps_applied", 3)]),
            }
        )

    gains = pd.DataFrame(rows)
    if gains.empty:
        return gains

    def q(value):
        return lambda data: data.quantile(value)

    return (
        gains.groupby(["profile", "K"], as_index=False)
        .agg(
            samples=("gain_3_vs_0_pct", "count"),
            gain_1_vs_0_pct_mean=("gain_1_vs_0_pct", "mean"),
            gain_2_vs_0_pct_mean=("gain_2_vs_0_pct", "mean"),
            gain_3_vs_0_pct_mean=("gain_3_vs_0_pct", "mean"),
            gain_3_vs_0_pct_p95=("gain_3_vs_0_pct", q(0.95)),
            gain_3_vs_2_pct_mean=("gain_3_vs_2_pct", "mean"),
            gain_3_vs_2_pct_p50=("gain_3_vs_2_pct", q(0.50)),
            gain_3_vs_2_pct_p95=("gain_3_vs_2_pct", q(0.95)),
            gain_3_vs_2_abs_mean=("gain_3_vs_2_abs", "mean"),
            runtime_ratio_3_vs_2_p50=("runtime_ratio_3_vs_2", q(0.50)),
            runtime_ratio_3_vs_0_p50=("runtime_ratio_3_vs_0", q(0.50)),
            swaps_applied_3_mean=("swaps_applied_3", "mean"),
        )
        .sort_values(["profile", "K"])
    )


def _build_cyclic_best_3swap_combined_summary(combined):
    case_cols = ["profile", "generator_seed", "sample", "K"]
    case_best = combined.groupby(case_cols)["u_g"].transform("max")
    tolerance = np.maximum(1e-9, 1e-12 * np.abs(case_best))
    combined = combined.copy()
    combined["best_observed_tie"] = np.abs(combined["u_g"] - case_best) <= tolerance

    def q(value):
        return lambda data: data.quantile(value)

    return (
        combined.groupby(["profile", "off_pct", "active_pct", "K", "method"], as_index=False)
        .agg(
            samples=("u_g", "count"),
            u_g_mean=("u_g", "mean"),
            fraction_best_combined_mean=("fraction_best_combined", "mean"),
            fraction_best_combined_p05=("fraction_best_combined", q(0.05)),
            best_observed_tie_rate=("best_observed_tie", "mean"),
            elapsed_p50=("elapsed_seconds", "median"),
        )
        .sort_values(
            ["profile", "K", "fraction_best_combined_mean"],
            ascending=[True, True, False],
        )
    )


def _write_cyclic_best_3swap_plots(runs, combined, summary, improvement, combined_summary, out_dir):
    use_agg_backend()
    import matplotlib.pyplot as plt

    _plot_seed_cdf(
        runs,
        "u_g",
        "U_G",
        "Cyclic best raw U_G CDF by swap depth",
        out_dir / "cyclic_best_3swap_raw_u_g_cdf.png",
    )
    _plot_seed_cdf(
        combined,
        "fraction_best_combined",
        "U_G / best observed combined",
        "Combined fraction-best CDF with cyclic 3-swap",
        out_dir / "cyclic_best_3swap_fraction_best_combined_cdf.png",
    )
    _plot_cyclic_best_3swap_gain_cdf(
        runs,
        out_dir / "cyclic_best_3swap_gain_cdf.png",
    )
    _plot_cyclic_best_3swap_mean_fraction(
        summary,
        out_dir / "cyclic_best_3swap_mean_fraction_vs_depth.png",
    )
    _plot_cyclic_best_3swap_runtime(
        summary,
        out_dir / "cyclic_best_3swap_runtime_vs_depth.png",
    )
    _plot_cyclic_best_3swap_tie_rate(
        combined_summary,
        out_dir / "cyclic_best_3swap_best_observed_tie_rate.png",
    )
    plt.close("all")


def _plot_cyclic_best_3swap_gain_cdf(runs, out_path):
    use_agg_backend()
    import matplotlib.pyplot as plt

    pivot = runs.pivot_table(
        index=["profile", "generator_seed", "sample", "K"],
        columns="max_swaps",
        values="u_g",
        aggfunc="first",
    )
    rows = []
    for key, row in pivot.iterrows():
        profile, generator_seed, sample, K = key
        del generator_seed, sample
        if not all(value in row.index for value in (0, 2, 3)):
            continue
        u0 = float(row[0])
        u2 = float(row[2])
        u3 = float(row[3])
        if u0 <= 0.0 or u2 <= 0.0:
            continue
        rows.append(
            {
                "profile": profile,
                "K": int(K),
                "gain_kind": "3swap vs 0swap",
                "gain_pct": 100.0 * (u3 / u0 - 1.0),
            }
        )
        rows.append(
            {
                "profile": profile,
                "K": int(K),
                "gain_kind": "3swap vs 2swap",
                "gain_pct": 100.0 * (u3 / u2 - 1.0),
            }
        )
    gains = pd.DataFrame(rows)
    if gains.empty:
        return

    K_values = sorted(gains["K"].unique())
    fig, axes = plt.subplots(1, len(K_values), figsize=(7.0 * len(K_values), 4.4), squeeze=False)
    colors = {"3swap vs 0swap": "#4C78A8", "3swap vs 2swap": "#F58518"}
    for col, K in enumerate(K_values):
        ax = axes[0, col]
        data = gains[gains["K"] == K]
        for gain_kind in ("3swap vs 0swap", "3swap vs 2swap"):
            values = data[data["gain_kind"] == gain_kind]["gain_pct"]
            x_values, y_values = _empirical_cdf(values)
            if len(x_values) == 0:
                continue
            ax.step(
                y_values,
                x_values,
                where="post",
                linewidth=1.4,
                color=colors[gain_kind],
                label=gain_kind,
            )
        ax.set_title(f"Cyclic 3-swap gain CDF, K={int(K)}")
        ax.set_xlabel("Cumulative fraction of examples")
        if col == 0:
            ax.set_ylabel("U_G gain, %")
        ax.set_xlim(0.0, 1.02)
        ax.grid(True, alpha=0.25)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.995, 0.5), fontsize=8)
    fig.tight_layout(rect=(0.0, 0.0, 0.84, 0.94))
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_cyclic_best_3swap_mean_fraction(summary, out_path):
    if summary.empty:
        return
    use_agg_backend()
    import matplotlib.pyplot as plt

    K_values = sorted(summary["K"].unique())
    fig, axes = plt.subplots(1, len(K_values), figsize=(6.5 * len(K_values), 4.4), squeeze=False)
    for col, K in enumerate(K_values):
        ax = axes[0, col]
        data = summary[summary["K"] == K].sort_values("max_swaps")
        ax.plot(
            data["max_swaps"],
            data["fraction_best_combined_mean"],
            marker="o",
            linewidth=1.5,
            color="#4C78A8",
        )
        off_pct = float(data["off_pct"].iloc[0]) if not data.empty else np.nan
        ax.set_title(f"Cyclic quality vs swaps, K={int(K)} ({off_pct:g}% off)")
        ax.set_xlabel("Max greedy swaps")
        if col == 0:
            ax.set_ylabel("Mean U_G / best observed combined")
        ax.set_xticks([0, 1, 2, 3])
        ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_cyclic_best_3swap_runtime(summary, out_path):
    if summary.empty:
        return
    use_agg_backend()
    import matplotlib.pyplot as plt

    K_values = sorted(summary["K"].unique())
    fig, axes = plt.subplots(1, len(K_values), figsize=(6.5 * len(K_values), 4.4), squeeze=False)
    for col, K in enumerate(K_values):
        ax = axes[0, col]
        data = summary[summary["K"] == K].sort_values("max_swaps")
        ax.plot(
            data["max_swaps"],
            data["elapsed_p50"],
            marker="o",
            linewidth=1.5,
            color="#E45756",
        )
        off_pct = float(data["off_pct"].iloc[0]) if not data.empty else np.nan
        ax.set_title(f"Cyclic runtime vs swaps, K={int(K)} ({off_pct:g}% off)")
        ax.set_xlabel("Max greedy swaps")
        if col == 0:
            ax.set_ylabel("Median runtime, seconds")
        ax.set_xticks([0, 1, 2, 3])
        ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_cyclic_best_3swap_tie_rate(combined_summary, out_path):
    if combined_summary.empty:
        return
    use_agg_backend()
    import matplotlib.pyplot as plt

    K_values = sorted(combined_summary["K"].unique())
    fig, axes = plt.subplots(1, len(K_values), figsize=(7.5 * len(K_values), 4.8), squeeze=False)
    colors = _method_colors(list(CYCLIC_3SWAP_COMPARISON_METHODS))
    for col, K in enumerate(K_values):
        ax = axes[0, col]
        data = combined_summary[combined_summary["K"] == K].copy()
        data = data.sort_values(
            ["best_observed_tie_rate", "fraction_best_combined_mean"],
            ascending=False,
        )
        positions = np.arange(len(data))
        ax.barh(
            positions,
            data["best_observed_tie_rate"],
            color=[colors.get(method, "#777777") for method in data["method"]],
        )
        ax.set_yticks(positions)
        ax.set_yticklabels(data["method"], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlim(0.0, 1.0)
        ax.set_xlabel("Best-observed tie rate")
        ax.set_title(f"Best observed tie rate, K={int(K)}")
        ax.grid(True, axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _write_cyclic_best_3swap_report(runs, summary, improvement, combined_summary, out_path):
    lines = [
        "# Cyclic Best-T 3-Swap Analysis",
        "",
        "This focused follow-up reruns only `cyclic_best` with `0/1/2/3` greedy raw-`U_G` swaps.",
        "Non-cyclic comparison methods are loaded from the configured baseline U_G swap result directory.",
        "",
        "## K Semantics",
        "",
        "`K` is the number of active/kept antennas selected by the solver.",
        "`K_off = N - K` is the number of disabled antennas.",
        "For `N=1000`, `K=750` means `25% off`; `K=500` means `50% off`.",
        "",
        "## Clear Answer",
        "",
        *_cyclic_best_3swap_answer_lines(summary, improvement, combined_summary),
        "",
        "## Cyclic Swap Depth Summary",
        "",
        _markdown_table(
            summary,
            [
                ("profile", "profile"),
                ("K", "K"),
                ("method", "method"),
                ("u_g_mean", "mean U_G"),
                ("fraction_best_combined_mean", "mean fraction combined best"),
                ("fraction_best_combined_p05", "p05 fraction"),
                ("elapsed_p50", "median runtime"),
                ("swaps_applied_mean", "mean swaps applied"),
            ],
        ),
        "",
        "## Third-Swap Marginal Gain",
        "",
        _markdown_table(
            improvement,
            [
                ("profile", "profile"),
                ("K", "K"),
                ("gain_1_vs_0_pct_mean", "mean +1 vs 0 %"),
                ("gain_2_vs_0_pct_mean", "mean +2 vs 0 %"),
                ("gain_3_vs_0_pct_mean", "mean +3 vs 0 %"),
                ("gain_3_vs_2_pct_mean", "mean +3 vs +2 %"),
                ("gain_3_vs_2_pct_p95", "p95 +3 vs +2 %"),
                ("runtime_ratio_3_vs_2_p50", "median runtime x, 3/2"),
                ("runtime_ratio_3_vs_0_p50", "median runtime x, 3/0"),
            ],
        ),
        "",
        "## Combined Comparison",
        "",
        "The table compares the previous major heuristics with the new cyclic 3-swap result.",
        "",
    ]

    for (profile, K), chunk in combined_summary.groupby(["profile", "K"], sort=True):
        ranked = chunk.sort_values(
            ["fraction_best_combined_mean", "best_observed_tie_rate"],
            ascending=False,
        )
        off_pct = float(ranked["off_pct"].iloc[0])
        lines.extend(
            [
                f"### {profile}, K={int(K)} ({off_pct:g}% off)",
                "",
                _markdown_table(
                    ranked,
                    [
                        ("method", "method"),
                        ("u_g_mean", "mean U_G"),
                        ("fraction_best_combined_mean", "mean fraction best"),
                        ("fraction_best_combined_p05", "p05 fraction"),
                        ("best_observed_tie_rate", "best tie rate"),
                        ("elapsed_p50", "median runtime"),
                    ],
                ),
                "",
            ]
        )

    lines.extend(
        [
            "## Output Figures",
            "",
            "- `cyclic_best_3swap_raw_u_g_cdf.png`",
            "- `cyclic_best_3swap_fraction_best_combined_cdf.png`",
            "- `cyclic_best_3swap_gain_cdf.png`",
            "- `cyclic_best_3swap_mean_fraction_vs_depth.png`",
            "- `cyclic_best_3swap_runtime_vs_depth.png`",
            "- `cyclic_best_3swap_best_observed_tie_rate.png`",
        ]
    )
    out_path.write_text("\n".join(lines), encoding="utf-8")


def _cyclic_best_3swap_answer_lines(summary, improvement, combined_summary):
    lines = []
    for (profile, K), chunk in summary.groupby(["profile", "K"], sort=True):
        off_pct = float(chunk["off_pct"].iloc[0])
        row2 = chunk[chunk["max_swaps"] == 2].iloc[0]
        row3 = chunk[chunk["max_swaps"] == 3].iloc[0]
        gain = improvement[(improvement["profile"] == profile) & (improvement["K"] == K)]
        gain_row = gain.iloc[0] if not gain.empty else None
        combined = combined_summary[
            (combined_summary["profile"] == profile) & (combined_summary["K"] == K)
        ].sort_values("fraction_best_combined_mean", ascending=False)
        best = combined.iloc[0] if not combined.empty else None
        lines.extend([f"### {profile}, K={int(K)} ({off_pct:g}% off)", ""])
        lines.append(
            f"- `cyclic_best+3swap` mean fraction of combined best: "
            f"`{_format_report_float(row3['fraction_best_combined_mean'])}`."
        )
        lines.append(
            f"- `cyclic_best+2swap` mean fraction of combined best: "
            f"`{_format_report_float(row2['fraction_best_combined_mean'])}`."
        )
        if gain_row is not None:
            lines.append(
                f"- Third-swap marginal gain over 2 swaps: mean "
                f"`{_format_report_float(gain_row['gain_3_vs_2_pct_mean'])}%`, "
                f"p95 `{_format_report_float(gain_row['gain_3_vs_2_pct_p95'])}%`."
            )
            lines.append(
                f"- Runtime cost: median 3-swap runtime is "
                f"`{_format_report_float(gain_row['runtime_ratio_3_vs_2_p50'])}x` "
                f"the 2-swap runtime."
            )
        if best is not None:
            lines.append(
                f"- Best combined method by mean fraction is `{best['method']}` "
                f"at `{_format_report_float(best['fraction_best_combined_mean'])}`."
            )
        lines.append("")
    lines.append(
        "Conclusion: the 3rd swap is a measurable but very small quality polish. It should be used only when best quality matters more than runtime. The previous conclusion remains: cyclic best-T is the primary heuristic, and local swaps are an incremental refinement."
    )
    return lines


def _archive_csv_files(out_dir, csv_names):
    archive_path = out_dir / "csv_data.tar.gz"
    with tarfile.open(archive_path, "w:gz") as archive:
        for name in csv_names:
            path = out_dir / name
            if not path.exists():
                continue
            archive.add(path, arcname=name)
            path.unlink()


UNIFIED_LOCAL_SWAP_FAMILIES = (
    "StrongWeak",
    "H3Threshold-T0.05N-Gen",
    "H3Threshold-CyclicBestT-Gen",
    "CapWindowFull-Gen",
    "CapSubmod-Gen",
)


def _unified_l_values(args):
    values = args.L_values if args.L_values is not None else [args.L]
    values = [int(value) for value in values]
    if any(value <= 0 for value in values):
        raise ValueError("All L values must be positive.")
    return values


def _unified_sigma_values(args):
    values = args.sigmas if args.sigmas is not None else [args.sigma]
    values = [float(value) for value in values]
    if any(value <= 0.0 for value in values):
        raise ValueError("All sigma values must be positive.")
    return values


def run_unified_local_swap_comparison(args):
    if int(args.workers) != 1:
        raise ValueError(
            "--unified-local-swap-comparison is intentionally sequential; use --workers 1."
        )

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    runs_path = out_dir / "unified_local_swap_runs.csv"
    existing, completed = (
        _load_existing_unified_runs(runs_path) if args.resume else (pd.DataFrame(), set())
    )
    rows = existing.to_dict("records") if not existing.empty else []
    l_values = _unified_l_values(args)
    sigma_values = _unified_sigma_values(args)
    total_runs = (
        len(l_values)
        * len(sigma_values)
        * len(args.data_profiles)
        * len(args.generator_seeds)
        * int(args.samples)
        * len(args.off_cases)
        * len(UNIFIED_LOCAL_SWAP_FAMILIES)
        * 2
    )
    progress = (
        tqdm(
            total=total_runs,
            initial=min(len(completed), total_runs),
            unit="run",
            dynamic_ncols=True,
        )
        if tqdm is not None
        else None
    )
    new_since_checkpoint = 0
    total_new = 0

    try:
        for L in l_values:
            for profile in args.data_profiles:
                for generator_seed in args.generator_seeds:
                    rng = np.random.RandomState(int(generator_seed))
                    for sample in range(int(args.samples)):
                        V = generate_v_profile_from_rng(rng, args.N, L, profile=profile)
                        for off_case in args.off_cases:
                            K = int(off_case["K"])
                            off_pct = float(off_case["off_pct"])
                            for sigma in sigma_values:
                                case_base = {
                                    "data_profile": str(profile),
                                    "generator_seed": int(generator_seed),
                                    "sample": int(sample),
                                    "N": int(args.N),
                                    "L": int(L),
                                    "K": K,
                                    "K_active": K,
                                    "K_off": int(args.N - K),
                                    "off_pct": off_pct,
                                    "active_pct": float(100.0 - off_pct),
                                    "sigma": float(sigma),
                                    "P": float(args.P),
                                }
                                seed_specs = _unified_seed_specs_for_case(V, K, args, case_base)
                                for seed in seed_specs:
                                    for max_swaps in (0, 1):
                                        method = f"{seed['seed_family']}+{int(max_swaps)}swap"
                                        key = _unified_run_key(case_base, method)
                                        if key in completed:
                                            continue
                                        row = _unified_local_swap_row(
                                            V,
                                            case_base,
                                            seed,
                                            max_swaps,
                                            include_selection_details=not args.compact_runs,
                                        )
                                        rows.append(row)
                                        completed.add(key)
                                        total_new += 1
                                        new_since_checkpoint += 1
                                        message = (
                                            f"L={L}, sigma={sigma:g}, profile={profile}, "
                                            f"seed={generator_seed}, sample={sample}, "
                                            f"off={off_pct:g}%, method={row['method']}"
                                        )
                                        if progress is not None:
                                            progress.set_postfix_str(message)
                                            progress.update(1)
                                        else:
                                            done = min(len(completed), total_runs)
                                            print(
                                                f"[{done}/{total_runs}] {message}",
                                                flush=True,
                                            )
                                        if (
                                            args.checkpoint_every
                                            and new_since_checkpoint >= args.checkpoint_every
                                        ):
                                            atomic_write_csv(pd.DataFrame(rows), runs_path)
                                            new_since_checkpoint = 0
    finally:
        if progress is not None:
            progress.close()

    runs = _attach_unified_best_observed(pd.DataFrame(rows))
    summary = _build_unified_local_swap_summary(runs)
    overall = _build_unified_local_swap_overall_summary(summary)

    atomic_write_csv(runs, runs_path)
    atomic_write_csv(summary, out_dir / "unified_local_swap_summary.csv")
    atomic_write_csv(overall, out_dir / "unified_local_swap_overall_summary.csv")
    _write_unified_local_swap_report(runs, summary, overall, out_dir, args)
    print(
        f"Completed {total_new} new unified local-swap runs; "
        f"{len(rows)} total rows in {runs_path}.",
        flush=True,
    )


def _load_existing_unified_runs(path):
    if not path.exists():
        return pd.DataFrame(), set()
    runs = pd.read_csv(path)
    runs = runs.drop(
        columns=[
            "best_observed_u_g",
            "fraction_best_observed_u_g",
            "is_best_observed",
        ],
        errors="ignore",
    )
    completed = {_unified_run_key(row, row["method"]) for row in runs.to_dict("records")}
    return runs, completed


def _unified_run_key(case_or_row, method):
    return (
        str(case_or_row["data_profile"]),
        int(case_or_row["generator_seed"]),
        int(case_or_row["sample"]),
        int(case_or_row["N"]),
        int(case_or_row["L"]),
        int(case_or_row["K"]),
        _canonical_float(case_or_row["off_pct"]),
        _canonical_float(case_or_row["sigma"]),
        _canonical_float(case_or_row["P"]),
        str(method),
    )


def _canonical_float(value):
    return f"{float(value):.12g}"


def _unified_seed_specs_for_case(V, K, args, case_base):
    N = int(args.N)
    sigma = float(case_base["sigma"])
    P = float(case_base["P"])
    row_power = np.sum(np.abs(V) ** 2, axis=1).real
    order = np.argsort(row_power)[::-1]
    radius = _candidate_radius(K, None)
    specs = []

    def add_spec(
        seed_family,
        x,
        seed_position,
        candidate_kind,
        seed_elapsed_seconds,
        seed_candidate_count=np.nan,
        candidate_pool=None,
        candidate_pool_kind="rank_neighborhood",
    ):
        specs.append(
            {
                "seed_family": seed_family,
                "seed_position": seed_position,
                "seed_candidate_count": seed_candidate_count,
                "candidate_kind": candidate_kind,
                "candidate_radius": radius,
                "candidate_pool": candidate_pool,
                "candidate_pool_kind": candidate_pool_kind,
                "seed_elapsed_seconds": float(seed_elapsed_seconds),
                "x": np.asarray(x, dtype=int),
            }
        )

    T_005 = int(np.clip(round(0.05 * N), 0, max(0, N - K)))
    started_at = time.perf_counter()
    h3_t005_x = solve_h3(
        V,
        K,
        target_obj="gen",
        sigma=sigma,
        P=P,
        t_tests=(T_005,),
        include_phase_nulling=False,
    )
    h3_t005_elapsed = time.perf_counter() - started_at
    h3_t005_position = contiguous_threshold_window_T(V, np.flatnonzero(h3_t005_x))
    h3_t005_pool = _linear_pool_or_none(order, h3_t005_x, h3_t005_position, K, radius)
    add_spec(
        "H3Threshold-T0.05N-Gen",
        h3_t005_x,
        T_005 if h3_t005_position is None else int(h3_t005_position),
        "h3_threshold_T0p05N",
        h3_t005_elapsed,
        seed_candidate_count=1 + int(0 < T_005 <= N - K),
        candidate_pool=h3_t005_pool,
        candidate_pool_kind="linear_boundary" if h3_t005_pool is not None else "rank_neighborhood",
    )

    started_at = time.perf_counter()
    cyclic = best_cyclic_threshold_window(V, K, sigma=sigma, P=P)
    cyclic_elapsed = time.perf_counter() - started_at
    cyclic_pool = _cyclic_boundary_add_pool(order, cyclic["x"], int(cyclic["T"]), K, radius)
    add_spec(
        "H3Threshold-CyclicBestT-Gen",
        cyclic["x"],
        int(cyclic["T"]),
        "cyclic_threshold_window",
        cyclic_elapsed,
        seed_candidate_count=int(cyclic["candidate_count"]),
        candidate_pool=cyclic_pool,
        candidate_pool_kind="cyclic_boundary",
    )

    random_state = solver_random_state(
        case_base["generator_seed"],
        case_base["sample"],
        K,
        UNIFIED_LOCAL_SWAP_FAMILIES.index("CapWindowFull-Gen"),
    )
    started_at = time.perf_counter()
    cap_window_full_x = solve_cap_window_full_gen(
        V,
        K,
        sigma=sigma,
        P=P,
        random_state=random_state,
    )
    cap_window_full_elapsed = time.perf_counter() - started_at
    cap_window_full_position = contiguous_threshold_window_T(
        V,
        np.flatnonzero(cap_window_full_x),
    )
    cap_window_full_pool = _linear_pool_or_none(
        order,
        cap_window_full_x,
        cap_window_full_position,
        K,
        radius,
    )
    add_spec(
        "CapWindowFull-Gen",
        cap_window_full_x,
        cap_window_full_position,
        "cap_window_full",
        cap_window_full_elapsed,
        candidate_pool=cap_window_full_pool,
        candidate_pool_kind="linear_boundary"
        if cap_window_full_pool is not None
        else "rank_neighborhood",
    )

    random_state = solver_random_state(
        case_base["generator_seed"],
        case_base["sample"],
        K,
        UNIFIED_LOCAL_SWAP_FAMILIES.index("CapSubmod-Gen"),
    )
    started_at = time.perf_counter()
    cap_submod_x = cap_submodular_seed_gen(
        V,
        K,
        sigma=sigma,
        P=P,
        random_state=random_state,
    )
    add_spec(
        "CapSubmod-Gen",
        cap_submod_x,
        np.nan,
        "cap_submodular",
        time.perf_counter() - started_at,
    )

    started_at = time.perf_counter()
    strong_x = solve_h3_strong_weak(V, K, sigma=sigma, P=P)
    strong_elapsed = time.perf_counter() - started_at
    off_count = int(args.N - K)
    strong_position = off_count - off_count // 2
    strong_pool = _linear_boundary_add_pool(order, strong_x, strong_position, K, radius)
    add_spec(
        "StrongWeak",
        strong_x,
        strong_position,
        "strong_weak",
        strong_elapsed,
        seed_candidate_count=1,
        candidate_pool=strong_pool,
        candidate_pool_kind="linear_boundary",
    )
    return specs


def _linear_pool_or_none(order, x, seed_position, K, radius):
    if seed_position is None or not np.isfinite(seed_position):
        return None
    return _linear_boundary_add_pool(order, x, int(seed_position), K, radius)


def _unified_local_swap_row(
    V,
    case_base,
    seed,
    max_swaps,
    include_selection_details=True,
):
    result = refine_selection_by_swaps(
        V,
        seed["x"],
        max_swaps=int(max_swaps),
        sigma=float(case_base["sigma"]),
        P=float(case_base["P"]),
        candidate_radius=seed["candidate_radius"],
        candidate_pool=seed["candidate_pool"],
        seed_position=0 if pd.isna(seed["seed_position"]) else int(seed["seed_position"]),
    )
    local_elapsed = float(result["elapsed_seconds"])
    total_elapsed = float(seed["seed_elapsed_seconds"]) + local_elapsed
    row = {
        **case_base,
        "method": f"{seed['seed_family']}+{int(max_swaps)}swap",
        "seed_family": seed["seed_family"],
        "max_swaps": int(max_swaps),
        "seed_position": seed["seed_position"],
        "seed_candidate_count": seed["seed_candidate_count"],
        "seed_elapsed_seconds": float(seed["seed_elapsed_seconds"]),
        "candidate_kind": seed["candidate_kind"]
        if int(max_swaps) == 0
        else f"{seed['candidate_kind']}_local_swap",
        "candidate_pool_kind": seed["candidate_pool_kind"],
        "candidate_radius": int(result["candidate_radius"]),
        "add_candidate_count": int(result["add_candidate_count"]),
        "evaluated_swap_count": int(result["evaluated_swap_count"]),
        "candidate_count": int(result["candidate_count"]),
        "swaps_applied": int(result["swaps_applied"]),
        "active_count": int(np.sum(result["x"])),
        "initial_u_bf": float(result["initial_u_bf"]),
        "initial_u_i": float(result["initial_u_i"]),
        "initial_u_g": float(result["initial_u_g"]),
        "u_bf": float(result["u_bf"]),
        "u_i": float(result["u_i"]),
        "u_g": float(result["u_g"]),
        "u_g_db": 10.0 * np.log10(max(float(result["u_g"]), np.finfo(float).tiny)),
        "local_elapsed_seconds": local_elapsed,
        "total_elapsed_seconds": total_elapsed,
    }
    if include_selection_details:
        row.update(
            {
                "initial_subset": _subset_to_string(result["initial_subset"]),
                "subset": _subset_to_string(result["subset"]),
                "swap_history": result["swap_history"],
            }
        )
    return row


def _attach_unified_best_observed(runs):
    if runs.empty:
        return runs
    case_keys = [
        "data_profile",
        "generator_seed",
        "sample",
        "N",
        "L",
        "K",
        "off_pct",
        "sigma",
        "P",
    ]
    runs = runs.copy()
    runs["best_observed_u_g"] = runs.groupby(case_keys)["u_g"].transform("max")
    runs["fraction_best_observed_u_g"] = (
        runs["u_g"].astype(float) / runs["best_observed_u_g"].astype(float)
    )
    runs["is_best_observed"] = (
        runs["u_g"].astype(float) >= runs["best_observed_u_g"].astype(float) - 1e-9
    )
    return runs


def _build_unified_local_swap_summary(runs):
    if runs.empty:
        return pd.DataFrame()

    def q(value):
        return lambda data: data.quantile(value)

    return (
        runs.groupby(
            [
                "data_profile",
                "N",
                "L",
                "off_pct",
                "K",
                "sigma",
                "P",
                "seed_family",
                "max_swaps",
                "method",
            ],
            as_index=False,
        )
        .agg(
            cases=("u_g", "count"),
            u_g_mean=("u_g", "mean"),
            u_g_p05=("u_g", q(0.05)),
            u_g_p50=("u_g", q(0.50)),
            u_g_p95=("u_g", q(0.95)),
            u_g_db_mean=("u_g_db", "mean"),
            u_g_db_p05=("u_g_db", q(0.05)),
            u_g_db_p50=("u_g_db", q(0.50)),
            u_g_db_p95=("u_g_db", q(0.95)),
            fraction_best_observed_mean=("fraction_best_observed_u_g", "mean"),
            fraction_best_observed_p05=("fraction_best_observed_u_g", q(0.05)),
            fraction_best_observed_p50=("fraction_best_observed_u_g", q(0.50)),
            best_observed_rate=("is_best_observed", "mean"),
            seed_elapsed_seconds_mean=("seed_elapsed_seconds", "mean"),
            local_elapsed_seconds_mean=("local_elapsed_seconds", "mean"),
            total_elapsed_seconds_mean=("total_elapsed_seconds", "mean"),
            add_candidate_count_mean=("add_candidate_count", "mean"),
            evaluated_swap_count_mean=("evaluated_swap_count", "mean"),
            swaps_applied_mean=("swaps_applied", "mean"),
        )
        .sort_values(
            ["data_profile", "L", "sigma", "off_pct", "max_swaps", "seed_family"]
        )
    )


def _build_unified_local_swap_overall_summary(summary):
    if summary.empty:
        return pd.DataFrame()
    rows = []
    for (seed_family, max_swaps, method), group in summary.groupby(
        ["seed_family", "max_swaps", "method"],
        sort=True,
    ):
        weights = group["cases"].astype(float)
        total = float(weights.sum())
        row = {
            "seed_family": seed_family,
            "max_swaps": int(max_swaps),
            "method": method,
            "cases": int(total),
        }
        for col in (
            "u_g_mean",
            "u_g_db_mean",
            "fraction_best_observed_mean",
            "best_observed_rate",
            "seed_elapsed_seconds_mean",
            "local_elapsed_seconds_mean",
            "total_elapsed_seconds_mean",
            "add_candidate_count_mean",
            "evaluated_swap_count_mean",
            "swaps_applied_mean",
        ):
            row[col] = float(np.average(group[col].astype(float), weights=weights))
        rows.append(row)
    return pd.DataFrame(rows).sort_values(
        ["fraction_best_observed_mean", "best_observed_rate"],
        ascending=[False, False],
    )


def _write_unified_local_swap_report(runs, summary, overall, out_dir, args):
    l_values = _unified_l_values(args)
    sigma_values = _unified_sigma_values(args)
    lines = [
        "# Unified Local Swap Comparison",
        "",
        f"- N: {args.N}",
        f"- L values: {', '.join(str(value) for value in l_values)}",
        f"- sigma values: {', '.join(format_number_slug(value) for value in sigma_values)}",
        f"- off-pcts: {', '.join(format_number_slug(case['off_pct']) for case in args.off_cases)}",
        f"- samples: {args.samples}",
        f"- generator seeds: {', '.join(str(value) for value in args.generator_seeds)}",
        f"- data profiles: {', '.join(args.data_profiles)}",
        "- workers: 1 sequential run",
        f"- compact runs: {bool(args.compact_runs)}",
        "- local search: existing `algorithms.h3_threshold_local.refine_selection_by_swaps`.",
        "- local swap depths: `0` and `1`.",
        "- default candidate radius: `max(8, ceil(0.05K))`.",
        "",
        "## Seed Families",
        "",
        *[f"- `{name}`" for name in UNIFIED_LOCAL_SWAP_FAMILIES],
        "",
        "## Overall Summary",
        "",
        _markdown_table(
            overall,
            [
                ("method", "method"),
                ("cases", "cases"),
                ("fraction_best_observed_mean", "mean fraction best"),
                ("best_observed_rate", "best rate"),
                ("u_g_db_mean", "mean U_G dB"),
                ("total_elapsed_seconds_mean", "mean total seconds"),
                ("add_candidate_count_mean", "mean add candidates"),
                ("evaluated_swap_count_mean", "mean evaluated swaps"),
            ],
        ),
        "",
        "## Artifacts",
        "",
        "- `unified_local_swap_runs.csv`",
        "- `unified_local_swap_summary.csv`",
        "- `unified_local_swap_overall_summary.csv`",
    ]
    (Path(out_dir) / "unified_local_swap_report.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def _subset_to_string(subset):
    return " ".join(str(int(value)) for value in subset)


def _ensure_single_experiment_mode(args):
    modes = [
        args.unified_local_swap_comparison,
        args.ug_swap_seed_comparison,
        args.cyclic_best_3swap_analysis,
        args.threshold_local_exact_analysis,
        args.threshold_real_off_cyclic_local_exact_analysis,
        args.threshold_active_k_cyclic_local_exact_analysis,
        args.threshold_large_cyclic_local_analysis,
        args.threshold_large_cyclic_honest_local_analysis,
    ]
    if sum(bool(mode) for mode in modes) > 1:
        raise ValueError("Use only one experiment mode flag at a time.")


def _print_existing_outputs(out_dir, output_names, docs_path=None):
    print("Wrote:", flush=True)
    for output_name in output_names:
        output_path = out_dir / output_name
        if output_path.exists():
            print(f"  {output_path}", flush=True)
    if docs_path is not None and docs_path.exists():
        print(f"  {docs_path}", flush=True)


def _run_threshold_analysis_mode(args):
    if args.threshold_large_cyclic_honest_local_analysis:
        if args.K_values is None:
            raise ValueError(
                "--threshold-large-cyclic-honest-local-analysis requires --K-values."
            )
        docs_path = Path("docs/local_threshold_large_honest_cyclic_report.md")
        run_large_cyclic_honest_local_analysis(
            out_dir=args.out_dir,
            N=args.N,
            K_values=args.K_values,
            profiles=args.data_profiles,
            generator_seeds=args.generator_seeds,
            samples=args.samples,
            L=args.L,
            sigma=args.sigma,
            P=args.P,
            docs_path=docs_path,
        )
        _print_existing_outputs(
            args.out_dir,
            (
                "csv_data.tar.gz",
                "local_threshold_large_honest_report.md",
                "large_raw_u_g_cdf_by_K.png",
                "large_fraction_cyclic_seed_cdf_by_K.png",
                "large_fraction_best_observed_cdf_by_K.png",
                "large_mean_fraction_cyclic_seed_by_K.png",
                "large_best_cyclic_T_boxplot.png",
                "large_best_cyclic_T_over_N.png",
                "large_runtime_by_method.png",
            ),
            docs_path=docs_path,
        )
        return True

    if args.threshold_large_cyclic_local_analysis:
        if args.K_values is None:
            raise ValueError("--threshold-large-cyclic-local-analysis requires --K-values.")
        docs_path = Path("docs/local_threshold_large_cyclic_report.md")
        run_large_cyclic_local_analysis(
            out_dir=args.out_dir,
            N=args.N,
            K_values=args.K_values,
            profiles=args.data_profiles,
            generator_seeds=args.generator_seeds,
            samples=args.samples,
            L=args.L,
            sigma=args.sigma,
            P=args.P,
            docs_path=docs_path,
        )
        _print_existing_outputs(
            args.out_dir,
            (
                "csv_data.tar.gz",
                "local_threshold_large_report.md",
                "large_raw_u_g_cdf_by_K.png",
                "large_fraction_cyclic_seed_cdf_by_K.png",
                "large_fraction_best_observed_cdf_by_K.png",
                "large_mean_fraction_cyclic_seed_by_K.png",
                "large_best_cyclic_T_boxplot.png",
                "large_best_cyclic_T_over_N.png",
                "large_runtime_by_method.png",
            ),
            docs_path=docs_path,
        )
        return True

    if args.threshold_active_k_cyclic_local_exact_analysis:
        docs_path = Path("docs/local_threshold_active_k_cyclic_report.md")
        run_active_k_cyclic_local_exact_analysis(
            out_dir=args.out_dir,
            n_values=args.N_values,
            active_pcts=args.K_pcts,
            profiles=args.data_profiles,
            generator_seeds=args.generator_seeds,
            samples=args.samples,
            L=args.L,
            sigma=args.sigma,
            P=args.P,
            exact_source_dir=args.exact_source_dir,
            exact_time_limit=args.exact_time_limit,
            docs_path=docs_path,
        )
        _print_existing_outputs(
            args.out_dir,
            (
                "csv_data.tar.gz",
                "local_threshold_active_k_report.md",
                "active_k_raw_u_g_cdf_by_requested_active_pct.png",
                "active_k_fraction_exact_cdf.png",
                "active_k_mean_fraction_by_requested_active_pct.png",
                "active_k_exact_recovery_by_requested_active_pct.png",
                "active_k_best_cyclic_start_hist.png",
                "active_k_failure_diagnostics.png",
                "active_k_runtime_by_method.png",
            ),
            docs_path=docs_path,
        )
        return True

    if args.threshold_real_off_cyclic_local_exact_analysis:
        docs_path = Path("docs/local_threshold_real_off_cyclic_report.md")
        run_real_off_cyclic_local_exact_analysis(
            out_dir=args.out_dir,
            n_values=args.N_values,
            off_pcts=args.off_pcts,
            profiles=args.data_profiles,
            generator_seeds=args.generator_seeds,
            samples=args.samples,
            L=args.L,
            sigma=args.sigma,
            P=args.P,
            exact_source_dir=args.exact_source_dir,
            exact_time_limit=args.exact_time_limit,
            docs_path=docs_path,
        )
        _print_existing_outputs(
            args.out_dir,
            (
                "csv_data.tar.gz",
                "local_threshold_real_off_report.md",
                "real_off_raw_u_g_cdf_by_off_pct.png",
                "real_off_fraction_exact_cdf.png",
                "real_off_mean_fraction_by_off_pct.png",
                "real_off_exact_recovery_by_off_pct.png",
                "real_off_best_cyclic_start_hist.png",
                "real_off_failure_diagnostics.png",
                "real_off_runtime_by_method.png",
            ),
            docs_path=docs_path,
        )
        return True

    if args.threshold_local_exact_analysis:
        docs_path = Path("docs/local_threshold_exact_gauss_report.md")
        run_local_threshold_exact_analysis(
            exact_dir=args.exact_source_dir,
            out_dir=args.out_dir,
            docs_path=docs_path,
            n_values=args.N_values,
            k_pcts=args.K_pcts,
            profiles=args.data_profiles,
        )
        _print_existing_outputs(
            args.out_dir,
            (
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
            ),
            docs_path=docs_path,
        )
        return True

    return False


def main():
    args = parse_args()
    _ensure_single_experiment_mode(args)
    args.off_cases = build_off_cases(args)
    if args.cyclic_best_3swap_analysis:
        args.out_dir = args.out_dir or default_cyclic_3swap_out_dir(args)
        run_cyclic_best_3swap_analysis(args)
        return

    if args.ug_swap_seed_comparison:
        args.out_dir = args.out_dir or default_ug_swap_out_dir(args)
        run_ug_swap_seed_comparison(args)
        return

    args.out_dir = args.out_dir or default_out_dir(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.unified_local_swap_comparison:
        run_unified_local_swap_comparison(args)
        return

    if _run_threshold_analysis_mode(args):
        return

    all_algorithms = (
        REQUESTED_GEN_SOLVERS if args.solver_set == "requested-gen" else CDF_SOLVERS
    )
    algorithms = select_algorithms(all_algorithms, args)
    if args.workers > 1 and args.solver_set != "requested-gen":
        raise ValueError("--workers > 1 is supported only with --solver-set requested-gen.")

    runs_path = args.out_dir / "cdf_runs.csv"
    if args.plot_only:
        if not runs_path.exists():
            raise FileNotFoundError(f"No existing runs file: {runs_path}")
        runs = pd.read_csv(runs_path)
    else:
        if args.workers > 1:
            runs = run_benchmark_parallel(args, algorithms, runs_path)
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
        "pareto_speed_accuracy.png",
        "cdf_u_g_db_h3_submodular_gen.png",
        "cdf_runtime_seconds_h3_submodular_gen.png",
    ):
        output_path = args.out_dir / output_name
        if output_path.exists():
            print(f"  {output_path}", flush=True)


if __name__ == "__main__":
    main()
