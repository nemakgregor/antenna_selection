import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None

from utils.solver_sets import CDF_SOLVERS, REQUESTED_GEN_SOLVERS
from utils.data import generate_v_from_rng
from utils.evaluation import evaluate_solver
from utils.io import atomic_write_csv
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
    "H3ThresholdT123-Gen",
    "CapWindow-Gen",
    "CapWindowFull-Gen",
    "CapSubmod-Gen",
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
    parser.add_argument("--out-dir", type=Path, default=None)
    return parser.parse_args()


def default_out_dir(args):
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


def main():
    args = parse_args()
    args.off_cases = build_off_cases(args)
    args.out_dir = args.out_dir or default_out_dir(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)

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
