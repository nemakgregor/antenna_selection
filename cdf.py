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

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None

from algorithms import (
    calculate_objectives,
    check_constraints,
    solve_coutino_greedy,
    solve_frame_portfolio,
    solve_h1,
    solve_h2,
    solve_h3,
    solve_h3_fast,
    solve_h3_strong_weak,
    solve_miso_energy_greedy,
    solve_pareto_interference_greedy,
)


FRAME_FAST_KWARGS = {
    "max_refined_starts": 3,
    "max_passes": 2,
    "remove_limit": 60,
    "add_limit": 60,
    "lambdas": (),
}

DEFAULT_ALGORITHMS = (
    "H1",
    "H2",
    "H3",
    "FrameOnly-Gen",
    "S-threshold-Gen",
    "Coutino",
)

OUR_ALGORITHMS = ("FrameOnly-Gen", "S-threshold-Gen", "Coutino")
FOCUSED_H3_FRAMEONLY = ("H3", "FrameOnly-Gen")

ALGORITHM_COLORS = {
    "H1": "#0072B2",
    "H2": "#D55E00",
    "H3": "#009E73",
    "FrameOnly-Gen": "#6A3D9A",
    "Frame-Gen": "#CC79A7",
    "S-threshold-Gen": "#E69F00",
    "Coutino": "#000000",
}


def generate_v_from_rng(rng, N, L):
    V = rng.normal(size=(N, L)) + 1j * rng.normal(size=(N, L))
    column_norms = np.linalg.norm(V, axis=0)
    V /= column_norms
    antenna_max = np.max(np.linalg.norm(V, axis=1))
    V /= antenna_max
    return V


def build_algorithms():
    return (
        ("H1", lambda V, K, sigma, P, random_state: solve_h1(V, K, sigma=sigma, P=P)),
        ("H2", lambda V, K, sigma, P, random_state: solve_h2(V, K, sigma=sigma, P=P)),
        (
            "H3",
            lambda V, K, sigma, P, random_state: solve_h3_strong_weak(
                V, K, sigma=sigma, P=P
            ),
        ),
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
            "S-threshold-BF",
            lambda V, K, sigma, P, random_state: solve_h3(
                V, K, target_obj="bf", sigma=sigma, P=P
            ),
        ),
        (
            "S-threshold-Int",
            lambda V, K, sigma, P, random_state: solve_h3(
                V, K, target_obj="int", sigma=sigma, P=P
            ),
        ),
        (
            "S-threshold-Gen",
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
            "FrameOnly-BF",
            lambda V, K, sigma, P, random_state: solve_frame_portfolio(
                V,
                K,
                target_obj="bf",
                sigma=sigma,
                P=P,
                random_state=random_state,
                external_starts=False,
                **FRAME_FAST_KWARGS,
            ),
        ),
        (
            "FrameOnly-Int",
            lambda V, K, sigma, P, random_state: solve_frame_portfolio(
                V,
                K,
                target_obj="int",
                sigma=sigma,
                P=P,
                random_state=random_state,
                external_starts=False,
                **FRAME_FAST_KWARGS,
            ),
        ),
        (
            "FrameOnly-Gen",
            lambda V, K, sigma, P, random_state: solve_frame_portfolio(
                V,
                K,
                target_obj="gen",
                sigma=sigma,
                P=P,
                random_state=random_state,
                external_starts=False,
                **FRAME_FAST_KWARGS,
            ),
        ),
        (
            "N-H3-Fast",
            lambda V, K, sigma, P, random_state: solve_h3_fast(
                V, K, random_state=random_state
            ),
        ),
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
    parser.add_argument("--sigma", type=float, default=1.0)
    parser.add_argument("--P", type=float, default=1.0)
    parser.add_argument(
        "--algorithms",
        nargs="+",
        default=None,
        help="Optional subset of algorithm names. Default: H1 H2 H3 plus three best distinct families.",
    )
    parser.add_argument(
        "--all-algorithms",
        action="store_true",
        help="Run and plot every algorithm instead of the compact default subset.",
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
    return parser.parse_args()


def default_out_dir(args):
    if args.off_counts is not None:
        off_label = "count" + "_".join(str(value) for value in args.off_counts)
    else:
        off_label = "_".join(format_number(value) for value in args.off_pcts)
    seed_label = "_".join(str(value) for value in args.generator_seeds)
    return Path(
        f"results/cdf_N{args.N}_L{args.L}_off{off_label}_"
        f"seeds{seed_label}_{args.samples}samples"
    )


def build_off_cases(args):
    cases = []
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


def format_number(value):
    if float(value).is_integer():
        return str(int(value))
    return str(value).replace(".", "p")


def atomic_write_csv(df, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(path.name + ".tmp")
    df.to_csv(tmp_path, index=False)
    os.replace(tmp_path, path)


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
    if args.all_algorithms:
        return all_algorithms

    selected_names = args.algorithms or DEFAULT_ALGORITHMS
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
    started_at = time.perf_counter()
    with np.errstate(all="ignore"):
        x = solver(V, K, args.sigma, args.P, random_state)
        elapsed_seconds = time.perf_counter() - started_at
        valid, active_count = check_constraints(x, K)
        u_bf, u_i, u_g = calculate_objectives(V, x, sigma=args.sigma, P=args.P)

    if not valid or not np.isfinite([u_bf, u_i, u_g]).all():
        raise RuntimeError(f"Invalid result for {name}, off_pct={off_pct}, K={K}.")

    u_g_safe = max(float(u_g), np.finfo(float).tiny)
    return {
        "algorithm": name,
        "active_count": int(active_count),
        "u_bf": float(u_bf),
        "u_i": float(u_i),
        "u_g": float(u_g),
        "u_g_db": float(10.0 * np.log10(u_g_safe)),
        "elapsed_seconds": float(elapsed_seconds),
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


def empirical_cdf(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return values, values
    values = np.sort(values)
    probs = np.arange(1, len(values) + 1, dtype=float) / len(values)
    return np.r_[values[0], values], np.r_[0.0, probs]


def plot_cdf(runs, algorithms, value_col, ylabel, title, out_path, log_y=False):
    off_pcts = sorted(runs["off_pct"].unique())
    generator_seeds = sorted(runs["generator_seed"].unique())
    fig, axes = plt.subplots(
        len(generator_seeds),
        len(off_pcts),
        figsize=(7.0 * len(off_pcts), 4.3 * len(generator_seeds)),
        sharey=True,
        squeeze=False,
    )
    color_map = plt.get_cmap("tab20")
    colors = {
        name: ALGORITHM_COLORS.get(name, color_map(index % color_map.N))
        for index, (name, _) in enumerate(algorithms)
    }

    handles = []
    labels = []
    for row, generator_seed in enumerate(generator_seeds):
        for col, off_pct in enumerate(off_pcts):
            ax = axes[row, col]
            data = runs[
                (runs["generator_seed"] == generator_seed)
                & (runs["off_pct"] == off_pct)
            ]
            if data.empty:
                ax.set_visible(False)
                continue

            K = int(data["K"].iloc[0])
            for name, _ in algorithms:
                values = data[data["algorithm"] == name][value_col]
                x_values, y_values = empirical_cdf(values)
                if len(x_values) == 0:
                    continue
                line = ax.step(
                    y_values,
                    x_values,
                    where="post",
                    linewidth=1.8,
                    color=colors[name],
                    marker=".",
                    markersize=2.0,
                    label=name,
                )[0]
                if row == 0 and col == 0:
                    handles.append(line)
                    labels.append(name)

            ax.set_title(f"seed={int(generator_seed)}, {off_pct:g}% off, K={K}")
            ax.set_xlabel("Cumulative fraction of examples")
            if col == 0:
                ax.set_ylabel(ylabel)
            ax.set_xlim(0.0, 1.02)
            ax.grid(True, alpha=0.25)
            if log_y:
                ax.set_yscale("log")

    fig.suptitle(title)
    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.995, 0.5), fontsize=8)
    fig.tight_layout(rect=(0.0, 0.0, 0.84, 0.93))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_outputs(runs, algorithms, out_dir):
    selected = {name for name, _ in algorithms}
    runs = runs[runs["algorithm"].isin(selected)].copy()
    if runs.empty:
        raise RuntimeError("No runs are available for the selected algorithms.")

    summary = build_summary(runs)
    atomic_write_csv(summary, out_dir / "cdf_summary.csv")
    improvement = build_baseline_improvement(runs)
    atomic_write_csv(improvement, out_dir / "cdf_baseline_improvement.csv")
    write_improvement_report(improvement, out_dir / "cdf_baseline_improvement.md")
    write_our_vs_h123_report(improvement, out_dir)
    plot_cdf(
        runs,
        algorithms,
        "u_g_db",
        "10 lg(U_G), dB",
        "Cumulative distribution of general objective U_G",
        out_dir / "cdf_u_g_db.png",
    )
    plot_cdf(
        runs,
        algorithms,
        "elapsed_seconds",
        "Elapsed time, seconds",
        "Cumulative distribution of solver runtime",
        out_dir / "cdf_runtime_seconds.png",
        log_y=True,
    )
    if set(FOCUSED_H3_FRAMEONLY).issubset(selected):
        focused_algorithms = tuple(
            (name, solver)
            for name, solver in algorithms
            if name in FOCUSED_H3_FRAMEONLY
        )
        plot_cdf(
            runs,
            focused_algorithms,
            "u_g_db",
            "10 lg(U_G), dB",
            "Focused cumulative distribution: H3 vs FrameOnly-Gen",
            out_dir / "cdf_u_g_db_h3_frameonly.png",
        )
        plot_cdf(
            runs,
            focused_algorithms,
            "elapsed_seconds",
            "Elapsed time, seconds",
            "Focused solver runtime: H3 vs FrameOnly-Gen",
            out_dir / "cdf_runtime_seconds_h3_frameonly.png",
            log_y=True,
        )


def main():
    args = parse_args()
    args.off_cases = build_off_cases(args)
    args.out_dir = args.out_dir or default_out_dir(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    algorithms = select_algorithms(build_algorithms(), args)

    runs_path = args.out_dir / "cdf_runs.csv"
    if args.plot_only:
        if not runs_path.exists():
            raise FileNotFoundError(f"No existing runs file: {runs_path}")
        runs = pd.read_csv(runs_path)
    else:
        runs = run_benchmark(args, algorithms, runs_path)
    write_outputs(runs, algorithms, args.out_dir)
    print("Wrote:", flush=True)
    print(f"  {args.out_dir / 'cdf_runs.csv'}", flush=True)
    print(f"  {args.out_dir / 'cdf_summary.csv'}", flush=True)
    print(f"  {args.out_dir / 'cdf_baseline_improvement.csv'}", flush=True)
    print(f"  {args.out_dir / 'cdf_baseline_improvement.md'}", flush=True)
    print(f"  {args.out_dir / 'cdf_our_vs_h123.csv'}", flush=True)
    print(f"  {args.out_dir / 'cdf_our_vs_h123.md'}", flush=True)
    print(f"  {args.out_dir / 'cdf_u_g_db.png'}", flush=True)
    print(f"  {args.out_dir / 'cdf_runtime_seconds.png'}", flush=True)
    if (args.out_dir / "cdf_u_g_db_h3_frameonly.png").exists():
        print(f"  {args.out_dir / 'cdf_u_g_db_h3_frameonly.png'}", flush=True)
    if (args.out_dir / "cdf_runtime_seconds_h3_frameonly.png").exists():
        print(
            f"  {args.out_dir / 'cdf_runtime_seconds_h3_frameonly.png'}",
            flush=True,
        )


if __name__ == "__main__":
    main()
