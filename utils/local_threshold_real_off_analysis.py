import math
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

from algorithms.h3_strong_weak import solve_h3_strong_weak
from algorithms.h3_threshold_local import (
    best_cyclic_threshold_window,
    refine_selection_by_swaps,
    threshold_window_selection,
)
from utils.brute_force import (
    brute_force_exact_u_g,
    contiguous_threshold_window_T,
    parse_subset_string,
    subset_to_string,
)
from utils.data import generate_v_profile_from_rng
from utils.io import archive_csv_artifacts, atomic_write_csv


REAL_OFF_SEED_RULES = ("best_cyclic_window", "T_0p05N", "strong_weak")
REAL_OFF_MAX_SWAPS = (0, 1, 2)

REAL_OFF_LABELS = {
    "exact_best": "exact best",
    "best_cyclic_window_s0": "cyclic best",
    "best_cyclic_window_s1": "cyclic best + 1 swap",
    "best_cyclic_window_s2": "cyclic best + 2 swaps",
    "T_0p05N_s0": "T=0.05N",
    "T_0p05N_s1": "T=0.05N + 1 swap",
    "T_0p05N_s2": "T=0.05N + 2 swaps",
    "strong_weak_s0": "strong/weak",
    "strong_weak_s1": "strong/weak + 1 swap",
    "strong_weak_s2": "strong/weak + 2 swaps",
}

REAL_OFF_COLORS = {
    "exact_best": "#000000",
    "best_cyclic_window_s0": "#0072B2",
    "best_cyclic_window_s1": "#009E73",
    "best_cyclic_window_s2": "#56B4E9",
    "T_0p05N_s0": "#E69F00",
    "T_0p05N_s1": "#CC79A7",
    "T_0p05N_s2": "#D55E00",
    "strong_weak_s0": "#6A3D9A",
    "strong_weak_s1": "#8E24AA",
    "strong_weak_s2": "#4A148C",
}


def active_k_from_off_pct(N, off_pct):
    N = int(N)
    off_pct = float(off_pct)
    if not (0.0 <= off_pct < 100.0):
        raise ValueError("off_pct must satisfy 0 <= off_pct < 100.")
    K_off = int(round(N * off_pct / 100.0))
    K_active = N - K_off
    if not (0 <= K_active <= N):
        raise ValueError("off_pct produced invalid K_active.")
    return K_active, K_off


def active_k_from_active_pct(N, active_pct):
    N = int(N)
    active_pct = float(active_pct)
    if not (0.0 <= active_pct <= 100.0):
        raise ValueError("active_pct must satisfy 0 <= active_pct <= 100.")
    K_active = int(round(N * active_pct / 100.0))
    K_off = N - K_active
    if not (0 <= K_active <= N):
        raise ValueError("active_pct produced invalid K_active.")
    return K_active, K_off


def run_active_k_cyclic_local_exact_analysis(
    out_dir,
    n_values,
    active_pcts,
    profiles,
    generator_seeds,
    samples,
    L,
    sigma=1.0,
    P=1.0,
    exact_source_dir=None,
    exact_time_limit=120.0,
    max_swaps_values=REAL_OFF_MAX_SWAPS,
    archive_csv=True,
    docs_path=None,
):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    exact_source = _load_exact_source(exact_source_dir)

    exact_rows = []
    local_rows = []
    total_cases = (
        len(n_values)
        * len(profiles)
        * len(generator_seeds)
        * int(samples)
        * len(active_pcts)
    )
    case_no = 0
    for N in n_values:
        N = int(N)
        for profile in profiles:
            for generator_seed in generator_seeds:
                rng = np.random.RandomState(int(generator_seed))
                for sample in range(int(samples)):
                    V = generate_v_profile_from_rng(rng, N, int(L), str(profile))
                    for requested_active_pct in active_pcts:
                        case_no += 1
                        K_active, K_off = active_k_from_active_pct(
                            N,
                            requested_active_pct,
                        )
                        actual_active_pct = 100.0 * float(K_active) / float(N)
                        actual_off_pct = 100.0 * float(K_off) / float(N)
                        if case_no == 1 or case_no % 100 == 0:
                            print(
                                "active-K local exact "
                                f"{case_no}/{total_cases}: N={N}, L={int(L)}, "
                                f"profile={profile}, seed={int(generator_seed)}, "
                                f"sample={sample}, requested_active_pct={float(requested_active_pct):g}, "
                                f"K_active={K_active}, K_off={K_off}",
                                flush=True,
                            )
                        case_base = {
                            "data_profile": str(profile),
                            "generator_seed": int(generator_seed),
                            "sample": int(sample),
                            "N": N,
                            "L": int(L),
                            "K": int(K_active),
                            "K_active": int(K_active),
                            "K_off": int(K_off),
                            "requested_active_pct": float(requested_active_pct),
                            "off_pct": actual_off_pct,
                            "active_pct": actual_active_pct,
                            "sigma": float(sigma),
                            "P": float(P),
                        }
                        exact = _exact_for_case(
                            V,
                            case_base,
                            exact_source,
                            time_limit_seconds=exact_time_limit,
                        )
                        exact_rows.append({**case_base, **exact})
                        local_rows.extend(
                            _local_rows_for_case(
                                V,
                                case_base,
                                exact,
                                max_swaps_values=max_swaps_values,
                            )
                        )

    exact_runs = pd.DataFrame(exact_rows)
    local_runs = pd.DataFrame(local_rows)
    summary = build_real_off_local_summary(local_runs)
    diagnostics = build_real_off_local_diagnostics(local_runs)
    failure_cases = build_real_off_failure_cases(local_runs, diagnostics)

    atomic_write_csv(exact_runs, out_dir / "local_threshold_active_k_exact_runs.csv")
    atomic_write_csv(local_runs, out_dir / "local_threshold_active_k_runs.csv")
    atomic_write_csv(summary, out_dir / "local_threshold_active_k_summary.csv")
    atomic_write_csv(
        failure_cases,
        out_dir / "local_threshold_active_k_failure_cases.csv",
    )
    atomic_write_csv(
        diagnostics,
        out_dir / "local_threshold_active_k_diagnostics.csv",
    )

    write_active_k_local_plots(exact_runs, local_runs, summary, diagnostics, out_dir)
    write_active_k_local_report(
        exact_runs,
        local_runs,
        summary,
        diagnostics,
        failure_cases,
        out_dir,
        exact_source_dir=exact_source_dir,
        docs_path=docs_path,
    )
    archive_path = None
    if archive_csv:
        archive_path = archive_csv_artifacts(out_dir, remove_originals=True)

    return {
        "exact_runs": exact_runs,
        "local_runs": local_runs,
        "summary": summary,
        "diagnostics": diagnostics,
        "failure_cases": failure_cases,
        "archive_path": archive_path,
    }


def run_real_off_cyclic_local_exact_analysis(
    out_dir,
    n_values,
    off_pcts,
    profiles,
    generator_seeds,
    samples,
    L,
    sigma=1.0,
    P=1.0,
    exact_source_dir=None,
    exact_time_limit=120.0,
    max_swaps_values=REAL_OFF_MAX_SWAPS,
    archive_csv=True,
    docs_path=None,
):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    exact_source = _load_exact_source(exact_source_dir)

    exact_rows = []
    local_rows = []
    total_cases = (
        len(n_values)
        * len(profiles)
        * len(generator_seeds)
        * int(samples)
        * len(off_pcts)
    )
    case_no = 0
    for N in n_values:
        N = int(N)
        for profile in profiles:
            for generator_seed in generator_seeds:
                rng = np.random.RandomState(int(generator_seed))
                for sample in range(int(samples)):
                    V = generate_v_profile_from_rng(rng, N, int(L), str(profile))
                    for off_pct in off_pcts:
                        case_no += 1
                        K_active, K_off = active_k_from_off_pct(N, off_pct)
                        if case_no == 1 or case_no % 100 == 0:
                            print(
                                "real-off local exact "
                                f"{case_no}/{total_cases}: N={N}, L={int(L)}, "
                                f"profile={profile}, seed={int(generator_seed)}, "
                                f"sample={sample}, off_pct={float(off_pct):g}, "
                                f"K_active={K_active}, K_off={K_off}",
                                flush=True,
                            )
                        case_base = {
                            "data_profile": str(profile),
                            "generator_seed": int(generator_seed),
                            "sample": int(sample),
                            "N": N,
                            "L": int(L),
                            "K": int(K_active),
                            "K_active": int(K_active),
                            "K_off": int(K_off),
                            "off_pct": float(off_pct),
                            "active_pct": 100.0 * float(K_active) / float(N),
                            "sigma": float(sigma),
                            "P": float(P),
                        }
                        exact = _exact_for_case(
                            V,
                            case_base,
                            exact_source,
                            time_limit_seconds=exact_time_limit,
                        )
                        exact_rows.append({**case_base, **exact})
                        local_rows.extend(
                            _local_rows_for_case(
                                V,
                                case_base,
                                exact,
                                max_swaps_values=max_swaps_values,
                            )
                        )

    exact_runs = pd.DataFrame(exact_rows)
    local_runs = pd.DataFrame(local_rows)
    summary = build_real_off_local_summary(local_runs)
    diagnostics = build_real_off_local_diagnostics(local_runs)
    failure_cases = build_real_off_failure_cases(local_runs, diagnostics)

    atomic_write_csv(exact_runs, out_dir / "local_threshold_real_off_exact_runs.csv")
    atomic_write_csv(local_runs, out_dir / "local_threshold_real_off_runs.csv")
    atomic_write_csv(summary, out_dir / "local_threshold_real_off_summary.csv")
    atomic_write_csv(
        failure_cases,
        out_dir / "local_threshold_real_off_failure_cases.csv",
    )
    atomic_write_csv(
        diagnostics,
        out_dir / "local_threshold_real_off_diagnostics.csv",
    )

    write_real_off_local_plots(exact_runs, local_runs, summary, diagnostics, out_dir)
    write_real_off_local_report(
        exact_runs,
        local_runs,
        summary,
        diagnostics,
        failure_cases,
        out_dir,
        exact_source_dir=exact_source_dir,
        docs_path=docs_path,
    )
    archive_path = None
    if archive_csv:
        archive_path = archive_csv_artifacts(out_dir, remove_originals=True)

    return {
        "exact_runs": exact_runs,
        "local_runs": local_runs,
        "summary": summary,
        "diagnostics": diagnostics,
        "failure_cases": failure_cases,
        "archive_path": archive_path,
    }


def build_real_off_local_summary(local_runs):
    if local_runs.empty:
        return pd.DataFrame()

    def q(value):
        return lambda data: data.quantile(value)

    group_cols = [
        "data_profile",
        "N",
        "L",
        "K",
        "K_active",
        "K_off",
        "off_pct",
        "active_pct",
        "seed_rule",
        "max_swaps",
    ]
    if "requested_active_pct" in local_runs.columns:
        group_cols.insert(8, "requested_active_pct")

    return (
        local_runs.groupby(
            group_cols,
            as_index=False,
        )
        .agg(
            cases=("u_g", "count"),
            seed_position_mean=("seed_position", "mean"),
            seed_position_p50=("seed_position", q(0.50)),
            u_g_mean=("u_g", "mean"),
            u_g_p05=("u_g", q(0.05)),
            u_g_p50=("u_g", q(0.50)),
            u_g_p95=("u_g", q(0.95)),
            fraction_exact_mean=("fraction_exact_u_g", "mean"),
            fraction_exact_p05=("fraction_exact_u_g", q(0.05)),
            fraction_exact_p50=("fraction_exact_u_g", q(0.50)),
            fraction_exact_p95=("fraction_exact_u_g", q(0.95)),
            gap_to_exact_pct_mean=("gap_to_exact_pct", "mean"),
            exact_rate=("is_exact_best", "mean"),
            near_99_rate=("near_99", "mean"),
            swaps_applied_mean=("swaps_applied", "mean"),
            add_candidate_count_mean=("add_candidate_count", "mean"),
            evaluated_swap_count_mean=("evaluated_swap_count", "mean"),
            elapsed_seconds_mean=("elapsed_seconds", "mean"),
            elapsed_seconds_p95=("elapsed_seconds", q(0.95)),
        )
        .sort_values(["data_profile", "N", "off_pct", "seed_rule", "max_swaps"])
    )


def build_real_off_local_diagnostics(local_runs):
    rows = []
    for _, row in local_runs.iterrows():
        exact_subset = parse_subset_string(row["exact_subset"])
        local_subset = parse_subset_string(row["subset"])
        exact_set = set(exact_subset)
        local_set = set(local_subset)
        K = int(row["K"])
        overlap = len(exact_set & local_set)
        exact_features = _features_from_serialized(row, "exact")
        local_features = _features_from_serialized(row, "local")
        rows.append(
            {
                **_case_columns(row),
                "seed_rule": str(row["seed_rule"]),
                "max_swaps": int(row["max_swaps"]),
                "seed_position": int(row["seed_position"]),
                "fraction_exact_u_g": float(row["fraction_exact_u_g"]),
                "gap_to_exact_pct": float(row["gap_to_exact_pct"]),
                "is_exact_best": bool(row["is_exact_best"]),
                "overlap_count": int(overlap),
                "overlap_fraction": overlap / max(K, 1),
                "swap_distance_to_exact": int(K - overlap),
                "exact_max_row_power": exact_features["max_row_power"],
                "local_max_row_power": local_features["max_row_power"],
                "local_exact_max_power_ratio": _safe_ratio(
                    local_features["max_row_power"],
                    exact_features["max_row_power"],
                ),
                "exact_scaled_trace": exact_features["scaled_trace"],
                "local_scaled_trace": local_features["scaled_trace"],
                "local_exact_scaled_trace_ratio": _safe_ratio(
                    local_features["scaled_trace"],
                    exact_features["scaled_trace"],
                ),
                "exact_eig_min": exact_features["eig_min"],
                "local_eig_min": local_features["eig_min"],
                "local_exact_eig_min_ratio": _safe_ratio(
                    local_features["eig_min"],
                    exact_features["eig_min"],
                ),
                "exact_eig_balance": exact_features["eig_balance"],
                "local_eig_balance": local_features["eig_balance"],
                "local_exact_eig_balance_delta": (
                    local_features["eig_balance"] - exact_features["eig_balance"]
                ),
            }
        )
    return pd.DataFrame(rows)


def build_real_off_failure_cases(local_runs, diagnostics):
    if local_runs.empty or diagnostics.empty:
        return pd.DataFrame()
    primary = diagnostics[
        (diagnostics["seed_rule"] == "best_cyclic_window")
        & (diagnostics["max_swaps"] == 2)
        & (diagnostics["fraction_exact_u_g"] < 1.0 - 1e-9)
    ].copy()
    if primary.empty:
        return primary
    local_cols = [
        "data_profile",
        "generator_seed",
        "sample",
        "N",
        "L",
        "K",
        "K_active",
        "K_off",
        *([] if "requested_active_pct" not in local_runs.columns else ["requested_active_pct"]),
        "off_pct",
        "active_pct",
        "seed_rule",
        "max_swaps",
        "seed_position",
        "u_g",
        "exact_u_g",
        "fraction_exact_u_g",
        "subset",
        "exact_subset",
        "swap_history",
        "swaps_applied",
    ]
    merged = primary.merge(
        local_runs[local_cols],
        on=[
            "data_profile",
            "generator_seed",
            "sample",
            "N",
            "L",
            "K",
            "K_active",
            "K_off",
            *([] if "requested_active_pct" not in local_runs.columns else ["requested_active_pct"]),
            "off_pct",
            "active_pct",
            "seed_rule",
            "max_swaps",
            "seed_position",
            "fraction_exact_u_g",
        ],
        how="left",
    )
    return merged.sort_values("fraction_exact_u_g").head(100)


def write_real_off_local_plots(exact_runs, local_runs, summary, diagnostics, out_dir):
    if exact_runs.empty or local_runs.empty:
        return

    from utils.plotting import use_agg_backend

    use_agg_backend()
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    off_pcts = sorted(exact_runs["off_pct"].unique())
    rule_specs = _rule_specs()

    cols = min(2, max(1, len(off_pcts)))
    rows = int(math.ceil(len(off_pcts) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(6.2 * cols, 4.2 * rows), squeeze=False)
    for ax, off_pct in zip(axes.ravel(), off_pcts):
        exact_chunk = exact_runs[exact_runs["off_pct"] == off_pct]
        values, cumulative = _empirical_cdf(exact_chunk["exact_u_g"])
        if len(values):
            ax.step(
                cumulative,
                values,
                where="post",
                color=REAL_OFF_COLORS["exact_best"],
                label=REAL_OFF_LABELS["exact_best"],
            )
        local_pct = local_runs[local_runs["off_pct"] == off_pct]
        for rule_key, seed_rule, max_swaps in rule_specs:
            chunk = local_pct[
                (local_pct["seed_rule"] == seed_rule)
                & (local_pct["max_swaps"] == max_swaps)
            ]
            values, cumulative = _empirical_cdf(chunk["u_g"])
            if len(values):
                ax.step(
                    cumulative,
                    values,
                    where="post",
                    color=REAL_OFF_COLORS.get(rule_key),
                    label=REAL_OFF_LABELS[rule_key],
                    linewidth=1.1,
                )
        active_mean = exact_chunk["active_pct"].astype(float).mean()
        k_mean = exact_chunk["K_active"].astype(float).mean()
        off_mean = exact_chunk["K_off"].astype(float).mean()
        ax.set_title(
            f"{off_pct:g}% off: mean K_active={k_mean:.1f}, K_off={off_mean:.1f}, "
            f"active={active_mean:.1f}%"
        )
        ax.set_xlabel("cumulative fraction")
        ax.set_ylabel("raw U_G")
        ax.set_yscale("log")
        ax.grid(True, alpha=0.25)
    for ax in axes.ravel()[len(off_pcts) :]:
        ax.axis("off")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Real-off local search raw U_G CDF")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_dir / "real_off_raw_u_g_cdf_by_off_pct.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    for rule_key, seed_rule, max_swaps in rule_specs:
        chunk = local_runs[
            (local_runs["seed_rule"] == seed_rule)
            & (local_runs["max_swaps"] == max_swaps)
        ]
        values, cumulative = _empirical_cdf(chunk["fraction_exact_u_g"])
        if len(values):
            ax.step(
                cumulative,
                values,
                where="post",
                color=REAL_OFF_COLORS.get(rule_key),
                label=REAL_OFF_LABELS[rule_key],
                linewidth=1.2,
            )
    ax.set_xlabel("cumulative fraction")
    ax.set_ylabel("U_G(method) / exact U_G")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.03)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7)
    ax.set_title("Real-off local search fraction of exact U_G")
    fig.tight_layout()
    fig.savefig(out_dir / "real_off_fraction_exact_cdf.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    global_by_off = _global_by_off_pct(summary)
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for rule_key, seed_rule, max_swaps in rule_specs:
        chunk = global_by_off[
            (global_by_off["seed_rule"] == seed_rule)
            & (global_by_off["max_swaps"] == max_swaps)
        ].sort_values("off_pct")
        if chunk.empty:
            continue
        ax.plot(
            chunk["off_pct"],
            chunk["fraction_exact_mean"],
            marker="o",
            color=REAL_OFF_COLORS.get(rule_key),
            label=REAL_OFF_LABELS[rule_key],
        )
    ax.set_xlabel("percent of antennas turned off")
    ax.set_ylabel("mean U_G / exact U_G")
    ax.set_ylim(0.0, 1.03)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7)
    ax.set_title("Mean quality by real off percentage")
    fig.tight_layout()
    fig.savefig(out_dir / "real_off_mean_fraction_by_off_pct.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for rule_key, seed_rule, max_swaps in rule_specs:
        chunk = global_by_off[
            (global_by_off["seed_rule"] == seed_rule)
            & (global_by_off["max_swaps"] == max_swaps)
        ].sort_values("off_pct")
        if chunk.empty:
            continue
        ax.plot(
            chunk["off_pct"],
            chunk["exact_rate"],
            marker="o",
            color=REAL_OFF_COLORS.get(rule_key),
            label=REAL_OFF_LABELS[rule_key],
        )
    ax.set_xlabel("percent of antennas turned off")
    ax.set_ylabel("exact recovery rate")
    ax.set_ylim(0.0, 1.03)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7)
    ax.set_title("Exact recovery by real off percentage")
    fig.tight_layout()
    fig.savefig(out_dir / "real_off_exact_recovery_by_off_pct.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    cyclic = local_runs[
        (local_runs["seed_rule"] == "best_cyclic_window")
        & (local_runs["max_swaps"] == 0)
    ].copy()
    if not cyclic.empty:
        fig, axes = plt.subplots(1, len(off_pcts), figsize=(5.2 * len(off_pcts), 4.0), squeeze=False)
        for ax, off_pct in zip(axes.ravel(), off_pcts):
            chunk = cyclic[cyclic["off_pct"] == off_pct]
            ax.hist(chunk["seed_position"], bins=range(0, int(chunk["N"].max()) + 2), color="#0072B2", alpha=0.85)
            ax.set_title(f"{off_pct:g}% off")
            ax.set_xlabel("best cyclic start")
            ax.set_ylabel("count")
            ax.grid(True, alpha=0.25)
        fig.suptitle("Best tested cyclic start distribution")
        fig.tight_layout(rect=(0, 0, 1, 0.93))
        fig.savefig(out_dir / "real_off_best_cyclic_start_hist.png", dpi=180, bbox_inches="tight")
        plt.close(fig)

    _plot_real_off_failure_diagnostics(diagnostics, out_dir)
    _plot_real_off_runtime(summary, out_dir)


def write_active_k_local_plots(exact_runs, local_runs, summary, diagnostics, out_dir):
    if exact_runs.empty or local_runs.empty:
        return

    from utils.plotting import use_agg_backend

    use_agg_backend()
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    active_pcts = sorted(exact_runs["requested_active_pct"].unique())
    rule_specs = _rule_specs()

    cols = 4
    rows = int(math.ceil(len(active_pcts) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(5.0 * cols, 3.7 * rows), squeeze=False)
    for ax, active_pct in zip(axes.ravel(), active_pcts):
        exact_chunk = exact_runs[exact_runs["requested_active_pct"] == active_pct]
        values, cumulative = _empirical_cdf(exact_chunk["exact_u_g"])
        if len(values):
            ax.step(cumulative, values, where="post", color=REAL_OFF_COLORS["exact_best"], label=REAL_OFF_LABELS["exact_best"])
        local_pct = local_runs[local_runs["requested_active_pct"] == active_pct]
        for rule_key, seed_rule, max_swaps in rule_specs:
            chunk = local_pct[
                (local_pct["seed_rule"] == seed_rule)
                & (local_pct["max_swaps"] == max_swaps)
            ]
            values, cumulative = _empirical_cdf(chunk["u_g"])
            if len(values):
                ax.step(
                    cumulative,
                    values,
                    where="post",
                    color=REAL_OFF_COLORS.get(rule_key),
                    label=REAL_OFF_LABELS[rule_key],
                    linewidth=1.35 if max_swaps == 2 else 1.0,
                )
        actual_mean = exact_chunk["active_pct"].astype(float).mean()
        k_mean = exact_chunk["K_active"].astype(float).mean()
        ax.set_title(f"{active_pct:g}% active: actual {actual_mean:.1f}%, K mean {k_mean:.1f}")
        ax.set_xlabel("cumulative fraction")
        ax.set_ylabel("raw U_G")
        ax.set_yscale("log")
        ax.grid(True, alpha=0.25)
    for ax in axes.ravel()[len(active_pcts) :]:
        ax.axis("off")
    axes[0, 0].legend(fontsize=6)
    fig.suptitle("Active-K sweep raw U_G CDF")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_dir / "active_k_raw_u_g_cdf_by_requested_active_pct.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    for rule_key, seed_rule, max_swaps in rule_specs:
        chunk = local_runs[
            (local_runs["seed_rule"] == seed_rule)
            & (local_runs["max_swaps"] == max_swaps)
        ]
        values, cumulative = _empirical_cdf(chunk["fraction_exact_u_g"])
        if len(values):
            ax.step(
                cumulative,
                values,
                where="post",
                color=REAL_OFF_COLORS.get(rule_key),
                label=REAL_OFF_LABELS[rule_key],
                linewidth=1.2,
            )
    ax.set_xlabel("cumulative fraction")
    ax.set_ylabel("U_G(method) / exact U_G")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.03)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7)
    ax.set_title("Active-K sweep fraction of exact U_G")
    fig.tight_layout()
    fig.savefig(out_dir / "active_k_fraction_exact_cdf.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    by_active = _global_by_requested_active_pct(summary)
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for rule_key, seed_rule, max_swaps in rule_specs:
        chunk = by_active[
            (by_active["seed_rule"] == seed_rule)
            & (by_active["max_swaps"] == max_swaps)
        ].sort_values("requested_active_pct")
        if chunk.empty:
            continue
        ax.plot(
            chunk["requested_active_pct"],
            chunk["fraction_exact_mean"],
            marker="o",
            color=REAL_OFF_COLORS.get(rule_key),
            label=REAL_OFF_LABELS[rule_key],
        )
    ax.set_xlabel("requested active K percentage")
    ax.set_ylabel("mean U_G / exact U_G")
    ax.set_ylim(0.0, 1.03)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7)
    ax.set_title("Mean quality by requested active K percentage")
    fig.tight_layout()
    fig.savefig(out_dir / "active_k_mean_fraction_by_requested_active_pct.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for rule_key, seed_rule, max_swaps in rule_specs:
        chunk = by_active[
            (by_active["seed_rule"] == seed_rule)
            & (by_active["max_swaps"] == max_swaps)
        ].sort_values("requested_active_pct")
        if chunk.empty:
            continue
        ax.plot(
            chunk["requested_active_pct"],
            chunk["exact_rate"],
            marker="o",
            color=REAL_OFF_COLORS.get(rule_key),
            label=REAL_OFF_LABELS[rule_key],
        )
    ax.set_xlabel("requested active K percentage")
    ax.set_ylabel("exact recovery rate")
    ax.set_ylim(0.0, 1.03)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7)
    ax.set_title("Exact recovery by requested active K percentage")
    fig.tight_layout()
    fig.savefig(out_dir / "active_k_exact_recovery_by_requested_active_pct.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    _plot_active_k_best_cyclic_t(local_runs, out_dir)
    _plot_active_k_runtime(summary, out_dir)
    _plot_real_off_failure_diagnostics(diagnostics, out_dir)


def write_real_off_local_report(
    exact_runs,
    local_runs,
    summary,
    diagnostics,
    failure_cases,
    out_dir,
    exact_source_dir=None,
    docs_path=None,
):
    global_summary = _global_summary(local_runs)
    lines = [
        "# Real-Off Cyclic Threshold Local Search Study",
        "",
        "This report uses the real task condition: turn off a requested percentage of antennas while keeping the rest active.",
        "",
        "## K Semantics",
        "",
        "- `K_off = round(N * off_pct / 100)` is the number of disabled antennas.",
        "- `K_active = N - K_off` is the number of selected/kept antennas.",
        "- The solver variable `K` equals `K_active`; it does not mean the number of antennas turned off.",
        "- Therefore `25% off` means `K_active = 0.75N`, and `50% off` means `K_active = 0.50N`.",
        "",
        "## Setup",
        "",
        f"- Exact source: `{exact_source_dir}`" if exact_source_dir is not None else "- Exact source: computed in this run when needed.",
        f"- Profiles: {', '.join(str(value) for value in sorted(local_runs['data_profile'].unique())) if not local_runs.empty else ''}",
        f"- N values: {', '.join(str(int(value)) for value in sorted(local_runs['N'].unique())) if not local_runs.empty else ''}",
        f"- L values: {', '.join(str(int(value)) for value in sorted(local_runs['L'].unique())) if not local_runs.empty else ''}",
        f"- Off percentages: {', '.join(_format_pct(value) for value in sorted(local_runs['off_pct'].unique())) if not local_runs.empty else ''}",
        "- Seeds compared: best tested cyclic window, `T=round(0.05N)`, and strong/weak.",
        "- Local search: greedy remove-one/add-one refinement by raw `U_G`, with 0, 1, or 2 swaps.",
        "",
        "## Local Swap Scheme And Cost",
        "",
        "- A `1-swap` removes one currently active antenna and adds one currently inactive antenna, preserving exact `K_active`.",
        "- One local-search pass evaluates every `(remove, add)` pair from the current active set and the configured add-candidate pool, applies only the single pair with the largest positive `U_G` improvement, and stops if no pair improves `U_G`.",
        "- `2-swap` means two greedy passes, not an exhaustive simultaneous two-pair exchange.",
        "- In this exact cyclic analysis the add-candidate pool is all inactive antennas, so one pass evaluates `K_active * (N - K_active)` swaps.",
        "- With `S` greedy passes, `A` add candidates, and `L` streams, local refinement costs `O(S * K_active * A * L^3)` time after row Gram matrices are built; here `A = N - K_active`, so it is `O(S * K_active * (N - K_active) * L^3)`.",
        "- Seed costs are separate: `T=0.05N` and strong/weak need one sorted window, while best cyclic scans `N` cyclic windows before local refinement.",
        "- Extra working space is `O(N * L^2 + K_active + A + L^2)` for row Gram matrices, active/add sets, and the current Gram matrix.",
        "",
        "## Direct Answer",
        "",
    ]

    if global_summary.empty:
        lines.append("- No rows were produced.")
    else:
        for seed_rule in REAL_OFF_SEED_RULES:
            chunk = global_summary[global_summary["seed_rule"] == seed_rule].sort_values("max_swaps")
            if chunk.empty:
                continue
            for _, row in chunk.iterrows():
                label = REAL_OFF_LABELS[f"{seed_rule}_s{int(row['max_swaps'])}"]
                lines.append(
                    f"- `{label}`: mean fraction exact `{row['fraction_exact_mean']:.4f}`, "
                    f"p05 `{row['fraction_exact_p05']:.4f}`, exact recovery `{row['exact_rate']:.1%}`, "
                    f"mean swaps applied `{row['swaps_applied_mean']:.3f}`."
                )
        cyclic = global_summary[global_summary["seed_rule"] == "best_cyclic_window"]
        zero = cyclic[cyclic["max_swaps"] == 0]
        two = cyclic[cyclic["max_swaps"] == 2]
        if not zero.empty and not two.empty:
            lines.append(
                f"- Two swaps improve cyclic-threshold mean exact fraction by "
                f"`{100.0 * (two.iloc[0]['fraction_exact_mean'] - zero.iloc[0]['fraction_exact_mean']):.3f}` percentage points."
            )

    primary_diag = diagnostics[
        (diagnostics["seed_rule"] == "best_cyclic_window")
        & (diagnostics["max_swaps"] == 2)
    ]
    misses = primary_diag[primary_diag["fraction_exact_u_g"] < 1.0 - 1e-9]
    if not primary_diag.empty:
        if misses.empty:
            lines.append("- After two swaps from cyclic best, every analyzed case recovered exact `U_G`.")
        else:
            lines.append(
                f"- After two swaps from cyclic best, remaining misses have mean overlap "
                f"`{misses['overlap_fraction'].mean():.3f}` with exact and mean swap distance "
                f"`{misses['swap_distance_to_exact'].mean():.3f}` rows."
            )

    lines.extend(
        [
            "",
            "## Global Summary",
            "",
            "| seed | swaps | cases | mean fraction exact | p05 | p50 | p95 | exact rate | near-99 rate | runtime s |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for _, row in global_summary.sort_values(["seed_rule", "max_swaps"]).iterrows():
        label = REAL_OFF_LABELS.get(
            f"{row['seed_rule']}_s{int(row['max_swaps'])}",
            str(row["seed_rule"]),
        )
        lines.append(
            "| "
            + " | ".join(
                [
                    label,
                    str(int(row["max_swaps"])),
                    str(int(row["cases"])),
                    f"{row['fraction_exact_mean']:.4f}",
                    f"{row['fraction_exact_p05']:.4f}",
                    f"{row['fraction_exact_p50']:.4f}",
                    f"{row['fraction_exact_p95']:.4f}",
                    f"{row['exact_rate']:.1%}",
                    f"{row['near_99_rate']:.1%}",
                    _format_float(row["elapsed_seconds_mean"]),
                ]
            )
            + " |"
        )

    by_off = _global_by_off_pct(summary)
    lines.extend(
        [
            "",
            "## By Real Off Percentage",
            "",
            "| off % | K semantics | seed | swaps | mean fraction exact | p05 | exact rate |",
            "|---:|---|---|---:|---:|---:|---:|",
        ]
    )
    for _, row in by_off.sort_values(["off_pct", "seed_rule", "max_swaps"]).iterrows():
        label = REAL_OFF_LABELS.get(
            f"{row['seed_rule']}_s{int(row['max_swaps'])}",
            str(row["seed_rule"]),
        )
        lines.append(
            "| "
            + " | ".join(
                [
                    _format_pct(row["off_pct"]),
                    f"K_active mean {row['K_active_mean']:.1f}, K_off mean {row['K_off_mean']:.1f}",
                    label,
                    str(int(row["max_swaps"])),
                    f"{row['fraction_exact_mean']:.4f}",
                    f"{row['fraction_exact_p05']:.4f}",
                    f"{row['exact_rate']:.1%}",
                ]
            )
            + " |"
        )

    if not failure_cases.empty:
        lines.extend(
            [
                "",
                "## Worst Remaining Cyclic + 2-Swap Cases",
                "",
                "| N | K_active | K_off | off % | sample | start | fraction exact | overlap | swap distance | exact subset | local subset |",
                "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|",
            ]
        )
        for _, row in failure_cases.head(10).iterrows():
            lines.append(
                "| "
                + " | ".join(
                    [
                        str(int(row["N"])),
                        str(int(row["K_active"])),
                        str(int(row["K_off"])),
                        _format_pct(row["off_pct"]),
                        str(int(row["sample"])),
                        str(int(row["seed_position"])),
                        f"{row['fraction_exact_u_g']:.4f}",
                        f"{row['overlap_fraction']:.3f}",
                        str(int(row["swap_distance_to_exact"])),
                        str(row["exact_subset"]),
                        str(row["subset"]),
                    ]
                )
                + " |"
            )

    exact_completed = exact_runs["exact_completed"].astype(bool).mean() if not exact_runs.empty else 0.0
    loaded_rate = exact_runs["exact_source"].astype(str).eq("loaded").mean() if not exact_runs.empty else 0.0
    lines.extend(
        [
            "",
            "## Notes",
            "",
            f"- Exact enumeration completed for `{exact_completed:.1%}` of cases.",
            f"- Exact rows loaded from previous artifacts for `{loaded_rate:.1%}` of cases; missing real-off cases were recomputed.",
            "- Strong/weak under real-off semantics disables `K_off` antennas from the weakest and strongest row-power tails, then keeps the middle `K_active` antennas.",
            "- Historical reports that say `25% active` are not the real `25% off` task. `25% active` means `75% off`.",
            "",
            "## Plots",
            "",
            "![Raw U_G CDF](real_off_raw_u_g_cdf_by_off_pct.png)",
            "",
            "![Fraction exact CDF](real_off_fraction_exact_cdf.png)",
            "",
            "![Mean fraction by off percentage](real_off_mean_fraction_by_off_pct.png)",
            "",
            "![Exact recovery by off percentage](real_off_exact_recovery_by_off_pct.png)",
            "",
            "![Best cyclic start](real_off_best_cyclic_start_hist.png)",
            "",
            "![Failure diagnostics](real_off_failure_diagnostics.png)",
            "",
            "![Runtime](real_off_runtime_by_method.png)",
            "",
            "## Artifacts",
            "",
            "- Detailed CSVs are packed in `csv_data.tar.gz` after the run.",
            "- Main report: `local_threshold_real_off_report.md`.",
        ]
    )
    text = "\n".join(lines) + "\n"
    report_path = Path(out_dir) / "local_threshold_real_off_report.md"
    report_path.write_text(text, encoding="utf-8")
    if docs_path is not None:
        docs_path = Path(docs_path)
        docs_path.parent.mkdir(parents=True, exist_ok=True)
        docs_path.write_text(text, encoding="utf-8")


def write_active_k_local_report(
    exact_runs,
    local_runs,
    summary,
    diagnostics,
    failure_cases,
    out_dir,
    exact_source_dir=None,
    docs_path=None,
):
    global_summary = _global_summary(local_runs)
    by_active = _global_by_requested_active_pct(summary)
    lines = [
        "# Active-K Cyclic Threshold Local Search Study",
        "",
        "This exploratory report sweeps requested active antenna fractions directly.",
        "",
        "## K Semantics",
        "",
        "- `K_active = round(N * requested_active_pct / 100)` is the number of selected/kept antennas.",
        "- `K_off = N - K_active` is the number of disabled antennas.",
        "- The solver variable `K` equals `K_active`; it does not mean the number of antennas turned off.",
        "",
        "## Setup",
        "",
        f"- Exact source: `{exact_source_dir}`" if exact_source_dir is not None else "- Exact source: computed in this run when needed.",
        f"- Profiles: {', '.join(str(value) for value in sorted(local_runs['data_profile'].unique())) if not local_runs.empty else ''}",
        f"- N values: {', '.join(str(int(value)) for value in sorted(local_runs['N'].unique())) if not local_runs.empty else ''}",
        f"- L values: {', '.join(str(int(value)) for value in sorted(local_runs['L'].unique())) if not local_runs.empty else ''}",
        f"- Requested active K percentages: {', '.join(_format_pct(value) for value in sorted(local_runs['requested_active_pct'].unique())) if not local_runs.empty else ''}",
        "- Seeds compared: best tested cyclic window, `T=round(0.05N)`, and strong/weak.",
        "- Local search: greedy remove-one/add-one refinement by raw `U_G`, with 0, 1, or 2 swaps.",
        "",
        "## Local Swap Scheme And Cost",
        "",
        "- A `1-swap` removes one currently active antenna and adds one currently inactive antenna, preserving exact `K_active`.",
        "- One local-search pass evaluates every `(remove, add)` pair from the current active set and the configured add-candidate pool, applies only the single pair with the largest positive `U_G` improvement, and stops if no pair improves `U_G`.",
        "- `2-swap` means two greedy passes, not an exhaustive simultaneous two-pair exchange.",
        "- In this exact cyclic analysis the add-candidate pool is all inactive antennas, so one pass evaluates `K_active * (N - K_active)` swaps.",
        "- With `S` greedy passes, `A` add candidates, and `L` streams, local refinement costs `O(S * K_active * A * L^3)` time after row Gram matrices are built; here `A = N - K_active`, so it is `O(S * K_active * (N - K_active) * L^3)`.",
        "- Seed costs are separate: `T=0.05N` and strong/weak need one sorted window, while best cyclic scans `N` cyclic windows before local refinement.",
        "- Extra working space is `O(N * L^2 + K_active + A + L^2)` for row Gram matrices, active/add sets, and the current Gram matrix.",
        "",
        "## Direct Answer",
        "",
    ]
    if global_summary.empty:
        lines.append("- No rows were produced.")
    else:
        for seed_rule in REAL_OFF_SEED_RULES:
            chunk = global_summary[global_summary["seed_rule"] == seed_rule].sort_values("max_swaps")
            if chunk.empty:
                continue
            for _, row in chunk.iterrows():
                label = REAL_OFF_LABELS[f"{seed_rule}_s{int(row['max_swaps'])}"]
                lines.append(
                    f"- `{label}`: mean fraction exact `{row['fraction_exact_mean']:.4f}`, "
                    f"p05 `{row['fraction_exact_p05']:.4f}`, exact recovery `{row['exact_rate']:.1%}`, "
                    f"mean swaps applied `{row['swaps_applied_mean']:.3f}`."
                )

    best_t_summary = _best_cyclic_t_summary(local_runs)
    if not best_t_summary.empty:
        lines.extend(
            [
                "",
                "## Best Cyclic T Usually",
                "",
                "`T` is the cyclic start position in descending row-power order. `T=0` starts at the strongest row; larger `T` shifts the active window toward weaker rows.",
                "",
                "| requested active % | actual active % mean | T p05 | T median | T p95 | T/N median | T/K median |",
                "|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for _, row in best_t_summary.iterrows():
            lines.append(
                "| "
                + " | ".join(
                    [
                        _format_pct(row["requested_active_pct"]),
                        f"{row['active_pct_mean']:.2f}",
                        f"{row['T_p05']:.1f}",
                        f"{row['T_p50']:.1f}",
                        f"{row['T_p95']:.1f}",
                        f"{row['T_over_N_p50']:.3f}",
                        f"{row['T_over_K_p50']:.3f}",
                    ]
                )
                + " |"
            )

    lines.extend(
        [
            "",
            "## By Requested Active K Percentage",
            "",
            "| requested active % | actual active % mean | K_active mean | K_off mean | seed | swaps | mean fraction exact | p05 | exact rate |",
            "|---:|---:|---:|---:|---|---:|---:|---:|---:|",
        ]
    )
    for _, row in by_active.sort_values(["requested_active_pct", "seed_rule", "max_swaps"]).iterrows():
        label = REAL_OFF_LABELS.get(
            f"{row['seed_rule']}_s{int(row['max_swaps'])}",
            str(row["seed_rule"]),
        )
        lines.append(
            "| "
            + " | ".join(
                [
                    _format_pct(row["requested_active_pct"]),
                    f"{row['active_pct_mean']:.2f}",
                    f"{row['K_active_mean']:.1f}",
                    f"{row['K_off_mean']:.1f}",
                    label,
                    str(int(row["max_swaps"])),
                    f"{row['fraction_exact_mean']:.4f}",
                    f"{row['fraction_exact_p05']:.4f}",
                    f"{row['exact_rate']:.1%}",
                ]
            )
            + " |"
        )

    exact_completed = exact_runs["exact_completed"].astype(bool).mean() if not exact_runs.empty else 0.0
    loaded_rate = exact_runs["exact_source"].astype(str).eq("loaded").mean() if not exact_runs.empty else 0.0
    lines.extend(
        [
            "",
            "## Notes",
            "",
            f"- Exact enumeration completed for `{exact_completed:.1%}` of cases.",
            f"- Exact rows loaded from previous artifacts for `{loaded_rate:.1%}` of cases; missing active-K cases were recomputed.",
            "- Strong/weak is one non-cyclic sorted row-power window: it disables weakest and strongest tails and keeps the middle rows.",
            "- The cyclic best seed includes the strong/weak window as one possible cyclic or non-cyclic window, but chooses by tested `U_G`.",
            "",
            "## Plots",
            "",
            "![Raw U_G CDF](active_k_raw_u_g_cdf_by_requested_active_pct.png)",
            "",
            "![Fraction exact CDF](active_k_fraction_exact_cdf.png)",
            "",
            "![Mean fraction by requested active percentage](active_k_mean_fraction_by_requested_active_pct.png)",
            "",
            "![Exact recovery by requested active percentage](active_k_exact_recovery_by_requested_active_pct.png)",
            "",
            "![Best cyclic T boxplot](active_k_best_cyclic_T_boxplot.png)",
            "",
            "![Best cyclic T over N](active_k_best_cyclic_T_over_N.png)",
            "",
            "![Failure diagnostics](real_off_failure_diagnostics.png)",
            "",
            "![Runtime](active_k_runtime_by_method.png)",
            "",
            "## Artifacts",
            "",
            "- Detailed CSVs are packed in `csv_data.tar.gz` after the run.",
            "- Main report: `local_threshold_active_k_report.md`.",
        ]
    )
    text = "\n".join(lines) + "\n"
    report_path = Path(out_dir) / "local_threshold_active_k_report.md"
    report_path.write_text(text, encoding="utf-8")
    if docs_path is not None:
        docs_path = Path(docs_path)
        docs_path.parent.mkdir(parents=True, exist_ok=True)
        docs_path.write_text(text, encoding="utf-8")


def read_csv_artifact(base_dir, name):
    base_dir = Path(base_dir)
    path = base_dir / name
    if path.exists():
        return pd.read_csv(path)

    archive_path = base_dir / "csv_data.tar.gz"
    if not archive_path.exists():
        raise FileNotFoundError(f"Missing {name} and {archive_path}")
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive.getmembers():
            if Path(member.name).name == name:
                handle = archive.extractfile(member)
                if handle is None:
                    break
                return pd.read_csv(handle)
    raise FileNotFoundError(f"{name} not found in {archive_path}")


def _load_exact_source(exact_source_dir):
    if exact_source_dir is None:
        return {}
    exact_source_dir = Path(exact_source_dir)
    if not exact_source_dir.exists():
        return {}
    try:
        exact_runs = read_csv_artifact(exact_source_dir, "exact_runs.csv")
    except FileNotFoundError:
        return {}

    source = {}
    for _, row in exact_runs.iterrows():
        if not _bool(row.get("exact_completed", False)):
            continue
        key = _exact_key(row)
        source[key] = row.to_dict()
    return source


def _exact_for_case(V, case_base, exact_source, time_limit_seconds=120.0):
    key = _exact_key(case_base)
    if key in exact_source:
        row = exact_source[key]
        return {
            "exact_source": "loaded",
            "exact_completed": True,
            "exact_timed_out": _bool(row.get("exact_timed_out", False)),
            "exact_candidate_count": int(row["exact_candidate_count"]),
            "exact_evaluated_count": int(row["exact_evaluated_count"]),
            "exact_elapsed_seconds": float(row["exact_elapsed_seconds"]),
            "exact_subset": str(row["exact_subset"]),
            "exact_window_T": row.get("exact_window_T", np.nan),
            "exact_is_threshold_window": _bool(row.get("exact_is_threshold_window", False)),
            "exact_u_bf": float(row["exact_u_bf"]),
            "exact_u_i": float(row["exact_u_i"]),
            "exact_u_g": float(row["exact_u_g"]),
            "exact_u_g_db": float(row["exact_u_g_db"]),
        }

    exact = brute_force_exact_u_g(
        V,
        int(case_base["K_active"]),
        sigma=float(case_base["sigma"]),
        P=float(case_base["P"]),
        time_limit_seconds=time_limit_seconds,
    )
    exact_window_T = (
        contiguous_threshold_window_T(V, exact["subset"])
        if bool(exact["completed"])
        else None
    )
    exact_u_g = float(exact["u_g"])
    return {
        "exact_source": "computed",
        "exact_completed": bool(exact["completed"]),
        "exact_timed_out": bool(exact["timed_out"]),
        "exact_candidate_count": int(exact["candidate_count"]),
        "exact_evaluated_count": int(exact["evaluated_count"]),
        "exact_elapsed_seconds": float(exact["elapsed_seconds"]),
        "exact_subset": subset_to_string(exact["subset"]),
        "exact_window_T": np.nan if exact_window_T is None else int(exact_window_T),
        "exact_is_threshold_window": bool(exact_window_T is not None and exact["completed"]),
        "exact_u_bf": float(exact["u_bf"]),
        "exact_u_i": float(exact["u_i"]),
        "exact_u_g": exact_u_g,
        "exact_u_g_db": 10.0 * np.log10(max(exact_u_g, np.finfo(float).tiny)),
    }


def _local_rows_for_case(V, case_base, exact, max_swaps_values):
    seed_specs = _seed_specs_for_case(V, case_base)
    rows = []
    exact_u_g = float(exact["exact_u_g"])
    exact_completed = _bool(exact["exact_completed"])
    for seed in seed_specs:
        for max_swaps in max_swaps_values:
            result = refine_selection_by_swaps(
                V,
                seed["x"],
                max_swaps=int(max_swaps),
                sigma=float(case_base["sigma"]),
                P=float(case_base["P"]),
                seed_position=int(seed["seed_position"]),
            )
            fraction = (
                float(result["u_g"]) / exact_u_g
                if exact_completed and exact_u_g > 0.0
                else np.nan
            )
            exact_subset = parse_subset_string(exact["exact_subset"])
            local_subset = result["subset"]
            exact_features = _subset_features(V, exact_subset)
            local_features = _subset_features(V, local_subset)
            rows.append(
                {
                    **case_base,
                    "exact_u_g": exact_u_g if exact_completed else np.nan,
                    "exact_u_bf": float(exact["exact_u_bf"]),
                    "exact_u_i": float(exact["exact_u_i"]),
                    "exact_subset": str(exact["exact_subset"]),
                    "seed_rule": seed["seed_rule"],
                    "seed_label": seed["seed_label"],
                    "seed_position": int(seed["seed_position"]),
                    "seed_candidate_count": int(seed["candidate_count"]),
                    "max_swaps": int(max_swaps),
                    "active_count": int(np.sum(result["x"])),
                    "candidate_kind": seed["candidate_kind"]
                    if int(max_swaps) == 0
                    else f"{seed['candidate_kind']}_local_swap",
                    "candidate_count": int(result["candidate_count"]),
                    "candidate_radius": int(result["candidate_radius"]),
                    "add_candidate_count": int(result["add_candidate_count"]),
                    "evaluated_swap_count": int(result["evaluated_swap_count"]),
                    "swaps_applied": int(result["swaps_applied"]),
                    "initial_subset": subset_to_string(result["initial_subset"]),
                    "subset": subset_to_string(local_subset),
                    "swap_history": result["swap_history"],
                    "initial_u_bf": float(result["initial_u_bf"]),
                    "initial_u_i": float(result["initial_u_i"]),
                    "initial_u_g": float(result["initial_u_g"]),
                    "u_bf": float(result["u_bf"]),
                    "u_i": float(result["u_i"]),
                    "u_g": float(result["u_g"]),
                    "u_g_db": 10.0
                    * np.log10(max(float(result["u_g"]), np.finfo(float).tiny)),
                    "fraction_exact_u_g": fraction,
                    "gap_to_exact_pct": 100.0 * (1.0 - fraction)
                    if np.isfinite(fraction)
                    else np.nan,
                    "is_exact_best": bool(fraction >= 1.0 - 1e-9)
                    if np.isfinite(fraction)
                    else False,
                    "near_99": bool(fraction >= 0.99)
                    if np.isfinite(fraction)
                    else False,
                    "elapsed_seconds": float(result["elapsed_seconds"]),
                    **_feature_columns("exact", exact_features),
                    **_feature_columns("local", local_features),
                }
            )
    return rows


def _seed_specs_for_case(V, case_base):
    N = int(case_base["N"])
    K = int(case_base["K_active"])
    sigma = float(case_base["sigma"])
    P = float(case_base["P"])

    cyclic = best_cyclic_threshold_window(V, K, sigma=sigma, P=P)
    T_005 = int(np.clip(round(0.05 * N), 0, max(0, N - K)))
    formula_x = threshold_window_selection(V, K, T_005)
    strong_x = solve_h3_strong_weak(V, K, sigma=sigma, P=P)
    off_count = N - K
    strong_position = off_count - off_count // 2
    return [
        {
            "seed_rule": "best_cyclic_window",
            "seed_label": REAL_OFF_LABELS["best_cyclic_window_s0"],
            "seed_position": int(cyclic["T"]),
            "candidate_count": int(cyclic["candidate_count"]),
            "candidate_kind": "cyclic_threshold_window",
            "x": cyclic["x"],
        },
        {
            "seed_rule": "T_0p05N",
            "seed_label": REAL_OFF_LABELS["T_0p05N_s0"],
            "seed_position": int(T_005),
            "candidate_count": 1,
            "candidate_kind": "threshold_window_T_0p05N",
            "x": formula_x,
        },
        {
            "seed_rule": "strong_weak",
            "seed_label": REAL_OFF_LABELS["strong_weak_s0"],
            "seed_position": int(strong_position),
            "candidate_count": 1,
            "candidate_kind": "strong_weak",
            "x": strong_x,
        },
    ]


def _exact_key(row):
    return (
        str(row["data_profile"]),
        int(row["generator_seed"]),
        int(row["sample"]),
        int(row["N"]),
        int(row["L"]),
        int(row["K"]),
        round(float(row["sigma"]), 12),
        round(float(row["P"]), 12),
    )


def _case_columns(row):
    columns = {
        "data_profile": str(row["data_profile"]),
        "generator_seed": int(row["generator_seed"]),
        "sample": int(row["sample"]),
        "N": int(row["N"]),
        "L": int(row["L"]),
        "K": int(row["K"]),
        "K_active": int(row["K_active"]),
        "K_off": int(row["K_off"]),
        "off_pct": float(row["off_pct"]),
        "active_pct": float(row["active_pct"]),
        "sigma": float(row["sigma"]),
        "P": float(row["P"]),
    }
    if "requested_active_pct" in row:
        columns["requested_active_pct"] = float(row["requested_active_pct"])
    return columns


def _subset_features(V, subset):
    subset = tuple(int(value) for value in subset)
    if not subset:
        return {
            "max_row_power": 0.0,
            "scaled_trace": 0.0,
            "eig_min": 0.0,
            "eig_max": 0.0,
            "eig_balance": 0.0,
        }
    row_powers = np.sum(np.abs(V) ** 2, axis=1).real
    gram = V[list(subset)].conj().T @ V[list(subset)]
    max_row_power = float(np.max(row_powers[list(subset)]))
    z2 = 1.0 / max_row_power if max_row_power > 0 else 0.0
    gram_sq = gram @ gram.conj().T
    scaled = z2 * gram_sq
    eigvals = np.linalg.eigvalsh(scaled)
    eig_min = float(np.min(np.real(eigvals)))
    eig_max = float(np.max(np.real(eigvals)))
    return {
        "max_row_power": max_row_power,
        "scaled_trace": float(np.real(np.trace(scaled))),
        "eig_min": eig_min,
        "eig_max": eig_max,
        "eig_balance": _safe_ratio(eig_min, eig_max),
    }


def _feature_columns(prefix, features):
    return {f"{prefix}_{key}": float(value) for key, value in features.items()}


def _features_from_serialized(row, prefix):
    return {
        "max_row_power": float(row[f"{prefix}_max_row_power"]),
        "scaled_trace": float(row[f"{prefix}_scaled_trace"]),
        "eig_min": float(row[f"{prefix}_eig_min"]),
        "eig_max": float(row[f"{prefix}_eig_max"]),
        "eig_balance": float(row[f"{prefix}_eig_balance"]),
    }


def _global_summary(local_runs):
    if local_runs.empty:
        return pd.DataFrame()

    def q(value):
        return lambda data: data.quantile(value)

    return (
        local_runs.groupby(["seed_rule", "max_swaps"], as_index=False)
        .agg(
            cases=("u_g", "count"),
            fraction_exact_mean=("fraction_exact_u_g", "mean"),
            fraction_exact_p05=("fraction_exact_u_g", q(0.05)),
            fraction_exact_p50=("fraction_exact_u_g", q(0.50)),
            fraction_exact_p95=("fraction_exact_u_g", q(0.95)),
            exact_rate=("is_exact_best", "mean"),
            near_99_rate=("near_99", "mean"),
            swaps_applied_mean=("swaps_applied", "mean"),
            elapsed_seconds_mean=("elapsed_seconds", "mean"),
        )
        .sort_values(["seed_rule", "max_swaps"])
    )


def _global_by_off_pct(summary):
    if summary.empty:
        return pd.DataFrame()
    weighted = []
    group_cols = ["off_pct", "seed_rule", "max_swaps"]
    for keys, group in summary.groupby(group_cols):
        cases = group["cases"].astype(float)
        total = float(cases.sum())
        if total <= 0:
            continue
        row = dict(zip(group_cols, keys))
        for col in (
            "fraction_exact_mean",
            "exact_rate",
            "near_99_rate",
            "swaps_applied_mean",
            "elapsed_seconds_mean",
        ):
            row[col] = float(np.average(group[col].astype(float), weights=cases))
        row["cases"] = int(total)
        row["fraction_exact_p05"] = float(group["fraction_exact_p05"].min())
        row["K_active_mean"] = float(np.average(group["K_active"].astype(float), weights=cases))
        row["K_off_mean"] = float(np.average(group["K_off"].astype(float), weights=cases))
        weighted.append(row)
    return pd.DataFrame(weighted)


def _global_by_requested_active_pct(summary):
    if summary.empty:
        return pd.DataFrame()
    if "requested_active_pct" not in summary.columns:
        return pd.DataFrame()
    weighted = []
    group_cols = ["requested_active_pct", "seed_rule", "max_swaps"]
    for keys, group in summary.groupby(group_cols):
        cases = group["cases"].astype(float)
        total = float(cases.sum())
        if total <= 0:
            continue
        row = dict(zip(group_cols, keys))
        for col in (
            "fraction_exact_mean",
            "exact_rate",
            "near_99_rate",
            "swaps_applied_mean",
            "elapsed_seconds_mean",
        ):
            row[col] = float(np.average(group[col].astype(float), weights=cases))
        row["cases"] = int(total)
        row["fraction_exact_p05"] = float(group["fraction_exact_p05"].min())
        row["active_pct_mean"] = float(np.average(group["active_pct"].astype(float), weights=cases))
        row["K_active_mean"] = float(np.average(group["K_active"].astype(float), weights=cases))
        row["K_off_mean"] = float(np.average(group["K_off"].astype(float), weights=cases))
        weighted.append(row)
    return pd.DataFrame(weighted)


def _plot_real_off_failure_diagnostics(diagnostics, out_dir):
    from utils.plotting import use_agg_backend

    use_agg_backend()
    import matplotlib.pyplot as plt

    primary = diagnostics[diagnostics["seed_rule"] == "best_cyclic_window"]
    if primary.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4), squeeze=False)
    colors = {0: "#0072B2", 1: "#009E73", 2: "#56B4E9"}
    for max_swaps, chunk in primary.groupby("max_swaps"):
        axes[0, 0].scatter(
            chunk["overlap_fraction"],
            chunk["fraction_exact_u_g"],
            s=16,
            alpha=0.55,
            color=colors.get(int(max_swaps)),
            label=f"{int(max_swaps)} swaps",
        )
        axes[0, 1].scatter(
            chunk["local_exact_eig_balance_delta"],
            chunk["fraction_exact_u_g"],
            s=16,
            alpha=0.55,
            color=colors.get(int(max_swaps)),
            label=f"{int(max_swaps)} swaps",
        )
    axes[0, 0].set_xlabel("overlap with exact subset")
    axes[0, 0].set_ylabel("U_G / exact U_G")
    axes[0, 0].grid(True, alpha=0.25)
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].set_xlabel("local minus exact eigen-balance")
    axes[0, 1].set_ylabel("U_G / exact U_G")
    axes[0, 1].grid(True, alpha=0.25)
    fig.suptitle("Cyclic local failure diagnostics")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(Path(out_dir) / "real_off_failure_diagnostics.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_real_off_runtime(summary, out_dir):
    from utils.plotting import use_agg_backend

    use_agg_backend()
    import matplotlib.pyplot as plt

    if summary.empty:
        return
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for seed_rule in REAL_OFF_SEED_RULES:
        for max_swaps in REAL_OFF_MAX_SWAPS:
            chunk = summary[
                (summary["seed_rule"] == seed_rule)
                & (summary["max_swaps"] == max_swaps)
            ]
            if chunk.empty:
                continue
            grouped = chunk.groupby("N", as_index=False)["elapsed_seconds_mean"].mean()
            key = f"{seed_rule}_s{max_swaps}"
            ax.plot(
                grouped["N"],
                grouped["elapsed_seconds_mean"],
                marker="o",
                color=REAL_OFF_COLORS.get(key),
                label=REAL_OFF_LABELS.get(key, key),
            )
    ax.set_xlabel("N")
    ax.set_ylabel("mean runtime seconds")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7)
    ax.set_title("Real-off local search runtime")
    fig.tight_layout()
    fig.savefig(Path(out_dir) / "real_off_runtime_by_method.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_active_k_best_cyclic_t(local_runs, out_dir):
    from utils.plotting import use_agg_backend

    use_agg_backend()
    import matplotlib.pyplot as plt

    best = local_runs[
        (local_runs["seed_rule"] == "best_cyclic_window")
        & (local_runs["max_swaps"] == 0)
    ].copy()
    if best.empty:
        return
    best["T_over_N"] = best["seed_position"].astype(float) / best["N"].astype(float)
    active_pcts = sorted(best["requested_active_pct"].unique())

    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    data = [
        best[best["requested_active_pct"] == active_pct]["seed_position"].astype(float)
        for active_pct in active_pcts
    ]
    ax.boxplot(data, tick_labels=[_format_pct(value) for value in active_pcts], showmeans=True)
    ax.set_xlabel("requested active K percentage")
    ax.set_ylabel("best cyclic T")
    ax.grid(True, axis="y", alpha=0.25)
    ax.set_title("Best cyclic T distribution by requested active percentage")
    fig.tight_layout()
    fig.savefig(Path(out_dir) / "active_k_best_cyclic_T_boxplot.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    summary = _best_cyclic_t_summary(local_runs)
    if summary.empty:
        return
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    x_values = summary["requested_active_pct"].astype(float).to_numpy()
    ax.plot(
        x_values,
        summary["T_over_N_p50"].astype(float),
        marker="o",
        color=REAL_OFF_COLORS["best_cyclic_window_s0"],
        label="median best T/N",
    )
    ax.fill_between(
        x_values,
        summary["T_over_N_p05"].astype(float),
        summary["T_over_N_p95"].astype(float),
        color=REAL_OFF_COLORS["best_cyclic_window_s0"],
        alpha=0.18,
        label="p05..p95",
    )
    ax.set_xlabel("requested active K percentage")
    ax.set_ylabel("best cyclic T / N")
    ax.set_ylim(0.0, 1.0)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    ax.set_title("Best cyclic T/N by requested active percentage")
    fig.tight_layout()
    fig.savefig(Path(out_dir) / "active_k_best_cyclic_T_over_N.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_active_k_runtime(summary, out_dir):
    from utils.plotting import use_agg_backend

    use_agg_backend()
    import matplotlib.pyplot as plt

    if summary.empty:
        return
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for seed_rule in REAL_OFF_SEED_RULES:
        for max_swaps in REAL_OFF_MAX_SWAPS:
            chunk = summary[
                (summary["seed_rule"] == seed_rule)
                & (summary["max_swaps"] == max_swaps)
            ]
            if chunk.empty:
                continue
            grouped = chunk.groupby("N", as_index=False)["elapsed_seconds_mean"].mean()
            key = f"{seed_rule}_s{max_swaps}"
            ax.plot(
                grouped["N"],
                grouped["elapsed_seconds_mean"],
                marker="o",
                color=REAL_OFF_COLORS.get(key),
                label=REAL_OFF_LABELS.get(key, key),
            )
    ax.set_xlabel("N")
    ax.set_ylabel("mean local-refinement runtime seconds")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7)
    ax.set_title("Active-K local refinement runtime")
    fig.tight_layout()
    fig.savefig(Path(out_dir) / "active_k_runtime_by_method.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def _best_cyclic_t_summary(local_runs):
    if local_runs.empty or "requested_active_pct" not in local_runs.columns:
        return pd.DataFrame()
    best = local_runs[
        (local_runs["seed_rule"] == "best_cyclic_window")
        & (local_runs["max_swaps"] == 0)
    ].copy()
    if best.empty:
        return pd.DataFrame()
    best["T"] = best["seed_position"].astype(float)
    best["T_over_N"] = best["T"] / best["N"].astype(float)
    best["T_over_K"] = best["T"] / best["K_active"].replace(0, np.nan).astype(float)

    def q(value):
        return lambda data: data.quantile(value)

    return (
        best.groupby("requested_active_pct", as_index=False)
        .agg(
            cases=("T", "count"),
            active_pct_mean=("active_pct", "mean"),
            K_active_mean=("K_active", "mean"),
            T_p05=("T", q(0.05)),
            T_p50=("T", q(0.50)),
            T_p95=("T", q(0.95)),
            T_over_N_p05=("T_over_N", q(0.05)),
            T_over_N_p50=("T_over_N", q(0.50)),
            T_over_N_p95=("T_over_N", q(0.95)),
            T_over_K_p50=("T_over_K", q(0.50)),
        )
        .sort_values("requested_active_pct")
    )


def _rule_specs():
    return [
        (f"{seed_rule}_s{max_swaps}", seed_rule, max_swaps)
        for seed_rule in REAL_OFF_SEED_RULES
        for max_swaps in REAL_OFF_MAX_SWAPS
    ]


def _empirical_cdf(values):
    values = pd.Series(values, dtype=float).replace([np.inf, -np.inf], np.nan).dropna()
    values = np.sort(values.to_numpy())
    if len(values) == 0:
        return values, values
    cumulative = np.arange(1, len(values) + 1, dtype=float) / float(len(values))
    return values, cumulative


def _safe_ratio(numerator, denominator):
    numerator = float(numerator)
    denominator = float(denominator)
    if abs(denominator) <= np.finfo(float).eps:
        return np.nan
    return numerator / denominator


def _bool(value):
    if pd.isna(value):
        return False
    if isinstance(value, str):
        return value.strip().lower() in ("true", "1", "yes")
    return bool(value)


def _format_pct(value):
    value = float(value)
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:g}"


def _format_float(value):
    if pd.isna(value):
        return ""
    value = float(value)
    if abs(value) >= 1e6 or (0 < abs(value) < 1e-3):
        return f"{value:.3e}"
    return f"{value:.4f}"
