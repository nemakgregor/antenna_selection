from pathlib import Path

import numpy as np
import pandas as pd

from algorithms.h3_threshold_local import evaluate_threshold_local_rules
from utils.brute_force import parse_subset_string
from utils.data import generate_v_profile_from_rng
from utils.io import atomic_write_csv


LOCAL_SEED_RULES = ("best_tested_T", "T_0", "T_0p025N", "T_0p05N")
LOCAL_MAX_SWAPS = (0, 1, 2)

LOCAL_LABELS = {
    "exact_best": "exact best",
    "best_tested_T_s0": "best tested T",
    "best_tested_T_s1": "best tested T + 1 swap",
    "best_tested_T_s2": "best tested T + 2 swaps",
    "T_0_s0": "T=0",
    "T_0p025N_s0": "T=0.025N",
    "T_0p05N_s0": "T=0.05N",
    "T_0p05N_s1": "T=0.05N + 1 swap",
    "T_0p05N_s2": "T=0.05N + 2 swaps",
}


def run_local_threshold_exact_analysis(
    exact_dir,
    out_dir,
    docs_path=None,
    n_values=(8, 12, 16, 20),
    k_pcts=(25, 30, 35, 40, 45, 50),
    profiles=("gaussian",),
    seed_rules=LOCAL_SEED_RULES,
    max_swaps_values=LOCAL_MAX_SWAPS,
    candidate_radius=None,
):
    exact_dir = Path(exact_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    exact_runs = _load_filtered_exact_runs(exact_dir, n_values, k_pcts, profiles)
    matrix_cache = _reconstruct_matrices(exact_runs)
    local_runs = _collect_local_runs(
        exact_runs,
        matrix_cache,
        seed_rules=seed_rules,
        max_swaps_values=max_swaps_values,
        candidate_radius=candidate_radius,
    )
    summary = build_local_threshold_summary(local_runs)
    diagnostics = build_local_threshold_diagnostics(local_runs, matrix_cache)
    failure_cases = build_local_threshold_failure_cases(local_runs, diagnostics)

    atomic_write_csv(local_runs, out_dir / "local_threshold_runs.csv")
    atomic_write_csv(summary, out_dir / "local_threshold_summary.csv")
    atomic_write_csv(failure_cases, out_dir / "local_threshold_failure_cases.csv")
    atomic_write_csv(diagnostics, out_dir / "local_threshold_diagnostics.csv")

    write_local_threshold_plots(exact_runs, local_runs, summary, diagnostics, out_dir)
    write_local_threshold_report(
        exact_runs,
        local_runs,
        summary,
        diagnostics,
        failure_cases,
        out_dir,
        exact_dir,
        docs_path=docs_path,
    )

    return {
        "exact_runs": exact_runs,
        "local_runs": local_runs,
        "summary": summary,
        "diagnostics": diagnostics,
        "failure_cases": failure_cases,
    }


def build_local_threshold_summary(local_runs):
    if local_runs.empty:
        return pd.DataFrame()

    def q(value):
        return lambda data: data.quantile(value)

    return (
        local_runs.groupby(
            [
                "data_profile",
                "N",
                "L",
                "K",
                "active_pct",
                "seed_rule",
                "max_swaps",
            ],
            as_index=False,
        )
        .agg(
            cases=("u_g", "count"),
            seed_T_mean=("seed_T", "mean"),
            seed_T_p50=("seed_T", q(0.50)),
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
        .sort_values(["data_profile", "N", "K", "active_pct", "seed_rule", "max_swaps"])
    )


def build_local_threshold_diagnostics(local_runs, matrix_cache):
    rows = []
    for _, row in local_runs.iterrows():
        key = _matrix_key(row)
        V = matrix_cache[key]
        row_powers = np.sum(np.abs(V) ** 2, axis=1).real
        order = np.argsort(row_powers)[::-1]
        ranks = np.empty(len(order), dtype=int)
        ranks[order] = np.arange(len(order))

        exact_subset = parse_subset_string(row["exact_subset"])
        local_subset = parse_subset_string(row["subset"])
        exact_set = set(exact_subset)
        local_set = set(local_subset)
        overlap = len(exact_set & local_set)
        K = int(row["K"])
        exact_ranks = sorted(int(ranks[index]) for index in exact_subset)
        local_ranks = sorted(int(ranks[index]) for index in local_subset)
        if len(exact_ranks) == len(local_ranks) and exact_ranks:
            rank_abs_gap = np.abs(np.asarray(exact_ranks) - np.asarray(local_ranks))
            rank_gap_mean = float(np.mean(rank_abs_gap))
            rank_gap_max = float(np.max(rank_abs_gap))
        else:
            rank_gap_mean = np.nan
            rank_gap_max = np.nan

        exact_features = _subset_features(V, exact_subset)
        local_features = _subset_features(V, local_subset)
        rows.append(
            {
                **_case_columns(row),
                "seed_rule": row["seed_rule"],
                "max_swaps": int(row["max_swaps"]),
                "seed_T": int(row["seed_T"]),
                "fraction_exact_u_g": float(row["fraction_exact_u_g"]),
                "gap_to_exact_pct": float(row["gap_to_exact_pct"]),
                "is_exact_best": bool(row["is_exact_best"]),
                "overlap_count": int(overlap),
                "overlap_fraction": overlap / max(K, 1),
                "swap_distance_to_exact": int(K - overlap),
                "rank_gap_mean": rank_gap_mean,
                "rank_gap_max": rank_gap_max,
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


def build_local_threshold_failure_cases(local_runs, diagnostics):
    if local_runs.empty or diagnostics.empty:
        return pd.DataFrame()
    primary = diagnostics[
        (diagnostics["seed_rule"] == "best_tested_T")
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
        "active_pct",
        "seed_rule",
        "max_swaps",
        "seed_T",
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
            "active_pct",
            "seed_rule",
            "max_swaps",
            "seed_T",
            "fraction_exact_u_g",
        ],
        how="left",
    )
    return merged.sort_values("fraction_exact_u_g").head(100)


def write_local_threshold_plots(exact_runs, local_runs, summary, diagnostics, out_dir):
    if exact_runs.empty or local_runs.empty:
        return

    from utils.plotting import use_agg_backend

    use_agg_backend()
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    active_pcts = sorted(exact_runs["active_pct"].unique())
    raw_rules = [
        ("exact_best", None, None),
        ("best_tested_T_s0", "best_tested_T", 0),
        ("best_tested_T_s1", "best_tested_T", 1),
        ("best_tested_T_s2", "best_tested_T", 2),
        ("T_0p05N_s0", "T_0p05N", 0),
    ]
    colors = {
        "exact_best": "#000000",
        "best_tested_T_s0": "#0072B2",
        "best_tested_T_s1": "#009E73",
        "best_tested_T_s2": "#56B4E9",
        "T_0p05N_s0": "#E69F00",
        "T_0p05N_s1": "#CC79A7",
        "T_0p05N_s2": "#D55E00",
        "T_0_s0": "#999999",
        "T_0p025N_s0": "#8B5A2B",
    }

    cols = min(3, max(1, len(active_pcts)))
    rows = int(np.ceil(len(active_pcts) / cols))
    fig, axes = plt.subplots(
        rows,
        cols,
        figsize=(5.2 * cols, 4.0 * rows),
        squeeze=False,
    )
    for ax, active_pct in zip(axes.ravel(), active_pcts):
        exact_chunk = exact_runs[exact_runs["active_pct"] == active_pct]
        values, cumulative = _empirical_cdf(exact_chunk["exact_u_g"])
        if len(values):
            ax.step(
                cumulative,
                values,
                where="post",
                color=colors["exact_best"],
                label=LOCAL_LABELS["exact_best"],
            )
        local_pct = local_runs[local_runs["active_pct"] == active_pct]
        for rule_key, seed_rule, max_swaps in raw_rules[1:]:
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
                    color=colors.get(rule_key),
                    label=LOCAL_LABELS[rule_key],
                )
        actual_pct = (
            100.0
            * exact_chunk["K"].astype(float)
            / exact_chunk["N"].astype(float)
        ).mean()
        ax.set_title(f"requested K={active_pct:g}% (actual mean {actual_pct:.2f}%)")
        ax.set_xlabel("cumulative fraction")
        ax.set_ylabel("raw U_G")
        ax.set_yscale("log")
        ax.grid(True, alpha=0.25)
    for ax in axes.ravel()[len(active_pcts) :]:
        ax.axis("off")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Local threshold raw U_G CDF by active K percentage")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_dir / "local_threshold_raw_u_g_cdf_by_k_pct.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fraction_rules = [
        ("best_tested_T_s0", "best_tested_T", 0),
        ("best_tested_T_s1", "best_tested_T", 1),
        ("best_tested_T_s2", "best_tested_T", 2),
        ("T_0p05N_s0", "T_0p05N", 0),
        ("T_0p05N_s1", "T_0p05N", 1),
        ("T_0p05N_s2", "T_0p05N", 2),
    ]
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    for rule_key, seed_rule, max_swaps in fraction_rules:
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
                color=colors.get(rule_key),
                label=LOCAL_LABELS[rule_key],
            )
    ax.set_xlabel("cumulative fraction")
    ax.set_ylabel("U_G(local) / exact U_G")
    ax.set_ylim(0.0, 1.03)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    ax.set_title("Local threshold fraction of exact U_G")
    fig.tight_layout()
    fig.savefig(out_dir / "local_threshold_fraction_exact_cdf.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    selected_summary = _global_by_active_pct(summary)
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    for rule_key, seed_rule, max_swaps in fraction_rules:
        chunk = selected_summary[
            (selected_summary["seed_rule"] == seed_rule)
            & (selected_summary["max_swaps"] == max_swaps)
        ].sort_values("active_pct")
        if chunk.empty:
            continue
        ax.plot(
            chunk["active_pct"],
            chunk["fraction_exact_mean"],
            marker="o",
            color=colors.get(rule_key),
            label=LOCAL_LABELS[rule_key],
        )
    ax.set_xlabel("requested active K percentage")
    ax.set_ylabel("mean U_G / exact U_G")
    ax.set_ylim(0.0, 1.03)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    ax.set_title("Mean local-search quality by active K percentage")
    fig.tight_layout()
    fig.savefig(out_dir / "local_threshold_mean_fraction_by_active_pct.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    for rule_key, seed_rule, max_swaps in fraction_rules[:3]:
        chunk = selected_summary[
            (selected_summary["seed_rule"] == seed_rule)
            & (selected_summary["max_swaps"] == max_swaps)
        ].sort_values("active_pct")
        if chunk.empty:
            continue
        ax.plot(
            chunk["active_pct"],
            chunk["exact_rate"],
            marker="o",
            color=colors.get(rule_key),
            label=LOCAL_LABELS[rule_key],
        )
    ax.set_xlabel("requested active K percentage")
    ax.set_ylabel("exact recovery rate")
    ax.set_ylim(0.0, 1.03)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    ax.set_title("Best-tested-T local search exact recovery")
    fig.tight_layout()
    fig.savefig(out_dir / "local_threshold_exact_rate_by_active_pct.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    _plot_seed_dependence(summary, out_dir, colors)
    _plot_failure_diagnostics(diagnostics, out_dir)
    _plot_runtime(summary, out_dir, colors)


def write_local_threshold_report(
    exact_runs,
    local_runs,
    summary,
    diagnostics,
    failure_cases,
    out_dir,
    exact_dir,
    docs_path=None,
):
    lines = [
        "# Local Threshold Exact Gaussian Study",
        "",
        "This experiment starts from a pure row-power threshold window and then applies greedy one-swap or two-swap local search by `U_G`.",
        "It uses saved exact Gaussian cases and reconstructs matrices from `(profile, seed, sample)`; brute force is not rerun here.",
        "",
        "## Setup",
        "",
        f"- Exact source: `{exact_dir}`",
        f"- Profiles: {', '.join(str(value) for value in sorted(local_runs['data_profile'].unique())) if not local_runs.empty else ''}",
        f"- N values: {', '.join(str(int(value)) for value in sorted(local_runs['N'].unique())) if not local_runs.empty else ''}",
        f"- L values: {', '.join(str(int(value)) for value in sorted(local_runs['L'].unique())) if not local_runs.empty else ''}",
        f"- Requested active K percentages: {', '.join(_format_pct(value) for value in sorted(local_runs['active_pct'].unique())) if not local_runs.empty else ''}",
        "- Seed rules: `best_tested_T`, `T=0`, `T=0.025N`, `T=0.05N`.",
        "- Local search: all active rows are removable; add candidates are inactive rows near the window boundary.",
        "",
        "## Direct Answer",
        "",
    ]

    global_summary = _global_summary(local_runs)
    primary = global_summary[global_summary["seed_rule"] == "best_tested_T"].sort_values(
        "max_swaps"
    )
    if not primary.empty:
        for _, row in primary.iterrows():
            lines.append(
                f"- `best_tested_T + {int(row['max_swaps'])} swaps`: "
                f"mean fraction exact `{row['fraction_exact_mean']:.4f}`, "
                f"p05 `{row['fraction_exact_p05']:.4f}`, exact recovery `{row['exact_rate']:.1%}`, "
                f"mean swaps applied `{row['swaps_applied_mean']:.3f}`."
            )
        zero = primary[primary["max_swaps"] == 0]
        two = primary[primary["max_swaps"] == 2]
        if not zero.empty and not two.empty:
            lines.append(
                f"- Two swaps improve the headline threshold-window approach by "
                f"`{100.0 * (two.iloc[0]['fraction_exact_mean'] - zero.iloc[0]['fraction_exact_mean']):.3f}` percentage points in mean exact fraction."
            )
    else:
        lines.append("- No local rows were produced.")

    t005 = global_summary[global_summary["seed_rule"] == "T_0p05N"].sort_values(
        "max_swaps"
    )
    if not t005.empty:
        best_t005 = t005.iloc[-1]
        lines.append(
            f"- `T=0.05N + 2 swaps` reaches mean fraction exact `{best_t005['fraction_exact_mean']:.4f}` with exact recovery `{best_t005['exact_rate']:.1%}`."
        )

    primary_diag = diagnostics[
        (diagnostics["seed_rule"] == "best_tested_T")
        & (diagnostics["max_swaps"] == 2)
    ]
    misses = primary_diag[primary_diag["fraction_exact_u_g"] < 1.0 - 1e-9]
    if not primary_diag.empty:
        lines.append(
            f"- After two swaps from `best_tested_T`, remaining misses have mean overlap `{misses['overlap_fraction'].mean():.3f}` with exact "
            f"and mean swap distance `{misses['swap_distance_to_exact'].mean():.3f}` rows."
            if not misses.empty
            else "- After two swaps from `best_tested_T`, every analyzed case recovered exact `U_G`."
        )

    lines.extend(
        [
            "",
            "## Global Rule Summary",
            "",
            "| seed rule | swaps allowed | cases | mean fraction exact | p05 | p50 | p95 | exact rate | near-99 rate | mean runtime s |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for _, row in global_summary.sort_values(
        ["seed_rule", "max_swaps"], ascending=[True, True]
    ).iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["seed_rule"]),
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

    active_summary = _global_by_active_pct(summary)
    lines.extend(
        [
            "",
            "## Best Tested T Local Search By K%",
            "",
            "| requested K% | swaps | mean fraction exact | p05 | exact rate | mean swaps applied |",
            "|---:|---:|---:|---:|---:|---:|",
        ]
    )
    best_active = active_summary[active_summary["seed_rule"] == "best_tested_T"]
    for _, row in best_active.sort_values(["active_pct", "max_swaps"]).iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    _format_pct(row["active_pct"]),
                    str(int(row["max_swaps"])),
                    f"{row['fraction_exact_mean']:.4f}",
                    f"{row['fraction_exact_p05']:.4f}",
                    f"{row['exact_rate']:.1%}",
                    f"{row['swaps_applied_mean']:.3f}",
                ]
            )
            + " |"
        )

    if not failure_cases.empty:
        lines.extend(
            [
                "",
                "## Worst Remaining Best-T + 2-Swap Cases",
                "",
                "| N | K | requested K% | sample | T | fraction exact | overlap | swap distance | exact subset | local subset |",
                "|---:|---:|---:|---:|---:|---:|---:|---:|---|---|",
            ]
        )
        display = failure_cases.head(10)
        for _, row in display.iterrows():
            lines.append(
                "| "
                + " | ".join(
                    [
                        str(int(row["N"])),
                        str(int(row["K"])),
                        _format_pct(row["active_pct"]),
                        str(int(row["sample"])),
                        str(int(row["seed_T"])),
                        f"{row['fraction_exact_u_g']:.4f}",
                        f"{row['overlap_fraction']:.3f}",
                        str(int(row["swap_distance_to_exact"])),
                        str(row["exact_subset"]),
                        str(row["subset"]),
                    ]
                )
                + " |"
            )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- If `best_tested_T + 1/2 swaps` closes most remaining exact gap, the threshold-window approach is structurally strong: sorting by row power finds almost the right neighborhood, and local swaps repair non-contiguous exact subsets.",
            "- Remaining misses mean the exact subset is outside the fixed boundary neighborhood, or needs a coordinated multi-row replacement that greedy one-at-a-time swaps cannot see.",
            "- `T=0.05N + local search` measures whether the practical formula can benefit from the same repair step, while `best_tested_T + local search` measures the approach itself independent of formula choice.",
            "",
            "## Plots",
            "",
            "![Raw U_G CDF by K percentage](local_threshold_raw_u_g_cdf_by_k_pct.png)",
            "",
            "![Fraction exact CDF](local_threshold_fraction_exact_cdf.png)",
            "",
            "![Mean fraction by active percentage](local_threshold_mean_fraction_by_active_pct.png)",
            "",
            "![Exact recovery rate](local_threshold_exact_rate_by_active_pct.png)",
            "",
            "![Seed threshold dependence](local_threshold_seed_dependence.png)",
            "",
            "![Failure diagnostics](local_threshold_failure_diagnostics.png)",
            "",
            "![Runtime](local_threshold_runtime_by_N_K.png)",
            "",
            "## Artifacts",
            "",
            "- `local_threshold_runs.csv`",
            "- `local_threshold_summary.csv`",
            "- `local_threshold_failure_cases.csv`",
            "- `local_threshold_diagnostics.csv`",
        ]
    )

    text = "\n".join(lines) + "\n"
    report_path = Path(out_dir) / "local_threshold_exact_gauss_report.md"
    report_path.write_text(text, encoding="utf-8")
    if docs_path is not None:
        docs_path = Path(docs_path)
        docs_path.parent.mkdir(parents=True, exist_ok=True)
        docs_path.write_text(text, encoding="utf-8")


def _load_filtered_exact_runs(exact_dir, n_values, k_pcts, profiles):
    exact_path = Path(exact_dir) / "exact_runs.csv"
    if not exact_path.exists():
        raise FileNotFoundError(f"No exact_runs.csv found at {exact_path}")
    exact_runs = pd.read_csv(exact_path)
    exact_runs = exact_runs[_boolean_series(exact_runs["exact_completed"])].copy()
    exact_runs = exact_runs[exact_runs["N"].astype(int).isin([int(value) for value in n_values])]
    exact_runs = exact_runs[
        exact_runs["data_profile"].astype(str).isin([str(value) for value in profiles])
    ]
    pct_values = {round(float(value), 8) for value in k_pcts}
    exact_runs = exact_runs[
        exact_runs["active_pct"].astype(float).round(8).isin(pct_values)
    ].copy()
    if exact_runs.empty:
        raise ValueError("No exact rows remain after local-analysis filters.")
    return exact_runs.sort_values(
        ["data_profile", "N", "L", "generator_seed", "sample", "active_pct", "K"]
    ).reset_index(drop=True)


def _reconstruct_matrices(exact_runs):
    cache = {}
    group_cols = ["data_profile", "generator_seed", "N", "L"]
    for (profile, generator_seed, N, L), group in exact_runs.groupby(group_cols):
        samples = set(int(value) for value in group["sample"].unique())
        max_sample = max(samples)
        rng = np.random.RandomState(int(generator_seed))
        for sample in range(max_sample + 1):
            V = generate_v_profile_from_rng(rng, int(N), int(L), str(profile))
            if sample in samples:
                cache[(str(profile), int(generator_seed), int(sample), int(N), int(L))] = V
    return cache


def _collect_local_runs(
    exact_runs,
    matrix_cache,
    seed_rules=LOCAL_SEED_RULES,
    max_swaps_values=LOCAL_MAX_SWAPS,
    candidate_radius=None,
):
    rows = []
    for _, exact_row in exact_runs.iterrows():
        key = _matrix_key(exact_row)
        V = matrix_cache[key]
        K = int(exact_row["K"])
        rules = _seed_rules_for_row(exact_row, seed_rules)
        local_rows = evaluate_threshold_local_rules(
            V,
            K,
            rules,
            max_swaps_values=max_swaps_values,
            sigma=float(exact_row["sigma"]),
            P=float(exact_row["P"]),
            candidate_radius=candidate_radius,
        )
        exact_u_g = float(exact_row["exact_u_g"])
        for local in local_rows:
            fraction = local["u_g"] / exact_u_g if exact_u_g > 0 else np.nan
            rows.append(
                {
                    **_case_columns(exact_row),
                    "exact_u_g": exact_u_g,
                    "exact_u_bf": float(exact_row["exact_u_bf"]),
                    "exact_u_i": float(exact_row["exact_u_i"]),
                    "exact_subset": str(exact_row["exact_subset"]),
                    "exact_window_T": exact_row["exact_window_T"],
                    "exact_is_threshold_window": bool(exact_row["exact_is_threshold_window"]),
                    **local,
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
                }
            )
    return pd.DataFrame(rows)


def _seed_rules_for_row(row, seed_rules):
    N = int(row["N"])
    K = int(row["K"])
    values = []
    for seed_rule in seed_rules:
        if seed_rule == "best_tested_T":
            T = int(row["best_tested_T"])
        elif seed_rule == "T_0":
            T = 0
        elif seed_rule == "T_0p025N":
            T = int(round(0.025 * N))
        elif seed_rule == "T_0p05N":
            T = int(round(0.05 * N))
        else:
            raise ValueError(f"Unknown local seed rule: {seed_rule}")
        values.append({"seed_rule": seed_rule, "T": int(np.clip(T, 0, K))})
    return values


def _case_columns(row):
    return {
        "data_profile": str(row["data_profile"]),
        "generator_seed": int(row["generator_seed"]),
        "sample": int(row["sample"]),
        "N": int(row["N"]),
        "L": int(row["L"]),
        "K": int(row["K"]),
        "off_pct": float(row["off_pct"]),
        "active_pct": float(row["active_pct"]),
        "sigma": float(row["sigma"]),
        "P": float(row["P"]),
    }


def _matrix_key(row):
    return (
        str(row["data_profile"]),
        int(row["generator_seed"]),
        int(row["sample"]),
        int(row["N"]),
        int(row["L"]),
    )


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


def _global_by_active_pct(summary):
    if summary.empty:
        return pd.DataFrame()
    weighted = []
    group_cols = ["active_pct", "seed_rule", "max_swaps"]
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
            "seed_T_mean",
        ):
            row[col] = float(np.average(group[col].astype(float), weights=cases))
        row["cases"] = int(total)
        row["fraction_exact_p05"] = float(group["fraction_exact_p05"].min())
        weighted.append(row)
    return pd.DataFrame(weighted)


def _plot_seed_dependence(summary, out_dir, colors):
    from utils.plotting import use_agg_backend

    use_agg_backend()
    import matplotlib.pyplot as plt

    if summary.empty:
        return
    seed_summary = summary[summary["max_swaps"] == 0].copy()
    seed_summary["seed_T_over_N"] = seed_summary["seed_T_mean"] / seed_summary["N"]
    n_values = sorted(seed_summary["N"].unique())
    cols = min(2, max(1, len(n_values)))
    rows = int(np.ceil(len(n_values) / cols))
    fig, axes = plt.subplots(
        rows,
        cols,
        figsize=(5.4 * cols, 4.0 * rows),
        squeeze=False,
    )
    for ax, N in zip(axes.ravel(), n_values):
        chunk_N = seed_summary[seed_summary["N"] == N]
        for seed_rule in LOCAL_SEED_RULES:
            chunk = chunk_N[chunk_N["seed_rule"] == seed_rule].sort_values("active_pct")
            if chunk.empty:
                continue
            ax.plot(
                chunk["active_pct"],
                chunk["seed_T_over_N"],
                marker="o",
                color=colors.get(f"{seed_rule}_s0"),
                label=seed_rule,
            )
        ax.set_title(f"N={int(N)}")
        ax.set_xlabel("requested active K percentage")
        ax.set_ylabel("seed T / N")
        ax.grid(True, alpha=0.25)
    for ax in axes.ravel()[len(n_values) :]:
        ax.axis("off")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Local seed threshold dependence")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(Path(out_dir) / "local_threshold_seed_dependence.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_failure_diagnostics(diagnostics, out_dir):
    from utils.plotting import use_agg_backend

    use_agg_backend()
    import matplotlib.pyplot as plt

    if diagnostics.empty:
        return
    primary = diagnostics[diagnostics["seed_rule"] == "best_tested_T"]
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
    fig.suptitle("Best-tested-T local failure diagnostics")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(Path(out_dir) / "local_threshold_failure_diagnostics.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_runtime(summary, out_dir, colors):
    from utils.plotting import use_agg_backend

    use_agg_backend()
    import matplotlib.pyplot as plt

    primary = summary[summary["seed_rule"] == "best_tested_T"]
    if primary.empty:
        return
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    for max_swaps in sorted(primary["max_swaps"].unique()):
        chunk = primary[primary["max_swaps"] == max_swaps]
        grouped = chunk.groupby("N", as_index=False)["elapsed_seconds_mean"].mean()
        ax.plot(
            grouped["N"],
            grouped["elapsed_seconds_mean"],
            marker="o",
            color=colors.get(f"best_tested_T_s{int(max_swaps)}"),
            label=f"best T + {int(max_swaps)} swaps",
        )
    ax.set_xlabel("N")
    ax.set_ylabel("mean runtime seconds")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    ax.set_title("Local threshold runtime by N")
    fig.tight_layout()
    fig.savefig(Path(out_dir) / "local_threshold_runtime_by_N_K.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def _empirical_cdf(values):
    values = pd.Series(values, dtype=float).replace([np.inf, -np.inf], np.nan).dropna()
    values = np.sort(values.to_numpy())
    if len(values) == 0:
        return values, values
    cumulative = np.arange(1, len(values) + 1, dtype=float) / float(len(values))
    return values, cumulative


def _boolean_series(series):
    if series.dtype == object:
        return series.astype(str).str.lower().isin(("true", "1", "yes"))
    return series.astype(bool)


def _safe_ratio(numerator, denominator):
    numerator = float(numerator)
    denominator = float(denominator)
    if abs(denominator) <= np.finfo(float).eps:
        return np.nan
    return numerator / denominator


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
