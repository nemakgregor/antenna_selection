import math
import tarfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

from algorithms.h3_strong_weak import solve_h3_strong_weak
from algorithms.h3_threshold_local import (
    best_cyclic_threshold_window,
    refine_selection_by_swaps,
    threshold_window_selection,
)
from utils.data import generate_v_profile_from_rng
from utils.io import archive_csv_artifacts, atomic_write_csv


LARGE_SEED_RULES = ("best_cyclic_window", "T_0p05N", "strong_weak")
LARGE_MAX_SWAPS = (0, 1, 2)
HONEST_MAX_SWAPS = (0, 1)
BOUNDED_REFERENCE_DIR = Path(
    "results/local_threshold_large_gauss_L2_N1000_K500_750_cyclic_s100"
)

LARGE_LABELS = {
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

LARGE_COLORS = {
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


def run_large_cyclic_local_analysis(
    out_dir,
    N,
    K_values,
    profiles,
    generator_seeds,
    samples,
    L,
    sigma=1.0,
    P=1.0,
    max_swaps_values=LARGE_MAX_SWAPS,
    candidate_radius=None,
    honest_all_inactive=False,
    output_prefix="local_threshold_large",
    report_title="Large-N Cyclic Threshold Local Search Study",
    bounded_reference_dir=None,
    archive_csv=True,
    docs_path=None,
):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    total_cases = (
        len(profiles)
        * len(generator_seeds)
        * int(samples)
        * len(K_values)
    )
    case_no = 0
    N = int(N)
    L = int(L)
    for profile in profiles:
        for generator_seed in generator_seeds:
            rng = np.random.RandomState(int(generator_seed))
            for sample in range(int(samples)):
                V = generate_v_profile_from_rng(rng, N, L, str(profile))
                for K in K_values:
                    K = int(K)
                    if not (0 <= K <= N):
                        raise ValueError(f"K must satisfy 0 <= K <= N; got K={K}, N={N}.")
                    case_no += 1
                    K_off = N - K
                    case_base = {
                        "data_profile": str(profile),
                        "generator_seed": int(generator_seed),
                        "sample": int(sample),
                        "N": N,
                        "L": L,
                        "K": K,
                        "K_active": K,
                        "K_off": K_off,
                        "off_pct": 100.0 * float(K_off) / float(N) if N else 0.0,
                        "active_pct": 100.0 * float(K) / float(N) if N else 0.0,
                        "sigma": float(sigma),
                        "P": float(P),
                    }
                    if case_no == 1 or case_no % 25 == 0:
                        print(
                            "large cyclic local "
                            f"{case_no}/{total_cases}: N={N}, L={L}, profile={profile}, "
                            f"seed={int(generator_seed)}, sample={sample}, "
                            f"K_active={K}, K_off={K_off}",
                            flush=True,
                        )
                    rows.extend(
                        _large_local_rows_for_case(
                            V,
                            case_base,
                            max_swaps_values=max_swaps_values,
                            candidate_radius=candidate_radius,
                            honest_all_inactive=honest_all_inactive,
                        )
                    )

    runs = pd.DataFrame(rows)
    runs = _attach_large_references(runs)
    summary = build_large_local_summary(runs)
    best_t_summary = build_large_best_t_summary(runs)

    atomic_write_csv(runs, out_dir / f"{output_prefix}_runs.csv")
    atomic_write_csv(summary, out_dir / f"{output_prefix}_summary.csv")
    atomic_write_csv(best_t_summary, out_dir / f"{output_prefix}_best_t_summary.csv")

    write_large_local_plots(runs, summary, best_t_summary, out_dir, max_swaps_values=max_swaps_values)
    write_large_local_report(
        runs,
        summary,
        best_t_summary,
        out_dir,
        candidate_radius=candidate_radius,
        honest_all_inactive=honest_all_inactive,
        max_swaps_values=max_swaps_values,
        output_prefix=output_prefix,
        report_title=report_title,
        bounded_reference_dir=bounded_reference_dir,
        docs_path=docs_path,
    )

    archive_path = None
    if archive_csv:
        archive_path = archive_csv_artifacts(out_dir, remove_originals=True)

    return {
        "runs": runs,
        "summary": summary,
        "best_t_summary": best_t_summary,
        "archive_path": archive_path,
    }


def run_large_cyclic_honest_local_analysis(
    out_dir,
    N,
    K_values,
    profiles,
    generator_seeds,
    samples,
    L,
    sigma=1.0,
    P=1.0,
    archive_csv=True,
    docs_path=None,
    bounded_reference_dir=BOUNDED_REFERENCE_DIR,
):
    return run_large_cyclic_local_analysis(
        out_dir=out_dir,
        N=N,
        K_values=K_values,
        profiles=profiles,
        generator_seeds=generator_seeds,
        samples=samples,
        L=L,
        sigma=sigma,
        P=P,
        max_swaps_values=HONEST_MAX_SWAPS,
        candidate_radius=None,
        honest_all_inactive=True,
        output_prefix="local_threshold_large_honest",
        report_title="Honest All-Inactive One-Swap Large-N Study",
        bounded_reference_dir=bounded_reference_dir,
        archive_csv=archive_csv,
        docs_path=docs_path,
    )


def build_large_local_summary(runs):
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
                "K",
                "K_active",
                "K_off",
                "off_pct",
                "active_pct",
                "seed_rule",
                "max_swaps",
            ],
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
            fraction_best_cyclic_seed_mean=("fraction_best_cyclic_seed_u_g", "mean"),
            fraction_best_cyclic_seed_p05=("fraction_best_cyclic_seed_u_g", q(0.05)),
            fraction_best_cyclic_seed_p50=("fraction_best_cyclic_seed_u_g", q(0.50)),
            fraction_best_cyclic_seed_p95=("fraction_best_cyclic_seed_u_g", q(0.95)),
            fraction_best_observed_mean=("fraction_best_observed_u_g", "mean"),
            fraction_best_observed_p05=("fraction_best_observed_u_g", q(0.05)),
            fraction_best_observed_p50=("fraction_best_observed_u_g", q(0.50)),
            fraction_best_observed_p95=("fraction_best_observed_u_g", q(0.95)),
            best_observed_rate=("is_best_observed", "mean"),
            swaps_applied_mean=("swaps_applied", "mean"),
            candidate_radius_mean=("candidate_radius", "mean"),
            add_candidate_count_mean=("add_candidate_count", "mean"),
            evaluated_swap_count_mean=("evaluated_swap_count", "mean"),
            seed_elapsed_seconds_mean=("seed_elapsed_seconds", "mean"),
            local_elapsed_seconds_mean=("local_elapsed_seconds", "mean"),
            total_elapsed_seconds_mean=("total_elapsed_seconds", "mean"),
            total_elapsed_seconds_p95=("total_elapsed_seconds", q(0.95)),
        )
        .sort_values(["data_profile", "K", "seed_rule", "max_swaps"])
    )


def build_large_best_t_summary(runs):
    if runs.empty:
        return pd.DataFrame()
    best = runs[
        (runs["seed_rule"] == "best_cyclic_window")
        & (runs["max_swaps"] == 0)
    ].copy()
    if best.empty:
        return pd.DataFrame()
    best["T"] = best["seed_position"].astype(float)
    best["T_over_N"] = best["T"] / best["N"].astype(float)
    best["T_over_K"] = best["T"] / best["K_active"].replace(0, np.nan).astype(float)

    def q(value):
        return lambda data: data.quantile(value)

    return (
        best.groupby(["data_profile", "N", "L", "K", "K_active", "K_off", "off_pct", "active_pct"], as_index=False)
        .agg(
            cases=("T", "count"),
            T_mean=("T", "mean"),
            T_std=("T", "std"),
            T_p05=("T", q(0.05)),
            T_p50=("T", q(0.50)),
            T_p95=("T", q(0.95)),
            T_over_N_mean=("T_over_N", "mean"),
            T_over_N_p05=("T_over_N", q(0.05)),
            T_over_N_p50=("T_over_N", q(0.50)),
            T_over_N_p95=("T_over_N", q(0.95)),
            T_over_K_p50=("T_over_K", q(0.50)),
        )
        .sort_values(["data_profile", "K"])
    )


def write_large_local_plots(runs, summary, best_t_summary, out_dir, max_swaps_values=LARGE_MAX_SWAPS):
    if runs.empty:
        return

    from utils.plotting import use_agg_backend

    use_agg_backend()
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    K_values = sorted(runs["K_active"].unique())
    rule_specs = _large_rule_specs(max_swaps_values)

    cols = min(2, max(1, len(K_values)))
    plot_rows = int(math.ceil(len(K_values) / cols))
    fig, axes = plt.subplots(plot_rows, cols, figsize=(6.3 * cols, 4.3 * plot_rows), squeeze=False)
    for ax, K in zip(axes.ravel(), K_values):
        chunk_k = runs[runs["K_active"] == K]
        for rule_key, seed_rule, max_swaps in rule_specs:
            chunk = chunk_k[
                (chunk_k["seed_rule"] == seed_rule)
                & (chunk_k["max_swaps"] == max_swaps)
            ]
            values, cumulative = _empirical_cdf(chunk["u_g"])
            if len(values):
                ax.step(
                    cumulative,
                    values,
                    where="post",
                    color=LARGE_COLORS.get(rule_key),
                    label=LARGE_LABELS[rule_key],
                    linewidth=1.35 if max_swaps == 2 else 1.0,
                )
        off_pct = chunk_k["off_pct"].astype(float).mean()
        ax.set_title(f"K_active={int(K)} ({off_pct:.1f}% off)")
        ax.set_xlabel("cumulative fraction")
        ax.set_ylabel("raw U_G")
        ax.set_yscale("log")
        ax.grid(True, alpha=0.25)
    for ax in axes.ravel()[len(K_values) :]:
        ax.axis("off")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Large-N local threshold raw U_G CDF")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_dir / "large_raw_u_g_cdf_by_K.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(plot_rows, cols, figsize=(6.3 * cols, 4.3 * plot_rows), squeeze=False)
    for ax, K in zip(axes.ravel(), K_values):
        chunk_k = runs[runs["K_active"] == K]
        for rule_key, seed_rule, max_swaps in rule_specs:
            chunk = chunk_k[
                (chunk_k["seed_rule"] == seed_rule)
                & (chunk_k["max_swaps"] == max_swaps)
            ]
            values, cumulative = _empirical_cdf(chunk["fraction_best_cyclic_seed_u_g"])
            if len(values):
                ax.step(
                    cumulative,
                    values,
                    where="post",
                    color=LARGE_COLORS.get(rule_key),
                    label=LARGE_LABELS[rule_key],
                    linewidth=1.35 if max_swaps == 2 else 1.0,
                )
        ax.axhline(1.0, color="#666666", linewidth=0.9, linestyle="--")
        ax.set_title(f"K_active={int(K)}")
        ax.set_xlabel("cumulative fraction")
        ax.set_ylabel("U_G / cyclic seed U_G")
        ax.grid(True, alpha=0.25)
    for ax in axes.ravel()[len(K_values) :]:
        ax.axis("off")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Large-N local threshold fraction of cyclic seed")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_dir / "large_fraction_cyclic_seed_cdf_by_K.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(plot_rows, cols, figsize=(6.3 * cols, 4.3 * plot_rows), squeeze=False)
    for ax, K in zip(axes.ravel(), K_values):
        chunk_k = runs[runs["K_active"] == K]
        for rule_key, seed_rule, max_swaps in rule_specs:
            chunk = chunk_k[
                (chunk_k["seed_rule"] == seed_rule)
                & (chunk_k["max_swaps"] == max_swaps)
            ]
            values, cumulative = _empirical_cdf(chunk["fraction_best_observed_u_g"])
            if len(values):
                ax.step(
                    cumulative,
                    values,
                    where="post",
                    color=LARGE_COLORS.get(rule_key),
                    label=LARGE_LABELS[rule_key],
                    linewidth=1.35 if max_swaps == 2 else 1.0,
                )
        ax.set_title(f"K_active={int(K)}")
        ax.set_xlabel("cumulative fraction")
        ax.set_ylabel("U_G / best observed local U_G")
        ax.set_ylim(0.0, 1.03)
        ax.grid(True, alpha=0.25)
    for ax in axes.ravel()[len(K_values) :]:
        ax.axis("off")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Large-N local threshold fraction of best observed method")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_dir / "large_fraction_best_observed_cdf_by_K.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    for rule_key, seed_rule, max_swaps in rule_specs:
        chunk = summary[
            (summary["seed_rule"] == seed_rule)
            & (summary["max_swaps"] == max_swaps)
        ].sort_values("K_active")
        if chunk.empty:
            continue
        ax.plot(
            chunk["K_active"],
            chunk["fraction_best_cyclic_seed_mean"],
            marker="o",
            color=LARGE_COLORS.get(rule_key),
            label=LARGE_LABELS[rule_key],
        )
    ax.axhline(1.0, color="#666666", linewidth=0.9, linestyle="--")
    ax.set_xlabel("K_active")
    ax.set_ylabel("mean U_G / cyclic seed U_G")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7)
    ax.set_title("Mean quality relative to cyclic seed")
    fig.tight_layout()
    fig.savefig(out_dir / "large_mean_fraction_cyclic_seed_by_K.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    _plot_large_best_cyclic_t(runs, best_t_summary, out_dir)
    _plot_large_runtime(summary, out_dir, max_swaps_values=max_swaps_values)


def write_large_local_report(
    runs,
    summary,
    best_t_summary,
    out_dir,
    candidate_radius=None,
    honest_all_inactive=False,
    max_swaps_values=LARGE_MAX_SWAPS,
    output_prefix="local_threshold_large",
    report_title="Large-N Cyclic Threshold Local Search Study",
    bounded_reference_dir=None,
    docs_path=None,
):
    global_summary = _large_global_summary(summary)
    bounded_summary = _load_bounded_reference_summary(bounded_reference_dir)
    lines = [
        f"# {report_title}",
        "",
        "This report studies the threshold-window approach without brute-force exact enumeration.",
        "",
        "## K Semantics",
        "",
        "- `K` means active/kept antennas selected by the solver.",
        "- `K_off = N - K` is the number of disabled antennas.",
        "- For `N=1000`, `K=750` means `25% off`; `K=500` means `50% off`.",
        "",
        "## Setup",
        "",
        f"- Profiles: {', '.join(str(value) for value in sorted(runs['data_profile'].unique())) if not runs.empty else ''}",
        f"- N values: {', '.join(str(int(value)) for value in sorted(runs['N'].unique())) if not runs.empty else ''}",
        f"- L values: {', '.join(str(int(value)) for value in sorted(runs['L'].unique())) if not runs.empty else ''}",
        f"- K_active values: {', '.join(str(int(value)) for value in sorted(runs['K_active'].unique())) if not runs.empty else ''}",
        "- Seeds compared: best tested cyclic window, `T=round(0.05N)`, and strong/weak.",
        f"- Local search: greedy remove-one/add-one refinement by raw `U_G`, with swaps {', '.join(str(int(value)) for value in max_swaps_values)}.",
        "",
        "## Local Swap Scheme And Cost",
        "",
        "- A `1-swap` removes one currently active antenna and adds one currently inactive antenna, preserving exact `K`.",
        "- One local-search pass evaluates every `(remove, add)` pair from the current active set and the configured add-candidate pool, applies only the single pair with the largest positive `U_G` improvement, and stops if no pair improves `U_G`.",
        _local_search_scope_line(honest_all_inactive, candidate_radius),
        _local_search_asymptotic_line(honest_all_inactive),
        "- Best cyclic seed construction scans `N` cyclic windows before local refinement; `T=0.05N` and strong/weak construct one sorted-window seed.",
        _local_search_space_line(honest_all_inactive),
        "",
        "## Direct Answer",
        "",
    ]

    if runs.empty:
        lines.append("- No rows were produced.")
    else:
        best_observed = global_summary.sort_values("fraction_best_observed_mean", ascending=False).head(1)
        if not best_observed.empty:
            row = best_observed.iloc[0]
            key = f"{row['seed_rule']}_s{int(row['max_swaps'])}"
            lines.append(
                f"- Best observed tested method by mean `U_G`: `{LARGE_LABELS.get(key, key)}` "
                f"with mean fraction of best observed `{row['fraction_best_observed_mean']:.4f}`."
            )
        cyclic = global_summary[global_summary["seed_rule"] == "best_cyclic_window"].sort_values("max_swaps")
        zero = cyclic[cyclic["max_swaps"] == 0]
        best_swap = cyclic[cyclic["max_swaps"] == max(int(value) for value in max_swaps_values)]
        if not zero.empty and not best_swap.empty:
            swap_count = int(best_swap.iloc[0]["max_swaps"])
            gain = best_swap.iloc[0]["fraction_best_cyclic_seed_mean"] - zero.iloc[0]["fraction_best_cyclic_seed_mean"]
            lines.append(
                f"- Cyclic best + {swap_count} swap changes mean `U_G` by `{100.0 * gain:.3f}`% of the cyclic seed baseline."
            )
        t005 = global_summary[global_summary["seed_rule"] == "T_0p05N"].sort_values("max_swaps")
        if not t005.empty:
            best_t005 = t005.iloc[-1]
            lines.append(
                f"- `T=0.05N + {int(best_t005['max_swaps'])} swap` reaches mean fraction of cyclic seed "
                f"`{best_t005['fraction_best_cyclic_seed_mean']:.4f}`."
            )

    if honest_all_inactive:
        lines.extend(_honest_bounded_comparison_lines(summary, bounded_summary))

    if not best_t_summary.empty:
        lines.extend(
            [
                "",
                "## Best Cyclic T",
                "",
                "`T` is the cyclic start position in descending row-power order. `T=0` starts at the strongest row; larger `T` shifts the selected cyclic window.",
                "",
                "| K_active | off % | cases | T mean | T std | T p05 | T median | T p95 | T/N mean | T/N median |",
                "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for _, row in best_t_summary.iterrows():
            lines.append(
                "| "
                + " | ".join(
                    [
                        str(int(row["K_active"])),
                        _format_pct(row["off_pct"]),
                        str(int(row["cases"])),
                        _format_float(row["T_mean"]),
                        _format_float(row["T_std"]),
                        _format_float(row["T_p05"]),
                        _format_float(row["T_p50"]),
                        _format_float(row["T_p95"]),
                        f"{row['T_over_N_mean']:.4f}",
                        f"{row['T_over_N_p50']:.4f}",
                    ]
                )
                + " |"
            )

    lines.extend(
        [
            "",
            "## Method Summary",
            "",
            "| K_active | seed | swaps | cases | mean U_G | p05 U_G | mean frac cyclic seed | p05 frac cyclic seed | mean frac best observed | mean total runtime s | mean swaps applied |",
            "|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for _, row in summary.sort_values(["K_active", "seed_rule", "max_swaps"]).iterrows():
        key = f"{row['seed_rule']}_s{int(row['max_swaps'])}"
        lines.append(
            "| "
            + " | ".join(
                [
                    str(int(row["K_active"])),
                    LARGE_LABELS.get(key, key),
                    str(int(row["max_swaps"])),
                    str(int(row["cases"])),
                    _format_float(row["u_g_mean"]),
                    _format_float(row["u_g_p05"]),
                    f"{row['fraction_best_cyclic_seed_mean']:.4f}",
                    f"{row['fraction_best_cyclic_seed_p05']:.4f}",
                    f"{row['fraction_best_observed_mean']:.4f}",
                    _format_float(row["total_elapsed_seconds_mean"]),
                    f"{row['swaps_applied_mean']:.3f}",
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Plots",
            "",
            "![Raw U_G CDF](large_raw_u_g_cdf_by_K.png)",
            "",
            "![Fraction cyclic seed CDF](large_fraction_cyclic_seed_cdf_by_K.png)",
            "",
            "![Fraction best observed CDF](large_fraction_best_observed_cdf_by_K.png)",
            "",
            "![Mean fraction cyclic seed](large_mean_fraction_cyclic_seed_by_K.png)",
            "",
            "![Best cyclic T boxplot](large_best_cyclic_T_boxplot.png)",
            "",
            "![Best cyclic T over N](large_best_cyclic_T_over_N.png)",
            "",
            "![Runtime](large_runtime_by_method.png)",
            "",
            "## Artifacts",
            "",
            "- Detailed CSVs are packed in `csv_data.tar.gz` after the run.",
            f"- Main report: `{output_prefix}_report.md`.",
        ]
    )

    text = "\n".join(lines) + "\n"
    report_path = Path(out_dir) / f"{output_prefix}_report.md"
    report_path.write_text(text, encoding="utf-8")
    if docs_path is not None:
        docs_path = Path(docs_path)
        docs_path.parent.mkdir(parents=True, exist_ok=True)
        docs_path.write_text(text, encoding="utf-8")


def _large_local_rows_for_case(
    V,
    case_base,
    max_swaps_values,
    candidate_radius=None,
    honest_all_inactive=False,
):
    seed_specs = _large_seed_specs_for_case(
        V,
        case_base,
        candidate_radius,
        honest_all_inactive=honest_all_inactive,
    )
    rows = []
    for seed in seed_specs:
        for max_swaps in max_swaps_values:
            result = refine_selection_by_swaps(
                V,
                seed["x"],
                max_swaps=int(max_swaps),
                sigma=float(case_base["sigma"]),
                P=float(case_base["P"]),
                candidate_radius=seed["candidate_radius"],
                candidate_pool=seed["candidate_pool"],
                seed_position=int(seed["seed_position"]),
            )
            local_elapsed = float(result["elapsed_seconds"])
            total_elapsed = float(seed["seed_elapsed_seconds"]) + local_elapsed
            rows.append(
                {
                    **case_base,
                    "seed_rule": seed["seed_rule"],
                    "seed_label": seed["seed_label"],
                    "seed_position": int(seed["seed_position"]),
                    "seed_candidate_count": int(seed["candidate_count"]),
                    "seed_elapsed_seconds": float(seed["seed_elapsed_seconds"]),
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
                    "initial_subset": _subset_to_string(result["initial_subset"]),
                    "subset": _subset_to_string(result["subset"]),
                    "swap_history": result["swap_history"],
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
            )
    return rows


def _large_global_summary(summary):
    if summary.empty:
        return pd.DataFrame()
    rows = []
    for (seed_rule, max_swaps), group in summary.groupby(["seed_rule", "max_swaps"]):
        weights = group["cases"].astype(float)
        total = float(weights.sum())
        if total <= 0:
            continue
        row = {
            "seed_rule": seed_rule,
            "max_swaps": int(max_swaps),
            "cases": int(total),
        }
        for col in (
            "u_g_mean",
            "fraction_best_cyclic_seed_mean",
            "fraction_best_observed_mean",
            "best_observed_rate",
            "swaps_applied_mean",
            "total_elapsed_seconds_mean",
        ):
            row[col] = float(np.average(group[col].astype(float), weights=weights))
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["seed_rule", "max_swaps"])


def _load_bounded_reference_summary(bounded_reference_dir):
    if bounded_reference_dir is None:
        return pd.DataFrame()
    base_dir = Path(bounded_reference_dir)
    csv_name = "local_threshold_large_summary.csv"
    csv_path = base_dir / csv_name
    if csv_path.exists():
        return pd.read_csv(csv_path)
    archive_path = base_dir / "csv_data.tar.gz"
    if not archive_path.exists():
        return pd.DataFrame()
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive.getmembers():
            if Path(member.name).name == csv_name:
                handle = archive.extractfile(member)
                if handle is not None:
                    return pd.read_csv(handle)
    return pd.DataFrame()


def _honest_bounded_comparison_lines(honest_summary, bounded_summary):
    lines = [
        "",
        "## Honest Versus Previous Bounded Local Search",
        "",
    ]
    if bounded_summary.empty:
        lines.append("- Previous bounded-neighborhood summary was not found, so no direct comparison table was generated.")
        return lines

    rows = []
    for _, honest in honest_summary[honest_summary["max_swaps"] == 1].iterrows():
        match = bounded_summary[
            (bounded_summary["data_profile"].astype(str) == str(honest["data_profile"]))
            & (bounded_summary["N"].astype(int) == int(honest["N"]))
            & (bounded_summary["L"].astype(int) == int(honest["L"]))
            & (bounded_summary["K_active"].astype(int) == int(honest["K_active"]))
            & (bounded_summary["seed_rule"].astype(str) == str(honest["seed_rule"]))
            & (bounded_summary["max_swaps"].astype(int) == 1)
        ]
        if match.empty:
            continue
        bounded = match.iloc[0]
        rows.append(
            {
                "K_active": int(honest["K_active"]),
                "seed_rule": str(honest["seed_rule"]),
                "honest_fraction": float(honest["fraction_best_cyclic_seed_mean"]),
                "bounded_fraction": float(bounded["fraction_best_cyclic_seed_mean"]),
                "delta_fraction": float(honest["fraction_best_cyclic_seed_mean"])
                - float(bounded["fraction_best_cyclic_seed_mean"]),
                "honest_runtime": float(honest["total_elapsed_seconds_mean"]),
                "bounded_runtime": float(bounded["total_elapsed_seconds_mean"]),
                "runtime_ratio": _safe_ratio(
                    float(honest["total_elapsed_seconds_mean"]),
                    float(bounded["total_elapsed_seconds_mean"]),
                ),
                "honest_pairs": float(honest["evaluated_swap_count_mean"]),
                "bounded_pairs": float(bounded["evaluated_swap_count_mean"]),
            }
        )

    if not rows:
        lines.append("- Previous bounded-neighborhood summary did not contain matching `max_swaps=1` rows.")
        return lines

    comparison = pd.DataFrame(rows)
    cyclic = comparison[comparison["seed_rule"] == "best_cyclic_window"]
    t005 = comparison[comparison["seed_rule"] == "T_0p05N"]
    strong = comparison[comparison["seed_rule"] == "strong_weak"]
    if not cyclic.empty:
        lines.append(
            f"- `cyclic best + 1 honest swap` changes mean fraction versus bounded by "
            f"`{cyclic['delta_fraction'].mean():.6f}` on average across K values."
        )
    if not t005.empty:
        lines.append(
            f"- `T=0.05N + 1 honest swap` changes mean fraction versus bounded by "
            f"`{t005['delta_fraction'].mean():.6f}` on average across K values."
        )
    if not strong.empty:
        lines.append(
            f"- `strong/weak + 1 honest swap` changes mean fraction versus bounded by "
            f"`{strong['delta_fraction'].mean():.6f}` on average across K values."
        )
    lines.extend(
        [
            "",
            "| K_active | seed | bounded frac cyclic seed | honest frac cyclic seed | delta | bounded runtime s | honest runtime s | runtime ratio | bounded pairs | honest pairs |",
            "|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for _, row in comparison.sort_values(["K_active", "seed_rule"]).iterrows():
        key = f"{row['seed_rule']}_s1"
        lines.append(
            "| "
            + " | ".join(
                [
                    str(int(row["K_active"])),
                    LARGE_LABELS.get(key, key),
                    f"{row['bounded_fraction']:.4f}",
                    f"{row['honest_fraction']:.4f}",
                    f"{row['delta_fraction']:.6f}",
                    _format_float(row["bounded_runtime"]),
                    _format_float(row["honest_runtime"]),
                    _format_float(row["runtime_ratio"]),
                    _format_float(row["bounded_pairs"]),
                    _format_float(row["honest_pairs"]),
                ]
            )
            + " |"
        )
    return lines


def _large_seed_specs_for_case(V, case_base, candidate_radius=None, honest_all_inactive=False):
    N = int(case_base["N"])
    K = int(case_base["K_active"])
    sigma = float(case_base["sigma"])
    P = float(case_base["P"])
    radius = None if honest_all_inactive else _candidate_radius(K, candidate_radius)
    if honest_all_inactive:
        order = None
    else:
        row_powers = np.sum(np.abs(V) ** 2, axis=1).real
        order = np.argsort(row_powers)[::-1]

    started_at = time.perf_counter()
    cyclic = best_cyclic_threshold_window(V, K, sigma=sigma, P=P)
    cyclic_pool = None if honest_all_inactive else _cyclic_boundary_add_pool(order, cyclic["x"], int(cyclic["T"]), K, radius)
    cyclic_elapsed = time.perf_counter() - started_at

    started_at = time.perf_counter()
    T_005 = int(np.clip(round(0.05 * N), 0, max(0, N - K)))
    formula_x = threshold_window_selection(V, K, T_005)
    formula_pool = None if honest_all_inactive else _linear_boundary_add_pool(order, formula_x, T_005, K, radius)
    formula_elapsed = time.perf_counter() - started_at

    started_at = time.perf_counter()
    strong_x = solve_h3_strong_weak(V, K, sigma=sigma, P=P)
    off_count = N - K
    strong_position = off_count - off_count // 2
    strong_pool = None if honest_all_inactive else _linear_boundary_add_pool(order, strong_x, strong_position, K, radius)
    strong_elapsed = time.perf_counter() - started_at

    return [
        {
            "seed_rule": "best_cyclic_window",
            "seed_label": LARGE_LABELS["best_cyclic_window_s0"],
            "seed_position": int(cyclic["T"]),
            "candidate_count": int(cyclic["candidate_count"]),
            "candidate_kind": "cyclic_threshold_window",
            "candidate_radius": radius,
            "candidate_pool": cyclic_pool,
            "seed_elapsed_seconds": cyclic_elapsed,
            "x": cyclic["x"],
        },
        {
            "seed_rule": "T_0p05N",
            "seed_label": LARGE_LABELS["T_0p05N_s0"],
            "seed_position": int(T_005),
            "candidate_count": 1,
            "candidate_kind": "threshold_window_T_0p05N",
            "candidate_radius": radius,
            "candidate_pool": formula_pool,
            "seed_elapsed_seconds": formula_elapsed,
            "x": formula_x,
        },
        {
            "seed_rule": "strong_weak",
            "seed_label": LARGE_LABELS["strong_weak_s0"],
            "seed_position": int(strong_position),
            "candidate_count": 1,
            "candidate_kind": "strong_weak",
            "candidate_radius": radius,
            "candidate_pool": strong_pool,
            "seed_elapsed_seconds": strong_elapsed,
            "x": strong_x,
        },
    ]


def _attach_large_references(runs):
    if runs.empty:
        return runs
    case_keys = [
        "data_profile",
        "generator_seed",
        "sample",
        "N",
        "L",
        "K",
        "K_active",
        "K_off",
        "sigma",
        "P",
    ]
    ref = runs[
        (runs["seed_rule"] == "best_cyclic_window")
        & (runs["max_swaps"] == 0)
    ][case_keys + ["u_g"]].rename(columns={"u_g": "best_cyclic_seed_u_g"})
    best = (
        runs.groupby(case_keys, as_index=False)["u_g"]
        .max()
        .rename(columns={"u_g": "best_observed_u_g"})
    )
    merged = runs.merge(ref, on=case_keys, how="left").merge(best, on=case_keys, how="left")
    merged["fraction_best_cyclic_seed_u_g"] = (
        merged["u_g"].astype(float) / merged["best_cyclic_seed_u_g"].astype(float)
    )
    merged["fraction_best_observed_u_g"] = (
        merged["u_g"].astype(float) / merged["best_observed_u_g"].astype(float)
    )
    merged["is_best_observed"] = (
        merged["u_g"].astype(float) >= merged["best_observed_u_g"].astype(float) - 1e-9
    )
    return merged


def _plot_large_best_cyclic_t(runs, best_t_summary, out_dir):
    if runs.empty:
        return

    from utils.plotting import use_agg_backend

    use_agg_backend()
    import matplotlib.pyplot as plt

    best = runs[
        (runs["seed_rule"] == "best_cyclic_window")
        & (runs["max_swaps"] == 0)
    ].copy()
    if best.empty:
        return
    K_values = sorted(best["K_active"].unique())

    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    data = [best[best["K_active"] == K]["seed_position"].astype(float) for K in K_values]
    ax.boxplot(data, tick_labels=[str(int(value)) for value in K_values], showmeans=True)
    ax.set_xlabel("K_active")
    ax.set_ylabel("best cyclic T")
    ax.grid(True, axis="y", alpha=0.25)
    ax.set_title("Best cyclic T distribution")
    fig.tight_layout()
    fig.savefig(Path(out_dir) / "large_best_cyclic_T_boxplot.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, len(K_values), figsize=(5.4 * len(K_values), 4.0), squeeze=False)
    for ax, K in zip(axes.ravel(), K_values):
        chunk = best[best["K_active"] == K]
        ax.hist(chunk["seed_position"], bins=30, color=LARGE_COLORS["best_cyclic_window_s0"], alpha=0.86)
        ax.set_title(f"K_active={int(K)}")
        ax.set_xlabel("best cyclic T")
        ax.set_ylabel("count")
        ax.grid(True, alpha=0.25)
    fig.suptitle("Best cyclic T histogram")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(Path(out_dir) / "large_best_cyclic_T_hist.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    if best_t_summary.empty:
        return
    fig, ax = plt.subplots(figsize=(7.8, 4.6))
    ax.plot(
        best_t_summary["K_active"],
        best_t_summary["T_over_N_p50"],
        marker="o",
        color=LARGE_COLORS["best_cyclic_window_s0"],
        label="median best T/N",
    )
    ax.fill_between(
        best_t_summary["K_active"].astype(float).to_numpy(),
        best_t_summary["T_over_N_p05"].astype(float).to_numpy(),
        best_t_summary["T_over_N_p95"].astype(float).to_numpy(),
        color=LARGE_COLORS["best_cyclic_window_s0"],
        alpha=0.18,
        label="p05..p95",
    )
    ax.set_xlabel("K_active")
    ax.set_ylabel("best cyclic T / N")
    ax.set_ylim(0.0, 1.0)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    ax.set_title("Best cyclic T/N by K")
    fig.tight_layout()
    fig.savefig(Path(out_dir) / "large_best_cyclic_T_over_N.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_large_runtime(summary, out_dir, max_swaps_values=LARGE_MAX_SWAPS):
    if summary.empty:
        return

    from utils.plotting import use_agg_backend

    use_agg_backend()
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    for rule_key, seed_rule, max_swaps in _large_rule_specs(max_swaps_values):
        chunk = summary[
            (summary["seed_rule"] == seed_rule)
            & (summary["max_swaps"] == max_swaps)
        ].sort_values("K_active")
        if chunk.empty:
            continue
        ax.plot(
            chunk["K_active"],
            chunk["total_elapsed_seconds_mean"],
            marker="o",
            color=LARGE_COLORS.get(rule_key),
            label=LARGE_LABELS[rule_key],
        )
    ax.set_xlabel("K_active")
    ax.set_ylabel("mean total runtime seconds")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7)
    ax.set_title("Large-N local threshold runtime")
    fig.tight_layout()
    fig.savefig(Path(out_dir) / "large_runtime_by_method.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def _candidate_radius(K, candidate_radius):
    if candidate_radius is None:
        return int(max(8, math.ceil(0.05 * int(K))))
    return int(max(0, candidate_radius))


def _candidate_radius_text(candidate_radius):
    if candidate_radius is None:
        return "default max(8, ceil(0.05K))"
    return str(int(candidate_radius))


def _local_search_scope_line(honest_all_inactive, candidate_radius):
    if honest_all_inactive:
        return "- Honest mode uses all inactive antennas as add candidates, so `A = N - K`."
    return (
        "- The add-candidate pool is fixed at seed time and is restricted to inactive "
        f"antennas near the sorted-window boundaries; candidate radius: `{_candidate_radius_text(candidate_radius)}`."
    )


def _local_search_asymptotic_line(honest_all_inactive):
    if honest_all_inactive:
        return "- One honest pass costs `O(K * (N-K) * L^3)` time after row Gram matrices are built."
    return "- With `S` greedy passes, `A` add candidates, and `L` streams, local refinement costs `O(S * K * A * L^3)` time after row Gram matrices are built."


def _local_search_space_line(honest_all_inactive):
    if honest_all_inactive:
        return "- Extra working space is `O(N * L^2 + K + (N-K) + L^2)` for row Gram matrices, active/add sets, and the current Gram matrix."
    return "- Extra working space is `O(N * L^2 + K + A + L^2)` for row Gram matrices, active/add sets, and the current Gram matrix."


def _cyclic_boundary_add_pool(order, x, start, K, radius):
    active = {int(value) for value in np.flatnonzero(np.asarray(x, dtype=int))}
    N = len(order)
    if int(K) >= N or int(radius) <= 0:
        return tuple()
    ranks = []
    start = int(start) % N
    end = (start + int(K) - 1) % N
    for offset in range(1, int(radius) + 1):
        ranks.append((start - offset) % N)
        ranks.append((end + offset) % N)
    return _pool_from_ranks(order, active, ranks)


def _linear_boundary_add_pool(order, x, start, K, radius):
    active = {int(value) for value in np.flatnonzero(np.asarray(x, dtype=int))}
    if int(radius) <= 0:
        return tuple()
    N = len(order)
    start = int(start)
    stop = int(start) + int(K)
    ranks = list(range(max(0, start - int(radius)), min(N, start + int(radius))))
    ranks.extend(range(max(0, stop - int(radius)), min(N, stop + int(radius))))
    return _pool_from_ranks(order, active, ranks)


def _pool_from_ranks(order, active, ranks):
    pool = []
    seen = set()
    N = len(order)
    for rank in ranks:
        if not (0 <= int(rank) < N):
            continue
        index = int(order[int(rank)])
        if index not in active and index not in seen:
            pool.append(index)
            seen.add(index)
    return tuple(pool)


def _large_rule_specs(max_swaps_values=LARGE_MAX_SWAPS):
    return [
        (f"{seed_rule}_s{max_swaps}", seed_rule, max_swaps)
        for seed_rule in LARGE_SEED_RULES
        for max_swaps in max_swaps_values
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


def _subset_to_string(subset):
    return " ".join(str(int(value)) for value in subset)


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
