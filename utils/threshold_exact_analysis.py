from pathlib import Path

import numpy as np
import pandas as pd


EXACT_K_RULES = (
    ("exact_best", "exact best"),
    ("best_tested_T", "best tested T"),
    ("T_0p05N", "T=0.05N"),
    ("strong_weak", "strong/weak H3"),
)

HISTORICAL_ACTIVE_K_NOTE = (
    "> Historical K semantics note: this report uses active-K semantics. "
    "Here `K` is the number of selected/kept antennas, not the number turned off. "
    "A `25% active` or `K=0.25N` case means `75% off`, not the real `25% off` task. "
    "For real off-percent experiments, `25% off => K_active=0.75N` and "
    "`50% off => K_active=0.50N`."
)


def write_threshold_exact_k_analysis(exact_runs, formula_runs, out_dir):
    out_dir = Path(out_dir)
    exact_runs = exact_runs.copy()
    formula_runs = formula_runs.copy()
    if exact_runs.empty or formula_runs.empty:
        return {}

    rule_runs = build_exact_k_rule_runs(exact_runs, formula_runs)
    rule_summary = build_exact_k_rule_summary(rule_runs)
    best_t_summary = build_exact_k_best_t_summary(exact_runs)

    rule_runs.to_csv(out_dir / "exact_k_pct_rule_runs.csv", index=False)
    rule_summary.to_csv(out_dir / "exact_k_pct_rule_summary.csv", index=False)
    best_t_summary.to_csv(out_dir / "exact_k_pct_best_t_summary.csv", index=False)

    write_exact_k_plots(rule_runs, best_t_summary, out_dir)
    write_exact_k_report(rule_summary, best_t_summary, out_dir)

    return {
        "rule_runs": rule_runs,
        "rule_summary": rule_summary,
        "best_t_summary": best_t_summary,
    }


def build_exact_k_rule_runs(exact_runs, formula_runs):
    exact_completed = exact_runs[_boolean_series(exact_runs["exact_completed"])].copy()
    exact_completed["formula"] = "exact_best"
    exact_completed["formula_label"] = "exact best"
    exact_completed["formula_T"] = exact_completed["exact_window_T"]
    exact_completed["u_g"] = exact_completed["exact_u_g"]
    exact_completed["fraction_exact_u_g"] = 1.0
    exact_completed["gap_to_exact_pct"] = 0.0
    exact_completed["is_exact_best"] = True

    keep_cols = [
        "data_profile",
        "generator_seed",
        "sample",
        "N",
        "L",
        "K",
        "off_pct",
        "active_pct",
        "sigma",
        "P",
        "formula",
        "formula_label",
        "formula_T",
        "u_g",
        "fraction_exact_u_g",
        "gap_to_exact_pct",
        "is_exact_best",
    ]
    exact_part = exact_completed[keep_cols].copy()

    formula_part = formula_runs[
        formula_runs["formula"].isin(["best_tested_T", "T_0p05N", "strong_weak"])
    ][keep_cols].copy()
    combined = pd.concat([exact_part, formula_part], ignore_index=True)
    combined["active_pct_actual"] = 100.0 * combined["K"].astype(float) / combined[
        "N"
    ].astype(float)
    combined["T_over_N"] = combined["formula_T"].astype(float) / combined["N"].astype(
        float
    )
    combined["T_over_K"] = combined["formula_T"].astype(float) / combined["K"].astype(
        float
    )
    label_map = dict(EXACT_K_RULES)
    combined["formula_label"] = combined["formula"].map(label_map).fillna(
        combined["formula_label"]
    )
    return combined.sort_values(
        ["active_pct", "N", "K", "formula", "generator_seed", "sample"]
    )


def build_exact_k_rule_summary(rule_runs):
    def q(value):
        return lambda data: data.quantile(value)

    return (
        rule_runs.groupby(["active_pct", "formula", "formula_label"], as_index=False)
        .agg(
            cases=("u_g", "count"),
            active_pct_actual_mean=("active_pct_actual", "mean"),
            N_min=("N", "min"),
            N_max=("N", "max"),
            K_min=("K", "min"),
            K_max=("K", "max"),
            T_mean=("formula_T", "mean"),
            T_p50=("formula_T", q(0.50)),
            T_over_N_mean=("T_over_N", "mean"),
            T_over_K_mean=("T_over_K", "mean"),
            u_g_mean=("u_g", "mean"),
            u_g_p05=("u_g", q(0.05)),
            u_g_p50=("u_g", q(0.50)),
            u_g_p95=("u_g", q(0.95)),
            fraction_exact_mean=("fraction_exact_u_g", "mean"),
            fraction_exact_p05=("fraction_exact_u_g", q(0.05)),
            fraction_exact_p50=("fraction_exact_u_g", q(0.50)),
            fraction_exact_p95=("fraction_exact_u_g", q(0.95)),
            exact_rate=("is_exact_best", "mean"),
        )
        .sort_values(["active_pct", "formula"])
    )


def build_exact_k_best_t_summary(exact_runs):
    def q(value):
        return lambda data: data.quantile(value)

    frame = exact_runs.copy()
    frame["active_pct_actual"] = 100.0 * frame["K"].astype(float) / frame["N"].astype(
        float
    )
    frame["best_tested_T_over_N"] = frame["best_tested_T"].astype(float) / frame[
        "N"
    ].astype(float)
    frame["best_tested_T_over_K"] = frame["best_tested_T"].astype(float) / frame[
        "K"
    ].astype(float)
    frame["T_0p05N_over_N"] = 0.05
    frame["T_0p05N_over_K"] = 0.05 * frame["N"].astype(float) / frame["K"].astype(
        float
    )
    return (
        frame.groupby(["active_pct"], as_index=False)
        .agg(
            cases=("best_tested_T", "count"),
            active_pct_actual_mean=("active_pct_actual", "mean"),
            best_T_mean=("best_tested_T", "mean"),
            best_T_p50=("best_tested_T", q(0.50)),
            best_T_p95=("best_tested_T", q(0.95)),
            best_T_over_N_mean=("best_tested_T_over_N", "mean"),
            best_T_over_K_mean=("best_tested_T_over_K", "mean"),
            T_0p05N_over_N_mean=("T_0p05N_over_N", "mean"),
            T_0p05N_over_K_mean=("T_0p05N_over_K", "mean"),
            threshold_fraction_exact_mean=("best_tested_fraction_exact_u_g", "mean"),
            threshold_fraction_exact_p05=("best_tested_fraction_exact_u_g", q(0.05)),
            threshold_exact_rate=(
                "best_tested_fraction_exact_u_g",
                lambda values: (values >= 1.0 - 1e-9).mean(),
            ),
            exact_is_threshold_window_rate=("exact_is_threshold_window", "mean"),
        )
        .sort_values("active_pct")
    )


def write_exact_k_plots(rule_runs, best_t_summary, out_dir):
    from utils.plotting import use_agg_backend

    use_agg_backend()
    import matplotlib.pyplot as plt

    rule_order = [rule for rule, _label in EXACT_K_RULES]
    label_map = dict(EXACT_K_RULES)
    active_pcts = sorted(rule_runs["active_pct"].unique())
    colors = {
        "exact_best": "#000000",
        "best_tested_T": "#0072B2",
        "T_0p05N": "#E69F00",
        "strong_weak": "#D55E00",
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
        chunk_pct = rule_runs[rule_runs["active_pct"] == active_pct]
        actual_pct = (
            float(chunk_pct["active_pct_actual"].mean()) if not chunk_pct.empty else active_pct
        )
        for rule in rule_order:
            chunk = chunk_pct[chunk_pct["formula"] == rule]
            x_values, y_values = _empirical_cdf(chunk["u_g"])
            if len(x_values):
                ax.step(
                    y_values,
                    x_values,
                    where="post",
                    color=colors.get(rule),
                    label=label_map[rule],
                )
        ax.set_title(f"requested K={active_pct:g}% (actual mean {actual_pct:.2f}%)")
        ax.set_xlabel("cumulative fraction")
        ax.set_ylabel("raw U_G")
        ax.set_yscale("log")
        ax.grid(True, alpha=0.25)
    for ax in axes.ravel()[len(active_pcts) :]:
        ax.axis("off")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Raw U_G CDF by active K percentage")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_dir / "exact_k_pct_raw_u_g_cdf.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fraction_rules = ["best_tested_T", "T_0p05N", "strong_weak"]
    pct_colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(active_pcts)))
    color_by_pct = dict(zip(active_pcts, pct_colors))
    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.4), squeeze=False)
    for ax, rule in zip(axes.ravel(), fraction_rules):
        chunk_rule = rule_runs[rule_runs["formula"] == rule]
        for active_pct in active_pcts:
            chunk = chunk_rule[chunk_rule["active_pct"] == active_pct]
            x_values, y_values = _empirical_cdf(chunk["fraction_exact_u_g"])
            if len(x_values):
                ax.step(
                    y_values,
                    x_values,
                    where="post",
                    color=color_by_pct[active_pct],
                    label=f"{active_pct:g}%",
                )
        ax.set_title(label_map[rule])
        ax.set_xlabel("cumulative fraction")
        ax.set_ylabel("U_G / exact U_G")
        ax.set_ylim(0.0, 1.05)
        ax.grid(True, alpha=0.25)
    axes[0, 0].legend(title="requested K/N", fontsize=8)
    fig.suptitle("Fraction of exact U_G CDF by active K percentage")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(
        out_dir / "exact_k_pct_fraction_exact_cdf.png",
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    summary = build_exact_k_rule_summary(rule_runs)
    for rule in ["exact_best", "best_tested_T", "T_0p05N", "strong_weak"]:
        chunk = summary[summary["formula"] == rule].sort_values("active_pct")
        if chunk.empty:
            continue
        ax.plot(
            chunk["active_pct"],
            chunk["fraction_exact_mean"],
            marker="o",
            label=label_map[rule],
        )
    ax.set_xlabel("requested active K percentage")
    ax.set_ylabel("mean U_G / exact U_G")
    ax.set_ylim(0.0, 1.04)
    ax.grid(True, alpha=0.25)
    ax.legend()
    ax.set_title("Mean fraction of exact U_G by active K percentage")
    fig.tight_layout()
    fig.savefig(
        out_dir / "exact_k_pct_fraction_by_active_pct.png",
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4), squeeze=False)
    ax = axes[0, 0]
    ax.plot(
        best_t_summary["active_pct"],
        best_t_summary["best_T_over_N_mean"],
        marker="o",
        label="best tested T/N",
    )
    ax.plot(
        best_t_summary["active_pct"],
        best_t_summary["T_0p05N_over_N_mean"],
        marker="o",
        label="T=0.05N",
    )
    ax.set_xlabel("requested active K percentage")
    ax.set_ylabel("T / N")
    ax.grid(True, alpha=0.25)
    ax.legend()
    ax.set_title("Threshold scale vs N")

    ax = axes[0, 1]
    ax.plot(
        best_t_summary["active_pct"],
        best_t_summary["best_T_over_K_mean"],
        marker="o",
        label="best tested T/K",
    )
    ax.plot(
        best_t_summary["active_pct"],
        best_t_summary["T_0p05N_over_K_mean"],
        marker="o",
        label="T=0.05N as T/K",
    )
    ax.set_xlabel("requested active K percentage")
    ax.set_ylabel("T / K")
    ax.grid(True, alpha=0.25)
    ax.legend()
    ax.set_title("Threshold scale vs K")
    fig.suptitle("Best tested T compared with T=0.05N")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(
        out_dir / "exact_k_pct_best_T_dependence.png",
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(fig)


def write_exact_k_report(rule_summary, best_t_summary, out_dir):
    lines = [
        "# Exact Study: K-Percentage Dependence",
        "",
        HISTORICAL_ACTIVE_K_NOTE,
        "",
        "This report compares exact best, best tested threshold, `T=0.05N`, and strong/weak H3 across requested active K percentages.",
        "For small `N`, different requested percentages can round to the same integer `K`; the tables include the mean actual `100*K/N` percentage.",
        "",
        "## Direct Answer",
        "",
    ]

    best = rule_summary[rule_summary["formula"] == "best_tested_T"]
    t005 = rule_summary[rule_summary["formula"] == "T_0p05N"]
    strong = rule_summary[rule_summary["formula"] == "strong_weak"]
    if not best.empty:
        lines.append(
            f"- Best tested threshold stays close to exact across K percentages: mean fraction range `{best['fraction_exact_mean'].min():.4f}..{best['fraction_exact_mean'].max():.4f}`."
        )
    if not t005.empty:
        lines.append(
            f"- `T=0.05N` is weaker than best tested T: mean fraction range `{t005['fraction_exact_mean'].min():.4f}..{t005['fraction_exact_mean'].max():.4f}`."
        )
    if not strong.empty:
        lines.append(
            f"- Strong/weak H3 remains far below exact on this exact small-N Gaussian grid: mean fraction range `{strong['fraction_exact_mean'].min():.4f}..{strong['fraction_exact_mean'].max():.4f}`."
        )
    if not best_t_summary.empty:
        lines.append(
            f"- Best tested `T/N` is small: mean range `{best_t_summary['best_T_over_N_mean'].min():.4f}..{best_t_summary['best_T_over_N_mean'].max():.4f}`, while `T=0.05N` is fixed at `0.0500`."
        )

    lines.extend(
        [
            "",
            "## Rule Quality By Requested K Percentage",
            "",
            "| requested K% | actual K% mean | rule | mean fraction exact | p05 fraction | exact rate | mean raw U_G |",
            "|---:|---:|---|---:|---:|---:|---:|",
        ]
    )
    for _, row in rule_summary.sort_values(["active_pct", "formula"]).iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    f"{row['active_pct']:.0f}",
                    f"{row['active_pct_actual_mean']:.2f}",
                    str(row["formula_label"]),
                    f"{row['fraction_exact_mean']:.4f}",
                    f"{row['fraction_exact_p05']:.4f}",
                    f"{row['exact_rate']:.1%}",
                    _format_float(row["u_g_mean"]),
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Best T Dependence",
            "",
            "| requested K% | actual K% mean | best T mean | best T p50 | best T/N mean | best T/K mean | T=0.05N as T/K | exact-window rate |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for _, row in best_t_summary.iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    f"{row['active_pct']:.0f}",
                    f"{row['active_pct_actual_mean']:.2f}",
                    f"{row['best_T_mean']:.3f}",
                    f"{row['best_T_p50']:.3f}",
                    f"{row['best_T_over_N_mean']:.4f}",
                    f"{row['best_T_over_K_mean']:.4f}",
                    f"{row['T_0p05N_over_K_mean']:.4f}",
                    f"{row['exact_is_threshold_window_rate']:.1%}",
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Plots",
            "",
            "![Raw U_G CDF by K percentage](exact_k_pct_raw_u_g_cdf.png)",
            "",
            "![Fraction of exact CDF by K percentage](exact_k_pct_fraction_exact_cdf.png)",
            "",
            "![Fraction by active percentage](exact_k_pct_fraction_by_active_pct.png)",
            "",
            "![Best T dependence](exact_k_pct_best_T_dependence.png)",
            "",
            "## Artifacts",
            "",
            "- `exact_k_pct_rule_runs.csv`",
            "- `exact_k_pct_rule_summary.csv`",
            "- `exact_k_pct_best_t_summary.csv`",
        ]
    )
    (out_dir / "threshold_exact_k_pct_analysis.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


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


def _format_float(value):
    if pd.isna(value):
        return ""
    value = float(value)
    if abs(value) >= 1e6 or (0 < abs(value) < 1e-3):
        return f"{value:.3e}"
    return f"{value:.3f}"
