from utils.plotting import use_agg_backend

use_agg_backend()
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ALGORITHM_COLORS = {
    "H1": "#0072B2",
    "H2": "#D55E00",
    "H3": "#009E73",
    "FrameOnly-Gen": "#9467BD",
    "CapWindow-Gen": "#1F77B4",
    "CapWindowFull-Gen": "#005F73",
    "Frame-Gen": "#E377C2",
    "S-threshold-Gen": "#E69F00",
    "H3ThresholdT123-Gen": "#B07AA1",
    "CapSubmod-Gen": "#17BECF",
    "CapSubmodPort-Gen": "#D62728",
    "Coutino": "#000000",
    "BackwardTrueGreedy": "#4D4D4D",
    "CoutinoSchur-Gen": "#7F3C8D",
    "ThreshDOpt-Gen": "#11A579",
    "ThreshWLogdet-Gen": "#3969AC",
    "ThreshDOptSwap-Gen": "#F2B701",
}


ALGORITHM_LINESTYLES = {
    "H3ThresholdT123-Gen": "-.",
    "CapSubmod-Gen": "-",
    "CapSubmodPort-Gen": "-.",
}


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
                    linewidth=1.1,
                    color=colors[name],
                    linestyle=ALGORITHM_LINESTYLES.get(name, "-"),
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


def pareto_front(points):
    points = points.sort_values(["elapsed_p50", "u_g_db_mean"], ascending=[True, False])
    front = []
    best_accuracy = -np.inf
    for _, row in points.iterrows():
        accuracy = float(row["u_g_db_mean"])
        if accuracy > best_accuracy + 1e-12:
            front.append(row)
            best_accuracy = accuracy
    return front


def _load_trace_upper_bounds(out_dir):
    upper_bound_path = out_dir / "cap_window_full_trace_upper_bound.csv"
    if not upper_bound_path.exists():
        return None

    upper_bounds = pd.read_csv(upper_bound_path)
    required = {"off_pct", "K", "trace_cap_upper_bound_db"}
    if not required.issubset(upper_bounds.columns):
        return None

    return (
        upper_bounds.groupby(["off_pct", "K"], as_index=False)
        .agg(trace_cap_upper_bound_db=("trace_cap_upper_bound_db", "mean"))
    )


def plot_speed_accuracy_pareto(runs, algorithms, out_path, upper_bounds=None):
    selected = [name for name, _ in algorithms]
    summary = (
        runs[runs["algorithm"].isin(selected)]
        .groupby(["off_pct", "K", "algorithm"], as_index=False)
        .agg(
            u_g_db_mean=("u_g_db", "mean"),
            elapsed_p50=("elapsed_seconds", "median"),
        )
    )
    off_pcts = sorted(summary["off_pct"].unique())
    fig, axes = plt.subplots(
        1,
        len(off_pcts),
        figsize=(7.2 * len(off_pcts), 5.0),
        squeeze=False,
    )
    color_map = plt.get_cmap("tab20")
    colors = {
        name: ALGORITHM_COLORS.get(name, color_map(index % color_map.N))
        for index, name in enumerate(selected)
    }

    for col, off_pct in enumerate(off_pcts):
        ax = axes[0, col]
        data = summary[summary["off_pct"] == off_pct].copy()
        if data.empty:
            ax.set_visible(False)
            continue

        K = int(data["K"].iloc[0])
        upper_bound_value = None
        if upper_bounds is not None:
            bound_row = upper_bounds[
                (upper_bounds["off_pct"] == off_pct)
                & (upper_bounds["K"] == K)
            ]
            if not bound_row.empty:
                upper_bound_value = float(bound_row["trace_cap_upper_bound_db"].iloc[0])

        front = pareto_front(data)
        front_names = {row["algorithm"] for row in front}
        for _, row in data.iterrows():
            name = row["algorithm"]
            on_front = name in front_names
            ax.scatter(
                row["elapsed_p50"],
                row["u_g_db_mean"],
                s=62 if on_front else 42,
                color=colors[name],
                edgecolor="black",
                linewidth=0.8 if on_front else 0.4,
                zorder=4 if on_front else 3,
                label=name if col == 0 else None,
            )
        if front:
            front_df = (
                data[data["algorithm"].isin([row["algorithm"] for row in front])]
                .sort_values("elapsed_p50")
            )
            ax.plot(
                front_df["elapsed_p50"],
                front_df["u_g_db_mean"],
                color="#222222",
                linewidth=1.5,
                marker="o",
                markersize=3.0,
                zorder=2,
            )
            for _, row in front_df.iterrows():
                ax.annotate(
                    row["algorithm"],
                    (row["elapsed_p50"], row["u_g_db_mean"]),
                    xytext=(5, 5),
                    textcoords="offset points",
                    fontsize=8,
                    weight="bold",
                )

        if upper_bound_value is not None and np.isfinite(upper_bound_value):
            ax.axhline(
                upper_bound_value,
                color="#111111",
                linestyle=(0, (4, 3)),
                linewidth=1.1,
                alpha=0.75,
                label="trace/cap upper bound" if col == 0 else None,
                zorder=1,
            )

        ax.set_xscale("log")
        ax.set_title(f"{off_pct:g}% off, K={K}")
        ax.set_xlabel("Median runtime, seconds (log scale)")
        if col == 0:
            ax.set_ylabel("Mean 10 lg(U_G), dB")
        ax.grid(True, alpha=0.25)

    fig.suptitle("Speed-accuracy Pareto front")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            loc="center left",
            bbox_to_anchor=(0.995, 0.5),
            fontsize=7,
        )
    fig.tight_layout(rect=(0.0, 0.0, 0.86, 0.93))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_algorithm_comparison_plots(runs, algorithms, out_dir, focused_names):
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
    plot_speed_accuracy_pareto(
        runs,
        algorithms,
        out_dir / "pareto_speed_accuracy.png",
        upper_bounds=_load_trace_upper_bounds(out_dir),
    )
    selected = {name for name, _ in algorithms}
    if set(focused_names).issubset(selected):
        focused_algorithms = tuple(
            (name, solver)
            for name, solver in algorithms
            if name in focused_names
        )
        plot_cdf(
            runs,
            focused_algorithms,
            "u_g_db",
            "10 lg(U_G), dB",
            "Focused cumulative distribution: H3 vs cap-aware candidates",
            out_dir / "cdf_u_g_db_h3_submodular_gen.png",
        )
        plot_cdf(
            runs,
            focused_algorithms,
            "elapsed_seconds",
            "Elapsed time, seconds",
            "Focused solver runtime: H3 vs cap-aware candidates",
            out_dir / "cdf_runtime_seconds_h3_submodular_gen.png",
            log_y=True,
        )
