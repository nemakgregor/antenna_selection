from utils.plotting import use_agg_backend

use_agg_backend()
import matplotlib.pyplot as plt
import numpy as np


ALGORITHM_COLORS = {
    "H1": "#0072B2",
    "H2": "#D55E00",
    "H3": "#009E73",
    "FrameOnly-Gen": "#6A3D9A",
    "CapWindow-Gen": "#1F78B4",
    "Frame-Gen": "#CC79A7",
    "S-threshold-Gen": "#E69F00",
    "Coutino": "#000000",
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
            "Focused cumulative distribution: H3 vs FrameOnly-Gen vs CapWindow-Gen",
            out_dir / "cdf_u_g_db_h3_frameonly_capwindow.png",
        )
        plot_cdf(
            runs,
            focused_algorithms,
            "elapsed_seconds",
            "Elapsed time, seconds",
            "Focused solver runtime: H3 vs FrameOnly-Gen vs CapWindow-Gen",
            out_dir / "cdf_runtime_seconds_h3_frameonly_capwindow.png",
            log_y=True,
        )


def write_threshold_exploration_plots(runs, out_dir, top_n=8):
    plot_threshold_cdf_grid(
        runs,
        "u_g",
        "U_G",
        "Threshold exploration: raw general objective",
        out_dir / "threshold_cdf_u_g.png",
        top_n=top_n,
    )
    plot_threshold_cdf_grid(
        runs,
        "u_g_db",
        "10 lg(U_G), dB",
        "Threshold exploration: general objective in dB",
        out_dir / "threshold_cdf_u_g_db.png",
        top_n=top_n,
    )
    plot_threshold_cdf_grid(
        runs,
        "u_g_vs_best_T",
        "U_G(T) / best_T U_G",
        "Threshold exploration: normalized objective gap",
        out_dir / "threshold_cdf_u_g_vs_best.png",
        top_n=top_n,
    )


def plot_threshold_cdf_grid(runs, value_col, ylabel, title, out_path, top_n=8):
    profiles = sorted(runs["data_profile"].unique())
    off_pcts = sorted(runs["off_pct"].unique())
    fig, axes = plt.subplots(
        len(profiles),
        len(off_pcts),
        figsize=(6.4 * len(off_pcts), 3.9 * len(profiles)),
        sharey=value_col == "u_g_vs_best_T",
        squeeze=False,
    )
    color_map = plt.get_cmap("tab20")
    handles = []
    labels = []

    for row, profile in enumerate(profiles):
        for col, off_pct in enumerate(off_pcts):
            ax = axes[row, col]
            data = runs[
                (runs["data_profile"] == profile)
                & (runs["off_pct"] == off_pct)
            ]
            if data.empty:
                ax.set_visible(False)
                continue

            K = int(data["K"].iloc[0])
            top_thresholds = (
                data.groupby("T")["u_g_vs_best_T"]
                .mean()
                .sort_values(ascending=False)
                .head(top_n)
                .index.tolist()
            )
            if 0 in set(data["T"]):
                top_thresholds = [0, *[T for T in top_thresholds if T != 0]]
                top_thresholds = top_thresholds[:top_n]

            for index, T in enumerate(top_thresholds):
                values = data[data["T"] == T][value_col]
                x_values, y_values = empirical_cdf(values)
                if len(x_values) == 0:
                    continue
                label = f"T={int(T)}"
                line = ax.step(
                    y_values,
                    x_values,
                    where="post",
                    linewidth=1.6,
                    color=color_map(index % color_map.N),
                    label=label,
                )[0]
                if row == 0 and col == 0:
                    handles.append(line)
                    labels.append(label)

            ax.set_title(f"{profile}, {off_pct:g}% off, K={K}")
            ax.set_xlabel("Cumulative fraction of examples")
            if col == 0:
                ax.set_ylabel(ylabel)
            ax.set_xlim(0.0, 1.02)
            ax.grid(True, alpha=0.25)

    fig.suptitle(title)
    if handles:
        fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.995, 0.5), fontsize=8)
    fig.tight_layout(rect=(0.0, 0.0, 0.84, 0.94))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
