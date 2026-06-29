from utils.plotting import use_agg_backend

use_agg_backend()
import matplotlib.pyplot as plt
import numpy as np


SIGMA_COLORS = {
    "H1": "#1f77b4",
    "H2": "#ff7f0e",
    "Coutino": "#2ca02c",
    "BackwardTrueGreedy": "#2ca02c",
    "CoutinoSchur-Gen": "#7a3e9d",
    "MISO-EE": "#4f6d7a",
    "Pareto-H2": "#9467bd",
    "H3-threshold-BF": "#8c564b",
    "H3-threshold-Int": "#e377c2",
    "H3-threshold-Gen": "#7f7f7f",
    "Frame-BF": "#17a398",
    "Frame-Int": "#b56576",
    "Frame-Gen": "#3366cc",
    "CapWindow-Gen": "#1f78b4",
    "CapSubmod-Gen": "#007a87",
    "ThreshDOpt-Gen": "#0099c6",
    "ThreshWLogdet-Gen": "#6c8e00",
    "ThreshDOptSwap-Gen": "#d55e00",
    "H3-Fast": "#bcbd22",
}

SIGMA_MARKERS = {
    "H1": "o",
    "H2": "s",
    "Coutino": "^",
    "BackwardTrueGreedy": "^",
    "CoutinoSchur-Gen": "P",
    "MISO-EE": "D",
    "Pareto-H2": "P",
    "H3-threshold-BF": "v",
    "H3-threshold-Int": "X",
    "H3-threshold-Gen": "*",
    "Frame-BF": ">",
    "Frame-Int": "<",
    "Frame-Gen": "8",
    "CapWindow-Gen": "d",
    "CapSubmod-Gen": "H",
    "ThreshDOpt-Gen": "*",
    "ThreshWLogdet-Gen": "p",
    "ThreshDOptSwap-Gen": "X",
    "H3-Fast": "h",
}


def plot_sigma_metrics(summary, out_dir, args, algorithms, K, data_k, sigma_values):
    k_label = data_k["K_label"].iloc[0]
    sigma_to_x = {sigma: idx for idx, sigma in enumerate(sigma_values)}
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), constrained_layout=True)
    panels = [
        ("u_bf_mean", "BF gain, higher is better"),
        ("u_i_mean", "Interference, lower is better"),
        ("u_g_mean", "U_G = det(V_eq V_eq* + sigma I), higher is better"),
    ]

    for ax, (column, title) in zip(axes.flat, panels):
        for heuristic, _ in algorithms:
            curve = data_k[data_k["heuristic"] == heuristic].sort_values("sigma")
            ax.plot(
                curve["sigma"].map(sigma_to_x),
                curve[column],
                marker=SIGMA_MARKERS.get(heuristic, "o"),
                linewidth=1.8,
                markersize=4,
                color=SIGMA_COLORS.get(heuristic),
                label=heuristic,
            )
        ax.set_title(title)
        ax.set_xlabel("sigma")
        ax.set_xticks(range(len(sigma_values)), [f"{sigma:g}" for sigma in sigma_values])
        ax.tick_params(axis="x", labelrotation=35)
        ax.grid(True, alpha=0.25)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))

    axes[0].legend(loc="best", fontsize=8)
    fig.suptitle(
        f"Raw objective values, N={args.N}, L={args.L}, {k_label}, samples={args.samples}"
    )
    fig.savefig(out_dir / f"sigma_sweep_K{K}.png", dpi=180)
    plt.close(fig)


def plot_sigma_winners(winners, out_dir, args, algorithms, K, data_k, sigma_values):
    k_label = data_k["K_label"].iloc[0]
    wins_k = winners[(winners["K"] == K) & (winners["metric"] == "u_g")]
    fig, ax = plt.subplots(figsize=(14, 4.8), constrained_layout=True)
    bottom = np.zeros(len(sigma_values), dtype=float)
    for heuristic, _ in algorithms:
        shares = []
        for sigma in sigma_values:
            row = wins_k[
                (wins_k["sigma"] == sigma) & (wins_k["heuristic"] == heuristic)
            ]
            shares.append(float(row["win_fraction"].iloc[0]) if len(row) else 0.0)
        ax.bar(
            range(len(sigma_values)),
            shares,
            bottom=bottom,
            color=SIGMA_COLORS.get(heuristic),
            label=heuristic,
        )
        bottom += np.asarray(shares)
    ax.set_title(f"Per-sample U_G winner share, N={args.N}, L={args.L}, {k_label}")
    ax.set_xlabel("sigma")
    ax.set_ylabel("split win share")
    ax.set_xticks(range(len(sigma_values)), [f"{sigma:g}" for sigma in sigma_values])
    ax.tick_params(axis="x", labelrotation=35)
    ax.set_ylim(0.0, 1.0)
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8)
    fig.savefig(out_dir / f"sigma_winners_K{K}.png", dpi=180)
    plt.close(fig)


def plot_sigma_sweep(summary, winners, out_dir, args, algorithms):
    for K in sorted(summary["K"].unique()):
        data_k = summary[summary["K"] == K]
        sigma_values = sorted(data_k["sigma"].unique())
        plot_sigma_metrics(summary, out_dir, args, algorithms, K, data_k, sigma_values)
        plot_sigma_winners(winners, out_dir, args, algorithms, K, data_k, sigma_values)
