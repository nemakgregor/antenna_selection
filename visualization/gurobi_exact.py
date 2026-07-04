from utils.plotting import use_agg_backend

use_agg_backend()
import matplotlib.pyplot as plt
import numpy as np


GUROBI_COLORS = {
    "Gurobi optimum": "#111111",
    "Gurobi-BF": "#111111",
    "Gurobi-Int": "#d62728",
    "Gurobi-Gen": "#17becf",
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


def gurobi_algorithm_name(target_obj):
    return {
        "bf": "Gurobi-BF",
        "int": "Gurobi-Int",
        "gen": "Gurobi-Gen",
    }[target_obj]


def plot_gurobi_comparison(comparison, out_dir, args):
    plot_df = comparison.copy()
    plot_df["log10_u_g"] = np.log10(
        np.maximum(plot_df["u_g"], np.finfo(float).tiny)
    )
    fig, axes = plt.subplots(2, 2, figsize=(15, 9), constrained_layout=True)
    axes = axes.ravel()
    panels = [
        ("u_bf", "BF gain (higher is better)"),
        ("u_i", "Interference (lower is better)"),
        ("log10_u_g", "General objective log10(U_G), higher is better"),
        ("elapsed_seconds", "Runtime seconds, lower is better"),
    ]
    x = np.arange(len(plot_df))
    for ax, (column, title) in zip(axes, panels):
        values = plot_df[column]
        if column == "elapsed_seconds":
            values = np.maximum(values, np.finfo(float).tiny)
        ax.bar(
            x,
            values,
            color=[GUROBI_COLORS.get(name, "#777777") for name in plot_df["algorithm"]],
        )
        ax.set_title(title)
        ax.set_xticks(x, plot_df["algorithm"], rotation=30, ha="right")
        ax.grid(True, axis="y", alpha=0.25)
        if column == "elapsed_seconds":
            ax.set_yscale("log")
    fig.suptitle(
        f"Small exact benchmark, N={args.N}, L={args.L}, "
        f"K<={int(round(args.N * args.active_frac))}, objective={args.objective}"
    )
    fig.savefig(out_dir / "gurobi_small_comparison.png", dpi=180)
    plt.close(fig)


def plot_multi_objective_summary(summary, objectives, objective_specs, out_dir):
    fig, axes = plt.subplots(
        2,
        len(objectives),
        figsize=(6.2 * len(objectives), 9.0),
        constrained_layout=True,
    )
    if len(objectives) == 1:
        axes = np.asarray(axes).reshape(2, 1)

    for col, target_obj in enumerate(objectives):
        spec = objective_specs[target_obj]
        data = summary[summary["target_obj"] == target_obj].copy()
        target_gurobi = gurobi_algorithm_name(target_obj)
        data["gurobi_rank"] = np.select(
            [
                data["algorithm"].eq(target_gurobi),
                data["algorithm"].str.startswith("Gurobi-"),
            ],
            [0, 1],
            default=2,
        )
        if spec["direction"] == "min":
            data = data.sort_values(
                ["gurobi_rank", "target_value_mean", "elapsed_seconds_mean"]
            )
        else:
            data = data.sort_values(
                ["gurobi_rank", "target_value_mean", "elapsed_seconds_mean"],
                ascending=[True, False, True],
            )
        algorithms = data["algorithm"].tolist()
        x = np.arange(len(algorithms))
        bar_colors = [GUROBI_COLORS.get(name, "#777777") for name in algorithms]

        axes[0, col].bar(x, data["target_value_mean"], color=bar_colors)
        axes[0, col].set_title(
            f"{spec['label']} mean ({'lower' if spec['direction'] == 'min' else 'higher'} is better)"
        )
        axes[0, col].set_xticks(x, algorithms, rotation=35, ha="right")
        axes[0, col].grid(True, axis="y", alpha=0.25)
        if spec["metric"] in {"u_i", "u_g"}:
            axes[0, col].set_yscale("log")

        runtime = np.maximum(data["elapsed_seconds_mean"], np.finfo(float).tiny)
        axes[1, col].bar(x, runtime, color=bar_colors)
        axes[1, col].set_title(f"Runtime mean for objective={target_obj}")
        axes[1, col].set_xticks(x, algorithms, rotation=35, ha="right")
        axes[1, col].set_yscale("log")
        axes[1, col].grid(True, axis="y", alpha=0.25)

    fig.savefig(out_dir / "gurobi_multi_objective_summary.png", dpi=180)
    plt.close(fig)
