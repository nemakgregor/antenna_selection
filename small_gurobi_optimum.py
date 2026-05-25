import argparse
import itertools
import math
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

from algorithms import (
    calculate_objectives,
    check_constraints,
    solve_coutino_greedy,
    solve_h1,
    solve_h2,
    solve_miso_energy_greedy,
    solve_pareto_interference_greedy,
)
from motor_challenge_1205 import generate_V


HEURISTICS = (
    ("H1", lambda V, K, sigma, P: solve_h1(V, K, sigma=sigma, P=P)),
    ("H2", lambda V, K, sigma, P: solve_h2(V, K, sigma=sigma, P=P)),
    ("Coutino", lambda V, K, sigma, P: solve_coutino_greedy(V, K, sigma=sigma, P=P)),
    (
        "MISO-EE",
        lambda V, K, sigma, P: solve_miso_energy_greedy(
            V, K, sigma=sigma, P=P, target_margin=0.05
        ),
    ),
    (
        "Pareto-H2",
        lambda V, K, sigma, P: solve_pareto_interference_greedy(
            V, K, sigma=sigma, P=P
        ),
    ),
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Small exact Gurobi benchmark for the antenna-selection task. "
            "It enumerates all feasible subsets, then lets gurobipy select "
            "the optimal subset through a binary linear model."
        )
    )
    parser.add_argument("--N", type=int, default=10, help="Small number of antennas.")
    parser.add_argument("--L", type=int, default=2, help="Number of layers.")
    parser.add_argument(
        "--active-frac",
        type=float,
        default=0.5,
        help="K/N active upper bound. Default 0.5 means at least 50%% off.",
    )
    parser.add_argument(
        "--min-active",
        type=int,
        default=None,
        help="Minimum active antennas during exact search. Default: L.",
    )
    parser.add_argument(
        "--exact-k",
        action="store_true",
        help="Optimize among subsets with exactly K active antennas instead of <=K.",
    )
    parser.add_argument(
        "--objective",
        choices=[
            "general",
            "bf",
            "interference",
            "bf_protected_interference",
            "bf_under_h2_interference",
            "balanced",
        ],
        default="bf_protected_interference",
        help="Exact objective to optimize.",
    )
    parser.add_argument(
        "--bf-floor",
        type=float,
        default=0.85,
        help=(
            "For bf_protected_interference: require BF gain >= bf-floor times "
            "the best feasible BF gain."
        ),
    )
    parser.add_argument(
        "--interference-weight",
        type=float,
        default=1.0,
        help="For balanced: score = normalized BF - weight * normalized interference.",
    )
    parser.add_argument(
        "--interference-ceiling-factor",
        type=float,
        default=2.0,
        help=(
            "For bf_under_h2_interference: require interference <= this factor "
            "times H2 interference."
        ),
    )
    parser.add_argument(
        "--max-candidates",
        type=int,
        default=200000,
        help="Safety limit for enumerated subsets.",
    )
    parser.add_argument("--sigma", type=float, default=1.0)
    parser.add_argument("--P", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--out-dir", type=Path, default=Path("results/gurobi_small")
    )
    return parser.parse_args()


def subset_to_x(N, subset):
    x = np.zeros(N, dtype=int)
    x[list(subset)] = 1
    return x


def count_candidates(N, sizes):
    return sum(math.comb(N, size) for size in sizes)


def enumerate_subsets(V, K, min_active, sigma, P, exact_k, max_candidates):
    N = V.shape[0]
    sizes = [K] if exact_k else list(range(min_active, K + 1))
    candidate_count = count_candidates(N, sizes)
    if candidate_count > max_candidates:
        raise ValueError(
            f"Exact search would enumerate {candidate_count} subsets; "
            f"increase --max-candidates or reduce N/K."
        )

    rows = []
    for size in sizes:
        for subset in itertools.combinations(range(N), size):
            x = subset_to_x(N, subset)
            u_bf, u_i, u_g = calculate_objectives(V, x, sigma=sigma, P=P)
            rows.append(
                {
                    "candidate_id": len(rows),
                    "subset": subset,
                    "active_count": size,
                    "u_bf": u_bf,
                    "u_i": u_i,
                    "u_g": u_g,
                }
            )
    return rows


def solve_exact_with_gurobi(
    rows, objective, bf_floor, interference_weight, interference_ceiling=None
):
    try:
        import gurobipy as gp
        from gurobipy import GRB
    except ImportError as exc:
        raise RuntimeError("gurobipy is not installed in the active Python env.") from exc

    model = gp.Model("antenna_subset_exact")
    model.Params.OutputFlag = 0
    y = model.addVars(len(rows), vtype=GRB.BINARY, name="subset")
    model.addConstr(gp.quicksum(y[i] for i in range(len(rows))) == 1)

    best_bf = max(row["u_bf"] for row in rows)
    max_interference = max(row["u_i"] for row in rows)
    max_interference = max(max_interference, np.finfo(float).eps)

    if objective == "general":
        sense = GRB.MAXIMIZE
        coeffs = [row["u_g"] for row in rows]
    elif objective == "bf":
        sense = GRB.MAXIMIZE
        coeffs = [row["u_bf"] for row in rows]
    elif objective == "interference":
        sense = GRB.MINIMIZE
        coeffs = [row["u_i"] for row in rows]
    elif objective == "bf_protected_interference":
        sense = GRB.MINIMIZE
        coeffs = [row["u_i"] for row in rows]
        model.addConstr(
            gp.quicksum(rows[i]["u_bf"] * y[i] for i in range(len(rows)))
            >= bf_floor * best_bf
        )
    elif objective == "bf_under_h2_interference":
        if interference_ceiling is None:
            raise ValueError("interference_ceiling is required for bf_under_h2_interference")
        sense = GRB.MAXIMIZE
        coeffs = [row["u_bf"] for row in rows]
        model.addConstr(
            gp.quicksum(rows[i]["u_i"] * y[i] for i in range(len(rows)))
            <= interference_ceiling
        )
    else:
        sense = GRB.MAXIMIZE
        coeffs = [
            (row["u_bf"] / best_bf)
            - interference_weight * (row["u_i"] / max_interference)
            for row in rows
        ]

    model.setObjective(
        gp.quicksum(coeffs[i] * y[i] for i in range(len(rows))), sense
    )
    model.optimize()
    if model.Status != GRB.OPTIMAL:
        raise RuntimeError(f"Gurobi did not find an optimum; status={model.Status}")

    selected_id = max(range(len(rows)), key=lambda i: y[i].X)
    selected = dict(rows[selected_id])
    selected["objective_value"] = model.ObjVal
    selected["bf_floor_reference"] = best_bf
    selected["bf_floor_required"] = (
        bf_floor * best_bf if objective == "bf_protected_interference" else np.nan
    )
    selected["interference_ceiling"] = (
        interference_ceiling if objective == "bf_under_h2_interference" else np.nan
    )
    return selected


def target_value(row, objective, best_bf=None, max_interference=None, weight=1.0):
    if objective == "general":
        return row["u_g"]
    if objective == "bf":
        return row["u_bf"]
    if objective == "bf_under_h2_interference":
        return row["u_bf"]
    if objective in {"interference", "bf_protected_interference"}:
        return row["u_i"]
    best_bf = max(best_bf, np.finfo(float).eps)
    max_interference = max(max_interference, np.finfo(float).eps)
    return row["u_bf"] / best_bf - weight * row["u_i"] / max_interference


def direct_exact_optimum(
    rows, objective, bf_floor, interference_weight, interference_ceiling=None
):
    best_bf = max(row["u_bf"] for row in rows)
    max_interference = max(row["u_i"] for row in rows)

    feasible = rows
    if objective == "bf_protected_interference":
        feasible = [row for row in rows if row["u_bf"] >= bf_floor * best_bf]
    elif objective == "bf_under_h2_interference":
        if interference_ceiling is None:
            raise ValueError("interference_ceiling is required for bf_under_h2_interference")
        feasible = [row for row in rows if row["u_i"] <= interference_ceiling]
    if not feasible:
        raise RuntimeError("No feasible subset remains after objective-side constraints.")

    kwargs = {
        "best_bf": best_bf,
        "max_interference": max_interference,
        "weight": interference_weight,
    }
    if objective in {"interference", "bf_protected_interference"}:
        return min(feasible, key=lambda row: target_value(row, objective, **kwargs))
    return max(feasible, key=lambda row: target_value(row, objective, **kwargs))


def objective_gap(candidate, optimum, objective, best_bf, max_interference, weight):
    cand_value = target_value(candidate, objective, best_bf, max_interference, weight)
    opt_value = target_value(optimum, objective, best_bf, max_interference, weight)
    denom = max(abs(opt_value), np.finfo(float).eps)
    if objective in {"interference", "bf_protected_interference"}:
        return 100.0 * (cand_value - opt_value) / denom
    return 100.0 * (opt_value - cand_value) / denom


def compare_heuristics(
    V,
    K,
    sigma,
    P,
    optimum,
    objective,
    bf_floor,
    rows,
    interference_weight,
    interference_ceiling=None,
):
    best_bf = max(row["u_bf"] for row in rows)
    max_interference = max(row["u_i"] for row in rows)
    comparison = [
        {
            "algorithm": "Gurobi optimum",
            "active_count": optimum["active_count"],
            "valid": True,
            "u_bf": optimum["u_bf"],
            "u_i": optimum["u_i"],
            "u_g": optimum["u_g"],
            "meets_bf_floor": True,
            "meets_interference_ceiling": True,
            "objective_gap_pct": 0.0,
            "subset": " ".join(map(str, optimum["subset"])),
            "elapsed_seconds": 0.0,
        }
    ]

    for name, solver in HEURISTICS:
        with np.errstate(all="ignore"):
            started_at = time.perf_counter()
            x = solver(V, K, sigma, P)
            elapsed_seconds = time.perf_counter() - started_at
            valid, active_count = check_constraints(x, K)
            u_bf, u_i, u_g = calculate_objectives(V, x, sigma=sigma, P=P)
        row = {
            "u_bf": u_bf,
            "u_i": u_i,
            "u_g": u_g,
        }
        meets_bf_floor = (
            u_bf >= bf_floor * best_bf
            if objective == "bf_protected_interference"
            else True
        )
        meets_interference_ceiling = (
            u_i <= interference_ceiling
            if objective == "bf_under_h2_interference"
            else True
        )
        meets_extra_constraint = meets_bf_floor and meets_interference_ceiling
        gap = (
            objective_gap(
                row,
                optimum,
                objective,
                best_bf,
                max_interference,
                weight=interference_weight,
            )
            if meets_extra_constraint
            else np.nan
        )
        comparison.append(
            {
                "algorithm": name,
                "active_count": active_count,
                "valid": valid,
                "u_bf": u_bf,
                "u_i": u_i,
                "u_g": u_g,
                "meets_bf_floor": meets_bf_floor,
                "meets_interference_ceiling": meets_interference_ceiling,
                "objective_gap_pct": gap,
                "subset": " ".join(map(str, np.flatnonzero(x))),
                "elapsed_seconds": elapsed_seconds,
            }
        )
    return pd.DataFrame(comparison)


def plot_comparison(comparison, out_dir, args):
    plot_df = comparison.copy()
    plot_df["log10_u_g"] = np.log10(
        np.maximum(plot_df["u_g"], np.finfo(float).tiny)
    )
    colors = {
        "Gurobi optimum": "#111111",
        "H1": "#1f77b4",
        "H2": "#ff7f0e",
        "Coutino": "#2ca02c",
        "MISO-EE": "#4f6d7a",
        "Pareto-H2": "#9467bd",
    }
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    panels = [
        ("u_bf", "BF gain (higher is better)"),
        ("u_i", "Interference (lower is better)"),
        ("log10_u_g", "General objective log10(U_G), higher is better"),
    ]
    x = np.arange(len(plot_df))
    for ax, (column, title) in zip(axes, panels):
        ax.bar(
            x,
            plot_df[column],
            color=[colors.get(name, "#777777") for name in plot_df["algorithm"]],
        )
        ax.set_title(title)
        ax.set_xticks(x, plot_df["algorithm"], rotation=30, ha="right")
        ax.grid(True, axis="y", alpha=0.25)
    fig.suptitle(
        f"Small exact benchmark, N={args.N}, L={args.L}, "
        f"K<={int(round(args.N * args.active_frac))}, objective={args.objective}"
    )
    fig.savefig(out_dir / "gurobi_small_comparison.png", dpi=180)
    plt.close(fig)


def write_report(comparison, optimum, candidate_count, out_dir, args):
    lines = [
        "# Small Gurobi Optimum",
        "",
        f"- N: {args.N}",
        f"- L: {args.L}",
        f"- K max: {int(round(args.N * args.active_frac))}",
        f"- Minimum active antennas: {args.min_active if args.min_active is not None else args.L}",
        f"- Exact K only: {args.exact_k}",
        f"- Objective: {args.objective}",
        f"- BF floor: {args.bf_floor}",
        f"- H2 interference ceiling factor: {args.interference_ceiling_factor}",
        f"- Enumerated subsets: {candidate_count}",
        f"- Direct scan verified: {bool(optimum.get('direct_scan_verified', False))}",
        f"- Optimum subset: {' '.join(map(str, optimum['subset']))}",
        "",
        "| algorithm | active | BF gain | interference | U_G | meets BF floor | meets int ceiling | objective gap |",
        "|:---|---:|---:|---:|---:|:---:|:---:|---:|",
    ]
    for _, row in comparison.iterrows():
        gap = (
            ""
            if pd.isna(row["objective_gap_pct"])
            else f"{row['objective_gap_pct']:.2f}%"
        )
        lines.append(
            "| "
            + " | ".join(
                [
                    row["algorithm"],
                    str(int(row["active_count"])),
                    f"{row['u_bf']:.4f}",
                    f"{row['u_i']:.4f}",
                    f"{row['u_g']:.4e}",
                    str(bool(row["meets_bf_floor"])),
                    str(bool(row["meets_interference_ceiling"])),
                    gap,
                ]
            )
            + " |"
        )
    (out_dir / "gurobi_small_report.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    K = int(round(args.N * args.active_frac))
    min_active = args.L if args.min_active is None else args.min_active
    if not (0 <= min_active <= K <= args.N):
        raise ValueError(f"Require 0 <= min_active <= K <= N, got {min_active}, {K}, {args.N}")

    np.random.seed(args.seed)
    V = generate_V(args.N, args.L)

    print(
        f"Enumerating exact candidates for N={args.N}, L={args.L}, "
        f"K<={K}, min_active={min_active}",
        flush=True,
    )
    rows = enumerate_subsets(
        V,
        K,
        min_active,
        args.sigma,
        args.P,
        args.exact_k,
        args.max_candidates,
    )
    with np.errstate(all="ignore"):
        x_h2 = solve_h2(V, K, sigma=args.sigma, P=args.P)
        _, h2_interference, _ = calculate_objectives(
            V, x_h2, sigma=args.sigma, P=args.P
        )
    interference_ceiling = args.interference_ceiling_factor * h2_interference
    print(f"Solving Gurobi model with {len(rows)} subset variables", flush=True)
    optimum = solve_exact_with_gurobi(
        rows,
        args.objective,
        args.bf_floor,
        args.interference_weight,
        interference_ceiling,
    )
    direct_optimum = direct_exact_optimum(
        rows,
        args.objective,
        args.bf_floor,
        args.interference_weight,
        interference_ceiling,
    )
    best_bf = max(row["u_bf"] for row in rows)
    max_interference = max(row["u_i"] for row in rows)
    optimum_value = target_value(
        optimum,
        args.objective,
        best_bf=best_bf,
        max_interference=max_interference,
        weight=args.interference_weight,
    )
    direct_value = target_value(
        direct_optimum,
        args.objective,
        best_bf=best_bf,
        max_interference=max_interference,
        weight=args.interference_weight,
    )
    optimum["direct_scan_objective_value"] = direct_value
    optimum["direct_scan_verified"] = bool(np.isclose(optimum_value, direct_value))
    if not optimum["direct_scan_verified"]:
        raise RuntimeError("Gurobi optimum did not match direct enumeration scan.")
    comparison = compare_heuristics(
        V,
        K,
        args.sigma,
        args.P,
        optimum,
        args.objective,
        args.bf_floor,
        rows,
        args.interference_weight,
        interference_ceiling,
    )

    comparison.to_csv(args.out_dir / "gurobi_small_comparison.csv", index=False)
    pd.DataFrame([{**optimum, "subset": " ".join(map(str, optimum["subset"]))}]).to_csv(
        args.out_dir / "gurobi_small_optimum.csv", index=False
    )
    plot_comparison(comparison, args.out_dir, args)
    write_report(comparison, optimum, len(rows), args.out_dir, args)

    print("\nBest exact subset:")
    print(
        f"  active={optimum['active_count']} subset={' '.join(map(str, optimum['subset']))}"
    )
    print(
        f"  BF={optimum['u_bf']:.4f}  Int={optimum['u_i']:.4f}  U_G={optimum['u_g']:.4e}"
    )
    print("\nSaved:")
    for path in [
        "gurobi_small_optimum.csv",
        "gurobi_small_comparison.csv",
        "gurobi_small_report.md",
        "gurobi_small_comparison.png",
    ]:
        print(f"  {args.out_dir / path}")


if __name__ == "__main__":
    main()
