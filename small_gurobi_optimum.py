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
    solve_frame_portfolio,
    solve_h1,
    solve_h2,
    solve_h3,
    solve_h3_fast,
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
    (
        "H3-threshold-BF",
        lambda V, K, sigma, P: solve_h3(
            V, K, target_obj="bf", sigma=sigma, P=P
        ),
    ),
    (
        "H3-threshold-Int",
        lambda V, K, sigma, P: solve_h3(
            V, K, target_obj="int", sigma=sigma, P=P
        ),
    ),
    (
        "H3-threshold-Gen",
        lambda V, K, sigma, P: solve_h3(
            V, K, target_obj="gen", sigma=sigma, P=P
        ),
    ),
    (
        "Frame-BF",
        lambda V, K, sigma, P: solve_frame_portfolio(
            V,
            K,
            target_obj="bf",
            sigma=sigma,
            P=P,
            random_state=0,
            max_refined_starts=12,
            max_passes=20,
            remove_limit=None,
            add_limit=None,
            random_restarts=20,
        ),
    ),
    (
        "Frame-Int",
        lambda V, K, sigma, P: solve_frame_portfolio(
            V,
            K,
            target_obj="int",
            sigma=sigma,
            P=P,
            random_state=0,
            max_refined_starts=12,
            max_passes=20,
            remove_limit=None,
            add_limit=None,
            random_restarts=20,
        ),
    ),
    (
        "Frame-Gen",
        lambda V, K, sigma, P: solve_frame_portfolio(
            V,
            K,
            target_obj="gen",
            sigma=sigma,
            P=P,
            random_state=0,
            max_refined_starts=12,
            max_passes=20,
            remove_limit=None,
            add_limit=None,
            random_restarts=20,
        ),
    ),
    ("H3-Fast", lambda V, K, sigma, P: solve_h3_fast(V, K, random_state=0)),
)

STANDARD_OBJECTIVES = {
    "bf": {
        "exact_objective": "bf",
        "metric": "u_bf",
        "direction": "max",
        "label": "BF gain",
    },
    "int": {
        "exact_objective": "interference",
        "metric": "u_i",
        "direction": "min",
        "label": "Interference",
    },
    "gen": {
        "exact_objective": "general",
        "metric": "u_g",
        "direction": "max",
        "label": "General objective",
    },
}


def gurobi_algorithm_name(target_obj):
    return {
        "bf": "Gurobi-BF",
        "int": "Gurobi-Int",
        "gen": "Gurobi-Gen",
    }[target_obj]


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
        "--samples",
        type=int,
        default=1,
        help=(
            "Number of random matrices to average. If greater than 1 and "
            "--objectives is omitted, the script runs bf/int/gen objectives."
        ),
    )
    parser.add_argument(
        "--objectives",
        nargs="+",
        choices=sorted(STANDARD_OBJECTIVES),
        default=None,
        help=(
            "Run several exact Gurobi objectives. Use bf/int/gen to mirror "
            "h3_threshold.py."
        ),
    )
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
    exact_timing=None,
    gurobi_solutions=None,
):
    exact_timing = exact_timing or {}
    best_bf = max(row["u_bf"] for row in rows)
    max_interference = max(row["u_i"] for row in rows)
    comparison = []

    if gurobi_solutions is None:
        gurobi_solutions = [
            {
                "algorithm": "Gurobi optimum",
                "optimum": optimum,
                "timing": exact_timing,
            }
        ]

    for solution in gurobi_solutions:
        solution_optimum = solution["optimum"]
        solution_timing = solution.get("timing", {})
        solution_row = {
            "u_bf": solution_optimum["u_bf"],
            "u_i": solution_optimum["u_i"],
            "u_g": solution_optimum["u_g"],
        }
        meets_bf_floor = (
            solution_row["u_bf"] >= bf_floor * best_bf
            if objective == "bf_protected_interference"
            else True
        )
        meets_interference_ceiling = (
            solution_row["u_i"] <= interference_ceiling
            if objective == "bf_under_h2_interference"
            else True
        )
        meets_extra_constraint = meets_bf_floor and meets_interference_ceiling
        comparison.append(
            {
                "algorithm": solution["algorithm"],
                "active_count": solution_optimum["active_count"],
                "valid": True,
                "u_bf": solution_optimum["u_bf"],
                "u_i": solution_optimum["u_i"],
                "u_g": solution_optimum["u_g"],
                "meets_bf_floor": meets_bf_floor,
                "meets_interference_ceiling": meets_interference_ceiling,
                "objective_gap_pct": (
                    objective_gap(
                        solution_row,
                        optimum,
                        objective,
                        best_bf,
                        max_interference,
                        weight=interference_weight,
                    )
                    if meets_extra_constraint
                    else np.nan
                ),
                "subset": " ".join(map(str, solution_optimum["subset"])),
                "elapsed_seconds": solution_timing.get("total_seconds", np.nan),
                "enumeration_seconds": solution_timing.get("enumeration_seconds", np.nan),
                "gurobi_model_seconds": solution_timing.get("gurobi_model_seconds", np.nan),
                "direct_scan_seconds": solution_timing.get("direct_scan_seconds", np.nan),
            }
        )

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
                "enumeration_seconds": np.nan,
                "gurobi_model_seconds": np.nan,
                "direct_scan_seconds": np.nan,
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
        "Gurobi-BF": "#111111",
        "Gurobi-Int": "#d62728",
        "Gurobi-Gen": "#17becf",
        "H1": "#1f77b4",
        "H2": "#ff7f0e",
        "Coutino": "#2ca02c",
        "MISO-EE": "#4f6d7a",
        "Pareto-H2": "#9467bd",
        "H3-threshold-BF": "#8c564b",
        "H3-threshold-Int": "#e377c2",
        "H3-threshold-Gen": "#7f7f7f",
        "Frame-BF": "#17a398",
        "Frame-Int": "#b56576",
        "Frame-Gen": "#3366cc",
        "H3-Fast": "#bcbd22",
    }
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
            color=[colors.get(name, "#777777") for name in plot_df["algorithm"]],
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
        f"- Algorithms: {', '.join(comparison['algorithm'].tolist())}",
        "",
        "## Timing",
        "",
        f"- Enumeration: {comparison.loc[0, 'enumeration_seconds']:.6f}s",
        f"- Gurobi model build/solve: {comparison.loc[0, 'gurobi_model_seconds']:.6f}s",
        f"- Direct scan verification: {comparison.loc[0, 'direct_scan_seconds']:.6f}s",
        f"- Gurobi total compared below: {comparison.loc[0, 'elapsed_seconds']:.6f}s",
        "",
        "| algorithm | active | BF gain | interference | U_G | time, s | meets BF floor | meets int ceiling | objective gap |",
        "|:---|---:|---:|---:|---:|---:|:---:|:---:|---:|",
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
                    f"{row['elapsed_seconds']:.6f}",
                    str(bool(row["meets_bf_floor"])),
                    str(bool(row["meets_interference_ceiling"])),
                    gap,
                ]
            )
            + " |"
        )
    (out_dir / "gurobi_small_report.md").write_text("\n".join(lines), encoding="utf-8")


def verify_exact_optimum(rows, optimum, objective, bf_floor, interference_weight, interference_ceiling):
    direct_started_at = time.perf_counter()
    direct_optimum = direct_exact_optimum(
        rows,
        objective,
        bf_floor,
        interference_weight,
        interference_ceiling,
    )
    direct_scan_seconds = time.perf_counter() - direct_started_at

    best_bf = max(row["u_bf"] for row in rows)
    max_interference = max(row["u_i"] for row in rows)
    optimum_value = target_value(
        optimum,
        objective,
        best_bf=best_bf,
        max_interference=max_interference,
        weight=interference_weight,
    )
    direct_value = target_value(
        direct_optimum,
        objective,
        best_bf=best_bf,
        max_interference=max_interference,
        weight=interference_weight,
    )
    optimum["direct_scan_objective_value"] = direct_value
    optimum["direct_scan_verified"] = bool(np.isclose(optimum_value, direct_value))
    if not optimum["direct_scan_verified"]:
        raise RuntimeError("Gurobi optimum did not match direct enumeration scan.")
    return direct_scan_seconds


def run_multi_objective_benchmark(args):
    if args.samples <= 0:
        raise ValueError("--samples must be positive.")

    objectives = list(dict.fromkeys(args.objectives or ["bf", "int", "gen"]))
    K = int(round(args.N * args.active_frac))
    min_active = args.L if args.min_active is None else args.min_active
    if not (0 <= min_active <= K <= args.N):
        raise ValueError(f"Require 0 <= min_active <= K <= N, got {min_active}, {K}, {args.N}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    run_frames = []
    optima = []
    total_cases = args.samples * len(objectives)

    for sample_idx in range(args.samples):
        seed = args.seed + sample_idx
        np.random.seed(seed)
        V = generate_V(args.N, args.L)

        print(
            f"[sample {sample_idx + 1}/{args.samples}] Enumerating N={args.N}, "
            f"L={args.L}, K<={K}, seed={seed}",
            flush=True,
        )
        enumeration_started_at = time.perf_counter()
        rows = enumerate_subsets(
            V,
            K,
            min_active,
            args.sigma,
            args.P,
            args.exact_k,
            args.max_candidates,
        )
        enumeration_seconds = time.perf_counter() - enumeration_started_at

        with np.errstate(all="ignore"):
            x_h2 = solve_h2(V, K, sigma=args.sigma, P=args.P)
            _, h2_interference, _ = calculate_objectives(
                V, x_h2, sigma=args.sigma, P=args.P
            )
        interference_ceiling = args.interference_ceiling_factor * h2_interference
        sample_gurobi_solutions = []

        for objective_pos, target_obj in enumerate(objectives, start=1):
            case_no = sample_idx * len(objectives) + objective_pos
            spec = STANDARD_OBJECTIVES[target_obj]
            objective = spec["exact_objective"]
            print(
                f"  [{case_no}/{total_cases}] Solving Gurobi objective={target_obj} "
                f"with {len(rows)} subset variables",
                flush=True,
            )
            gurobi_started_at = time.perf_counter()
            optimum = solve_exact_with_gurobi(
                rows,
                objective,
                args.bf_floor,
                args.interference_weight,
                interference_ceiling,
            )
            gurobi_model_seconds = time.perf_counter() - gurobi_started_at
            direct_scan_seconds = verify_exact_optimum(
                rows,
                optimum,
                objective,
                args.bf_floor,
                args.interference_weight,
                interference_ceiling,
            )
            exact_timing = {
                "enumeration_seconds": enumeration_seconds,
                "gurobi_model_seconds": gurobi_model_seconds,
                "direct_scan_seconds": direct_scan_seconds,
                "total_seconds": enumeration_seconds + gurobi_model_seconds,
            }
            sample_gurobi_solutions.append(
                {
                    "target_obj": target_obj,
                    "algorithm": gurobi_algorithm_name(target_obj),
                    "optimum": optimum,
                    "timing": exact_timing,
                }
            )

            optima.append(
                {
                    "target_obj": target_obj,
                    "algorithm": gurobi_algorithm_name(target_obj),
                    "exact_objective": objective,
                    "sample": sample_idx,
                    "seed": seed,
                    "N": args.N,
                    "L": args.L,
                    "K": K,
                    "candidate_count": len(rows),
                    "active_count": optimum["active_count"],
                    "u_bf": optimum["u_bf"],
                    "u_i": optimum["u_i"],
                    "u_g": optimum["u_g"],
                    "objective_value": optimum["objective_value"],
                    "direct_scan_verified": optimum["direct_scan_verified"],
                    "subset": " ".join(map(str, optimum["subset"])),
                    **exact_timing,
                }
            )

        for target_obj in objectives:
            spec = STANDARD_OBJECTIVES[target_obj]
            objective = spec["exact_objective"]
            target_solution = next(
                solution
                for solution in sample_gurobi_solutions
                if solution["target_obj"] == target_obj
            )
            comparison = compare_heuristics(
                V,
                K,
                args.sigma,
                args.P,
                target_solution["optimum"],
                objective,
                args.bf_floor,
                rows,
                args.interference_weight,
                interference_ceiling,
                exact_timing=target_solution["timing"],
                gurobi_solutions=sample_gurobi_solutions,
            )
            comparison.insert(0, "target_obj", target_obj)
            comparison.insert(1, "exact_objective", objective)
            comparison.insert(2, "sample", sample_idx)
            comparison.insert(3, "seed", seed)
            comparison.insert(4, "N", args.N)
            comparison.insert(5, "L", args.L)
            comparison.insert(6, "K", K)
            comparison["target_metric"] = spec["metric"]
            comparison["target_direction"] = spec["direction"]
            comparison["target_value"] = comparison[spec["metric"]]
            comparison["candidate_count"] = len(rows)
            run_frames.append(comparison)

    runs = pd.concat(run_frames, ignore_index=True)
    summary = build_multi_summary(runs)
    wins = build_multi_wins(runs)
    optima_df = pd.DataFrame(optima)

    runs.to_csv(args.out_dir / "gurobi_multi_objective_runs.csv", index=False)
    summary.to_csv(args.out_dir / "gurobi_multi_objective_summary.csv", index=False)
    wins.to_csv(args.out_dir / "gurobi_multi_objective_wins.csv", index=False)
    optima_df.to_csv(args.out_dir / "gurobi_multi_objective_optima.csv", index=False)
    plot_multi_objective_summary(summary, objectives, args.out_dir)
    write_multi_objective_report(summary, wins, objectives, args.out_dir, args)

    print(f"wrote multi-objective results to {args.out_dir}")


def build_multi_summary(runs):
    group_cols = ["target_obj", "exact_objective", "target_metric", "target_direction", "algorithm"]
    return (
        runs.groupby(group_cols, as_index=False)
        .agg(
            active_count_mean=("active_count", "mean"),
            active_count_std=("active_count", "std"),
            u_bf_mean=("u_bf", "mean"),
            u_bf_std=("u_bf", "std"),
            u_i_mean=("u_i", "mean"),
            u_i_std=("u_i", "std"),
            u_g_mean=("u_g", "mean"),
            u_g_std=("u_g", "std"),
            target_value_mean=("target_value", "mean"),
            target_value_std=("target_value", "std"),
            objective_gap_pct_mean=("objective_gap_pct", "mean"),
            objective_gap_pct_std=("objective_gap_pct", "std"),
            elapsed_seconds_mean=("elapsed_seconds", "mean"),
            elapsed_seconds_std=("elapsed_seconds", "std"),
            enumeration_seconds_mean=("enumeration_seconds", "mean"),
            gurobi_model_seconds_mean=("gurobi_model_seconds", "mean"),
            direct_scan_seconds_mean=("direct_scan_seconds", "mean"),
            samples=("sample", "count"),
        )
        .fillna(0.0)
    )


def build_multi_wins(runs):
    win_rows = []
    for keys, chunk in runs.groupby(["target_obj", "sample"]):
        target_obj, sample = keys
        spec = STANDARD_OBJECTIVES[target_obj]
        values = chunk.set_index("algorithm")[spec["metric"]]
        if spec["direction"] == "min":
            best_value = values.min()
        else:
            best_value = values.max()
        winners = values[np.isclose(values, best_value)].index.tolist()
        share = 1.0 / len(winners)
        for winner in winners:
            win_rows.append(
                {
                    "target_obj": target_obj,
                    "sample": sample,
                    "algorithm": winner,
                    "win_share": share,
                    "winner_hit": 1.0,
                }
            )
    wins = pd.DataFrame(win_rows)
    return (
        wins.groupby(["target_obj", "algorithm"], as_index=False)
        .agg(win_share=("win_share", "sum"), winner_hits=("winner_hit", "sum"))
        .merge(
            runs.groupby("target_obj", as_index=False)["sample"]
            .nunique()
            .rename(columns={"sample": "samples"}),
            on="target_obj",
            how="left",
        )
        .assign(
            win_fraction=lambda df: df["win_share"] / df["samples"],
            winner_rate=lambda df: df["winner_hits"] / df["samples"],
        )
    )


def plot_multi_objective_summary(summary, objectives, out_dir):
    fig, axes = plt.subplots(
        2,
        len(objectives),
        figsize=(6.2 * len(objectives), 9.0),
        constrained_layout=True,
    )
    if len(objectives) == 1:
        axes = np.asarray(axes).reshape(2, 1)

    colors = {
        "Gurobi optimum": "#111111",
        "Gurobi-BF": "#111111",
        "Gurobi-Int": "#d62728",
        "Gurobi-Gen": "#17becf",
        "H1": "#1f77b4",
        "H2": "#ff7f0e",
        "Coutino": "#2ca02c",
        "MISO-EE": "#4f6d7a",
        "Pareto-H2": "#9467bd",
        "H3-threshold-BF": "#8c564b",
        "H3-threshold-Int": "#e377c2",
        "H3-threshold-Gen": "#7f7f7f",
        "Frame-BF": "#17a398",
        "Frame-Int": "#b56576",
        "Frame-Gen": "#3366cc",
        "H3-Fast": "#bcbd22",
    }

    for col, target_obj in enumerate(objectives):
        spec = STANDARD_OBJECTIVES[target_obj]
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
        bar_colors = [colors.get(name, "#777777") for name in algorithms]

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


def format_number(value):
    if pd.isna(value):
        return ""
    if abs(value) >= 1e5 or (0 < abs(value) < 1e-3):
        return f"{value:.4e}"
    return f"{value:.4f}"


def write_multi_objective_report(summary, wins, objectives, out_dir, args):
    lines = [
        "# Multi-Objective Small Gurobi Benchmark",
        "",
        f"- N: {args.N}",
        f"- L: {args.L}",
        f"- K max: {int(round(args.N * args.active_frac))}",
        f"- Minimum active antennas: {args.min_active if args.min_active is not None else args.L}",
        f"- Exact K only: {args.exact_k}",
        f"- Samples: {args.samples}",
        f"- Seed range: {args.seed}..{args.seed + args.samples - 1}",
        f"- Objectives: {', '.join(objectives)}",
        "",
        "`Gurobi-BF`, `Gurobi-Int`, and `Gurobi-Gen` are separate exact solutions evaluated in every objective section.",
        "Runtime for each Gurobi variant includes subset enumeration plus that variant's Gurobi model build/solve; model-only timing is in the CSV.",
        "",
    ]

    for target_obj in objectives:
        spec = STANDARD_OBJECTIVES[target_obj]
        data = summary[summary["target_obj"] == target_obj].copy()
        win_data = wins[wins["target_obj"] == target_obj].set_index("algorithm")
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

        lines.extend(
            [
                f"## Objective `{target_obj}`",
                "",
                f"Target metric: `{spec['metric']}` ({'minimize' if spec['direction'] == 'min' else 'maximize'}).",
                "",
                "| algorithm | target mean | gap mean | BF mean | Int mean | U_G mean | time mean, s | winner rate | split win share |",
                "|:---|---:|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for _, row in data.iterrows():
            winner_rate = (
                win_data.loc[row["algorithm"], "winner_rate"]
                if row["algorithm"] in win_data.index
                else 0.0
            )
            win_fraction = (
                win_data.loc[row["algorithm"], "win_fraction"]
                if row["algorithm"] in win_data.index
                else 0.0
            )
            lines.append(
                "| "
                + " | ".join(
                    [
                        row["algorithm"],
                        format_number(row["target_value_mean"]),
                        f"{row['objective_gap_pct_mean']:.2f}%",
                        format_number(row["u_bf_mean"]),
                        format_number(row["u_i_mean"]),
                        format_number(row["u_g_mean"]),
                        f"{row['elapsed_seconds_mean']:.6f}",
                        f"{winner_rate:.2f}",
                        f"{win_fraction:.2f}",
                    ]
                )
                + " |"
            )
        lines.append("")

    (out_dir / "gurobi_multi_objective_report.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def main():
    args = parse_args()
    if args.objectives is not None or args.samples > 1:
        run_multi_objective_benchmark(args)
        return

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
    enumeration_started_at = time.perf_counter()
    rows = enumerate_subsets(
        V,
        K,
        min_active,
        args.sigma,
        args.P,
        args.exact_k,
        args.max_candidates,
    )
    enumeration_seconds = time.perf_counter() - enumeration_started_at
    with np.errstate(all="ignore"):
        x_h2 = solve_h2(V, K, sigma=args.sigma, P=args.P)
        _, h2_interference, _ = calculate_objectives(
            V, x_h2, sigma=args.sigma, P=args.P
        )
    interference_ceiling = args.interference_ceiling_factor * h2_interference
    print(f"Solving Gurobi model with {len(rows)} subset variables", flush=True)
    gurobi_started_at = time.perf_counter()
    optimum = solve_exact_with_gurobi(
        rows,
        args.objective,
        args.bf_floor,
        args.interference_weight,
        interference_ceiling,
    )
    gurobi_model_seconds = time.perf_counter() - gurobi_started_at
    direct_scan_started_at = time.perf_counter()
    direct_optimum = direct_exact_optimum(
        rows,
        args.objective,
        args.bf_floor,
        args.interference_weight,
        interference_ceiling,
    )
    direct_scan_seconds = time.perf_counter() - direct_scan_started_at
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
    exact_timing = {
        "enumeration_seconds": enumeration_seconds,
        "gurobi_model_seconds": gurobi_model_seconds,
        "direct_scan_seconds": direct_scan_seconds,
        "total_seconds": enumeration_seconds + gurobi_model_seconds,
    }
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
        exact_timing=exact_timing,
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
    print(
        "  timing: "
        f"enumeration={enumeration_seconds:.4f}s, "
        f"gurobi_model={gurobi_model_seconds:.4f}s, "
        f"direct_scan={direct_scan_seconds:.4f}s"
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
