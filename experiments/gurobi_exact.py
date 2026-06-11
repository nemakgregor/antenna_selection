import argparse
import itertools
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd

from algorithms import (
    calculate_objectives,
    solve_h2,
)
from utils.solver_sets import SMALL_GUROBI_HEURISTICS as HEURISTICS
from utils.data import generate_V
from utils.evaluation import evaluate_solver
from utils.reporting import format_number, split_win_shares, write_markdown
from visualization.gurobi_exact import (
    plot_gurobi_comparison,
    plot_multi_objective_summary,
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


def extra_constraints_status(row, objective, bf_floor, best_bf, interference_ceiling):
    meets_bf_floor = (
        row["u_bf"] >= bf_floor * best_bf
        if objective == "bf_protected_interference"
        else True
    )
    meets_interference_ceiling = (
        row["u_i"] <= interference_ceiling
        if objective == "bf_under_h2_interference"
        else True
    )
    return meets_bf_floor, meets_interference_ceiling


def gap_or_nan(
    row,
    optimum,
    objective,
    best_bf,
    max_interference,
    interference_weight,
    meets_bf_floor,
    meets_interference_ceiling,
):
    if not (meets_bf_floor and meets_interference_ceiling):
        return np.nan
    return objective_gap(
        row,
        optimum,
        objective,
        best_bf,
        max_interference,
        weight=interference_weight,
    )


def gurobi_solution_comparison_row(
    solution,
    optimum,
    objective,
    bf_floor,
    best_bf,
    max_interference,
    interference_weight,
    interference_ceiling,
):
    solution_optimum = solution["optimum"]
    solution_timing = solution.get("timing", {})
    row = {
        "u_bf": solution_optimum["u_bf"],
        "u_i": solution_optimum["u_i"],
        "u_g": solution_optimum["u_g"],
    }
    meets_bf_floor, meets_interference_ceiling = extra_constraints_status(
        row, objective, bf_floor, best_bf, interference_ceiling
    )
    return {
        "algorithm": solution["algorithm"],
        "active_count": solution_optimum["active_count"],
        "valid": True,
        **row,
        "meets_bf_floor": meets_bf_floor,
        "meets_interference_ceiling": meets_interference_ceiling,
        "objective_gap_pct": gap_or_nan(
            row,
            optimum,
            objective,
            best_bf,
            max_interference,
            interference_weight,
            meets_bf_floor,
            meets_interference_ceiling,
        ),
        "subset": " ".join(map(str, solution_optimum["subset"])),
        "elapsed_seconds": solution_timing.get("total_seconds", np.nan),
        "enumeration_seconds": solution_timing.get("enumeration_seconds", np.nan),
        "gurobi_model_seconds": solution_timing.get("gurobi_model_seconds", np.nan),
        "direct_scan_seconds": solution_timing.get("direct_scan_seconds", np.nan),
    }


def heuristic_comparison_row(
    name,
    solver,
    V,
    K,
    sigma,
    P,
    optimum,
    objective,
    bf_floor,
    best_bf,
    max_interference,
    interference_weight,
    interference_ceiling,
):
    x, result = evaluate_solver(name, solver, V, K, sigma, P)
    row = {
        "u_bf": result["u_bf"],
        "u_i": result["u_i"],
        "u_g": result["u_g"],
    }
    meets_bf_floor, meets_interference_ceiling = extra_constraints_status(
        row, objective, bf_floor, best_bf, interference_ceiling
    )
    return {
        "algorithm": name,
        "active_count": result["active_count"],
        "valid": result["valid"],
        **row,
        "meets_bf_floor": meets_bf_floor,
        "meets_interference_ceiling": meets_interference_ceiling,
        "objective_gap_pct": gap_or_nan(
            row,
            optimum,
            objective,
            best_bf,
            max_interference,
            interference_weight,
            meets_bf_floor,
            meets_interference_ceiling,
        ),
        "subset": " ".join(map(str, np.flatnonzero(x))),
        "elapsed_seconds": result["elapsed_seconds"],
        "enumeration_seconds": np.nan,
        "gurobi_model_seconds": np.nan,
        "direct_scan_seconds": np.nan,
    }


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
    if gurobi_solutions is None:
        gurobi_solutions = [
            {
                "algorithm": "Gurobi optimum",
                "optimum": optimum,
                "timing": exact_timing,
            }
        ]

    comparison = [
        gurobi_solution_comparison_row(
            solution,
            optimum,
            objective,
            bf_floor,
            best_bf,
            max_interference,
            interference_weight,
            interference_ceiling,
        )
        for solution in gurobi_solutions
    ]
    comparison.extend(
        heuristic_comparison_row(
            name,
            solver,
            V,
            K,
            sigma,
            P,
            optimum,
            objective,
            bf_floor,
            best_bf,
            max_interference,
            interference_weight,
            interference_ceiling,
        )
        for name, solver in HEURISTICS
    )
    return pd.DataFrame(comparison)


def write_report(comparison, optimum, candidate_count, out_dir, args):
    best_u_g = comparison.loc[comparison["u_g"].idxmax()]
    fastest = comparison.loc[comparison["elapsed_seconds"].idxmin()]
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
        f"- Enumeration: {format_number(comparison.loc[0, 'enumeration_seconds'], precision=6)}s",
        f"- Gurobi model build/solve: {format_number(comparison.loc[0, 'gurobi_model_seconds'], precision=6)}s",
        f"- Direct scan verification: {format_number(comparison.loc[0, 'direct_scan_seconds'], precision=6)}s",
        "",
        "## Summary",
        "",
        f"- Best `U_G`: {best_u_g['algorithm']} (`{format_number(best_u_g['u_g'])}`).",
        f"- Fastest: {fastest['algorithm']} (`{format_number(fastest['elapsed_seconds'], precision=6)}s`).",
        "",
        "Full per-algorithm values are in `gurobi_small_comparison.csv`.",
    ]
    write_markdown(out_dir / "gurobi_small_report.md", lines)


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


def enumerate_multi_sample(args, K, min_active, sample_idx):
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
    return V, seed, rows, enumeration_seconds


def solve_multi_target(args, rows, target_obj, enumeration_seconds, interference_ceiling):
    objective = STANDARD_OBJECTIVES[target_obj]["exact_objective"]
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
    timing = {
        "enumeration_seconds": enumeration_seconds,
        "gurobi_model_seconds": gurobi_model_seconds,
        "direct_scan_seconds": direct_scan_seconds,
        "total_seconds": enumeration_seconds + gurobi_model_seconds,
    }
    return {
        "target_obj": target_obj,
        "algorithm": gurobi_algorithm_name(target_obj),
        "optimum": optimum,
        "timing": timing,
    }


def optima_row(args, target_obj, sample_idx, seed, K, candidate_count, solution):
    optimum = solution["optimum"]
    return {
        "target_obj": target_obj,
        "algorithm": gurobi_algorithm_name(target_obj),
        "exact_objective": STANDARD_OBJECTIVES[target_obj]["exact_objective"],
        "sample": sample_idx,
        "seed": seed,
        "N": args.N,
        "L": args.L,
        "K": K,
        "candidate_count": candidate_count,
        "active_count": optimum["active_count"],
        "u_bf": optimum["u_bf"],
        "u_i": optimum["u_i"],
        "u_g": optimum["u_g"],
        "objective_value": optimum["objective_value"],
        "direct_scan_verified": optimum["direct_scan_verified"],
        "subset": " ".join(map(str, optimum["subset"])),
        **solution["timing"],
    }


def solve_multi_targets(args, objectives, rows, enumeration_seconds, interference_ceiling, sample_idx, total_cases):
    solutions = []
    for objective_pos, target_obj in enumerate(objectives, start=1):
        case_no = sample_idx * len(objectives) + objective_pos
        print(
            f"  [{case_no}/{total_cases}] Solving Gurobi objective={target_obj} "
            f"with {len(rows)} subset variables",
            flush=True,
        )
        solutions.append(
            solve_multi_target(
                args,
                rows,
                target_obj,
                enumeration_seconds,
                interference_ceiling,
            )
        )
    return solutions


def multi_comparison_frame(
    args,
    target_obj,
    sample_idx,
    seed,
    K,
    V,
    rows,
    interference_ceiling,
    solutions,
):
    spec = STANDARD_OBJECTIVES[target_obj]
    target_solution = next(
        solution for solution in solutions if solution["target_obj"] == target_obj
    )
    comparison = compare_heuristics(
        V,
        K,
        args.sigma,
        args.P,
        target_solution["optimum"],
        spec["exact_objective"],
        args.bf_floor,
        rows,
        args.interference_weight,
        interference_ceiling,
        exact_timing=target_solution["timing"],
        gurobi_solutions=solutions,
    )
    comparison.insert(0, "target_obj", target_obj)
    comparison.insert(1, "exact_objective", spec["exact_objective"])
    comparison.insert(2, "sample", sample_idx)
    comparison.insert(3, "seed", seed)
    comparison.insert(4, "N", args.N)
    comparison.insert(5, "L", args.L)
    comparison.insert(6, "K", K)
    comparison["target_metric"] = spec["metric"]
    comparison["target_direction"] = spec["direction"]
    comparison["target_value"] = comparison[spec["metric"]]
    comparison["candidate_count"] = len(rows)
    return comparison


def write_multi_outputs(runs, optima, objectives, args):
    summary = build_multi_summary(runs)
    wins = build_multi_wins(runs)
    optima_df = pd.DataFrame(optima)

    runs.to_csv(args.out_dir / "gurobi_multi_objective_runs.csv", index=False)
    summary.to_csv(args.out_dir / "gurobi_multi_objective_summary.csv", index=False)
    wins.to_csv(args.out_dir / "gurobi_multi_objective_wins.csv", index=False)
    optima_df.to_csv(args.out_dir / "gurobi_multi_objective_optima.csv", index=False)
    plot_multi_objective_summary(summary, objectives, STANDARD_OBJECTIVES, args.out_dir)
    write_multi_objective_report(summary, wins, objectives, args.out_dir, args)

    print(f"wrote multi-objective results to {args.out_dir}")


def run_multi_objective_benchmark(args):
    if args.samples <= 0:
        raise ValueError("--samples must be positive.")

    objectives = list(dict.fromkeys(args.objectives or ["bf", "int", "gen"]))
    K = int(round(args.N * args.active_frac))
    min_active = args.L if args.min_active is None else args.min_active
    validate_exact_args(args, K, min_active)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    run_frames = []
    optima = []
    total_cases = args.samples * len(objectives)

    for sample_idx in range(args.samples):
        V, seed, rows, enumeration_seconds = enumerate_multi_sample(
            args, K, min_active, sample_idx
        )
        interference_ceiling = h2_interference_ceiling(V, K, args)
        solutions = solve_multi_targets(
            args,
            objectives,
            rows,
            enumeration_seconds,
            interference_ceiling,
            sample_idx,
            total_cases,
        )
        optima.extend(
            optima_row(args, solution["target_obj"], sample_idx, seed, K, len(rows), solution)
            for solution in solutions
        )
        run_frames.extend(
            multi_comparison_frame(
                args,
                target_obj,
                sample_idx,
                seed,
                K,
                V,
                rows,
                interference_ceiling,
                solutions,
            )
            for target_obj in objectives
        )

    runs = pd.concat(run_frames, ignore_index=True)
    write_multi_outputs(runs, optima, objectives, args)


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
    def metric_resolver(case, chunk):
        target_obj = case["target_obj"]
        spec = STANDARD_OBJECTIVES[target_obj]
        return [(spec["metric"], spec["label"], spec["direction"])]

    return split_win_shares(
        runs,
        case_cols=["target_obj", "sample"],
        winner_group_cols=["target_obj"],
        sample_group_cols=["target_obj"],
        algorithm_col="algorithm",
        metric_resolver=metric_resolver,
    )


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
        "Full per-algorithm values are in `gurobi_multi_objective_summary.csv`.",
        "",
        "## Objective Leaders",
        "",
    ]

    for target_obj in objectives:
        spec = STANDARD_OBJECTIVES[target_obj]
        data = summary[summary["target_obj"] == target_obj].copy()
        leader = (
            data.loc[data["target_value_mean"].idxmin()]
            if spec["direction"] == "min"
            else data.loc[data["target_value_mean"].idxmax()]
        )
        exact = data[data["algorithm"] == gurobi_algorithm_name(target_obj)].iloc[0]
        win_row = wins[
            (wins["target_obj"] == target_obj)
            & (wins["algorithm"] == exact["algorithm"])
        ]
        win_rate = float(win_row["winner_rate"].iloc[0]) if len(win_row) else 0.0
        lines.append(
            f"- `{target_obj}` ({spec['metric']}): leader {leader['algorithm']} "
            f"with mean `{format_number(leader['target_value_mean'])}`; "
            f"exact {exact['algorithm']} winner rate `{win_rate:.2f}`."
        )

    lines.extend(
        [
            "",
            "## Artifacts",
            "",
            "- `gurobi_multi_objective_runs.csv`: per-run comparisons.",
            "- `gurobi_multi_objective_summary.csv`: aggregate metrics.",
            "- `gurobi_multi_objective_wins.csv`: split winner shares.",
            "- `gurobi_multi_objective_optima.csv`: exact optima.",
            "- `gurobi_multi_objective_summary.png`: summary plot.",
        ]
    )
    write_markdown(out_dir / "gurobi_multi_objective_report.md", lines)


def validate_exact_args(args, K, min_active):
    if not (0 <= min_active <= K <= args.N):
        raise ValueError(
            f"Require 0 <= min_active <= K <= N, got {min_active}, {K}, {args.N}"
        )


def generate_exact_case(args):
    K = int(round(args.N * args.active_frac))
    min_active = args.L if args.min_active is None else args.min_active
    validate_exact_args(args, K, min_active)

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
    return V, K, rows, enumeration_seconds


def h2_interference_ceiling(V, K, args):
    with np.errstate(all="ignore"):
        x_h2 = solve_h2(V, K, sigma=args.sigma, P=args.P)
        _, h2_interference, _ = calculate_objectives(
            V, x_h2, sigma=args.sigma, P=args.P
        )
    return args.interference_ceiling_factor * h2_interference


def solve_verified_exact(rows, args, interference_ceiling):
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
    return optimum, gurobi_model_seconds, direct_scan_seconds


def write_single_objective_outputs(comparison, optimum, candidate_count, args):
    comparison.to_csv(args.out_dir / "gurobi_small_comparison.csv", index=False)
    pd.DataFrame([{**optimum, "subset": " ".join(map(str, optimum["subset"]))}]).to_csv(
        args.out_dir / "gurobi_small_optimum.csv", index=False
    )
    plot_gurobi_comparison(comparison, args.out_dir, args)
    write_report(comparison, optimum, candidate_count, args.out_dir, args)


def print_single_objective_summary(optimum, exact_timing, args):
    print("\nBest exact subset:")
    print(
        f"  active={optimum['active_count']} subset={' '.join(map(str, optimum['subset']))}"
    )
    print(
        f"  BF={optimum['u_bf']:.4f}  Int={optimum['u_i']:.4f}  U_G={optimum['u_g']:.4e}"
    )
    print(
        "  timing: "
        f"enumeration={exact_timing['enumeration_seconds']:.4f}s, "
        f"gurobi_model={exact_timing['gurobi_model_seconds']:.4f}s, "
        f"direct_scan={exact_timing['direct_scan_seconds']:.4f}s"
    )
    print("\nSaved:")
    for path in [
        "gurobi_small_optimum.csv",
        "gurobi_small_comparison.csv",
        "gurobi_small_report.md",
        "gurobi_small_comparison.png",
    ]:
        print(f"  {args.out_dir / path}")


def run_single_objective_benchmark(args):
    args.out_dir.mkdir(parents=True, exist_ok=True)
    V, K, rows, enumeration_seconds = generate_exact_case(args)
    interference_ceiling = h2_interference_ceiling(V, K, args)
    optimum, gurobi_model_seconds, direct_scan_seconds = solve_verified_exact(
        rows, args, interference_ceiling
    )
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
    write_single_objective_outputs(comparison, optimum, len(rows), args)
    print_single_objective_summary(optimum, exact_timing, args)


def main():
    args = parse_args()
    if args.objectives is not None or args.samples > 1:
        run_multi_objective_benchmark(args)
        return

    run_single_objective_benchmark(args)


if __name__ == "__main__":
    main()
