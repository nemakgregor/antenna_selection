import time

import numpy as np

from algorithms import calculate_objectives, check_constraints


def evaluate_solver(name, solver, V, K, sigma=1.0, P=1.0, random_state=None):
    started_at = time.perf_counter()
    with np.errstate(all="ignore"):
        x = solver(V, K, sigma, P, random_state)
        elapsed_seconds = time.perf_counter() - started_at
        valid, active_count = check_constraints(x, K)
        u_bf, u_i, u_g = calculate_objectives(V, x, sigma=sigma, P=P)

    metrics = {
        "valid": bool(valid),
        "active_count": int(active_count),
        "u_bf": float(u_bf),
        "u_i": float(u_i),
        "u_g": float(u_g),
        "elapsed_seconds": float(elapsed_seconds),
    }
    if not metrics["valid"] or not np.isfinite([u_bf, u_i, u_g]).all():
        raise RuntimeError(f"Invalid result for {name}.")
    return x, metrics


def evaluate_algorithms(V, K, sigma=1.0, P=1.0, random_state=None, solvers=None):
    if solvers is None:
        raise ValueError("evaluate_algorithms requires an explicit solvers list.")

    results = {}

    for name, solver in solvers:
        _, metrics = evaluate_solver(name, solver, V, K, sigma, P, random_state)
        results[name] = {
            "valid": metrics["valid"],
            "active_count": metrics["active_count"],
            "u_bf": metrics["u_bf"],
            "u_i": metrics["u_i"],
            "u_g": metrics["u_g"],
        }

    return results
