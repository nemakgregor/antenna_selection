import time

import numpy as np

from algorithms import calculate_objectives, check_constraints


def evaluate_solver(name, solver, V, K, sigma=1.0, P=1.0, random_state=None):
    started_at = time.perf_counter()
    with np.errstate(all="ignore"):
        kwargs = {"sigma": sigma, "P": P}
        if random_state is not None:
            kwargs["random_state"] = random_state
        x = solver(V, K, **kwargs)
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
