import math
import time

import numpy as np

from .common import objective_from_gram


def cyclic_threshold_window_selection(V, K, start):
    V, K = _validate_matrix_and_k(V, K)
    N = V.shape[0]
    x = np.zeros(N, dtype=int)
    if K == 0:
        return x
    if K == N:
        x[:] = 1
        return x

    start = int(start) % N
    row_powers = np.sum(np.abs(V) ** 2, axis=1).real
    order = np.argsort(row_powers)[::-1]
    ranks = (start + np.arange(K)) % N
    x[order[ranks]] = 1
    return x


def best_cyclic_threshold_window(V, K, sigma=1.0, P=1.0):
    started_at = time.perf_counter()
    V, K = _validate_matrix_and_k(V, K)
    N, L = V.shape
    if K == 0:
        x = np.zeros(N, dtype=int)
        return _seed_result(x, 0, 1, (0.0, 0.0, 0.0), started_at)

    row_powers = np.sum(np.abs(V) ** 2, axis=1).real
    order = np.argsort(row_powers)[::-1]
    row_grams = V[:, :, None].conj() * V[:, None, :]

    starts = [0] if K == N else range(N)
    best = None
    candidate_count = 0
    for start in starts:
        ranks = (int(start) + np.arange(K)) % N
        subset = tuple(sorted(int(value) for value in order[ranks]))
        gram = row_grams[list(subset)].sum(axis=0)
        max_row_power = float(np.max(row_powers[list(subset)]))
        values = objective_from_gram(
            gram,
            max_row_power,
            L,
            sigma=sigma,
            P=P,
        )
        candidate_count += 1
        if best is None or values[2] > best["values"][2] + 1e-12:
            best = {"start": int(start), "subset": subset, "values": values}

    x = np.zeros(N, dtype=int)
    x[list(best["subset"])] = 1
    return _seed_result(x, best["start"], candidate_count, best["values"], started_at)


def threshold_window_selection(V, K, T):
    V, K, T = _validate_and_clip_threshold(V, K, T)
    N = V.shape[0]
    x = np.zeros(N, dtype=int)
    if K == 0:
        return x

    row_powers = np.sum(np.abs(V) ** 2, axis=1).real
    order = np.argsort(row_powers)[::-1]
    x[order[T : T + K]] = 1
    return x


def refine_threshold_by_swaps(
    V,
    K,
    T,
    max_swaps,
    sigma=1.0,
    P=1.0,
    candidate_radius=None,
):
    started_at = time.perf_counter()
    V, K, T = _validate_and_clip_threshold(V, K, T)
    max_swaps = int(max_swaps)
    if max_swaps < 0:
        raise ValueError("max_swaps must be non-negative.")

    x0 = threshold_window_selection(V, K, T)
    row_powers = np.sum(np.abs(V) ** 2, axis=1).real
    order = np.argsort(row_powers)[::-1]
    initial_subset = tuple(int(value) for value in np.flatnonzero(x0))
    active = set(initial_subset)
    radius = _candidate_radius(K, candidate_radius)
    add_pool = _boundary_add_pool(order, active, T, K, radius)
    return _refine_selection_by_swaps(
        V,
        x0,
        max_swaps=max_swaps,
        sigma=sigma,
        P=P,
        candidate_radius=radius,
        candidate_pool=add_pool,
        seed_position=T,
        started_at=started_at,
    )


def refine_selection_by_swaps(
    V,
    x0,
    max_swaps,
    sigma=1.0,
    P=1.0,
    candidate_radius=None,
    candidate_pool=None,
    seed_position=0,
):
    return _refine_selection_by_swaps(
        V,
        x0,
        max_swaps=max_swaps,
        sigma=sigma,
        P=P,
        candidate_radius=candidate_radius,
        candidate_pool=candidate_pool,
        seed_position=seed_position,
        started_at=time.perf_counter(),
    )


def _refine_selection_by_swaps(
    V,
    x0,
    max_swaps,
    sigma=1.0,
    P=1.0,
    candidate_radius=None,
    candidate_pool=None,
    seed_position=0,
    started_at=None,
):
    if started_at is None:
        started_at = time.perf_counter()
    V = np.asarray(V)
    if V.ndim != 2:
        raise ValueError("V must be a 2D complex matrix of shape (N, L).")
    if not np.iscomplexobj(V):
        V = V.astype(np.complex128)
    x0 = np.asarray(x0, dtype=int)
    if x0.ndim != 1 or x0.shape[0] != V.shape[0]:
        raise ValueError("x0 must be a 1D selection vector with length N.")
    if not np.isin(x0, [0, 1]).all():
        raise ValueError("x0 must be a binary selection vector.")

    max_swaps = int(max_swaps)
    if max_swaps < 0:
        raise ValueError("max_swaps must be non-negative.")

    N, L = V.shape
    K = int(np.sum(x0))
    initial_subset = tuple(int(value) for value in np.flatnonzero(x0))
    if K == 0:
        return _result(
            np.zeros(N, dtype=int),
            initial_subset,
            initial_subset,
            seed_position,
            max_swaps,
            0 if candidate_radius is None else int(candidate_radius),
            0,
            0,
            0,
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
            [],
            started_at,
        )

    row_powers = np.sum(np.abs(V) ** 2, axis=1).real
    order = np.argsort(row_powers)[::-1]
    row_grams = V[:, :, None].conj() * V[:, None, :]

    active = set(initial_subset)
    current_gram = row_grams[list(active)].sum(axis=0)
    current_max_power = float(max(row_powers[list(active)]))
    current_values = objective_from_gram(
        current_gram,
        current_max_power,
        L,
        sigma=sigma,
        P=P,
    )
    initial_values = current_values

    if candidate_pool is None:
        radius = -1 if candidate_radius is None else int(max(0, candidate_radius))
        add_pool = _rank_neighborhood_add_pool(order, active, radius)
    else:
        radius = -1 if candidate_radius is None else int(max(0, candidate_radius))
        add_pool = _dedupe_pool(candidate_pool, active)
    evaluated_swap_count = 0
    swap_history = []

    for _pass_no in range(max_swaps):
        active_idx = tuple(sorted(active))
        inactive_add = tuple(index for index in add_pool if index not in active)
        if not active_idx or not inactive_add:
            break

        max_idx, max_power, second_power = _active_power_leaders(
            active_idx,
            row_powers,
        )
        best_swap = None
        best_values = current_values
        best_gram = None
        best_max_power = current_max_power

        for remove_idx in active_idx:
            gram_without = current_gram - row_grams[remove_idx]
            max_after_remove = second_power if remove_idx == max_idx else max_power
            for add_idx in inactive_add:
                candidate_gram = gram_without + row_grams[add_idx]
                candidate_max_power = float(max(max_after_remove, row_powers[add_idx]))
                values = objective_from_gram(
                    candidate_gram,
                    candidate_max_power,
                    L,
                    sigma=sigma,
                    P=P,
                )
                evaluated_swap_count += 1
                if values[2] > best_values[2] + 1e-12:
                    best_values = values
                    best_swap = (int(remove_idx), int(add_idx))
                    best_gram = candidate_gram
                    best_max_power = candidate_max_power

        if best_swap is None:
            break

        remove_idx, add_idx = best_swap
        active.remove(remove_idx)
        active.add(add_idx)
        current_gram = best_gram
        current_max_power = best_max_power
        current_values = best_values
        swap_history.append(
            {
                "remove": remove_idx,
                "add": add_idx,
                "u_g": float(current_values[2]),
            }
        )

    final_subset = tuple(sorted(active))
    x = np.zeros(N, dtype=int)
    x[list(final_subset)] = 1
    return _result(
        x,
        initial_subset,
        final_subset,
        seed_position,
        max_swaps,
        radius,
        len(add_pool),
        evaluated_swap_count,
        len(swap_history),
        initial_values,
        current_values,
        swap_history,
        started_at,
    )


def evaluate_threshold_local_rules(
    V,
    K,
    seed_rules,
    max_swaps_values=(0, 1, 2),
    sigma=1.0,
    P=1.0,
    candidate_radius=None,
):
    rows = []
    for seed_rule in seed_rules:
        label, T = _seed_rule_parts(seed_rule)
        for max_swaps in max_swaps_values:
            result = refine_threshold_by_swaps(
                V,
                K,
                T,
                max_swaps=max_swaps,
                sigma=sigma,
                P=P,
                candidate_radius=candidate_radius,
            )
            rows.append(
                {
                    "seed_rule": label,
                    "seed_T": int(result["T"]),
                    "max_swaps": int(max_swaps),
                    "active_count": int(np.sum(result["x"])),
                    "candidate_kind": "threshold_window"
                    if int(max_swaps) == 0
                    else "threshold_window_local_swap",
                    "candidate_count": int(result["candidate_count"]),
                    "candidate_radius": int(result["candidate_radius"]),
                    "add_candidate_count": int(result["add_candidate_count"]),
                    "evaluated_swap_count": int(result["evaluated_swap_count"]),
                    "swaps_applied": int(result["swaps_applied"]),
                    "initial_subset": _subset_to_string(result["initial_subset"]),
                    "subset": _subset_to_string(result["subset"]),
                    "swap_history": result["swap_history"],
                    "initial_u_bf": float(result["initial_u_bf"]),
                    "initial_u_i": float(result["initial_u_i"]),
                    "initial_u_g": float(result["initial_u_g"]),
                    "u_bf": float(result["u_bf"]),
                    "u_i": float(result["u_i"]),
                    "u_g": float(result["u_g"]),
                    "u_g_db": 10.0
                    * np.log10(max(float(result["u_g"]), np.finfo(float).tiny)),
                    "score": float(result["u_g"]),
                    "elapsed_seconds": float(result["elapsed_seconds"]),
                }
            )
    return rows


def _validate_and_clip_threshold(V, K, T):
    V, K = _validate_matrix_and_k(V, K)
    N = V.shape[0]
    max_T = max(0, N - K)
    T = int(np.clip(int(round(T)), 0, max_T))
    return V, K, T


def _validate_matrix_and_k(V, K):
    V = np.asarray(V)
    if V.ndim != 2:
        raise ValueError("V must be a 2D complex matrix of shape (N, L).")
    if not np.iscomplexobj(V):
        V = V.astype(np.complex128)

    N = V.shape[0]
    K = int(K)
    if not (0 <= K <= N):
        raise ValueError("K must satisfy 0 <= K <= N.")
    return V, K


def _candidate_radius(K, candidate_radius):
    if candidate_radius is None:
        return int(max(8, math.ceil(0.05 * int(K))))
    return int(max(0, candidate_radius))


def _boundary_add_pool(order, active, T, K, radius):
    N = len(order)
    start = max(0, int(T) - int(radius))
    stop = min(N, int(T) + int(K) + int(radius))
    add_pool = []
    seen = set()
    for rank in range(start, stop):
        index = int(order[rank])
        if index not in active and index not in seen:
            add_pool.append(index)
            seen.add(index)
    return tuple(add_pool)


def _rank_neighborhood_add_pool(order, active, radius):
    if radius < 0:
        return tuple(int(index) for index in order if int(index) not in active)
    ranks = {int(rank) for rank, index in enumerate(order) if int(index) in active}
    pool = []
    seen = set()
    N = len(order)
    for rank in sorted(ranks):
        start = max(0, rank - radius)
        stop = min(N, rank + radius + 1)
        for candidate_rank in range(start, stop):
            index = int(order[candidate_rank])
            if index not in active and index not in seen:
                pool.append(index)
                seen.add(index)
    return tuple(pool)


def _dedupe_pool(candidate_pool, active):
    pool = []
    seen = set()
    for value in candidate_pool:
        index = int(value)
        if index not in active and index not in seen:
            pool.append(index)
            seen.add(index)
    return tuple(pool)


def _active_power_leaders(active_idx, row_powers):
    powers = [(int(index), float(row_powers[int(index)])) for index in active_idx]
    powers.sort(key=lambda item: item[1], reverse=True)
    max_idx, max_power = powers[0]
    second_power = powers[1][1] if len(powers) > 1 else 0.0
    return max_idx, max_power, second_power


def _seed_rule_parts(seed_rule):
    if isinstance(seed_rule, dict):
        return str(seed_rule["seed_rule"]), int(seed_rule["T"])
    if isinstance(seed_rule, (tuple, list)) and len(seed_rule) == 2:
        return str(seed_rule[0]), int(seed_rule[1])
    raise ValueError("seed_rules must contain dicts or (label, T) pairs.")


def _subset_to_string(subset):
    return " ".join(str(int(value)) for value in subset)


def _swap_history_to_string(history):
    return ";".join(
        f"{int(item['remove'])}>{int(item['add'])}@{float(item['u_g']):.12g}"
        for item in history
    )


def _seed_result(x, seed_position, candidate_count, values, started_at):
    subset = tuple(int(value) for value in np.flatnonzero(x))
    return {
        "x": x,
        "T": int(seed_position),
        "initial_subset": subset,
        "subset": subset,
        "candidate_count": int(candidate_count),
        "u_bf": float(values[0]),
        "u_i": float(values[1]),
        "u_g": float(values[2]),
        "u_g_db": 10.0 * np.log10(max(float(values[2]), np.finfo(float).tiny)),
        "elapsed_seconds": float(time.perf_counter() - started_at),
    }


def _result(
    x,
    initial_subset,
    final_subset,
    T,
    max_swaps,
    candidate_radius,
    add_candidate_count,
    evaluated_swap_count,
    swaps_applied,
    initial_values,
    final_values,
    swap_history,
    started_at,
):
    return {
        "x": x,
        "T": int(T),
        "max_swaps": int(max_swaps),
        "initial_subset": tuple(initial_subset),
        "subset": tuple(final_subset),
        "candidate_count": 1 + int(evaluated_swap_count),
        "candidate_radius": int(candidate_radius),
        "add_candidate_count": int(add_candidate_count),
        "evaluated_swap_count": int(evaluated_swap_count),
        "swaps_applied": int(swaps_applied),
        "swap_history": _swap_history_to_string(swap_history),
        "initial_u_bf": float(initial_values[0]),
        "initial_u_i": float(initial_values[1]),
        "initial_u_g": float(initial_values[2]),
        "u_bf": float(final_values[0]),
        "u_i": float(final_values[1]),
        "u_g": float(final_values[2]),
        "score": float(final_values[2]),
        "elapsed_seconds": float(time.perf_counter() - started_at),
    }
