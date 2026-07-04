import time

import numpy as np

from .common import objective_from_gram


def threshold_window_selection(V, K, T):
    V, K = _validate_matrix_and_k(V, K)
    N = V.shape[0]
    T = int(np.clip(int(round(T)), 0, max(0, N - K)))
    x = np.zeros(N, dtype=int)
    if K == 0:
        return x
    if K == N:
        x[:] = 1
        return x

    row_power = np.sum(np.abs(V) ** 2, axis=1).real
    order = np.argsort(row_power)[::-1]
    x[order[T : T + K]] = 1
    return x


def cyclic_threshold_window_selection(V, K, start):
    V, K = _validate_matrix_and_k(V, K)
    N = V.shape[0]
    x = np.zeros(N, dtype=int)
    if K == 0:
        return x
    if K == N:
        x[:] = 1
        return x

    row_power = np.sum(np.abs(V) ** 2, axis=1).real
    order = np.argsort(row_power)[::-1]
    ranks = (int(start) + np.arange(K)) % N
    x[order[ranks]] = 1
    return x


def best_cyclic_threshold_window(V, K, sigma=1.0, P=1.0):
    started_at = time.perf_counter()
    V, K = _validate_matrix_and_k(V, K)
    N, L = V.shape
    if K == 0:
        x = np.zeros(N, dtype=int)
        return _seed_result(x, 0, 1, (0.0, 0.0, 0.0), started_at)
    if K == N:
        x = np.ones(N, dtype=int)
        row_power = np.sum(np.abs(V) ** 2, axis=1).real
        row_grams = V.conj()[:, :, None] * V[:, None, :]
        values = objective_from_gram(
            np.sum(row_grams, axis=0),
            float(np.max(row_power)),
            L,
            sigma=sigma,
            P=P,
        )
        return _seed_result(x, 0, 1, values, started_at)

    row_power = np.sum(np.abs(V) ** 2, axis=1).real
    order = np.argsort(row_power)[::-1]
    row_grams = V.conj()[:, :, None] * V[:, None, :]
    ordered_grams = row_grams[order]
    prefix = np.concatenate(
        [
            np.zeros((1, L, L), dtype=complex),
            np.cumsum(np.concatenate([ordered_grams, ordered_grams], axis=0), axis=0),
        ],
        axis=0,
    )

    best_start = 0
    best_values = None
    best_score = -np.inf
    top_power = float(row_power[order[0]])

    for start in range(N):
        gram = prefix[start + K] - prefix[start]
        max_row_power = top_power if start + K > N else float(row_power[order[start]])
        values = objective_from_gram(
            gram,
            max_row_power,
            L,
            sigma=sigma,
            P=P,
        )
        if values[2] > best_score + 1e-12:
            best_score = values[2]
            best_values = values
            best_start = int(start)

    x = cyclic_threshold_window_selection(V, K, best_start)
    return _seed_result(x, best_start, N, best_values, started_at)


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


def _seed_result(x, seed_position, candidate_count, values, started_at):
    return {
        "x": np.asarray(x, dtype=int),
        "T": int(seed_position),
        "candidate_count": int(candidate_count),
        "u_bf": float(values[0]),
        "u_i": float(values[1]),
        "u_g": float(values[2]),
        "u_g_db": float(10.0 * np.log10(max(float(values[2]), np.finfo(float).tiny))),
        "elapsed_seconds": float(time.perf_counter() - started_at),
    }
