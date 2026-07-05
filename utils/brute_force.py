import itertools
import math
import time

import numpy as np

from algorithms.common import objective_from_gram


def selection_vector_from_subset(N, subset):
    x = np.zeros(int(N), dtype=int)
    if subset:
        x[np.asarray(subset, dtype=int)] = 1
    return x


def subset_to_string(subset):
    return " ".join(str(int(value)) for value in subset)


def parse_subset_string(value):
    if value is None:
        return tuple()
    try:
        if np.isnan(value):
            return tuple()
    except TypeError:
        pass
    text = str(value).strip()
    if not text:
        return tuple()
    return tuple(int(part) for part in text.split())


def brute_force_exact_u_g(
    V,
    K,
    sigma=1.0,
    P=1.0,
    time_limit_seconds=None,
    timeout_check_interval=1024,
):
    V = np.asarray(V)
    if V.ndim != 2:
        raise ValueError("V must be a 2D complex matrix.")
    if not np.iscomplexobj(V):
        V = V.astype(np.complex128)

    N, L = V.shape
    K = int(K)
    if not (0 <= K <= N):
        raise ValueError("K must satisfy 0 <= K <= N.")

    started_at = time.perf_counter()
    candidate_count = math.comb(N, K)
    if K == 0:
        return {
            "x": np.zeros(N, dtype=int),
            "subset": tuple(),
            "u_bf": 0.0,
            "u_i": 0.0,
            "u_g": 0.0,
            "candidate_count": int(candidate_count),
            "evaluated_count": 1,
            "elapsed_seconds": time.perf_counter() - started_at,
            "timed_out": False,
            "completed": True,
        }

    row_grams = V[:, :, None].conj() * V[:, None, :]
    row_powers = np.sum(np.abs(V) ** 2, axis=1).real
    best_subset = None
    best_values = (np.nan, np.nan, -np.inf)
    evaluated_count = 0
    timed_out = False
    batch_size = max(1, int(timeout_check_interval))
    combinations_iter = itertools.combinations(range(N), K)
    eye = np.eye(L)

    while True:
        if (
            time_limit_seconds is not None
            and evaluated_count > 0
            and time.perf_counter() - started_at > float(time_limit_seconds)
        ):
            timed_out = True
            break

        batch = list(itertools.islice(combinations_iter, batch_size))
        if not batch:
            break

        idx = np.asarray(batch, dtype=np.int64)
        grams = row_grams[idx].sum(axis=1)
        max_row_powers = row_powers[idx].max(axis=1)
        valid = max_row_powers > 0
        z2 = np.zeros(len(batch), dtype=float)
        z2[valid] = float(P) / max_row_powers[valid]

        gram_sq = grams @ np.swapaxes(grams.conj(), 1, 2)
        traces = np.real(np.trace(gram_sq, axis1=1, axis2=2))
        diag_values = np.diagonal(gram_sq, axis1=1, axis2=2)
        offdiag_energy = (
            np.sum(np.abs(gram_sq) ** 2, axis=(1, 2))
            - np.sum(np.abs(diag_values) ** 2, axis=1)
        )
        u_bf_values = z2 * traces
        u_i_values = (z2**2) * offdiag_energy
        det_matrices = z2[:, None, None] * gram_sq + float(sigma) * eye[None, :, :]
        u_g_values = np.real(np.linalg.det(det_matrices))
        u_bf_values[~valid] = 0.0
        u_i_values[~valid] = 0.0
        u_g_values[~valid] = 0.0

        batch_best_pos = int(np.argmax(u_g_values))
        batch_best_u_g = float(u_g_values[batch_best_pos])
        if batch_best_u_g > best_values[2]:
            best_subset = tuple(int(value) for value in batch[batch_best_pos])
            best_values = (
                float(u_bf_values[batch_best_pos]),
                float(u_i_values[batch_best_pos]),
                batch_best_u_g,
            )
        evaluated_count += len(batch)

    completed = (not timed_out) and evaluated_count == candidate_count
    elapsed_seconds = time.perf_counter() - started_at
    if best_subset is None:
        best_subset = tuple()
        best_values = (np.nan, np.nan, np.nan)

    return {
        "x": selection_vector_from_subset(N, best_subset),
        "subset": best_subset,
        "u_bf": float(best_values[0]),
        "u_i": float(best_values[1]),
        "u_g": float(best_values[2]),
        "candidate_count": int(candidate_count),
        "evaluated_count": int(evaluated_count),
        "elapsed_seconds": float(elapsed_seconds),
        "timed_out": bool(timed_out),
        "completed": bool(completed),
    }


def contiguous_threshold_window_T(V, subset):
    V = np.asarray(V)
    N = V.shape[0]
    subset = tuple(sorted(int(value) for value in subset))
    K = len(subset)
    if K == 0:
        return 0
    if K > N:
        return None

    p_n = np.sum(np.abs(V) ** 2, axis=1).real
    idx_desc = np.argsort(p_n)[::-1]
    rank_by_index = np.empty(N, dtype=int)
    rank_by_index[idx_desc] = np.arange(N)
    ranks = sorted(int(rank_by_index[index]) for index in subset)
    start = ranks[0]
    if ranks == list(range(start, start + K)) and start <= N - K:
        return int(start)
    return None


def threshold_window_subset_string(V, K, T):
    V = np.asarray(V)
    p_n = np.sum(np.abs(V) ** 2, axis=1).real
    idx_desc = np.argsort(p_n)[::-1]
    active_idx = idx_desc[int(T) : int(T) + int(K)]
    return subset_to_string(tuple(sorted(int(value) for value in active_idx)))
