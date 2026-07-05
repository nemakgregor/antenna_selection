import time

import numpy as np

from .common import calculate_objectives, objective_from_gram


LEGACY_RAW_THRESHOLDS = (1, 2, 5, 10, 25, 50)
PERCENT_THRESHOLDS = (0.005, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20)
QUANTILE_THRESHOLDS = (0.80, 0.90, 0.95, 0.99)


def legacy_thresholds(N, K):
    return _dedupe_thresholds([*LEGACY_RAW_THRESHOLDS, int(N) // 10], N, K)


def threshold_power_window(V, K, T):
    V, K, T = _validate_inputs(V, K, T)
    N = V.shape[0]
    p_n = np.sum(np.abs(V) ** 2, axis=1).real
    idx_desc = np.argsort(p_n)[::-1]

    x = np.zeros(N, dtype=int)
    x[idx_desc[T : T + K]] = 1
    return x


def dense_thresholds(K):
    K = int(K)
    if K < 0:
        raise ValueError("K must be non-negative.")
    return list(range(K + 1))


def evaluate_power_window_thresholds(V, K, thresholds, sigma=1.0, P=1.0):
    """
    Evaluate shifted power windows for a dense threshold experiment.

    This isolates the threshold parameter: for each T, select rows
    idx_desc[T:T+K], where idx_desc sorts rows by descending row power.
    """

    V = np.asarray(V)
    if V.ndim != 2:
        raise ValueError("V must be a 2D complex matrix of shape (N, L).")
    if not np.iscomplexobj(V):
        V = V.astype(np.complex128)

    N, L = V.shape
    K = int(K)
    if not (0 <= K <= N):
        raise ValueError("K must satisfy 0 <= K <= N.")

    thresholds = [int(T) for T in thresholds]
    max_T = N - K
    invalid = [T for T in thresholds if not (0 <= T <= max_T)]
    if invalid:
        raise ValueError(
            f"All thresholds must satisfy 0 <= T <= N-K={max_T}; "
            f"got {invalid[:5]}."
        )

    if K == 0:
        rows = []
        for T in thresholds:
            started_at = time.perf_counter()
            rows.append(
                _threshold_row(
                    T,
                    0,
                    "power_window",
                    0.0,
                    0.0,
                    0.0,
                    time.perf_counter() - started_at,
                )
            )
        return rows

    p_n = np.sum(np.abs(V) ** 2, axis=1).real
    idx_desc = np.argsort(p_n)[::-1]
    V_sorted = V[idx_desc]
    p_sorted = p_n[idx_desc]

    row_grams = V_sorted[:, :, None].conj() * V_sorted[:, None, :]
    prefix = np.zeros((N + 1, L, L), dtype=np.complex128)
    prefix[1:] = np.cumsum(row_grams, axis=0)

    rows = []
    for T in thresholds:
        started_at = time.perf_counter()
        gram = prefix[T + K] - prefix[T]
        max_row_power = float(p_sorted[T]) if K > 0 else 0.0
        u_bf, u_i, u_g = objective_from_gram(
            gram,
            max_row_power,
            L,
            sigma=sigma,
            P=P,
        )
        rows.append(
            _threshold_row(
                T,
                K,
                "power_window",
                u_bf,
                u_i,
                u_g,
                time.perf_counter() - started_at,
            )
        )
    return rows


def row_power_distribution_metrics(V):
    V = np.asarray(V)
    p_n = np.sum(np.abs(V) ** 2, axis=1).real
    eps = np.finfo(float).eps
    total = float(np.sum(p_n))
    mean = float(np.mean(p_n))
    std = float(np.std(p_n))
    centered = p_n - mean
    skew = float(np.mean(centered**3) / max(std**3, eps)) if std > 0 else 0.0

    p50 = float(np.quantile(p_n, 0.50))
    p80 = float(np.quantile(p_n, 0.80))
    p90 = float(np.quantile(p_n, 0.90))
    p95 = float(np.quantile(p_n, 0.95))
    p99 = float(np.quantile(p_n, 0.99))
    pmax = float(np.max(p_n)) if len(p_n) else 0.0

    sorted_log = np.log(np.maximum(np.sort(p_n)[::-1], eps))
    if len(sorted_log) > 1:
        gaps = sorted_log[:-1] - sorted_log[1:]
        gap_idx = int(np.argmax(gaps))
        gap_value = float(gaps[gap_idx])
        gap_rank = gap_idx + 1
    else:
        gap_value = 0.0
        gap_rank = 0

    return {
        "row_power_mean": mean,
        "row_power_std": std,
        "row_power_cv": std / max(mean, eps),
        "row_power_skew3": skew,
        "row_power_p50": p50,
        "row_power_p80": p80,
        "row_power_p90": p90,
        "row_power_p95": p95,
        "row_power_p99": p99,
        "row_power_max": pmax,
        "row_power_p95_p50": p95 / max(p50, eps),
        "row_power_p99_p50": p99 / max(p50, eps),
        "row_power_max_p95": pmax / max(p95, eps),
        "log_power_gap_max": gap_value,
        "log_power_gap_rank": gap_rank,
        "log_power_gap_rank_pct": gap_rank / max(len(p_n), 1),
        "tail_mass_p80": _tail_mass(p_n, p80, total),
        "tail_mass_p90": _tail_mass(p_n, p90, total),
        "tail_mass_p95": _tail_mass(p_n, p95, total),
        "tail_mass_p99": _tail_mass(p_n, p99, total),
    }


def build_threshold_grid(V, K):
    V = np.asarray(V)
    N = V.shape[0]
    K = int(K)
    max_T = max(0, N - K)

    records = []

    def add(T, source):
        T = int(np.clip(int(round(T)), 0, max_T))
        records.append({"T": T, "threshold_source": source})

    add(0, "top_power")
    for T in legacy_thresholds(N, K):
        add(T, "legacy")
    for pct in PERCENT_THRESHOLDS:
        T = max(1, int(round(max_T * pct))) if max_T > 0 else 0
        add(T, f"off_pct_{100.0 * pct:g}")

    p_n = np.sum(np.abs(V) ** 2, axis=1).real
    eps = np.finfo(float).eps
    sorted_log = np.log(np.maximum(np.sort(p_n)[::-1], eps))
    if len(sorted_log) > 1:
        gaps = sorted_log[:-1] - sorted_log[1:]
        valid_count = min(max_T, len(gaps))
        if valid_count > 0:
            top_gap_indices = np.argsort(gaps[:valid_count])[::-1][:3]
            for rank, gap_idx in enumerate(top_gap_indices, start=1):
                add(int(gap_idx) + 1, f"gap_top{rank}")

    for quantile in QUANTILE_THRESHOLDS:
        cutoff = float(np.quantile(p_n, quantile))
        add(int(np.sum(p_n >= cutoff)), f"quantile_{quantile:g}")

    mean = float(np.mean(p_n))
    std = float(np.std(p_n))
    for scale in (1.0, 1.5, 2.0):
        add(int(np.sum(p_n >= mean + scale * std)), f"mean_plus_{scale:g}std")

    return _merge_threshold_records(records, N, K)


def evaluate_threshold_T(V, K, T, target_obj="gen", sigma=1.0, P=1.0):
    started_at = time.perf_counter()
    V, K, T = _validate_inputs(V, K, T)

    if target_obj not in {"bf", "int", "gen"}:
        raise ValueError("target_obj must be one of {'bf', 'int', 'gen'}.")

    N, L = V.shape
    if K == 0:
        x = np.zeros(N, dtype=int)
        u_bf, u_i, u_g = calculate_objectives(V, x, sigma=sigma, P=P)
        return _result(x, "empty", 1, u_bf, u_i, u_g, target_obj, started_at)
    if K == N:
        x = np.ones(N, dtype=int)
        u_bf, u_i, u_g = calculate_objectives(V, x, sigma=sigma, P=P)
        return _result(x, "all", 1, u_bf, u_i, u_g, target_obj, started_at)

    p_n = np.sum(np.abs(V) ** 2, axis=1).real
    idx_desc = np.argsort(p_n)[::-1]
    row_interference = V.conj()[:, :, None] * V[:, None, :]
    diag = np.arange(L)
    row_interference[:, diag, diag] = 0

    def phase_nulling_subset(pool_indices, lock_index, target_drop_count):
        pool_indices = np.asarray(pool_indices, dtype=int)
        target_active_count = len(pool_indices) - int(target_drop_count)
        target_active_count = int(np.clip(target_active_count, 0, len(pool_indices)))

        x = np.zeros(N, dtype=int)
        x[pool_indices] = 1
        active_count = len(pool_indices)
        S = np.sum(row_interference[pool_indices], axis=0)
        locked = (
            np.array([], dtype=int)
            if lock_index is None
            else np.asarray(lock_index, dtype=int)
        )

        while active_count > target_active_count:
            active_idx = np.where(x == 1)[0]
            if locked.size:
                active_idx = np.setdiff1d(active_idx, locked, assume_unique=False)
            if len(active_idx) == 0:
                break

            residuals = S[None, :, :] - row_interference[active_idx]
            interfs = np.sum(np.abs(residuals) ** 2, axis=(1, 2))
            best_n = active_idx[int(np.argmin(interfs))]

            x[best_n] = 0
            S -= row_interference[best_n]
            active_count -= 1

        return x

    candidates = []
    if target_obj in {"bf", "gen"}:
        x_power = np.zeros(N, dtype=int)
        x_power[idx_desc[T : T + K]] = 1
        candidates.append(("power_window", x_power))

    if target_obj in {"int", "gen"}:
        pool = idx_desc[T:]
        if len(pool) > K:
            lock_anchor = None if T == 0 else [int(pool[0])]
            candidates.append(
                (
                    "tail_phase_null",
                    phase_nulling_subset(
                        pool,
                        lock_index=lock_anchor,
                        target_drop_count=len(pool) - K,
                    ),
                )
            )

    if target_obj == "gen":
        buffer_size = min(30, len(idx_desc[T:]) - K)
        if buffer_size > 0:
            pool = idx_desc[T : T + K + buffer_size]
            lock_anchor = int(pool[0])
            candidates.append(
                (
                    "buffered_phase_null",
                    phase_nulling_subset(
                        pool,
                        lock_index=[lock_anchor],
                        target_drop_count=buffer_size,
                    ),
                )
            )

    if not candidates:
        raise RuntimeError(f"No threshold candidates generated for T={T}.")

    best_x = None
    best_kind = None
    best_values = None
    best_score = -np.inf

    for kind, x in candidates:
        u_bf, u_i, u_g = calculate_objectives(V, x, sigma=sigma, P=P)
        if target_obj == "bf":
            score = u_bf
        elif target_obj == "int":
            score = -u_i
        else:
            score = u_g

        if score > best_score:
            best_score = score
            best_x = x.copy()
            best_kind = kind
            best_values = (u_bf, u_i, u_g)

    return _result(
        best_x,
        best_kind,
        len(candidates),
        *best_values,
        target_obj,
        started_at,
    )


def _validate_inputs(V, K, T):
    V = np.asarray(V)
    if V.ndim != 2:
        raise ValueError("V must be a 2D complex matrix of shape (N, L).")
    if not np.iscomplexobj(V):
        V = V.astype(np.complex128)

    N = V.shape[0]
    K = int(K)
    T = int(T)
    if not (0 <= K <= N):
        raise ValueError("K must satisfy 0 <= K <= N.")
    if not (0 <= T <= N - K):
        raise ValueError("T must satisfy 0 <= T <= N - K.")
    return V, K, T


def _dedupe_thresholds(values, N, K):
    max_T = max(0, int(N) - int(K))
    result = []
    seen = set()
    for value in values:
        T = int(value)
        if 0 < T <= max_T and T not in seen:
            seen.add(T)
            result.append(T)
    return result


def _merge_threshold_records(records, N, K):
    max_T = max(0, int(N) - int(K))
    by_T = {}
    for record in records:
        T = int(np.clip(record["T"], 0, max_T))
        by_T.setdefault(T, []).append(record["threshold_source"])
    return [
        {"T": T, "threshold_source": "+".join(dict.fromkeys(sources))}
        for T, sources in sorted(by_T.items())
    ]


def _tail_mass(values, cutoff, total):
    if total <= 0:
        return 0.0
    return float(np.sum(values[values >= cutoff]) / total)


def _result(x, kind, candidate_count, u_bf, u_i, u_g, target_obj, started_at):
    if target_obj == "bf":
        score = u_bf
    elif target_obj == "int":
        score = -u_i
    else:
        score = u_g

    return {
        "x": x,
        "candidate_kind": kind,
        "candidate_count": int(candidate_count),
        "u_bf": float(u_bf),
        "u_i": float(u_i),
        "u_g": float(u_g),
        "score": float(score),
        "elapsed_seconds": float(time.perf_counter() - started_at),
    }


def _threshold_row(T, active_count, candidate_kind, u_bf, u_i, u_g, elapsed_seconds):
    return {
        "T": int(T),
        "candidate_kind": candidate_kind,
        "candidate_count": 1,
        "active_count": int(active_count),
        "u_bf": float(u_bf),
        "u_i": float(u_i),
        "u_g": float(u_g),
        "u_g_db": float(10.0 * np.log10(max(u_g, np.finfo(float).tiny))),
        "score": float(u_g),
        "elapsed_seconds": float(elapsed_seconds),
    }
