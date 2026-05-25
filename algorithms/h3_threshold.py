import numpy as np

from .common import calculate_objectives


def solve_h3(V, K, target_obj="gen", sigma=1.0, P=1.0):
    """
    H3: standalone threshold-swept projection.

    Runs the same threshold candidate family under three possible target
    objectives. The selected vector is still evaluated later on all metrics.
    """

    if target_obj not in {"bf", "int", "gen"}:
        raise ValueError("target_obj must be one of {'bf', 'int', 'gen'}.")

    V = np.asarray(V)
    if V.ndim != 2:
        raise ValueError("V must be a 2D complex matrix of shape (N, L).")
    if not np.iscomplexobj(V):
        V = V.astype(np.complex128)

    N, L = V.shape
    K = int(K)

    if not (0 <= K <= N):
        raise ValueError("K must satisfy 0 <= K <= N.")
    if K == 0:
        return np.zeros(N, dtype=int)
    if K == N:
        return np.ones(N, dtype=int)

    p_n = np.sum(np.abs(V) ** 2, axis=1).real
    idx_desc = np.argsort(p_n)[::-1]
    row_interference = V.conj()[:, :, None] * V[:, None, :]
    diag = np.arange(L)
    row_interference[:, diag, diag] = 0

    candidates = []

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

    x_power_bound = np.zeros(N, dtype=int)
    x_power_bound[idx_desc[:K]] = 1
    candidates.append(x_power_bound)

    if N > K:
        candidates.append(
            phase_nulling_subset(
                np.arange(N),
                lock_index=None,
                target_drop_count=N - K,
            )
        )

    raw_T_tests = [1, 2, 5, 10, 25, 50, N // 10]
    T_tests = []
    seen = set()
    for T in raw_T_tests:
        T = int(T)
        if 0 < T <= N - K and T not in seen:
            T_tests.append(T)
            seen.add(T)

    for T in T_tests:
        if target_obj in {"bf", "gen"}:
            x_cand = np.zeros(N, dtype=int)
            x_cand[idx_desc[T : T + K]] = 1
            candidates.append(x_cand)

        if target_obj in {"int", "gen"}:
            pool = idx_desc[T:]
            if len(pool) > K:
                lock_anchor = int(pool[0])
                candidates.append(
                    phase_nulling_subset(
                        pool,
                        lock_index=[lock_anchor],
                        target_drop_count=len(pool) - K,
                    )
                )

        if target_obj == "gen":
            buffer_size = min(30, len(idx_desc[T:]) - K)
            if buffer_size > 0:
                pool = idx_desc[T : T + K + buffer_size]
                lock_anchor = int(pool[0])
                candidates.append(
                    phase_nulling_subset(
                        pool,
                        lock_index=[lock_anchor],
                        target_drop_count=buffer_size,
                    )
                )

    best_x = None
    best_score = -np.inf

    for x in candidates:
        bf, ui, gen = calculate_objectives(V, x, sigma=sigma, P=P)
        if target_obj == "bf":
            score = bf
        elif target_obj == "int":
            score = -ui
        else:
            score = gen

        if score > best_score:
            best_score = score
            best_x = x.copy()

    return best_x
