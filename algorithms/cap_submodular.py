import heapq

import numpy as np


DEFAULT_SCAN_SIZE = 17
DEFAULT_PORTFOLIO_WINDOW_SCAN_SIZE = 129
DEFAULT_PORTFOLIO_REFINE_CAPS = 1


def solve_cap_submodular_gen(
    V,
    K,
    sigma=1.0,
    P=1.0,
    random_state=None,
    scan_size=DEFAULT_SCAN_SIZE,
):
    """
    Cap-aware submodular relaxation for the general objective.

    For a fixed row-power cap tau, replacing the active max row power by tau
    gives the fixed-cap log objective

        epsilon_tau = sqrt(sigma * tau / P),
        h_tau(A) = log det(sigma I + (P / tau) G_A^2).

    If lambda_i are the eigenvalues of G_A, then h_tau(A) differs by at most
    L log 2 from L log(sigma) + 2 f_tau(A), where

        f_tau(A) = log det(I + G_A / epsilon_tau).

    The solver greedily maximizes f_tau for a compact set of caps, evaluates
    each candidate by the true objective, and returns the best true-scoring
    candidate.
    """

    del random_state

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
    if sigma <= 0.0:
        raise ValueError("sigma must be positive for the cap-submodular relaxation.")
    if P <= 0.0:
        raise ValueError("P must be positive for the cap-submodular relaxation.")

    row_power = np.sum(np.abs(V) ** 2, axis=1).real
    row_grams = V.conj()[:, :, None] * V[:, None, :]
    eye = np.eye(L, dtype=complex)

    best_active = None
    best_score = -np.inf

    for tau in _candidate_thresholds(row_power, K, scan_size):
        if tau <= 0.0:
            continue
        pool = np.flatnonzero(row_power <= tau + max(1e-12, abs(tau) * 1e-12))
        if len(pool) < K:
            continue

        epsilon = np.sqrt(float(sigma) * float(tau) / float(P))
        active = _greedy_logdet_cap(V, pool, K, epsilon, eye)
        score = _true_log_general_score(row_grams, row_power, active, sigma, P, eye)
        if score > best_score:
            best_score = score
            best_active = active

    if best_active is None:
        best_active = np.zeros(N, dtype=bool)
        best_active[np.argsort(row_power)[-K:]] = True
    return best_active.astype(int)


def solve_cap_submodular_lazy_gen(
    V,
    K,
    sigma=1.0,
    P=1.0,
    random_state=None,
    scan_size=DEFAULT_SCAN_SIZE,
):
    """
    Lazy-greedy cap-aware submodular relaxation for the general objective.

    This is the same fixed-cap greedy objective as solve_cap_submodular_gen, but
    it uses the accelerated greedy rule of Minoux: because f_tau is submodular,
    previously computed marginal gains are valid upper bounds after later
    insertions.  A heap of those upper bounds avoids recomputing every remaining
    antenna at every greedy step while preserving the exact greedy choice.
    """

    del random_state

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
    if sigma <= 0.0:
        raise ValueError("sigma must be positive for the cap-submodular relaxation.")
    if P <= 0.0:
        raise ValueError("P must be positive for the cap-submodular relaxation.")

    row_power = np.sum(np.abs(V) ** 2, axis=1).real
    row_grams = V.conj()[:, :, None] * V[:, None, :]
    eye = np.eye(L, dtype=complex)

    best_active = None
    best_score = -np.inf

    for tau in _candidate_thresholds(row_power, K, scan_size):
        if tau <= 0.0:
            continue
        pool = np.flatnonzero(row_power <= tau + max(1e-12, abs(tau) * 1e-12))
        if len(pool) < K:
            continue

        epsilon = np.sqrt(float(sigma) * float(tau) / float(P))
        active = _lazy_greedy_logdet_cap(V, pool, K, epsilon, eye)
        score = _true_log_general_score(row_grams, row_power, active, sigma, P, eye)
        if score > best_score:
            best_score = score
            best_active = active

    if best_active is None:
        best_active = np.zeros(N, dtype=bool)
        best_active[np.argsort(row_power)[-K:]] = True
    return best_active.astype(int)


def solve_cap_submodular_portfolio_gen(
    V,
    K,
    sigma=1.0,
    P=1.0,
    random_state=None,
    window_scan_size=DEFAULT_PORTFOLIO_WINDOW_SCAN_SIZE,
    refine_cap_count=DEFAULT_PORTFOLIO_REFINE_CAPS,
):
    """
    Fast cap-aware portfolio of monotone submodular fixed-cap objectives.

    For each fixed cap tau, the modular objective

        m_tau(A) = sum_{n in A} ||v_n||_2^2

    is a monotone submodular function. Its greedy maximizer under |A| = K is
    the top-power K-set within the cap, i.e. the cap-window candidate.  The best
    true-scoring modular caps are then refined with lazy greedy on the fixed-cap
    logdet submodular relaxation used by solve_cap_submodular_gen.
    """

    del random_state

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
    if sigma <= 0.0:
        raise ValueError("sigma must be positive for the cap-submodular relaxation.")
    if P <= 0.0:
        raise ValueError("P must be positive for the cap-submodular relaxation.")

    row_power = np.sum(np.abs(V) ** 2, axis=1).real
    row_grams = V.conj()[:, :, None] * V[:, None, :]
    eye = np.eye(L, dtype=complex)

    best_active = None
    best_score = -np.inf
    cap_scores = []
    seen = set()

    def consider(active, tau):
        nonlocal best_active, best_score
        key = np.packbits(active).tobytes()
        if key in seen:
            return None
        seen.add(key)

        score = _true_log_general_score(row_grams, row_power, active, sigma, P, eye)
        if score > best_score:
            best_score = score
            best_active = active
        if tau is not None:
            cap_scores.append((score, float(tau)))
        return score

    for tau in _candidate_thresholds(row_power, K, window_scan_size):
        active = _top_power_under_cap(row_power, K, tau)
        if active is not None:
            consider(active, tau)

    cap_scores.sort(key=lambda item: item[0], reverse=True)
    refine_caps = []
    seen_caps = set()
    for _, tau in cap_scores[: max(1, int(refine_cap_count))]:
        key = round(float(tau), 15)
        if key in seen_caps:
            continue
        seen_caps.add(key)
        refine_caps.append(float(tau))

    for tau in refine_caps:
        if tau <= 0.0:
            continue
        pool = np.flatnonzero(row_power <= tau + max(1e-12, abs(tau) * 1e-12))
        if len(pool) < K:
            continue

        epsilon = np.sqrt(float(sigma) * float(tau) / float(P))
        active = _lazy_greedy_logdet_cap(V, pool, K, epsilon, eye)
        consider(active, None)

    if best_active is None:
        best_active = np.zeros(N, dtype=bool)
        best_active[np.argsort(row_power)[-K:]] = True
    return best_active.astype(int)


def _candidate_thresholds(row_power, K, scan_size):
    N = len(row_power)
    max_drop = N - K
    if max_drop <= 0:
        return np.asarray([np.max(row_power)], dtype=float)

    sorted_power = np.sort(row_power)
    if scan_size is None:
        return np.unique(sorted_power[K - 1 :]).astype(float)

    off_count = max_drop
    weak_drop = off_count // 2
    strong_drop = off_count - weak_drop

    raw_drops = [
        0,
        1,
        2,
        3,
        5,
        8,
        10,
        15,
        20,
        25,
        40,
        50,
        75,
        100,
        N // 200,
        N // 100,
        N // 50,
        N // 25,
        N // 20,
        N // 10,
        strong_drop,
        max_drop,
    ]
    for radius in (3, 5, 10, 20, 30, 50, 75, 100):
        raw_drops.extend((strong_drop - radius, strong_drop + radius))

    grid_count = min(max_drop + 1, max(1, int(scan_size)))
    raw_drops.extend(np.linspace(0, max_drop, grid_count, dtype=int).tolist())

    thresholds = []
    seen = set()
    for drop in raw_drops:
        drop = int(np.clip(drop, 0, max_drop))
        threshold = float(sorted_power[N - drop - 1])
        key = round(threshold, 15)
        if key in seen:
            continue
        seen.add(key)
        thresholds.append(threshold)

    return np.asarray(sorted(thresholds), dtype=float)


def _greedy_logdet_cap(V, pool, K, epsilon, eye):
    N = V.shape[0]
    if len(pool) <= K:
        active = np.zeros(N, dtype=bool)
        active[np.asarray(pool[:K], dtype=int)] = True
        return active

    V_pool = V[pool]
    selected = np.zeros(len(pool), dtype=bool)
    inv_matrix = (1.0 / float(epsilon)) * eye.copy()

    for _ in range(K):
        leverage = np.real(
            np.einsum("ni,ij,nj->n", V_pool, inv_matrix, V_pool.conj(), optimize=True)
        )
        scores = np.log1p(np.maximum(leverage, 0.0))
        scores[selected] = -np.inf
        pos = int(np.argmax(scores))
        selected[pos] = True

        row = V_pool[pos]
        vec = row.conj()
        inv_vec = inv_matrix @ vec
        denom = max(1.0 + float(np.real(row @ inv_vec)), 1e-12)
        inv_matrix -= np.outer(inv_vec, inv_vec.conj()) / denom
        inv_matrix = 0.5 * (inv_matrix + inv_matrix.conj().T)

    active = np.zeros(N, dtype=bool)
    active[pool[selected]] = True
    return active


def _lazy_greedy_logdet_cap(V, pool, K, epsilon, eye):
    N = V.shape[0]
    if len(pool) <= K:
        active = np.zeros(N, dtype=bool)
        active[np.asarray(pool[:K], dtype=int)] = True
        return active

    V_pool = V[pool]
    selected = np.zeros(len(pool), dtype=bool)
    inv_matrix = (1.0 / float(epsilon)) * eye.copy()
    leverage = np.real(
        np.einsum("ni,ij,nj->n", V_pool, inv_matrix, V_pool.conj(), optimize=True)
    )
    gains = np.log1p(np.maximum(leverage, 0.0))
    heap = [(-float(gain), int(pos)) for pos, gain in enumerate(gains)]
    heapq.heapify(heap)

    for _ in range(K):
        while heap:
            _, pos = heapq.heappop(heap)
            if selected[pos]:
                continue

            row = V_pool[pos]
            gain = _single_logdet_gain(row, inv_matrix)
            while heap and selected[heap[0][1]]:
                heapq.heappop(heap)
            next_upper = -heap[0][0] if heap else -np.inf
            if gain >= next_upper - 1e-12:
                selected[pos] = True
                vec = row.conj()
                inv_vec = inv_matrix @ vec
                denom = max(1.0 + float(np.real(row @ inv_vec)), 1e-12)
                inv_matrix -= np.outer(inv_vec, inv_vec.conj()) / denom
                inv_matrix = 0.5 * (inv_matrix + inv_matrix.conj().T)
                break

            heapq.heappush(heap, (-float(gain), pos))

    active = np.zeros(N, dtype=bool)
    active[pool[selected]] = True
    return active


def _top_power_under_cap(row_power, K, tau):
    if tau <= 0.0:
        return None

    tolerance = max(1e-12, abs(float(tau)) * 1e-12)
    pool = np.flatnonzero(row_power <= float(tau) + tolerance)
    if len(pool) < K:
        return None

    if len(pool) == K:
        chosen = pool
    else:
        local = np.argpartition(row_power[pool], -K)[-K:]
        chosen = pool[local]

    active = np.zeros(len(row_power), dtype=bool)
    active[chosen] = True
    return active


def _single_logdet_gain(row, inv_matrix):
    vec = row.conj()
    leverage = float(np.real(row @ (inv_matrix @ vec)))
    return float(np.log1p(max(leverage, 0.0)))


def _true_log_general_score(row_grams, row_power, active, sigma, P, eye):
    if not np.any(active):
        return -np.inf

    max_power = float(np.max(row_power[active]))
    if max_power <= 0.0:
        return -np.inf

    gram = np.sum(row_grams[active], axis=0)
    gram_sq = gram @ gram.conj().T
    matrix = float(sigma) * eye + (float(P) / max_power) * gram_sq
    sign, logdet = np.linalg.slogdet(matrix)
    if sign <= 0.0 or not np.isfinite(logdet):
        return -np.inf
    return float(np.real(logdet))
