import numpy as np


def solve_h3_fast(
    V,
    K,
    *,
    mode="auto",
    use_refinement=True,
    refinement_iter=500,
    beam_size=100,
    eps=1e-10,
    random_state=None,
):
    """
    Fast antenna selection from heur.py.

    This heuristic always returns exactly K active antennas. For L=2 it uses
    angularly balanced phase selection; for L>2 it uses discrepancy-style
    dependent rounding, optionally followed by local swap refinement.
    """

    rng = np.random.default_rng(random_state)
    V = np.asarray(V)

    if V.ndim != 2:
        raise ValueError("V must be a 2D complex matrix of shape (N, L).")
    if not np.iscomplexobj(V):
        V = V.astype(np.complex128)

    N, L = V.shape

    if not (0 <= K <= N):
        raise ValueError("K must satisfy 0 <= K <= N.")
    if K == 0:
        return np.zeros(N, dtype=int)
    if K == N:
        return np.ones(N, dtype=int)

    if mode not in {"auto", "l2_phase", "general"}:
        raise ValueError("mode must be one of {'auto', 'l2_phase', 'general'}.")

    if mode == "auto":
        mode = "l2_phase" if L == 2 else "general"

    if mode == "l2_phase":
        if L != 2:
            raise ValueError("mode='l2_phase' requires V.shape[1] == 2.")

        x = _antenna_selection_l2_phase_fast(V, K, eps=eps)

        if use_refinement:
            A = _offdiag_features(V)
            x = _local_swap_refinement_fast(
                A,
                x,
                max_iter=refinement_iter,
                beam_size=beam_size,
                eps=eps,
            )

        return x

    A = _offdiag_features(V)
    x = _dependent_rounding_fast(A, K, eps=eps, rng=rng)

    if use_refinement:
        x = _local_swap_refinement_fast(
            A,
            x,
            max_iter=refinement_iter,
            beam_size=beam_size,
            eps=eps,
        )

    return x


def _antenna_selection_l2_phase_fast(V, K, eps=1e-10):
    N, L = V.shape
    if L != 2:
        raise ValueError("L2 phase algorithm requires L == 2.")

    g = np.conj(V[:, 0]) * V[:, 1]
    abs_g = np.abs(g)

    if np.all(abs_g <= eps):
        x = np.zeros(N, dtype=int)
        x[:K] = 1
        return x

    theta = np.angle(g)
    order = np.argsort(theta)
    selected_sorted_positions = _evenly_spaced_positions(N, K)
    selected_indices = order[selected_sorted_positions]

    x = np.zeros(N, dtype=int)
    x[selected_indices] = 1

    return _repair_cardinality_by_interference(g, x, K)


def _evenly_spaced_positions(N, K):
    if K == 0:
        return np.array([], dtype=int)
    if K == N:
        return np.arange(N, dtype=int)

    raw = (np.arange(K) + 0.5) * N / K
    pos = np.floor(raw).astype(int)
    pos = np.clip(pos, 0, N - 1)

    used = set()
    result = []

    for p in pos:
        if p not in used:
            used.add(p)
            result.append(p)
            continue

        left = p - 1
        right = p + 1
        while True:
            if left >= 0 and left not in used:
                used.add(left)
                result.append(left)
                break
            if right < N and right not in used:
                used.add(right)
                result.append(right)
                break
            left -= 1
            right += 1

    return np.array(result, dtype=int)


def _repair_cardinality_by_interference(g, x, K):
    x = x.copy()
    current = np.sum(g[x == 1])

    while x.sum() < K:
        inactive = np.where(x == 0)[0]
        vals = np.abs(current + g[inactive]) ** 2
        j = inactive[np.argmin(vals)]
        x[j] = 1
        current += g[j]

    while x.sum() > K:
        active = np.where(x == 1)[0]
        vals = np.abs(current - g[active]) ** 2
        i = active[np.argmin(vals)]
        x[i] = 0
        current -= g[i]

    return x


def _offdiag_features(V):
    N, L = V.shape

    if L < 2:
        return np.zeros((N, 0), dtype=float)

    features = []

    for p in range(L):
        for q in range(p + 1, L):
            g = np.conj(V[:, p]) * V[:, q]
            features.append(np.real(g))
            features.append(np.imag(g))

    return np.column_stack(features).astype(float)


def _dependent_rounding_fast(A, K, eps=1e-10, rng=None):
    if rng is None:
        rng = np.random.default_rng()

    N, d = A.shape

    if K == 0:
        return np.zeros(N, dtype=int)
    if K == N:
        return np.ones(N, dtype=int)

    if d == 0:
        x = np.zeros(N, dtype=int)
        x[:K] = 1
        return x

    rho = K / N
    y = np.full(N, rho, dtype=float)

    while True:
        F = np.where((y > eps) & (y < 1.0 - eps))[0]

        if len(F) <= d + 1:
            break

        M = np.vstack([np.ones(len(F)), A[F].T])
        alpha = _random_null_vector_projected(M, eps=eps, rng=rng)

        if alpha is None:
            alpha = _null_vector_svd(M, eps=eps)

        if alpha is None:
            break

        t_min, t_max = _max_feasible_interval(y[F], alpha, eps=eps)

        if not np.isfinite(t_min) or not np.isfinite(t_max):
            break

        t = t_max if abs(t_max) >= abs(t_min) else t_min

        if abs(t) <= eps:
            break

        y[F] += t * alpha
        y[y < eps] = 0.0
        y[y > 1.0 - eps] = 1.0

    return _final_rounding_greedy(A, y, K, eps=eps)


def _random_null_vector_projected(M, eps=1e-10, rng=None, max_tries=8):
    if rng is None:
        rng = np.random.default_rng()

    r_dim, c_dim = M.shape

    if c_dim <= r_dim:
        return None

    with np.errstate(all="ignore"):
        G = M @ M.T
        scale = max(1.0, np.linalg.norm(G))
        G_reg = G + eps * scale * np.eye(r_dim)
        M_norm = max(1.0, np.linalg.norm(M))

    if not np.isfinite(G_reg).all() or not np.isfinite(M_norm):
        return None

    for _ in range(max_tries):
        r = rng.standard_normal(c_dim)
        with np.errstate(all="ignore"):
            Mr = M @ r

        if not np.isfinite(Mr).all():
            return None

        try:
            coeff = np.linalg.solve(G_reg, Mr)
        except np.linalg.LinAlgError:
            return None

        with np.errstate(all="ignore"):
            alpha = r - M.T @ coeff
            norm = np.linalg.norm(alpha)

        if not np.isfinite(norm) or norm <= eps:
            continue

        alpha /= norm
        with np.errstate(all="ignore"):
            residual = np.linalg.norm(M @ alpha)

        if np.isfinite(residual) and residual <= 1e-6 * M_norm:
            return alpha

    return None


def _null_vector_svd(M, eps=1e-10):
    try:
        _, _, vh = np.linalg.svd(M, full_matrices=True)
    except np.linalg.LinAlgError:
        return None

    alpha = vh[-1, :]
    norm = np.linalg.norm(alpha)

    if norm <= eps:
        return None

    alpha = alpha / norm

    with np.errstate(all="ignore"):
        residual = np.linalg.norm(M @ alpha)
        scale = max(1.0, np.linalg.norm(M))

    if not np.isfinite(residual) or residual > 1e-6 * scale:
        return None

    return alpha


def _max_feasible_interval(yF, alpha, eps=1e-10):
    t_min = -np.inf
    t_max = np.inf

    for y_i, a_i in zip(yF, alpha):
        if abs(a_i) <= eps:
            continue

        if a_i > 0:
            t_max = min(t_max, (1.0 - y_i) / a_i)
            t_min = max(t_min, -y_i / a_i)
        else:
            t_max = min(t_max, -y_i / a_i)
            t_min = max(t_min, (1.0 - y_i) / a_i)

    return t_min, t_max


def _final_rounding_greedy(A, y, K, eps=1e-10):
    N, _ = A.shape
    x = np.zeros(N, dtype=int)

    ones = np.where(y >= 1.0 - eps)[0]
    frac = np.where((y > eps) & (y < 1.0 - eps))[0]

    x[ones] = 1
    need = int(K - x.sum())

    if need == 0:
        return x

    if need < 0:
        with np.errstate(all="ignore"):
            current = A.T @ x
        for _ in range(-need):
            active = np.where(x == 1)[0]
            vals = np.empty(len(active))

            for t, i in enumerate(active):
                with np.errstate(all="ignore"):
                    cand = current - A[i]
                    vals[t] = cand @ cand

            i = active[np.argmin(vals)]
            x[i] = 0
            current -= A[i]

        return x

    if need > len(frac):
        frac = np.where(x == 0)[0]

    with np.errstate(all="ignore"):
        current = A.T @ x
    candidates = frac.copy()

    for _ in range(need):
        if len(candidates) == 0:
            break

        with np.errstate(all="ignore"):
            vals = np.sum((A[candidates] + current[None, :]) ** 2, axis=1)
        pos = int(np.argmin(vals))
        j = candidates[pos]

        x[j] = 1
        current += A[j]

        candidates = np.delete(candidates, pos)

    if x.sum() != K:
        x = _repair_binary_cardinality_by_features(A, x, K)

    return x


def _repair_binary_cardinality_by_features(A, x, K):
    x = x.copy()
    with np.errstate(all="ignore"):
        current = A.T @ x

    while x.sum() < K:
        inactive = np.where(x == 0)[0]
        with np.errstate(all="ignore"):
            vals = np.sum((A[inactive] + current[None, :]) ** 2, axis=1)
        j = inactive[np.argmin(vals)]
        x[j] = 1
        current += A[j]

    while x.sum() > K:
        active = np.where(x == 1)[0]
        with np.errstate(all="ignore"):
            vals = np.sum((current[None, :] - A[active]) ** 2, axis=1)
        i = active[np.argmin(vals)]
        x[i] = 0
        current -= A[i]

    return x


def _local_swap_refinement_fast(A, x, max_iter=500, beam_size=100, eps=1e-10):
    x = x.copy().astype(int)
    _, d = A.shape

    if d == 0:
        return x

    beam_size = max(1, int(beam_size))
    active_mask = x.astype(bool)

    with np.errstate(all="ignore"):
        current = A.T @ x
        current_val = current @ current

    for _ in range(max_iter):
        with np.errstate(all="ignore"):
            scores = A @ current

        active_idx = np.where(active_mask)[0]
        inactive_idx = np.where(~active_mask)[0]

        if len(active_idx) == 0 or len(inactive_idx) == 0:
            break

        b_rem = min(beam_size, len(active_idx))
        b_add = min(beam_size, len(inactive_idx))

        if b_rem < len(active_idx):
            rem_part = np.argpartition(-scores[active_idx], b_rem - 1)[:b_rem]
            rem_candidates = active_idx[rem_part]
        else:
            rem_candidates = active_idx

        if b_add < len(inactive_idx):
            add_part = np.argpartition(scores[inactive_idx], b_add - 1)[:b_add]
            add_candidates = inactive_idx[add_part]
        else:
            add_candidates = inactive_idx

        best_i = None
        best_j = None
        best_val = current_val
        best_current = None

        for i in rem_candidates:
            base = current - A[i]
            with np.errstate(all="ignore"):
                candidates = base[None, :] + A[add_candidates]
                vals = np.sum(candidates * candidates, axis=1)

            pos = int(np.argmin(vals))
            val = vals[pos]

            if val + eps < best_val:
                best_val = val
                best_i = i
                best_j = add_candidates[pos]
                best_current = candidates[pos].copy()

        if best_i is None:
            break

        active_mask[best_i] = False
        active_mask[best_j] = True

        x[best_i] = 0
        x[best_j] = 1

        current = best_current
        current_val = best_val

    return x
