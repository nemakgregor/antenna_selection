import numpy as np


DEFAULT_EPS_VALUES = (1e-6, 1e-3, 1e-1)
DEFAULT_LAMBDAS = (1.0,)
DEFAULT_THRESHOLD_SCAN_SIZE = 9


def solve_thresholded_logdet_greedy(
    V,
    K,
    sigma=1.0,
    P=1.0,
    random_state=None,
    eps_values=DEFAULT_EPS_VALUES,
    lambdas=DEFAULT_LAMBDAS,
    threshold_scan_size=DEFAULT_THRESHOLD_SCAN_SIZE,
    include_h3_candidate=True,
    swap_max_passes=0,
    swap_remove_limit=None,
    swap_add_limit=None,
):
    """
    Thresholded D-optimal greedy for the general objective.

    For each upper row-power threshold, the solver greedily maximizes the
    monotone submodular surrogate

        lambda * log det(eps I + sum v_n^* v_n)
        + (1 - lambda) * sum ||v_n||^2

    over the threshold-feasible rows. It then scores each candidate by the
    true log det(sigma I + z(S)^2 G_S^2) objective and returns the best set.
    """

    del random_state

    context = _ThresholdLogdetContext(V, K, sigma=sigma, P=P)
    if context.K == 0:
        return np.zeros(context.N, dtype=int)
    if context.K == context.N:
        return np.ones(context.N, dtype=int)

    eps_values = _positive_tuple(eps_values, "eps_values")
    lambdas = _lambda_tuple(lambdas)
    thresholds = _candidate_thresholds(
        context.row_power,
        context.K,
        threshold_scan_size,
    )

    best_active = None
    best_score = -np.inf
    seen = set()

    def consider(indices):
        nonlocal best_active, best_score
        active = np.zeros(context.N, dtype=bool)
        active[np.asarray(indices, dtype=int)] = True
        if int(active.sum()) != context.K:
            active = _repair_active(context, active)

        key = np.packbits(active).tobytes()
        if key in seen:
            return
        seen.add(key)

        score = _true_log_general_score(context, active)
        if score > best_score:
            best_score = score
            best_active = active

    consider(_top_power_indices(context))
    if include_h3_candidate:
        consider(_h3_power_band_indices(context))

    for threshold in thresholds:
        tolerance = max(1e-12, abs(float(threshold)) * 1e-12)
        pool = np.flatnonzero(context.row_power <= float(threshold) + tolerance)
        if len(pool) < context.K:
            continue

        for eps in eps_values:
            for lam in lambdas:
                indices = _greedy_logdet_trace_indices(context, pool, eps, lam)
                consider(indices)

    if best_active is None:
        best_active = _repair_active(context, np.zeros(context.N, dtype=bool))

    if int(swap_max_passes) > 0:
        polished = _true_swap_polish(
            context,
            best_active,
            max_passes=swap_max_passes,
            remove_limit=swap_remove_limit,
            add_limit=swap_add_limit,
        )
        polished_score = _true_log_general_score(context, polished)
        if polished_score > best_score:
            best_active = polished

    return best_active.astype(int)


class _ThresholdLogdetContext:
    def __init__(self, V, K, sigma=1.0, P=1.0):
        V = np.asarray(V)
        if V.ndim != 2:
            raise ValueError("V must be a 2D complex matrix of shape (N, L).")
        if not np.iscomplexobj(V):
            V = V.astype(np.complex128)

        self.V = V
        self.N, self.L = V.shape
        self.K = int(K)
        if not (0 <= self.K <= self.N):
            raise ValueError("K must satisfy 0 <= K <= N.")

        self.sigma = float(sigma)
        self.P = float(P)
        self.eye = np.eye(self.L, dtype=complex)
        self.row_power = np.sum(np.abs(V) ** 2, axis=1).real
        self.row_grams = V.conj()[:, :, None] * V[:, None, :]


def _positive_tuple(values, name):
    result = tuple(float(value) for value in values)
    if not result or any(value <= 0.0 for value in result):
        raise ValueError(f"{name} must contain at least one positive value.")
    return result


def _lambda_tuple(values):
    result = tuple(float(value) for value in values)
    if not result or any(value < 0.0 or value > 1.0 for value in result):
        raise ValueError("lambdas must contain values in [0, 1].")
    return result


def _candidate_thresholds(row_power, K, scan_size):
    N = len(row_power)
    max_drop = N - K
    if max_drop <= 0:
        return np.asarray([np.max(row_power)], dtype=float)

    off_count = max_drop
    weak_drop = off_count // 2
    strong_drop = off_count - weak_drop

    raw = [
        0,
        1,
        2,
        5,
        10,
        25,
        50,
        100,
        N // 200,
        N // 100,
        N // 50,
        N // 20,
        N // 10,
        strong_drop,
        max_drop,
    ]
    for radius in (5, 10, 25, 50, 100):
        raw.extend((strong_drop - radius, strong_drop + radius))

    grid_count = min(max_drop + 1, max(1, int(scan_size)))
    raw.extend(np.linspace(0, max_drop, grid_count, dtype=int).tolist())

    sorted_power = np.sort(row_power)
    thresholds = []
    seen_drops = set()
    seen_thresholds = set()
    for drop in raw:
        drop = int(np.clip(drop, 0, max_drop))
        if drop in seen_drops:
            continue
        seen_drops.add(drop)

        threshold = float(sorted_power[N - drop - 1])
        key = round(threshold, 15)
        if key in seen_thresholds:
            continue
        seen_thresholds.add(key)
        thresholds.append(threshold)

    return np.asarray(sorted(thresholds), dtype=float)


def _greedy_logdet_trace_indices(context, pool, eps, lam):
    if len(pool) <= context.K:
        return np.asarray(pool[: context.K], dtype=int)

    V_pool = context.V[pool]
    row_power = context.row_power[pool]
    selected = np.zeros(len(pool), dtype=bool)
    inv_matrix = (1.0 / float(eps)) * context.eye.copy()

    for _ in range(context.K):
        leverage = _row_quadratic_scores(V_pool, inv_matrix)
        leverage = np.maximum(leverage, 0.0)

        if lam >= 1.0:
            scores = np.log1p(leverage)
        elif lam <= 0.0:
            scores = row_power.copy()
        else:
            scores = lam * np.log1p(leverage) + (1.0 - lam) * row_power

        scores[selected] = -np.inf
        pos = int(np.argmax(scores))
        selected[pos] = True

        row = V_pool[pos]
        vec = row.conj()
        inv_vec = inv_matrix @ vec
        denom = max(1.0 + float(np.real(row @ inv_vec)), 1e-12)
        inv_matrix = inv_matrix - np.outer(inv_vec, inv_vec.conj()) / denom
        inv_matrix = 0.5 * (inv_matrix + inv_matrix.conj().T)

    return np.asarray(pool[selected], dtype=int)


def _true_log_general_score(context, active):
    if not np.any(active):
        if context.sigma > 0.0:
            return float(context.L * np.log(context.sigma))
        return -np.inf

    max_power = float(np.max(context.row_power[active]))
    if max_power <= 0.0:
        if context.sigma > 0.0:
            return float(context.L * np.log(context.sigma))
        return -np.inf

    gram = np.sum(context.row_grams[active], axis=0)
    return _true_log_general_score_from_gram(context, gram, max_power)


def _true_log_general_score_from_gram(context, gram, max_power):
    if max_power <= 0.0:
        if context.sigma > 0.0:
            return float(context.L * np.log(context.sigma))
        return -np.inf

    gram_sq = gram @ gram.conj().T
    matrix = context.sigma * context.eye + (context.P / max_power) * gram_sq
    sign, logdet = np.linalg.slogdet(matrix)
    if sign <= 0.0 or not np.isfinite(logdet):
        return -np.inf
    return float(np.real(logdet))


def _true_log_general_scores_from_grams(context, grams, max_power):
    scores = np.full(len(grams), -np.inf, dtype=float)
    valid = max_power > 0.0
    if not np.any(valid):
        if context.sigma > 0.0:
            scores[~valid] = float(context.L * np.log(context.sigma))
        return scores

    grams_valid = grams[valid]
    z2 = context.P / max_power[valid]
    gram_sq = grams_valid @ np.swapaxes(grams_valid.conj(), 1, 2)
    matrices = z2[:, None, None] * gram_sq + context.sigma * context.eye[None, :, :]
    signs, logdets = np.linalg.slogdet(matrices)
    scores[valid] = np.where(signs > 0.0, np.real(logdets), -np.inf)
    return scores


def _true_swap_polish(context, active, max_passes, remove_limit, add_limit):
    active = np.asarray(active, dtype=bool).copy()
    gram = np.sum(context.row_grams[active], axis=0)
    max_power = float(np.max(context.row_power[active])) if np.any(active) else 0.0
    base_score = _true_log_general_score_from_gram(context, gram, max_power)

    for _ in range(max(0, int(max_passes))):
        rem_order, add_order = _swap_candidate_orders(
            context,
            active,
            gram,
            remove_limit,
            add_limit,
        )
        if len(rem_order) == 0 or len(add_order) == 0:
            break

        sorted_active = np.flatnonzero(active)
        sorted_active = sorted_active[np.argsort(context.row_power[sorted_active])]
        max1_i = int(sorted_active[-1])
        max1 = float(context.row_power[max1_i])
        max2 = (
            float(context.row_power[sorted_active[-2]])
            if len(sorted_active) >= 2
            else max1
        )

        best_i, best_j, best_score = _best_true_swap(
            context,
            gram,
            base_score,
            rem_order,
            add_order,
            max1_i,
            max1,
            max2,
        )
        if best_i is None:
            break

        active[best_i] = False
        active[best_j] = True
        gram = gram - context.row_grams[best_i] + context.row_grams[best_j]
        base_score = best_score

    return active


def _swap_candidate_orders(context, active, gram, remove_limit, add_limit):
    active_idx = np.flatnonzero(active)
    inactive_idx = np.flatnonzero(~active)
    if len(active_idx) == 0 or len(inactive_idx) == 0:
        return np.array([], dtype=int), np.array([], dtype=int)

    if remove_limit is None and add_limit is None:
        return active_idx, inactive_idx

    max_power = float(np.max(context.row_power[active]))
    grad = _true_logdet_gradient(context, gram, max_power)
    grad_scores = _row_quadratic_scores(context.V, grad)

    rem_order = _interleave_limited(
        [
            active_idx[np.argsort(grad_scores[active_idx])],
            active_idx[np.argsort(-context.row_power[active_idx])],
            active_idx[np.argsort(context.row_power[active_idx])],
        ],
        remove_limit,
    )
    add_order = _interleave_limited(
        [
            inactive_idx[np.argsort(-grad_scores[inactive_idx])],
            inactive_idx[np.argsort(context.row_power[inactive_idx])],
            inactive_idx[np.argsort(-context.row_power[inactive_idx])],
        ],
        add_limit,
    )
    return rem_order, add_order


def _best_true_swap(
    context,
    gram,
    base_score,
    rem_order,
    add_order,
    max1_i,
    max1,
    max2,
):
    best_score = base_score
    best_i = None
    best_j = None
    add_grams = context.row_grams[add_order]
    add_power = context.row_power[add_order]
    max_pairs_per_batch = _max_swap_batch_size(context.L)
    add_count = len(add_order)
    remove_chunk_size = max(1, max_pairs_per_batch // max(1, add_count))

    for start in range(0, len(rem_order), remove_chunk_size):
        rem_chunk = rem_order[start : start + remove_chunk_size]
        remove_grams = context.row_grams[rem_chunk]
        candidate_grams = (
            gram[None, None, :, :]
            - remove_grams[:, None, :, :]
            + add_grams[None, :, :, :]
        )
        max_without = np.where(rem_chunk == max1_i, max2, max1)
        max_power = np.maximum(max_without[:, None], add_power[None, :])
        scores = _true_log_general_scores_from_grams(
            context,
            candidate_grams.reshape(-1, context.L, context.L),
            max_power.ravel(),
        ).reshape(len(rem_chunk), add_count)

        pos = int(np.argmax(scores))
        score = float(scores.ravel()[pos])
        if score > best_score + 1e-12:
            rem_pos, add_pos = np.unravel_index(pos, scores.shape)
            best_score = score
            best_i = int(rem_chunk[rem_pos])
            best_j = int(add_order[add_pos])

    return best_i, best_j, best_score


def _max_swap_batch_size(L):
    return max(1024, min(25000, 2_000_000 // max(1, L * L)))


def _true_logdet_gradient(context, gram, max_power):
    if max_power <= 0.0:
        return np.zeros_like(gram)

    scale = context.P / max_power
    matrix = scale * (gram @ gram.conj().T) + context.sigma * context.eye
    try:
        gram_times_inv = np.linalg.solve(matrix.T, gram.T).T
    except np.linalg.LinAlgError:
        gram_times_inv = gram @ np.linalg.pinv(matrix)
    grad = 2.0 * scale * gram_times_inv
    return 0.5 * (grad + grad.conj().T)


def _top_power_indices(context):
    return np.argsort(context.row_power)[-context.K :]


def _h3_power_band_indices(context):
    off_count = context.N - context.K
    weak_drop = off_count // 2
    strong_drop = off_count - weak_drop
    power_order = np.argsort(context.row_power)

    active = np.ones(context.N, dtype=bool)
    if weak_drop:
        active[power_order[:weak_drop]] = False
    if strong_drop:
        active[power_order[context.N - strong_drop :]] = False
    return np.flatnonzero(active)


def _repair_active(context, active):
    active = np.asarray(active, dtype=bool).copy()
    count = int(active.sum())
    if count == context.K:
        return active

    if count > context.K:
        selected = np.flatnonzero(active)
        drop = selected[np.argsort(context.row_power[selected])[: count - context.K]]
        active[drop] = False
        return active

    missing = context.K - count
    inactive = np.flatnonzero(~active)
    add = inactive[np.argsort(context.row_power[inactive])[-missing:]]
    active[add] = True
    return active


def _row_quadratic_scores(V, matrix):
    return np.real(np.einsum("ni,ij,nj->n", V, matrix, V.conj(), optimize=True))


def _interleave_limited(order_lists, limit):
    if limit is None:
        limit = sum(len(values) for values in order_lists)
    limit = max(0, int(limit))
    result = []
    seen = set()
    max_len = max((len(values) for values in order_lists), default=0)
    for pos in range(max_len):
        for values in order_lists:
            if pos >= len(values):
                continue
            value = int(values[pos])
            if value in seen:
                continue
            seen.add(value)
            result.append(value)
            if len(result) >= limit:
                return np.asarray(result, dtype=int)
    return np.asarray(result, dtype=int)
