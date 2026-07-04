import numpy as np


DEFAULT_MAX_PASSES = 3
DEFAULT_REMOVE_LIMIT = 64
DEFAULT_ADD_LIMIT = 128
DEFAULT_BOUNDARY_LIMIT = 128


def refine_general_1swap(
    V,
    x,
    K=None,
    sigma=1.0,
    P=1.0,
    *,
    max_passes=DEFAULT_MAX_PASSES,
    remove_limit=DEFAULT_REMOVE_LIMIT,
    add_limit=DEFAULT_ADD_LIMIT,
    boundary_limit=DEFAULT_BOUNDARY_LIMIT,
    min_gain=1e-12,
):
    """
    Cap-aware one-swap hill climbing for the general objective.

    The refinement scores the true log U_G objective for every considered
    remove/add pair, but keeps N=1000 runs fast by restricting the neighborhood
    to high-signal candidates: low marginal active rows, cap-heavy active rows,
    high marginal inactive rows, and rows near the current power-band boundary.
    """

    context = _GeneralSwapContext(V, int(np.sum(x)) if K is None else K, sigma, P)
    if context.K == 0:
        return np.zeros(context.N, dtype=int)
    if context.K == context.N:
        return np.ones(context.N, dtype=int)

    active = _repair_active(context, np.asarray(x, dtype=bool))
    gram = _build_gram(context, active)
    max_power = _active_max_power(context, active)
    base_score = _log_general_score_from_gram(context, gram, max_power)

    for _ in range(max(0, int(max_passes))):
        rem_order, add_order = _candidate_orders(
            context,
            active,
            gram,
            remove_limit=remove_limit,
            add_limit=add_limit,
            boundary_limit=boundary_limit,
        )
        if len(rem_order) == 0 or len(add_order) == 0:
            break

        max1_i, max1, max2 = _top_active_powers(context, active)
        best_i, best_j, best_score = _best_swap(
            context,
            gram,
            base_score,
            rem_order,
            add_order,
            max1_i,
            max1,
            max2,
            min_gain,
        )
        if best_i is None:
            break

        active[best_i] = False
        active[best_j] = True
        gram = gram - context.row_grams[best_i] + context.row_grams[best_j]
        base_score = best_score

    return active.astype(int)


class _GeneralSwapContext:
    def __init__(self, V, K, sigma, P):
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
        self.power_order = np.argsort(self.row_power)
        self.power_rank = np.empty(self.N, dtype=int)
        self.power_rank[self.power_order] = np.arange(self.N)


def _candidate_orders(
    context,
    active,
    gram,
    remove_limit,
    add_limit,
    boundary_limit,
):
    active_idx = np.flatnonzero(active)
    inactive_idx = np.flatnonzero(~active)
    if len(active_idx) == 0 or len(inactive_idx) == 0:
        return np.array([], dtype=int), np.array([], dtype=int)

    max_power = _active_max_power(context, active)
    grad = _logdet_gradient(context, gram, max_power)
    grad_scores = np.nan_to_num(
        _row_quadratic_scores(context.V, grad),
        nan=-np.inf,
        posinf=np.inf,
        neginf=-np.inf,
    )

    weak_align, strong_align = _principal_alignments(context, gram)
    boundary_rem, boundary_add = _power_boundary_orders(
        context,
        active,
        boundary_limit,
    )

    rem_order = _interleave_limited(
        [
            active_idx[np.argsort(grad_scores[active_idx])],
            active_idx[np.argsort(-context.row_power[active_idx])],
            active_idx[np.argsort(context.row_power[active_idx])],
            active_idx[
                np.argsort(-(strong_align[active_idx] - weak_align[active_idx]))
            ],
            boundary_rem,
        ],
        remove_limit,
    )
    add_order = _interleave_limited(
        [
            inactive_idx[np.argsort(-grad_scores[inactive_idx])],
            inactive_idx[np.argsort(-weak_align[inactive_idx])],
            inactive_idx[np.argsort(context.row_power[inactive_idx])],
            inactive_idx[np.argsort(-context.row_power[inactive_idx])],
            boundary_add,
        ],
        add_limit,
    )
    return rem_order, add_order


def _best_swap(
    context,
    gram,
    base_score,
    rem_order,
    add_order,
    max1_i,
    max1,
    max2,
    min_gain,
):
    best_score = base_score
    best_i = None
    best_j = None
    add_grams = context.row_grams[add_order]
    add_power = context.row_power[add_order]
    add_count = len(add_order)
    max_pairs_per_batch = _max_swap_batch_size(context.L)
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
        scores = _log_general_scores_from_grams(
            context,
            candidate_grams.reshape(-1, context.L, context.L),
            max_power.ravel(),
        ).reshape(len(rem_chunk), add_count)

        pos = int(np.argmax(scores))
        score = float(scores.ravel()[pos])
        if score > best_score + float(min_gain):
            rem_pos, add_pos = np.unravel_index(pos, scores.shape)
            best_score = score
            best_i = int(rem_chunk[rem_pos])
            best_j = int(add_order[add_pos])

    return best_i, best_j, best_score


def _log_general_score_from_gram(context, gram, max_power):
    return float(
        _log_general_scores_from_grams(
            context,
            gram[None, :, :],
            np.asarray([max_power], dtype=float),
        )[0]
    )


def _log_general_scores_from_grams(context, grams, max_power):
    scores = np.full(len(grams), -np.inf, dtype=float)
    valid = max_power > 0.0
    if np.any(~valid) and context.sigma > 0.0:
        scores[~valid] = float(context.L * np.log(context.sigma))
    if not np.any(valid):
        return scores

    if context.L == 2:
        scores[valid] = _log_general_scores_l2(
            context,
            grams[valid],
            max_power[valid],
        )
        return scores

    grams_valid = grams[valid]
    z2 = context.P / max_power[valid]
    gram_sq = grams_valid @ np.swapaxes(grams_valid.conj(), 1, 2)
    matrices = z2[:, None, None] * gram_sq + context.sigma * context.eye[None, :, :]
    signs, logdets = np.linalg.slogdet(matrices)
    scores[valid] = np.where(signs > 0.0, np.real(logdets), -np.inf)
    return scores


def _log_general_scores_l2(context, grams, max_power):
    scores = np.full(len(grams), -np.inf, dtype=float)
    z2 = context.P / max_power

    a = np.real(grams[:, 0, 0])
    d = np.real(grams[:, 1, 1])
    b_abs_sq = np.abs(grams[:, 0, 1]) ** 2

    sq00 = a * a + b_abs_sq
    sq11 = d * d + b_abs_sq
    offdiag_abs_sq = b_abs_sq * (a + d) ** 2

    det = (
        (context.sigma + z2 * sq00) * (context.sigma + z2 * sq11)
        - (z2**2) * offdiag_abs_sq
    )
    good = (det > 0.0) & np.isfinite(det)
    scores[good] = np.log(det[good])
    return scores


def _logdet_gradient(context, gram, max_power):
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


def _principal_alignments(context, gram):
    if context.L == 1:
        align = np.abs(context.V[:, 0]) ** 2
        return align, align

    hermitian = 0.5 * (gram + gram.conj().T)
    eigvals, eigvecs = np.linalg.eigh(hermitian)
    weak_vec = eigvecs[:, int(np.argmin(eigvals))]
    strong_vec = eigvecs[:, int(np.argmax(eigvals))]
    weak_align = np.abs(context.V @ weak_vec) ** 2
    strong_align = np.abs(context.V @ strong_vec) ** 2
    return weak_align, strong_align


def _power_boundary_orders(context, active, boundary_limit):
    active_idx = np.flatnonzero(active)
    inactive_idx = np.flatnonzero(~active)
    if len(active_idx) == 0 or len(inactive_idx) == 0:
        return np.array([], dtype=int), np.array([], dtype=int)

    active_ranks = context.power_rank[active_idx]
    low_rank = int(np.min(active_ranks))
    high_rank = int(np.max(active_ranks))

    active_low = active_idx[np.argsort(active_ranks)]
    active_high = active_idx[np.argsort(-active_ranks)]
    remove_order = _interleave_limited([active_low, active_high], boundary_limit)

    inactive_ranks = context.power_rank[inactive_idx]
    distance = np.minimum(
        np.abs(inactive_ranks - low_rank),
        np.abs(inactive_ranks - high_rank),
    )
    order = np.lexsort((-context.row_power[inactive_idx], distance))
    add_order = inactive_idx[order]
    if boundary_limit is not None:
        add_order = add_order[: max(0, int(boundary_limit))]
    return remove_order, add_order


def _top_active_powers(context, active):
    active_idx = np.flatnonzero(active)
    if len(active_idx) == 0:
        return -1, 0.0, 0.0
    sorted_active = active_idx[np.argsort(context.row_power[active_idx])]
    max1_i = int(sorted_active[-1])
    max1 = float(context.row_power[max1_i])
    max2 = float(context.row_power[sorted_active[-2]]) if len(sorted_active) >= 2 else 0.0
    return max1_i, max1, max2


def _active_max_power(context, active):
    if not np.any(active):
        return 0.0
    return float(np.max(context.row_power[active]))


def _build_gram(context, active):
    if not np.any(active):
        return np.zeros((context.L, context.L), dtype=complex)
    return np.sum(context.row_grams[active], axis=0)


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


def _max_swap_batch_size(L):
    return max(1024, min(25000, 2_000_000 // max(1, L * L)))


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
