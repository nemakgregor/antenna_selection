import numpy as np

from .h3_fast import solve_h3_fast
from .h3_threshold import solve_h3


DEFAULT_LAMBDAS = (1e-2, 1e-1, 1.0)


def _solve_frame_variant(
    V,
    K,
    target_obj,
    sigma,
    P,
    random_state,
    external_starts,
    kwargs,
):
    kwargs = dict(kwargs)
    if "external_starts" in kwargs and kwargs["external_starts"] != external_starts:
        raise ValueError("external_starts is fixed by the selected Frame wrapper.")
    kwargs["external_starts"] = external_starts
    return solve_frame_portfolio(
        V,
        K,
        target_obj=target_obj,
        sigma=sigma,
        P=P,
        random_state=random_state,
        **kwargs,
    )


def solve_frame_general(
    V,
    K,
    sigma=1.0,
    P=1.0,
    random_state=None,
    **kwargs,
):
    """Frame-aware portfolio tuned for the exact general determinant objective."""

    return _solve_frame_variant(
        V,
        K,
        "gen",
        sigma,
        P,
        random_state,
        True,
        kwargs,
    )


def solve_frame_bf(
    V,
    K,
    sigma=1.0,
    P=1.0,
    random_state=None,
    **kwargs,
):
    """Frame-aware portfolio tuned for exact BF gain."""

    return _solve_frame_variant(
        V,
        K,
        "bf",
        sigma,
        P,
        random_state,
        True,
        kwargs,
    )


def solve_frame_interference(
    V,
    K,
    sigma=1.0,
    P=1.0,
    random_state=None,
    **kwargs,
):
    """Frame-aware portfolio tuned for exact interference minimization."""

    return _solve_frame_variant(
        V,
        K,
        "int",
        sigma,
        P,
        random_state,
        True,
        kwargs,
    )


def solve_frame_only_general(
    V,
    K,
    sigma=1.0,
    P=1.0,
    random_state=None,
    **kwargs,
):
    """Self-contained frame algorithm tuned for the exact general objective."""

    return _solve_frame_variant(
        V,
        K,
        "gen",
        sigma,
        P,
        random_state,
        False,
        kwargs,
    )


def solve_frame_only_bf(
    V,
    K,
    sigma=1.0,
    P=1.0,
    random_state=None,
    **kwargs,
):
    """Self-contained frame algorithm tuned for exact BF gain."""

    return _solve_frame_variant(
        V,
        K,
        "bf",
        sigma,
        P,
        random_state,
        False,
        kwargs,
    )


def solve_frame_only_interference(
    V,
    K,
    sigma=1.0,
    P=1.0,
    random_state=None,
    **kwargs,
):
    """Self-contained frame algorithm tuned for exact interference minimization."""

    return _solve_frame_variant(
        V,
        K,
        "int",
        sigma,
        P,
        random_state,
        False,
        kwargs,
    )


def solve_frame_portfolio(
    V,
    K,
    *,
    target_obj="gen",
    sigma=1.0,
    P=1.0,
    random_state=None,
    max_refined_starts=6,
    max_passes=4,
    remove_limit=120,
    add_limit=120,
    lambdas=DEFAULT_LAMBDAS,
    random_restarts=0,
    use_h3_fast_int_start=True,
    h3_fast_refinement_iter=None,
    h3_fast_beam_size=None,
    external_starts=True,
):
    """
    Multi-start frame subset selection.

    The starts are cap-aware: they deliberately try prefixes after dropping the
    strongest rows, because the exact objectives scale by z(S), and z(S) is set
    by the largest selected row power. The selected starts are then improved by
    exact vectorized one-swap local search for the requested target objective.
    """

    if target_obj not in {"bf", "int", "gen"}:
        raise ValueError("target_obj must be one of {'bf', 'int', 'gen'}.")

    context = _FrameContext(V, K, sigma=sigma, P=P)
    if context.K == 0:
        return np.zeros(context.N, dtype=int)
    if context.K == context.N:
        return np.ones(context.N, dtype=int)

    rng = np.random.default_rng(0 if random_state is None else int(random_state))
    starts = _build_starts(
        context,
        target_obj=target_obj,
        rng=rng,
        lambdas=lambdas,
        random_restarts=random_restarts,
        use_h3_fast_int_start=use_h3_fast_int_start,
        h3_fast_refinement_iter=h3_fast_refinement_iter,
        h3_fast_beam_size=h3_fast_beam_size,
        external_starts=external_starts,
    )

    scored_starts = []
    seen = set()
    for active in starts:
        active = _repair_active(context, active)
        key = np.packbits(active).tobytes()
        if key in seen:
            continue
        seen.add(key)
        score = _score_active(context, active, target_obj)
        if np.isfinite(score):
            scored_starts.append((score, active))

    if not scored_starts:
        return _active_to_x(_top_power_start(context))

    scored_starts.sort(key=lambda item: item[0], reverse=True)
    best_score, best_active = scored_starts[0][0], scored_starts[0][1].copy()

    for _, active in scored_starts[: max(1, int(max_refined_starts))]:
        improved = _swap_local_search(
            context,
            active,
            target_obj=target_obj,
            max_passes=max_passes,
            remove_limit=remove_limit,
            add_limit=add_limit,
        )
        score = _score_active(context, improved, target_obj)
        if score > best_score:
            best_score = score
            best_active = improved

    return _active_to_x(best_active)


class _FrameContext:
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
        self.sorted_power = np.argsort(self.row_power)[::-1]
        self._frame_features = None
        self._offdiag_features = None

    @property
    def frame_features(self):
        if self._frame_features is None:
            self._frame_features = _hermitian_features(self.row_grams, include_diag=True)
        return self._frame_features

    @property
    def offdiag_features(self):
        if self._offdiag_features is None:
            self._offdiag_features = _hermitian_features(
                self.row_grams,
                include_diag=False,
            )
        return self._offdiag_features


def _build_starts(
    context,
    target_obj,
    rng,
    lambdas,
    random_restarts,
    use_h3_fast_int_start,
    h3_fast_refinement_iter,
    h3_fast_beam_size,
    external_starts,
):
    starts = []
    sorted_idx = context.sorted_power
    offsets = _cap_offsets(context.N, context.K)
    buffers = _buffer_sizes(context.N, context.K)

    starts.append(_top_power_start(context))

    for offset in offsets:
        if offset + context.K <= context.N:
            starts.append(_indices_to_active(context, sorted_idx[offset : offset + context.K]))

    for offset in offsets[:5]:
        for buffer in buffers:
            if buffer <= 0:
                continue
            starts.append(_frame_pool_start(context, offset, buffer))
            if target_obj == "int" and not external_starts:
                starts.append(_offdiag_pool_start(context, offset, buffer))
            if target_obj in {"bf", "gen"}:
                for lam in lambdas:
                    starts.append(_dopt_pool_start(context, offset, buffer, lam))

    if target_obj == "int":
        if not external_starts:
            for offset in offsets[:4]:
                starts.append(_offdiag_greedy_add_start(context, offset))
            starts.append(_phase_balanced_start(context))

        if external_starts:
            starts.append(
                solve_h3(
                    context.V,
                    context.K,
                    target_obj="int",
                    sigma=context.sigma,
                    P=context.P,
                ).astype(bool)
            )
        if external_starts and use_h3_fast_int_start:
            if h3_fast_refinement_iter is None:
                refinement_iter = 250 if context.N <= 2000 else 120
            else:
                refinement_iter = int(h3_fast_refinement_iter)
            if h3_fast_beam_size is None:
                beam_size = 80 if context.N <= 2000 else 50
            else:
                beam_size = int(h3_fast_beam_size)
            starts.append(
                solve_h3_fast(
                    context.V,
                    context.K,
                    refinement_iter=refinement_iter,
                    beam_size=beam_size,
                    random_state=int(rng.integers(0, 2**31 - 1)),
                ).astype(bool)
            )

    for offset in offsets[:4]:
        pool = sorted_idx[offset:]
        if len(pool) < context.K:
            continue
        for _ in range(int(random_restarts)):
            starts.append(_indices_to_active(context, rng.choice(pool, context.K, replace=False)))

    return starts


def _cap_offsets(N, K):
    max_offset = N - K
    raw = [0, 1, 2, 5, 10, 25, 50, N // 20, N // 10]
    if max_offset >= 100:
        raw.append(100)
    if max_offset >= 150 and N <= 2000:
        raw.append(150)
    result = []
    seen = set()
    for value in raw:
        value = int(value)
        if 0 <= value <= max_offset and value not in seen:
            seen.add(value)
            result.append(value)
    return result


def _buffer_sizes(N, K):
    slack = N - K
    raw = [20, 50, 100] if N <= 2000 else [20, 50]
    result = []
    seen = set()
    for value in raw:
        value = min(int(value), slack)
        if value > 0 and value not in seen:
            seen.add(value)
            result.append(value)
    return result


def _top_power_start(context):
    return _indices_to_active(context, context.sorted_power[: context.K])


def _frame_pool_start(context, offset, buffer):
    pool = context.sorted_power[offset : offset + context.K + buffer]
    if len(pool) <= context.K:
        return _indices_to_active(context, pool[: context.K])

    drop_count = len(pool) - context.K
    features = context.frame_features[pool]
    norms = np.sum(features * features, axis=1)
    residual = (drop_count / len(pool)) * np.sum(features, axis=0)
    available = np.ones(len(pool), dtype=bool)

    for _ in range(drop_count):
        scores = features @ residual - 0.5 * norms
        scores[~available] = -np.inf
        pos = int(np.argmax(scores))
        if not np.isfinite(scores[pos]):
            break
        available[pos] = False
        residual -= features[pos]

    return _indices_to_active(context, pool[available])


def _offdiag_pool_start(context, offset, buffer):
    pool = context.sorted_power[offset : offset + context.K + buffer]
    if len(pool) <= context.K:
        return _indices_to_active(context, pool[: context.K])

    features = context.offdiag_features[pool]
    if features.shape[1] == 0:
        return _indices_to_active(context, pool[: context.K])

    drop_count = len(pool) - context.K
    norms = np.sum(features * features, axis=1)
    current = np.sum(features, axis=0)
    available = np.ones(len(pool), dtype=bool)

    for _ in range(drop_count):
        scores = features @ current - 0.5 * norms
        scores[~available] = -np.inf
        pos = int(np.argmax(scores))
        if not np.isfinite(scores[pos]):
            break
        available[pos] = False
        current -= features[pos]

    return _indices_to_active(context, pool[available])


def _offdiag_greedy_add_start(context, offset):
    pool = context.sorted_power[offset:]
    if len(pool) <= context.K:
        return _indices_to_active(context, pool[: context.K])

    features = context.offdiag_features[pool]
    if features.shape[1] == 0:
        return _indices_to_active(context, pool[: context.K])

    selected = np.zeros(len(pool), dtype=bool)
    current = np.zeros(features.shape[1], dtype=float)
    norms = np.sum(features * features, axis=1)

    for _ in range(context.K):
        scores = features @ current + 0.5 * norms
        scores[selected] = np.inf
        pos = int(np.argmin(scores))
        if not np.isfinite(scores[pos]):
            break
        selected[pos] = True
        current += features[pos]

    return _indices_to_active(context, pool[selected])


def _phase_balanced_start(context):
    if context.L != 2:
        return _top_power_start(context)

    g = context.V[:, 0].conj() * context.V[:, 1]
    if np.all(np.abs(g) <= 1e-12):
        return _top_power_start(context)

    order = np.argsort(np.angle(g))
    positions = np.floor((np.arange(context.K) + 0.5) * context.N / context.K).astype(int)
    positions = np.clip(positions, 0, context.N - 1)
    selected = order[_unique_nearby_positions(positions, context.N)]
    active = _indices_to_active(context, selected)
    return _repair_active_by_offdiag(context, active)


def _unique_nearby_positions(positions, N):
    used = set()
    result = []
    for pos in positions:
        pos = int(pos)
        if pos not in used:
            used.add(pos)
            result.append(pos)
            continue

        left = pos - 1
        right = pos + 1
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

    return np.asarray(result, dtype=int)


def _repair_active_by_offdiag(context, active):
    active = _repair_active(context, active)
    features = context.offdiag_features
    if features.shape[1] == 0:
        return active

    current = np.sum(features[active], axis=0)
    while active.sum() < context.K:
        inactive_idx = np.flatnonzero(~active)
        scores = np.sum((current[None, :] + features[inactive_idx]) ** 2, axis=1)
        antenna = int(inactive_idx[np.argmin(scores)])
        active[antenna] = True
        current += features[antenna]

    while active.sum() > context.K:
        active_idx = np.flatnonzero(active)
        scores = np.sum((current[None, :] - features[active_idx]) ** 2, axis=1)
        antenna = int(active_idx[np.argmin(scores)])
        active[antenna] = False
        current -= features[antenna]

    return active


def _dopt_pool_start(context, offset, buffer, lam):
    pool = context.sorted_power[offset : offset + context.K + buffer]
    if len(pool) <= context.K:
        return _indices_to_active(context, pool[: context.K])

    active_local = np.ones(len(pool), dtype=bool)
    gram = np.sum(context.row_grams[pool], axis=0)
    matrix = gram + float(lam) * context.eye
    try:
        inv_matrix = np.linalg.inv(matrix)
    except np.linalg.LinAlgError:
        inv_matrix = np.linalg.pinv(matrix)

    for _ in range(len(pool) - context.K):
        live = np.flatnonzero(active_local)
        scores = _row_quadratic_scores(context.V[pool[live]], inv_matrix)
        pos = int(live[np.argmin(scores)])
        antenna = int(pool[pos])
        a = context.V[antenna].conj()
        inv_a = inv_matrix @ a
        leverage = float(np.real(context.V[antenna] @ inv_a))
        denom = max(1.0 - leverage, 1e-10)
        inv_matrix = inv_matrix + np.outer(inv_a, inv_a.conj()) / denom
        active_local[pos] = False

    return _indices_to_active(context, pool[active_local])


def _swap_local_search(
    context,
    active,
    target_obj,
    max_passes,
    remove_limit,
    add_limit,
):
    active = _repair_active(context, active)
    gram = _build_gram(context, active)
    base_score = _score_from_gram(context, gram, np.max(context.row_power[active]), target_obj)

    for _ in range(max(0, int(max_passes))):
        rem_order, add_order = _candidate_orders(
            context,
            active,
            gram,
            target_obj=target_obj,
            remove_limit=remove_limit,
            add_limit=add_limit,
        )
        if len(rem_order) == 0 or len(add_order) == 0:
            break

        active_idx = np.flatnonzero(active)
        sorted_active = active_idx[np.argsort(context.row_power[active_idx])]
        max1_i = int(sorted_active[-1])
        max1 = float(context.row_power[max1_i])
        max2 = (
            float(context.row_power[sorted_active[-2]])
            if len(sorted_active) >= 2
            else max1
        )

        best_i, best_j, best_score = _best_swap(
            context,
            gram,
            base_score,
            rem_order,
            add_order,
            target_obj,
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


def _best_swap(
    context,
    gram,
    base_score,
    rem_order,
    add_order,
    target_obj,
    max1_i,
    max1,
    max2,
):
    if target_obj == "bf":
        return _best_bf_swap(
            context,
            gram,
            base_score,
            rem_order,
            add_order,
            max1_i,
            max1,
            max2,
        )

    if not _use_batched_matrix_swaps(context, target_obj, rem_order, add_order):
        return _best_matrix_swap_loop(
            context,
            gram,
            base_score,
            rem_order,
            add_order,
            target_obj,
            max1_i,
            max1,
            max2,
        )

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
        scores = _score_many_from_grams(
            context,
            candidate_grams.reshape(-1, context.L, context.L),
            max_power.ravel(),
            target_obj,
        ).reshape(len(rem_chunk), add_count)

        pos = int(np.argmax(scores))
        score = float(scores.ravel()[pos])
        if score > best_score + 1e-12:
            rem_pos, add_pos = np.unravel_index(pos, scores.shape)
            best_score = score
            best_i = int(rem_chunk[rem_pos])
            best_j = int(add_order[add_pos])

    return best_i, best_j, best_score


def _best_matrix_swap_loop(
    context,
    gram,
    base_score,
    rem_order,
    add_order,
    target_obj,
    max1_i,
    max1,
    max2,
):
    best_score = base_score
    best_i = None
    best_j = None
    add_grams = context.row_grams[add_order]
    add_power = context.row_power[add_order]

    for antenna in rem_order:
        gram_minus = gram - context.row_grams[antenna]
        max_without_i = max2 if int(antenna) == max1_i else max1
        max_power = np.maximum(max_without_i, add_power)
        candidate_grams = gram_minus[None, :, :] + add_grams
        scores = _score_many_from_grams(
            context,
            candidate_grams,
            max_power,
            target_obj,
        )
        pos = int(np.argmax(scores))
        score = float(scores[pos])

        if score > best_score + 1e-12:
            best_score = score
            best_i = int(antenna)
            best_j = int(add_order[pos])

    return best_i, best_j, best_score


def _best_bf_swap(
    context,
    gram,
    base_score,
    rem_order,
    add_order,
    max1_i,
    max1,
    max2,
):
    q = _row_quadratic_scores(context.V, gram)
    base_trace = float(np.real(np.trace(gram @ gram.conj().T)))
    rem_power = context.row_power[rem_order]
    add_power = context.row_power[add_order]
    max_without = np.where(rem_order == max1_i, max2, max1)
    max_power = np.maximum(max_without[:, None], add_power[None, :])

    cross = np.abs(context.V[rem_order] @ context.V[add_order].conj().T) ** 2
    trace_values = (
        base_trace
        + 2.0 * (q[add_order][None, :] - q[rem_order][:, None])
        + add_power[None, :] ** 2
        + rem_power[:, None] ** 2
        - 2.0 * cross
    )
    scores = np.full(trace_values.shape, -np.inf, dtype=float)
    valid = max_power > 0
    scores[valid] = (context.P / max_power[valid]) * trace_values[valid]

    pos = int(np.argmax(scores))
    score = float(scores.ravel()[pos])
    if score <= base_score + 1e-12:
        return None, None, base_score

    rem_pos, add_pos = np.unravel_index(pos, scores.shape)
    return int(rem_order[rem_pos]), int(add_order[add_pos]), score


def _max_swap_batch_size(L):
    return max(1024, min(25000, 2_000_000 // max(1, L * L)))


def _use_batched_matrix_swaps(context, target_obj, rem_order, add_order):
    if target_obj == "int":
        return False
    return context.L <= 4 and len(rem_order) * len(add_order) <= _max_swap_batch_size(
        context.L
    )


def _candidate_orders(context, active, gram, target_obj, remove_limit, add_limit):
    active_idx = np.flatnonzero(active)
    inactive_idx = np.flatnonzero(~active)
    if len(active_idx) == 0 or len(inactive_idx) == 0:
        return np.array([], dtype=int), np.array([], dtype=int)

    if target_obj == "gen":
        grad = _logdet_gradient(context, gram, np.max(context.row_power[active]))
        grad_scores = _row_quadratic_scores(context.V, grad)
        eigvals, eigvecs = np.linalg.eigh(0.5 * (gram + gram.conj().T))
        weak_vec = eigvecs[:, int(np.argmin(eigvals))]
        strong_vec = eigvecs[:, int(np.argmax(eigvals))]
        weak_align = np.abs(context.V @ weak_vec) ** 2
        strong_align = np.abs(context.V @ strong_vec) ** 2
        rem_order = _interleave_limited(
            [
                active_idx[np.argsort(grad_scores[active_idx])],
                active_idx[np.argsort(context.row_power[active_idx])],
                active_idx[np.argsort(-context.row_power[active_idx])],
                active_idx[np.argsort(-(strong_align[active_idx] - weak_align[active_idx]))],
            ],
            remove_limit,
        )
        add_order = _interleave_limited(
            [
                inactive_idx[np.argsort(-grad_scores[inactive_idx])],
                inactive_idx[np.argsort(-weak_align[inactive_idx])],
                inactive_idx[np.argsort(-context.row_power[inactive_idx])],
                inactive_idx[np.argsort(context.row_power[inactive_idx])],
            ],
            add_limit,
        )
        return rem_order, add_order

    if target_obj == "bf":
        q = _row_quadratic_scores(context.V, gram)
        gain = 2.0 * q + context.row_power**2
        loss = 2.0 * q - context.row_power**2
        rem_order = _interleave_limited(
            [
                active_idx[np.argsort(loss[active_idx])],
                active_idx[np.argsort(context.row_power[active_idx])],
                active_idx[np.argsort(-context.row_power[active_idx])],
            ],
            remove_limit,
        )
        add_order = _interleave_limited(
            [
                inactive_idx[np.argsort(-gain[inactive_idx])],
                inactive_idx[np.argsort(-context.row_power[inactive_idx])],
                inactive_idx[np.argsort(context.row_power[inactive_idx])],
            ],
            add_limit,
        )
        return rem_order, add_order

    features = context.offdiag_features
    if features.shape[1] == 0:
        rem_order = active_idx[np.argsort(context.row_power[active_idx])]
        add_order = inactive_idx[np.argsort(-context.row_power[inactive_idx])]
        return rem_order[:remove_limit], add_order[:add_limit]

    current = _hermitian_features(gram[None, :, :], include_diag=False)[0]
    dot = features @ current
    norms = np.sum(features * features, axis=1)
    rem_score = dot - 0.5 * norms
    add_score = dot + 0.5 * norms
    rem_order = _interleave_limited(
        [
            active_idx[np.argsort(-rem_score[active_idx])],
            active_idx[np.argsort(-context.row_power[active_idx])],
            active_idx[np.argsort(context.row_power[active_idx])],
        ],
        remove_limit,
    )
    add_order = _interleave_limited(
        [
            inactive_idx[np.argsort(add_score[inactive_idx])],
            inactive_idx[np.argsort(context.row_power[inactive_idx])],
            inactive_idx[np.argsort(-context.row_power[inactive_idx])],
        ],
        add_limit,
    )
    return rem_order, add_order


def _logdet_gradient(context, gram, max_row_power):
    if max_row_power <= 0:
        return np.zeros_like(gram)
    scale = context.P / max_row_power
    matrix = scale * (gram @ gram.conj().T) + context.sigma * context.eye
    try:
        gram_times_inv = np.linalg.solve(matrix.T, gram.T).T
    except np.linalg.LinAlgError:
        gram_times_inv = gram @ np.linalg.pinv(matrix)
    grad = 2.0 * scale * gram_times_inv
    return 0.5 * (grad + grad.conj().T)


def _score_active(context, active, target_obj):
    gram = _build_gram(context, active)
    return _score_from_gram(
        context,
        gram,
        np.max(context.row_power[active]) if np.any(active) else 0.0,
        target_obj,
    )


def _score_from_gram(context, gram, max_row_power, target_obj):
    if max_row_power <= 0:
        if target_obj == "gen" and context.sigma > 0:
            return float(context.L * np.log(context.sigma))
        return -np.inf
    return float(
        _score_many_from_grams(
            context,
            gram[None, :, :],
            np.asarray([max_row_power], dtype=float),
            target_obj,
        )[0]
    )


def _score_many_from_grams(context, grams, max_row_power, target_obj):
    scores = np.full(len(grams), -np.inf, dtype=float)
    valid = max_row_power > 0
    if not np.any(valid):
        return scores

    grams_valid = grams[valid]
    z2 = context.P / max_row_power[valid]
    gram_sq = grams_valid @ np.swapaxes(grams_valid.conj(), 1, 2)

    if target_obj == "bf":
        scores[valid] = z2 * np.real(np.trace(gram_sq, axis1=1, axis2=2))
        return scores

    offdiag = np.sum(np.abs(gram_sq) ** 2, axis=(1, 2)) - np.sum(
        np.abs(np.diagonal(gram_sq, axis1=1, axis2=2)) ** 2,
        axis=1,
    )
    if target_obj == "int":
        scores[valid] = -(z2**2) * offdiag
        return scores

    matrices = z2[:, None, None] * gram_sq + context.sigma * context.eye[None, :, :]
    signs, logdets = np.linalg.slogdet(matrices)
    scores[valid] = np.where(signs > 0, np.real(logdets), -np.inf)
    return scores


def _build_gram(context, active):
    idx = np.flatnonzero(active)
    if len(idx) == 0:
        return np.zeros((context.L, context.L), dtype=complex)
    return context.V[idx].conj().T @ context.V[idx]


def _row_quadratic_scores(V, matrix):
    return np.real(np.einsum("ni,ij,nj->n", V, matrix, V.conj(), optimize=True))


def _hermitian_features(matrices, include_diag):
    matrices = np.asarray(matrices)
    _, L, _ = matrices.shape
    features = []
    if include_diag:
        for i in range(L):
            features.append(np.real(matrices[:, i, i]))
    scale = np.sqrt(2.0)
    for i in range(L):
        for j in range(i + 1, L):
            values = matrices[:, i, j]
            features.append(scale * np.real(values))
            features.append(scale * np.imag(values))
    if not features:
        return np.zeros((matrices.shape[0], 0), dtype=float)
    return np.column_stack(features).astype(float, copy=False)


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


def _indices_to_active(context, indices):
    active = np.zeros(context.N, dtype=bool)
    active[np.asarray(indices, dtype=int)] = True
    return _repair_active(context, active)


def _active_to_x(active):
    return active.astype(int)


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
