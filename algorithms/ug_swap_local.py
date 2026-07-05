import time

import numpy as np

from .common import objective_from_gram


def refine_selection_by_ug_swaps(
    V,
    x0,
    max_swaps,
    sigma=1.0,
    P=1.0,
    *,
    K=None,
    candidate_pool=None,
    min_gain=1e-12,
):
    return refine_selection_by_ug_swaps_steps(
        V,
        x0,
        max_swaps_values=(int(max_swaps),),
        sigma=sigma,
        P=P,
        K=K,
        candidate_pool=candidate_pool,
        min_gain=min_gain,
    )[int(max_swaps)]


def refine_selection_by_ug_swaps_steps(
    V,
    x0,
    max_swaps_values=(0, 1, 2),
    sigma=1.0,
    P=1.0,
    *,
    K=None,
    candidate_pool=None,
    min_gain=1e-12,
):
    started_at = time.perf_counter()
    context = _SwapContext(V, x0, sigma=sigma, P=P, K=K)
    requested = sorted({int(value) for value in max_swaps_values})
    if not requested:
        return {}
    if requested[0] < 0:
        raise ValueError("max_swaps values must be non-negative.")

    max_requested = requested[-1]
    active = context.active.copy()
    current_gram = _build_gram(context, active)
    current_max_power = _active_max_power(context, active)
    current_values = objective_from_gram(
        current_gram,
        current_max_power,
        context.L,
        sigma=context.sigma,
        P=context.P,
    )
    initial_values = current_values
    add_pool = _prepare_add_pool(context, active, candidate_pool)
    evaluated_swap_count = 0
    swap_history = []

    snapshots = {
        0: _result(
            context,
            active,
            max_swaps=0,
            add_candidate_count=len(add_pool),
            evaluated_swap_count=0,
            swaps_applied=0,
            initial_values=initial_values,
            final_values=current_values,
            swap_history=swap_history,
            started_at=started_at,
        )
    }

    for pass_no in range(1, max_requested + 1):
        best = _best_ug_swap(
            context,
            active,
            current_gram,
            current_values[2],
            add_pool,
            min_gain=min_gain,
        )
        evaluated_swap_count += best["evaluated_swap_count"]
        if best["remove"] is None:
            for missing in range(pass_no, max_requested + 1):
                snapshots[missing] = _result(
                    context,
                    active,
                    max_swaps=missing,
                    add_candidate_count=len(add_pool),
                    evaluated_swap_count=evaluated_swap_count,
                    swaps_applied=len(swap_history),
                    initial_values=initial_values,
                    final_values=current_values,
                    swap_history=swap_history,
                    started_at=started_at,
                )
            break

        remove_idx = int(best["remove"])
        add_idx = int(best["add"])
        active[remove_idx] = False
        active[add_idx] = True
        current_gram = best["gram"]
        current_values = best["values"]
        swap_history.append(
            {
                "remove": remove_idx,
                "add": add_idx,
                "u_g": float(current_values[2]),
            }
        )
        snapshots[pass_no] = _result(
            context,
            active,
            max_swaps=pass_no,
            add_candidate_count=len(add_pool),
            evaluated_swap_count=evaluated_swap_count,
            swaps_applied=len(swap_history),
            initial_values=initial_values,
            final_values=current_values,
            swap_history=swap_history,
            started_at=started_at,
        )

    for value in requested:
        if value not in snapshots:
            snapshots[value] = snapshots[max(key for key in snapshots if key < value)]
    return {value: snapshots[value] for value in requested}


class _SwapContext:
    def __init__(self, V, x0, sigma=1.0, P=1.0, K=None):
        V = np.asarray(V)
        if V.ndim != 2:
            raise ValueError("V must be a 2D complex matrix of shape (N, L).")
        if not np.iscomplexobj(V):
            V = V.astype(np.complex128)

        x0 = np.asarray(x0, dtype=int)
        if x0.ndim != 1 or x0.shape[0] != V.shape[0]:
            raise ValueError("x0 must be a 1D selection vector with length N.")
        if not np.isin(x0, [0, 1]).all():
            raise ValueError("x0 must be binary.")

        self.V = V
        self.N, self.L = V.shape
        self.K = int(np.sum(x0) if K is None else K)
        if not (0 <= self.K <= self.N):
            raise ValueError("K must satisfy 0 <= K <= N.")
        self.sigma = float(sigma)
        self.P = float(P)
        self.eye = np.eye(self.L, dtype=complex)
        self.row_power = np.sum(np.abs(V) ** 2, axis=1).real
        self.row_grams = V.conj()[:, :, None] * V[:, None, :]
        self.active = _repair_active(self, x0.astype(bool))
        self.initial_x = self.active.astype(int)


def _prepare_add_pool(context, active, candidate_pool):
    if candidate_pool is None:
        return np.flatnonzero(~active).astype(int)

    pool = []
    seen = set()
    for value in candidate_pool:
        index = int(value)
        if not (0 <= index < context.N):
            continue
        if active[index] or index in seen:
            continue
        seen.add(index)
        pool.append(index)
    return np.asarray(pool, dtype=int)


def _best_ug_swap(context, active, gram, base_u_g, add_pool, min_gain):
    active_idx = np.flatnonzero(active)
    add_order = np.asarray([idx for idx in add_pool if not active[int(idx)]], dtype=int)
    if len(active_idx) == 0 or len(add_order) == 0:
        return _empty_swap()

    max1_i, max1, max2 = _top_active_powers(context, active)
    add_grams = context.row_grams[add_order]
    add_power = context.row_power[add_order]
    max_pairs_per_batch = _max_swap_batch_size(context.L)
    remove_chunk_size = max(1, max_pairs_per_batch // max(1, len(add_order)))

    best_score = float(base_u_g)
    best_remove = None
    best_add = None
    best_gram = None
    best_values = None
    evaluated = 0

    for start in range(0, len(active_idx), remove_chunk_size):
        rem_chunk = active_idx[start : start + remove_chunk_size]
        remove_grams = context.row_grams[rem_chunk]
        candidate_grams = (
            gram[None, None, :, :]
            - remove_grams[:, None, :, :]
            + add_grams[None, :, :, :]
        )
        max_without = np.where(rem_chunk == max1_i, max2, max1)
        max_power = np.maximum(max_without[:, None], add_power[None, :])
        flat_grams = candidate_grams.reshape(-1, context.L, context.L)
        flat_power = max_power.ravel()
        scores = _u_g_scores_from_grams(context, flat_grams, flat_power).reshape(
            len(rem_chunk),
            len(add_order),
        )
        evaluated += int(scores.size)

        pos = int(np.argmax(scores))
        score = float(scores.ravel()[pos])
        if score > best_score + float(min_gain):
            rem_pos, add_pos = np.unravel_index(pos, scores.shape)
            best_score = score
            best_remove = int(rem_chunk[rem_pos])
            best_add = int(add_order[add_pos])
            best_gram = flat_grams[pos].copy()
            best_values = objective_from_gram(
                best_gram,
                float(flat_power[pos]),
                context.L,
                sigma=context.sigma,
                P=context.P,
            )

    if best_remove is None:
        return _empty_swap(evaluated)
    return {
        "remove": best_remove,
        "add": best_add,
        "gram": best_gram,
        "values": best_values,
        "evaluated_swap_count": evaluated,
    }


def _u_g_scores_from_grams(context, grams, max_power):
    grams = np.asarray(grams)
    max_power = np.asarray(max_power, dtype=float)
    scores = np.zeros(len(grams), dtype=float)
    valid = max_power > 0.0
    if not np.any(valid):
        return scores

    if context.L == 2:
        scores[valid] = _u_g_scores_l2(context, grams[valid], max_power[valid])
        return scores

    z2 = context.P / max_power[valid]
    gram_sq = grams[valid] @ np.swapaxes(grams[valid].conj(), 1, 2)
    matrices = z2[:, None, None] * gram_sq + context.sigma * context.eye[None, :, :]
    dets = np.real(np.linalg.det(matrices))
    scores[valid] = np.where(np.isfinite(dets), dets, -np.inf)
    return scores


def _u_g_scores_l2(context, grams, max_power):
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
    return np.where(np.isfinite(det), np.real(det), -np.inf)


def _result(
    context,
    active,
    max_swaps,
    add_candidate_count,
    evaluated_swap_count,
    swaps_applied,
    initial_values,
    final_values,
    swap_history,
    started_at,
):
    x = active.astype(int)
    return {
        "x": x,
        "active_count": int(np.sum(x)),
        "max_swaps": int(max_swaps),
        "add_candidate_count": int(add_candidate_count),
        "evaluated_swap_count": int(evaluated_swap_count),
        "swaps_applied": int(swaps_applied),
        "initial_u_bf": float(initial_values[0]),
        "initial_u_i": float(initial_values[1]),
        "initial_u_g": float(initial_values[2]),
        "u_bf": float(final_values[0]),
        "u_i": float(final_values[1]),
        "u_g": float(final_values[2]),
        "u_g_db": float(10.0 * np.log10(max(float(final_values[2]), np.finfo(float).tiny))),
        "swap_history": _swap_history_to_string(swap_history),
        "elapsed_seconds": float(time.perf_counter() - started_at),
    }


def _swap_history_to_string(history):
    return ";".join(
        f"{int(item['remove'])}>{int(item['add'])}@{float(item['u_g']):.12g}"
        for item in history
    )


def _empty_swap(evaluated=0):
    return {
        "remove": None,
        "add": None,
        "gram": None,
        "values": None,
        "evaluated_swap_count": int(evaluated),
    }


def _repair_active(context, active):
    active = np.asarray(active, dtype=bool).copy()
    count = int(np.sum(active))
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


def _build_gram(context, active):
    if not np.any(active):
        return np.zeros((context.L, context.L), dtype=complex)
    return np.sum(context.row_grams[active], axis=0)


def _active_max_power(context, active):
    if not np.any(active):
        return 0.0
    return float(np.max(context.row_power[active]))


def _top_active_powers(context, active):
    active_idx = np.flatnonzero(active)
    if len(active_idx) == 0:
        return -1, 0.0, 0.0
    sorted_active = active_idx[np.argsort(context.row_power[active_idx])]
    max1_i = int(sorted_active[-1])
    max1 = float(context.row_power[max1_i])
    max2 = float(context.row_power[sorted_active[-2]]) if len(sorted_active) >= 2 else 0.0
    return max1_i, max1, max2


def _max_swap_batch_size(L):
    return max(1024, min(25000, 2_000_000 // max(1, int(L) * int(L))))
