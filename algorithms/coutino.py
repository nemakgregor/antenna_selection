import numpy as np

from .common import default_min_active


def solve_coutino_greedy(V, K, sigma=1.0, P=1.0, min_active=None):
    """
    Coutino-style greedy deletion adapted to the task constraint sum(x) <= K.

    Greedily delete the antenna whose removal leaves the largest
    log det(V_eq V_eq^* + sigma I), continue below K, and return the feasible
    set on that path with the best general objective.
    """

    N, L = V.shape
    K, min_count = default_min_active(V, K, min_active)
    if K == 0:
        return np.zeros(N, dtype=int)

    row_power = np.sum(np.abs(V) ** 2, axis=1).real
    row_grams = V.conj()[:, :, None] * V[:, None, :]
    gram = np.sum(row_grams, axis=0)

    active = np.ones(N, dtype=bool)
    eye = np.eye(L, dtype=complex)
    best_x = None
    best_log_g = -np.inf

    for active_count in range(N, min_count - 1, -1):
        if active_count <= K:
            active_power = row_power[active]
            max_power = np.max(active_power) if active_count > 0 else 0.0
            z2 = P / max_power if max_power > 0 else 0.0
            gram_sq = gram @ gram.conj().T
            sign, log_g = np.linalg.slogdet(sigma * eye + z2 * gram_sq)
            log_g = log_g.real if sign > 0.0 and np.isfinite(log_g) else -np.inf
            if log_g > best_log_g:
                best_log_g = log_g
                best_x = active.astype(int).copy()

        if active_count == min_count:
            break

        active_idx = np.flatnonzero(active)
        candidate_grams = gram[None, :, :] - row_grams[active_idx]

        active_power = row_power[active_idx]
        current_max_power = np.max(active_power)
        z2 = np.zeros(len(active_idx), dtype=float)
        if current_max_power > 0:
            z2.fill(P / current_max_power)

            max_positions = np.flatnonzero(active_power == current_max_power)
            if len(max_positions) == 1 and len(active_idx) > 1:
                max_pos = max_positions[0]
                second_max_power = np.max(np.delete(active_power, max_pos))
                if second_max_power > 0:
                    z2[max_pos] = P / second_max_power

        candidate_sq = candidate_grams @ candidate_grams.conj().transpose(0, 2, 1)
        objective_mats = sigma * eye + z2[:, None, None] * candidate_sq
        signs, log_abs_dets = np.linalg.slogdet(objective_mats)
        log_abs_dets = np.where(
            np.isfinite(log_abs_dets) & (signs > 0.0),
            log_abs_dets.real,
            -np.inf,
        )

        antenna_to_delete = active_idx[int(np.argmax(log_abs_dets))]
        active[antenna_to_delete] = False
        gram -= row_grams[antenna_to_delete]

    return best_x if best_x is not None else active.astype(int)
