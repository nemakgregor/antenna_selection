import numpy as np

from .common import calculate_objectives
from .coutino import solve_coutino_greedy
from .h1 import solve_h1
from .h2 import solve_h2


def solve_miso_energy_greedy(
    V, K, sigma=1.0, P=1.0, target_margin=0.05, min_active=None
):
    """
    Energy-aware MISO antenna selection under the task constraint sum(x) <= K.

    Deletes antennas by preserving the task's general objective, then returns
    the smallest active set on that greedy path that still beats the best of
    H1/H2 by target_margin on U_G.
    """

    N, L = V.shape
    K = int(np.clip(K, 0, N))
    min_active = L if min_active is None else int(min_active)
    min_active = int(np.clip(min_active, 1, max(1, N)))

    if K <= 0:
        return np.zeros(N, dtype=int)
    if min_active > K:
        min_active = K

    with np.errstate(all="ignore"):
        _, _, h1_g = calculate_objectives(V, solve_h1(V, K), sigma=sigma, P=P)
        _, _, h2_g = calculate_objectives(V, solve_h2(V, K), sigma=sigma, P=P)
    baseline_g = max(h1_g, h2_g)
    target_log_g = np.log(max(baseline_g * (1.0 + target_margin), np.finfo(float).tiny))

    row_power = np.sum(np.abs(V) ** 2, axis=1).real
    row_grams = V.conj()[:, :, None] * V[:, None, :]
    gram = np.sum(row_grams, axis=0)

    active = np.ones(N, dtype=bool)
    eye = np.eye(L, dtype=complex)
    best_log_g = -np.inf
    best_x = None
    smallest_target_x = None

    for active_count in range(N, min_active - 1, -1):
        active_power = row_power[active]
        max_power = np.max(active_power) if active_count > 0 else 0.0
        z2 = P / max_power if max_power > 0 else 0.0
        gram_sq = gram @ gram.conj().T
        sign, log_g = np.linalg.slogdet(sigma * eye + z2 * gram_sq)
        log_g = log_g.real if sign > 0.0 and np.isfinite(log_g) else -np.inf

        if active_count <= K:
            if log_g > best_log_g:
                best_log_g = log_g
                best_x = active.astype(int).copy()
            if log_g >= target_log_g:
                smallest_target_x = active.astype(int).copy()

        if active_count == min_active:
            break

        active_idx = np.flatnonzero(active)
        candidate_grams = gram[None, :, :] - row_grams[active_idx]

        candidate_power = row_power[active_idx]
        current_max_power = np.max(candidate_power)
        z2_candidates = np.zeros(len(active_idx), dtype=float)
        if current_max_power > 0:
            z2_candidates.fill(P / current_max_power)

            max_positions = np.flatnonzero(candidate_power == current_max_power)
            if len(max_positions) == 1 and len(active_idx) > 1:
                max_pos = max_positions[0]
                second_max_power = np.max(np.delete(candidate_power, max_pos))
                if second_max_power > 0:
                    z2_candidates[max_pos] = P / second_max_power

        candidate_sq = candidate_grams @ candidate_grams.conj().transpose(0, 2, 1)
        objective_mats = sigma * eye + z2_candidates[:, None, None] * candidate_sq
        signs, log_abs_dets = np.linalg.slogdet(objective_mats)
        log_abs_dets = np.where(
            np.isfinite(log_abs_dets) & (signs > 0.0),
            log_abs_dets.real,
            -np.inf,
        )

        antenna_to_delete = active_idx[int(np.argmax(log_abs_dets))]
        active[antenna_to_delete] = False
        gram -= row_grams[antenna_to_delete]

    if smallest_target_x is not None:
        return smallest_target_x
    if best_x is not None:
        return best_x
    return solve_coutino_greedy(V, K, sigma=sigma, P=P)
