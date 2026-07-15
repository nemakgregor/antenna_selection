import numpy as np

from .common import calculate_objectives, default_min_active, objective_from_gram
from .h1 import solve_h1
from .h2 import solve_h2


def solve_pareto_interference_greedy(
    V,
    K,
    sigma=1.0,
    P=1.0,
    step_bf_keep=0.95,
    h2_bf_floor=1.0,
    h12_bf_floor=0.85,
    h2_interference_ceiling=2.0,
    min_active=None,
    random_state=None,
):
    """
    Interference-first greedy with BF protection.

    Follows a deletion path that minimizes interference only among candidates
    whose BF gain remains close to the best BF candidate at the current step.
    """

    del random_state
    N, L = V.shape
    K, min_count = default_min_active(V, K, min_active)
    if K == 0:
        return np.zeros(N, dtype=int)

    with np.errstate(all="ignore"):
        x_h1 = solve_h1(V, K, sigma=sigma, P=P, min_active=min_active)
        h1_bf, h1_i, h1_g = calculate_objectives(V, x_h1, sigma=sigma, P=P)
        x_h2 = solve_h2(V, K, sigma=sigma, P=P, min_active=min_active)
        h2_bf, h2_i, h2_g = calculate_objectives(V, x_h2, sigma=sigma, P=P)

    row_power = np.sum(np.abs(V) ** 2, axis=1).real
    row_grams = V.conj()[:, :, None] * V[:, None, :]
    gram = np.sum(row_grams, axis=0)
    active = np.ones(N, dtype=bool)
    off_diag = ~np.eye(L, dtype=bool)
    records = []

    def record_state():
        active_count = int(np.sum(active))
        if active_count <= K:
            active_power = row_power[active]
            max_power = np.max(active_power) if active_count > 0 else 0.0
            u_bf, u_i, u_g = objective_from_gram(
                gram, max_power, L, sigma=sigma, P=P
            )
            records.append((u_bf, u_i, u_g, active.astype(int).copy()))

    for active_count in range(N, min_count - 1, -1):
        record_state()
        if active_count == min_count:
            break

        active_idx = np.flatnonzero(active)
        candidate_grams = gram[None, :, :] - row_grams[active_idx]

        candidate_power = row_power[active_idx]
        current_max_power = np.max(candidate_power)
        z2 = np.zeros(len(active_idx), dtype=float)
        if current_max_power > 0:
            z2.fill(P / current_max_power)

            max_positions = np.flatnonzero(candidate_power == current_max_power)
            if len(max_positions) == 1 and len(active_idx) > 1:
                max_pos = max_positions[0]
                second_max_power = np.max(np.delete(candidate_power, max_pos))
                if second_max_power > 0:
                    z2[max_pos] = P / second_max_power

        candidate_sq = candidate_grams @ candidate_grams.conj().transpose(0, 2, 1)
        candidate_bf = z2 * np.trace(candidate_sq, axis1=1, axis2=2).real
        candidate_i = (z2**2) * np.sum(
            np.abs(candidate_sq[:, off_diag]) ** 2, axis=1
        )

        bf_allowed = candidate_bf >= step_bf_keep * np.max(candidate_bf)
        scores = np.where(bf_allowed, candidate_i, np.inf)
        if not np.isfinite(scores).any():
            scores = candidate_i

        antenna_to_delete = active_idx[int(np.argmin(scores))]
        active[antenna_to_delete] = False
        gram -= row_grams[antenna_to_delete]

    records.append((h1_bf, h1_i, h1_g, x_h1.copy()))
    records.append((h2_bf, h2_i, h2_g, x_h2.copy()))
    bf_floor = max(h2_bf_floor * h2_bf, h12_bf_floor * max(h1_bf, h2_bf))
    interference_ceiling = h2_interference_ceiling * max(h2_i, np.finfo(float).eps)
    feasible = [
        record
        for record in records
        if record[0] >= bf_floor and record[1] <= interference_ceiling
    ]
    if feasible:
        return min(feasible, key=lambda record: record[1])[3]

    if h2_bf >= h2_bf_floor * h2_bf:
        return x_h2

    eps = np.finfo(float).eps
    return min(
        records,
        key=lambda record: np.log1p(record[1] / max(h2_i, eps))
        + max(0.0, bf_floor - record[0]) / max(h2_bf, eps),
    )[3]
