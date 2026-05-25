import numpy as np

from .common import default_min_active, objective_from_gram


def solve_h1(V, K, sigma=1.0, P=1.0, min_active=None):
    """
    H1 under the task constraint sum(x) <= K.

    Sort antennas by power and scan feasible prefixes of the strongest
    antennas. Return the feasible H1 set with the best general objective.
    """

    N, L = V.shape
    K, min_count = default_min_active(V, K, min_active)
    if K == 0:
        return np.zeros(N, dtype=int)

    row_power = np.sum(np.abs(V) ** 2, axis=1).real
    strongest_first = np.argsort(row_power)[::-1]
    gram = np.zeros((L, L), dtype=complex)
    best_x = None
    best_u_g = -np.inf

    for active_count, antenna in enumerate(strongest_first[:K], start=1):
        v = V[antenna, :]
        gram += v.conj()[:, None] * v[None, :]
        if active_count < min_count:
            continue

        _, _, u_g = objective_from_gram(
            gram, row_power[strongest_first[0]], L, sigma=sigma, P=P
        )
        if u_g > best_u_g:
            best_u_g = u_g
            x = np.zeros(N, dtype=int)
            x[strongest_first[:active_count]] = 1
            best_x = x

    return best_x if best_x is not None else np.zeros(N, dtype=int)
