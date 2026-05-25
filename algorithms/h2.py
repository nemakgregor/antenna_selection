import numpy as np

from .common import calculate_objectives, default_min_active


def solve_h2(V, K, sigma=1.0, P=1.0, min_active=None):
    """
    H2 under the task constraint sum(x) <= K.

    Follow iterative interference-based deletion, but continue below K and
    return the feasible set on that path with the best general objective.
    Uses the exact formula for L=2, and a generalized matrix approach for L>2.
    """

    N, L = V.shape
    K, min_count = default_min_active(V, K, min_active)
    if K == 0:
        return np.zeros(N, dtype=int)

    x = np.ones(N, dtype=int)
    active_count = N
    best_x = None
    best_u_g = -np.inf

    def record_if_feasible():
        nonlocal best_x, best_u_g
        if active_count <= K:
            _, _, u_g = calculate_objectives(V, x, sigma=sigma, P=P)
            if u_g > best_u_g:
                best_u_g = u_g
                best_x = x.copy()

    if L == 2:
        contributions = V[:, 0].conj() * V[:, 1]
        S = np.sum(contributions)

        record_if_feasible()
        while active_count > min_count:
            active_idx = np.where(x == 1)[0]
            residuals = np.abs(S - contributions[active_idx])
            best_n = active_idx[int(np.argmin(residuals))]

            x[best_n] = 0
            S -= contributions[best_n]
            active_count -= 1
            record_if_feasible()
    else:
        with np.errstate(all="ignore"):
            S = V.conj().T @ V
        np.fill_diagonal(S, 0)

        C = np.zeros((N, L, L), dtype=complex)
        for n in range(N):
            v_n = V[n, :]
            C_n = v_n.conj()[:, None] * v_n[None, :]
            np.fill_diagonal(C_n, 0)
            C[n] = C_n

        record_if_feasible()
        while active_count > min_count:
            active_idx = np.where(x == 1)[0]
            residuals = S - C[active_idx]
            interfs = np.sum(np.abs(residuals) ** 2, axis=(1, 2))
            best_n = active_idx[int(np.argmin(interfs))]

            x[best_n] = 0
            S -= C[best_n]
            active_count -= 1
            record_if_feasible()

    return best_x if best_x is not None else x
