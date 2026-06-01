import numpy as np

from .common import default_min_active


def solve_h2(V, K, sigma=1.0, P=1.0, min_active=None):
    """
    H2 from the task statement.

    Follow iterative interference-based deletion until exactly K antennas
    remain active.
    Uses the exact formula for L=2, and a generalized matrix approach for L>2.
    """

    N, L = V.shape
    K, _ = default_min_active(V, K, min_active)
    if K == 0:
        return np.zeros(N, dtype=int)
    if K == N:
        return np.ones(N, dtype=int)

    x = np.ones(N, dtype=int)
    active_count = N

    if L == 2:
        contributions = V[:, 0].conj() * V[:, 1]
        S = np.sum(contributions)

        while active_count > K:
            active_idx = np.where(x == 1)[0]
            residuals = np.abs(S - contributions[active_idx])
            best_n = active_idx[int(np.argmin(residuals))]

            x[best_n] = 0
            S -= contributions[best_n]
            active_count -= 1
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

        while active_count > K:
            active_idx = np.where(x == 1)[0]
            residuals = S - C[active_idx]
            interfs = np.sum(np.abs(residuals) ** 2, axis=(1, 2))
            best_n = active_idx[int(np.argmin(interfs))]

            x[best_n] = 0
            S -= C[best_n]
            active_count -= 1

    return x
