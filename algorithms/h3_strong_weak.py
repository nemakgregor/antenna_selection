import numpy as np

from .common import default_min_active


def solve_h3_strong_weak(V, K, sigma=1.0, P=1.0, min_active=None, random_state=None):
    """
    H3 split-off heuristic under the task constraint sum(x) <= K.

    Disable half of the antennas to turn off from the weakest rows and half
    from the strongest rows, leaving the middle-power rows active.
    """

    del random_state
    N, _ = V.shape
    K, _ = default_min_active(V, K, min_active)
    if K == 0:
        return np.zeros(N, dtype=int)
    if K == N:
        return np.ones(N, dtype=int)

    off_count = N - K
    weak_drop = off_count // 2
    strong_drop = off_count - weak_drop

    row_power = np.sum(np.abs(V) ** 2, axis=1).real
    power_order = np.argsort(row_power)

    active = np.ones(N, dtype=bool)
    if weak_drop:
        active[power_order[:weak_drop]] = False
    if strong_drop:
        active[power_order[N - strong_drop :]] = False

    return active.astype(int)
