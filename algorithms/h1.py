import numpy as np

from .common import default_min_active


def solve_h1(V, K, sigma=1.0, P=1.0, min_active=None):
    """
    H1 from the task statement.

    Sort antennas by row power and switch off the weakest rows until exactly
    K antennas remain active.
    """

    N, _ = V.shape
    K, _ = default_min_active(V, K, min_active)
    if K == 0:
        return np.zeros(N, dtype=int)
    if K == N:
        return np.ones(N, dtype=int)

    row_power = np.sum(np.abs(V) ** 2, axis=1).real
    weakest_first = np.argsort(row_power)
    off_count = N - K

    x = np.ones(N, dtype=int)
    x[weakest_first[:off_count]] = 0
    return x
