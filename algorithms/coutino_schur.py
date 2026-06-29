import numpy as np


def solve_coutino_schur_greedy(
    V,
    K,
    sigma=1.0,
    P=1.0,
    random_state=None,
    jitter=1e-6,
    loading_fraction=0.5,
):
    """
    Coutino-Schur submodular greedy candidate.

    This implements the Schur-complement log-det surrogate from Coutino,
    Chepuri, and Leus for each stream column of V and sums the stream
    surrogates. The greedy step maximizes the exact marginal gain of
    this monotone submodular surrogate under the uniform cardinality
    constraint |A| = K.

    The guarantee is for the surrogate optimized here, not for the true
    antenna-selection objective with the max-row-power normalization.
    """

    del P, random_state

    V = np.asarray(V)
    if V.ndim != 2:
        raise ValueError("V must be a 2D complex matrix of shape (N, L).")
    if not np.iscomplexobj(V):
        V = V.astype(np.complex128)

    N, L = V.shape
    K = int(K)
    if not (0 <= K <= N):
        raise ValueError("K must satisfy 0 <= K <= N.")
    if K == 0:
        return np.zeros(N, dtype=int)
    if K == N:
        return np.ones(N, dtype=int)

    states = [
        _CoutinoStreamState(
            V,
            target_col=target_col,
            sigma=sigma,
            jitter=jitter,
            loading_fraction=loading_fraction,
        )
        for target_col in range(L)
    ]

    selected = np.zeros(N, dtype=bool)
    candidates = np.arange(N)

    for _ in range(K):
        remaining = candidates[~selected]
        gains = np.zeros(len(remaining), dtype=float)
        for state in states:
            gains += state.marginal_gains(remaining)

        best = int(remaining[int(np.argmax(gains))])
        selected[best] = True
        for state in states:
            state.add(best)

    return selected.astype(int)


class _CoutinoStreamState:
    def __init__(self, V, target_col, sigma, jitter, loading_fraction):
        self.N = V.shape[0]
        self.loading = _diagonal_loading(sigma, loading_fraction)
        self.add_weight = 1.0 / self.loading
        self.inv_matrix = _regularized_schur_inverse(
            V,
            target_col=target_col,
            sigma=sigma,
            loading=self.loading,
            jitter=jitter,
        )

    def marginal_gains(self, indices):
        diag_values = np.real(np.diag(self.inv_matrix)[: self.N][indices])
        diag_values = np.maximum(diag_values, 0.0)
        return np.log1p(self.add_weight * diag_values)

    def add(self, index):
        column = self.inv_matrix[:, index].copy()
        denom = 1.0 + self.add_weight * np.real(column[index])
        denom = max(float(denom), 1e-12)
        self.inv_matrix -= (self.add_weight / denom) * np.outer(column, column.conj())
        self.inv_matrix = 0.5 * (self.inv_matrix + self.inv_matrix.conj().T)


def _diagonal_loading(sigma, loading_fraction):
    noise_floor = max(float(sigma), 1e-6)
    fraction = float(np.clip(loading_fraction, 1e-3, 0.999))
    return max(fraction * noise_floor, 1e-9)


def _regularized_schur_inverse(V, target_col, sigma, loading, jitter):
    N, L = V.shape
    target = V[:, target_col]
    interferers = np.delete(V, target_col, axis=1)
    residual_noise = max(float(sigma) - loading, 1e-6)

    S_inv = _low_rank_loaded_inverse(interferers, residual_noise)
    steering = S_inv @ target
    lower = np.vdot(target, steering).real

    dim = N + 1
    matrix = np.zeros((dim, dim), dtype=complex)
    matrix[:N, :N] = S_inv
    matrix[:N, N] = steering
    matrix[N, :N] = steering.conj()
    matrix[N, N] = lower

    scale = max(1.0, float(np.linalg.norm(matrix, ord=np.inf)))
    matrix += max(float(jitter), 1e-12) * scale * np.eye(dim, dtype=complex)
    matrix = 0.5 * (matrix + matrix.conj().T)

    try:
        return np.linalg.inv(matrix)
    except np.linalg.LinAlgError:
        return np.linalg.pinv(matrix)


def _low_rank_loaded_inverse(U, diagonal):
    N = U.shape[0]
    if U.shape[1] == 0:
        return (1.0 / diagonal) * np.eye(N, dtype=complex)

    gram = U.conj().T @ U
    middle = np.eye(U.shape[1], dtype=complex) + gram / diagonal
    middle_inv = np.linalg.inv(middle)
    return (
        (1.0 / diagonal) * np.eye(N, dtype=complex)
        - (1.0 / diagonal**2) * U @ middle_inv @ U.conj().T
    )
