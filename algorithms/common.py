import numpy as np


def calculate_objectives(V, x, sigma=1.0, P=1.0):
    """
    Calculate the three objectives based on the formal problem statement.
    """

    _, L = V.shape
    active_idx = np.where(x == 1)[0]

    if len(active_idx) == 0:
        return 0.0, 0.0, 0.0

    row_norms_sq = np.sum(np.abs(V[active_idx, :]) ** 2, axis=1)
    z = np.sqrt(P / np.max(row_norms_sq)) if np.max(row_norms_sq) > 0 else 0

    W = np.zeros_like(V, dtype=complex)
    W[active_idx, :] = z * V[active_idx, :]

    with np.errstate(all="ignore"):
        V_eq = V.conj().T @ W
        V_eq_sq = V_eq @ V_eq.conj().T

    u_bf = np.real(np.trace(V_eq_sq))
    u_i = np.sum(np.abs(V_eq_sq) ** 2) - np.sum(
        np.abs(np.diag(V_eq_sq)) ** 2
    )
    u_g = np.real(np.linalg.det(V_eq_sq + sigma * np.eye(L)))

    return u_bf, u_i, u_g


def check_constraints(x, K):
    """
    Check if the cardinality constraint is met (sum(x_n) <= K).
    """

    num_active = np.sum(x)
    return num_active <= K, num_active


def objective_from_gram(gram, max_row_power, L, sigma=1.0, P=1.0):
    if max_row_power <= 0:
        return 0.0, 0.0, 0.0

    z2 = P / max_row_power
    gram_sq = gram @ gram.conj().T
    u_bf = z2 * np.real(np.trace(gram_sq))
    u_i = (z2**2) * (
        np.sum(np.abs(gram_sq) ** 2)
        - np.sum(np.abs(np.diag(gram_sq)) ** 2)
    )
    u_g = np.real(np.linalg.det(z2 * gram_sq + sigma * np.eye(L)))
    return float(u_bf), float(u_i), float(u_g)


def default_min_active(V, K, min_active=None):
    N, L = V.shape
    K = int(np.clip(K, 0, N))
    if K == 0:
        return K, 0

    min_count = L if min_active is None else int(min_active)
    min_count = int(np.clip(min_count, 1, K))
    return K, min_count
