import numpy as np


def normalize_channel_matrix(V):
    V = np.asarray(V, dtype=complex).copy()
    column_norms = np.linalg.norm(V, axis=0)
    V /= column_norms
    antenna_max = np.max(np.linalg.norm(V, axis=1))
    V /= antenna_max
    return V


def generate_v_from_rng(rng, N, L):
    V = rng.normal(size=(N, L)) + 1j * rng.normal(size=(N, L))
    return normalize_channel_matrix(V)


def generate_V(N, L):
    return generate_v_from_rng(np.random, N, L)
