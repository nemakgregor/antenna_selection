import numpy as np


DATA_PROFILES = (
    "gaussian",
    "rayleigh",
    "rician",
    "nakagami",
    "lognormal",
    "thin_tail",
    "twdp",
)

TWDP_K = 8.0
TWDP_DELTA = 0.9
TWDP_OMEGA = 1.0


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


def _uniform_phase(rng, shape):
    return np.exp(1j * rng.uniform(0.0, 2.0 * np.pi, size=shape))


def _twdp_specular_amplitudes(twdp_K=TWDP_K, twdp_Delta=TWDP_DELTA, omega=TWDP_OMEGA):
    specular_power = omega * twdp_K / (twdp_K + 1.0)
    balance = np.sqrt(max(0.0, 1.0 - twdp_Delta**2))
    v1 = np.sqrt(0.5 * specular_power * (1.0 + balance))
    v2 = np.sqrt(0.5 * specular_power * (1.0 - balance))
    diffuse_sigma = np.sqrt(omega / (2.0 * (twdp_K + 1.0)))
    return v1, v2, diffuse_sigma


def generate_v_profile_from_rng(rng, N, L, profile="gaussian"):
    """
    Generate normalized channel matrices for threshold-distribution experiments.

    The default `gaussian` profile is intentionally identical to the existing
    generator. Other profiles are synthetic fading/tail stress cases, not a
    replacement for measured channel data.
    """

    if profile == "gaussian":
        return generate_v_from_rng(rng, N, L)

    if profile == "rayleigh":
        amplitudes = rng.rayleigh(scale=1.0, size=(N, L))
        V = amplitudes * _uniform_phase(rng, (N, L))
    elif profile == "rician":
        k_factor = 4.0
        scatter = rng.normal(size=(N, L)) + 1j * rng.normal(size=(N, L))
        los_phase = _uniform_phase(rng, (1, L))
        los = np.ones((N, L), dtype=complex) * los_phase
        V = (
            np.sqrt(k_factor / (k_factor + 1.0)) * los
            + np.sqrt(1.0 / (k_factor + 1.0)) * scatter
        )
    elif profile == "nakagami":
        m_shape = 2.0
        omega = 1.0
        amplitudes_sq = rng.gamma(
            shape=m_shape,
            scale=omega / m_shape,
            size=(N, L),
        )
        V = np.sqrt(amplitudes_sq) * _uniform_phase(rng, (N, L))
    elif profile == "lognormal":
        shadow = rng.lognormal(mean=0.0, sigma=0.6, size=(N, 1))
        scatter = rng.normal(size=(N, L)) + 1j * rng.normal(size=(N, L))
        V = shadow * scatter
    elif profile == "thin_tail":
        row_scale = 0.2 + 0.8 * rng.beta(a=2.0, b=8.0, size=(N, 1))
        scatter = rng.normal(size=(N, L)) + 1j * rng.normal(size=(N, L))
        V = row_scale * scatter
    elif profile == "twdp":
        v1, v2, diffuse_sigma = _twdp_specular_amplitudes()
        diffuse = diffuse_sigma * (
            rng.normal(size=(N, L)) + 1j * rng.normal(size=(N, L))
        )
        V = (
            v1 * _uniform_phase(rng, (N, L))
            + v2 * _uniform_phase(rng, (N, L))
            + diffuse
        )
    else:
        raise ValueError(
            f"Unknown data profile {profile!r}; expected one of {DATA_PROFILES}."
        )

    return normalize_channel_matrix(V)


def generate_V_profile(N, L, profile="gaussian"):
    return generate_v_profile_from_rng(np.random, N, L, profile=profile)


def generate_V(N, L):
    return generate_v_from_rng(np.random, N, L)
