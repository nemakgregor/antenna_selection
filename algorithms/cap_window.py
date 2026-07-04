import numpy as np


DEFAULT_SCAN_SIZE = 129


def solve_cap_window_gen(
    V,
    K,
    sigma=1.0,
    P=1.0,
    random_state=None,
    scan_size=DEFAULT_SCAN_SIZE,
):
    """
    Fast cap-window scan for the general objective.

    The method sorts antennas by row power, evaluates a compact set of
    contiguous power windows with the exact log-det objective, and returns the
    best window. The H3 strong/weak split is always one of the windows.
    """

    del random_state

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

    row_power = np.sum(np.abs(V) ** 2, axis=1).real
    order = np.argsort(row_power)[::-1]
    offsets = _candidate_offsets(N, K, scan_size)
    ordered_power = row_power[order]

    if L == 1:
        best_offset = _best_scalar_window_offset(
            V[order, 0],
            ordered_power,
            offsets,
            K,
            sigma,
            P,
        )
    else:
        best_offset = _best_matrix_window_offset(
            V[order],
            ordered_power,
            offsets,
            K,
            sigma,
            P,
        )

    x = np.zeros(N, dtype=int)
    x[order[best_offset : best_offset + K]] = 1
    return x


def solve_cap_window_full_gen(V, K, sigma=1.0, P=1.0, random_state=None):
    """
    Exhaustive cap-window scan for the general objective.

    This is the same cap-window family as solve_cap_window_gen, but it evaluates
    every contiguous power window of length K. It is still fast for N=1000,
    L=2, and gives the exact best candidate inside this H3-style window family.
    """

    V = np.asarray(V)
    if V.ndim != 2:
        raise ValueError("V must be a 2D complex matrix of shape (N, L).")
    N = V.shape[0]
    K = int(K)
    scan_size = max(1, N - K + 1)
    return solve_cap_window_gen(
        V,
        K,
        sigma=sigma,
        P=P,
        random_state=random_state,
        scan_size=scan_size,
    )


def _candidate_offsets(N, K, scan_size):
    max_offset = N - K
    off_count = max_offset
    weak_drop = off_count // 2
    h3_offset = off_count - weak_drop

    raw = [
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        8,
        10,
        12,
        15,
        20,
        25,
        30,
        40,
        50,
        75,
        100,
        N // 200,
        N // 100,
        N // 50,
        N // 25,
        N // 20,
        N // 10,
        h3_offset,
        max_offset,
    ]
    for radius in (3, 5, 10, 20, 30, 50, 75, 100, 150, 200):
        raw.extend((h3_offset - radius, h3_offset + radius))

    grid_count = min(max(1, int(scan_size)), max_offset + 1)
    raw.extend(np.linspace(0, max_offset, grid_count, dtype=int).tolist())

    offsets = []
    seen = set()
    for offset in raw:
        offset = int(np.clip(offset, 0, max_offset))
        if offset in seen:
            continue
        seen.add(offset)
        offsets.append(offset)
    return np.asarray(sorted(offsets), dtype=int)


def _best_scalar_window_offset(values, ordered_power, offsets, K, sigma, P):
    values_sq = np.abs(values) ** 2
    prefix = np.r_[0.0, np.cumsum(values_sq)]
    grams = prefix[offsets + K] - prefix[offsets]
    caps = ordered_power[offsets]
    scores = np.full(len(offsets), -np.inf, dtype=float)
    valid = caps > 0.0
    scores[valid] = np.log(float(sigma) + (float(P) / caps[valid]) * grams[valid] ** 2)
    return int(offsets[int(np.argmax(scores))])


def _best_matrix_window_offset(V_ordered, ordered_power, offsets, K, sigma, P):
    _, L = V_ordered.shape
    row_grams = V_ordered.conj()[:, :, None] * V_ordered[:, None, :]
    prefix = np.concatenate(
        [np.zeros((1, L, L), dtype=complex), np.cumsum(row_grams, axis=0)],
        axis=0,
    )
    grams = prefix[offsets + K] - prefix[offsets]
    caps = ordered_power[offsets]

    scores = np.full(len(offsets), -np.inf, dtype=float)
    valid = caps > 0.0
    if np.any(valid):
        if L == 2:
            scores[valid] = _logdet_scores_l2(grams[valid], caps[valid], sigma, P)
        else:
            gram_sq = grams[valid] @ np.swapaxes(grams[valid].conj(), 1, 2)
            matrices = (
                (float(P) / caps[valid])[:, None, None] * gram_sq
                + float(sigma) * np.eye(L, dtype=complex)[None, :, :]
            )
            signs, logdets = np.linalg.slogdet(matrices)
            scores[valid] = np.where(signs > 0, np.real(logdets), -np.inf)

    return int(offsets[int(np.argmax(scores))])


def _logdet_scores_l2(grams, caps, sigma, P):
    a = np.real(grams[:, 0, 0])
    d = np.real(grams[:, 1, 1])
    b_abs_sq = np.abs(grams[:, 0, 1]) ** 2

    sq00 = a * a + b_abs_sq
    sq11 = d * d + b_abs_sq
    off_abs_sq = b_abs_sq * (a + d) ** 2
    scale = float(P) / caps

    dets = (
        (float(sigma) + scale * sq00)
        * (float(sigma) + scale * sq11)
        - (scale**2) * off_abs_sq
    )
    return np.where(dets > 0.0, np.log(dets), -np.inf)
