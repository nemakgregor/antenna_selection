import numpy as np


DEFAULT_SHORTLIST_SIZE = 8
DEFAULT_MAX_ROUNDS = 250
DEFAULT_IMPROVEMENT_TOL = 1e-12


def r2_window(V, K, sigma=1.0, P=1.0, return_info=False):
    """
    Select the best contiguous power window using the R2 spike-background rule.

    The window score estimates the spectrum of the selected Gram matrix from
    two moments: T = tr(G) and Q = ||G||_F^2.  The selected set has exactly K
    active antennas unless K is zero.
    """

    V, K, sigma, P = _validate_inputs(V, K, sigma, P)
    N, L = V.shape
    if K == 0:
        return _format_result(np.zeros(N, dtype=int), return_info, {"score": 0.0})
    if K == N:
        x = np.ones(N, dtype=int)
        row_power = _row_power(V)
        gram = _gram_from_rows(V)
        info = {
            "score": _log_general_from_gram(
                gram,
                float(np.max(row_power)) if N else 0.0,
                sigma,
                P,
                np.eye(L, dtype=complex),
            ),
            "offset": 0,
            "cap": float(np.max(row_power)) if N else 0.0,
            "T": float(np.sum(row_power)),
            "Q": float(np.sum(np.abs(gram) ** 2)),
        }
        return _format_result(x, return_info, info)

    row_power = _row_power(V)
    order = np.argsort(row_power, kind="mergesort")[::-1]
    V_ordered = V[order]
    ordered_power = row_power[order]

    row_grams = V_ordered.conj()[:, :, None] * V_ordered[:, None, :]
    gram_prefix = np.concatenate(
        [np.zeros((1, L, L), dtype=complex), np.cumsum(row_grams, axis=0)],
        axis=0,
    )
    power_prefix = np.r_[0.0, np.cumsum(ordered_power)]

    best_score = -np.inf
    best_offset = 0
    best_info = {}
    for offset in range(N - K + 1):
        cap = float(ordered_power[offset])
        total_power = float(power_prefix[offset + K] - power_prefix[offset])
        gram = gram_prefix[offset + K] - gram_prefix[offset]
        second_moment = float(np.real(np.sum(np.abs(gram) ** 2)))
        score, a, b = _r2_score_from_moments(
            total_power,
            second_moment,
            cap,
            L,
            sigma,
            P,
        )
        if score > best_score:
            best_score = score
            best_offset = offset
            best_info = {
                "score": float(score),
                "offset": int(offset),
                "cap": cap,
                "T": total_power,
                "Q": second_moment,
                "a": float(a),
                "b": float(b),
            }

    x = np.zeros(N, dtype=int)
    x[order[best_offset : best_offset + K]] = 1
    return _format_result(x, return_info, best_info)


def solve_r2_delta_gen(
    V,
    K,
    sigma=1.0,
    P=1.0,
    random_state=None,
    shortlist_size=DEFAULT_SHORTLIST_SIZE,
    max_rounds=DEFAULT_MAX_ROUNDS,
    improvement_tol=DEFAULT_IMPROVEMENT_TOL,
    return_info=False,
):
    """
    Analytical R2+Delta solver for the general objective U_G.

    The algorithm starts from the R2 power window and then applies monotone
    pair swaps.  Swap candidates are shortlisted with the closed-form
    derivative of log det(sigma I + c G^2); every accepted swap is checked with
    the exact objective.
    """

    del random_state

    V, K, sigma, P = _validate_inputs(V, K, sigma, P)
    N, L = V.shape
    if K == 0:
        return _format_result(
            np.zeros(N, dtype=int),
            return_info,
            {"initial_score": 0.0, "final_score": 0.0, "rounds": 0},
        )
    if K == N:
        return _format_result(
            np.ones(N, dtype=int),
            return_info,
            {"initial_score": None, "final_score": None, "rounds": 0},
        )

    shortlist_size = max(1, int(shortlist_size))
    max_rounds = max(0, int(max_rounds))
    improvement_tol = float(improvement_tol)

    seed = r2_window(V, K, sigma=sigma, P=P, return_info=True)
    active = seed["x"].astype(bool)

    row_power = _row_power(V)
    row_grams = V.conj()[:, :, None] * V[:, None, :]
    eye = np.eye(L, dtype=complex)

    gram = np.sum(row_grams[active], axis=0)
    cap = float(np.max(row_power[active]))
    current_score = _log_general_from_gram(gram, cap, sigma, P, eye)
    initial_score = current_score
    if not np.isfinite(current_score):
        x = active.astype(int)
        return _format_result(
            x,
            return_info,
            {
                "initial_score": current_score,
                "final_score": current_score,
                "rounds": 0,
                "window": seed,
            },
        )

    rounds = 0
    for _ in range(max_rounds):
        if cap <= 0.0:
            break

        c = float(P) / cap
        eigvals, eigvecs = np.linalg.eigh(_hermitian_part(gram))
        eigvals = np.maximum(np.real(eigvals), 0.0)
        weights = eigvals / (float(sigma) + c * eigvals**2)
        projections = np.abs(V @ eigvecs) ** 2
        usefulness = np.real(projections @ weights)

        active_idx = np.flatnonzero(active)
        inactive_idx = np.flatnonzero(~active)
        if len(active_idx) == 0 or len(inactive_idx) == 0:
            break

        remove_count = min(shortlist_size, len(active_idx))
        remove_local = np.argsort(usefulness[active_idx], kind="mergesort")[
            :remove_count
        ]
        remove_candidates = active_idx[remove_local]

        cap_tol = max(1e-12, abs(cap) * 1e-12)
        eligible_add = inactive_idx[row_power[inactive_idx] <= cap + cap_tol]
        if len(eligible_add) == 0:
            break
        add_count = min(shortlist_size, len(eligible_add))
        add_local = np.argsort(usefulness[eligible_add], kind="mergesort")[
            -add_count:
        ][::-1]
        add_candidates = eligible_add[add_local]

        best_pair = None
        best_gram = None
        best_cap = None
        best_score = current_score
        for remove_idx in remove_candidates:
            base_active = active.copy()
            base_active[remove_idx] = False
            base_gram = gram - row_grams[remove_idx]
            for add_idx in add_candidates:
                if base_active[add_idx]:
                    continue
                trial_active = base_active.copy()
                trial_active[add_idx] = True
                trial_cap = float(np.max(row_power[trial_active]))
                trial_gram = base_gram + row_grams[add_idx]
                trial_score = _log_general_from_gram(
                    trial_gram,
                    trial_cap,
                    sigma,
                    P,
                    eye,
                )
                if trial_score > best_score + improvement_tol:
                    best_score = trial_score
                    best_pair = (int(remove_idx), int(add_idx))
                    best_gram = trial_gram
                    best_cap = trial_cap

        if best_pair is None:
            break

        remove_idx, add_idx = best_pair
        active[remove_idx] = False
        active[add_idx] = True
        gram = _hermitian_part(best_gram)
        cap = float(best_cap)
        current_score = float(best_score)
        rounds += 1

    x = active.astype(int)
    info = {
        "initial_score": float(initial_score),
        "final_score": float(current_score),
        "rounds": int(rounds),
        "window": seed,
    }
    return _format_result(x, return_info, info)


def _validate_inputs(V, K, sigma, P):
    V = np.asarray(V)
    if V.ndim != 2:
        raise ValueError("V must be a 2D complex matrix of shape (N, L).")
    if not np.iscomplexobj(V):
        V = V.astype(np.complex128)
    else:
        V = V.astype(np.complex128, copy=False)

    N = V.shape[0]
    K = int(K)
    if not (0 <= K <= N):
        raise ValueError("K must satisfy 0 <= K <= N.")
    sigma = float(sigma)
    P = float(P)
    if sigma <= 0.0:
        raise ValueError("sigma must be positive.")
    if P <= 0.0:
        raise ValueError("P must be positive.")
    return V, K, sigma, P


def _format_result(x, return_info, info):
    if not return_info:
        return x
    result = dict(info)
    result["x"] = x
    return result


def _row_power(V):
    return np.real(np.sum(np.abs(V) ** 2, axis=1))


def _gram_from_rows(V):
    if V.shape[0] == 0:
        return np.zeros((V.shape[1], V.shape[1]), dtype=complex)
    return V.conj().T @ V


def _r2_score_from_moments(total_power, second_moment, cap, L, sigma, P):
    if cap <= 0.0:
        return float(L * np.log(float(sigma))), 0.0, 0.0

    c = float(P) / cap
    if L == 1:
        value = max(float(total_power), 0.0)
        return float(np.log(float(sigma) + c * value**2)), value, 0.0

    total_power = max(float(total_power), 0.0)
    second_moment = max(float(second_moment), 0.0)
    discriminant = max((L - 1) * (L * second_moment - total_power**2), 0.0)
    spike = (total_power + np.sqrt(discriminant)) / L
    background = (total_power - spike) / (L - 1)
    background = max(float(background), 0.0)
    score = np.log(float(sigma) + c * spike**2) + (L - 1) * np.log(
        float(sigma) + c * background**2
    )
    return float(score), float(spike), float(background)


def _log_general_from_gram(gram, max_row_power, sigma, P, eye):
    if max_row_power <= 0.0:
        return float(eye.shape[0] * np.log(float(sigma)))

    gram = _hermitian_part(gram)
    gram_sq = gram @ gram.conj().T
    matrix = float(sigma) * eye + (float(P) / float(max_row_power)) * gram_sq
    matrix = _hermitian_part(matrix)
    sign, logdet = np.linalg.slogdet(matrix)
    if sign <= 0.0 or not np.isfinite(logdet):
        return -np.inf
    return float(np.real(logdet))


def _hermitian_part(matrix):
    return 0.5 * (matrix + matrix.conj().T)
