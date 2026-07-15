from functools import partial

import numpy as np

from algorithms import (
    solve_cap_submodular_gen,
    solve_cap_submodular_portfolio_gen,
    solve_cap_window_full_gen,
    solve_cap_window_gen,
    solve_coutino_schur_greedy,
    solve_frame_portfolio,
    solve_h1,
    solve_h2,
    solve_h3,
    solve_h3_fast,
    solve_h3_strong_weak,
    solve_miso_energy_greedy,
    solve_pareto_interference_greedy,
    refine_general_1swap,
    solve_r2_delta_gen,
    solve_thresholded_logdet_greedy,
    solve_true_backward_greedy,
)


FAST_FRAME_KWARGS = {
    "max_refined_starts": 3,
    "max_passes": 2,
    "remove_limit": 60,
    "add_limit": 60,
    "lambdas": (),
}

EXACT_FRAME_KWARGS = {
    "max_refined_starts": 12,
    "max_passes": 20,
    "remove_limit": None,
    "add_limit": None,
    "random_restarts": 20,
}

FAST_THRESH_LOGDET_KWARGS = {
    "eps_values": (1e-6, 1e-2),
    "lambdas": (1.0,),
    "threshold_scan_size": 5,
}

WEIGHTED_THRESH_LOGDET_KWARGS = {
    "eps_values": (1e-6, 1e-2),
    "lambdas": (0.9, 0.7, 0.5),
    "threshold_scan_size": 5,
}

SWAP_THRESH_LOGDET_KWARGS = {
    "eps_values": (1e-6, 1e-2),
    "lambdas": (1.0,),
    "threshold_scan_size": 5,
    "swap_max_passes": 1,
}

FAST_1SWAP_LS_KWARGS = {
    "max_passes": 3,
    "remove_limit": 64,
    "add_limit": 128,
    "boundary_limit": 128,
}


def _round_half_up(value):
    return int(np.floor(float(value) + 0.5))


def _h3_t1(N, K, L):
    del K, L
    return _round_half_up(0.05 * N)


def _h3_t2(N, K, L):
    del N, L
    return _round_half_up(0.15 * K)


def _h3_t3(N, K, L):
    return _round_half_up(0.125 * N * L / (L + 2))


def _h3_threshold_with_ts(V, K, sigma, P, t_values):
    return solve_h3(
        V,
        K,
        target_obj="gen",
        sigma=sigma,
        P=P,
        t_tests=tuple(t_values),
        include_phase_nulling=False,
    )


def _h3_threshold_t123_gen(V, K, sigma, P, random_state=None):
    del random_state
    N, L = V.shape
    K = int(K)
    return _h3_threshold_with_ts(
        V,
        K,
        sigma,
        P,
        (
            _h3_t1(N, K, L),
            _h3_t2(N, K, L),
            _h3_t3(N, K, L),
        ),
    )


def _general_1swap_from(base_solver, V, K, sigma, P, random_state=None):
    x = base_solver(V, K, sigma=sigma, P=P, random_state=random_state)
    return refine_general_1swap(
        V,
        x,
        K,
        sigma=sigma,
        P=P,
        **FAST_1SWAP_LS_KWARGS,
    )


def _h3_strong_weak_1swap(V, K, sigma, P, random_state=None):
    return _general_1swap_from(solve_h3_strong_weak, V, K, sigma, P, random_state)


def _h3_threshold_t123_gen_1swap(V, K, sigma, P, random_state=None):
    return _general_1swap_from(_h3_threshold_t123_gen, V, K, sigma, P, random_state)


def _cap_window_1swap(V, K, sigma, P, random_state=None):
    return _general_1swap_from(solve_cap_window_gen, V, K, sigma, P, random_state)


def _cap_window_full_1swap(V, K, sigma, P, random_state=None):
    return _general_1swap_from(solve_cap_window_full_gen, V, K, sigma, P, random_state)


def _cap_submodular_1swap(V, K, sigma, P, random_state=None):
    return _general_1swap_from(solve_cap_submodular_gen, V, K, sigma, P, random_state)


def _threshold_specs(prefix):
    return (
        (f"{prefix}-BF", partial(solve_h3, target_obj="bf")),
        (f"{prefix}-Int", partial(solve_h3, target_obj="int")),
        (f"{prefix}-Gen", partial(solve_h3, target_obj="gen")),
    )


def _frame_specs(frame_kwargs, fixed_random_state=None):
    common = dict(frame_kwargs)
    if fixed_random_state is not None:
        common["random_state"] = fixed_random_state
    return (
        ("Frame-BF", partial(solve_frame_portfolio, target_obj="bf", **common)),
        ("Frame-Int", partial(solve_frame_portfolio, target_obj="int", **common)),
        ("Frame-Gen", partial(solve_frame_portfolio, target_obj="gen", **common)),
    )


def _frame_only_specs(frame_kwargs, fixed_random_state=None):
    common = dict(frame_kwargs)
    common["external_starts"] = False
    if fixed_random_state is not None:
        common["random_state"] = fixed_random_state
    return (
        (
            "FrameOnly-BF",
            partial(solve_frame_portfolio, target_obj="bf", **common),
        ),
        (
            "FrameOnly-Int",
            partial(solve_frame_portfolio, target_obj="int", **common),
        ),
        (
            "FrameOnly-Gen",
            partial(solve_frame_portfolio, target_obj="gen", **common),
        ),
    )


def _submodular_general_specs(include_swap=True):
    specs = [
        ("ThreshDOpt-Gen", partial(solve_thresholded_logdet_greedy, **FAST_THRESH_LOGDET_KWARGS)),
        (
            "ThreshWLogdet-Gen",
            partial(solve_thresholded_logdet_greedy, **WEIGHTED_THRESH_LOGDET_KWARGS),
        ),
    ]
    if include_swap:
        specs.append(
            (
                "ThreshDOptSwap-Gen",
                partial(solve_thresholded_logdet_greedy, **SWAP_THRESH_LOGDET_KWARGS),
            )
        )
    return tuple(specs)


BASELINE_SOLVERS = (
    ("H1", solve_h1),
    ("H2", solve_h2),
    ("TrueBackwardGreedy", solve_true_backward_greedy),
    ("CoutinoSchur-Gen", solve_coutino_schur_greedy),
    ("MISO-EE", partial(solve_miso_energy_greedy, target_margin=0.05)),
    ("Pareto-H2", solve_pareto_interference_greedy),
)

SIGMA_SWEEP_SOLVERS = (
    ("H1", solve_h1),
    ("H2", solve_h2),
    ("TrueBackwardGreedy", solve_true_backward_greedy),
    ("MISO-EE", partial(solve_miso_energy_greedy, target_margin=0.05)),
    ("Pareto-H2", solve_pareto_interference_greedy),
    *_threshold_specs("H3-threshold"),
    *_frame_specs(FAST_FRAME_KWARGS),
    ("CapWindow-Gen", solve_cap_window_gen),
    ("R2Delta-Gen", solve_r2_delta_gen),
    ("CapSubmod-Gen", solve_cap_submodular_gen),
    ("CapSubmodPort-Gen", solve_cap_submodular_portfolio_gen),
    *_submodular_general_specs(include_swap=False),
    ("H3-Fast", solve_h3_fast),
)

CDF_SOLVERS = (
    ("H1", solve_h1),
    ("H2", solve_h2),
    ("H3", solve_h3_strong_weak),
    ("TrueBackwardGreedy", solve_true_backward_greedy),
    ("CoutinoSchur-Gen", solve_coutino_schur_greedy),
    ("MISO-EE", partial(solve_miso_energy_greedy, target_margin=0.05)),
    ("Pareto-H2", solve_pareto_interference_greedy),
    *_threshold_specs("S-threshold"),
    *_frame_specs(FAST_FRAME_KWARGS),
    *_frame_only_specs(FAST_FRAME_KWARGS),
    ("CapWindow-Gen", solve_cap_window_gen),
    ("R2Delta-Gen", solve_r2_delta_gen),
    ("CapSubmod-Gen", solve_cap_submodular_gen),
    ("CapSubmodPort-Gen", solve_cap_submodular_portfolio_gen),
    *_submodular_general_specs(),
    ("N-H3-Fast", solve_h3_fast),
)

REQUESTED_GEN_SOLVERS = (
    ("H3", solve_h3_strong_weak),
    ("H3-1SwapLS-Gen", _h3_strong_weak_1swap),
    ("H3ThresholdT123-Gen", _h3_threshold_t123_gen),
    ("H3ThresholdT123-1SwapLS-Gen", _h3_threshold_t123_gen_1swap),
    ("CapWindow-Gen", solve_cap_window_gen),
    ("CapWindow-1SwapLS-Gen", _cap_window_1swap),
    ("CapWindowFull-Gen", solve_cap_window_full_gen),
    ("CapWindowFull-1SwapLS-Gen", _cap_window_full_1swap),
    ("R2Delta-Gen", solve_r2_delta_gen),
    ("CapSubmod-Gen", solve_cap_submodular_gen),
    ("CapSubmod-1SwapLS-Gen", _cap_submodular_1swap),
)

SMALL_GUROBI_HEURISTICS = (
    *BASELINE_SOLVERS,
    *_threshold_specs("H3-threshold"),
    *_frame_specs(EXACT_FRAME_KWARGS, fixed_random_state=0),
    ("CapWindow-Gen", solve_cap_window_gen),
    ("R2Delta-Gen", solve_r2_delta_gen),
    ("CapSubmod-Gen", solve_cap_submodular_gen),
    ("CapSubmodPort-Gen", solve_cap_submodular_portfolio_gen),
    *_submodular_general_specs(),
    ("H3-Fast", partial(solve_h3_fast, random_state=0)),
)
