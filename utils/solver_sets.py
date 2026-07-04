import numpy as np

from algorithms import (
    solve_backward_true_greedy,
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
    solve_thresholded_logdet_greedy,
)


FAST_FRAME_KWARGS = {
    "max_refined_starts": 3,
    "max_passes": 2,
    "remove_limit": 60,
    "add_limit": 60,
    "lambdas": (),
}

TEST_FRAME_KWARGS = {
    "max_refined_starts": 3,
    "max_passes": 2,
    "remove_limit": 50,
    "add_limit": 50,
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


def _h1(V, K, sigma, P, random_state=None):
    return solve_h1(V, K, sigma=sigma, P=P)


def _h2(V, K, sigma, P, random_state=None):
    return solve_h2(V, K, sigma=sigma, P=P)


def _h3_strong_weak(V, K, sigma, P, random_state=None):
    return solve_h3_strong_weak(V, K, sigma=sigma, P=P)


def _backward_true_greedy(V, K, sigma, P, random_state=None):
    return solve_backward_true_greedy(V, K, sigma=sigma, P=P)


def _coutino_schur(V, K, sigma, P, random_state=None):
    return solve_coutino_schur_greedy(
        V,
        K,
        sigma=sigma,
        P=P,
        random_state=random_state,
    )


def _miso_energy(V, K, sigma, P, random_state=None):
    return solve_miso_energy_greedy(V, K, sigma=sigma, P=P, target_margin=0.05)


def _pareto_h2(V, K, sigma, P, random_state=None):
    return solve_pareto_interference_greedy(V, K, sigma=sigma, P=P)


def _threshold_solver(target_obj):
    def solver(V, K, sigma, P, random_state=None):
        return solve_h3(V, K, target_obj=target_obj, sigma=sigma, P=P)

    return solver


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


def _frame_solver(target_obj, frame_kwargs, external_starts=True, fixed_random_state=None):
    def solver(V, K, sigma, P, random_state=None):
        seed = fixed_random_state if fixed_random_state is not None else random_state
        return solve_frame_portfolio(
            V,
            K,
            target_obj=target_obj,
            sigma=sigma,
            P=P,
            random_state=seed,
            external_starts=external_starts,
            **frame_kwargs,
        )

    return solver


def _cap_window(V, K, sigma, P, random_state=None):
    return solve_cap_window_gen(V, K, sigma=sigma, P=P, random_state=random_state)


def _cap_window_full(V, K, sigma, P, random_state=None):
    return solve_cap_window_full_gen(V, K, sigma=sigma, P=P, random_state=random_state)


def _cap_submodular(V, K, sigma, P, random_state=None):
    return solve_cap_submodular_gen(
        V,
        K,
        sigma=sigma,
        P=P,
        random_state=random_state,
    )


def _cap_submodular_portfolio(V, K, sigma, P, random_state=None):
    return solve_cap_submodular_portfolio_gen(
        V,
        K,
        sigma=sigma,
        P=P,
        random_state=random_state,
    )


def _general_1swap_from(base_solver, V, K, sigma, P, random_state=None):
    x = base_solver(V, K, sigma, P, random_state=random_state)
    return refine_general_1swap(
        V,
        x,
        K,
        sigma=sigma,
        P=P,
        **FAST_1SWAP_LS_KWARGS,
    )


def _h3_strong_weak_1swap(V, K, sigma, P, random_state=None):
    return _general_1swap_from(_h3_strong_weak, V, K, sigma, P, random_state)


def _h3_threshold_t123_gen_1swap(V, K, sigma, P, random_state=None):
    return _general_1swap_from(_h3_threshold_t123_gen, V, K, sigma, P, random_state)


def _cap_window_1swap(V, K, sigma, P, random_state=None):
    return _general_1swap_from(_cap_window, V, K, sigma, P, random_state)


def _cap_window_full_1swap(V, K, sigma, P, random_state=None):
    return _general_1swap_from(_cap_window_full, V, K, sigma, P, random_state)


def _cap_submodular_1swap(V, K, sigma, P, random_state=None):
    return _general_1swap_from(_cap_submodular, V, K, sigma, P, random_state)


def _thresholded_logdet_solver(solver_kwargs):
    def solver(V, K, sigma, P, random_state=None):
        return solve_thresholded_logdet_greedy(
            V,
            K,
            sigma=sigma,
            P=P,
            random_state=random_state,
            **solver_kwargs,
        )

    return solver


def _h3_fast(fixed_random_state=None):
    def solver(V, K, sigma, P, random_state=None):
        seed = fixed_random_state if fixed_random_state is not None else random_state
        return solve_h3_fast(V, K, random_state=seed)

    return solver


def _threshold_specs(prefix):
    return (
        (f"{prefix}-BF", _threshold_solver("bf")),
        (f"{prefix}-Int", _threshold_solver("int")),
        (f"{prefix}-Gen", _threshold_solver("gen")),
    )


def _frame_specs(frame_kwargs, fixed_random_state=None):
    return (
        ("Frame-BF", _frame_solver("bf", frame_kwargs, fixed_random_state=fixed_random_state)),
        ("Frame-Int", _frame_solver("int", frame_kwargs, fixed_random_state=fixed_random_state)),
        ("Frame-Gen", _frame_solver("gen", frame_kwargs, fixed_random_state=fixed_random_state)),
    )


def _frame_only_specs(frame_kwargs, fixed_random_state=None):
    return (
        (
            "FrameOnly-BF",
            _frame_solver(
                "bf",
                frame_kwargs,
                external_starts=False,
                fixed_random_state=fixed_random_state,
            ),
        ),
        (
            "FrameOnly-Int",
            _frame_solver(
                "int",
                frame_kwargs,
                external_starts=False,
                fixed_random_state=fixed_random_state,
            ),
        ),
        (
            "FrameOnly-Gen",
            _frame_solver(
                "gen",
                frame_kwargs,
                external_starts=False,
                fixed_random_state=fixed_random_state,
            ),
        ),
    )


def _submodular_general_specs(include_swap=True):
    specs = [
        ("ThreshDOpt-Gen", _thresholded_logdet_solver(FAST_THRESH_LOGDET_KWARGS)),
        (
            "ThreshWLogdet-Gen",
            _thresholded_logdet_solver(WEIGHTED_THRESH_LOGDET_KWARGS),
        ),
    ]
    if include_swap:
        specs.append(
            (
                "ThreshDOptSwap-Gen",
                _thresholded_logdet_solver(SWAP_THRESH_LOGDET_KWARGS),
            )
        )
    return tuple(specs)


BASELINE_SOLVERS = (
    ("H1", _h1),
    ("H2", _h2),
    ("BackwardTrueGreedy", _backward_true_greedy),
    ("CoutinoSchur-Gen", _coutino_schur),
    ("MISO-EE", _miso_energy),
    ("Pareto-H2", _pareto_h2),
)

GRID_SOLVERS = (
    *BASELINE_SOLVERS,
    *_threshold_specs("S-threshold"),
    *_frame_specs(FAST_FRAME_KWARGS),
    *_frame_only_specs(FAST_FRAME_KWARGS),
    ("CapWindow-Gen", _cap_window),
    ("CapSubmod-Gen", _cap_submodular),
    ("CapSubmodPort-Gen", _cap_submodular_portfolio),
    *_submodular_general_specs(),
    ("N-H3-Fast", _h3_fast()),
)

SIGMA_SWEEP_SOLVERS = (
    ("H1", _h1),
    ("H2", _h2),
    ("BackwardTrueGreedy", _backward_true_greedy),
    ("MISO-EE", _miso_energy),
    ("Pareto-H2", _pareto_h2),
    *_threshold_specs("H3-threshold"),
    *_frame_specs(FAST_FRAME_KWARGS),
    ("CapWindow-Gen", _cap_window),
    ("CapSubmod-Gen", _cap_submodular),
    ("CapSubmodPort-Gen", _cap_submodular_portfolio),
    *_submodular_general_specs(include_swap=False),
    ("H3-Fast", _h3_fast()),
)

CDF_SOLVERS = (
    ("H1", _h1),
    ("H2", _h2),
    ("H3", _h3_strong_weak),
    ("BackwardTrueGreedy", _backward_true_greedy),
    ("CoutinoSchur-Gen", _coutino_schur),
    ("MISO-EE", _miso_energy),
    ("Pareto-H2", _pareto_h2),
    *_threshold_specs("S-threshold"),
    *_frame_specs(FAST_FRAME_KWARGS),
    *_frame_only_specs(FAST_FRAME_KWARGS),
    ("CapWindow-Gen", _cap_window),
    ("CapSubmod-Gen", _cap_submodular),
    ("CapSubmodPort-Gen", _cap_submodular_portfolio),
    *_submodular_general_specs(),
    ("N-H3-Fast", _h3_fast()),
)

REQUESTED_GEN_SOLVERS = (
    ("H3", _h3_strong_weak),
    ("H3-1SwapLS-Gen", _h3_strong_weak_1swap),
    ("H3ThresholdT123-Gen", _h3_threshold_t123_gen),
    ("H3ThresholdT123-1SwapLS-Gen", _h3_threshold_t123_gen_1swap),
    ("CapWindow-Gen", _cap_window),
    ("CapWindow-1SwapLS-Gen", _cap_window_1swap),
    ("CapWindowFull-Gen", _cap_window_full),
    ("CapWindowFull-1SwapLS-Gen", _cap_window_full_1swap),
    ("CapSubmod-Gen", _cap_submodular),
    ("CapSubmod-1SwapLS-Gen", _cap_submodular_1swap),
)

MOTOR_SOLVERS = (
    ("H1", _h1),
    ("H2", _h2),
    ("H3", _h3_strong_weak),
    ("H3-Fast", _h3_fast()),
    *_threshold_specs("H3"),
    *_frame_specs(TEST_FRAME_KWARGS),
    ("CapWindow-Gen", _cap_window),
    ("CapSubmod-Gen", _cap_submodular),
    ("CapSubmodPort-Gen", _cap_submodular_portfolio),
    *_submodular_general_specs(),
    ("BackwardTrueGreedy", _backward_true_greedy),
    ("CoutinoSchur-Gen", _coutino_schur),
    ("MISO-EE", _miso_energy),
    ("Pareto-H2", _pareto_h2),
)

SMALL_GUROBI_HEURISTICS = (
    *BASELINE_SOLVERS,
    *_threshold_specs("H3-threshold"),
    *_frame_specs(EXACT_FRAME_KWARGS, fixed_random_state=0),
    ("CapWindow-Gen", _cap_window),
    ("CapSubmod-Gen", _cap_submodular),
    ("CapSubmodPort-Gen", _cap_submodular_portfolio),
    *_submodular_general_specs(),
    ("H3-Fast", _h3_fast(fixed_random_state=0)),
)
