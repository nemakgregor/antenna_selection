from algorithms import (
    solve_cap_window_gen,
    solve_coutino_greedy,
    solve_frame_portfolio,
    solve_h1,
    solve_h2,
    solve_h3,
    solve_h3_fast,
    solve_h3_strong_weak,
    solve_miso_energy_greedy,
    solve_pareto_interference_greedy,
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


def _h1(V, K, sigma, P, random_state=None):
    return solve_h1(V, K, sigma=sigma, P=P)


def _h2(V, K, sigma, P, random_state=None):
    return solve_h2(V, K, sigma=sigma, P=P)


def _h3_strong_weak(V, K, sigma, P, random_state=None):
    return solve_h3_strong_weak(V, K, sigma=sigma, P=P)


def _coutino(V, K, sigma, P, random_state=None):
    return solve_coutino_greedy(V, K, sigma=sigma, P=P)


def _miso_energy(V, K, sigma, P, random_state=None):
    return solve_miso_energy_greedy(V, K, sigma=sigma, P=P, target_margin=0.05)


def _pareto_h2(V, K, sigma, P, random_state=None):
    return solve_pareto_interference_greedy(V, K, sigma=sigma, P=P)


def _threshold_solver(target_obj):
    def solver(V, K, sigma, P, random_state=None):
        return solve_h3(V, K, target_obj=target_obj, sigma=sigma, P=P)

    return solver


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


BASELINE_SOLVERS = (
    ("H1", _h1),
    ("H2", _h2),
    ("Coutino", _coutino),
    ("MISO-EE", _miso_energy),
    ("Pareto-H2", _pareto_h2),
)

GRID_SOLVERS = (
    *BASELINE_SOLVERS,
    *_threshold_specs("S-threshold"),
    *_frame_specs(FAST_FRAME_KWARGS),
    *_frame_only_specs(FAST_FRAME_KWARGS),
    ("CapWindow-Gen", _cap_window),
    ("N-H3-Fast", _h3_fast()),
)

SIGMA_SWEEP_SOLVERS = (
    *BASELINE_SOLVERS,
    *_threshold_specs("H3-threshold"),
    *_frame_specs(FAST_FRAME_KWARGS),
    ("CapWindow-Gen", _cap_window),
    ("H3-Fast", _h3_fast()),
)

CDF_SOLVERS = (
    ("H1", _h1),
    ("H2", _h2),
    ("H3", _h3_strong_weak),
    ("Coutino", _coutino),
    ("MISO-EE", _miso_energy),
    ("Pareto-H2", _pareto_h2),
    *_threshold_specs("S-threshold"),
    *_frame_specs(FAST_FRAME_KWARGS),
    *_frame_only_specs(FAST_FRAME_KWARGS),
    ("CapWindow-Gen", _cap_window),
    ("N-H3-Fast", _h3_fast()),
)

MOTOR_SOLVERS = (
    ("H1", _h1),
    ("H2", _h2),
    ("H3", _h3_strong_weak),
    ("H3-Fast", _h3_fast()),
    *_threshold_specs("H3"),
    *_frame_specs(TEST_FRAME_KWARGS),
    ("CapWindow-Gen", _cap_window),
    ("Coutino", _coutino),
    ("MISO-EE", _miso_energy),
    ("Pareto-H2", _pareto_h2),
)

SMALL_GUROBI_HEURISTICS = (
    *BASELINE_SOLVERS,
    *_threshold_specs("H3-threshold"),
    *_frame_specs(EXACT_FRAME_KWARGS, fixed_random_state=0),
    ("CapWindow-Gen", _cap_window),
    ("H3-Fast", _h3_fast(fixed_random_state=0)),
)
