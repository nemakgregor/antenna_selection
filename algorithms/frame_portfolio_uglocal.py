import numpy as np

from .common import calculate_objectives
from .frame_portfolio import (
    DEFAULT_LAMBDAS,
    _FrameContext,
    _active_to_x,
    _build_starts,
    _repair_active,
    _score_active,
    _top_power_start,
    solve_frame_portfolio,
)
from .ug_swap_local import refine_selection_by_ug_swaps


def frame_portfolio_seed_gen(
    V,
    K,
    sigma=1.0,
    P=1.0,
    *,
    random_state=None,
    max_refined_starts=6,
    lambdas=DEFAULT_LAMBDAS,
    random_restarts=0,
    external_starts=True,
    include_power_band_start=True,
):
    return solve_frame_portfolio(
        V,
        K,
        sigma=sigma,
        P=P,
        target_obj="gen",
        random_state=random_state,
        max_refined_starts=max_refined_starts,
        max_passes=0,
        remove_limit=0,
        add_limit=0,
        lambdas=lambdas,
        random_restarts=random_restarts,
        external_starts=external_starts,
        include_power_band_start=include_power_band_start,
    )


def solve_frame_portfolio_uglocal_gen(
    V,
    K,
    sigma=1.0,
    P=1.0,
    *,
    random_state=None,
    max_swaps=1,
    max_refined_starts=6,
    lambdas=DEFAULT_LAMBDAS,
    random_restarts=0,
    external_starts=True,
    include_power_band_start=True,
    h3_fast_refinement_iter=None,
    h3_fast_beam_size=None,
    power_band_dopt_lambdas=(1e-2,),
    power_band_dopt_buffers=(50, 100),
    power_band_dopt_offset_margin=50,
):
    context = _FrameContext(V, K, sigma=sigma, P=P)
    if context.K == 0:
        return np.zeros(context.N, dtype=int)
    if context.K == context.N:
        return np.ones(context.N, dtype=int)

    rng = np.random.default_rng(0 if random_state is None else int(random_state))
    starts = _build_starts(
        context,
        target_obj="gen",
        rng=rng,
        lambdas=lambdas,
        random_restarts=random_restarts,
        use_h3_fast_int_start=True,
        h3_fast_refinement_iter=h3_fast_refinement_iter,
        h3_fast_beam_size=h3_fast_beam_size,
        external_starts=external_starts,
        include_power_band_start=include_power_band_start,
        power_band_dopt_lambdas=power_band_dopt_lambdas,
        power_band_dopt_buffers=power_band_dopt_buffers,
        power_band_dopt_offset_margin=power_band_dopt_offset_margin,
    )

    scored_starts = []
    seen = set()
    for active in starts:
        active = _repair_active(context, active)
        key = np.packbits(active).tobytes()
        if key in seen:
            continue
        seen.add(key)
        score = _score_active(context, active, "gen")
        if np.isfinite(score):
            scored_starts.append((score, active))

    if not scored_starts:
        return _active_to_x(_top_power_start(context))

    scored_starts.sort(key=lambda item: item[0], reverse=True)
    best_x = _active_to_x(scored_starts[0][1])
    best_score = calculate_objectives(V, best_x, sigma=sigma, P=P)[2]

    for _, active in scored_starts[: max(1, int(max_refined_starts))]:
        x0 = _active_to_x(active)
        improved = refine_selection_by_ug_swaps(
            V,
            x0,
            max_swaps=max_swaps,
            sigma=sigma,
            P=P,
            K=K,
        )["x"]
        score = calculate_objectives(V, improved, sigma=sigma, P=P)[2]
        if score > best_score + 1e-12:
            best_score = score
            best_x = improved

    return best_x
