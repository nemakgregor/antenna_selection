from .cap_submodular import solve_cap_submodular_gen
from .ug_swap_local import refine_selection_by_ug_swaps


def cap_submodular_seed_gen(
    V,
    K,
    sigma=1.0,
    P=1.0,
    *,
    random_state=None,
):
    return solve_cap_submodular_gen(
        V,
        K,
        sigma=sigma,
        P=P,
        random_state=random_state,
    )


def solve_cap_submodular_uglocal_gen(
    V,
    K,
    sigma=1.0,
    P=1.0,
    *,
    random_state=None,
    max_swaps=1,
):
    x0 = cap_submodular_seed_gen(
        V,
        K,
        sigma=sigma,
        P=P,
        random_state=random_state,
    )
    return refine_selection_by_ug_swaps(
        V,
        x0,
        max_swaps=max_swaps,
        sigma=sigma,
        P=P,
        K=K,
    )["x"]
