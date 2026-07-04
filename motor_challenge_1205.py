import argparse
import sys
import unittest

import numpy as np

from algorithms import (
    calculate_objectives,
    check_constraints,
    solve_cap_submodular_gen,
    solve_cap_submodular_portfolio_gen,
    solve_cap_window_full_gen,
    solve_coutino_greedy,
    solve_coutino_schur_greedy,
    solve_cap_window_gen,
    solve_frame_portfolio,
    solve_h1,
    solve_h2,
    solve_h3,
    solve_h3_fast,
    solve_h3_strong_weak,
    solve_miso_energy_greedy,
    solve_pareto_interference_greedy,
    solve_thresholded_logdet_greedy,
)
from utils.solver_sets import MOTOR_SOLVERS
from utils.data import generate_V
from utils.evaluation import evaluate_algorithms as evaluate_solver_set


def evaluate_algorithms(V, K, sigma=1.0, P=1.0, random_state=None):
    return evaluate_solver_set(
        V,
        K,
        sigma=sigma,
        P=P,
        random_state=random_state,
        solvers=MOTOR_SOLVERS,
    )


# ==========================================
# 2. Unit Tests
# ==========================================

class TestAntennaSelection(unittest.TestCase):
    def setUp(self):
        np.random.seed(42)
        self.N, self.K = 100, 50
        self.V_L2 = generate_V(self.N, 2)
        self.V_L4 = generate_V(self.N, 4)

    def test_generate_data(self):
        self.assertEqual(self.V_L2.shape, (self.N, 2))
        self.assertTrue(np.iscomplexobj(self.V_L2))

    def test_h1_logic(self):
        x = solve_h1(self.V_L2, self.K)
        is_valid, num_active = check_constraints(x, self.K)
        self.assertTrue(is_valid)
        self.assertEqual(num_active, self.K)

        row_power = np.sum(np.abs(self.V_L2) ** 2, axis=1).real
        weakest_first = np.argsort(row_power)
        expected_off = weakest_first[: self.N - self.K]
        actual_off = np.flatnonzero(x == 0)
        self.assertTrue(np.array_equal(np.sort(actual_off), np.sort(expected_off)))

    def test_h2_logic_L2(self):
        x = solve_h2(self.V_L2, self.K)
        is_valid, num_active = check_constraints(x, self.K)
        self.assertTrue(is_valid)
        self.assertEqual(num_active, self.K)

    def test_h2_logic_generalized_L4(self):
        x = solve_h2(self.V_L4, self.K)
        is_valid, num_active = check_constraints(x, self.K)
        self.assertTrue(is_valid)
        self.assertEqual(num_active, self.K)

    def test_coutino_greedy_logic(self):
        x = solve_coutino_greedy(self.V_L4, self.K)
        is_valid, num_active = check_constraints(x, self.K)
        self.assertTrue(is_valid)
        self.assertLessEqual(num_active, self.K)
        self.assertGreaterEqual(num_active, self.V_L4.shape[1])

    def test_coutino_schur_logic(self):
        x = solve_coutino_schur_greedy(
            self.V_L4,
            self.K,
            sigma=1.0,
            P=1.0,
        )
        is_valid, num_active = check_constraints(x, self.K)
        self.assertTrue(is_valid)
        self.assertEqual(num_active, self.K)

    def test_miso_energy_greedy_logic(self):
        x = solve_miso_energy_greedy(self.V_L4, self.K)
        is_valid, num_active = check_constraints(x, self.K)
        self.assertTrue(is_valid)
        self.assertGreaterEqual(num_active, self.V_L4.shape[1])

    def test_pareto_interference_greedy_logic(self):
        x = solve_pareto_interference_greedy(self.V_L4, self.K)
        is_valid, num_active = check_constraints(x, self.K)
        self.assertTrue(is_valid)
        self.assertGreaterEqual(num_active, self.V_L4.shape[1])

    def test_h3_fast_logic(self):
        x_l2 = solve_h3_fast(
            self.V_L2,
            self.K,
            refinement_iter=50,
            beam_size=20,
            random_state=42,
        )
        x_l4 = solve_h3_fast(
            self.V_L4,
            self.K,
            refinement_iter=50,
            beam_size=20,
            random_state=42,
        )

        for x in (x_l2, x_l4):
            is_valid, num_active = check_constraints(x, self.K)
            self.assertTrue(is_valid)
            self.assertEqual(num_active, self.K)

    def test_h3_strong_weak_logic(self):
        x = solve_h3_strong_weak(self.V_L4, self.K)
        is_valid, num_active = check_constraints(x, self.K)
        self.assertTrue(is_valid)
        self.assertEqual(num_active, self.K)

        row_power = np.sum(np.abs(self.V_L4) ** 2, axis=1).real
        power_order = np.argsort(row_power)
        off_count = self.N - self.K
        weak_drop = off_count // 2
        strong_drop = off_count - weak_drop
        expected_off = np.r_[power_order[:weak_drop], power_order[self.N - strong_drop :]]
        actual_off = np.flatnonzero(x == 0)
        self.assertTrue(np.array_equal(np.sort(actual_off), np.sort(expected_off)))

    def test_h3_threshold_objective_modes(self):
        for target_obj in ("bf", "int", "gen"):
            x = solve_h3(
                self.V_L4,
                self.K,
                target_obj=target_obj,
                sigma=1.0,
                P=1.0,
            )
            is_valid, num_active = check_constraints(x, self.K)
            self.assertTrue(is_valid)
            self.assertEqual(num_active, self.K)

    def test_frame_portfolio_logic(self):
        variants = (
            ("bf", True),
            ("int", True),
            ("gen", True),
            ("bf", False),
            ("int", False),
            ("gen", False),
        )
        for target_obj, external_starts in variants:
            x = solve_frame_portfolio(
                self.V_L4,
                self.K,
                target_obj=target_obj,
                external_starts=external_starts,
                random_state=42,
                max_refined_starts=2,
                max_passes=1,
                remove_limit=20,
                add_limit=20,
                lambdas=(),
            )
            is_valid, num_active = check_constraints(x, self.K)
            self.assertTrue(is_valid)
            self.assertEqual(num_active, self.K)

    def test_cap_window_gen_logic(self):
        x = solve_cap_window_gen(
            self.V_L4,
            self.K,
            sigma=1.0,
            P=1.0,
            random_state=42,
        )
        is_valid, num_active = check_constraints(x, self.K)
        self.assertTrue(is_valid)
        self.assertEqual(num_active, self.K)

    def test_cap_window_full_gen_logic(self):
        x = solve_cap_window_full_gen(
            self.V_L4,
            self.K,
            sigma=1.0,
            P=1.0,
            random_state=42,
        )
        is_valid, num_active = check_constraints(x, self.K)
        self.assertTrue(is_valid)
        self.assertEqual(num_active, self.K)

    def test_cap_submodular_gen_logic(self):
        x = solve_cap_submodular_gen(
            self.V_L4,
            self.K,
            sigma=1.0,
            P=1.0,
            random_state=42,
        )
        is_valid, num_active = check_constraints(x, self.K)
        self.assertTrue(is_valid)
        self.assertEqual(num_active, self.K)

    def test_cap_submodular_portfolio_gen_logic(self):
        x = solve_cap_submodular_portfolio_gen(
            self.V_L4,
            self.K,
            sigma=1.0,
            P=1.0,
            random_state=42,
        )
        is_valid, num_active = check_constraints(x, self.K)
        self.assertTrue(is_valid)
        self.assertEqual(num_active, self.K)

    def test_thresholded_logdet_logic(self):
        x = solve_thresholded_logdet_greedy(
            self.V_L4,
            self.K,
            sigma=1.0,
            P=1.0,
            random_state=42,
            threshold_scan_size=3,
            eps_values=(1e-6,),
            lambdas=(1.0,),
            swap_max_passes=1,
        )
        is_valid, num_active = check_constraints(x, self.K)
        self.assertTrue(is_valid)
        self.assertEqual(num_active, self.K)

    def test_frame_portfolio_gen_includes_h3_start(self):
        h3_x = solve_h3_strong_weak(self.V_L4, self.K, sigma=1.0, P=1.0)
        frame_x = solve_frame_portfolio(
            self.V_L4,
            self.K,
            target_obj="gen",
            external_starts=False,
            random_state=42,
            max_refined_starts=2,
            max_passes=1,
            remove_limit=20,
            add_limit=20,
            lambdas=(),
        )

        h3_u_g = calculate_objectives(self.V_L4, h3_x, sigma=1.0, P=1.0)[2]
        frame_u_g = calculate_objectives(self.V_L4, frame_x, sigma=1.0, P=1.0)[2]
        self.assertGreaterEqual(frame_u_g, h3_u_g - 1e-9)

    def test_objectives_calculation(self):
        x = solve_h1(self.V_L2, self.K)
        u_bf, u_i, u_g = calculate_objectives(self.V_L2, x, sigma=1.0)
        self.assertIsInstance(u_bf, float)
        self.assertIsInstance(u_i, float)
        self.assertIsInstance(u_g, float)
        self.assertTrue(u_bf >= 0)
        self.assertTrue(u_i >= 0)


# ==========================================
# 3. Main Execution
# ==========================================

def parse_args():
    parser = argparse.ArgumentParser(description="Antenna Selection Optimization")
    parser.add_argument("--N", type=int, default=1000, help="Number of antennas (default: 1000)")
    parser.add_argument("--L", type=int, default=2, help="Number of layers/streams (default: 2)")
    parser.add_argument("--K", type=int, default=500, help="Max active antennas limit (default: 500)")
    parser.add_argument("--sigma", type=float, default=1.0, help="Noise power (default: 1.0)")
    parser.add_argument("--P", type=float, default=1.0, help="Given Power parameter (default: 1.0)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for generation (default: 42)")
    parser.add_argument("--samples", type=int, default=100, help="Random examples per case (default: 100)")
    parser.add_argument("--test", action="store_true", help="Run the unit tests and exit")
    return parser.parse_args()


def print_aggregate_results(aggregate, win_counts, samples):
    metric_specs = [
        ("u_bf", "BF Gain", "max"),
        ("u_i", "Interference", "min"),
        ("u_g", "General Obj", "max"),
    ]
    algorithm_names = list(aggregate)

    print("--- Average Objective Values ---")
    print(
        f"{'Algorithm':<16} | {'Valid':>7} | {'Active':>9} | "
        f"{'BF mean':>13} | {'Int mean':>13} | {'Gen mean':>13}"
    )
    print("-" * 78)

    for name in algorithm_names:
        row = aggregate[name]
        print(
            f"{name:<16} | {row['valid_count'] / samples:>6.1%} | "
            f"{row['active_count'] / samples:>9.2f} | "
            f"{row['u_bf'] / samples:>13.6g} | "
            f"{row['u_i'] / samples:>13.6g} | "
            f"{row['u_g'] / samples:>13.6g}"
        )

    print("\n--- Objective Comparison Over Samples ---")
    print(f"{'Metric':<22} | {'Best by mean':<16} | {'Per-sample wins'}")
    print("-" * 78)

    for metric, label, direction in metric_specs:
        means = {name: aggregate[name][metric] / samples for name in algorithm_names}
        best = max(means, key=means.get) if direction == "max" else min(means, key=means.get)
        wins = ", ".join(
            f"{name}:{win_counts[metric].get(name, 0)}"
            for name in algorithm_names
            if win_counts[metric].get(name, 0)
        )
        print(f"{label:<22} | {best:<16} | {wins}")


def main():
    args = parse_args()

    if args.test:
        sys.argv = [sys.argv[0]]
        unittest.main(exit=True)

    if args.N <= 0:
        raise ValueError("--N must be positive.")
    if args.L <= 0:
        raise ValueError("--L must be positive.")
    if not (0 <= args.K <= args.N):
        raise ValueError("--K must satisfy 0 <= K <= N.")
    if args.samples <= 0:
        raise ValueError("--samples must be positive.")

    np.random.seed(args.seed)
    print("--- Initialization ---")
    print(
        f"Generating {args.samples} V matrices "
        f"(N={args.N}, L={args.L}) | Seed: {args.seed} | Sigma: {args.sigma}"
    )
    print(f"Targeting K = {args.K} active antennas ({(args.K/args.N)*100:.1f}% active).\n")

    metric_specs = [
        ("u_bf", "max"),
        ("u_i", "min"),
        ("u_g", "max"),
    ]
    aggregate = {}
    win_counts = {metric: {} for metric, _ in metric_specs}
    progress_step = max(1, args.samples // 10)

    for sample_idx in range(args.samples):
        V_matrix = generate_V(args.N, args.L)
        sample_results = evaluate_algorithms(
            V_matrix,
            args.K,
            sigma=args.sigma,
            P=args.P,
            random_state=args.seed + sample_idx,
        )

        if not aggregate:
            aggregate = {
                name: {
                    "valid_count": 0,
                    "active_count": 0.0,
                    "u_bf": 0.0,
                    "u_i": 0.0,
                    "u_g": 0.0,
                }
                for name in sample_results
            }

        for name, result in sample_results.items():
            aggregate[name]["valid_count"] += int(result["valid"])
            aggregate[name]["active_count"] += result["active_count"]
            aggregate[name]["u_bf"] += result["u_bf"]
            aggregate[name]["u_i"] += result["u_i"]
            aggregate[name]["u_g"] += result["u_g"]

        for metric, direction in metric_specs:
            values = {
                name: result[metric]
                for name, result in sample_results.items()
                if result["valid"]
            }
            best = max(values, key=values.get) if direction == "max" else min(values, key=values.get)
            win_counts[metric][best] = win_counts[metric].get(best, 0) + 1

        if args.samples > 1 and (
            (sample_idx + 1) % progress_step == 0 or sample_idx + 1 == args.samples
        ):
            print(f"  completed {sample_idx + 1}/{args.samples}")

    print()
    print_aggregate_results(aggregate, win_counts, args.samples)


if __name__ == "__main__":
    main()
