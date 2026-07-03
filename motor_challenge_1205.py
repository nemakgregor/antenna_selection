import argparse
import itertools
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from algorithms import (
    calculate_objectives,
    check_constraints,
    solve_backward_true_greedy,
    solve_cap_submodular_gen,
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
    solve_thresholded_logdet_greedy,
)
from algorithms.h3_threshold_explore import (
    dense_thresholds,
    evaluate_power_window_thresholds,
    evaluate_threshold_T,
    legacy_thresholds,
    threshold_power_window,
)
from algorithms.h3_threshold_local import (
    evaluate_threshold_local_rules,
    refine_threshold_by_swaps,
    threshold_window_selection,
)
from utils.solver_sets import MOTOR_SOLVERS
from utils.data import DATA_PROFILES, generate_V, generate_v_profile_from_rng
from utils.evaluation import evaluate_algorithms as evaluate_solver_set
from experiments.algorithm_comparison import (
    build_active_pct_cases,
    build_full_sweep_best_thresholds,
    build_scaling_strong_weak_runs,
    collect_threshold_scaling_profile_runs,
    collect_threshold_exact_study_runs,
    threshold_formula_T,
)
from utils.brute_force import brute_force_exact_u_g, contiguous_threshold_window_T
from utils.local_threshold_analysis import run_local_threshold_exact_analysis


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

    def test_backward_true_greedy_logic(self):
        x = solve_backward_true_greedy(self.V_L4, self.K)
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

    def test_h3_threshold_explore_explicit_T(self):
        result = evaluate_threshold_T(
            self.V_L4,
            self.K,
            T=5,
            target_obj="gen",
            sigma=1.0,
            P=1.0,
        )
        x = result["x"]
        is_valid, num_active = check_constraints(x, self.K)
        self.assertTrue(is_valid)
        self.assertEqual(num_active, self.K)
        self.assertTrue(np.isfinite([result["u_bf"], result["u_i"], result["u_g"]]).all())

    def test_h3_threshold_explore_top_window_T0(self):
        x = threshold_power_window(self.V_L4, self.K, T=0)
        row_power = np.sum(np.abs(self.V_L4) ** 2, axis=1).real
        expected_on = np.argsort(row_power)[::-1][: self.K]
        actual_on = np.flatnonzero(x == 1)
        self.assertTrue(np.array_equal(np.sort(actual_on), np.sort(expected_on)))

    def test_h3_threshold_explore_legacy_thresholds(self):
        self.assertEqual(legacy_thresholds(100, 50), [1, 2, 5, 10, 25, 50])
        self.assertEqual(legacy_thresholds(100, 80), [1, 2, 5, 10])

    def test_h3_threshold_full_sweep_dense_thresholds(self):
        self.assertEqual(dense_thresholds(250), list(range(251)))
        self.assertEqual(dense_thresholds(500), list(range(501)))

    def test_h3_threshold_full_sweep_power_windows(self):
        rows = evaluate_power_window_thresholds(
            self.V_L2,
            self.K,
            dense_thresholds(self.K),
            sigma=1.0,
            P=1.0,
        )
        self.assertEqual(len(rows), self.K + 1)
        for row in rows:
            self.assertEqual(row["active_count"], self.K)
            self.assertEqual(row["candidate_kind"], "power_window")
            self.assertTrue(np.isfinite([row["u_bf"], row["u_i"], row["u_g"]]).all())

    def test_brute_force_exact_solver_returns_valid_exact_k(self):
        V = np.array(
            [
                [1.0 + 0j, 0.0 + 0j],
                [0.0 + 0j, 1.0 + 0j],
                [0.4 + 0j, 0.2 + 0j],
                [0.1 + 0j, 0.3 + 0j],
            ]
        )
        result = brute_force_exact_u_g(V, 2, sigma=1.0, P=1.0)
        self.assertTrue(result["completed"])
        self.assertFalse(result["timed_out"])
        self.assertEqual(int(np.sum(result["x"])), 2)
        self.assertEqual(len(result["subset"]), 2)
        self.assertEqual(result["candidate_count"], 6)
        self.assertEqual(result["evaluated_count"], 6)
        self.assertTrue(np.isfinite([result["u_bf"], result["u_i"], result["u_g"]]).all())

    def test_brute_force_exact_solver_matches_manual_enumeration(self):
        V = np.array(
            [
                [1.0 + 0j, 0.0 + 0j],
                [0.0 + 0j, 1.0 + 0j],
                [0.1 + 0j, 0.0 + 0j],
                [0.0 + 0j, 0.1 + 0j],
            ]
        )
        exact = brute_force_exact_u_g(V, 2, sigma=1.0, P=1.0)
        manual_best = -np.inf
        manual_subset = None
        for subset in itertools.combinations(range(4), 2):
            x = np.zeros(4, dtype=int)
            x[list(subset)] = 1
            u_g = calculate_objectives(V, x, sigma=1.0, P=1.0)[2]
            if u_g > manual_best:
                manual_best = u_g
                manual_subset = subset
        self.assertEqual(exact["subset"], manual_subset)
        self.assertAlmostEqual(exact["u_g"], manual_best)
        self.assertEqual(exact["subset"], (0, 1))

    def test_brute_force_exact_solver_timeout_guard_marks_partial(self):
        rng = np.random.RandomState(123)
        V = generate_v_profile_from_rng(rng, 10, 2, profile="gaussian")
        result = brute_force_exact_u_g(
            V,
            5,
            sigma=1.0,
            P=1.0,
            time_limit_seconds=0.0,
            timeout_check_interval=1,
        )
        self.assertFalse(result["completed"])
        self.assertTrue(result["timed_out"])
        self.assertGreater(result["evaluated_count"], 0)
        self.assertLess(result["evaluated_count"], result["candidate_count"])
        self.assertEqual(int(np.sum(result["x"])), 5)

    def test_contiguous_threshold_window_detector(self):
        V = np.array(
            [
                [2.0 + 0j, 0.0 + 0j],
                [0.0 + 0j, np.sqrt(3.0) + 0j],
                [np.sqrt(2.0) + 0j, 0.0 + 0j],
                [0.0 + 0j, 1.0 + 0j],
            ]
        )
        self.assertEqual(contiguous_threshold_window_T(V, (0, 1)), 0)
        self.assertIsNone(contiguous_threshold_window_T(V, (0, 2)))

    def test_threshold_scaling_active_pct_cases(self):
        cases = build_active_pct_cases(1000, [25, 50])
        self.assertEqual([case["K"] for case in cases], [250, 500])
        self.assertEqual([case["active_pct"] for case in cases], [25.0, 50.0])
        self.assertEqual([case["off_pct"] for case in cases], [75.0, 50.0])

    def test_threshold_scaling_formula_clipping(self):
        self.assertEqual(threshold_formula_T("T_0p10N", 5000, 4, 250), 250)
        self.assertEqual(threshold_formula_T("T_0p05K", 1000, 8, 250), 12)
        self.assertEqual(threshold_formula_T("legacy_T100", 1000, 4, 50), 50)

    def test_threshold_scaling_tiny_shard_rows(self):
        args = SimpleNamespace(
            K_pcts=[25, 50],
            generator_seeds=[42],
            samples=1,
            sigma=1.0,
            P=1.0,
        )
        runs = collect_threshold_scaling_profile_runs(args, 20, 4, "gaussian")
        self.assertEqual(set(runs["K"]), {5, 10})
        self.assertEqual(len(runs[runs["K"] == 5]), 6)
        self.assertEqual(len(runs[runs["K"] == 10]), 11)
        self.assertTrue((runs["active_count"] == runs["K"]).all())
        self.assertTrue(np.isfinite(runs[["u_bf", "u_i", "u_g"]].to_numpy()).all())
        self.assertTrue(np.isfinite(runs["fraction_best_tested_u_g"]).all())

    def test_threshold_scaling_strong_weak_rows(self):
        args = SimpleNamespace(
            K_pcts=[25, 50],
            generator_seeds=[42],
            samples=1,
            sigma=1.0,
            P=1.0,
        )
        runs = collect_threshold_scaling_profile_runs(args, 20, 4, "gaussian")
        best = build_full_sweep_best_thresholds(runs)
        strong_weak = build_scaling_strong_weak_runs(best)
        self.assertEqual(set(strong_weak["K"]), {5, 10})
        self.assertEqual(len(strong_weak), 2)
        self.assertTrue((strong_weak["formula"] == "strong_weak").all())
        self.assertTrue(np.isfinite(strong_weak["fraction_best_tested_u_g"]).all())
        self.assertTrue((strong_weak["fraction_best_tested_u_g"] >= 0.0).all())

    def test_threshold_exact_study_tiny_rows(self):
        args = SimpleNamespace(
            N_values=[8],
            L=2,
            K_pcts=[25, 50],
            data_profiles=["gaussian"],
            generator_seeds=[42],
            samples=1,
            sigma=1.0,
            P=1.0,
            exact_time_limit=5.0,
        )
        threshold_runs, exact_runs, formula_runs = collect_threshold_exact_study_runs(args)
        self.assertEqual(set(exact_runs["K"]), {2, 4})
        self.assertTrue(exact_runs["exact_completed"].astype(bool).all())
        self.assertTrue(np.isfinite(exact_runs["best_tested_fraction_exact_u_g"]).all())
        self.assertTrue((exact_runs["best_tested_fraction_exact_u_g"] <= 1.0 + 1e-9).all())
        self.assertEqual(set(threshold_runs["K"]), {2, 4})
        self.assertTrue((threshold_runs["active_count"] == threshold_runs["K"]).all())
        self.assertTrue(np.isfinite(formula_runs["fraction_exact_u_g"]).all())

    def test_h3_threshold_local_exact_k_and_zero_swaps(self):
        x = threshold_window_selection(self.V_L2, self.K, T=5)
        result = refine_threshold_by_swaps(
            self.V_L2,
            self.K,
            T=5,
            max_swaps=0,
            sigma=1.0,
            P=1.0,
        )
        self.assertTrue(np.array_equal(result["x"], x))
        self.assertEqual(int(np.sum(result["x"])), self.K)
        self.assertEqual(result["swaps_applied"], 0)
        self.assertTrue(np.isfinite([result["u_bf"], result["u_i"], result["u_g"]]).all())

    def test_h3_threshold_local_swaps_never_decrease(self):
        rows = evaluate_threshold_local_rules(
            self.V_L2,
            self.K,
            seed_rules=[("best_tested_T", 5)],
            max_swaps_values=(0, 1, 2),
            sigma=1.0,
            P=1.0,
        )
        by_swaps = {row["max_swaps"]: row for row in rows}
        self.assertGreaterEqual(by_swaps[1]["u_g"], by_swaps[0]["u_g"] - 1e-9)
        self.assertGreaterEqual(by_swaps[2]["u_g"], by_swaps[1]["u_g"] - 1e-9)
        for row in rows:
            self.assertEqual(row["active_count"], self.K)

    def test_h3_threshold_local_best_T_improves_or_equals(self):
        rows = evaluate_power_window_thresholds(
            self.V_L2,
            self.K,
            dense_thresholds(self.K),
            sigma=1.0,
            P=1.0,
        )
        best_T = max(rows, key=lambda row: row["u_g"])["T"]
        pure = refine_threshold_by_swaps(
            self.V_L2,
            self.K,
            T=best_T,
            max_swaps=0,
            sigma=1.0,
            P=1.0,
        )
        local = refine_threshold_by_swaps(
            self.V_L2,
            self.K,
            T=best_T,
            max_swaps=2,
            sigma=1.0,
            P=1.0,
        )
        self.assertGreaterEqual(local["u_g"], pure["u_g"] - 1e-9)

    def test_h3_threshold_local_constructed_one_swap_reaches_exact(self):
        V = np.array(
            [
                [np.sqrt(10.0) + 0j, 0.0 + 0j],
                [np.sqrt(9.0) + 0j, 0.0 + 0j],
                [0.0 + 0j, np.sqrt(8.0) + 0j],
                [0.0 + 0j, np.sqrt(0.1) + 0j],
            ]
        )
        exact = brute_force_exact_u_g(V, 2, sigma=1.0, P=1.0)
        pure = refine_threshold_by_swaps(V, 2, T=0, max_swaps=0, sigma=1.0, P=1.0)
        local = refine_threshold_by_swaps(V, 2, T=0, max_swaps=1, sigma=1.0, P=1.0)
        self.assertLess(pure["u_g"], exact["u_g"])
        self.assertEqual(local["subset"], exact["subset"])
        self.assertAlmostEqual(local["u_g"], exact["u_g"])

    def test_h3_threshold_local_boundary_clipping(self):
        for T in (0, self.K, self.N + 100):
            result = refine_threshold_by_swaps(
                self.V_L2,
                self.K,
                T=T,
                max_swaps=1,
                sigma=1.0,
                P=1.0,
            )
            self.assertEqual(int(np.sum(result["x"])), self.K)
            self.assertTrue(np.isfinite([result["u_bf"], result["u_i"], result["u_g"]]).all())

    def test_threshold_local_exact_analysis_tiny_output(self):
        args = SimpleNamespace(
            N_values=[8],
            L=2,
            K_pcts=[25, 50],
            data_profiles=["gaussian"],
            generator_seeds=[42],
            samples=1,
            sigma=1.0,
            P=1.0,
            exact_time_limit=5.0,
        )
        _threshold_runs, exact_runs, _formula_runs = collect_threshold_exact_study_runs(args)
        with tempfile.TemporaryDirectory() as tmpdir:
            exact_dir = Path(tmpdir) / "exact"
            out_dir = Path(tmpdir) / "local"
            exact_dir.mkdir()
            exact_runs.to_csv(exact_dir / "exact_runs.csv", index=False)
            outputs = run_local_threshold_exact_analysis(
                exact_dir,
                out_dir,
                n_values=[8],
                k_pcts=[25, 50],
                profiles=["gaussian"],
            )
            runs = outputs["local_runs"]
            self.assertFalse(runs.empty)
            self.assertEqual(set(runs["max_swaps"]), {0, 1, 2})
            self.assertTrue((runs["active_count"] == runs["K"]).all())
            self.assertTrue(np.isfinite(runs["fraction_exact_u_g"]).all())
            self.assertTrue((out_dir / "local_threshold_exact_gauss_report.md").exists())

    def test_distribution_profile_generators(self):
        self.assertIn("twdp", DATA_PROFILES)
        for profile in DATA_PROFILES:
            rng = np.random.RandomState(7)
            V = generate_v_profile_from_rng(rng, 40, 2, profile=profile)
            self.assertEqual(V.shape, (40, 2))
            self.assertTrue(np.iscomplexobj(V))
            self.assertTrue(np.isfinite(V.real).all())
            self.assertTrue(np.isfinite(V.imag).all())

    def test_twdp_distribution_profile_is_deterministic_and_nonzero(self):
        rng_a = np.random.RandomState(42)
        rng_b = np.random.RandomState(42)
        V_a = generate_v_profile_from_rng(rng_a, 40, 2, profile="twdp")
        V_b = generate_v_profile_from_rng(rng_b, 40, 2, profile="twdp")
        self.assertTrue(np.allclose(V_a, V_b))
        self.assertTrue(np.isfinite(V_a).all())
        self.assertGreater(np.linalg.norm(V_a), 0.0)
        self.assertGreater(np.min(np.linalg.norm(V_a, axis=0)), 0.0)

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

    def test_cap_submodular_gen_logic(self):
        x = solve_cap_submodular_gen(
            self.V_L4,
            self.K,
            sigma=1.0,
            P=1.0,
            random_state=42,
            scan_size=5,
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
            eps_values=(1e-6, 1e-2),
            threshold_scan_size=5,
        )
        is_valid, num_active = check_constraints(x, self.K)
        self.assertTrue(is_valid)
        self.assertEqual(num_active, self.K)

    def test_thresholded_weighted_logdet_logic(self):
        x = solve_thresholded_logdet_greedy(
            self.V_L4,
            self.K,
            sigma=1.0,
            P=1.0,
            eps_values=(1e-6,),
            lambdas=(1.0, 0.7),
            threshold_scan_size=3,
        )
        is_valid, num_active = check_constraints(x, self.K)
        self.assertTrue(is_valid)
        self.assertEqual(num_active, self.K)

    def test_thresholded_logdet_swap_logic(self):
        x = solve_thresholded_logdet_greedy(
            self.V_L4,
            self.K,
            sigma=1.0,
            P=1.0,
            eps_values=(1e-6,),
            threshold_scan_size=3,
            swap_max_passes=1,
            swap_remove_limit=20,
            swap_add_limit=20,
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
