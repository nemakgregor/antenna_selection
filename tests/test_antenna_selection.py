import itertools
import subprocess
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from algorithms import (
    best_cyclic_threshold_window,
    calculate_objectives,
    check_constraints,
    cyclic_threshold_window_selection,
    solve_cap_submodular_gen,
    solve_cap_submodular_portfolio_gen,
    solve_cap_window_full_gen,
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
    refine_general_1swap,
    refine_selection_by_ug_swaps,
    refine_selection_by_ug_swaps_steps,
    threshold_window_selection,
    solve_thresholded_logdet_greedy,
    solve_true_backward_greedy,
)
from algorithms.h3_threshold_local import refine_threshold_by_swaps
from utils.brute_force import brute_force_exact_u_g, contiguous_threshold_window_T
from utils.data import generate_V

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

    def test_true_backward_greedy_logic(self):
        x = solve_true_backward_greedy(self.V_L4, self.K)
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

    def test_general_1swap_local_search_nonworse(self):
        x = solve_h3_strong_weak(self.V_L4, self.K, sigma=1.0, P=1.0)
        refined = refine_general_1swap(
            self.V_L4,
            x,
            self.K,
            sigma=1.0,
            P=1.0,
            max_passes=2,
            remove_limit=20,
            add_limit=30,
            boundary_limit=30,
        )
        is_valid, num_active = check_constraints(refined, self.K)
        self.assertTrue(is_valid)
        self.assertEqual(num_active, self.K)

        base_u_g = calculate_objectives(self.V_L4, x, sigma=1.0, P=1.0)[2]
        refined_u_g = calculate_objectives(self.V_L4, refined, sigma=1.0, P=1.0)[2]
        self.assertGreaterEqual(refined_u_g, base_u_g - 1e-8)

    def test_cyclic_threshold_window_exact_k_with_wraparound(self):
        x = cyclic_threshold_window_selection(self.V_L2, self.K, start=self.N - 5)
        is_valid, num_active = check_constraints(x, self.K)
        self.assertTrue(is_valid)
        self.assertEqual(num_active, self.K)

        best = best_cyclic_threshold_window(self.V_L2, self.K, sigma=1.0, P=1.0)
        is_valid, num_active = check_constraints(best["x"], self.K)
        self.assertTrue(is_valid)
        self.assertEqual(num_active, self.K)
        self.assertGreaterEqual(best["T"], 0)
        self.assertLess(best["T"], self.N)

    def test_fixed_threshold_window_exact_k(self):
        T = int(np.floor(0.05 * self.N + 0.5))
        x = threshold_window_selection(self.V_L2, self.K, T)
        is_valid, num_active = check_constraints(x, self.K)
        self.assertTrue(is_valid)
        self.assertEqual(num_active, self.K)

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

        self.assertTrue(exact["completed"])
        self.assertFalse(exact["timed_out"])
        self.assertEqual(exact["subset"], manual_subset)
        self.assertAlmostEqual(exact["u_g"], manual_best)
        self.assertEqual(contiguous_threshold_window_T(V, exact["subset"]), 0)

    def test_threshold_local_helper_never_decreases_u_g(self):
        seed = refine_threshold_by_swaps(
            self.V_L2,
            self.K,
            T=5,
            max_swaps=0,
            sigma=1.0,
            P=1.0,
            candidate_radius=6,
        )
        refined = refine_threshold_by_swaps(
            self.V_L2,
            self.K,
            T=5,
            max_swaps=1,
            sigma=1.0,
            P=1.0,
            candidate_radius=6,
        )
        is_valid, num_active = check_constraints(refined["x"], self.K)
        self.assertTrue(is_valid)
        self.assertEqual(num_active, self.K)
        self.assertGreaterEqual(refined["u_g"], seed["u_g"] - 1e-8)

    def test_ug_swap_zero_one_two_three_nonworse(self):
        x0 = solve_h3_strong_weak(self.V_L2, self.K, sigma=1.0, P=1.0)
        zero = refine_selection_by_ug_swaps(
            self.V_L2,
            x0,
            max_swaps=0,
            sigma=1.0,
            P=1.0,
            K=self.K,
        )
        self.assertTrue(np.array_equal(zero["x"], x0))

        steps = refine_selection_by_ug_swaps_steps(
            self.V_L2,
            x0,
            max_swaps_values=(0, 1, 2, 3),
            sigma=1.0,
            P=1.0,
            K=self.K,
        )
        self.assertEqual(set(steps), {0, 1, 2, 3})
        self.assertGreaterEqual(steps[1]["u_g"], steps[0]["u_g"] - 1e-8)
        self.assertGreaterEqual(steps[2]["u_g"], steps[1]["u_g"] - 1e-8)
        self.assertGreaterEqual(steps[3]["u_g"], steps[2]["u_g"] - 1e-8)
        for result in steps.values():
            is_valid, num_active = check_constraints(result["x"], self.K)
            self.assertTrue(is_valid)
            self.assertEqual(num_active, self.K)

    def test_direct_frame_and_cap_seeds_refine_to_exact_k(self):
        frame_seed = solve_frame_portfolio(
            self.V_L2,
            self.K,
            sigma=1.0,
            P=1.0,
            target_obj="gen",
            random_state=42,
            max_refined_starts=2,
            max_passes=0,
            remove_limit=0,
            add_limit=0,
            lambdas=(),
        )
        cap_seed = solve_cap_submodular_gen(
            self.V_L2,
            self.K,
            sigma=1.0,
            P=1.0,
            random_state=42,
        )
        for seed in (frame_seed, cap_seed):
            x = refine_selection_by_ug_swaps(
                self.V_L2,
                seed,
                max_swaps=1,
                sigma=1.0,
                P=1.0,
                K=self.K,
            )["x"]
            is_valid, num_active = check_constraints(x, self.K)
            self.assertTrue(is_valid)
            self.assertEqual(num_active, self.K)

    def test_ug_swap_seed_comparison_smoke(self):
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp) / "ug_swap_smoke"
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "experiments.algorithm_comparison",
                    "--ug-swap-seed-comparison",
                    "--N",
                    "20",
                    "--L",
                    "2",
                    "--K-values",
                    "10",
                    "15",
                    "--samples",
                    "1",
                    "--generator-seeds",
                    "42",
                    "--data-profiles",
                    "gaussian",
                    "--sigma",
                    "1",
                    "--out-dir",
                    str(out_dir),
                ],
                check=True,
            )
            self.assertTrue((out_dir / "ug_swap_seed_summary.csv").exists())
            self.assertTrue((out_dir / "ug_swap_seed_report.md").exists())
            self.assertTrue((out_dir / "ug_swap_raw_u_g_cdf.png").exists())
            archive_path = out_dir / "csv_data.tar.gz"
            self.assertTrue(archive_path.exists())
            with tarfile.open(archive_path, "r:gz") as archive:
                self.assertIn("ug_swap_seed_runs.csv", archive.getnames())

    def test_unified_local_swap_comparison_smoke(self):
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp) / "unified_local_swap"
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "experiments.algorithm_comparison",
                    "--unified-local-swap-comparison",
                    "--N",
                    "20",
                    "--L",
                    "2",
                    "--off-pcts",
                    "50",
                    "--samples",
                    "1",
                    "--generator-seeds",
                    "42",
                    "--data-profiles",
                    "gaussian",
                    "--workers",
                    "1",
                    "--out-dir",
                    str(out_dir),
                ],
                check=True,
            )
            runs_path = out_dir / "unified_local_swap_runs.csv"
            summary_path = out_dir / "unified_local_swap_summary.csv"
            self.assertTrue(runs_path.exists())
            self.assertTrue(summary_path.exists())
            self.assertTrue((out_dir / "unified_local_swap_report.md").exists())

            runs = pd.read_csv(runs_path)
            self.assertEqual(len(runs), 10)
            self.assertEqual(set(runs["max_swaps"]), {0, 1})
            self.assertEqual(runs["method"].nunique(), 10)
            self.assertTrue((runs["active_count"] == runs["K"]).all())

    def test_cyclic_best_3swap_analysis_smoke_and_plot_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline_dir = root / "baseline"
            out_dir = root / "cyclic3"
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "experiments.algorithm_comparison",
                    "--ug-swap-seed-comparison",
                    "--N",
                    "20",
                    "--L",
                    "2",
                    "--K-values",
                    "10",
                    "15",
                    "--samples",
                    "1",
                    "--generator-seeds",
                    "42",
                    "--data-profiles",
                    "gaussian",
                    "--sigma",
                    "1",
                    "--out-dir",
                    str(baseline_dir),
                ],
                check=True,
            )
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "experiments.algorithm_comparison",
                    "--cyclic-best-3swap-analysis",
                    "--N",
                    "20",
                    "--L",
                    "2",
                    "--K-values",
                    "10",
                    "15",
                    "--samples",
                    "1",
                    "--generator-seeds",
                    "42",
                    "--data-profiles",
                    "gaussian",
                    "--sigma",
                    "1",
                    "--baseline-dir",
                    str(baseline_dir),
                    "--out-dir",
                    str(out_dir),
                ],
                check=True,
            )
            self.assertTrue((out_dir / "cyclic_best_3swap_summary.csv").exists())
            self.assertTrue((out_dir / "cyclic_best_3swap_improvement.csv").exists())
            self.assertTrue((out_dir / "combined_previous_vs_cyclic3_summary.csv").exists())
            self.assertTrue((out_dir / "cyclic_best_3swap_report.md").exists())
            self.assertTrue((out_dir / "cyclic_best_3swap_raw_u_g_cdf.png").exists())
            archive_path = out_dir / "csv_data.tar.gz"
            self.assertTrue(archive_path.exists())
            with tarfile.open(archive_path, "r:gz") as archive:
                self.assertIn("cyclic_best_3swap_runs.csv", archive.getnames())

            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "experiments.algorithm_comparison",
                    "--cyclic-best-3swap-analysis",
                    "--plot-only",
                    "--N",
                    "20",
                    "--L",
                    "2",
                    "--K-values",
                    "10",
                    "15",
                    "--samples",
                    "1",
                    "--generator-seeds",
                    "42",
                    "--data-profiles",
                    "gaussian",
                    "--sigma",
                    "1",
                    "--baseline-dir",
                    str(baseline_dir),
                    "--out-dir",
                    str(out_dir),
                ],
                check=True,
            )

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

if __name__ == "__main__":
    unittest.main()
