import os, sys, unittest
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import vepbench as vb


class TestBetweenGroup(unittest.TestCase):
    def test_between_is_in_unit_interval(self):
        d = vb.simulate_benchmark(n_genes=120, seed=11)
        b = vb.between_group_auroc(d, "gene", "label", "within_signal")
        self.assertTrue(0 <= b <= 1)

    def test_gene_prior_has_high_between_and_chance_within(self):
        d = vb.simulate_benchmark(n_genes=200, seed=3)
        w = vb.within_group_auroc(d, "gene", "label", "gene_prior")
        b = vb.between_group_auroc(d, "gene", "label", "gene_prior")
        self.assertAlmostEqual(w, 0.5, delta=0.05)
        self.assertGreater(b, 0.65)

    def test_within_signal_has_high_both(self):
        d = vb.simulate_benchmark(n_genes=200, seed=3)
        w = vb.within_group_auroc(d, "gene", "label", "within_signal")
        b = vb.between_group_auroc(d, "gene", "label", "within_signal")
        self.assertGreater(w, 0.7)
        self.assertGreater(b, 0.7)


class TestPairComposition(unittest.TestCase):
    def test_counts_add_up(self):
        d = vb.simulate_benchmark(n_genes=80, seed=7)
        pc = vb.pair_composition(d, "gene", "label", "random")
        self.assertAlmostEqual(pc["within_pairs"] + pc["between_pairs"], pc["total_pairs"])
        self.assertAlmostEqual(pc["within_share"] + pc["between_share"], 1.0)

    def test_share_small_for_many_groups(self):
        d = vb.simulate_benchmark(n_genes=500, seed=1)
        pc = vb.pair_composition(d, "gene", "label", "random")
        self.assertLess(pc["within_share"], 0.2)


class TestGapDiagnostic(unittest.TestCase):
    def test_within_signal_has_small_gap(self):
        d = vb.simulate_benchmark(n_genes=200, seed=3)
        g = vb.protein_mean_vs_residual_auroc(d, "gene", "label", "within_signal")
        self.assertLess(abs(g["gap_mean_minus_residual"]), 0.3)

    def test_gene_prior_has_large_gap(self):
        d = vb.simulate_benchmark(n_genes=200, seed=3)
        g = vb.protein_mean_vs_residual_auroc(d, "gene", "label", "gene_prior")
        self.assertGreater(g["gap_mean_minus_residual"], 0.15)

    def test_random_has_zero_gap(self):
        d = vb.simulate_benchmark(n_genes=200, seed=3)
        g = vb.protein_mean_vs_residual_auroc(d, "gene", "label", "random")
        self.assertLess(abs(g["gap_mean_minus_residual"]), 0.05)


if __name__ == "__main__":
    unittest.main()
