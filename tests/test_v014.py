import os, sys, itertools, tempfile, unittest
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import vepbench as vb
import importlib.util
spec = importlib.util.spec_from_file_location("rec", os.path.join(os.path.dirname(__file__), "..", "analysis", "reconcile_with_first_run.py"))
rec = importlib.util.module_from_spec(spec); spec.loader.exec_module(rec)


def brute(df, g, y, s, within):
    num = den = 0.0
    P = df[df[y] == 1]; N = df[df[y] == 0]
    for (_, p), (_, n) in itertools.product(P.iterrows(), N.iterrows()):
        if (p[g] == n[g]) != within:
            continue
        den += 1; num += 1.0 if p[s] > n[s] else (0.5 if p[s] == n[s] else 0.0)
    return num / den if den else float("nan")


class TestDecomposition(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(11)
        n = 150
        self.df = pd.DataFrame(dict(g=rng.integers(0, 12, n), y=rng.integers(0, 2, n)))
        self.df["s"] = np.round(rng.normal(self.df.y * 0.8 + 0.3 * self.df.g / 12, 1.0), 1)   # ties, some gene signal
        self.df.loc[self.df.g == 11, "y"] = 0                                                   # a single-class group

    def test_between_matches_brute_force(self):
        self.assertAlmostEqual(vb.between_group_auroc(self.df, "g", "y", "s"), brute(self.df, "g", "y", "s", False), places=12)

    def test_within_matches_brute_force(self):
        self.assertAlmostEqual(vb.within_group_auroc(self.df, "g", "y", "s"), brute(self.df, "g", "y", "s", True), places=12)

    def test_pair_counts_and_identity(self):
        c = vb.pair_composition(self.df, "g", "y", "s")
        P = self.df[self.df.y == 1]; N = self.df[self.df.y == 0]
        self.assertEqual(c["total_pairs"], len(P) * len(N))
        self.assertEqual(c["within_pairs"] + c["between_pairs"], c["total_pairs"])
        self.assertAlmostEqual(c["within_share"] + c["between_share"], 1.0, places=12)
        pooled = vb.auroc(self.df.y, self.df.s)
        w = vb.within_group_auroc(self.df, "g", "y", "s"); b = vb.between_group_auroc(self.df, "g", "y", "s")
        self.assertAlmostEqual((w * c["within_pairs"] + b * c["between_pairs"]) / c["total_pairs"], pooled, places=12)

    def test_missing_scores_are_handled_consistently(self):
        d = self.df.copy(); d.loc[d.index[::7], "s"] = np.nan
        c = vb.pair_composition(d, "g", "y", "s"); k = d.dropna(subset=["s"])
        self.assertEqual(c["total_pairs"], (k.y == 1).sum() * (k.y == 0).sum())
        self.assertAlmostEqual(vb.between_group_auroc(d, "g", "y", "s"), brute(k, "g", "y", "s", False), places=12)

    def test_no_between_pairs_is_nan(self):
        d = pd.DataFrame(dict(g=[1, 1, 2, 2], y=[1, 0, 1, 0], s=[.9, .1, .8, .2]))
        d = d[d.g == 1]
        self.assertTrue(np.isnan(vb.between_group_auroc(d, "g", "y", "s")))

    def test_within_share_is_tiny_with_many_small_groups(self):
        d = vb.simulate_benchmark(n_genes=300, seed=4)
        self.assertLess(vb.pair_composition(d, "gene", "label", "mixed")["within_share"], 0.02)

    def test_gap_separates_gene_tracking_from_within_signal(self):
        d = vb.simulate_benchmark(n_genes=300, seed=3)
        gap = {c: vb.between_minus_within(d, "gene", "label", c)["gap_between_minus_within"] for c in ["within_signal", "gene_prior", "mixed", "random"]}
        self.assertAlmostEqual(gap["within_signal"], 0.0, delta=0.03)
        self.assertAlmostEqual(gap["random"], 0.0, delta=0.03)
        self.assertGreater(gap["gene_prior"], 0.15)
        self.assertGreater(gap["mixed"], 0.05); self.assertLess(gap["mixed"], gap["gene_prior"])

    def test_negative_gap_when_scores_are_not_comparable_across_groups(self):
        rng = np.random.default_rng(2)
        rows = []
        for g in range(150):
            y = rng.integers(0, 2, 30); shift = rng.normal(0, 3)              # a group-specific scale shift, no gene information
            rows.append(pd.DataFrame(dict(g=g, y=y, s=1.5 * y + rng.normal(0, 1, 30) + shift)))
        d = pd.concat(rows)
        r = vb.between_minus_within(d, "g", "y", "s")
        self.assertLess(r["gap_between_minus_within"], -0.1)

    def test_protein_mean_gap_is_not_a_gene_reliance_diagnostic(self):
        """Documents why the protein-mean versus residual AUROC gap is NOT used: a predictor with no gene information
        still has a high protein-mean AUROC, because proteins with more pathogenic variants get higher mean scores."""
        d = vb.simulate_benchmark(n_genes=300, seed=3)
        d["mean_s"] = d.groupby("gene")["within_signal"].transform("mean")
        d["resid_s"] = d["within_signal"] - d["mean_s"]
        mean_auc, resid_auc = vb.auroc(d.label, d.mean_s), vb.auroc(d.label, d.resid_s)
        self.assertGreater(mean_auc - resid_auc, 0.1)                               # large 'gap' for a predictor built without gene information
        self.assertAlmostEqual(vb.between_minus_within(d, "gene", "label", "within_signal")["gap_between_minus_within"], 0.0, delta=0.03)

    def test_bootstrap_decomposition(self):
        d = vb.simulate_benchmark(n_genes=80, seed=6)
        pt = vb.between_minus_within(d, "gene", "label", "mixed")
        bs = vb.bootstrap_decomposition(d, "gene", "label", "mixed", n_boot=300, seed=0)
        self.assertEqual(list(bs.columns), ["pooled", "within", "between", "gap_between_minus_within"])
        lo, hi = vb.percentile_ci(bs["gap_between_minus_within"].to_numpy())
        self.assertTrue(lo <= pt["gap_between_minus_within"] <= hi)
        lo, hi = vb.percentile_ci(bs["pooled"].to_numpy()); self.assertTrue(lo <= vb.auroc(d.label, d.mixed) <= hi)

    def test_evaluate_reports_decomposition_columns(self):
        from vepbench.cli import evaluate
        d = vb.simulate_benchmark(n_genes=60, seed=9)
        r = evaluate(d, "gene", "label", ["within_signal", "gene_prior"], n_boot=20)["models"].set_index("model")
        for col in ("within_pairwise_auroc", "between_pairwise_auroc", "gap_between_minus_within", "within_pair_share",
                    "gap_between_minus_within_ci_low", "gap_between_minus_within_ci_high"):
            self.assertIn(col, r.columns)
        self.assertGreater(r.loc["gene_prior", "gap_between_minus_within"], r.loc["within_signal", "gap_between_minus_within"])


class TestReconcile(unittest.TestCase):
    def test_reconcile_flags_a_changed_model(self):
        from vepbench.cli import evaluate
        d = vb.simulate_benchmark(n_genes=150, seed=10)
        cols = ["within_signal", "gene_prior", "mixed", "random"]
        res = evaluate(d, "gene", "label", cols, n_boot=0)
        with tempfile.TemporaryDirectory() as t:
            pp = res["per_group"].rename(columns={"group": "group"}); pp.to_csv(os.path.join(t, "pp.csv"), index=False)
            old = res["models"]; old.to_csv(os.path.join(t, "old.csv"), index=False)
            new = old.copy(); new.to_csv(os.path.join(t, "new.csv"), index=False)
            r = rec.main(["--first-run-per-protein", os.path.join(t, "pp.csv"), "--first-run-models", os.path.join(t, "old.csv"),
                          "--new-models", os.path.join(t, "new.csv"), "--tol", "0.05"])
            self.assertEqual((r["flag"] == "CHECK").sum(), 0)
            new.loc[new.model == "mixed", "within_pairwise_auroc"] += 0.2
            new.to_csv(os.path.join(t, "new2.csv"), index=False)
            r2 = rec.main(["--first-run-per-protein", os.path.join(t, "pp.csv"), "--first-run-models", os.path.join(t, "old.csv"),
                           "--new-models", os.path.join(t, "new2.csv"), "--tol", "0.05"])
            self.assertEqual(r2.loc["mixed", "flag"], "CHECK"); self.assertEqual((r2["flag"] == "CHECK").sum(), 1)


if __name__ == "__main__":
    unittest.main()
