import os, sys, json, tempfile, subprocess, unittest
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from sklearn.metrics import roc_auc_score, average_precision_score
import vepbench as vb
from vepbench.metrics import roc_points


class TestMetrics(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(1)
        self.y = rng.integers(0, 2, 500)
        self.s = np.round(rng.normal(self.y * 0.8, 1.0), 1)   # rounding creates ties

    def test_auroc_matches_sklearn_with_ties(self):
        self.assertAlmostEqual(vb.auroc(self.y, self.s), roc_auc_score(self.y, self.s), places=12)

    def test_auroc_perfect_and_reversed_and_constant(self):
        y = np.array([0, 0, 1, 1])
        self.assertEqual(vb.auroc(y, [0.1, 0.2, 0.8, 0.9]), 1.0)
        self.assertEqual(vb.auroc(y, [0.9, 0.8, 0.2, 0.1]), 0.0)
        self.assertEqual(vb.auroc(y, [1, 1, 1, 1]), 0.5)

    def test_single_class_gives_nan(self):
        self.assertTrue(np.isnan(vb.auroc([1, 1, 1], [0.1, 0.2, 0.3])))
        self.assertTrue(np.isnan(vb.partial_auroc([1, 1, 1], [0.1, 0.2, 0.3])))   # undefined, not an error

    def test_missing_scores_are_dropped(self):
        y = np.array([0, 1, 0, 1, 1]); s = np.array([0.1, 0.9, np.nan, 0.8, np.nan])
        self.assertEqual(vb.auroc(y, s), 1.0)

    def test_bad_labels_raise(self):
        with self.assertRaises(ValueError):
            vb.auroc([0, 2, 1], [0.1, 0.2, 0.3])

    def test_average_precision_matches_sklearn(self):
        self.assertAlmostEqual(vb.average_precision(self.y, self.s),
                               average_precision_score(self.y, self.s), places=12)

    def test_partial_auroc_matches_sklearn(self):
        for m in (0.05, 0.1, 0.3):
            self.assertAlmostEqual(vb.partial_auroc(self.y, self.s, max_fpr=m),
                                   roc_auc_score(self.y, self.s, max_fpr=m), places=10)

    def test_roc_points_monotone(self):
        fpr, tpr = roc_points(self.y, self.s)
        self.assertTrue(np.all(np.diff(fpr) >= 0) and np.all(np.diff(tpr) >= 0))
        self.assertEqual((fpr[0], tpr[0], fpr[-1], tpr[-1]), (0, 0, 1, 1))

    def test_sensitivity_at_specificity_properties(self):
        for spec in (0.8, 0.9, 0.99):
            sens = vb.sensitivity_at_specificity(self.y, self.s, spec)
            thr_neg = np.quantile(self.s[self.y == 0], spec, method="higher")
            realised_spec = np.mean(self.s[self.y == 0] <= thr_neg)
            self.assertGreaterEqual(realised_spec, spec)
            self.assertTrue(0 <= sens <= 1)
        self.assertGreaterEqual(vb.sensitivity_at_specificity(self.y, self.s, 0.8),
                                vb.sensitivity_at_specificity(self.y, self.s, 0.99))

    def test_specificity_at_sensitivity_properties(self):
        sp = vb.specificity_at_sensitivity(self.y, self.s, 0.9)
        thr = np.quantile(self.s[self.y == 1], 0.1, method="lower")
        self.assertGreaterEqual(np.mean(self.s[self.y == 1] >= thr), 0.9)
        self.assertTrue(0 <= sp <= 1)


class TestGeneLevel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.df = vb.simulate_benchmark(n_genes=300, seed=3)
        cls.cols = ["within_signal", "gene_prior", "mixed", "random"]

    def test_per_group_counts_and_eligibility(self):
        pg = vb.per_group_auroc(self.df, "gene", "label", "within_signal", min_per_class=5)
        self.assertEqual(int(pg["n_pos"].sum() + pg["n_neg"].sum()), len(self.df))
        el = pg["auroc"].notna()
        self.assertTrue(((pg.loc[el, "n_pos"] >= 5) & (pg.loc[el, "n_neg"] >= 5)).all())
        self.assertTrue(((pg.loc[~el, "n_pos"] < 5) | (pg.loc[~el, "n_neg"] < 5)).all())

    def test_per_group_matches_manual(self):
        pg = vb.per_group_auroc(self.df, "gene", "label", "mixed", min_per_class=5).set_index("group")
        g = pg.index[pg["auroc"].notna()][0]
        sub = self.df[self.df.gene == g]
        self.assertAlmostEqual(pg.loc[g, "auroc"], roc_auc_score(sub.label, sub.mixed), places=12)

    def test_between_gene_predictor_inflated_when_pooled(self):
        agg = vb.compare_aggregations(self.df, "gene", "label", self.cols).set_index("model")
        self.assertGreater(agg.loc["gene_prior", "pooled"], 0.65)
        self.assertAlmostEqual(agg.loc["gene_prior", "pergroup"], 0.5, delta=0.03)
        self.assertAlmostEqual(agg.loc["random", "pooled"], 0.5, delta=0.03)
        self.assertGreater(agg.loc["within_signal", "pergroup"], 0.70)

    def test_gap_decomposition_and_selection_effect(self):
        agg = vb.compare_aggregations(self.df, "gene", "label", self.cols).set_index("model")
        lhs = agg["gap_selection"] + agg["gap_aggregation"]
        self.assertTrue(np.allclose(lhs, agg["gap_pooled_minus_pergroup"], atol=1e-12))
        g = agg.loc["gene_prior"]
        self.assertGreater(g["pooled"], g["pooled_eligible"])            # extreme genes inflate the all-variant pooled value
        self.assertGreater(g["pooled_eligible"], g["pergroup"])          # pooling across genes still inflates it
        # pooled_eligible equals a direct computation on the eligible genes
        pg = vb.per_group_auroc(self.df, "gene", "label", "gene_prior", 5)
        eligible = set(pg.loc[pg.auroc.notna(), "group"])
        sub = self.df[self.df.gene.isin(eligible)]
        self.assertAlmostEqual(g["pooled_eligible"], roc_auc_score(sub.label, sub.gene_prior), places=12)

    def test_rank_change_when_aggregation_changes(self):
        agg = vb.compare_aggregations(self.df, "gene", "label", self.cols).set_index("model")
        self.assertGreater(agg.loc["gene_prior", "rank_pergroup"], agg.loc["gene_prior", "rank_pooled"])
        cmp_ = vb.compare_rankings(agg["pooled"], agg["pergroup"])
        self.assertLess(cmp_["kendall"], 1.0)
        self.assertGreater(cmp_["n_rank_changed"], 0)

    def test_gene_identity_baseline_schemes(self):
        rep = vb.gene_identity_report(self.df, "gene", "label").set_index("scheme")
        self.assertGreater(rep.loc["kfold", "pooled"], 0.65)            # random split: gene identity works
        self.assertGreaterEqual(rep.loc["in_sample", "pooled"], rep.loc["kfold", "pooled"] - 1e-9)
        self.assertAlmostEqual(rep.loc["group_kfold", "pooled"], 0.5, delta=0.04)  # no gene information
        self.assertAlmostEqual(rep.loc["in_sample", "pergroup"], 0.5, places=12)   # constant within a gene
        self.assertAlmostEqual(rep.loc["group_kfold", "pergroup"], 0.5, places=12)
        self.assertTrue(np.isnan(rep.loc["kfold", "pergroup"]) and np.isnan(rep.loc["loo", "pergroup"]))
        # documented artefact: out-of-fold gene rates are biased BELOW 0.5 within genes
        self.assertLess(rep.loc["kfold", "pergroup_biased"], 0.45)
        self.assertLess(rep.loc["loo", "pergroup_biased"], 0.45)

    def test_gene_identity_unknown_scheme(self):
        with self.assertRaises(ValueError):
            vb.gene_identity_baseline(self.df, "gene", "label", scheme="nope")

    def test_gene_identity_is_reproducible(self):
        a = vb.gene_identity_baseline(self.df, "gene", "label", "kfold", seed=5)
        b = vb.gene_identity_baseline(self.df, "gene", "label", "kfold", seed=5)
        self.assertTrue(np.array_equal(a, b))


class TestBootstrap(unittest.TestCase):
    def test_percentile_ci(self):
        lo, hi = vb.percentile_ci(np.arange(101), alpha=0.1)
        self.assertAlmostEqual(lo, 5); self.assertAlmostEqual(hi, 95)

    def test_cluster_bootstrap_se_matches_replicate_sd(self):
        """Bootstrap SE of the mean per-gene AUROC should match the SD over fresh datasets."""
        reps = []
        for k in range(60):
            d = vb.simulate_benchmark(n_genes=60, seed=100 + k)
            pg = vb.per_group_auroc(d, "gene", "label", "within_signal")
            reps.append(vb.summarize_per_group(pg)["mean"])
        true_sd = np.std(reps, ddof=1)
        d = vb.simulate_benchmark(n_genes=60, seed=100)
        pg = vb.per_group_auroc(d, "gene", "label", "within_signal")
        boots = vb.bootstrap_pergroup_mean(pg["auroc"].to_numpy(), n_boot=2000, seed=1)
        ratio = np.std(boots, ddof=1) / true_sd
        self.assertTrue(0.7 < ratio < 1.4, ratio)

    def test_pooled_cluster_bootstrap_wider_than_naive(self):
        """Resampling genes should give a wider interval than resampling variants for a between-gene predictor."""
        from vepbench.bootstrap import bootstrap_pooled_auroc
        d = vb.simulate_benchmark(n_genes=80, seed=7)
        cl = bootstrap_pooled_auroc(d, "gene", "label", "gene_prior", n_boot=300, seed=0)
        rng = np.random.default_rng(0)
        y, s = d["label"].to_numpy(), d["gene_prior"].to_numpy()
        naive = [vb.auroc(y[i], s[i]) for i in (rng.integers(0, len(d), len(d)) for _ in range(300))]
        self.assertGreater(np.std(cl), 1.3 * np.std(naive))

    def test_paired_difference(self):
        d = vb.simulate_benchmark(n_genes=150, seed=9)
        pa = vb.per_group_auroc(d, "gene", "label", "within_signal")
        pr = vb.per_group_auroc(d, "gene", "label", "random")
        res = vb.paired_bootstrap_pergroup(pa, pr, n_boot=1000, seed=0)
        self.assertGreater(res["ci_low"], 0); self.assertLess(res["p_value"], 0.01)
        same = vb.paired_bootstrap_pergroup(pa, pa, n_boot=500, seed=0)
        self.assertEqual(same["mean_diff"], 0.0)
        self.assertEqual(same["ci_low"], 0.0); self.assertEqual(same["ci_high"], 0.0)

    def test_bootstrap_statistic_runs(self):
        d = vb.simulate_benchmark(n_genes=30, seed=2)
        out = vb.bootstrap_statistic(d, "gene", lambda x: vb.auroc(x["label"], x["mixed"]), n_boot=20, seed=0)
        self.assertEqual(out.shape, (20,)); self.assertFalse(np.isnan(out).any())


class TestRanking(unittest.TestCase):
    def test_identical_and_reversed(self):
        a = pd.Series([0.9, 0.8, 0.7, 0.6], index=list("wxyz"))
        same = vb.compare_rankings(a, a)
        self.assertEqual(same["kendall"], 1.0); self.assertEqual(same["n_rank_changed"], 0)
        rev = vb.compare_rankings(a, a[::-1].set_axis(a.index))
        self.assertEqual(rev["kendall"], -1.0); self.assertEqual(rev["discordant_fraction"], 1.0)

    def test_rank_models_ties(self):
        r = vb.rank_models(pd.Series([3.0, 3.0, 1.0], index=list("abc")))
        self.assertEqual(r["a"], 1.5); self.assertEqual(r["c"], 3.0)

    def test_bootstrap_rank_agreement_with_external(self):
        d = vb.simulate_benchmark(n_genes=120, seed=11)
        cols = ["within_signal", "gene_prior", "mixed", "random"]
        # independent 'assay' ranking that rewards within-gene signal
        ext = pd.Series({"within_signal": 0.9, "mixed": 0.7, "gene_prior": 0.5, "random": 0.4})
        samples, summ = vb.bootstrap_rank_agreement(d, "gene", "label", cols, n_boot=60, seed=0, external=ext)
        self.assertEqual(len(samples), 60)
        self.assertGreater(summ["external_gain_pergroup_minus_pooled"]["mean"], 0)
        self.assertIn("kendall_pooled_vs_pergroup", summ)


class TestIOAndCLI(unittest.TestCase):
    def test_load_table_label_mapping_and_flip(self):
        with tempfile.TemporaryDirectory() as t:
            p = os.path.join(t, "x.csv")
            pd.DataFrame(dict(g=["a", "a", "b"], lab=["P", "B", "P"], s1=[1., 2., 3.], s2=[1., 2., 3.])).to_csv(p, index=False)
            df = vb_io().load_table(p, "g", "lab", ["s1", "s2"], positive_label="P", flip=["s2"])
            self.assertEqual(df["lab"].tolist(), [1, 0, 1])
            self.assertEqual(df["s2"].tolist(), [-1., -2., -3.])
            with self.assertRaises(KeyError):
                vb_io().load_table(p, "g", "lab", ["nope"], positive_label="P")
            with self.assertRaises(ValueError):
                vb_io().load_table(p, "g", "lab", ["s1"])   # labels not 0/1 and no positive_label

    def test_merge_group_files(self):
        with tempfile.TemporaryDirectory() as t:
            L, S = os.path.join(t, "L"), os.path.join(t, "S"); os.makedirs(L); os.makedirs(S)
            pd.DataFrame(dict(mutant=["A1B", "A2C"], DMS_score_bin=[1, 0])).to_csv(os.path.join(L, "P1.csv"), index=False)
            pd.DataFrame(dict(mutant=["A1B", "A2C"], m1=[0.9, 0.1])).to_csv(os.path.join(S, "P1.csv"), index=False)
            pd.DataFrame(dict(mutant=["A1B"], DMS_score_bin=[1])).to_csv(os.path.join(L, "P2.csv"), index=False)
            df, skipped = vb_io().merge_group_files(L, S)
            self.assertEqual(skipped, ["P2"]); self.assertEqual(df["group"].unique().tolist(), ["P1"])
            self.assertEqual(df["label"].tolist(), [1, 0])

    def test_cli_end_to_end(self):
        with tempfile.TemporaryDirectory() as t:
            csv = os.path.join(t, "sim.csv"); out = os.path.join(t, "out")
            d = vb.simulate_benchmark(n_genes=60, seed=4); d.to_csv(csv, index=False)
            ext = os.path.join(t, "ext.csv")
            pd.DataFrame(dict(model=["within_signal", "gene_prior", "mixed", "random"], value=[.9, .5, .7, .4])).to_csv(ext, index=False)
            root = os.path.join(os.path.dirname(__file__), "..")
            r = subprocess.run([sys.executable, "-m", "vepbench", "evaluate", "--input", csv, "--group", "gene",
                                "--label", "label", "--scores", "within_signal,gene_prior,mixed,random",
                                "--n-boot", "30", "--external", ext, "--out", out],
                               cwd=root, capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            for f in ("models.csv", "per_group_auroc.csv", "aggregation_comparison.csv",
                      "gene_identity_baseline.csv", "summary.json", "report.md"):
                self.assertTrue(os.path.exists(os.path.join(out, f)), f)
            m = pd.read_csv(os.path.join(out, "models.csv")).set_index("model")
            self.assertGreater(m.loc["gene_prior", "pooled_auroc"], m.loc["gene_prior", "pergroup_mean_auroc"])
            self.assertTrue((m["pooled_ci_low"] <= m["pooled_auroc"]).all())
            with open(os.path.join(out, "summary.json")) as fh:
                s = json.load(fh)
            self.assertIn("ranking_vs_external", s)


def vb_io():
    import vepbench.io as m
    return m


if __name__ == "__main__":
    unittest.main()
