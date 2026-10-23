import importlib.util, json, os, sys, tempfile, unittest
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("run_clinical_analysis", os.path.join(HERE, "..", "analysis", "run_clinical_analysis.py"))
rca = importlib.util.module_from_spec(spec); spec.loader.exec_module(rca)


def build_mock(root, seed=0):
    """Mock folder layout with deliberate traps."""
    rng = np.random.default_rng(seed)
    L, S = os.path.join(root, "labels"), os.path.join(root, "scores")
    os.makedirs(L)
    models = ["ModelGood", "ModelLower", "ModelAmbiguous", "ModelBroken", "ModelPartial", "ModelBetween"]
    for m in models:
        os.makedirs(os.path.join(S, m))
    for i in range(40):
        prot = "P%02d" % i
        n = int(rng.integers(30, 61))
        prior = rng.beta(0.6, 0.6)
        y = (rng.random(n) < prior).astype(int)
        mut = ["A%dC" % k for k in range(1, n + 1)]
        lab = pd.DataFrame(dict(mutant=mut, mutated_sequence=["x"] * n, DMS_score_bin=y))
        if i % 2 == 0:
            lab.insert(0, "protein", prot)          # half the files carry a protein column, half rely on the file name
        lab.to_csv(os.path.join(L, prot + ".csv"), index=False)
        good = 1.2 * y + rng.standard_normal(n)
        base = dict(protein=prot, mutant=mut)
        pd.DataFrame({**base, "note": ["n"] * n, "ModelGood_score": good}).to_csv(os.path.join(S, "ModelGood", prot + ".csv"), index=False)
        pd.DataFrame({**base, "raw": -(1.2 * y + rng.standard_normal(n))}).to_csv(os.path.join(S, "ModelLower", prot + ".csv"), index=False)
        pd.DataFrame({**base, "a": 1.0 * y + rng.standard_normal(n), "b": rng.standard_normal(n)}).to_csv(os.path.join(S, "ModelAmbiguous", prot + ".csv"), index=False)
        pd.DataFrame(dict(protein=prot, position=np.arange(n), score=rng.standard_normal(n))).to_csv(os.path.join(S, "ModelBroken", prot + ".csv"), index=False)
        keep = rng.random(n) < 0.6
        part = pd.DataFrame({**base, "s": 1.2 * y + rng.standard_normal(n)})
        part[keep].to_csv(os.path.join(S, "ModelPartial", prot + ".csv"), index=False)
        pd.DataFrame({**base, "s": np.log(prior / (1 - prior)) + rng.standard_normal(n)}).to_csv(os.path.join(S, "ModelBetween", prot + ".csv"), index=False)
    # traps in the labels folder
    pd.DataFrame(dict(mutant=["A1C", "A2C"], DMS_score_bin=[0, 2])).to_csv(os.path.join(L, "BADLABEL.csv"), index=False)
    pd.DataFrame(dict(mutant=["A1C"], other=[1])).to_csv(os.path.join(L, "NOLABELCOL.csv"), index=False)
    return L, S


class TestClinicalScript(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        r = cls.tmp.name
        cls.L, cls.S = build_mock(r)
        prev = pd.DataFrame(dict(Model=["ModelGood", "GhostModel"], Original_Correlation=[0.4, 0.3],
                                 Absolute_Correlation=[0.4, 0.3], N_Mutations=[1000, 900]))
        cls.prev = os.path.join(r, "prev.csv"); prev.to_csv(cls.prev, index=False)
        cls.out1 = os.path.join(r, "run1")
        cls.sum1 = rca.main(["run", "--labels-dir", cls.L, "--scores-dir", cls.S, "--out", cls.out1,
                             "--n-boot", "20", "--previous-table", cls.prev])
        mp = pd.DataFrame(dict(model=["ModelAmbiguous", "ModelLower"], column=["a", "raw"],
                               orientation=["higher_pathogenic", "lower_pathogenic"]))
        cls.mp = os.path.join(r, "map.csv"); mp.to_csv(cls.mp, index=False)
        cls.out2 = os.path.join(r, "run2")
        cls.sum2 = rca.main(["run", "--labels-dir", cls.L, "--scores-dir", cls.S, "--out", cls.out2,
                             "--n-boot", "20", "--score-map", cls.mp])

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_nothing_is_skipped_silently(self):
        st = pd.read_csv(os.path.join(self.out1, "audit_models.csv")).set_index("model")
        self.assertEqual(st.loc["ModelAmbiguous", "status"], "skipped")
        self.assertIn("cannot choose the score column", st.loc["ModelAmbiguous", "message"])
        self.assertEqual(st.loc["ModelBroken", "status"], "skipped")
        self.assertIn("mutant", st.loc["ModelBroken", "message"])
        self.assertEqual(st.loc["ModelGood", "status"], "analysed")
        self.assertEqual(self.sum1["models_skipped"], 2)

    def test_skipped_models_absent_from_results(self):
        m = pd.read_csv(os.path.join(self.out1, "all_scored", "models.csv"))
        self.assertNotIn("ModelAmbiguous", m["model"].tolist()); self.assertNotIn("ModelBroken", m["model"].tolist())

    def test_bad_label_files_are_reported(self):
        la = pd.read_csv(os.path.join(self.out1, "audit_label_files.csv")).set_index("file")
        self.assertEqual(la.loc["BADLABEL.csv", "status"], "failed"); self.assertIn("not 0/1", la.loc["BADLABEL.csv", "message"])
        self.assertEqual(la.loc["NOLABELCOL.csv", "status"], "failed")
        self.assertEqual(self.sum1["label_files_failed"], 2)
        self.assertEqual(set(la.loc[la.status == "ok", "protein_source"]), {"column", "filename"})

    def test_orientation_warning_and_map(self):
        o1 = pd.read_csv(os.path.join(self.out1, "orientation_check.csv")).set_index("model")
        self.assertLess(o1.loc["ModelLower", "pooled_auroc"], 0.5); self.assertIsInstance(o1.loc["ModelLower", "warning"], str)
        self.assertIn("orientation", o1.loc["ModelLower", "warning"])
        o2 = pd.read_csv(os.path.join(self.out2, "orientation_check.csv")).set_index("model")
        self.assertGreater(o2.loc["ModelLower", "pooled_auroc"], 0.5); self.assertTrue(pd.isna(o2.loc["ModelLower", "warning"]))

    def test_map_resolves_ambiguous_model(self):
        st = pd.read_csv(os.path.join(self.out2, "audit_models.csv")).set_index("model")
        self.assertEqual(st.loc["ModelAmbiguous", "status"], "analysed"); self.assertEqual(st.loc["ModelAmbiguous", "score_column"], "a")
        self.assertEqual(self.sum2["models_skipped"], 1)   # only ModelBroken remains

    def test_common_set_differs_from_all_scored(self):
        a = pd.read_csv(os.path.join(self.out1, "all_scored", "models.csv")).set_index("model")
        c = pd.read_csv(os.path.join(self.out1, "common_variants", "models.csv")).set_index("model")
        self.assertAlmostEqual(a.loc["ModelPartial", "coverage"], 0.6, delta=0.05)
        self.assertAlmostEqual(c.loc["ModelPartial", "coverage"], 1.0)
        self.assertLess(self.sum1["common_variants"], self.sum1["variants"])
        self.assertEqual(c.loc["ModelGood", "n_scored"], self.sum1["common_variants"])

    def test_pooled_versus_perprotein_gap(self):
        agg = pd.read_csv(os.path.join(self.out1, "all_scored", "aggregation_comparison.csv")).set_index("model")
        self.assertGreater(agg.loc["ModelBetween", "pooled"] - agg.loc["ModelBetween", "pergroup"], 0.1)
        self.assertGreater(agg.loc["ModelGood", "pergroup"], 0.65)

    def test_legacy_spearman_matches_scipy_and_previous_table_reconciled(self):
        leg = pd.read_csv(os.path.join(self.out1, "legacy_pooled_spearman.csv")).set_index("model")
        parts = []
        for f in sorted(os.listdir(os.path.join(self.S, "ModelGood"))):
            s = pd.read_csv(os.path.join(self.S, "ModelGood", f)); lab = pd.read_csv(os.path.join(self.L, f))
            parts.append(pd.DataFrame(dict(y=lab["DMS_score_bin"].to_numpy(), s=s["ModelGood_score"].to_numpy())))
        d = pd.concat(parts)
        self.assertAlmostEqual(leg.loc["ModelGood", "original_correlation"], round(spearmanr(d.y, d.s).statistic, 4), places=4)
        self.assertAlmostEqual(leg.loc["ModelGood", "absolute_correlation"], abs(leg.loc["ModelGood", "original_correlation"]), places=6)
        c = pd.read_csv(os.path.join(self.out1, "comparison_with_previous_table.csv")).set_index("model_name")
        self.assertEqual(c.loc["GhostModel", "note"], "in previous table, not analysed now")
        self.assertIn("analysed now", c.loc["ModelLower", "note"])

    def test_outputs_exist(self):
        for f in ("RUN_SUMMARY.json", "audit_label_files.csv", "audit_score_files.csv", "audit_models.csv",
                  "legacy_pooled_spearman.csv", "orientation_check.csv"):
            self.assertTrue(os.path.exists(os.path.join(self.out1, f)), f)
        for f in ("models.csv", "per_protein_auroc.csv", "aggregation_comparison.csv", "gene_identity_baseline.csv", "summary.json"):
            self.assertTrue(os.path.exists(os.path.join(self.out1, "all_scored", f)), f)

    def test_inspect_builds_template(self):
        out = os.path.join(self.tmp.name, "insp")
        rca.main(["inspect", "--labels-dir", self.L, "--scores-dir", self.S, "--out", out])
        t = pd.read_csv(os.path.join(out, "score_map_template.csv")).set_index("model")
        self.assertEqual(t.loc["ModelGood", "column"], "ModelGood_score")
        self.assertTrue(pd.isna(t.loc["ModelAmbiguous", "column"]))   # left blank, never guessed


if __name__ == "__main__":
    unittest.main()
