#!/usr/bin/env python3
"""Regression check: do new within-protein results agree with an earlier run's per-protein table?

An earlier run's per_protein_auroc.csv (columns group, n_pos, n_neg, auroc, model) gives, for each model, a pair-weighted
mean of per-protein AUROCs over the proteins that were eligible then. The new pairwise within-protein AUROC uses every protein
with both classes, so the two should agree closely (they differ only by the small share of pairs in ineligible proteins).
A model whose values differ by more than --tol is flagged: investigate before trusting either run.
Also compares pooled AUROC.

  python analysis/reconcile_with_first_run.py --first-run-per-protein old/all_scored/per_protein_auroc.csv \
      --first-run-models old/all_scored/models.csv --new-models new/all_scored/models.csv --tol 0.01
"""
import argparse
import pandas as pd


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--first-run-per-protein", required=True)
    ap.add_argument("--first-run-models", required=True)
    ap.add_argument("--new-models", required=True)
    ap.add_argument("--tol", type=float, default=0.01)
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    pp = pd.read_csv(a.first_run_per_protein)
    old = pd.read_csv(a.first_run_models).set_index("model")
    new = pd.read_csv(a.new_models).set_index("model")
    rows = []
    for k, t in pp.groupby("model"):
        pairs = (t["n_pos"] * t["n_neg"]).astype(float)
        el = t["auroc"].notna()
        old_within = float((t["auroc"][el] * pairs[el]).sum() / pairs[el].sum())
        rows.append(dict(model=k, old_within_pair_weighted=old_within))
    r = pd.DataFrame(rows).set_index("model")
    r["new_within_pairwise"] = new["within_pairwise_auroc"]
    r["within_difference"] = r["new_within_pairwise"] - r["old_within_pair_weighted"]
    r["old_pooled"] = old["pooled_auroc"]
    r["new_pooled"] = new["pooled_auroc"]
    r["pooled_difference"] = r["new_pooled"] - r["old_pooled"]
    r["flag"] = ((r["within_difference"].abs() > a.tol) | (r["pooled_difference"].abs() > 1e-6)).map({True: "CHECK", False: ""})
    missing = sorted(set(old.index) ^ set(new.index))
    print(r.round(4).to_string())
    print("\nmodels flagged: %d of %d; models in only one run: %s" % ((r["flag"] == "CHECK").sum(), len(r), missing or "none"))
    if a.out:
        r.to_csv(a.out)
    return r


if __name__ == "__main__":
    main()
