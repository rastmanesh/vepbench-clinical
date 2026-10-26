"""Command-line interface: ``python -m vepbench evaluate ...``"""
import argparse
import json
import os
import sys
import numpy as np
import pandas as pd

from .io import load_table
from .metrics import auroc, average_precision, sensitivity_at_specificity, partial_auroc
from .genelevel import per_group_auroc, summarize_per_group, compare_aggregations, gene_identity_report
from .bootstrap import (bootstrap_decomposition, bootstrap_pergroup_mean, bootstrap_pooled_auroc, bootstrap_within_group_auroc,
                        percentile_ci)
from .genelevel import within_group_auroc, between_minus_within
from .ranking import compare_rankings, bootstrap_rank_agreement
from . import __version__


def evaluate(df, group, label, scores, min_per_class=5, n_boot=200, seed=0,
             specificity=0.95, max_fpr=0.05, external=None, n_splits=5):
    """Run the full evaluation; returns a dict of tables and summaries."""
    rows, pg_all = [], []
    for c in scores:
        sub = df[[group, label, c]].dropna(subset=[c])
        pooled = auroc(sub[label].to_numpy(), sub[c].to_numpy())
        pg = per_group_auroc(df, group, label, c, min_per_class)
        sm = summarize_per_group(pg)
        pg = pg.assign(model=c)
        pg_all.append(pg)
        dec = between_minus_within(df, group, label, c)
        elig = set(pg.loc[pg["auroc"].notna(), "group"])
        se = sub[sub[group].isin(elig)]
        row = dict(model=c, n_scored=int(len(sub)), coverage=float(len(sub) / len(df)),
                   pooled_auroc=pooled,
                   within_pairwise_auroc=dec["within_auroc"], between_pairwise_auroc=dec["between_auroc"],
                   gap_between_minus_within=dec["gap_between_minus_within"], within_pair_share=dec["within_share"],
                   pooled_auroc_eligible_genes=(auroc(se[label].to_numpy(), se[c].to_numpy()) if elig else float("nan")), pergroup_mean_auroc=sm["mean"],
                   pergroup_weighted_auroc=sm["weighted_mean"], pergroup_median_auroc=sm["median"],
                   n_eligible_groups=sm["n_eligible"],
                   average_precision=average_precision(sub[label].to_numpy(), sub[c].to_numpy()),
                   sensitivity_at_specificity=sensitivity_at_specificity(
                       sub[label].to_numpy(), sub[c].to_numpy(), specificity),
                   partial_auroc=partial_auroc(sub[label].to_numpy(), sub[c].to_numpy(), max_fpr))
        if n_boot > 0:
            bp = bootstrap_pooled_auroc(df, group, label, c, n_boot, seed)
            row["pooled_ci_low"], row["pooled_ci_high"] = percentile_ci(bp)
            bg = bootstrap_pergroup_mean(pg["auroc"].to_numpy(), n_boot, seed)
            row["pergroup_ci_low"], row["pergroup_ci_high"] = percentile_ci(bg)
            bd = bootstrap_decomposition(df, group, label, c, n_boot, seed)
            row["within_pairwise_ci_low"], row["within_pairwise_ci_high"] = percentile_ci(bd["within"].to_numpy())
            row["between_pairwise_ci_low"], row["between_pairwise_ci_high"] = percentile_ci(bd["between"].to_numpy())
            row["gap_between_minus_within_ci_low"], row["gap_between_minus_within_ci_high"] = percentile_ci(bd["gap_between_minus_within"].to_numpy())
        rows.append(row)
    models = pd.DataFrame(rows)
    agg = compare_aggregations(df, group, label, scores, min_per_class)
    baseline = gene_identity_report(df, group, label, n_splits=n_splits, seed=seed,
                                    min_per_class=min_per_class)
    ranking = compare_rankings(agg.set_index("model")["pooled"], agg.set_index("model")["pergroup"])
    ai = agg.set_index("model")
    result = dict(models=models, per_group=pd.concat(pg_all, ignore_index=True),
                  aggregation=agg, gene_identity=baseline, ranking_pooled_vs_pergroup=ranking,
                  ranking_pooled_eligible_vs_pergroup=compare_rankings(ai["pooled_eligible"], ai["pergroup"]),
                  ranking_pooled_vs_within_pairwise=compare_rankings(ai["pooled"], ai["within_pairwise"]))
    if n_boot > 0 and len(scores) > 1:
        samples, summ = bootstrap_rank_agreement(df, group, label, scores, n_boot=n_boot, seed=seed,
                                                 min_per_class=min_per_class, external=external)
        result["ranking_bootstrap"] = summ
    if external is not None:
        ext = pd.Series(external)
        result["ranking_vs_external"] = dict(
            pooled=compare_rankings(agg.set_index("model")["pooled"], ext),
            pergroup=compare_rankings(agg.set_index("model")["pergroup"], ext))
    return result


def _table(df):
    """Markdown table if the optional 'tabulate' package is installed, plain text otherwise."""
    try:
        return df.to_markdown(index=False)
    except ImportError:
        return df.to_string(index=False)


def _report(result, args, n_rows, n_groups):
    m = result["models"]
    lines = ["# vepbench evaluation report", "",
             "vepbench version %s; seed %d; bootstrap resamples %d (genes resampled, not variants)."
             % (__version__, args.seed, args.n_boot),
             "Variants: %d in %d groups. Labels: 1 = pathogenic, 0 = benign. Minimum variants per class per group: %d." %
             (n_rows, n_groups, args.min_per_class), "",
             "## Predictors", "", _table(m.round(3)),
             "", "## Gene-identity baseline (a predictor that knows only the gene)", "",
             result["gene_identity"].round(3).to_string(index=False),
             "", "## Pooled versus per-group ranking", "", json.dumps(result["ranking_pooled_vs_pergroup"], indent=2)]
    if "ranking_bootstrap" in result:
        lines += ["", "## Bootstrap agreement", "", json.dumps(result["ranking_bootstrap"], indent=2)]
    if "ranking_vs_external" in result:
        lines += ["", "## Agreement with the external ranking", "", json.dumps(result["ranking_vs_external"], indent=2)]
    return "\n".join(lines) + "\n"


def main(argv=None):
    p = argparse.ArgumentParser(prog="vepbench", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("evaluate", help="evaluate predictors on a labelled, grouped table")
    e.add_argument("--input", required=True, help="CSV/TSV/Parquet file, one row per variant")
    e.add_argument("--group", required=True, help="column naming the gene or protein")
    e.add_argument("--label", required=True, help="column with 1 = pathogenic, 0 = benign")
    e.add_argument("--scores", required=True, help="comma-separated predictor columns (higher = more pathogenic)")
    e.add_argument("--positive-label", default=None, help="value in the label column meaning pathogenic")
    e.add_argument("--flip", default="", help="comma-separated predictors where LOWER means pathogenic")
    e.add_argument("--min-per-class", type=int, default=5, help="minimum variants of each class for a group to count (default 5)")
    e.add_argument("--n-boot", type=int, default=200, help="bootstrap resamples; 0 skips intervals (default 200)")
    e.add_argument("--seed", type=int, default=0, help="random seed (default 0)")
    e.add_argument("--specificity", type=float, default=0.95, help="specificity for sensitivity-at-specificity (default 0.95)")
    e.add_argument("--max-fpr", type=float, default=0.05, help="false-positive-rate limit for partial AUROC (default 0.05)")
    e.add_argument("--n-splits", type=int, default=5, help="folds for the gene-identity baseline (default 5)")
    e.add_argument("--external", default=None, help="CSV with columns model,value: an independent ranking (for example from deep mutational scanning)")
    e.add_argument("--out", required=True, help="output directory")
    a = p.parse_args(argv)
    scores = [s for s in a.scores.split(",") if s]
    flip = [s for s in a.flip.split(",") if s]
    df = load_table(a.input, a.group, a.label, scores, a.positive_label, flip)
    ext = None
    if a.external:
        t = pd.read_csv(a.external)
        ext = pd.Series(t["value"].to_numpy(), index=t["model"].to_numpy())
    res = evaluate(df, a.group, a.label, scores, a.min_per_class, a.n_boot, a.seed, a.specificity,
                   a.max_fpr, ext, a.n_splits)
    os.makedirs(a.out, exist_ok=True)
    res["models"].to_csv(os.path.join(a.out, "models.csv"), index=False)
    res["per_group"].to_csv(os.path.join(a.out, "per_group_auroc.csv"), index=False)
    res["aggregation"].to_csv(os.path.join(a.out, "aggregation_comparison.csv"), index=False)
    res["gene_identity"].to_csv(os.path.join(a.out, "gene_identity_baseline.csv"), index=False)
    summary = {k: v for k, v in res.items() if isinstance(v, dict)}
    with open(os.path.join(a.out, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    with open(os.path.join(a.out, "report.md"), "w") as fh:
        fh.write(_report(res, a, len(df), df[a.group].nunique()))
    print("wrote results to", a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
