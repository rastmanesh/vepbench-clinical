"""Pair-level diagnostics for gene-aware evaluation.

Three quantities are computed per model:
  * within_group_auroc        -- share of correctly ordered pairs that share a group
  * between_group_auroc       -- share of correctly ordered pairs that do not
  * pair_composition          -- total pairs, within, between, their shares

Plus the gene-reliance gap diagnostic:
  * protein_mean_auroc        -- AUROC of the per-group mean score
  * protein_residual_auroc    -- AUROC of the within-group residual
  * gap_mean_minus_residual   -- how much gene-level tracking dominates

The gap is the diagnostic suggested by hostile review: a purely variant-level
predictor has mean and residual AUROC close to each other and close to the
pooled value; a purely gene-identity predictor has high mean AUROC and residual
AUROC near 0.5.
"""
import numpy as np
import pandas as pd
from scipy.stats import rankdata

from .metrics import auroc
from .genelevel import _group_u_and_pairs


def between_group_auroc(df, group, label, score, min_per_class=1):
    """Pairwise between-group AUROC.

    Share of correctly ordered (pathogenic, benign) pairs, counting only pairs
    where the two variants are in DIFFERENT groups. Ties count one half.
    Returns NaN if there are fewer than two groups with both classes.
    """
    d = df[[group, label, score]].dropna(subset=[score]).copy()
    d["_g"] = pd.factorize(d[group])[0]
    y = d[label].to_numpy().astype(int)
    s = d[score].to_numpy(float)
    g = d["_g"].to_numpy()

    ng = int(g.max()) + 1
    n_pos = np.bincount(g, weights=(y == 1), minlength=ng).astype(int)
    n_neg = np.bincount(g, weights=(y == 0), minlength=ng).astype(int)

    eligible = (n_pos > 0) & (n_neg > 0)
    if eligible.sum() < 2:
        return float("nan")

    u_within, p_within = _group_u_and_pairs(d, group, label, score, min_per_class=1)
    within_correct = float(u_within.sum())
    within_total = float(p_within.sum())

    n_pos_t = int(n_pos.sum())
    n_neg_t = int(n_neg.sum())
    total_pairs = float(n_pos_t * n_neg_t)
    between_total = total_pairs - within_total
    if between_total <= 0:
        return float("nan")

    r = rankdata(s)
    pooled_correct = float(r[y == 1].sum() - n_pos_t * (n_pos_t + 1) / 2.0)
    between_correct = pooled_correct - within_correct
    return float(between_correct / between_total)


def pair_composition(df, group, label, score):
    """Return a dict with total_pairs, within_pairs, between_pairs and shares."""
    d = df[[group, label, score]].dropna(subset=[score])
    n_pos_t = int((d[label] == 1).sum())
    n_neg_t = int((d[label] == 0).sum())
    total = float(n_pos_t * n_neg_t)
    within = 0.0
    for _, sub in d.groupby(group):
        within += float(int((sub[label] == 1).sum()) * int((sub[label] == 0).sum()))
    between = total - within
    return dict(
        total_pairs=total,
        within_pairs=within,
        between_pairs=between,
        within_share=float(within / total) if total > 0 else float("nan"),
        between_share=float(between / total) if total > 0 else float("nan"),
    )


def protein_mean_vs_residual_auroc(df, group, label, score):
    """Gene-reliance gap diagnostic.

    protein_mean_auroc:      AUROC of the per-group mean score.
    protein_residual_auroc:  AUROC of the within-group residual (score minus group mean).
    gap:                     protein_mean_auroc - protein_residual_auroc.

    A purely variant-level predictor has mean and residual AUROC close; a
    purely gene-identity predictor has high mean and residual near 0.5.
    """
    d = df[[group, label, score]].dropna(subset=[score]).copy()
    gm = d.groupby(group)[score].transform("mean")
    d["_mean"] = gm
    d["_resid"] = d[score] - gm
    y = d[label].to_numpy().astype(int)

    mean_au = auroc(y, d["_mean"].to_numpy())
    resid_au = auroc(y, d["_resid"].to_numpy())
    return dict(
        protein_mean_auroc=float(mean_au),
        protein_residual_auroc=float(resid_au),
        gap_mean_minus_residual=float(mean_au - resid_au) if not (np.isnan(mean_au) or np.isnan(resid_au)) else float("nan"),
    )