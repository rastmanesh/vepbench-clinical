"""Ranking comparison: do conclusions about which predictor is best survive a change of design?"""
import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr
from .metrics import auroc
from .genelevel import per_group_auroc


def rank_models(values, higher_is_better=True):
    """Rank a Series of per-model values (1 = best; ties share the average rank)."""
    s = pd.Series(values)
    return s.rank(ascending=not higher_is_better, method="average")


def compare_rankings(a, b):
    """Agreement between two Series of per-model values indexed by model.

    Returns Spearman and Kendall correlations, the fraction of model pairs
    ordered differently, the largest absolute rank change, and the number of
    models whose rank changes by at least one position.
    """
    a = pd.Series(a)
    b = pd.Series(b)
    common = a.index.intersection(b.index)
    a, b = a.loc[common], b.loc[common]
    ra, rb = rank_models(a), rank_models(b)
    tau = kendalltau(a, b).statistic
    rho = spearmanr(a, b).statistic
    n = len(common)
    disc = 0
    tot = 0
    va, vb = a.to_numpy(), b.to_numpy()
    for i in range(n):
        for j in range(i + 1, n):
            da, db = va[i] - va[j], vb[i] - vb[j]
            if da == 0 or db == 0:
                continue
            tot += 1
            disc += (da > 0) != (db > 0)
    change = (ra - rb).abs()
    return dict(n_models=int(n), spearman=float(rho), kendall=float(tau),
                discordant_fraction=float(disc / tot) if tot else float("nan"),
                max_abs_rank_change=float(change.max()),
                n_rank_changed=int((change >= 1).sum()))


def bootstrap_rank_agreement(df, group, label, score_cols, n_boot=200, seed=0,
                             min_per_class=5, external=None):
    """Bootstrap (resampling groups) the agreement between pooled and per-group rankings.

    For each resample of groups, predictors are ranked by pooled AUROC and by the
    mean per-group AUROC, and Kendall's tau between the two rankings is recorded.
    If ``external`` (a Series of per-model values, for example a ranking from an
    independent assay-based benchmark) is given, the Spearman correlation of each
    ranking with it is recorded as well.
    Returns (samples DataFrame, summary dict).
    """
    codes, _ = pd.factorize(df[group])
    ng = int(codes.max()) + 1
    order = np.argsort(codes, kind="stable")
    bounds = np.r_[0, np.cumsum(np.bincount(codes))]
    idx = [order[bounds[i]:bounds[i + 1]] for i in range(ng)]
    y = df[label].to_numpy()
    scores = {c: df[c].to_numpy(float) for c in score_cols}
    pg = {}
    for c in score_cols:
        t = per_group_auroc(df.assign(_c=codes), "_c", label, c, min_per_class)
        arr = np.full(ng, np.nan)
        arr[t["group"].to_numpy()] = t["auroc"].to_numpy()
        pg[c] = arr
    rng = np.random.default_rng(seed)
    rows = []
    ext = None if external is None else pd.Series(external).reindex(score_cols)
    for _ in range(n_boot):
        pick = rng.integers(0, ng, size=ng)
        rws = np.concatenate([idx[i] for i in pick])
        pooled = pd.Series({c: auroc(y[rws], scores[c][rws]) for c in score_cols})
        pergrp = pd.Series({c: np.nanmean(pg[c][pick]) if np.isfinite(pg[c][pick]).any() else np.nan
                            for c in score_cols})
        rec = {"kendall_pooled_vs_pergroup": kendalltau(pooled, pergrp).statistic}
        if ext is not None:
            rec["spearman_pooled_vs_external"] = spearmanr(pooled, ext).statistic
            rec["spearman_pergroup_vs_external"] = spearmanr(pergrp, ext).statistic
            rec["external_gain_pergroup_minus_pooled"] = (rec["spearman_pergroup_vs_external"]
                                                          - rec["spearman_pooled_vs_external"])
        rows.append(rec)
    samples = pd.DataFrame(rows)
    from .bootstrap import percentile_ci
    summary = {}
    for c in samples.columns:
        lo, hi = percentile_ci(samples[c].to_numpy())
        summary[c] = dict(mean=float(np.nanmean(samples[c])), ci_low=lo, ci_high=hi)
    return samples, summary
