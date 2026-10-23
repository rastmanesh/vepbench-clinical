"""Gene-aware evaluation: per-group metrics, aggregation contrasts, gene-identity baselines."""
import numpy as np
import pandas as pd
from .metrics import auroc


def per_group_auroc(df, group, label, score, min_per_class=5):
    """AUROC within each group (gene or protein).

    Groups with fewer than ``min_per_class`` variants of either class get NaN.
    Returns a DataFrame with columns group, n_pos, n_neg, auroc.
    """
    rows = []
    for g, sub in df.groupby(group, sort=False):
        s = sub[score].to_numpy(float)
        y = sub[label].to_numpy()
        ok = ~np.isnan(s)
        y, s = y[ok], s[ok]
        n1 = int(y.sum())
        n0 = int(y.size - n1)
        a = auroc(y, s) if (n1 >= min_per_class and n0 >= min_per_class) else float("nan")
        rows.append((g, n1, n0, a))
    return pd.DataFrame(rows, columns=["group", "n_pos", "n_neg", "auroc"])


def summarize_per_group(pg):
    """Summaries of a per_group_auroc table (eligible groups only)."""
    ok = pg["auroc"].notna()
    e = pg[ok]
    if e.empty:
        return dict(n_groups=int(len(pg)), n_eligible=0, mean=float("nan"),
                    weighted_mean=float("nan"), median=float("nan"))
    w = (e["n_pos"] * e["n_neg"] / (e["n_pos"] + e["n_neg"])).to_numpy(float)
    return dict(n_groups=int(len(pg)), n_eligible=int(ok.sum()),
                mean=float(e["auroc"].mean()),
                weighted_mean=float(np.average(e["auroc"], weights=w)),
                median=float(e["auroc"].median()))


def compare_aggregations(df, group, label, score_cols, min_per_class=5):
    """Pooled versus per-group AUROC for each predictor, with ranks.

    pooled: one AUROC over all scored variants.
    pooled_eligible: the same pooled AUROC restricted to the groups that qualify for the
        per-group summary (at least ``min_per_class`` variants of each class).
    pergroup: mean of per-group AUROCs (eligible groups only).
    The per-group summary cannot use groups that are nearly all one class, so comparing it with the
    all-variant pooled value mixes two effects: a different aggregation AND a different set of
    variants. ``gap_selection`` (pooled minus pooled_eligible) is the part due to the variant set;
    ``gap_aggregation`` (pooled_eligible minus pergroup) is the part due to pooling across groups.
    """
    rows = []
    for c in score_cols:
        sub = df[[group, label, c]].dropna(subset=[c])
        pooled = auroc(sub[label].to_numpy(), sub[c].to_numpy())
        pg = per_group_auroc(df, group, label, c, min_per_class)
        sm = summarize_per_group(pg)
        elig = pg.loc[pg["auroc"].notna(), "group"]
        se = sub[sub[group].isin(set(elig))]
        pooled_e = auroc(se[label].to_numpy(), se[c].to_numpy()) if len(elig) else float("nan")
        rows.append(dict(model=c, pooled=pooled, pooled_eligible=pooled_e, pergroup=sm["mean"],
                         pergroup_weighted=sm["weighted_mean"],
                         n_eligible_groups=sm["n_eligible"],
                         coverage=float(df[c].notna().mean())))
    out = pd.DataFrame(rows)
    out["gap_pooled_minus_pergroup"] = out["pooled"] - out["pergroup"]
    out["gap_selection"] = out["pooled"] - out["pooled_eligible"]
    out["gap_aggregation"] = out["pooled_eligible"] - out["pergroup"]
    out["rank_pooled"] = out["pooled"].rank(ascending=False, method="average")
    out["rank_pooled_eligible"] = out["pooled_eligible"].rank(ascending=False, method="average")
    out["rank_pergroup"] = out["pergroup"].rank(ascending=False, method="average")
    out["rank_change"] = out["rank_pooled"] - out["rank_pergroup"]
    return out


def _smoothed_rate(k, n, prior, a):
    return (k + a * prior) / (n + a)


def gene_identity_baseline(df, group, label, scheme="kfold", n_splits=5, seed=0, alpha=1.0):
    """Out-of-sample 'predictor' that knows only the gene.

    Each variant is scored by the smoothed fraction of pathogenic variants in its
    gene, estimated without that variant's own fold. Schemes:
      in_sample    gene rate from all variants (leaky; an upper bound).
      kfold        random variant-level K-fold (the usual random split).
      loo          leave-one-variant-out: the variant's own label is excluded.
    Caution: with kfold and loo, excluding the variant's own label (or fold) from the
    gene rate makes the score negatively related to the label *within* a gene, so the
    within-gene AUROC of these baselines is biased below 0.5 and is not a chance
    anchor. Use the pooled AUROC of these schemes, and the per-group AUROC only for
    in_sample and group_kfold (where scores are constant within a gene, giving 0.5).
      group_kfold  K-fold over genes; unseen genes get the training-set rate,
                   so the score carries no gene information (a sanity anchor).
    Returns a float array of scores aligned to ``df``.
    """
    y = df[label].to_numpy().astype(int)
    g, _ = pd.factorize(df[group])
    ng = int(g.max()) + 1
    n = y.size
    rng = np.random.default_rng(seed)
    if scheme == "in_sample":
        k = np.bincount(g, weights=y, minlength=ng)
        m = np.bincount(g, minlength=ng)
        return _smoothed_rate(k[g], m[g], y.mean(), alpha)
    if scheme == "loo":
        k = np.bincount(g, weights=y, minlength=ng)
        m = np.bincount(g, minlength=ng)
        prior = (y.sum() - y) / (n - 1)
        return _smoothed_rate(k[g] - y, m[g] - 1, prior, alpha)
    if scheme not in ("kfold", "group_kfold"):
        raise ValueError("unknown scheme: %s" % scheme)
    if scheme == "kfold":
        fold = rng.integers(0, n_splits, size=n)
    else:
        gfold = rng.integers(0, n_splits, size=ng)
        fold = gfold[g]
    scores = np.empty(n, dtype=float)
    for f in range(n_splits):
        test = fold == f
        train = ~test
        if not test.any():
            continue
        prior = y[train].mean()
        k = np.bincount(g[train], weights=y[train], minlength=ng)
        m = np.bincount(g[train], minlength=ng)
        scores[test] = _smoothed_rate(k[g[test]], m[g[test]], prior, alpha)
    return scores


def gene_identity_report(df, group, label, schemes=("in_sample", "kfold", "loo", "group_kfold"),
                         n_splits=5, seed=0, alpha=1.0, min_per_class=5):
    """Pooled AUROC of the gene-identity baseline under each scheme.

    The per-group AUROC is reported only where it is meaningful (in_sample and
    group_kfold, where scores are constant within a gene); for kfold and loo it is
    NaN because out-of-fold gene rates are biased within genes (see
    ``gene_identity_baseline``). ``pergroup_biased`` gives the biased value for
    transparency.
    """
    rows = []
    tmp = df[[group, label]].copy()
    for sch in schemes:
        tmp["_s"] = gene_identity_baseline(df, group, label, sch, n_splits, seed, alpha)
        pg = summarize_per_group(per_group_auroc(tmp, group, label, "_s", min_per_class))["mean"]
        valid = sch in ("in_sample", "group_kfold")
        rows.append(dict(scheme=sch, pooled=auroc(tmp[label].to_numpy(), tmp["_s"].to_numpy()),
                         pergroup=pg if valid else float("nan"), pergroup_biased=pg))
    return pd.DataFrame(rows)
