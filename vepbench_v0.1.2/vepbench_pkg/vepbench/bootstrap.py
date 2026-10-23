"""Cluster-aware resampling: genes (or proteins) are resampled, not variants."""
import numpy as np
import pandas as pd


def percentile_ci(samples, alpha=0.05):
    """Percentile interval of bootstrap samples (NaNs ignored)."""
    x = np.asarray(samples, dtype=float)
    x = x[~np.isnan(x)]
    if x.size == 0:
        return (float("nan"), float("nan"))
    lo, hi = np.quantile(x, [alpha / 2.0, 1.0 - alpha / 2.0])
    return (float(lo), float(hi))


def _group_indices(df, group):
    codes, _ = pd.factorize(df[group])
    order = np.argsort(codes, kind="stable")
    counts = np.bincount(codes)
    bounds = np.r_[0, np.cumsum(counts)]
    return [order[bounds[i]:bounds[i + 1]] for i in range(len(counts))]


def bootstrap_statistic(df, group, stat_fn, n_boot=1000, seed=0):
    """Resample groups with replacement and recompute ``stat_fn(sub_df)``.

    Use this for statistics that need all variants (for example the pooled AUROC).
    Returns an array of ``n_boot`` bootstrap values.
    """
    idx = _group_indices(df, group)
    ng = len(idx)
    rng = np.random.default_rng(seed)
    out = np.empty(n_boot)
    for b in range(n_boot):
        pick = rng.integers(0, ng, size=ng)
        rows = np.concatenate([idx[i] for i in pick])
        out[b] = stat_fn(df.iloc[rows])
    return out


def bootstrap_pergroup_mean(values, n_boot=1000, seed=0, weights=None):
    """Bootstrap the (weighted) mean of per-group values; NaN groups are dropped."""
    v = np.asarray(values, dtype=float)
    w = np.ones_like(v) if weights is None else np.asarray(weights, dtype=float)
    ok = ~np.isnan(v)
    v, w = v[ok], w[ok]
    if v.size == 0:
        return np.full(n_boot, np.nan)
    rng = np.random.default_rng(seed)
    pick = rng.integers(0, v.size, size=(n_boot, v.size))
    return (v[pick] * w[pick]).sum(axis=1) / w[pick].sum(axis=1)


def paired_bootstrap_pergroup(pg_a, pg_b, n_boot=2000, seed=0):
    """Paired comparison of two predictors on the same groups.

    ``pg_a`` and ``pg_b`` are per_group_auroc tables. Only groups eligible in both
    are used. Returns the mean difference (a minus b), a percentile interval,
    and a two-sided bootstrap p-value.
    """
    m = pg_a[["group", "auroc"]].merge(pg_b[["group", "auroc"]], on="group", suffixes=("_a", "_b"))
    m = m.dropna()
    if m.empty:
        return dict(n_groups=0, mean_diff=float("nan"), ci_low=float("nan"),
                    ci_high=float("nan"), p_value=float("nan"))
    d = (m["auroc_a"] - m["auroc_b"]).to_numpy()
    boots = bootstrap_pergroup_mean(d, n_boot=n_boot, seed=seed)
    lo, hi = percentile_ci(boots)
    p = 2.0 * min(np.mean(boots <= 0), np.mean(boots >= 0))
    return dict(n_groups=int(d.size), mean_diff=float(d.mean()), ci_low=lo, ci_high=hi,
                p_value=float(min(1.0, p)))


def bootstrap_pooled_auroc(df, group, label, score, n_boot=200, seed=0):
    """Cluster bootstrap of the pooled AUROC (fast array version).

    The same seed gives the same sequence of resampled groups for every predictor,
    so intervals and differences for different predictors are paired.
    """
    from .metrics import auroc
    codes, _ = pd.factorize(df[group])
    order = np.argsort(codes, kind="stable")
    bounds = np.r_[0, np.cumsum(np.bincount(codes))]
    ng = bounds.size - 1
    idx = [order[bounds[i]:bounds[i + 1]] for i in range(ng)]
    y = df[label].to_numpy()
    s = df[score].to_numpy(float)
    rng = np.random.default_rng(seed)
    out = np.empty(n_boot)
    for b in range(n_boot):
        pick = rng.integers(0, ng, size=ng)
        rows = np.concatenate([idx[i] for i in pick])
        out[b] = auroc(y[rows], s[rows])
    return out
