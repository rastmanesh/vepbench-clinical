"""Threshold-free and threshold-based metrics for binary labels.

Convention: labels are 0 (benign) and 1 (pathogenic); a higher score means
more likely pathogenic. Pairs with a missing (NaN) score are dropped.
"""
import numpy as np
from scipy.stats import rankdata


def _clean(y, s):
    y = np.asarray(y)
    s = np.asarray(s, dtype=float)
    if y.shape != s.shape or y.ndim != 1:
        raise ValueError("labels and scores must be one-dimensional and the same length")
    keep = ~np.isnan(s)
    y, s = y[keep], s[keep]
    if y.size and not np.isin(y, (0, 1)).all():
        raise ValueError("labels must be 0 or 1")
    return y.astype(int), s


def auroc(y, s):
    """Area under the ROC curve (Mann-Whitney U); ties count one half.

    Returns NaN if either class is absent after dropping missing scores.
    """
    y, s = _clean(y, s)
    n1 = int(y.sum())
    n0 = y.size - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    ranks = rankdata(s)
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def roc_points(y, s):
    """(fpr, tpr) at each distinct threshold, starting at (0, 0)."""
    y, s = _clean(y, s)
    n1 = int(y.sum())
    n0 = y.size - n1
    if n1 == 0 or n0 == 0:
        raise ValueError("both classes are required")
    order = np.argsort(-s, kind="stable")
    ys, ss = y[order], s[order]
    tps = np.cumsum(ys)
    fps = np.cumsum(1 - ys)
    last = np.r_[np.flatnonzero(np.diff(ss)), ys.size - 1]
    tpr = np.r_[0.0, tps[last] / n1]
    fpr = np.r_[0.0, fps[last] / n0]
    return fpr, tpr


def partial_auroc(y, s, max_fpr=0.05, standardize=True):
    """Area under the ROC curve for FPR <= max_fpr.

    With ``standardize=True`` the McClish correction is applied, so 0.5 is
    chance and 1.0 is perfect (this matches scikit-learn's ``max_fpr``).
    Otherwise the raw area is returned.
    """
    if not 0 < max_fpr <= 1:
        raise ValueError("max_fpr must be in (0, 1]")
    try:
        fpr, tpr = roc_points(y, s)
    except ValueError:                      # one class missing: undefined, not an error
        return float("nan")
    stop = np.searchsorted(fpr, max_fpr, side="right")
    if stop < fpr.size:
        t_at = np.interp(max_fpr, fpr[stop - 1:stop + 1], tpr[stop - 1:stop + 1])
        fpr = np.r_[fpr[:stop], max_fpr]
        tpr = np.r_[tpr[:stop], t_at]
    integ = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    area = float(integ(tpr, fpr))
    if not standardize:
        return area
    lo = 0.5 * max_fpr ** 2
    hi = max_fpr
    return float(0.5 * (1.0 + (area - lo) / (hi - lo)))


def average_precision(y, s):
    """Average precision (area under the precision-recall step curve)."""
    y, s = _clean(y, s)
    n1 = int(y.sum())
    if n1 == 0:
        return float("nan")
    order = np.argsort(-s, kind="stable")
    ys, ss = y[order], s[order]
    tps = np.cumsum(ys)
    last = np.r_[np.flatnonzero(np.diff(ss)), ys.size - 1]
    tp = tps[last]
    precision = tp / (last + 1.0)
    recall = tp / n1
    recall_prev = np.r_[0.0, recall[:-1]]
    return float(np.sum((recall - recall_prev) * precision))


def sensitivity_at_specificity(y, s, specificity=0.95):
    """Sensitivity at the most permissive threshold with specificity >= target."""
    y, s = _clean(y, s)
    pos, neg = s[y == 1], s[y == 0]
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    thr = np.quantile(neg, specificity, method="higher")
    return float(np.mean(pos > thr))


def specificity_at_sensitivity(y, s, sensitivity=0.95):
    """Specificity at the most permissive threshold with sensitivity >= target."""
    y, s = _clean(y, s)
    pos, neg = s[y == 1], s[y == 0]
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    thr = np.quantile(pos, 1.0 - sensitivity, method="lower")
    return float(np.mean(neg < thr))
