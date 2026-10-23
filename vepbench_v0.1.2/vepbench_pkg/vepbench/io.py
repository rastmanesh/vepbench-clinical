"""Reading benchmark tables.

IMPORTANT: ``merge_group_files`` assumes one labels file and one scores file per
protein, sharing a key column. The layout of the real ProteinGym release was NOT
tested here; check file names and column names against your download and adjust
the arguments.
"""
import os
import numpy as np
import pandas as pd


def load_table(path, group, label, scores, positive_label=None, flip=()):
    """Read a CSV/TSV (or Parquet if pyarrow is installed) benchmark table.

    group, label: column names. scores: list of predictor columns.
    positive_label: value in the label column that means pathogenic (if the labels
    are not already 0/1). flip: predictors where a LOWER score means pathogenic;
    their scores are negated.
    """
    ext = os.path.splitext(path)[1].lower()
    if ext == ".parquet":
        df = pd.read_parquet(path)
    else:
        df = pd.read_csv(path, sep="\t" if ext in (".tsv", ".txt") else ",")
    missing = [c for c in [group, label] + list(scores) if c not in df.columns]
    if missing:
        raise KeyError("columns not found: %s" % missing)
    if positive_label is not None:
        df[label] = (df[label] == positive_label).astype(int)
    if not np.isin(df[label].dropna().unique(), (0, 1)).all():
        raise ValueError("label column must be 0/1 (or give positive_label)")
    df = df.dropna(subset=[label]).copy()
    df[label] = df[label].astype(int)
    for c in flip:
        df[c] = -df[c]
    return df[[group, label] + list(scores)]


def merge_group_files(label_dir, score_dir, label_col="DMS_score_bin", key="mutant", ext=".csv"):
    """Merge per-protein label files with per-protein score files.

    Files are matched by name (without extension). Returns one DataFrame with a
    ``group`` column. Proteins without a matching scores file are skipped, and the
    skipped names are returned as the second value.
    """
    frames, skipped = [], []
    for fn in sorted(os.listdir(label_dir)):
        if not fn.endswith(ext):
            continue
        stem = fn[: -len(ext)]
        sp = os.path.join(score_dir, fn)
        if not os.path.exists(sp):
            skipped.append(stem)
            continue
        lab = pd.read_csv(os.path.join(label_dir, fn), usecols=[key, label_col])
        sc = pd.read_csv(sp)
        m = lab.merge(sc, on=key, how="inner")
        m.insert(0, "group", stem)
        frames.append(m.rename(columns={label_col: "label"}))
    if not frames:
        raise ValueError("no matching files found")
    return pd.concat(frames, ignore_index=True), skipped
