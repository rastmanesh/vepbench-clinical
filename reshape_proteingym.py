import zipfile
import pandas as pd
from pathlib import Path

LABELS_ZIP = Path(r"D:\Bio\clinical_ProteinGym_substitutions.zip")
SCORES_DIR = Path(r"D:\Bio\zero_shot_clinical_substitutions_scores")
LABELS_OUT = Path(r"D:\Bio\labels_reshaped")
SCORES_OUT = Path(r"D:\Bio\scores_reshaped")

IDENTIFIER_COLS = ["Unnamed: 0", "protein", "protein_sequence",
                   "mutant", "mutated_sequence", "DMS_bin_score"]
LABEL_MAP = {"Pathogenic": 1, "Benign": 0}

LABELS_OUT.mkdir(parents=True, exist_ok=True)
SCORES_OUT.mkdir(parents=True, exist_ok=True)

# --- 1. Labels: from the zip, one per-protein CSV ---
z = zipfile.ZipFile(LABELS_ZIP)
n_labels = 0
for name in z.namelist():
    df = pd.read_csv(z.open(name))
    if "mutant" not in df.columns or "DMS_bin_score" not in df.columns:
        continue
    out = pd.DataFrame({
        "mutant": df["mutant"].astype(str),
        "DMS_bin_score": df["DMS_bin_score"].map(LABEL_MAP)
    }).dropna(subset=["DMS_bin_score"])
    out["DMS_bin_score"] = out["DMS_bin_score"].astype(int)
    out.to_csv(LABELS_OUT / name, index=False)
    n_labels += 1
print("label files written:", n_labels)

# --- 2. Scores: wide per-protein file -> one folder per model ---
first = next(SCORES_DIR.glob("*.csv"))
cols = list(pd.read_csv(first, nrows=0).columns)
model_cols = [c for c in cols if c not in IDENTIFIER_COLS]
print("models:", len(model_cols))

written = {m: 0 for m in model_cols}
for f in sorted(SCORES_DIR.glob("*.csv")):
    df = pd.read_csv(f)
    for m in model_cols:
        if m not in df.columns:
            continue
        sub = df[["mutant", m]].rename(columns={m: "score"}).dropna(subset=["score"])
        if len(sub) == 0:
            continue
        outdir = SCORES_OUT / m
        outdir.mkdir(exist_ok=True)
        sub.to_csv(outdir / f.name, index=False)
        written[m] += 1
print("score files per model:")
for m, n in written.items():
    print(f"  {m}: {n}")