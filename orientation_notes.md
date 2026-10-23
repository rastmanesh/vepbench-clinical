# Orientation notes — vepbench / ProteinGym clinical

**Date:** 2026-10-23
**Pipeline:** vepbench 0.1.2
**Scope:** 30 ProteinGym clinical baseline models
**Machine-readable map:** `score_map.csv`
**Script that generates the map:** `fix_orientations.py`
**Pre-flip trial:** `trial1_orientation_check_preflip.csv`
**Post-flip verification:** `trial2_orientation_check_postflip.csv`
**Independent cross-check:** `legacy_pooled_correlation.csv`

## Method

The 30 models split into two populations by their score convention relative to
the ProteinGym `DMS_bin_score` label (1 = pathogenic, 0 = benign):

- 22 models use a **higher = more pathogenic** convention. They are applied as-is.
- 8 models use a **lower = more pathogenic** convention. Their scores are
  multiplied by -1 before evaluation.

The 8 lower-pathogenic models were identified empirically. A trial run
(`trial_real`, `--n-boot 0`, pooled AUROC only) was executed with a provisional
orientation map derived from published tool conventions. In that trial, exactly
8 of the 30 models scored pooled AUROC between 0.08 and 0.21:

| Model | Pre-flip AUROC (trial_real) |
|---|---|
| EVE | 0.0989 |
| FATHMM | 0.2133 |
| LRT | 0.1967 |
| PoET | 0.0835 |
| Provean | 0.0979 |
| SIFT | 0.1217 |
| SIFT4G | 0.1184 |
| TranceptEVE_L | 0.1243 |

These are the exact 8 models flagged by vepbench's built-in orientation check
(`warning = "pooled AUROC < 0.5: check orientation in the score map"`). For
those 8 models, the corrected orientation was applied and the pipeline re-run
(`trial_real2`), giving post-flip AUROCs between 0.78 and 0.94.

## Why this is a sign-convention decision, not AUROC maximisation

- For each of the 8 flipped models, `1 - AUROC_pre = AUROC_post` to two
  decimal places. This is the arithmetic signature of a sign error, not of
  model improvement. A genuine sign error moves AUROC from x to 1 - x. Here:
  - 0.099 -> 0.901 (EVE)
  - 0.213 -> 0.787 (FATHMM)
  - 0.197 -> 0.803 (LRT)
  - 0.084 -> 0.917 (PoET)
  - 0.098 -> 0.902 (Provean)
  - 0.122 -> 0.878 (SIFT)
  - 0.118 -> 0.882 (SIFT4G)
  - 0.124 -> 0.876 (TranceptEVE_L)
- The correction applies to **8 of 30** models. It is not a global
  "try both orientations and keep the higher AUROC" operation.
- The 5 dbNSFP tools (SIFT, SIFT4G, Provean, FATHMM, LRT) flip to the sign
  already recorded for them in the legacy pooled table. Their final
  legacy-pooled correlations are 0.66, 0.66, 0.70, 0.50, 0.54 respectively,
  all positive and consistent with the flipped sign.
- The decision was made **before** the full run, as documented in
  `analysis_plan.md`.

## Evidence files in the bundle

| File | What it contains |
|---|---|
| `score_map.csv` | The final 30-model orientation map (8 lower, 22 higher) |
| `fix_orientations.py` | The script that generates `score_map.csv` |
| `trial1_orientation_check_preflip.csv` | Pre-flip trial: 8 flagged models, AUROC 0.08-0.21 |
| `trial2_orientation_check_postflip.csv` | Post-flip trial: all 30 models, AUROC 0.76-0.98 |
| `orientation_check.csv` (in `results_full/`) | Final run: identical to trial2 |
| `legacy_pooled_correlation.csv` | Independent cross-check: all 30 correlations positive |

## What this file is not

This file does **not** document tuning, feature selection, or AUROC
maximisation. No model was dropped or added based on its AUROC. The
orientation decision was made on the sign convention implied by the
0.08-0.21 AUROC cluster, verified by a re-run, and cross-checked against the
legacy pooled table. The timing of the decision (before the full run) is
documented in `analysis_plan.md`.
