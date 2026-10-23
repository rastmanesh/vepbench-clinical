# Analysis plan (fill in and commit BEFORE the full run; do not edit afterwards)

Date written: 2026-10-23

Data: ProteinGym version 1.3 (Zenodo 10.5281/zenodo.15293562); files and SHA-256 in hashes.txt
  clinical_ProteinGym_substitutions.zip  (3,356,697 bytes; 2,525 per-protein label files)
  zero_shot_clinical_substitutions_scores.zip  (2,525 wide per-protein score files, 30 model columns)
  Preprocessing: labels zip extracted in-memory; each score file split into
  scores_reshaped/<model>/<protein>.csv with columns mutant, score.
  Label column DMS_bin_score mapped Pathogenic->1, Benign->0.
  No rows were filtered; no variants were dropped by hand.

Code: vepbench version 0.1.2 (editable install from vepbench_v0.1.2/vepbench_pkg)
  Two local patches applied before any real-data run:
   (1) analysis/run_clinical_analysis.py lines 120 and 166:
       os.path.relpath(...) -> os.path.basename(...)
       (cross-drive ValueError on Windows; both calls were informational
        "file"/"example_file" audit fields)
  39/39 unit tests pass after the patches.

Models analysed (and any excluded in advance, with the reason):
  All 30 ProteinGym clinical baselines present in the score files are analysed;
  none excluded in advance. Sparse coverage (e.g. MPC 1820/2525 proteins,
  MutPred 2307/2525) is a property of the source data and is reported via the
  coverage column, not used to exclude a model.

Orientation of each model and its documentation source: (separate notes file)
  Provisional orientation was set from published tool conventions. A trial run
  (--n-boot 0, pooled AUROC only) showed 8 ProteinGym-provided score columns with
  the opposite sign: EVE, FATHMM, LRT, PoET, Provean, SIFT, SIFT4G, TranceptEVE_L.
  These 8 were flipped. The decision is a sign-convention decision, not an
  AUROC-maximisation decision: each flipped column moved from AUROC ~0.1 to
  ~0.8-0.9 (i.e. to 1 - AUROC), the signature of a sign error. The five dbNSFP
  tools (SIFT, SIFT4G, Provean, FATHMM, LRT) agree with the legacy pooled table's
  signs, which is an independent cross-check. Final orientation map is in
  score_map.csv, with per-model sources in orientation_notes.md.

Parameters: min-per-class = 5 ; bootstrap resamples = 1000 ; seed = 0 ;
  specificity = 0.95 ; max FPR = 0.05

Primary question:
  Do the 30 ProteinGym clinical baselines discriminate pathogenic from benign
  variants (a) pooled across genes, (b) within genes, and (c) beyond what gene
  identity alone predicts?

Primary outcome (one number, with its interval):
  Pooled AUROC over all variants with a non-missing score for that model, with a
  percentile 95% bootstrap interval over genes (bootstrap resamples=1000, seed=0).
  The 30 per-model pooled AUROCs and their intervals are the primary table.

Rule for calling a ranking change material (write the threshold now, say how it was chosen):
  A ranking change between the pooled-AUROC and per-gene-AUROC model orderings is
  material if either (a) Kendall's tau < 0.6, or (b) the top-two models by pooled
  AUROC differ from the top-two by per-gene AUROC. Chosen now, before the full run:
  (a) is reachable with 30 models (it was not with 4 models in the simulated test),
  and (b) is the threshold that matters to a reader who asks "which model is best".
  The rule was fixed after the pooled-only trial run and before the full run with
  intervals, and that timing is reported in the Methods.

Secondary outcomes:
  Per-gene AUROC (mean, weighted mean, median) restricted to genes with
  >= 5 positives and >= 5 negatives (80-gene scale on the simulated set; the count
  on the real data is reported in audit/ n_eligible_groups); average precision;
  sensitivity at 95% specificity; partial AUROC (FPR <= 5%, McClish-standardized);
  gene-identity baseline under in-sample, k-fold, leave-one-out and group k-fold
  schemes; pooled-vs-per-gene gap and its decomposition; rank-change table.

How missing scores and models scored on different variants are handled:
  Per-model analysis uses all variants with a non-missing score for that model
  (coverage reported in audit_models.csv and models.csv). The common-variant
  analysis additionally restricts to variants scored by every model
  (common_variants = 16,626 in the trial); it is reported alongside, not instead of,
  the per-model analysis. No imputation is performed.

Planned reporting: all models and all outcomes will be reported, including null or
unfavourable results; any deviation from this plan will be stated.

Signed by (author): [Reza Rastmanesh]