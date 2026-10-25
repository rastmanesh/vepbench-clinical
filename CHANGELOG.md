# Changelog

## 0.1.4
- New: `between_group_auroc` — pairwise between-group AUROC (complement of `within_group_auroc` on pair counts).
- New: `pair_composition` — total / within / between pair counts and shares per model.
- New: `protein_mean_vs_residual_auroc` — gene-reliance gap diagnostic (protein-mean AUROC minus protein-residual AUROC). A pure variant-level predictor has a small gap; a pure gene-identity predictor has a large gap.
- New columns in `evaluate()` output: `between_group_auroc`, `protein_mean_auroc`, `protein_residual_auroc`, `gap_mean_minus_residual`, `total_pairs`, `within_pairs`, `between_pairs`, `within_share`, `between_share`.
- Docs: `pair_diagnostics.py` module docstring explains that between-group pairs are NOT gene-identity comparisons — a common misinterpretation.

## 0.1.3
- **Fix (critical): orientation map corrected.** `score_map.csv` in v0.1.3 has 4 `lower_pathogenic` models (ESM1b, EVE, GEMME, TranceptEVE_L), not 8. The 10/23/2026 10:10 AM edit of the map (which marked FATHMM, LRT, PoET, Provean, SIFT, SIFT4G as `lower_pathogenic` in addition to the 4) was empirically wrong: raw AUROC on the full 62,727-variant set is > 0.5 for all 6, meaning they were already correctly oriented. The corrected map reproduces v0.1.2's published pooled AUROC for all 30 models to 6 decimal places. SHA-256: CE333240CB453BF5222A2A051FA9CE68789F8B02A3946713A6B2338C3E9CFE79.
- **Fix: provenance.** `score_map.csv` timestamped 10/23/2026 10:10 AM — after the v0.1.2 full run (10/23/2026 01:13 AM) — was flagged by hostile review. The map is now frozen before the reported run and hashed into the v0.1.3 release.
- New: `within_group_auroc` (pairwise within-group AUROC using all groups with both classes), its cluster bootstrap, and columns `within_pairwise_auroc` / `rank_within_pairwise`.
- New: `split_half_rank_reliability` (noise ceiling for per-group rankings); the clinical script reports it in each `summary.json`.
- New: `ranking_pooled_eligible_vs_pergroup` and `ranking_pooled_vs_within_pairwise` in the evaluation summary.
- Fix: audit paths no longer raise `ValueError` when the data are on a different Windows drive from the working directory (previously patched by hand with `os.path.basename`).
- Fix: SciPy results are read by index, not `.statistic`, so older SciPy works; report tables fall back to plain text when `tabulate` is missing; minimum Python lowered to 3.8 (checked for 3.8 grammar only; tested on Python 3.12).

## 0.1.2
- Pooled AUROC on the eligible genes and the split of the pooled-minus-per-group gap into selection and aggregation parts; single-class partial AUROC returns NaN.

## 0.1.1 and 0.1.0
- Initial release and the audited clinical analysis script.


