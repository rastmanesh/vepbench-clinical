# Changelog

## 0.1.4
- New: `pair_composition`, `between_group_auroc`, `between_minus_within`, `bootstrap_decomposition`; the models table gains `within_pair_share`, `between_pairwise_auroc`, `gap_between_minus_within` and their intervals.
- New: `analysis/reconcile_with_first_run.py`, a regression check against an earlier run.
- Decision: no protein-mean versus residual AUROC output (see README); a test documents why.

## 0.1.3
- New: `within_group_auroc` (pairwise within-group AUROC using all groups with both classes), its cluster bootstrap, and columns `within_pairwise_auroc` / `rank_within_pairwise`.
- New: `split_half_rank_reliability` (noise ceiling for per-group rankings); the clinical script reports it in each `summary.json`.
- New: `ranking_pooled_eligible_vs_pergroup` and `ranking_pooled_vs_within_pairwise` in the evaluation summary.
- Fix: audit paths no longer raise `ValueError` when the data are on a different Windows drive from the working directory (previously patched by hand with `os.path.basename`).
- Fix: SciPy results are read by index, not `.statistic`, so older SciPy works; report tables fall back to plain text when `tabulate` is missing; minimum Python lowered to 3.8 (checked for 3.8 grammar only; tested on Python 3.12).

## 0.1.2
- Pooled AUROC on the eligible genes and the split of the pooled-minus-per-group gap into selection and aggregation parts; single-class partial AUROC returns NaN.

## 0.1.1 and 0.1.0
- Initial release and the audited clinical analysis script.
