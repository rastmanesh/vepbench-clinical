[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23589773.svg)](https://doi.org/10.5281/zenodo.23589773)
# vepbench (version 0.1.2)

Gene-aware evaluation of variant effect predictors on labelled missense variants. When variants
cluster in genes, a predictor can score well by learning which genes tend to harbour pathogenic
variants, not which variants are damaging. `vepbench` reports pooled and within-gene performance side by
side, adds a gene-identity baseline, puts bootstrap intervals on both by resampling genes (not variants),
and measures how much a predictor ranking changes when the evaluation design changes.

**Status:** Validated on the ProteinGym v1.3 clinical set (2,525 proteins, 30 models).

## Install

```
pip install -e .            # needs numpy, pandas, scipy
pip install -e ".[test]"    # adds scikit-learn, used only by the tests
```

## Quick start

```
sh examples/run_example.sh          # bundled simulated data, no real biological data
python -m unittest discover -s tests
```

Command line:

```
python -m vepbench evaluate --input table.csv --group gene --label label \
    --scores modelA,modelB,modelC --n-boot 200 --out results
```

Python:

```python
import vepbench as vb
df = vb.simulate_benchmark(n_genes=200, seed=0)
agg = vb.compare_aggregations(df, "gene", "label", ["within_signal", "gene_prior", "mixed", "random"])
base = vb.gene_identity_report(df, "gene", "label")
```

## Input

One row per variant, in a CSV, TSV or Parquet file (Parquet needs pyarrow): a group column (gene or protein), a label
column (1 = pathogenic, 0 = benign; or any value given with `--positive-label`) and one numeric column per
predictor. Higher scores must mean more pathogenic; list predictors where lower means pathogenic with `--flip`.
Missing scores are allowed: variants a predictor could not score are dropped for that predictor only, and the fraction scored is reported as `coverage`.

## Output files (in `--out`)

`models.csv` (per predictor: pooled AUROC, mean per-group AUROC, weighted mean, median, average precision,
sensitivity at fixed specificity, partial AUROC, coverage and bootstrap intervals), `per_group_auroc.csv`,
`aggregation_comparison.csv` (pooled versus per-group values and ranks), `gene_identity_baseline.csv`,
`summary.json` and `report.md`.

## Parameters

| Option | Default | Purpose and effect |
|---|---|---|
| `--input` | required | Variant table. |
| `--group` | required | Column naming the gene or protein; groups are the resampling units. |
| `--label` | required | Binary label column. |
| `--scores` | required | Comma-separated predictor columns. |
| `--positive-label` | none | Label value meaning pathogenic when labels are not 0/1. |
| `--flip` | none | Predictors where a lower score means pathogenic (scores are negated). |
| `--min-per-class` | 5 | A group enters the per-group AUROC only if it has at least this many variants of each class. Raising it keeps fewer, more precise groups; results can depend on it, so report it. |
| `--n-boot` | 200 | Gene-level bootstrap resamples; 0 skips intervals and the ranking bootstrap. More resamples give smoother intervals and take longer. |
| `--seed` | 0 | Random seed for the bootstrap and the baseline folds. |
| `--specificity` | 0.95 | Specificity at which sensitivity is reported. |
| `--max-fpr` | 0.05 | False-positive-rate limit for the standardized partial AUROC (0.5 is chance). |
| `--n-splits` | 5 | Folds for the gene-identity baseline. |
| `--external` | none | CSV with columns `model,value`: an independent per-predictor ranking (for example from deep mutational scanning). Reports how well pooled and per-group rankings each agree with it. |
| `--out` | required | Output directory. |

## Definitions

- **Pooled AUROC:** one AUROC over all scored variants (Mann-Whitney form; ties count one half).
- **Pooled AUROC on eligible genes:** the pooled AUROC restricted to the genes that qualify for the per-group summary. Comparing the all-variant pooled value with the per-group mean mixes two effects, so the report splits the gap into `gap_selection` (different variants) and `gap_aggregation` (pooling across genes).
- **Per-group AUROC:** AUROC within each eligible gene; the summary is the unweighted mean (also a weighted mean, with weights n_pos*n_neg/(n_pos+n_neg), and the median).
- **Gene-identity baseline:** each variant is scored by the smoothed fraction of pathogenic variants in its gene (additive smoothing toward the overall rate, alpha = 1). Schemes: `in_sample` (leaky upper bound), `kfold` (random variant-level split), `loo` (leave one variant out), `group_kfold` (split by gene; carries no gene information).
- **Bootstrap:** percentile intervals from resampling whole genes with replacement; the same resamples are used for every predictor, so comparisons are paired.
- **Ranking agreement:** Spearman and Kendall correlations, fraction of predictor pairs ordered differently, and rank changes, between pooled and per-group rankings (and an optional external ranking).

### A caution about the baseline

With `kfold` and `loo`, leaving out a variant's own label (or fold) from the gene rate makes the score
negatively related to the label *within* a gene. The within-gene AUROC of these baselines is therefore
biased below 0.5 (exactly 0 for `loo` in the example). `gene_identity_report` gives a per-group value only for `in_sample`
and `group_kfold` (0.5 by construction); the biased values are shown separately as `pergroup_biased`.
Use the pooled values to judge how much of a benchmark gene identity can explain.

## Clinical benchmarks laid out as in ProteinGym

`analysis/run_clinical_analysis.py` runs the evaluation on a folder of per-protein label files plus one sub-folder of
score files per model. **It has been tested only on a mock folder structure (tests/test_clinical_script.py), never on the real files.**
Check the layout against your download: the labels folder holds one CSV per protein with a mutant column and a 0/1 label column
(default `DMS_score_bin`); the scores folder holds one sub-folder per model, named for the model, each with CSVs that have a mutant column and one score column.
If the protein column is missing, the file name is used as the protein.

```
python analysis/run_clinical_analysis.py inspect --labels-dir LABELS --scores-dir SCORES --out inspect_out
# open inspect_out/score_map_template.csv; fill 'column' where blank; set 'orientation' to lower_pathogenic where needed
python analysis/run_clinical_analysis.py run --labels-dir LABELS --scores-dir SCORES --score-map score_map.csv \
    --previous-table Final_Paper_Table.csv --out results
```

What it guarantees: nothing is skipped silently (`audit_label_files.csv`, `audit_score_files.csv`, `audit_models.csv` give a status and message
for every file and model); the score column is never guessed (a model without a mapped column is skipped unless its files have exactly one numeric
column besides protein and mutant); models with a pooled AUROC below 0.5 are listed in `orientation_check.csv`; results are given for each model on all
variants it scored (`all_scored/`) and on the variants scored by every model (`common_variants/`); `legacy_pooled_spearman.csv` reproduces the earlier pooled
Spearman table and `comparison_with_previous_table.csv` lists differences in N and correlation and any model present in only one table.
Read `RUN_SUMMARY.json` first: it counts proteins, variants, models analysed and models skipped.

Note that ProteinGym's clinical files are one per protein and its leaderboard aggregates by protein, so the per-protein results reproduce a design the benchmark
already uses; the pooled comparison and the gene-identity baseline are what add information. Many models in the clinical files may be supervised tools trained on
ClinVar-type data; do not call them zero-shot, and do not interpret their scores on ClinVar-derived labels as independent evidence.

## Simulated example data

Everything in `examples/` is simulated. The four example predictors (`within_signal`, `gene_prior`, `mixed`, `random`) are built from the simulation's own labels and gene-level pathogenic fractions, so results on them show that the software behaves as designed and say nothing about real predictors or real variants. Never report them as findings.

## Limitations

- Binary labels and one score per variant only. No calibration or evidence-strength mapping, no homology-aware splitting, no structural or conservation stratification, no review-status handling (all planned or left to the user).
- Intervals are percentile bootstrap intervals; with few genes they can be unreliable.
- The per-group summary depends on `--min-per-class` and on how many genes qualify.
- Tested on simulated data and against scikit-learn for AUROC, average precision and partial AUROC; not yet validated on real benchmarks. The ProteinGym-style script is untested on real files.
- Runtime in a single test environment (about 95,000 variants, 2,500 genes, 12 predictors): about 15 s without intervals and about 60 s with 100 resamples; roughly linear in predictors and resamples.
- The name may collide with Ensembl's Variant Effect Predictor ("VEP"); consider renaming before release.

## Release checklist (for the author)

Fill in the copyright holder in LICENSE; place the code in a public repository; archive the release with a DOI (for example Zenodo) and put that URL in the manuscript's Availability and implementation section; confirm nothing named here already exists on PyPI or GitHub; read and run the code yourself and describe any AI help as in AI_ASSISTANCE.md.

## Licence

MIT (see LICENSE).
