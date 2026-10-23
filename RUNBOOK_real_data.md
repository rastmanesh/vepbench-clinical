# Runbook: running vepbench 0.1.2 on the real ProteinGym clinical files

Written for Windows PowerShell (Linux/macOS equivalents are the same commands without `.exe`).
Nothing in this runbook has been run on the real files. Stop and ask if anything looks different from what is described.

## 1. Install and test the code (about 5 minutes)

```
cd vepbench_pkg
pip install -e .
pip install scikit-learn          # used only by the tests
python -m unittest discover -s tests
```
Expect: `Ran 39 tests ... OK`. If not, stop and send the message.

## 2. The data: ProteinGym version 1.3 (the latest at the time of writing)

Two files are needed (names and sizes from the ProteinGym README):

| File | What it holds | Size unzipped |
|---|---|---|
| `clinical_ProteinGym_substitutions.zip` | clinical labels, one file per protein | 123 MB |
| `zero_shot_clinical_substitutions_scores.zip` | clinical model scores | 0.9 GB |

If these are already on your work computer from an earlier run, check that they are version 1.3 (same names and sizes) and re-download if unsure.
Download (PowerShell):
```
curl.exe -o clinical_ProteinGym_substitutions.zip https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/clinical_ProteinGym_substitutions.zip
curl.exe -o zero_shot_clinical_substitutions_scores.zip https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/zero_shot_clinical_substitutions_scores.zip
Expand-Archive clinical_ProteinGym_substitutions.zip -DestinationPath labels_raw
Expand-Archive zero_shot_clinical_substitutions_scores.zip -DestinationPath scores_raw
Get-FileHash *.zip -Algorithm SHA256 | Out-File hashes.txt
```
Keep `hashes.txt`. Cite ProteinGym (NeurIPS 2023) and the version, with its Zenodo record DOI 10.5281/zenodo.15293562.

## 3. Look before you run

Per the ProteinGym README each label file has the columns `mutant`, `mutated_sequence` and `DMS_score_bin` (1 = pathogenic, 0 = benign) and no continuous `DMS_score` column.
List the folders and open one label file and one score file. If the layout is not "labels in one folder, one sub-folder per model with the score files inside",
stop and send me: the folder listing to two levels, and the first line of one label file and one score file.

## 4. Build the score map

```
python analysis/run_clinical_analysis.py inspect --labels-dir <LABELS_FOLDER> --scores-dir <SCORES_FOLDER> --out inspect_out
```
Open `inspect_out/score_map_template.csv`. For every model fill in `column` (the score column) and `orientation`.
**Set orientation from the model's documented convention, never from its results.** Write the source for each model in a notes file.
Likelihood-based scores (protein language models, alignment models) usually mean higher = fitter = more benign, so those are normally `lower_pathogenic`;
for dbNSFP-derived tools the convention differs by tool (SIFT, PROVEAN, FATHMM and LRT are lower = pathogenic; most others are higher = pathogenic).
Your earlier table's signs agree with this for those five tools and are a useful cross-check.

## 5. Trial run (fast, no intervals)

```
python analysis/run_clinical_analysis.py run --labels-dir <LABELS_FOLDER> --scores-dir <SCORES_FOLDER> --score-map score_map.csv --previous-table Final_Paper_Table.csv --n-boot 0 --out trial
```
Read, in this order: `RUN_SUMMARY.json`, `audit_models.csv`, `audit_label_files.csv`, `audit_score_files.csv`, `orientation_check.csv`, `comparison_with_previous_table.csv`.
Explain every skipped model and failed file (the earlier script dropped EVE and ESM1b without saying why).
A model with pooled AUROC below 0.5 may have the wrong orientation: check its documentation. Do not flip it just to improve its result.

## 6. Write the analysis plan BEFORE the full run

Copy `analysis_plan_TEMPLATE.md` to `analysis_plan.md`, fill it in, and commit it with its date. Change nothing after you have seen full results; report any deviation.

## 7. Full run

```
python analysis/run_clinical_analysis.py run --labels-dir <LABELS_FOLDER> --scores-dir <SCORES_FOLDER> --score-map score_map.csv --previous-table Final_Paper_Table.csv --n-boot 1000 --seed 0 --out results_full
```
Runtime is a guess: in my test environment about 0.04 s per resample per model at roughly 95,000 variants, so expect tens of minutes to a couple of hours. Time a run with `--n-boot 200` first.
Optional: `--external dms_by_model.csv` (columns `model,value`, for example mean Spearman on the ProteinGym DMS benchmark, model names matching the score map).

## 8. What to send back (no raw ProteinGym files are needed)

`RUN_SUMMARY.json`; `audit_models.csv`; `audit_label_files.csv`; `audit_score_files.csv` (or only its failed rows if it is long); `orientation_check.csv`;
`comparison_with_previous_table.csv`; from `all_scored/` and `common_variants/`: `models.csv`, `aggregation_comparison.csv`, `gene_identity_baseline.csv`, `summary.json`;
your `analysis_plan.md`; the orientation notes; `hashes.txt`; the exact commands you ran; the vepbench version.

## 9. If something goes wrong

- Coverage near zero or `variants` far below the label total: protein names do not match between label and score files. See `protein_source` in `audit_label_files.csv`; use `--protein-col` if the column has another name.
- "cannot choose the score column": fill `column` in the score map for that model.
- Models as columns inside per-protein files (a flat folder): not supported; send me the header line and the listing.
- Out-of-memory: tell me the machine's memory and the error.

## 10. How to read the results

Many ProteinGym clinical baselines come from dbNSFP 4.4a (per the ProteinGym README), and several are supervised tools trained on ClinVar-type data, so do not call them zero-shot
and do not treat their scores on ClinVar-derived labels as independent evidence. The ProteinGym leaderboard already averages AUC by protein; what this analysis adds is the pooled
comparison on the same proteins, the gene-identity baseline, and how much rankings move.
