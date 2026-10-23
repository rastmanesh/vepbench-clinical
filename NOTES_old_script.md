# Problems found in the first analysis script (run_final_analysis.py) and how the new script differs

Found by reading the script in the repository rastmanesh/ProteinGym-Clinical-Benchmark:

1. Every file and model error is swallowed (`except Exception: pass`), so a model that fails to parse disappears from the table.
   EVE and ESM1b are handled specially in the script but absent from the result table; the script cannot say why.
   -> New script: every file and model gets a status and message in an audit file.
2. For models other than EVE and ESM1b the score column is taken to be the last column of the file.
   -> New script: the column comes from a score map, or the file must have exactly one numeric column besides protein and mutant.
3. The label column is the first column whose name contains "dms"; this may not be the binary label.
   -> New script: the label column is named (default DMS_score_bin) and must contain only 0 and 1.
4. The metric is a pooled Spearman correlation with the binary label, reported as an absolute value, which hides wrong orientation or a wrongly parsed column.
   -> New script: signed values, an orientation check, and the vepbench evaluation (pooled and per-protein AUROC, gene-cluster bootstrap intervals, gene-identity baseline).
5. Models are scored on different variant subsets (N ranges widely) without saying so.
   -> New script: results on all scored variants and on the variants scored by every model, with coverage reported.
6. Paths are hard-coded to a local Windows drive, and the README has one line.
   -> New script: paths are arguments; README documents the steps.
