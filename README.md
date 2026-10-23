# vepbench-clinical
Gene-aware benchmarking pipeline for ProteinGym clinical variant effect predictors. Author: Reza Rastmanesh
# vepbench-clinical

**A gene-aware benchmarking pipeline for ProteinGym clinical variant effect predictors.**

This repository contains the code and analysis pipeline for benchmarking 30 ProteinGym clinical baseline models. It evaluates pooled variant-level performance versus per-gene performance, controls for gene-identity confounding, and audits predictor orientation.

## Pipeline Overview
1. **Data ingestion**: Extracts per-protein scores and labels from ProteinGym v1.3 clinical set.
2. **Orientation check**: Audits 30 models for sign inversion (8 models require inversion).
3. **Evaluation**: Computes pooled AUROC, per-gene AUROC (minimum 5 positive and 5 negative variants per gene), and bootstrap confidence intervals.
4. **Gene-identity baseline**: Evaluates a target-encoding baseline using k-fold, leave-one-out (LOO), and group k-fold validation.

## Installation
```bash
pip install -e .
