"""vepbench: gene-aware evaluation of variant effect predictors.

Functions for evaluating predictors on labelled missense variants when variants
are clustered in genes or proteins: per-gene and pooled AUROC, gene-identity
baselines, cluster-aware bootstrap intervals, and ranking-stability analysis.
"""
from .metrics import (auroc, average_precision, partial_auroc,
                      sensitivity_at_specificity, specificity_at_sensitivity)
from .genelevel import (within_group_auroc, per_group_auroc, summarize_per_group,
                        compare_aggregations, gene_identity_baseline,
                        gene_identity_report)
from .bootstrap import (bootstrap_within_group_auroc, percentile_ci, bootstrap_statistic,
                        bootstrap_pergroup_mean, paired_bootstrap_pergroup)
from .ranking import rank_models, compare_rankings, bootstrap_rank_agreement, split_half_rank_reliability
from .simulate import simulate_benchmark
from .pair_diagnostics import (between_group_auroc, pair_composition,
                               protein_mean_vs_residual_auroc)

__version__ = "0.1.4"

__all__ = [
    "auroc", "average_precision", "partial_auroc",
    "sensitivity_at_specificity", "specificity_at_sensitivity",
    "within_group_auroc", "per_group_auroc", "summarize_per_group",
    "compare_aggregations", "gene_identity_baseline", "gene_identity_report",
    "bootstrap_within_group_auroc", "percentile_ci", "bootstrap_statistic",
    "bootstrap_pergroup_mean", "paired_bootstrap_pergroup",
    "rank_models", "compare_rankings", "bootstrap_rank_agreement",
    "split_half_rank_reliability", "simulate_benchmark",
    "between_group_auroc", "pair_composition", "protein_mean_vs_residual_auroc",
]