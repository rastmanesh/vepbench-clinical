# vepbench evaluation report

vepbench version 0.1.0; seed 0; bootstrap resamples 200 (genes resampled, not variants).
Variants: 6264 in 120 groups. Labels: 1 = pathogenic, 0 = benign. Minimum variants per class per group: 5.

## Predictors

| model         |   n_scored |   coverage |   pooled_auroc |   pergroup_mean_auroc |   pergroup_weighted_auroc |   pergroup_median_auroc |   n_eligible_groups |   average_precision |   sensitivity_at_specificity |   partial_auroc |   pooled_ci_low |   pooled_ci_high |   pergroup_ci_low |   pergroup_ci_high |
|:--------------|-----------:|-----------:|---------------:|----------------------:|--------------------------:|------------------------:|--------------------:|--------------------:|-----------------------------:|----------------:|----------------:|-----------------:|------------------:|-------------------:|
| within_signal |       6264 |          1 |          0.771 |                 0.762 |                     0.765 |                   0.774 |                  80 |               0.758 |                        0.286 |           0.577 |           0.761 |            0.783 |             0.744 |              0.779 |
| gene_prior    |       6264 |          1 |          0.765 |                 0.522 |                     0.517 |                   0.52  |                  80 |               0.753 |                        0.266 |           0.579 |           0.724 |            0.799 |             0.502 |              0.548 |
| mixed         |       6264 |          1 |          0.852 |                 0.713 |                     0.713 |                   0.709 |                  80 |               0.846 |                        0.453 |           0.647 |           0.83  |            0.871 |             0.695 |              0.732 |
| random        |       6264 |          1 |          0.493 |                 0.484 |                     0.486 |                   0.497 |                  80 |               0.48  |                        0.037 |           0.497 |           0.48  |            0.506 |             0.468 |              0.499 |

## Gene-identity baseline (a predictor that knows only the gene)

     scheme  pooled  pergroup  pergroup_biased
  in_sample   0.885       0.5            0.500
      kfold   0.873       NaN            0.334
        loo   0.873       NaN            0.000
group_kfold   0.474       0.5            0.500

## Pooled versus per-group ranking

{
  "n_models": 4,
  "spearman": 0.7999999999999999,
  "kendall": 0.6666666666666669,
  "discordant_fraction": 0.16666666666666666,
  "max_abs_rank_change": 1.0,
  "n_rank_changed": 2
}

## Bootstrap agreement

{
  "kendall_pooled_vs_pergroup": {
    "mean": 0.5533333333333335,
    "ci_low": 0.3333333333333334,
    "ci_high": 0.6666666666666669
  },
  "spearman_pooled_vs_external": {
    "mean": 0.6639999999999999,
    "ci_low": 0.39999999999999997,
    "ci_high": 0.7999999999999999
  },
  "spearman_pergroup_vs_external": {
    "mean": 1.0,
    "ci_low": 1.0,
    "ci_high": 1.0
  },
  "external_gain_pergroup_minus_pooled": {
    "mean": 0.3360000000000001,
    "ci_low": 0.20000000000000007,
    "ci_high": 0.6000000000000001
  }
}

## Agreement with the external ranking

{
  "pooled": {
    "n_models": 4,
    "spearman": 0.7999999999999999,
    "kendall": 0.6666666666666669,
    "discordant_fraction": 0.16666666666666666,
    "max_abs_rank_change": 1.0,
    "n_rank_changed": 2
  },
  "pergroup": {
    "n_models": 4,
    "spearman": 1.0,
    "kendall": 1.0,
    "discordant_fraction": 0.0,
    "max_abs_rank_change": 0.0,
    "n_rank_changed": 0
  }
}
