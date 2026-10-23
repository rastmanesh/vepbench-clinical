"""Simulated benchmark for tests, examples and the bundled test data.

The simulation has known structure: genes differ in their fraction of pathogenic
variants, and four artificial predictors carry between-gene signal, within-gene
signal, both, or neither. It contains no real biological data.
"""
import numpy as np
import pandas as pd


def simulate_benchmark(n_genes=200, min_variants=12, mean_extra=40, prior_a=0.6, prior_b=0.6,
                       effect=1.0, noise=1.0, seed=0):
    """Return a DataFrame with columns gene, variant, label and four predictor scores.

    within_signal: informative only inside genes.
    gene_prior:    informative only between genes (tracks each gene's pathogenic fraction).
    mixed:         half of each.
    random:        pure noise.
    """
    rng = np.random.default_rng(seed)
    sizes = min_variants + rng.poisson(mean_extra, size=n_genes)
    prior = rng.beta(prior_a, prior_b, size=n_genes)
    gene = np.repeat(np.arange(n_genes), sizes)
    p = prior[gene]
    label = (rng.random(gene.size) < p).astype(int)
    z_gene = (np.log(prior / (1 - prior)) - np.log(prior / (1 - prior)).mean()) / np.log(prior / (1 - prior)).std()
    within = effect * label + noise * rng.standard_normal(gene.size)
    between = z_gene[gene] + noise * rng.standard_normal(gene.size)
    mixed = 0.5 * within + 0.5 * between
    rand = rng.standard_normal(gene.size)
    variant = ["%s:v%d" % ("G%04d" % g, i) for i, g in enumerate(gene)]
    return pd.DataFrame(dict(gene=["G%04d" % g for g in gene], variant=variant, label=label,
                             within_signal=within, gene_prior=between, mixed=mixed, random=rand))
