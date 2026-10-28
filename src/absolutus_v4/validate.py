"""Validation statistics: held-out enrichment with degree matching, replication."""
from __future__ import annotations

import numpy as np
from scipy.stats import binomtest

from .network import degree_bins


def degree_matched_enrichment(selected, validation, deg, n_perm=1000, n_bins=10, rng=None):
    """Overlap of `selected` with a held-out `validation` gene set, against a null that
    draws random gene sets with the same degree-bin composition as `selected`.
    Returns (observed, null_mean, empirical one-sided p with +1 correction)."""
    rng = rng or np.random.default_rng(0)
    deg = np.asarray(deg)
    bins = degree_bins(deg, n_bins)
    val = np.zeros(len(deg), bool)
    val[list(validation)] = True
    sel = np.asarray(list(selected), int)
    obs = int(val[sel].sum())
    plan = [(np.where(bins == b)[0], int((bins[sel] == b).sum())) for b in np.unique(bins)]
    plan = [(idx, c) for idx, c in plan if c > 0]
    null = np.empty(n_perm)
    for j in range(n_perm):
        null[j] = sum(int(val[rng.choice(idx, c, replace=False)].sum()) for idx, c in plan)
    return obs, float(null.mean()), (1 + int((null >= obs).sum())) / (n_perm + 1)


def replication_concordance(beta_disc, beta_rep, module):
    """Fraction of module genes whose effect sign agrees between discovery and replication
    cohorts, with a one-sided binomial test against 0.5."""
    m = list(module)
    agree = int((np.sign(np.asarray(beta_disc)[m]) == np.sign(np.asarray(beta_rep)[m])).sum())
    return agree / len(m), float(binomtest(agree, len(m), 0.5, alternative="greater").pvalue)
