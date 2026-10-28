"""Donor-level (pseudobulk) differential expression with covariate adjustment.

The donor is the unit of replication. This is a plain OLS on log2-CPM with
classical t-tests: no empirical-Bayes moderation (limma) and no count model
(DESeq2). Use it as a transparent, tested baseline; compare against pydeseq2 /
limma-voom on real data before drawing conclusions.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


def log_cpm(counts, prior: float = 0.5) -> np.ndarray:
    counts = np.asarray(counts, float)
    lib = counts.sum(1, keepdims=True)
    return np.log2((counts + prior) / (lib + 2 * prior) * 1e6)


def size_factors(counts) -> np.ndarray:
    """Median-of-ratios size factors (DESeq2-style) using genes with all counts > 0.
    Robust to a minority of strongly changed genes, unlike plain library size."""
    counts = np.asarray(counts, float)
    ok = (counts > 0).all(0)
    if ok.sum() < 20:
        raise ValueError("too few genes with non-zero counts in every donor for median-of-ratios")
    logc = np.log(counts[:, ok])
    ref = logc.mean(0)
    return np.exp(np.median(logc - ref, axis=1))


def log_norm(counts, prior: float = 0.5) -> np.ndarray:
    """log2 of counts scaled by median-of-ratios size factors, on a CPM-like scale."""
    counts = np.asarray(counts, float)
    scale = size_factors(counts) * counts.sum(1).mean()
    return np.log2((counts + prior) / (scale[:, None] + 2 * prior) * 1e6)


def make_design(meta: pd.DataFrame, term: str = "disease", covariates=()) -> pd.DataFrame:
    cols = [term, *covariates]
    X = meta[cols].astype(float).copy()
    X.insert(0, "intercept", 1.0)
    return X


def fit_de(expr, design: pd.DataFrame, term: str = "disease", genes=None) -> pd.DataFrame:
    X = design.to_numpy(float)
    j = list(design.columns).index(term)
    n, p = X.shape
    if n <= p:
        raise ValueError("need more donors than model parameters")
    XtXi = np.linalg.inv(X.T @ X)
    B = XtXi @ X.T @ expr
    R = expr - X @ B
    s2 = (R ** 2).sum(0) / (n - p)
    se = np.sqrt(s2 * XtXi[j, j])
    t = B[j] / se
    pval = 2 * stats.t.sf(np.abs(t), n - p)
    fdr = multipletests(pval, method="fdr_bh")[1]
    out = pd.DataFrame({"beta": B[j], "se": se, "t": t, "p": pval, "fdr": fdr})
    if genes is not None:
        out.index = genes
    return out


def pseudobulk(X, donors, celltypes, min_cells: int = 10):
    """Sum raw counts per donor within each cell type.

    X: cells x genes (dense or scipy sparse). Returns {celltype: (counts DataFrame
    donors x genes, n_cells Series)}; donors with < min_cells cells are dropped."""
    donors = np.asarray(donors)
    celltypes = np.asarray(celltypes)
    out = {}
    for ct in np.unique(celltypes):
        rows, ids, ncells = [], [], []
        for d in np.unique(donors[celltypes == ct]):
            idx = np.where((celltypes == ct) & (donors == d))[0]
            if len(idx) < min_cells:
                continue
            sub = X[idx]
            s = np.asarray(sub.sum(0)).ravel() if sparse.issparse(sub) else sub.sum(0)
            rows.append(s)
            ids.append(d)
            ncells.append(len(idx))
        if rows:
            out[ct] = (pd.DataFrame(np.vstack(rows), index=ids), pd.Series(ncells, index=ids))
    return out
