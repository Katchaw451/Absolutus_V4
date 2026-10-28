"""Memory-bounded pseudobulk aggregation from streamed sparse (COO) chunks.

Pure functions (no Census dependency) so they can be tested exhaustively on synthetic data.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

OBS_COLS = ["soma_joinid", "donor_id", "disease", "sex", "development_stage", "assay"]
META_COLS = ("disease", "sex", "development_stage", "assay")


def _lit(s: str) -> str:
    if any(ch in s for ch in "'\"\\"):
        raise ValueError(f"unsafe character in filter literal: {s!r}")
    return f"'{s}'"


def value_filter(dataset_ids, tissue: str, labels) -> str:
    ids = ", ".join(_lit(i) for i in dataset_ids)
    labs = ", ".join(_lit(l) for l in labels)
    return f"dataset_id in [{ids}] and tissue == {_lit(tissue)} and cell_type in [{labs}]"


def group_cells(obs: pd.DataFrame, key: str = "donor_id", meta_cols=META_COLS):
    """One group per donor. Returns (donor table, group index per row of `obs`).
    Metadata that varies within a donor is joined with '; ' (and therefore visible)."""
    codes, uniques = pd.factorize(obs[key], sort=True)
    table = {key: np.asarray(uniques), "n_cells": np.bincount(codes, minlength=len(uniques))}
    for c in meta_cols:
        agg = obs[c].astype(str).groupby(codes).agg(lambda s: "; ".join(sorted(set(s))))
        table[c] = agg.to_numpy()
    return pd.DataFrame(table), codes


def accumulate_chunks(chunks, obs_joinids, group_idx, n_groups, var_joinids):
    """Sum COO chunks into an (n_groups x n_vars) CSR matrix.

    chunks: iterable of (obs soma_joinid, var soma_joinid, value) array triples.
    obs_joinids/group_idx: the query's cells and the group of each (aligned).
    var_joinids: defines the column order of the result (column i <-> var_joinids[i]).
    A chunk referring to a cell or gene outside the query raises: it would mean the
    query and the stream disagree."""
    obs_joinids = np.asarray(obs_joinids, np.int64)
    if len(obs_joinids) == 0:
        raise ValueError("empty query")
    o = np.argsort(obs_joinids, kind="stable")
    obs_sorted, grp_sorted = obs_joinids[o], np.asarray(group_idx)[o]
    var_joinids = np.asarray(var_joinids, np.int64)
    v = np.argsort(var_joinids, kind="stable")
    var_sorted, n_vars = var_joinids[v], len(var_joinids)
    total = sparse.csr_matrix((n_groups, n_vars), dtype=np.float64)
    for d0, d1, val in chunks:
        d0, d1 = np.asarray(d0, np.int64), np.asarray(d1, np.int64)
        if len(d0) == 0:
            continue
        pos = np.minimum(np.searchsorted(obs_sorted, d0), len(obs_sorted) - 1)
        if not np.array_equal(obs_sorted[pos], d0):
            raise ValueError("stream contains cells outside the query")
        vpos = np.minimum(np.searchsorted(var_sorted, d1), n_vars - 1)
        if not np.array_equal(var_sorted[vpos], d1):
            raise ValueError("stream contains genes outside the var table")
        part = sparse.coo_matrix((np.asarray(val, np.float64), (grp_sorted[pos], v[vpos])),
                                 shape=(n_groups, n_vars)).tocsr()
        total = total + part
    return total


def save_group(stem, counts, donors: pd.DataFrame) -> None:
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    sparse.save_npz(f"{stem}.counts.npz", counts)
    donors.to_csv(f"{stem}.donors.csv", index=False)


def load_group(stem):
    stem = Path(stem)
    return sparse.load_npz(f"{stem}.counts.npz").tocsr(), pd.read_csv(f"{stem}.donors.csv")


def group_cells_by_class(obs: pd.DataFrame, classes, key: str = "donor_id", meta_cols=META_COLS):
    """One group per (cell class, donor). Returns (table with a 'cell_class' column plus the
    donor columns of `group_cells`, group index per row of `obs`)."""
    combo = np.char.add(np.char.add(np.asarray(classes, dtype=str), "\t"), obs[key].astype(str).to_numpy())
    codes, uniques = pd.factorize(combo, sort=True)
    split = [u.split("\t", 1) for u in uniques]
    table = {"cell_class": np.array([s[0] for s in split]), key: np.array([s[1] for s in split]),
             "n_cells": np.bincount(codes, minlength=len(uniques))}
    for c in meta_cols:
        agg = obs[c].astype(str).groupby(codes).agg(lambda s: "; ".join(sorted(set(s))))
        table[c] = agg.to_numpy()
    return pd.DataFrame(table), codes


F32_EXACT = 2 ** 24


def add_entries(dense_flat, groups, cols, vals, n_vars: int) -> None:
    """dense_flat[group * n_vars + col] += val (vectorised; duplicates are summed)."""
    np.add.at(dense_flat, np.asarray(groups, np.int64) * n_vars + np.asarray(cols, np.int64),
              np.asarray(vals, dense_flat.dtype))


def check_exact(dense, dtype) -> None:
    """float32 accumulators are exact for integer counts below 2**24; refuse to continue otherwise."""
    if np.dtype(dtype) == np.float32 and dense.size and float(dense.max()) >= F32_EXACT:
        raise ValueError("a pseudobulk sum reached 2**24; rerun with dtype=float64")
