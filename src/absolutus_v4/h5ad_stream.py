"""Pseudobulk from a downloaded CELLxGENE source .h5ad by streaming the CSR matrix in row blocks.

Memory stays bounded (a block of a few million non-zeros plus a groups x genes accumulator) and the
read is sequential, so cost = file size / disk speed. Needs `h5py` and `anndata` (for element readers)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .pseudobulk_stream import META_COLS, add_entries, check_exact, group_cells_by_class


def _read_elem():
    try:
        from anndata.io import read_elem
    except ImportError:  # older anndata
        from anndata.experimental import read_elem
    return read_elem


def obs_columns(path) -> list[str]:
    import h5py
    with h5py.File(path, "r") as f:
        return [k for k in f["obs"].keys() if not k.startswith("_")]


def inspect(path) -> dict:
    """Layout summary used to decide how to read a file (no matrix data is loaded)."""
    import h5py
    with h5py.File(path, "r") as f:
        use_raw = "raw" in f and "X" in f["raw"]
        X = f["raw/X"] if use_raw else f["X"]
        enc = X.attrs.get("encoding-type", "dense array" if isinstance(X, h5py.Dataset) else "?")
        shape = tuple(X.attrs["shape"]) if "shape" in X.attrs else tuple(X.shape)
        return {"counts_in": "raw/X" if use_raw else "X", "encoding": str(enc), "shape": shape,
                "obs_columns": [k for k in f["obs"].keys() if not k.startswith("_")]}


def pseudobulk_h5ad(path, tissue, class_labels, *, donor_col="donor_id", meta_cols=META_COLS,
                    progress=None, block_entries=8_000_000, dtype=np.float32):
    """class_labels: {class: [cell_type labels]}. One sequential pass over the file.
    Returns (dense groups x genes accumulator, group table, gene table, info dict)."""
    import h5py
    read_elem = _read_elem()
    label_to_class = {l: c for c, ls in class_labels.items() for l in ls}
    with h5py.File(path, "r") as f:
        use_raw = "raw" in f and "X" in f["raw"]
        X = f["raw/X"] if use_raw else f["X"]
        if str(X.attrs.get("encoding-type", "")) != "csr_matrix":
            raise ValueError(f"matrix is {X.attrs.get('encoding-type', 'dense')}, expected csr_matrix")
        n_obs, n_vars = (int(v) for v in X.attrs["shape"])
        cols = ["tissue", "cell_type", donor_col, *meta_cols]
        obs = pd.DataFrame({c: np.asarray(read_elem(f["obs"][c])) for c in dict.fromkeys(cols)})
        cls = obs["cell_type"].astype(str).map(label_to_class)
        keep = ((obs["tissue"].astype(str) == tissue) & cls.notna()).to_numpy()
        if not keep.any():
            return None
        table, codes = group_cells_by_class(obs[keep].reset_index(drop=True), cls[keep].to_numpy(),
                                            key=donor_col, meta_cols=meta_cols)
        row_group = np.full(n_obs, -1, np.int64)
        row_group[np.where(keep)[0]] = codes
        var = read_elem(f["raw/var" if use_raw else "var"])
        genes = pd.DataFrame({"feature_id": var.index.astype(str)})
        genes["feature_name"] = var["feature_name"].astype(str).to_numpy() if "feature_name" in var else genes["feature_id"]
        assert len(genes) == n_vars, "var length does not match the matrix"
        indptr = X["indptr"][:].astype(np.int64)
        data_ds, idx_ds = X["data"], X["indices"]
        dense = np.zeros(len(table) * n_vars, dtype=dtype)
        total, a, checked = int(indptr[-1]), 0, False
        while a < n_obs:
            b = min(max(int(np.searchsorted(indptr, indptr[a] + block_entries, side="right")) - 1, a + 1), n_obs)
            rg = row_group[a:b]
            if (rg >= 0).any():
                lo, hi = int(indptr[a]), int(indptr[b])
                d, ix = data_ds[lo:hi], idx_ds[lo:hi]
                if not checked:
                    if not np.all(d == np.round(d)):
                        raise ValueError(f"{'raw/X' if use_raw else 'X'} holds non-integer values (normalised?), not counts")
                    checked = True
                g = rg[np.repeat(np.arange(b - a), np.diff(indptr[a:b + 1]))]
                m = g >= 0
                add_entries(dense, g[m], ix[m], d[m], n_vars)
            if progress:
                progress(indptr[b] / max(total, 1))
            a = b
        dense = dense.reshape(len(table), n_vars)
        check_exact(dense, dtype)
    return dense, table, genes, {"counts_in": "raw/X" if use_raw else "X", "n_obs": n_obs, "n_vars": n_vars,
                                 "cells_used": int(keep.sum())}


def _col(path, name) -> np.ndarray:
    import h5py
    read_elem = _read_elem()
    with h5py.File(path, "r") as f:
        return np.asarray(read_elem(f["obs"][name])).astype(str)


def observed_tissues(path) -> pd.Series:
    return pd.Series(_col(path, "tissue")).value_counts()


def observed_labels_h5ad(path, tissue) -> pd.Series:
    """cell_type -> n_cells within one tissue (metadata only)."""
    t, c = _col(path, "tissue"), _col(path, "cell_type")
    return pd.Series(c[t == tissue]).value_counts()
