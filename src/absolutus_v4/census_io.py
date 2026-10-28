"""CELLxGENE Census access. Needs network + `cellxgene-census`; imports are lazy so the rest of
the package (and CI) works without them. NOT covered by unit tests: validate with
`python scripts/build_pseudobulk.py --smoke`, which cross-checks against in-memory AnnData."""
from __future__ import annotations

import numpy as np
import pandas as pd

from scipy import sparse

from .pseudobulk_stream import OBS_COLS, accumulate_chunks, group_cells, group_cells_by_class, value_filter


def open_census(version: str):
    """open_soma with longer S3 timeouts (TileDB aborts a read that stalls for 3 s by default)."""
    import cellxgene_census
    cfg = {"vfs.s3.request_timeout_ms": 120000, "vfs.s3.connect_timeout_ms": 30000}
    try:
        return cellxgene_census.open_soma(census_version=version, tiledb_config=cfg)
    except TypeError:
        return cellxgene_census.open_soma(census_version=version)


def observed_labels(human, dataset_ids, tissue: str) -> pd.Series:
    """cell_type -> n_cells for the given datasets/tissue (metadata only)."""
    import collections
    ids = ", ".join(f"'{i}'" for i in dataset_ids)
    acc = collections.Counter()
    for tbl in human.obs.read(value_filter=f"dataset_id in [{ids}] and tissue == '{tissue}'",
                              column_names=["cell_type"]):
        for k, n in tbl.to_pandas().astype(str)["cell_type"].value_counts().items():
            acc[k] += int(n)
    return pd.Series(acc, dtype="int64").sort_values(ascending=False)


def build_group(human, dataset_ids, tissue: str, labels):
    """Donor-level summed raw counts for the cells of the given labels. Streams X in chunks.
    Returns (CSR donors x genes, donor table, gene table) or None if the query is empty."""
    import tiledbsoma as soma
    vf = value_filter(dataset_ids, tissue, labels)
    with human.axis_query(measurement_name="RNA", obs_query=soma.AxisQuery(value_filter=vf)) as q:
        obs = q.obs(column_names=OBS_COLS).concat().to_pandas()
        if len(obs) == 0:
            return None
        var = (q.var(column_names=["soma_joinid", "feature_id", "feature_name"]).concat().to_pandas()
               .sort_values("soma_joinid").reset_index(drop=True))
        donors, codes = group_cells(obs)
        stream = ((t["soma_dim_0"].to_numpy(), t["soma_dim_1"].to_numpy(), t["soma_data"].to_numpy())
                  for t in q.X("raw").tables())
        counts = accumulate_chunks(stream, obs["soma_joinid"].to_numpy(), codes, len(donors),
                                   var["soma_joinid"].to_numpy())
    return counts, donors, var[["feature_id", "feature_name"]]


def crosscheck_in_memory(census, dataset_ids, tissue: str, labels, counts, donors) -> float:
    """Max absolute difference between the streamed pseudobulk and one computed from an
    in-memory AnnData of the same cells (small queries only)."""
    import cellxgene_census
    from scipy import sparse
    ad = cellxgene_census.get_anndata(
        census, organism="Homo sapiens", X_name="raw",
        obs_value_filter=value_filter(dataset_ids, tissue, labels),
        obs_column_names=["soma_joinid", "donor_id"], var_column_names=["soma_joinid"])
    X = sparse.csr_matrix(ad.X)[:, np.argsort(ad.var["soma_joinid"].to_numpy(), kind="stable")]
    codes, uniq = pd.factorize(ad.obs["donor_id"], sort=True)
    if list(uniq) != list(donors["donor_id"]):
        raise ValueError("donor sets differ between streamed and in-memory data")
    G = sparse.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(uniq), len(codes)))
    return float(abs((G @ X) - counts).max())


def build_tissue(human, dataset_ids, tissue: str, class_labels: dict):
    """One streamed pass over all requested classes. Returns ({class: (CSR donors x genes, donor table)}, genes)
    or None if the query is empty. Only practical for small selections (see scripts/build_pseudobulk.py)."""
    import tiledbsoma as soma
    label_to_class = {l: c for c, ls in class_labels.items() for l in ls}
    vf = value_filter(dataset_ids, tissue, sorted(label_to_class))
    with human.axis_query(measurement_name="RNA", obs_query=soma.AxisQuery(value_filter=vf)) as q:
        obs = q.obs(column_names=OBS_COLS + ["cell_type"]).concat().to_pandas()
        if len(obs) == 0:
            return None
        classes = obs["cell_type"].astype(str).map(label_to_class)
        var = (q.var(column_names=["soma_joinid", "feature_id", "feature_name"]).concat().to_pandas()
               .sort_values("soma_joinid").reset_index(drop=True))
        table, codes = group_cells_by_class(obs, classes.to_numpy())
        stream = ((t["soma_dim_0"].to_numpy(), t["soma_dim_1"].to_numpy(), t["soma_data"].to_numpy())
                  for t in q.X("raw").tables())
        total = accumulate_chunks(stream, obs["soma_joinid"].to_numpy(), codes, len(table),
                                  var["soma_joinid"].to_numpy())
    out = {}
    for cls in class_labels:
        rows = (table["cell_class"].to_numpy() == cls).nonzero()[0]
        if len(rows):
            out[cls] = (sparse.csr_matrix(total[rows]), table.iloc[rows].drop(columns="cell_class").reset_index(drop=True))
    return out, var[["feature_id", "feature_name"]]
