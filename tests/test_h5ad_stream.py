import numpy as np
import pandas as pd
import pytest
from scipy import sparse

pytest.importorskip("anndata")
pytest.importorskip("h5py")
import anndata as ad  # noqa: E402

from absolutus_v4.h5ad_stream import inspect, pseudobulk_h5ad  # noqa: E402

CLASS_LABELS = {"astrocyte": ["astrocyte"], "excitatory": ["glutamatergic neuron", "L5 neuron"]}


def _make(tmp_path, with_raw=True, seed=0, n_cells=500, n_genes=80):
    rng = np.random.default_rng(seed)
    counts = sparse.random(n_cells, n_genes, density=0.15, random_state=seed, format="csr")
    counts.data = np.round(counts.data * 30) + 1
    obs = pd.DataFrame({
        "tissue": rng.choice(["prefrontal cortex", "primary motor cortex"], n_cells, p=[0.7, 0.3]),
        "cell_type": rng.choice(["astrocyte", "glutamatergic neuron", "L5 neuron", "mural cell"], n_cells),
        "donor_id": rng.choice([f"d{i}" for i in range(6)], n_cells),
        "disease": rng.choice(["normal", "Parkinson disease"], n_cells),
        "sex": "female", "development_stage": "70 year-old human stage", "assay": "10x 3' v3"}).astype("category")
    var = pd.DataFrame({"feature_name": [f"G{i}" for i in range(n_genes)]}, index=[f"ENSG{i:05d}" for i in range(n_genes)])
    if with_raw:
        a = ad.AnnData(X=sparse.csr_matrix(np.log1p(counts.toarray())), obs=obs, var=var)
        a.raw = ad.AnnData(X=counts, obs=obs, var=var)
    else:
        a = ad.AnnData(X=counts, obs=obs, var=var)
    path = tmp_path / "t.h5ad"
    a.write_h5ad(path)
    return path, counts, obs, var


def _truth(counts, obs, tissue):
    out = {}
    for cls, labels in CLASS_LABELS.items():
        m = ((obs["tissue"].astype(str) == tissue) & obs["cell_type"].astype(str).isin(labels)).to_numpy()
        for d in sorted(set(obs["donor_id"][m].astype(str))):
            mm = m & (obs["donor_id"].astype(str) == d).to_numpy()
            out[(cls, d)] = np.asarray(counts[mm].sum(0)).ravel()
    return out


@pytest.mark.parametrize("with_raw", [True, False])
@pytest.mark.parametrize("block", [50, 500, 8_000_000])
def test_matches_ground_truth_with_and_without_raw_and_any_block_size(tmp_path, with_raw, block):
    path, counts, obs, var = _make(tmp_path, with_raw)
    info0 = inspect(path)
    assert info0["counts_in"] == ("raw/X" if with_raw else "X") and info0["encoding"] == "csr_matrix"
    seen = []
    dense, table, genes, info = pseudobulk_h5ad(path, "prefrontal cortex", CLASS_LABELS,
                                                block_entries=block, progress=seen.append)
    truth = _truth(counts, obs, "prefrontal cortex")
    assert len(table) == len(truth) and info["counts_in"] == info0["counts_in"]
    for i, r in table.iterrows():
        assert np.array_equal(dense[i], truth[(r["cell_class"], r["donor_id"])])
        assert r["n_cells"] > 0
    assert list(genes["feature_id"]) == list(var.index) and list(genes["feature_name"]) == list(var["feature_name"])
    assert seen and abs(seen[-1] - 1.0) < 1e-9 and all(b >= a - 1e-12 for a, b in zip(seen, seen[1:]))


def test_other_tissue_and_missing_tissue(tmp_path):
    path, counts, obs, _ = _make(tmp_path)
    dense, table, _, _ = pseudobulk_h5ad(path, "primary motor cortex", CLASS_LABELS)
    truth = _truth(counts, obs, "primary motor cortex")
    assert all(np.array_equal(dense[i], truth[(r.cell_class, r.donor_id)]) for i, r in table.iterrows())
    assert pseudobulk_h5ad(path, "no such tissue", CLASS_LABELS) is None


def test_normalised_matrix_is_rejected(tmp_path):
    rng = np.random.default_rng(0)
    X = sparse.csr_matrix(rng.random((30, 10)) + 0.3)
    obs = pd.DataFrame({"tissue": "t", "cell_type": "astrocyte", "donor_id": "d", "disease": "x", "sex": "f",
                        "development_stage": "s", "assay": "a"}, index=[f"c{i}" for i in range(30)])
    path = tmp_path / "n.h5ad"
    ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(10)])).write_h5ad(path)
    with pytest.raises(ValueError, match="non-integer"):
        pseudobulk_h5ad(path, "t", {"astrocyte": ["astrocyte"]})


def test_observed_labels_and_tissues(tmp_path):
    from absolutus_v4.h5ad_stream import observed_labels_h5ad, observed_tissues
    path, _, obs, _ = _make(tmp_path)
    lab = observed_labels_h5ad(path, "prefrontal cortex")
    assert int(lab.sum()) == int((obs["tissue"].astype(str) == "prefrontal cortex").sum())
    assert set(observed_tissues(path).index) == {"prefrontal cortex", "primary motor cortex"}
