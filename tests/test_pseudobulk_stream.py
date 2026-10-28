import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from absolutus_v4.pseudobulk_stream import (accumulate_chunks, group_cells, load_group, save_group,
                                            value_filter)


def _synthetic(seed=0, n_cells=400, n_genes=60, n_donors=7):
    rng = np.random.default_rng(seed)
    X = sparse.random(n_cells, n_genes, density=0.2, random_state=seed, format="coo")
    X.data = np.round(X.data * 20) + 1
    donor = rng.integers(0, n_donors, n_cells)
    obs_ids = rng.choice(10_000_000, n_cells, replace=False)       # non-contiguous global ids
    var_ids = rng.choice(100_000, n_genes, replace=False)           # unsorted gene ids
    return X.tocsr(), donor, obs_ids, var_ids


def _dense_truth(X, donor, n_donors):
    G = sparse.csr_matrix((np.ones(len(donor)), (donor, np.arange(len(donor)))), shape=(n_donors, len(donor)))
    return (G @ X).toarray()


def _chunks(X, obs_ids, var_ids, size, rng, split_duplicates=False):
    c = X.tocoo()
    d0, d1, v = obs_ids[c.row], var_ids[c.col], c.data
    idx = rng.permutation(len(v))
    d0, d1, v = d0[idx], d1[idx], v[idx]
    if split_duplicates:  # the same (cell, gene) arriving as two partial values must be summed
        d0, d1, v = np.r_[d0, d0], np.r_[d1, d1], np.r_[v / 2, v / 2]
    for s in range(0, len(v), size):
        yield d0[s:s + size], d1[s:s + size], v[s:s + size]


def test_matches_dense_ground_truth_for_any_chunking_and_order():
    X, donor, obs_ids, var_ids = _synthetic()
    truth = _dense_truth(X, donor, 7)
    for size in (1, 37, 500, 10_000):
        out = accumulate_chunks(_chunks(X, obs_ids, var_ids, size, np.random.default_rng(1)),
                                obs_ids, donor, 7, var_ids)
        assert np.array_equal(out.toarray(), truth), size


def test_duplicate_entries_are_summed():
    X, donor, obs_ids, var_ids = _synthetic(seed=3)
    out = accumulate_chunks(_chunks(X, obs_ids, var_ids, 100, np.random.default_rng(2), split_duplicates=True),
                            obs_ids, donor, 7, var_ids)
    assert np.allclose(out.toarray(), _dense_truth(X, donor, 7))


def test_column_i_corresponds_to_var_joinids_i():
    X, donor, obs_ids, var_ids = _synthetic(seed=4)
    out = accumulate_chunks(_chunks(X, obs_ids, var_ids, 200, np.random.default_rng(3)), obs_ids, donor, 7, var_ids)
    j = 5
    expect = _dense_truth(X[:, [j]], donor, 7).ravel()
    assert np.array_equal(out.toarray()[:, j], expect)


def test_foreign_cells_or_genes_raise():
    X, donor, obs_ids, var_ids = _synthetic(seed=5)
    bad_cell = [(np.array([max(obs_ids) + 1]), np.array([var_ids[0]]), np.array([1.0]))]
    bad_gene = [(np.array([obs_ids[0]]), np.array([max(var_ids) + 1]), np.array([1.0]))]
    with pytest.raises(ValueError):
        accumulate_chunks(bad_cell, obs_ids, donor, 7, var_ids)
    with pytest.raises(ValueError):
        accumulate_chunks(bad_gene, obs_ids, donor, 7, var_ids)
    with pytest.raises(ValueError):
        accumulate_chunks([], [], [], 0, var_ids)


def test_group_cells_counts_and_visible_metadata_conflicts():
    obs = pd.DataFrame({"donor_id": ["b", "a", "a", "b", "c"], "disease": ["AD", "normal", "normal", "AD", "AD"],
                        "sex": ["F", "M", "M", "F", "F"],
                        "development_stage": ["80 year-old human stage"] * 5,
                        "assay": ["10x 3' v3", "10x 3' v3", "10x 3' v2", "10x 3' v3", "10x 3' v3"]})
    donors, codes = group_cells(obs)
    assert list(donors["donor_id"]) == ["a", "b", "c"] and list(donors["n_cells"]) == [2, 2, 1]
    assert list(codes) == [1, 0, 0, 1, 2]
    assert donors.loc[0, "assay"] == "10x 3' v2; 10x 3' v3"   # within-donor mixture is visible, not hidden


def test_value_filter_and_unsafe_literals():
    vf = value_filter(["id1", "id2"], "Brodmann (1909) area 9", ["astrocyte", "L5/6 NP"])
    assert vf == ("dataset_id in ['id1', 'id2'] and tissue == 'Brodmann (1909) area 9' "
                  "and cell_type in ['astrocyte', 'L5/6 NP']")
    with pytest.raises(ValueError):
        value_filter(["id1"], "o'brien", ["x"])


def test_save_and_load_roundtrip(tmp_path):
    X, donor, obs_ids, var_ids = _synthetic(seed=6)
    out = accumulate_chunks(_chunks(X, obs_ids, var_ids, 300, np.random.default_rng(4)), obs_ids, donor, 7, var_ids)
    donors = pd.DataFrame({"donor_id": list("abcdefg"), "n_cells": np.bincount(donor, minlength=7)})
    save_group(tmp_path / "x" / "astrocyte", out, donors)
    out2, donors2 = load_group(tmp_path / "x" / "astrocyte")
    assert np.array_equal(out.toarray(), out2.toarray()) and donors.equals(donors2)
