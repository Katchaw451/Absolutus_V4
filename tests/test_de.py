import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from absolutus_v4.de import fit_de, log_cpm, log_norm, make_design, pseudobulk
from absolutus_v4.sim import simulate_bulk_mixture, simulate_pseudobulk


def test_null_calibration():
    rng = np.random.default_rng(11)
    c, meta, _, _ = simulate_pseudobulk(60, 4000, rng, de_frac=0.0, age_frac=0.0)
    r = fit_de(log_norm(c), make_design(meta, "disease", ["age_z", "sex"]))
    assert 0.035 < (r.p < 0.05).mean() < 0.065


def test_power_and_fdr_with_symmetric_effects():
    tpr, fdp = [], []
    for seed in (1, 2, 3, 4):
        rng = np.random.default_rng(seed)
        c, meta, is_de, _ = simulate_pseudobulk(60, 3000, rng, lfc=0.7, age_frac=0.0, symmetric=True)
        r = fit_de(log_norm(c), make_design(meta, "disease", ["age_z", "sex"]))
        hit = (r.fdr < 0.05).to_numpy()
        tpr.append(hit[is_de].mean())
        fdp.append((hit & ~is_de).sum() / max(1, hit.sum()))
    assert np.mean(tpr) > 0.8
    assert np.mean(fdp) < 0.10


def test_age_confounding_is_removed_by_adjustment():
    rng = np.random.default_rng(1)
    c, meta, _, _ = simulate_pseudobulk(80, 3000, rng, de_frac=0.0, confound=2.0, age_frac=0.2, age_slope=0.4)
    ex = log_norm(c)
    unadj = fit_de(ex, make_design(meta, "disease", []))
    adj = fit_de(ex, make_design(meta, "disease", ["age_z"]))
    assert (unadj.fdr < 0.05).sum() > 100
    assert (adj.fdr < 0.05).sum() <= 5
    assert (adj.p < 0.05).mean() < 0.08


def test_cell_composition_artefact_is_removed_by_adjustment():
    rng = np.random.default_rng(1)
    c, meta, _, _ = simulate_bulk_mixture(80, 3000, rng)
    ex = log_norm(c)
    unadj = fit_de(ex, make_design(meta, "disease", []))
    adj = fit_de(ex, make_design(meta, "disease", ["neuron_prop"]))
    assert (unadj.fdr < 0.05).sum() > 200  # spurious DE with no intrinsic change
    assert (adj.fdr < 0.05).sum() <= 5


def test_median_of_ratios_beats_library_size_when_de_is_one_sided():
    lib, med = [], []
    for seed in (1, 2, 3):
        rng = np.random.default_rng(seed)
        c, meta, is_de, _ = simulate_pseudobulk(60, 3000, rng, lfc=0.7, age_frac=0.0, symmetric=False)
        d = make_design(meta, "disease", ["age_z", "sex"])
        lib.append((fit_de(log_cpm(c), d).p[~is_de] < 0.05).mean())
        med.append((fit_de(log_norm(c), d).p[~is_de] < 0.05).mean())
    assert np.mean(lib) > np.mean(med) + 0.02


def test_pseudobulk_sums_per_donor_and_drops_small_groups():
    X = np.arange(24, dtype=float).reshape(8, 3)
    donors = np.array(["a", "a", "a", "b", "b", "b", "c", "c"])
    cts = np.array(["n"] * 8)
    for M in (X, sparse.csr_matrix(X)):
        out = pseudobulk(M, donors, cts, min_cells=3)
        counts, n = out["n"]
        assert list(counts.index) == ["a", "b"] and n["a"] == 3
        assert np.allclose(counts.loc["a"], X[:3].sum(0))


def test_fit_de_rejects_underdetermined_design():
    meta = pd.DataFrame({"disease": [0.0, 1.0], "age_z": [0.1, 0.2]})
    with pytest.raises(ValueError):
        fit_de(np.zeros((2, 5)), make_design(meta, "disease", ["age_z"]))
