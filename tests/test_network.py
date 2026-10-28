import networkx as nx
import numpy as np
from scipy.stats import spearmanr

from absolutus_v4.network import (adjacency, degree_bins, diamond, greedy_connected,
                                  permute_within_bins, prizes, recovery_metrics, rwr,
                                  rwr_corrected, topk, transition)
from absolutus_v4.sim import degrees, make_ppi, plant_module, simulate_z


def _toy(seed=0, n=1000):
    g = make_ppi(n, 3, seed)
    A = adjacency(g)
    return g, A, transition(A), degrees(g)


def test_transition_is_column_stochastic_and_rwr_conserves_mass():
    _, A, W, _ = _toy()
    assert np.allclose(np.asarray(W.sum(axis=0)).ravel(), 1.0)
    p0 = np.random.default_rng(0).random(A.shape[0])
    p = rwr(W, p0, 0.5)
    assert abs(p.sum() - 1.0) < 1e-8
    assert np.allclose(rwr(W, p0, 1.0), p0 / p0.sum())


def test_permutation_preserves_values_within_degree_bins():
    _, _, _, deg = _toy()
    rng = np.random.default_rng(0)
    w = rng.random(len(deg))
    bins = degree_bins(deg, 10)
    out = permute_within_bins(w, bins, 5, rng)
    for j in range(5):
        for b in np.unique(bins):
            idx = bins == b
            assert np.allclose(np.sort(out[idx, j]), np.sort(w[idx]))


def test_degree_correction_removes_hub_bias_of_plain_rwr():
    _, A, W, deg = _toy()
    rng = np.random.default_rng(0)
    w = rng.random(len(deg)) * (rng.random(len(deg)) < 0.05)  # random prizes: no true signal
    obs, z = rwr_corrected(W, w, deg, 0.5, 200, rng)
    assert spearmanr(obs, deg)[0] > 0.3          # plain RWR favours hubs
    assert abs(spearmanr(z, deg)[0]) < 0.3       # corrected score does not


def test_diamond_and_greedy_return_valid_modules():
    g, A, _, deg = _toy()
    rng = np.random.default_rng(1)
    mod, _ = plant_module(g, 30, 5, rng)
    z = simulate_z(1000, mod, 3.0, rng, deg, 0.0)
    seeds = topk(np.abs(z), 6)
    d = diamond(A, seeds, 30)
    assert len(d) == 30 and len(set(d.tolist())) == 30 and set(seeds.tolist()) <= set(d.tolist())
    gr = greedy_connected(A, z, 30)
    assert len(gr) == 30 and nx.is_connected(g.subgraph(gr.tolist()))


def test_strong_signal_is_recovered_far_above_chance():
    prec = {"top": [], "rwr_corrected": [], "greedy": []}
    for seed in range(5):
        rng = np.random.default_rng(seed)
        g, A, W, deg = _toy(seed)
        mod, dec = plant_module(g, 30, 10, rng)
        z = simulate_z(1000, mod, 3.0, rng, deg, 0.0)
        w = prizes(z)
        sel = {"top": topk(np.abs(z), 30),
               "rwr_corrected": topk(rwr_corrected(W, w, deg, 0.5, 100, rng)[1], 30),
               "greedy": greedy_connected(A, z, 30)}
        for k, v in sel.items():
            prec[k].append(recovery_metrics(v, mod, dec, deg)["precision"])
    chance = 30 / 1000
    for k, v in prec.items():
        assert np.mean(v) > 0.5, k
        assert np.mean(v) > 10 * chance


def test_recovery_metrics():
    m = recovery_metrics([1, 2, 3, 4], [1, 2, 9, 10], [3, 4, 5])
    assert m["precision"] == 0.5 and m["recall"] == 0.5 and abs(m["decoy_rate"] - 2 / 3) < 1e-12
