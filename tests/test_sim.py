import networkx as nx
import numpy as np

from absolutus_v4.sim import degrees, make_ppi, plant_module, simulate_z


def test_graph_and_module_are_deterministic():
    out = []
    for _ in range(2):
        rng = np.random.default_rng(7)
        g = make_ppi(500, 3, 7)
        mod, dec = plant_module(g, 20, 5, rng)
        z = simulate_z(500, mod, 2.0, rng, degrees(g), 0.0)
        out.append((sorted(g.edges()), mod, dec, z))
    assert out[0][0] == out[1][0] and out[0][1] == out[1][1] and out[0][2] == out[1][2]
    assert np.array_equal(out[0][3], out[1][3])


def test_planted_module_is_connected_and_decoys_are_signal_free_hubs():
    rng = np.random.default_rng(3)
    g = make_ppi(800, 3, 3)
    mod, dec = plant_module(g, 25, 8, rng)
    deg = degrees(g)
    assert len(mod) == 25 and nx.is_connected(g.subgraph(mod))
    assert not set(mod) & set(dec)
    assert min(deg[dec]) >= np.percentile(deg, 95)  # decoys are genuine hubs
