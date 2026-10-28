import numpy as np

from absolutus_v4.crossdisease import conservative_combine
from absolutus_v4.sim import degrees, make_ppi, plant_module, simulate_two_diseases
from absolutus_v4.validate import degree_matched_enrichment, replication_concordance


def test_conservative_combine_rules():
    out = conservative_combine([3.0, -2.0, 1.0, 0.5], [2.0, 2.5, -1.0, 4.0])
    assert np.allclose(out, [2.0, 0.0, 0.0, 0.5])  # discordant sign -> 0, else signed minimum


def test_shared_module_stands_out_after_conservative_combination():
    rng = np.random.default_rng(2)
    g = make_ppi(1000, 3, 2)
    s = simulate_two_diseases(g, 30, 30, 3.0, rng)
    c = np.abs(conservative_combine(s["z_ad"], s["z_pd"]))
    assert c[s["shared"]].mean() > 2 * c[s["ad_only"]].mean()
    assert c[s["shared"]].mean() > 2 * c[s["pd_only"]].mean()


def test_enrichment_null_is_calibrated():
    rng = np.random.default_rng(5)
    g = make_ppi(500, 3, 5)
    deg = degrees(g)
    pv = []
    for _ in range(100):
        sel = rng.choice(500, 30, replace=False)
        val = rng.choice(500, 40, replace=False)
        pv.append(degree_matched_enrichment(sel, val, deg, n_perm=200, rng=rng)[2])
    assert np.mean(np.array(pv) < 0.05) < 0.12  # nominal 0.05


def test_enrichment_detects_true_overlap():
    rng = np.random.default_rng(6)
    g = make_ppi(500, 3, 6)
    mod, _ = plant_module(g, 30, 0, rng)
    val = set(mod[:15]) | set(rng.choice(500, 10, replace=False).tolist())
    obs, null_mean, p = degree_matched_enrichment(mod, val, degrees(g), n_perm=500, rng=rng)
    assert obs >= 15 and obs > 5 * max(null_mean, 0.1) and p < 0.01


def test_replication_concordance():
    rng = np.random.default_rng(0)
    b = rng.normal(size=200)
    frac, p = replication_concordance(b, b, range(50))
    assert frac == 1.0 and p < 1e-10
    frac, p = replication_concordance(b, rng.normal(size=200), range(50))
    assert p > 0.01
