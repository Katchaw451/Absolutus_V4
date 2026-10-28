"""Simulators with known ground truth, used to validate the pipeline.

Everything here is synthetic. Results obtained on these data validate the code;
they say nothing about Alzheimer's or Parkinson's disease.
"""
from __future__ import annotations

import networkx as nx
import numpy as np
import pandas as pd


def make_ppi(n_genes: int, m: int, seed: int) -> nx.Graph:
    """Scale-free toy PPI (Barabasi-Albert). Nodes are 0..n_genes-1."""
    return nx.barabasi_albert_graph(n_genes, m, seed=seed)


def degrees(g: nx.Graph) -> np.ndarray:
    return np.array([g.degree(i) for i in range(g.number_of_nodes())])


def plant_module(g, k, n_decoys, rng, exclude=()):
    """Plant a connected module of size k (grown from a non-hub start, preferring
    low-degree neighbours) and pick `n_decoys` highest-degree nodes outside it.
    Decoy hubs carry NO signal: they test whether a method confuses
    connectivity with evidence."""
    deg = dict(g.degree())
    excl = set(exclude)
    nodes = [v for v in g.nodes() if v not in excl]
    d = np.array([deg[v] for v in nodes])
    cutoff = np.percentile(d, 75)
    starts = [v for v, dv in zip(nodes, d) if dv <= cutoff]
    start = int(rng.choice(starts))
    module, mset = [start], {start}
    frontier = {v for v in g.neighbors(start) if v not in excl}
    while len(module) < k and frontier:
        fl = sorted(frontier)
        w = np.array([1.0 / deg[v] for v in fl])
        v = int(rng.choice(fl, p=w / w.sum()))
        module.append(v)
        mset.add(v)
        frontier |= {u for u in g.neighbors(v) if u not in excl}
        frontier -= mset
    ranked = sorted((v for v in g.nodes() if v not in mset and v not in excl), key=lambda v: -deg[v])
    return sorted(module), ranked[:n_decoys]


def simulate_z(n, module, effect, rng, degree=None, degree_bias=0.0):
    """Gene-level z-scores: N(0,1) noise, +/-effect on module genes (random sign).
    degree_bias > 0 inflates the noise sd of high-degree genes (toy study-bias)."""
    sd = np.ones(n)
    if degree is not None and degree_bias > 0:
        ld = np.log(np.asarray(degree, float))
        s = np.maximum((ld - ld.mean()) / ld.std(), 0.0)
        sd = 1.0 + degree_bias * s
    z = rng.standard_normal(n) * sd
    z[module] += effect * rng.choice([-1.0, 1.0], size=len(module))
    return z


def simulate_two_diseases(g, k_shared, k_specific, effect, rng):
    """Two diseases sharing one planted module plus one disease-specific module each."""
    n = g.number_of_nodes()
    shared, _ = plant_module(g, k_shared, 0, rng)
    ad_only, _ = plant_module(g, k_specific, 0, rng, exclude=shared)
    pd_only, _ = plant_module(g, k_specific, 0, rng, exclude=set(shared) | set(ad_only))
    signs = rng.choice([-1.0, 1.0], size=len(shared))
    z_ad, z_pd = rng.standard_normal(n), rng.standard_normal(n)
    z_ad[shared] += effect * signs
    z_pd[shared] += effect * signs
    z_ad[ad_only] += effect * rng.choice([-1.0, 1.0], size=len(ad_only))
    z_pd[pd_only] += effect * rng.choice([-1.0, 1.0], size=len(pd_only))
    return {"z_ad": z_ad, "z_pd": z_pd, "shared": shared, "ad_only": ad_only, "pd_only": pd_only}


def _nb(mu, phi, rng):
    return rng.negative_binomial(1.0 / phi, 1.0 / (1.0 + mu * phi))


def simulate_pseudobulk(n_donors, n_genes, rng, de_frac=0.1, lfc=0.7, confound=0.0,
                        age_frac=0.2, age_slope=0.3, phi=0.1, symmetric=False):
    """Donor-level counts (NB). Disease shifts `de_frac` of genes by `lfc` (natural log).
    Age shifts `age_frac` other genes; `confound` makes disease probability depend on age.
    symmetric=True gives DE genes random signs (up/down); False makes them all up, which
    biases any global normalisation (a known failure mode worth testing)."""
    age_z = rng.standard_normal(n_donors)
    p = 1.0 / (1.0 + np.exp(-confound * age_z))
    disease = (rng.random(n_donors) < p).astype(float)
    sex = rng.integers(0, 2, n_donors).astype(float)
    perm = rng.permutation(n_genes)
    n_de, n_age = int(de_frac * n_genes), int(age_frac * n_genes)
    is_de = np.zeros(n_genes, bool)
    is_de[perm[:n_de]] = True
    is_age = np.zeros(n_genes, bool)
    is_age[perm[n_de:n_de + n_age]] = True
    base = rng.normal(3.0, 1.5, n_genes)
    lib = rng.normal(0, 0.3, n_donors)
    sgn = rng.choice([-1.0, 1.0], size=n_genes) if symmetric else np.ones(n_genes)
    logmu = (base[None, :] + lib[:, None] + lfc * np.outer(disease, is_de * sgn)
             + age_slope * np.outer(age_z, is_age))
    counts = _nb(np.exp(logmu), phi, rng)
    meta = pd.DataFrame({"disease": disease, "age_z": age_z, "sex": sex})
    return counts, meta, is_de, is_age


def simulate_bulk_mixture(n_donors, n_genes, rng, neuron_shift=-0.15, marker_frac=0.2,
                          fold=8.0, phi=0.1, deconv_noise=0.03):
    """Bulk tissue = mixture of neurons / two glial types. Disease lowers the neuron
    fraction but NO gene is intrinsically differentially expressed: any DE call on
    marker genes without composition adjustment is a composition artefact."""
    disease = (np.arange(n_donors) % 2).astype(float)
    rng.shuffle(disease)
    f_n = np.clip(0.5 + neuron_shift * disease + rng.normal(0, 0.05, n_donors), 0.05, 0.95)
    rest = 1.0 - f_n
    a = rng.beta(5, 5, n_donors)
    props = np.column_stack([f_n, rest * a, rest * (1 - a)])
    base = rng.normal(3.0, 1.0, n_genes)
    prof = np.tile(np.exp(base), (3, 1))
    perm = rng.permutation(n_genes)
    m = int(marker_frac * n_genes)
    neuron_markers = np.zeros(n_genes, bool)
    neuron_markers[perm[:m]] = True
    glia_markers = np.zeros(n_genes, bool)
    glia_markers[perm[m:2 * m]] = True
    prof[0, neuron_markers] *= fold
    prof[1:, glia_markers] *= fold
    mu = (props @ prof) * np.exp(rng.normal(0, 0.3, n_donors))[:, None]
    counts = _nb(mu, phi, rng)
    meta = pd.DataFrame({"disease": disease, "neuron_prop": f_n + rng.normal(0, deconv_noise, n_donors)})
    return counts, meta, neuron_markers, glia_markers
