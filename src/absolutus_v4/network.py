"""Network module discovery with and without degree-bias correction."""
from __future__ import annotations

import networkx as nx
import numpy as np
from scipy import sparse
from scipy.stats import hypergeom


def adjacency(g: nx.Graph):
    n = g.number_of_nodes()
    return nx.to_scipy_sparse_array(g, nodelist=range(n), format="csr", dtype=float)


def transition(A):
    """Column-stochastic transition matrix W[i, j] = A[i, j] / deg(j)."""
    d = np.asarray(A.sum(axis=0)).ravel()
    d[d == 0] = 1.0
    c = A.tocoo()
    return sparse.csr_array((c.data / d[c.col], (c.row, c.col)), shape=A.shape)


def prizes(z, z0: float = 1.0):
    """Non-negative node prizes from z-scores: max(|z| - z0, 0)."""
    return np.maximum(np.abs(np.asarray(z, float)) - z0, 0.0)


def rwr(W, P0, restart: float = 0.5, tol: float = 1e-10, max_iter: int = 500):
    """Random walk with restart. P0: (n,) or (n, P); each column is normalised
    (a column with no mass becomes uniform)."""
    P0 = np.asarray(P0, float)
    single = P0.ndim == 1
    if single:
        P0 = P0[:, None]
    s = P0.sum(0, keepdims=True)
    P0 = np.where(s > 0, P0 / np.where(s > 0, s, 1.0), 1.0 / P0.shape[0])
    P = P0.copy()
    for _ in range(max_iter):
        Pn = (1 - restart) * (W @ P) + restart * P0
        done = np.abs(Pn - P).max() < tol
        P = Pn
        if done:
            break
    return P[:, 0] if single else P


def degree_bins(deg, n_bins: int = 10):
    qs = np.unique(np.quantile(deg, np.linspace(0, 1, n_bins + 1)[1:-1]))
    return np.digitize(deg, qs)


def permute_within_bins(w, bins, n_perm, rng):
    """n_perm permutations of w that only shuffle values among genes of the same degree bin."""
    n = len(w)
    out = np.empty((n, n_perm))
    groups = [np.where(bins == b)[0] for b in np.unique(bins)]
    for j in range(n_perm):
        for idx in groups:
            out[idx, j] = w[rng.permutation(idx)]
    return out


def rwr_corrected(W, w, deg, restart=0.5, n_perm=200, rng=None, n_bins=10):
    """RWR probabilities and their degree-bin-permutation z-scores.
    The z-score asks: is this gene reached more than expected for prizes of the
    same magnitude placed on genes of comparable degree?"""
    rng = rng or np.random.default_rng(0)
    obs = rwr(W, w, restart)
    null = rwr(W, permute_within_bins(w, degree_bins(deg, n_bins), n_perm, rng), restart)
    return obs, (obs - null.mean(1)) / (null.std(1) + 1e-12)


def topk(score, k):
    return np.argsort(-np.asarray(score), kind="stable")[:k]


def diamond(A, seeds, k_total):
    """DIAMOnD-style expansion (hypergeometric connectivity p-value). Returns k_total nodes."""
    n = A.shape[0]
    deg = np.asarray(A.sum(axis=1)).ravel()
    in_mod = np.zeros(n, bool)
    in_mod[list(seeds)] = True
    ks = np.asarray(A @ in_mod.astype(float)).ravel()
    module = [int(s) for s in seeds]
    while len(module) < k_total:
        cand = np.where((~in_mod) & (ks > 0))[0]
        if len(cand) == 0:
            break
        p = hypergeom.sf(ks[cand] - 1, n, len(module), deg[cand])
        best = int(cand[np.argmin(p)])
        module.append(best)
        in_mod[best] = True
        ks = ks + np.asarray(A[[best], :].toarray()).ravel()
    return np.array(module)


def greedy_connected(A, z, k, z0: float = 1.0, n_starts: int = 10):
    """Heinz-like greedy connected subgraph of size k maximising sum(|z| - z0)."""
    w = np.abs(np.asarray(z, float)) - z0
    ip, ix = A.indptr, A.indices
    best, best_val = None, -np.inf
    for s in np.argsort(-np.abs(z), kind="stable")[:n_starts]:
        S = [int(s)]
        inS = {int(s)}
        frontier = set(ix[ip[s]:ip[s + 1]].tolist())
        while len(S) < k and frontier:
            fl = np.array(sorted(frontier))
            v = int(fl[np.argmax(w[fl])])
            S.append(v)
            inS.add(v)
            frontier |= set(ix[ip[v]:ip[v + 1]].tolist())
            frontier -= inS
        val = float(w[S].sum())
        if val > best_val:
            best, best_val = S, val
    return np.array(best)


def recovery_metrics(selected, module, decoys, deg=None):
    sel, mod = set(int(v) for v in selected), set(int(v) for v in module)
    tp = len(sel & mod)
    out = {
        "precision": tp / max(1, len(sel)),
        "recall": tp / max(1, len(mod)),
        "decoy_rate": len(sel & set(int(v) for v in decoys)) / max(1, len(decoys)),
    }
    if deg is not None:
        out["mean_degree"] = float(np.mean(deg[list(sel)])) if sel else float("nan")
    return out
