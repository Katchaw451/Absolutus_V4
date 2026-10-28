#!/usr/bin/env python3
"""Build donor x cell-class pseudobulk (summed raw counts), per cohort and tissue.

  python scripts/fetch_h5ad.py --cohort sea_ad_dlpfc            # 1. download the source file(s)
  python scripts/build_pseudobulk.py --cohort sea_ad_dlpfc      # 2. stream them (default --source h5ad)
  python scripts/build_pseudobulk.py --list
  python scripts/build_pseudobulk.py --smoke                    # Census streaming self-test (small data only)

Outputs (gitignored): data/pseudobulk/<cohort>/<tissue>/<class>.counts.npz + .donors.csv, genes.csv and
mapping_report.csv. Finished classes are skipped (resumable). Replication/exploratory cohorts stay sealed until
the git tag `protocol-v1` exists. `--source census` reads from the Census over the network: only suitable for
small selections (reading a few cell classes across many donors was found to be impractically slow).
"""
import argparse
import collections
import re
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
from scipy import sparse

from absolutus_v4 import celltypes
from absolutus_v4.cohorts import CENSUS_VERSION, COHORTS, PROTOCOL_TAG, SEALED_ROLES, SMOKE
from absolutus_v4.pseudobulk_stream import save_group

OUT = Path("data/pseudobulk")
H5AD = Path("data/h5ad")
SIZE_ORDER = ["endothelial", "opc", "microglia", "astrocyte", "inhibitory", "excitatory", "oligodendrocyte"]


def tag_exists(tag: str) -> bool:
    r = subprocess.run(["git", "tag", "--list", tag], capture_output=True, text=True)
    return r.stdout.strip() == tag


def slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_").lower()


def make_progress(tag):
    t0, state = time.time(), {"last": 0.0}

    def cb(frac):
        now = time.time()
        if now - state["last"] < 30 and frac < 1.0:
            return
        state["last"] = now
        el = (now - t0) / 60
        eta = el * (1 - frac) / frac if frac > 0 else float("nan")
        print(f"[{tag}] {frac * 100:5.1f}% | {el:6.1f} min elapsed | ETA ~{eta:6.1f} min", flush=True)
    return cb


def guard(key):
    cfg = COHORTS[key]
    if cfg["role"] in SEALED_ROLES and not tag_exists(PROTOCOL_TAG):
        sys.exit(f"'{key}' is a {cfg['role']} cohort and stays sealed until the git tag '{PROTOCOL_TAG}' exists. "
                 "Freeze docs/PROTOCOL.md first (commit, then: git tag protocol-v1).")
    return cfg


def report_and_todo(obs, key, tissue, classes, rebuild):
    """Print/save the label mapping; return {class: labels} for classes not yet built."""
    report = pd.DataFrame({"cell_type": obs.index, "n_cells": obs.values})
    st = [celltypes.status(l) for l in report["cell_type"]]
    report["status"], report["class_or_reason"] = [s[0] for s in st], [s[1] for s in st]
    d = OUT / key / slug(tissue)
    d.mkdir(parents=True, exist_ok=True)
    report.to_csv(d / "mapping_report.csv", index=False)
    print(f"\n== {key} / {tissue}\n{report.to_string(index=False)}", flush=True)
    todo = {}
    for cls in [c for c in SIZE_ORDER if c in classes]:
        labels = celltypes.labels_for(cls, obs.index)
        if not labels:
            print(f"[{cls}] no mapped labels in this file; skipped")
        elif Path(f"{d / cls}.counts.npz").exists() and not rebuild:
            print(f"[{cls}] exists; skipped")
        else:
            todo[cls] = labels
    return d, todo


def save_classes(d, dense, table, genes, only):
    gp = d / "genes.csv"
    if gp.exists():
        if not pd.read_csv(gp)["feature_id"].equals(genes["feature_id"]):
            sys.exit(f"gene table differs from {gp}; aborting")
    else:
        genes.to_csv(gp, index=False)
    for cls in only:
        rows = (table["cell_class"].to_numpy() == cls).nonzero()[0]
        if len(rows) == 0:
            continue
        donors = table.iloc[rows].drop(columns="cell_class").reset_index(drop=True)
        save_group(d / cls, sparse.csr_matrix(dense[rows]), donors)
        print(f"[{cls}] saved: {len(donors)} donors, {int(donors['n_cells'].sum()):,} cells", flush=True)


def build_h5ad(key, classes, rebuild):
    from absolutus_v4.h5ad_stream import observed_labels_h5ad, observed_tissues, pseudobulk_h5ad
    cfg = guard(key)
    built = set()
    for did in cfg["dataset_ids"]:
        path = H5AD / f"{did}.h5ad"
        if not path.exists():
            sys.exit(f"{path} not found. Run: python scripts/fetch_h5ad.py --cohort {key}")
        for tissue in cfg["tissues"]:
            obs = observed_labels_h5ad(path, tissue)
            if len(obs) == 0:
                print(f"{did}: no cells for tissue {tissue!r}; tissues in file:\n{observed_tissues(path).head(15)}")
                continue
            d, todo = report_and_todo(obs, key, tissue, classes, rebuild)
            dup = built & set(todo)
            if dup:
                sys.exit(f"class(es) {sorted(dup)} appear in several files of '{key}'; merging is not implemented")
            if not todo:
                continue
            print(f"reading {path.name} for classes {list(todo)} (one sequential pass)", flush=True)
            res = pseudobulk_h5ad(path, tissue, todo, progress=make_progress(f"{key}/{slug(tissue)}"))
            dense, table, genes, info = res
            print(f"matrix read from {info['counts_in']}; {info['cells_used']:,} of {info['n_obs']:,} cells used", flush=True)
            save_classes(d, dense, table, genes, todo)
            built |= set(todo)
    return 0


def build_census(key, classes, rebuild):
    from absolutus_v4.census_io import build_tissue, observed_labels, open_census
    cfg = guard(key)
    with open_census(CENSUS_VERSION) as census:
        human = census["census_data"]["homo_sapiens"]
        for tissue in cfg["tissues"]:
            obs = observed_labels(human, cfg["dataset_ids"], tissue)
            d, todo = report_and_todo(obs, key, tissue, classes, rebuild)
            if not todo:
                continue
            res = build_tissue(human, cfg["dataset_ids"], tissue, todo)
            if res is None:
                continue
            per_class, genes = res
            import numpy as np
            parts = [(c, m, dn) for c, (m, dn) in per_class.items()]
            table = pd.concat([dn.assign(cell_class=c) for c, _, dn in parts], ignore_index=True)
            dense = np.vstack([m.toarray() for _, m, _ in parts])
            save_classes(d, dense, table, genes, todo)
    return 0


def smoke():
    from absolutus_v4.census_io import build_tissue, crosscheck_in_memory, observed_labels, open_census
    ids = SMOKE["dataset_ids"]
    with open_census(CENSUS_VERSION) as census:
        human = census["census_data"]["homo_sapiens"]
        tissues = collections.Counter()
        for tbl in human.obs.read(value_filter=f"dataset_id in ['{ids[0]}']", column_names=["tissue"]):
            for k, n in tbl.to_pandas().astype(str)["tissue"].value_counts().items():
                tissues[k] += int(n)
        tissue = tissues.most_common(1)[0][0]
        labels = observed_labels(human, ids, tissue)
        labels = labels[labels >= 300]
        label = labels.idxmin()
        print(f"smoke: tissue={tissue!r} cell_type={label!r} ({int(labels[label]):,} cells)", flush=True)
        t0 = time.time()
        per_class, _ = build_tissue(human, ids, tissue, {"smoke": [label]})
        counts, donors = per_class["smoke"]
        print(f"streamed in {time.time() - t0:.1f}s: {counts.shape[0]} donors x {counts.shape[1]} genes, "
              f"nnz={counts.nnz:,}, total counts={counts.sum():,.0f}")
        assert (counts.sum(axis=1) > 0).all(), "a donor has zero total counts"
        diff = crosscheck_in_memory(census, ids, tissue, [label], counts, donors)
        print(f"max abs difference vs in-memory AnnData: {diff}")
        print("SMOKE OK" if diff < 1e-6 else "SMOKE FAILED")
        return 0 if diff < 1e-6 else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--cohort", choices=sorted(COHORTS))
    ap.add_argument("--source", choices=["h5ad", "census"], default="h5ad")
    ap.add_argument("--classes", nargs="+", default=list(celltypes.PRIMARY_CLASSES), choices=celltypes.CLASSES)
    ap.add_argument("--rebuild", action="store_true")
    a = ap.parse_args()
    if a.list:
        print(f"protocol tag '{PROTOCOL_TAG}' present: {tag_exists(PROTOCOL_TAG)}")
        for k, c in COHORTS.items():
            print(f"{k:16s} {c['disease']}  {c['role']:12s} {c['tissues']}")
        return 0
    if a.smoke:
        return smoke()
    if a.cohort:
        return (build_h5ad if a.source == "h5ad" else build_census)(a.cohort, a.classes, a.rebuild)
    ap.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
