#!/usr/bin/env python3
"""Build donor x cell-class pseudobulk (summed raw counts) from CELLxGENE Census by streaming X.

  python scripts/build_pseudobulk.py --list
  python scripts/build_pseudobulk.py --smoke                       # validates the streaming code (no analysis data)
  python scripts/build_pseudobulk.py --cohort sea_ad_dlpfc [--classes astrocyte microglia]

Outputs (gitignored): data/pseudobulk/<cohort>/<tissue>/<class>.counts.npz + .donors.csv,
data/pseudobulk/genes.csv and a mapping_report.csv per tissue. Existing outputs are skipped (resumable).
Replication/exploratory cohorts stay sealed until the git tag `protocol-v1` exists.
"""
import argparse
import collections
import re
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

from absolutus_v4 import celltypes
from absolutus_v4.cohorts import CENSUS_VERSION, COHORTS, PROTOCOL_TAG, SEALED_ROLES, SMOKE
from absolutus_v4.pseudobulk_stream import save_group

OUT = Path("data/pseudobulk")
SIZE_ORDER = ["endothelial", "opc", "microglia", "astrocyte", "inhibitory", "excitatory", "oligodendrocyte"]


def tag_exists(tag: str) -> bool:
    r = subprocess.run(["git", "tag", "--list", tag], capture_output=True, text=True)
    return r.stdout.strip() == tag


def slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_").lower()


def smoke():
    import cellxgene_census
    from absolutus_v4.census_io import build_group, crosscheck_in_memory, observed_labels
    ids = SMOKE["dataset_ids"]
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
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
        counts, donors, var = build_group(human, ids, tissue, [label])
        print(f"streamed in {time.time() - t0:.1f}s: {counts.shape[0]} donors x {counts.shape[1]} genes, "
              f"nnz={counts.nnz:,}, total counts={counts.sum():,.0f}")
        assert (counts.sum(axis=1) > 0).all(), "a donor has zero total counts"
        diff = crosscheck_in_memory(census, ids, tissue, [label], counts, donors)
        print(f"max abs difference vs in-memory AnnData: {diff}")
        print("SMOKE OK" if diff < 1e-6 else "SMOKE FAILED")
        return 0 if diff < 1e-6 else 1


def build(key: str, classes, rebuild: bool):
    cfg = COHORTS[key]
    if cfg["role"] in SEALED_ROLES and not tag_exists(PROTOCOL_TAG):
        sys.exit(f"'{key}' is a {cfg['role']} cohort and stays sealed until the git tag '{PROTOCOL_TAG}' exists. "
                 "Freeze docs/PROTOCOL.md first (commit, then: git tag protocol-v1).")
    import cellxgene_census
    from absolutus_v4.census_io import build_group, observed_labels
    OUT.mkdir(parents=True, exist_ok=True)
    genes_path = OUT / "genes.csv"
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        human = census["census_data"]["homo_sapiens"]
        for tissue in cfg["tissues"]:
            obs = observed_labels(human, cfg["dataset_ids"], tissue)
            report = pd.DataFrame({"cell_type": obs.index, "n_cells": obs.values})
            st = [celltypes.status(l) for l in report["cell_type"]]
            report["status"] = [s[0] for s in st]
            report["class_or_reason"] = [s[1] for s in st]
            d = OUT / key / slug(tissue)
            d.mkdir(parents=True, exist_ok=True)
            report.to_csv(d / "mapping_report.csv", index=False)
            print(f"\n== {key} / {tissue}\n{report.to_string(index=False)}", flush=True)
            for cls in [c for c in SIZE_ORDER if c in classes]:
                labels = celltypes.labels_for(cls, obs.index)
                stem = d / cls
                if not labels:
                    print(f"[{cls}] no mapped labels in this cohort; skipped")
                    continue
                if Path(f"{stem}.counts.npz").exists() and not rebuild:
                    print(f"[{cls}] exists; skipped")
                    continue
                t0 = time.time()
                res = build_group(human, cfg["dataset_ids"], tissue, labels)
                if res is None:
                    print(f"[{cls}] empty query; skipped")
                    continue
                counts, donors, var = res
                if genes_path.exists():
                    old = pd.read_csv(genes_path)
                    if not old["feature_id"].equals(var["feature_id"]):
                        sys.exit("gene table differs from the one already on disk; aborting")
                else:
                    var.to_csv(genes_path, index=False)
                save_group(stem, counts, donors)
                print(f"[{cls}] {len(donors)} donors, {int(donors['n_cells'].sum()):,} cells, nnz={counts.nnz:,}, "
                      f"{(time.time() - t0) / 60:.1f} min", flush=True)
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--cohort", choices=sorted(COHORTS))
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
        return build(a.cohort, a.classes, a.rebuild)
    ap.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
