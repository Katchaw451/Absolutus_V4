#!/usr/bin/env python3
"""Locate, size-check and download the source .h5ad of each cohort dataset (resumable).

  python scripts/fetch_h5ad.py --sizes                       # sizes of every cohort's files (no download)
  python scripts/fetch_h5ad.py --cohort pd_atlas_pfc         # download missing files to data/h5ad/
Downloads go to data/h5ad/<dataset_id>.h5ad.part and are renamed when complete; rerun the same command to resume.
"""
import argparse
import subprocess
import sys
import urllib.request
from pathlib import Path

from absolutus_v4.cohorts import CENSUS_VERSION, COHORTS
from absolutus_v4.fetch import https_url, human_gb

OUT = Path("data/h5ad")


def locate(did):
    import cellxgene_census
    loc = cellxgene_census.get_source_h5ad_uri(did, census_version=CENSUS_VERSION)
    get = (lambda k: loc[k]) if isinstance(loc, dict) else (lambda k: getattr(loc, k))
    try:
        region = get("s3_region")
    except Exception:
        region = "us-west-2"
    return get("uri"), region


def head_size(url):
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=30) as r:
            n = r.headers.get("Content-Length")
            return int(n) if n else None
    except Exception as e:  # network/permission problems are reported, not fatal
        print(f"  HEAD failed: {e}")
        return None


def datasets(cohort):
    keys = [cohort] if cohort else list(COHORTS)
    seen = {}
    for k in keys:
        for did in COHORTS[k]["dataset_ids"]:
            seen.setdefault(did, k)
    return seen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", choices=sorted(COHORTS))
    ap.add_argument("--sizes", action="store_true")
    ap.add_argument("--use-census-downloader", action="store_true",
                    help="fallback: cellxgene_census.download_source_h5ad (not resumable)")
    a = ap.parse_args()
    if not a.sizes and not a.cohort:
        ap.print_help()
        return 1
    OUT.mkdir(parents=True, exist_ok=True)
    total = 0
    for did, key in datasets(a.cohort).items():
        uri, region = locate(did)
        url = https_url(uri, region)
        size = head_size(url)
        total += size or 0
        final = OUT / f"{did}.h5ad"
        state = "downloaded" if final.exists() else "missing"
        print(f"{key:16s} {did}  {human_gb(size):>10s}  {state}\n  {url}", flush=True)
        if a.sizes or final.exists():
            continue
        if a.use_census_downloader:
            import cellxgene_census
            cellxgene_census.download_source_h5ad(did, to_path=str(final), census_version=CENSUS_VERSION)
            continue
        part = Path(f"{final}.part")
        r = subprocess.run(["curl", "-L", "-C", "-", "--fail", "--retry", "30", "--retry-delay", "5",
                            "--retry-all-errors", "--connect-timeout", "30",
                            "--speed-limit", "102400", "--speed-time", "30",
                            "--progress-bar", "-o", str(part), url])
        if r.returncode != 0:
            sys.exit(f"download failed (curl exit {r.returncode}); rerun the same command to resume")
        if size is not None and part.stat().st_size != size:
            sys.exit(f"size mismatch: {part.stat().st_size} != {size}; rerun to resume")
        part.rename(final)
        print(f"  saved {final}", flush=True)
    print(f"total listed: {human_gb(total)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
