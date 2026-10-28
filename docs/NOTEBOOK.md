# Lab notebook

## 2026-10-28 (build session, sandbox without access to biological data hosts)
- Literature scan and partial data audit via web search (docs/GAP.md, docs/DATA_AUDIT.md). Result: the headline question is not novel; the project is repositioned as a pre-specified benchmark/replication.
- Implemented DE, network methods, cross-disease combination, validation statistics, simulators, benchmark. 20 tests pass.
- Finding during exploration: with plain library-size (CPM) normalisation, one-sided DE (10% of genes up by 0.7 log-units) inflated false positives in null genes (p<0.05 rate roughly 0.10-0.20 instead of 0.05). Median-of-ratios reduced it (about 0.09) but did not remove it; with symmetric effects both are calibrated (about 0.05). Kept median-of-ratios as default; documented as a limitation; tested.
- Benchmark (20 seeds x 4 settings, ~20 s): plain RWR picks decoy hubs 65-79% of the time, corrected 0.5-5.5%. Corrected RWR only clearly wins on precision where scores are degree-biased by construction. Twelve paired tests were run with no multiplicity correction; treat p-values as descriptive.
- Not done: any real data, STRING, GWAS validation sets, NetCore baseline, moderated DE, ILP formulation, push to GitHub (no credentials in the build environment).

## 2026-10-28 (real-data phase, PalmFace)
- Census metadata audits done (docs/DATA_AUDIT_CENSUS.md, docs/DATA_AUDIT_CELLTYPES.md); cohorts and roles drafted in docs/PROTOCOL.md (not frozen).
- Streaming self-test against the real Census passed (730 cells, 10 donors, max abs difference 0.0 vs in-memory AnnData, 74 s).
- Class-by-class streaming of PD-atlas prefrontal cortex from the Census failed: after roughly an hour not even the smallest class (OPC, 24,510 cells) had finished, then the job died with an S3 timeout (TileDB default request timeout is 3 s). Hypothesis, not measured: read amplification, because the cells of one class are scattered through the X array. Measured: accumulating 2e6 entries costs 53 ms (about 0.6 min for 1.25e9 entries), so compute is not the bottleneck.
- Decision: download each cohort's source .h5ad (resumable) and stream its CSR matrix sequentially with a progress bar (src/absolutus_v4/h5ad_stream.py). Validated on synthetic files: equality with ground truth for several block sizes, with and without `raw/X`, resume, non-integer guard. Real file sizes, download speed and the presence of integer counts in each real file are still unknown.
- Possible bonus to check: source files may carry author metadata (e.g. batch, PMI) that Census drops.
