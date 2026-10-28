# Lab notebook

## 2026-10-28 (build session, sandbox without access to biological data hosts)
- Literature scan and partial data audit via web search (docs/GAP.md, docs/DATA_AUDIT.md). Result: the headline question is not novel; the project is repositioned as a pre-specified benchmark/replication.
- Implemented DE, network methods, cross-disease combination, validation statistics, simulators, benchmark. 20 tests pass.
- Finding during exploration: with plain library-size (CPM) normalisation, one-sided DE (10% of genes up by 0.7 log-units) inflated false positives in null genes (p<0.05 rate roughly 0.10-0.20 instead of 0.05). Median-of-ratios reduced it (about 0.09) but did not remove it; with symmetric effects both are calibrated (about 0.05). Kept median-of-ratios as default; documented as a limitation; tested.
- Benchmark (20 seeds x 4 settings, ~20 s): plain RWR picks decoy hubs 65-79% of the time, corrected 0.5-5.5%. Corrected RWR only clearly wins on precision where scores are degree-biased by construction. Twelve paired tests were run with no multiplicity correction; treat p-values as descriptive.
- Not done: any real data, STRING, GWAS validation sets, NetCore baseline, moderated DE, ILP formulation, push to GitHub (no credentials in the build environment).
