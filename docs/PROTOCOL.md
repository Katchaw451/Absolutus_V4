# Analysis protocol (v1 DRAFT for review; NOT FROZEN)

Freeze = commit the reviewed text, then `git tag protocol-v1`. Until the tag exists,
`scripts/build_pseudobulk.py` refuses to read the replication/exploratory cohorts. Changes after the tag go in
`docs/DEVIATIONS.md`. Values marked PROPOSED are my suggestions and need your decision (and, ideally, a supervisor's review).

## 1. Questions
Q1. After donor-level, covariate-adjusted testing, which PPI modules are shared by AD and PD within the same cell class?
Q2. Does degree-bias-corrected propagation change the recovered modules, and do corrected modules enrich held-out genetic evidence more than uncorrected ones?

## 2. Cohorts and roles (details and numbers in docs/DATA_AUDIT.md)
- AD discovery: SEA-AD DLPFC (`sea_ad_dlpfc`). PD discovery: PD multi-region atlas, prefrontal cortex (`pd_atlas_pfc`). Same region in both diseases.
- AD replication: Layer-4 resilience dataset, Brodmann area 9 (`ad_l4_ba9`). PD replication, exploratory: Kamath substantia nigra (`pd_kamath_snpc`); a different region, so a weaker test; Lewy body dementia donors are excluded from the PD contrast.
- Disease labels exactly as in Census: SEA-AD "dementia" vs "normal"; "Parkinson disease" vs "normal"; "Alzheimer disease" vs "normal".
- Precondition for freezing: published donor provenance shows the cohorts above do not share donors (if they do, the overlapping cohort is dropped from replication).

## 3. Units and cell classes
Unit of replication = donor. Pseudobulk = summed raw counts per donor x broad class. Primary classes: excitatory, inhibitory,
astrocyte, microglia, oligodendrocyte, OPC (mapping fixed in `src/absolutus_v4/celltypes.py`; endothelial too sparse for the primary analysis).
Inclusion (PROPOSED): a donor enters a class analysis with >= 30 cells in that class; a class is analysed if each arm has >= 8 donors.

## 4. Differential expression
log2 expression after median-of-ratios size factors; OLS per gene with disease + age + sex (+ assay where a cohort has more than one);
age taken from `development_stage`, lower-bound values ("80 and over") used as 80 with a flag; BH-FDR 5%.
Pre-specified sensitivity: NB GLM (pydeseq2); if conclusions differ, both are reported. Known limitation: no PMI/RIN/batch in Census.
Cell-class proportions per donor are tested separately and never absorbed into DE (snRNA proportions are not tissue composition).

## 5. Network analysis
PPI: STRING v12.0 (version pinned), combined score >= 700 primary; 400 and 900 as sensitivity (PROPOSED). Module size k = 50 primary, 30 sensitivity (PROPOSED).
Methods: RWR; RWR + degree-bin permutation z-score (200 permutations, 10 degree bins); DIAMOnD; greedy connected subgraph; NetCore as an external comparator (to be implemented before freezing).
Cross-disease score per gene within a class: signed minimum of the AD and PD z (conservative combination).
Known behaviour from simulation: the permutation z-score slightly over-corrects (Spearman with degree about -0.19).

## 6. Held-out genetic evidence (sealed)
Gene sets from Open Targets genetic association scores for Alzheimer disease and Parkinson disease: top 200 genes each (PROPOSED), downloaded, hashed
and committed as a manifest BEFORE discovery results exist. Sensitivity: remove the APOE locus genes for AD. Access status: GraphQL API refused (403); bulk files pending.
These sets are never used to tune any parameter.

## 7. Endpoints
Primary (single test): for each of the 12 analyses (6 classes x AD, PD), the difference in -log10 enrichment p-value (corrected minus uncorrected modules,
degree-matched permutation null, 1000 permutations); one-sided Wilcoxon signed-rank over the 12 differences, alpha 0.05. The 12 analyses share donors and are not
independent; a donor bootstrap CI of the mean difference is reported alongside.
Secondary (BH within family): replication sign concordance (binomial test), donor-bootstrap module stability, cell-class specificity, comparison with top-N DE genes, enrichment sensitivity analyses.

## 8. Nulls and controls
Negative control: donor-label permutation through the whole pipeline (200 repeats, PROPOSED), expecting no enrichment excess. Degree-matched random modules.
Power: planted-module simulations on the real STRING topology.

## 9. Sealing and reporting
Replication cohorts and genetic sets are not read before the tag (enforced for expression data by the build script). Every comparison run is reported, including negative ones.
No claim of novelty (docs/GAP.md) and no clinical claim.
