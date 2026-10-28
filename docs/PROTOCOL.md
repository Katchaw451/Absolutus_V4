# Analysis protocol (DRAFT v0, NOT FROZEN)

Freeze by committing the final version and tagging `protocol-v1` after data QC but before any
look at the replication cohort or the held-out genetic gene sets. Changes after the tag go in
`docs/DEVIATIONS.md`. Items marked TBD depend on the data audit.

## Questions
Q1. After donor-level, covariate- and composition-adjusted testing, which PPI network modules are shared by AD and PD, and in which cell types?
Q2. Does degree-bias-corrected propagation (permutation z-score; NetCore as comparator) change the recovered modules, and do corrected modules replicate better and enrich held-out GWAS genes more than uncorrected ones?

## Data and splits (TBD after audit)
Discovery: one AD cohort and one PD cohort. Replication: a different AD cohort (and PD if available). Held-out genetics: AD and PD GWAS-supported gene sets from an open source, fixed and stored (hash) before discovery. None of the held-out data is used for tuning.

## Statistics
- Unit of replication: donor. Pseudobulk per donor x cell type, min cells per donor TBD.
- Model: log2 normalised expression ~ disease + age + sex + PMI + batch (+ estimated cell proportions in bulk). OLS baseline; moderated model (limma/pydeseq2) as cross-check. BH-FDR 5%.
- Composition shifts tested separately, never absorbed into DE.

## Network analysis
PPI: STRING, confidence threshold sweep (TBD values). Methods: RWR, RWR + degree-bin permutation z-score, DIAMOnD, greedy connected subgraph, NetCore. Module size k fixed per analysis in advance. Cross-disease score: signed minimum of AD and PD z.

## Endpoints
Primary: enrichment of held-out GWAS genes in corrected vs uncorrected modules (degree-matched permutation null). Secondary: replication sign concordance; donor-bootstrap stability; cell-type specificity; comparison with top-N DE genes.

## Nulls and controls
Donor-label permutation; degree-matched random modules; planted-module simulations on the real PPI topology (power).

## Stopping and reporting
Report every comparison run, including negative ones. No claim without a null model and a baseline.
