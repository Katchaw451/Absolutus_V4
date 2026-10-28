# Data audit (Phase 0, partial)

**Nothing below was downloaded.** The build environment could only reach PyPI and GitHub; STRING,
NCBI/GEO and CELLxGENE were unreachable (connection refused by the sandbox proxy). Entries come from
web-search snippets of the providers' pages and must be verified from a machine with network access.

| Resource | Role | What the search showed | Status |
|---|---|---|---|
| SEA-AD cohort (Allen Institute) | AD discovery or replication | snRNA-seq/snATAC-seq processed data from 84 aged donors on the AWS Registry of Open Data (Allen community post) | verify download, size, metadata |
| SEA-AD MTG reference | cell-type reference | 166,868 nuclei, 5 younger donors; h5ad of several GB listed on brain-map.org | not a disease cohort |
| CELLxGENE Census | programmatic access to many public datasets | Python/R API; source h5ad can be downloaded per dataset | verify AD/PD datasets present |
| Cross-dementia h5ad | AD/FTD/PSP vs controls | third-party README describes 432,555 cells, 3 regions, 9 cell types, on a CELLxGENE collection and Synapse (syn52074156) | verify; Synapse may need login |
| PD single-nucleus data | PD arm | AD/PD papers report using public PD snRNA-seq (hundreds of thousands of nuclei); accession not identified here | **unverified: biggest risk** |
| STRING | PPI network | not checked | verify version, licence, confidence thresholds |
| Open Targets / GWAS Catalog | held-out genetic validation | not checked | verify access, define AD/PD gene sets BEFORE discovery |
| ROSMAP / AMP-AD | large AD cohorts | controlled access | exclude unless access is already approved |

Planning notes: several datasets are multi-GB h5ad files; use backed mode, pseudobulk early, and
record actual RAM limits. If no PD single-nucleus cohort is usable, fall back to bulk PD with
cell-composition adjustment and state it as a limitation.
