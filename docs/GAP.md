# Related work and positioning

Searches run on 2026-10-28 (web search; abstracts and snippets only, full texts not reviewed).
This is a quick scan, not a systematic review; redo it properly before any submission.

## What already exists

Cross-disease, cell-type-resolved AD/PD analyses:
- Cell-type PLS-DA gene modules shared between AD and PD snRNA-seq: https://www.biorxiv.org/content/10.1101/2025.02.17.638647
- Comparative AD/PD single-cell cortical transcriptomics with pathway, network and cell-cell communication analysis: https://link.springer.com/article/10.1007/s12035-024-04419-7
- Large prefrontal-cortex snRNA-seq (AD, PD, both, controls) with gene regulatory network modules: https://scholarship.miami.edu/esploro/outputs/journalArticle/Transcriptomic-profiling-of-prefrontal-cortex-reveals/991031979505702976
- Cross-disease single-cell atlas of transcriptomic vulnerability: https://medrxiv.org/content/10.1101/2024.10.31.24316513

Degree / study bias in PPI networks and its correction:
- NetCore, network propagation using node coreness to correct degree bias: https://pmc.ncbi.nlm.nih.gov/articles/PMC7515737
- Study bias in PPI degree and a randomisation-based correction: https://pmc.ncbi.nlm.nih.gov/articles/PMC4523822/
- Study and technical bias as an explanation for PPI degree distributions (eLife): https://elifesciences.org/articles/99951.pdf

## Conclusion

1. "Do AD and PD share cell-type-specific modules?" is **not novel**.
2. "Correct PPI hub/degree bias in network propagation" is **not novel** (NetCore and others).
3. What remains defensible, and what V4 targets: a pre-registered, open, reproducible comparison of
   degree-bias-corrected versus uncorrected propagation (including NetCore as a baseline) applied to
   donor-level, composition-adjusted, cell-type-resolved effect sizes, with held-out GWAS enrichment
   and cross-cohort replication as the endpoints. That is a rigorous benchmark/replication
   contribution. Whether it adds anything beyond the above is an empirical question the project must answer.

## Not yet checked

QUBOLib/QOBLIB-style benchmark libraries (relevant only if the quantum-optimisation angle is revived),
DOMINO/heinz/DIAMOnD benchmark papers on neurodegeneration, and any 2026 preprints after the search date.
