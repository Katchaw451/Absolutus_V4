# Data audit (Phase 0, metadata only)

Source: CELLxGENE Census 2025-11-08, queried from PalmFace with `scripts/audit_census_obs.py` and
`scripts/audit_census_celltypes.py` (outputs: `docs/DATA_AUDIT_CENSUS.md`, `docs/DATA_AUDIT_CELLTYPES.md`).
Only metadata was read. Expression data of the sealed cohorts has not been opened.

## Cohorts (numbers from the audit outputs)

| Cohort | Disease labels (donors) | Region(s) | Assay | Planned role |
|---|---|---|---|---|
| SEA-AD, DLPFC | dementia 39, normal 44 (83 donors; 1,395,601 cells) | dorsolateral prefrontal cortex | 10x 3' v3 (1.31 M cells), 10x multiome (86 k) | AD discovery |
| PD multi-region atlas | Parkinson disease 73, normal 24 (97 donors overall; 90 donors / 503,414 cells in prefrontal cortex) | prefrontal cortex, M1, V1, medial globus pallidus, dorsal motor nucleus of vagus (no substantia nigra) | 10x 3' v3 | PD discovery (prefrontal cortex) |
| Layer-4 resilience (preprint) | AD 31, normal 15 (46 donors); 42 donors / 288,030 cells in Brodmann area 9 | BA9, BA7, BA17 | 10x v2, v3 and Drop-seq mixed | AD replication (BA9) |
| Kamath substantia nigra | PD 6, normal 8, Lewy body dementia 4 (18 donors; 387,483 cells in 7 cell-type datasets) | substantia nigra pars compacta | 10x 3' v3 | PD replication, exploratory (region differs; LBD excluded from the PD contrast) |
| Rexach cross-dementia | AD 10, Pick 9, PSP 11, normal 10 (40 donors; 432,555 cells) | insular cortex, BA4, V1 | 10x v2 (BA4, V1) and v3 (insular): assay is confounded with region | disease-specificity control, later |
| SEA-AD, MTG | 89 donors, 83 of them also in the DLPFC dataset | middle temporal gyrus | | not independent of SEA-AD DLPFC: not used for replication |
| Leng entorhinal / SFG | AD 7, normal 3 per region | | | too small; used only as a smoke-test dataset |

## Limitations that follow from the data

- Census carries no PMI, RIN or batch; available covariates are age (from `development_stage`, sometimes a lower bound such as "80 and over"), sex and assay.
- SEA-AD labels are "dementia" vs "normal", not neuropathologically defined AD.
- `donor_id` strings are only comparable within one study. Whether SEA-AD, the resilience dataset and the PD atlas share donors is not answerable from Census; it must be checked in the publications before freezing the protocol.
- The PD atlas has no substantia nigra, so cross-disease comparison is matched on prefrontal cortex; the Kamath replication is in a different region.
- "central nervous system macrophage" (PD atlas) is pooled with microglia (flagged in `celltypes.py`).

## Still to verify

- Open Targets: the GraphQL API refused PalmFace (HTTP 403); bulk files on the EBI FTP are the fallback (test pending).
- STRING version, licence and the Ensembl-to-symbol mapping.
- Per-region PD/normal donor split inside the prefrontal cortex of the PD atlas.
- Published donor provenance of each cohort (to confirm independence).
