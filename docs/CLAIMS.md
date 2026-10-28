# Claims and their status

| # | Claim | Evidence | Status |
|---|---|---|---|
| 1 | Simulators are deterministic; planted modules are connected; decoys are signal-free hubs | tests/test_sim.py | supported (code) |
| 2 | DE is calibrated under the null (p<0.05 rate in [0.035, 0.065]) and has high power with controlled FDR when effects are symmetric | tests/test_de.py | supported on simulated data only |
| 3 | Age confounding and cell-composition artefacts are removed by covariate adjustment | tests/test_de.py | supported on simulated data only |
| 4 | Median-of-ratios normalisation is less biased than library size when DE is one-sided; neither is calibrated in that case | tests/test_de.py, docs/NOTEBOOK.md | supported on simulated data only |
| 5 | Plain RWR is hub-biased; the degree-bin permutation z-score removes the correlation with degree (slightly over-corrects: Spearman about -0.19 in one run) | tests/test_network.py, results/ | supported on simulated data only |
| 6 | Corrected RWR improves module recovery when scores are degree-biased | docs/RESULTS_SIM.md | supported only under a toy generator that builds in the effect it corrects; **not evidence for real data** |
| 7 | Cross-disease conservative combination highlights shared over disease-specific modules | tests/test_crossdisease_validate.py | supported on simulated data only |
| 8 | AD and PD share cell-type-resolved network modules | none | untested; also not novel (docs/GAP.md) |
| 9 | Bias-corrected modules replicate better / enrich GWAS genes more on real data | none | untested |
| 10 | Novelty of the method | docs/GAP.md | not supported (prior art exists) |
