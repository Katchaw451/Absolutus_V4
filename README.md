# Absolutus V4

Study-bias-aware network module discovery for neurodegeneration transcriptomics
(Alzheimer's / Parkinson's). **Status: validation stage (v0.1.0).** The statistical and network
code is implemented and tested against simulated data with known ground truth. **No real-data
results exist yet and this repository makes no claim about Alzheimer's or Parkinson's biology.**

## What is here

| Component | Module | What it does |
|---|---|---|
| Donor-level DE | `de.py` | Pseudobulk aggregation, median-of-ratios normalisation, OLS with covariates (age, sex, estimated cell proportions), BH-FDR |
| Network methods | `network.py` | Random walk with restart (RWR), degree-bin permutation z-score (hub-bias correction), DIAMOnD-style expansion, greedy connected subgraph |
| Cross-disease | `crossdisease.py` | Conservative (signed-minimum) combination of two diseases' z-scores |
| Validation | `validate.py` | Degree-matched enrichment of held-out gene sets, cross-cohort replication concordance |
| Simulators | `sim.py` | Scale-free PPI with planted modules and signal-free decoy hubs; pseudobulk counts with age confounding; bulk cell-composition mixtures |
| Benchmark | `bench.py` | Planted-module recovery of every method, with paired tests |

20 pytest tests cover calibration of the DE model under the null, power/FDR, removal of age and
cell-composition artefacts by adjustment, the hub bias of plain RWR and its removal, and
null calibration of the enrichment test.

## Quickstart

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install -e ".[dev]"
pytest -q                 # 20 tests
make bench                # regenerates results/ and docs/RESULTS_SIM.md (~20 s)
```

## What the simulations show (and do not show)

On a synthetic scale-free network with a planted module and 10 signal-free decoy hubs
(`docs/RESULTS_SIM.md`, 20 seeds per setting):

- Plain RWR selects decoy hubs 65-79% of the time; the degree-bin permutation z-score
  selects them 0.5-5.5% of the time.
- Module recovery is hard for every method at these effect sizes (precision <= ~0.37).
- The corrected RWR clearly beats the alternatives only in the settings where the scores
  themselves are degree-biased (`degree_bias=0.5`). That bias is a toy construction chosen by
  the authors, so this is a check of the mechanism, **not evidence that the correction helps on real data**.

## Where this sits in the literature

Cross-disease, cell-type-resolved AD/PD analyses and degree-bias-aware network propagation
already exist (see `docs/GAP.md`). The defensible contribution of this project is a
pre-specified, fully reproducible head-to-head evaluation with donor-level statistics and
held-out genetic validation, not a new method. Novelty is not claimed.

## Roadmap

1. Real-data ingestion and QC (needs network access to CELLxGENE / GEO / STRING; see `docs/DATA_AUDIT.md`, `docs/PHASE2_PROMPT.md`)
2. Freeze `docs/PROTOCOL.md` (git tag `protocol-v1`) before touching replication data or held-out genetics
3. NetCore baseline and an ILP/greedy cardinality-constrained formulation
4. Moderated DE (limma/pydeseq2 cross-check), analysis, validation, report

## Disclosure

Initial scaffolding, code and documentation were written with Claude (Anthropic); every number in
`docs/RESULTS_SIM.md` comes from `make bench`.

## License

MIT
