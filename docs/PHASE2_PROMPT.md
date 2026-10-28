# Prompt to continue on a machine with network access (e.g. katchaw451@PalmFace)

Paste into Claude Code (or another agent) from the repository root.

```
You are continuing Absolutus_V4 in this repository. Read README.md, docs/GAP.md, docs/DATA_AUDIT.md,
docs/PROTOCOL.md, docs/CLAIMS.md and docs/NOTEBOOK.md first. Work only through bash; never print or commit credentials.

1. Verify network access (curl -I to CELLxGENE, GEO, STRING, Open Targets) and resources (free -h, nproc, df -h).
2. Complete docs/DATA_AUDIT.md by actually testing every download (headers + small sample). Pick >=2 independent AD
   cohorts, >=1 PD cohort (or document the bulk fallback), STRING, and an open genetics source. Exclude controlled-access data.
3. Write download scripts (data/ is gitignored; commit scripts and a SHA256 manifest), QC, and donor-level metadata tables, with tests.
4. Fill the TBD items in docs/PROTOCOL.md, ask me to review it, then commit and tag protocol-v1. Do not open the replication
   cohort or held-out genetics before that tag.
5. Add a NetCore baseline and a moderated DE cross-check (pydeseq2 or limma) with tests on the existing simulators.
6. Run the analysis exactly as the protocol states, report every comparison including negative ones, update docs/CLAIMS.md
   and docs/NOTEBOOK.md, regenerate everything with make targets, commit and push after each phase.
Never claim novelty or biological conclusions beyond what the protocol endpoints support.
```
