"""Cohort registry and sealing rules. Roles follow docs/PROTOCOL.md."""
from __future__ import annotations

COHORTS = {
    "sea_ad_dlpfc": dict(
        disease="AD", role="discovery",
        dataset_ids=["6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3"],
        tissues=["dorsolateral prefrontal cortex"]),
    "pd_atlas_pfc": dict(
        disease="PD", role="discovery",
        dataset_ids=["d3cb449b-c2b1-4b50-a7f1-21203535fe61"],
        tissues=["prefrontal cortex"]),
    "ad_l4_ba9": dict(
        disease="AD", role="replication",
        dataset_ids=["cff99df2-4904-44f7-9173-ff837f95606e"],
        tissues=["Brodmann (1909) area 9"]),
    "pd_kamath_snpc": dict(
        disease="PD", role="exploratory",
        dataset_ids=["2856d06c-0ff9-4e01-bfc9-202b74d0b60f", "251b1a7e-d050-4486-8d50-4c2619eb0f46",
                     "07760522-707a-4a1c-8891-dbd1226d6b27", "f8d8b443-bca6-4c3c-9042-669dfb7f8030",
                     "5829c7ba-697f-418e-8b98-d605b192dc48", "4dd1cd23-fc4d-4fd1-9709-602540f3ca6f",
                     "f5be4b96-f5a3-4c3d-84ac-6f69daf744d5"],
        tissues=["substantia nigra pars compacta"]),
}
# Used only to validate the streaming code against in-memory AnnData; not part of any analysis.
SMOKE = dict(dataset_ids=["2727d83a-0af0-443a-bff8-58dc7028289a"])
SEALED_ROLES = {"replication", "exploratory"}
PROTOCOL_TAG = "protocol-v1"
CENSUS_VERSION = "2025-11-08"
