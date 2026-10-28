import math

from absolutus_v4 import celltypes
from absolutus_v4.cohorts import COHORTS, SEALED_ROLES
from absolutus_v4.covariates import parse_age

SEA_AD_LABELS = [
    "L2/3-6 intratelencephalic projecting glutamatergic neuron", "oligodendrocyte",
    "pvalb GABAergic cortical interneuron", "VIP GABAergic cortical interneuron",
    "astrocyte of the cerebral cortex", "lamp5 GABAergic cortical interneuron",
    "sst GABAergic cortical interneuron", "microglial cell",
    "corticothalamic-projecting glutamatergic cortical neuron", "oligodendrocyte precursor cell",
    "sncg GABAergic cortical interneuron", "near-projecting glutamatergic cortical neuron",
    "L6b glutamatergic cortical neuron", "chandelier pvalb GABAergic cortical interneuron",
    "caudal ganglionic eminence derived cortical interneuron", "vascular leptomeningeal cell",
    "L5 extratelencephalic projecting glutamatergic cortical neuron", "cerebral cortex endothelial cell",
]


def test_every_observed_sea_ad_label_is_mapped_or_explicitly_excluded():
    for l in SEA_AD_LABELS:
        assert celltypes.status(l)[0] in ("class", "excluded"), l


def test_other_cohort_labels_from_the_audit():
    expect = {"glutamatergic neuron": "excitatory", "GABAergic neuron": "inhibitory",
              "inhibitory interneuron": "inhibitory", "astrocyte": "astrocyte",
              "microglial cell": "microglia", "central nervous system macrophage": "microglia",
              "oligodendrocyte precursor cell": "opc", "blood vessel endothelial cell": "endothelial",
              "endothelial cell of vascular tree": "endothelial"}
    for l, c in expect.items():
        assert celltypes.harmonize(l) == c
    assert celltypes.status("neuron")[0] == "excluded"
    for l in ("mural cell", "leukocyte", "central nervous system neuron", "ependymal cell"):
        assert celltypes.status(l)[0] == "excluded", l
    assert celltypes.status("dopaminergic neuron")[0] == "excluded"


def test_unknown_labels_are_reported_not_assigned():
    assert celltypes.harmonize("some new cell type") is None
    assert celltypes.status("some new cell type") == ("unmapped", None)


def test_classes_are_consistent_and_disjoint():
    assert set(celltypes.LOOKUP.values()) <= set(celltypes.CLASSES)
    assert not set(celltypes.LOOKUP) & set(celltypes.EXCLUDED)
    assert set(celltypes.PRIMARY_CLASSES) < set(celltypes.CLASSES)


def test_labels_for_returns_only_the_requested_class():
    obs = ["astrocyte", "oligodendrocyte", "microglial cell", "neuron"]
    assert celltypes.labels_for("astrocyte", obs) == ["astrocyte"]
    assert celltypes.labels_for("excitatory", obs) == []


def test_parse_age():
    assert parse_age("80 year-old human stage") == (80.0, False)
    assert parse_age("80 year-old and over human stage") == (80.0, True)
    assert parse_age("58-year-old human stage") == (58.0, False)
    a, c = parse_age("eighth decade human stage")
    assert math.isnan(a) and c is False


def test_cohort_registry_is_well_formed():
    for key, c in COHORTS.items():
        assert c["role"] in {"discovery", "replication", "exploratory"}, key
        assert c["disease"] in {"AD", "PD"} and c["tissues"], key
        for i in c["dataset_ids"]:
            assert len(i) == 36 and i.count("-") == 4, (key, i)
    assert {c["role"] for c in COHORTS.values()} & SEALED_ROLES
    # at least one discovery cohort per disease, so the primary comparison is possible
    assert {c["disease"] for c in COHORTS.values() if c["role"] == "discovery"} == {"AD", "PD"}
