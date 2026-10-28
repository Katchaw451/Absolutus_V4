"""Harmonisation of Census `cell_type` labels into broad classes shared across cohorts.

The mapping is fixed here, in version control, BEFORE any analysis. Labels not listed are
reported as `unmapped` and dropped (never silently assigned). Mapping choices that need
scrutiny are flagged in the comments and in docs/PROTOCOL.md.
"""
from __future__ import annotations

CLASSES = ("excitatory", "inhibitory", "astrocyte", "microglia", "oligodendrocyte", "opc", "endothelial")
PRIMARY_CLASSES = CLASSES[:6]  # endothelial is too sparse in several cohorts for the primary analysis

_MAP = {
    "excitatory": [
        "glutamatergic neuron",
        "L2/3-6 intratelencephalic projecting glutamatergic neuron",
        "corticothalamic-projecting glutamatergic cortical neuron",
        "near-projecting glutamatergic cortical neuron",
        "L6b glutamatergic cortical neuron",
        "L5 extratelencephalic projecting glutamatergic cortical neuron",
    ],
    "inhibitory": [
        "GABAergic neuron",
        "inhibitory interneuron",
        "pvalb GABAergic cortical interneuron",
        "VIP GABAergic cortical interneuron",
        "lamp5 GABAergic cortical interneuron",
        "sst GABAergic cortical interneuron",
        "sncg GABAergic cortical interneuron",
        "chandelier pvalb GABAergic cortical interneuron",
        "caudal ganglionic eminence derived cortical interneuron",
    ],
    "astrocyte": ["astrocyte", "astrocyte of the cerebral cortex"],
    # FLAG: the PD atlas labels its myeloid cells "central nervous system macrophage"; they are
    # pooled with microglia here. Treat as a sensitivity item (see protocol).
    "microglia": ["microglial cell", "central nervous system macrophage"],
    "oligodendrocyte": ["oligodendrocyte"],
    "opc": ["oligodendrocyte precursor cell"],
    "endothelial": [
        "cerebral cortex endothelial cell",
        "endothelial cell of vascular tree",
        "blood vessel endothelial cell",
        "endothelial cell",
    ],
}

EXCLUDED = {
    "vascular leptomeningeal cell": "non-parenchymal cell type",
    "pericyte": "mural cell, not in the primary classes",
    "T cell": "infiltrating immune cell",
    "ependymal cell": "rare, not in the primary classes",
    "macrophage": "border/peripheral macrophage (rare); not pooled with microglia",
    "neuron": "unspecified neuron type (Kamath); cannot be assigned to excitatory/inhibitory",
    "dopaminergic neuron": "analysed separately (substantia nigra specific)",
}

LOOKUP = {label: cls for cls, labels in _MAP.items() for label in labels}


def harmonize(label: str):
    """Broad class for a Census cell_type label, or None."""
    return LOOKUP.get(label)


def status(label: str):
    """('class', name) | ('excluded', reason) | ('unmapped', None)."""
    if label in LOOKUP:
        return "class", LOOKUP[label]
    if label in EXCLUDED:
        return "excluded", EXCLUDED[label]
    return "unmapped", None


def labels_for(cls: str, observed) -> list[str]:
    """Observed labels that belong to class `cls`, sorted."""
    return sorted(l for l in observed if LOOKUP.get(l) == cls)
