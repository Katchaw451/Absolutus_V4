"""Cross-disease combination of per-gene z-scores."""
from __future__ import annotations

import numpy as np


def conservative_combine(z_a, z_b):
    """Signed minimum: keeps min(|z_a|, |z_b|) when signs agree, else 0.
    A gene scores high only if it is strongly and concordantly altered in both."""
    z_a, z_b = np.asarray(z_a, float), np.asarray(z_b, float)
    same = np.sign(z_a) == np.sign(z_b)
    return np.where(same, np.minimum(np.abs(z_a), np.abs(z_b)) * np.sign(z_a), 0.0)
