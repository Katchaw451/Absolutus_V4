"""Donor-level covariate parsing."""
from __future__ import annotations

import math
import re

_AGE = re.compile(r"(\d+)[\s-]*year", re.IGNORECASE)


def parse_age(stage: str):
    """Age in years and whether it is a lower bound, from a Census `development_stage` string.

    '80 year-old human stage' -> (80.0, False); '80 year-old and over human stage' -> (80.0, True).
    Strings without an explicit number of years (e.g. 'eighth decade human stage') -> (nan, False)."""
    m = _AGE.search(str(stage))
    if not m:
        return math.nan, False
    return float(m.group(1)), "and over" in str(stage).lower()
