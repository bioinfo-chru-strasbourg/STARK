#!/usr/bin/env python

# Utility functions for the STARK API, not specific to any single router.


import os
from typing import Optional


def _safe_mtime(path: str) -> Optional[float]:
    try:
        return os.path.getmtime(path)
    except FileNotFoundError:
        return None


def is_empty(v):
    return v is None or (isinstance(v, str) and v.strip() == "")
