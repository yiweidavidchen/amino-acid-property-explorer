"""
export.py
=========

Write a SequenceProfile's per-residue data out to CSV.
"""

import csv
from pathlib import Path

from .aa_data import PROPERTY_FIELDS
from .profiler import SequenceProfile


def profile_to_csv(profile: SequenceProfile, out_path: str | Path) -> Path:
    out_path = Path(out_path)
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["position", "residue"] + PROPERTY_FIELDS + [f"{f}_smoothed" for f in PROPERTY_FIELDS]
        writer.writerow(header)
        for i, res in enumerate(profile.sequence):
            row = [i + 1, res]
            row += [profile.per_residue[f][i] for f in PROPERTY_FIELDS]
            row += [profile.smoothed[f][i] for f in PROPERTY_FIELDS]
            writer.writerow(row)
    return out_path
