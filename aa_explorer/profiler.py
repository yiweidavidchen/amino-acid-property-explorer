"""
profiler.py
============

Computes per-residue property arrays along a sequence, optional sliding-
window smoothing (the classic way hydropathy plots are made readable --
raw per-residue values are noisy), and whole-sequence composition
summaries for the radar plot.
"""

from dataclasses import dataclass

from .aa_data import AMINO_ACIDS, PROPERTY_FIELDS, get_property
from .sequence_io import validate_sequence


@dataclass
class SequenceProfile:
    sequence: str
    label: str
    per_residue: dict[str, list[float]]      # field -> value at each position
    smoothed: dict[str, list[float]]           # field -> sliding-window average
    composition_mean: dict[str, float]           # field -> mean over whole sequence
    window_size: int


def compute_profile(sequence: str, label: str = "sequence", window_size: int = 9) -> SequenceProfile:
    invalid = validate_sequence(sequence)
    if invalid:
        raise ValueError(
            f"Sequence contains non-standard amino acid code(s): {invalid}. "
            "Only the 20 standard one-letter codes are supported."
        )
    if len(sequence) == 0:
        raise ValueError("Sequence is empty.")
    if window_size < 1 or window_size % 2 == 0:
        raise ValueError("window_size must be a positive odd integer (e.g. 5, 9, 19).")

    per_residue: dict[str, list[float]] = {}
    for field in PROPERTY_FIELDS:
        per_residue[field] = [get_property(res, field) for res in sequence]

    smoothed: dict[str, list[float]] = {}
    half = window_size // 2
    n = len(sequence)
    for field, values in per_residue.items():
        smoothed_values = []
        for i in range(n):
            lo = max(0, i - half)
            hi = min(n, i + half + 1)
            window = values[lo:hi]
            smoothed_values.append(sum(window) / len(window))
        smoothed[field] = smoothed_values

    composition_mean = {
        field: sum(values) / len(values) for field, values in per_residue.items()
    }

    return SequenceProfile(
        sequence=sequence,
        label=label,
        per_residue=per_residue,
        smoothed=smoothed,
        composition_mean=composition_mean,
        window_size=window_size,
    )


def compare_point_mutation(sequence: str, position: int, mutant_residue: str) -> dict:
    """
    Report the property deltas for a single-position substitution.
    `position` is 1-indexed, matching how people normally refer to
    residue positions (E406D, not sequence[405]).
    """
    if position < 1 or position > len(sequence):
        raise ValueError(f"Position {position} out of range for sequence of length {len(sequence)}.")

    idx = position - 1
    wt_residue = sequence[idx]
    mutant_residue = mutant_residue.upper()

    if mutant_residue not in AMINO_ACIDS:
        raise ValueError(f"Unknown mutant residue code: '{mutant_residue}'")

    deltas = {}
    for field in PROPERTY_FIELDS:
        wt_val = get_property(wt_residue, field)
        mut_val = get_property(mutant_residue, field)
        deltas[field] = {
            "wt": wt_val,
            "mutant": mut_val,
            "delta": mut_val - wt_val,
        }

    return {
        "position": position,
        "wt_residue": wt_residue,
        "mutant_residue": mutant_residue,
        "label": f"{wt_residue}{position}{mutant_residue}",
        "deltas": deltas,
    }
