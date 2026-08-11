"""
plotting.py
============

Two plot types:

- property_line_plot: one property (e.g. hydropathy) plotted along the
  sequence, smoothed line + raw per-residue scatter underneath. This is
  the classic Kyte-Doolittle-style hydropathy plot, generalized to any
  property in aa_data.py.
- composition_radar_plot: mean value of each property across the whole
  sequence, radar/spider chart -- a "fingerprint" of the sequence's
  overall physicochemical character. Useful for comparing two sequences
  (e.g. WT vs mutant, or two homologs) at a glance.

Both save directly to a file path (PNG by default) rather than calling
plt.show(), since this is meant to run headless / from a script.
"""

import math

import matplotlib.pyplot as plt
import numpy as np

from .aa_data import PROPERTY_FIELDS, PROPERTY_LABELS
from .profiler import SequenceProfile


def property_line_plot(profile: SequenceProfile, field: str, out_path: str, highlight_position: int | None = None) -> str:
    if field not in PROPERTY_FIELDS:
        raise ValueError(f"Unknown property '{field}'. Choose from {PROPERTY_FIELDS}")

    positions = list(range(1, len(profile.sequence) + 1))
    raw = profile.per_residue[field]
    smooth = profile.smoothed[field]

    fig, ax = plt.subplots(figsize=(max(8, len(positions) * 0.05), 4))
    ax.scatter(positions, raw, s=8, alpha=0.25, color="#4C72B0", label="per-residue")
    ax.plot(positions, smooth, linewidth=2, color="#C44E52",
            label=f"smoothed (window={profile.window_size})")
    ax.axhline(0, color="grey", linewidth=0.8, linestyle="--", alpha=0.6)

    if highlight_position is not None:
        ax.axvline(highlight_position, color="black", linewidth=1, linestyle=":", alpha=0.7)
        ax.annotate(
            f"pos {highlight_position}",
            xy=(highlight_position, ax.get_ylim()[1]),
            xytext=(highlight_position, ax.get_ylim()[1]),
            ha="center", va="bottom", fontsize=8,
        )

    ax.set_xlabel("Residue position")
    ax.set_ylabel(PROPERTY_LABELS[field])
    ax.set_title(f"{PROPERTY_LABELS[field]} along {profile.label}")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def composition_radar_plot(profiles: list[SequenceProfile], out_path: str) -> str:
    """
    Radar plot comparing whole-sequence mean composition across one or
    more SequenceProfile objects. Values are min-max normalized per
    field across the provided profiles so axes are visually comparable
    despite very different natural ranges (e.g. volume ~60-240 vs charge -1 to 1).
    """
    if not profiles:
        raise ValueError("Need at least one profile to plot.")

    fields = PROPERTY_FIELDS
    n_fields = len(fields)

    # normalize each field to [0, 1] across the provided profiles
    normalized = {p.label: [] for p in profiles}
    for field in fields:
        values = [p.composition_mean[field] for p in profiles]
        lo, hi = min(values), max(values)
        span = (hi - lo) or 1.0  # avoid div-by-zero if all equal
        for p, v in zip(profiles, values):
            normalized[p.label].append((v - lo) / span)

    angles = [n / float(n_fields) * 2 * math.pi for n in range(n_fields)]
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))
    colors = plt.cm.tab10(np.linspace(0, 1, len(profiles)))

    for p, color in zip(profiles, colors):
        values = normalized[p.label] + normalized[p.label][:1]
        ax.plot(angles, values, linewidth=2, label=p.label, color=color)
        ax.fill(angles, values, alpha=0.1, color=color)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels([PROPERTY_LABELS[f] for f in fields], fontsize=8)
    ax.set_yticklabels([])
    ax.set_title("Sequence composition profile\n(normalized across compared sequences)", fontsize=10, pad=20)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.15), ncol=min(len(profiles), 3), fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path
