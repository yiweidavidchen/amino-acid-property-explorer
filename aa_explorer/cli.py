"""
cli.py
======

Usage:
    # Plot one property along a sequence
    python -m aa_explorer.cli profile --seq MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWELVMGDGERAFSTLTETIEAVWKGADFEIYETLKQR --field hydropathy -o hydropathy.png

    # Plot from a FASTA file (uses first record, or --record to pick by header)
    python -m aa_explorer.cli profile --fasta examples/example.fasta --field charge -o charge.png

    # Export full per-residue table to CSV
    python -m aa_explorer.cli export --fasta examples/example.fasta -o profile.csv

    # Compare composition of multiple sequences on one radar plot
    python -m aa_explorer.cli radar --fasta examples/example.fasta -o radar.png

    # Single point-mutation property delta report
    python -m aa_explorer.cli mutate --seq MKTAYIAK... --position 5 --to D
"""

import argparse
import json

from .aa_data import PROPERTY_FIELDS
from .export import profile_to_csv
from .plotting import composition_radar_plot, property_line_plot
from .profiler import compare_point_mutation, compute_profile
from .sequence_io import clean_sequence, read_fasta


def _load_sequences(args) -> dict[str, str]:
    if args.fasta:
        records = read_fasta(args.fasta)
        if getattr(args, "record", None):
            if args.record not in records:
                raise KeyError(f"Record '{args.record}' not found in {args.fasta}. "
                                f"Available: {list(records.keys())}")
            return {args.record: records[args.record]}
        return records
    elif args.seq:
        return {"sequence": clean_sequence(args.seq)}
    else:
        raise ValueError("Provide either --seq or --fasta")


def cmd_profile(args):
    sequences = _load_sequences(args)
    label, seq = next(iter(sequences.items()))
    profile = compute_profile(seq, label=label, window_size=args.window)
    out = property_line_plot(profile, args.field, args.output, highlight_position=args.highlight)
    print(f"Wrote {out}")


def cmd_export(args):
    sequences = _load_sequences(args)
    label, seq = next(iter(sequences.items()))
    profile = compute_profile(seq, label=label, window_size=args.window)
    out = profile_to_csv(profile, args.output)
    print(f"Wrote {out}")


def cmd_radar(args):
    sequences = _load_sequences(args)
    profiles = [compute_profile(seq, label=label, window_size=args.window)
                for label, seq in sequences.items()]
    out = composition_radar_plot(profiles, args.output)
    print(f"Wrote {out} ({len(profiles)} sequence(s) compared)")


def cmd_mutate(args):
    seq = clean_sequence(args.seq) if args.seq else next(iter(read_fasta(args.fasta).values()))
    result = compare_point_mutation(seq, args.position, args.to)
    print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser(prog="aa-explorer")
    sub = parser.add_subparsers(dest="command", required=True)

    def add_seq_args(p):
        p.add_argument("--seq", help="Raw sequence string")
        p.add_argument("--fasta", help="Path to FASTA file")
        p.add_argument("--record", help="Header of the FASTA record to use (default: first)")
        p.add_argument("--window", type=int, default=9, help="Smoothing window size (odd integer, default 9)")

    p_profile = sub.add_parser("profile", help="Plot one property along the sequence")
    add_seq_args(p_profile)
    p_profile.add_argument("--field", required=True, choices=PROPERTY_FIELDS)
    p_profile.add_argument("--highlight", type=int, help="1-indexed position to mark on the plot")
    p_profile.add_argument("-o", "--output", required=True)
    p_profile.set_defaults(func=cmd_profile)

    p_export = sub.add_parser("export", help="Export per-residue property table to CSV")
    add_seq_args(p_export)
    p_export.add_argument("-o", "--output", required=True)
    p_export.set_defaults(func=cmd_export)

    p_radar = sub.add_parser("radar", help="Composition radar plot (compares all records if --fasta has several)")
    add_seq_args(p_radar)
    p_radar.add_argument("-o", "--output", required=True)
    p_radar.set_defaults(func=cmd_radar)

    p_mutate = sub.add_parser("mutate", help="Report property deltas for a single point mutation")
    p_mutate.add_argument("--seq", help="Raw sequence string")
    p_mutate.add_argument("--fasta", help="Path to FASTA file (uses first record)")
    p_mutate.add_argument("--position", type=int, required=True, help="1-indexed residue position")
    p_mutate.add_argument("--to", required=True, help="Mutant residue one-letter code")
    p_mutate.set_defaults(func=cmd_mutate)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
