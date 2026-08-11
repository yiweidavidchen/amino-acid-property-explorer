"""
sequence_io.py
===============

Minimal FASTA reading and sequence validation. Deliberately doesn't
depend on Biopython -- this tool has zero heavy dependencies by design.
"""

from pathlib import Path

from .aa_data import AMINO_ACIDS

VALID_CODES = set(AMINO_ACIDS.keys())


def clean_sequence(raw: str) -> str:
    """Uppercase, strip whitespace/newlines from a raw sequence string."""
    return "".join(raw.split()).upper()


def validate_sequence(seq: str) -> list[str]:
    """Return a list of invalid characters found in the sequence (empty if valid)."""
    return sorted({c for c in seq if c not in VALID_CODES})


def read_fasta(path: str | Path) -> dict[str, str]:
    """
    Parse a (possibly multi-record) FASTA file into {header: sequence}.
    Header is the text after '>' on the record's first line, up to the
    first whitespace.
    """
    path = Path(path)
    records: dict[str, str] = {}
    header = None
    seq_chunks: list[str] = []

    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records[header] = clean_sequence("".join(seq_chunks))
            header = line[1:].split()[0] if len(line) > 1 else "unnamed"
            seq_chunks = []
        else:
            seq_chunks.append(line)

    if header is not None:
        records[header] = clean_sequence("".join(seq_chunks))

    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    return records
