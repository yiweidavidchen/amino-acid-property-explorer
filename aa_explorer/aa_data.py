"""
aa_data.py
==========

Reference table of standard amino acid physicochemical properties, keyed
by one-letter code.

Sources / conventions (standard textbook values -- if you need
publication-grade precision, cross-check against the primary references
below rather than trusting this table blind):

- hydropathy: Kyte & Doolittle (1982) hydropathy index.
- volume: Zamyatnin (1972) residue volume, cubic angstroms.
- charge: approximate formal charge at physiological pH (~7.4). Note this
  is a simplification -- His is ~10% protonated at pH 7.4 and is
  represented here as a fractional charge (0.1) rather than a strict 0/+1,
  to reflect that ambiguity rather than hide it.
- aromaticity: 1.0 for Phe/Trp/Tyr (full ring systems with UV absorbance
  and pi-stacking capacity), 0.0 otherwise. His is NOT counted as aromatic
  here despite having an imidazole ring, since it lacks the extended
  pi-system relevant to most "aromaticity" analyses (stacking, UV
  absorbance) -- adjust if your use case wants it included.
- mol_weight: residue mass (i.e. amino acid mass minus water lost in
  peptide bond formation), daltons.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class AAProperties:
    one_letter: str
    three_letter: str
    name: str
    hydropathy: float      # Kyte-Doolittle scale
    volume: float           # Zamyatnin, Angstrom^3
    charge: float             # approx formal charge at pH 7.4
    aromaticity: float         # 1.0 aromatic, 0.0 not
    mol_weight: float           # residue mass, Da


AMINO_ACIDS: dict[str, AAProperties] = {
    "A": AAProperties("A", "Ala", "Alanine",        1.8,  92.5,  0.0, 0.0,  71.08),
    "R": AAProperties("R", "Arg", "Arginine",       -4.5, 202.1,  1.0, 0.0, 156.19),
    "N": AAProperties("N", "Asn", "Asparagine",     -3.5, 135.2,  0.0, 0.0, 114.10),
    "D": AAProperties("D", "Asp", "Aspartate",      -3.5, 124.5, -1.0, 0.0, 115.09),
    "C": AAProperties("C", "Cys", "Cysteine",        2.5, 106.0,  0.0, 0.0, 103.14),
    "Q": AAProperties("Q", "Gln", "Glutamine",      -3.5, 161.1,  0.0, 0.0, 128.13),
    "E": AAProperties("E", "Glu", "Glutamate",      -3.5, 155.1, -1.0, 0.0, 129.12),
    "G": AAProperties("G", "Gly", "Glycine",        -0.4,  66.0,  0.0, 0.0,  57.05),
    "H": AAProperties("H", "His", "Histidine",      -3.2, 167.3,  0.1, 0.0, 137.14),
    "I": AAProperties("I", "Ile", "Isoleucine",      4.5, 168.8,  0.0, 0.0, 113.16),
    "L": AAProperties("L", "Leu", "Leucine",         3.8, 167.9,  0.0, 0.0, 113.16),
    "K": AAProperties("K", "Lys", "Lysine",         -3.9, 171.3,  1.0, 0.0, 128.17),
    "M": AAProperties("M", "Met", "Methionine",      1.9, 170.8,  0.0, 0.0, 131.19),
    "F": AAProperties("F", "Phe", "Phenylalanine",   2.8, 203.4,  0.0, 1.0, 147.18),
    "P": AAProperties("P", "Pro", "Proline",        -1.6, 129.0,  0.0, 0.0,  97.12),
    "S": AAProperties("S", "Ser", "Serine",         -0.8,  91.0,  0.0, 0.0,  87.08),
    "T": AAProperties("T", "Thr", "Threonine",      -0.7, 122.1,  0.0, 0.0, 101.10),
    "W": AAProperties("W", "Trp", "Tryptophan",     -0.9, 237.6,  0.0, 1.0, 186.21),
    "Y": AAProperties("Y", "Tyr", "Tyrosine",       -1.3, 203.6,  0.0, 1.0, 163.18),
    "V": AAProperties("V", "Val", "Valine",          4.2, 141.7,  0.0, 0.0,  99.13),
}

PROPERTY_FIELDS = ["hydropathy", "volume", "charge", "aromaticity", "mol_weight"]

PROPERTY_LABELS = {
    "hydropathy": "Hydropathy (Kyte-Doolittle)",
    "volume": "Volume (\u00c5\u00b3, Zamyatnin)",
    "charge": "Charge at pH 7.4",
    "aromaticity": "Aromaticity",
    "mol_weight": "Molecular weight (Da)",
}


def get_property(one_letter: str, field: str) -> float:
    code = one_letter.upper()
    if code not in AMINO_ACIDS:
        raise KeyError(f"Unknown amino acid code: '{one_letter}'")
    return getattr(AMINO_ACIDS[code], field)
