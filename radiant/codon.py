# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""Deterministic amino-acid to nucleotide back-translation.

A single fixed codon is used per amino acid (a frequent human codon) so that
back-translation is reproducible and independent of any external table. This
is used only to render CDR amino-acid sequences (for example from OAS) into
nucleotides when assembling scaffold-database entries and simulated reads; it
is not intended as a codon-optimization tool.
"""

from __future__ import annotations

# One frequent human codon per amino acid. Deterministic by construction.
CODON_TABLE: dict[str, str] = {
    "A": "GCC", "R": "CGG", "N": "AAC", "D": "GAC", "C": "TGC",
    "Q": "CAG", "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC",
    "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
    "*": "TAA",
}


def back_translate(aa: str) -> str:
    """Return a deterministic nucleotide sequence for an amino-acid string.

    Parameters
    ----------
    aa:
        Amino-acid sequence (one-letter codes, case-insensitive). Non-standard
        characters raise ``KeyError`` so that malformed input fails loudly
        rather than silently producing an incorrect sequence.

    Returns
    -------
    str
        Upper-case nucleotide sequence, three bases per residue.
    """
    aa = aa.strip().upper()
    return "".join(CODON_TABLE[residue] for residue in aa)
