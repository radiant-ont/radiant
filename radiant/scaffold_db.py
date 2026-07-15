# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""Assemble a germline-family scaffold database from OAS clones.

The output is a space-separated, eight-field scaffold database
(``name FR1 CDR1 FR2 CDR2 FR3 CDR3 FR4``) with one entry per germline V gene:
frameworks are taken from the germline reference and the CDR fields from a
representative clone of that gene. This is the same file format the workflow
uses for its production scaffolds; building one from public germline families
demonstrates that the workflow is not tied to any particular library.
"""

from __future__ import annotations

from .codon import back_translate
from .germline import GermlineFrameworks
from .oas import OasClone


def _representative_per_gene(clones: list[OasClone]) -> dict[str, OasClone]:
    """Pick one clone per gene, deterministically (most abundant first)."""
    best: dict[str, OasClone] = {}
    for clone in clones:
        current = best.get(clone.v_gene)
        key = (clone.redundancy, len(clone.cdr3_aa), clone.cdr3_aa)
        if current is None or key > (
            current.redundancy,
            len(current.cdr3_aa),
            current.cdr3_aa,
        ):
            best[clone.v_gene] = clone
    return best


def build_germline_scaffold_db(
    clones: list[OasClone],
    frameworks: dict[str, GermlineFrameworks],
) -> list[tuple[str, list[str]]]:
    """Return scaffold entries for every gene present in both inputs.

    Returns
    -------
    list of (name, fields)
        ``fields`` is ``[FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4]`` in nucleotides.
        Entries are sorted by gene name for a stable, diff-friendly output.
    """
    representatives = _representative_per_gene(clones)
    entries: list[tuple[str, list[str]]] = []
    for gene in sorted(representatives):
        fw = frameworks.get(gene)
        if fw is None:
            continue
        clone = representatives[gene]
        fields = [
            fw.fr1,
            back_translate(clone.cdr1_aa),
            fw.fr2,
            back_translate(clone.cdr2_aa),
            fw.fr3,
            back_translate(clone.cdr3_aa),
            fw.fr4,
        ]
        entries.append((f"germline_{gene}", fields))
    return entries


def write_scaffold_db(entries: list[tuple[str, list[str]]], path: str) -> None:
    """Write scaffold entries in the eight-field, space-separated format."""
    with open(path, "w") as handle:
        handle.write(
            "# germline-family scaffold database "
            "(name FR1 CDR1 FR2 CDR2 FR3 CDR3 FR4)\n"
        )
        for name, fields in entries:
            handle.write(" ".join([name, *fields]) + "\n")
