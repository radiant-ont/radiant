# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""Load germline framework nucleotides from a space-separated reference.

The reference format is one record per line:

    name FR1 CDR1 FR2 CDR2 FR3 CDR3 FR4

which is the same eight-field, space-separated layout the annotator consumes
for a scaffold database. IMGT alignment gaps (``-``) and any placeholder bases
are permitted in the reference and are removed on load. This module reads only
the framework regions (FR1-FR4); CDRs for a given clone come from the antibody
data being modelled (for example OAS), not from the germline reference.
"""

from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass(frozen=True)
class GermlineFrameworks:
    """Framework nucleotides (FR1-FR4) for one germline V/J scaffold."""

    name: str
    fr1: str
    fr2: str
    fr3: str
    fr4: str


def _clean(seq: str) -> str:
    """Strip alignment gaps and placeholder N runs; upper-case."""
    return re.sub(r"[^ACGT]", "", seq.upper())


def _gene_key(name: str) -> str:
    """Reduce a reference name to a bare gene symbol.

    ``human-IGHV3-23`` -> ``IGHV3-23``. Names without a recognisable IG gene
    are returned unchanged (upper-cased), so custom references still load.
    """
    match = re.search(r"IG[HKL][VDJ][0-9][0-9A-Za-z\-]*", name)
    return match.group(0).upper() if match else name.upper()


def load_germline_frameworks(path: str) -> dict[str, GermlineFrameworks]:
    """Read framework nucleotides keyed by bare gene symbol.

    Parameters
    ----------
    path:
        Path to a space-separated germline reference (see module docstring).

    Returns
    -------
    dict
        Mapping of gene symbol (for example ``IGHV3-23``) to
        :class:`GermlineFrameworks`. When several alleles of one gene are
        present the first encountered is kept, which is deterministic given a
        stable input ordering.
    """
    frameworks: dict[str, GermlineFrameworks] = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 8:
                continue
            name, fr1, _cdr1, fr2, _cdr2, fr3, _cdr3, fr4 = fields
            key = _gene_key(name)
            if key in frameworks:
                continue
            frameworks[key] = GermlineFrameworks(
                name=key,
                fr1=_clean(fr1),
                fr2=_clean(fr2),
                fr3=_clean(fr3),
                fr4=_clean(fr4),
            )
    return frameworks
