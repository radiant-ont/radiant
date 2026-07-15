# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""Read antibody clones from an OAS-style CSV subset.

Expected columns: ``v_call``, ``cdr1_aa``, ``cdr2_aa``, ``cdr3_aa`` and an
optional ``Redundancy`` count. Only the heavy chain is modelled here (the
subset is heavy-only), so clones are treated as single-domain, matching the
VHH/single-chain path of the workflow.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass


@dataclass(frozen=True)
class OasClone:
    """One heavy-chain clone: germline gene plus observed CDR amino acids."""

    v_gene: str
    cdr1_aa: str
    cdr2_aa: str
    cdr3_aa: str
    redundancy: int

    @property
    def cdr_key(self) -> tuple[str, str, str]:
        """Identity of the clone at the amino-acid CDR level."""
        return (self.cdr1_aa, self.cdr2_aa, self.cdr3_aa)


_VALID_AA = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


def _bare_gene(v_call: str) -> str:
    """``IGHV3-23*01`` -> ``IGHV3-23``."""
    return v_call.split("*", 1)[0].strip().upper()


def load_oas_clones(path: str) -> list[OasClone]:
    """Load clones from an OAS CSV, skipping rows with non-standard residues.

    Rows whose CDRs contain characters outside the twenty standard amino acids
    (for example ``X`` or ``*``) are skipped so that back-translation cannot
    fail downstream.
    """
    clones: list[OasClone] = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            cdr1 = row["cdr1_aa"].strip().upper()
            cdr2 = row["cdr2_aa"].strip().upper()
            cdr3 = row["cdr3_aa"].strip().upper()
            if not all(_VALID_AA.match(c) for c in (cdr1, cdr2, cdr3)):
                continue
            try:
                redundancy = int(row.get("Redundancy", "1") or "1")
            except ValueError:
                redundancy = 1
            clones.append(
                OasClone(
                    v_gene=_bare_gene(row["v_call"]),
                    cdr1_aa=cdr1,
                    cdr2_aa=cdr2,
                    cdr3_aa=cdr3,
                    redundancy=redundancy,
                )
            )
    return clones
