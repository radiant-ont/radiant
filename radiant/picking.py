# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""Lead selection on the re-annotated consensus sequences.

Keep functional only (no stop codons/frameshifts/truncations), group by 100%
HCDR3 AA identity keeping the most abundant per cluster, rank by frequency, and
optionally take the top N. Input is the annotated consensus table with a
cluster-size column.
"""

from __future__ import annotations

import re

import pandas as pd

ABERRANT_AA = re.compile(r"[*X#]")


def is_functional(row: pd.Series, single_chain: bool) -> bool:
    """Return True if an annotated consensus row is a clean functional sequence."""
    issue = row.get("sequence_issue", "")
    # A missing/blank sequence_issue means "no issue"; NaN must not be read as text.
    if not pd.isna(issue) and str(issue).strip():
        return False
    chain_cols = ["functional_1"] if single_chain else ["functional_1", "functional_2"]
    for col in chain_cols:
        if col in row and str(row[col]).lower() != "functional":
            return False
    merged = str(row.get("Merged_CDRs_AA", "") or "")
    if not merged or ABERRANT_AA.search(merged):
        return False
    return True


def select_functional(df: pd.DataFrame, single_chain: bool) -> pd.DataFrame:
    """Keep only functional consensus rows."""
    mask = df.apply(lambda r: is_functional(r, single_chain), axis=1)
    return df[mask].copy()


def dedup_by_hcdr3(
    df: pd.DataFrame,
    size_col: str = "cluster_size",
    hcdr3_col: str = "HCDR3_AA",
) -> pd.DataFrame:
    """Collapse to one row per unique HCDR3 (100% identity).

    Sequences sharing an identical HCDR3 amino-acid sequence form one functional
    cluster. The representative kept is the most-abundant member, and the
    cluster's population frequency (``hcdr3_reads``) is the total reads across
    every cluster sharing that HCDR3, the manuscript ranks clusters by the
    relative population frequency of the HCDR3 sequence.
    """
    hcdr3_reads = df.groupby(hcdr3_col)[size_col].sum()
    representative = (
        df.sort_values(size_col, ascending=False, kind="stable")
        .drop_duplicates(subset=hcdr3_col, keep="first")
        .copy()
    )
    representative["hcdr3_reads"] = representative[hcdr3_col].map(hcdr3_reads)
    return representative


def rank_top_n(
    df: pd.DataFrame,
    freq_col: str = "hcdr3_reads",
    top_n: int | None = None,
) -> pd.DataFrame:
    """Rank by HCDR3 population frequency; add ``pick_rank`` and ``rel_freq``."""
    ranked = df.sort_values(freq_col, ascending=False, kind="stable").reset_index(drop=True)
    total = ranked[freq_col].sum()
    ranked["rel_freq"] = ranked[freq_col] / total if total else 0.0
    ranked["pick_rank"] = range(1, len(ranked) + 1)
    if top_n is not None:
        ranked = ranked.head(top_n).copy()
    return ranked


def pick_leads(
    df: pd.DataFrame,
    single_chain: bool,
    top_n: int | None = None,
    size_col: str = "cluster_size",
    hcdr3_col: str = "HCDR3_AA",
) -> pd.DataFrame:
    """Run the full selection: functional -> HCDR3 grouping -> rank -> optional top N."""
    functional = select_functional(df, single_chain)
    deduped = dedup_by_hcdr3(functional, size_col=size_col, hcdr3_col=hcdr3_col)
    return rank_top_n(deduped, top_n=top_n)
