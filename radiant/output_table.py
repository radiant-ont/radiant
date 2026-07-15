# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""Assemble the operator-facing output table from the annotated consensus.

This reproduces the columns of the supplementary ONT output table: cluster rank
and size, full-length and per-chain sequences (AA and nucleotide), germline match
names, per-region frameworks and CDRs, votes, ROI and HCDR3 cluster labels,
liabilities, biophysical properties, and a confidence score.

The confidence score follows the legacy definition: a power law is fit to the
cluster-size distribution and each cluster is scored

    confidence = (1 - (cdf(rank) - cdf(1))) * count

with cdf the power-law cumulative distribution. The power-law exponent is
estimated by maximum likelihood (xmin = 1), so no external fitting package is
required.
"""

from __future__ import annotations

import math

import pandas as pd


def _powerlaw_alpha(sizes: list[int]) -> float:
    """MLE exponent of a power law with xmin = 1 over the cluster sizes."""
    xs = [s for s in sizes if s >= 1]
    if len(xs) < 2:
        return 2.0
    denom = sum(math.log(x) for x in xs)
    if denom <= 0:
        return 2.0
    return 1.0 + len(xs) / denom


def _confidence_scores(rank: pd.Series, count: pd.Series) -> pd.Series:
    """Power-law confidence: (1 - (cdf(rank) - cdf(1))) * count, xmin = 1."""
    alpha = _powerlaw_alpha([int(c) for c in count])
    exp = 1.0 - alpha  # 1 - cdf(x) = x^(1-alpha) for xmin = 1
    # cdf(1) = 0, so scaler = 1 - cdf(rank) = rank^(1-alpha)
    scaler = rank.astype(float).pow(exp).round(3)
    return (scaler * count.astype(float)).round(3)


# (display name, source column). Missing source columns are filled with "".
def _column_spec(single_chain: bool) -> list[tuple[str, str]]:
    spec: list[tuple[str, str]] = [
        ("Cluster Rank", "cluster_rank"),
        ("Full Length Sequence AA", "sequence_aa_1_2"),
    ]
    if not single_chain:
        spec += [
            ("VL Sequence AA", "sequence_aa_1"),
            ("VH Sequence AA", "sequence_aa_2"),
            ("VL Sequence Nucleotide", "sequence_1"),
            ("VH Sequence Nucleotide", "sequence_2"),
            ("VL Match Name", "match_name_1"),
            ("VH Match Name", "match_name_2"),
            ("VL/VH Match Name", "match_name_1_2"),
        ]
        chains = (("VL", "1"), ("VH", "2"))
    else:
        spec += [
            ("VH Sequence AA", "sequence_aa_1"),
            ("VH Sequence Nucleotide", "sequence_1"),
            ("VH Match Name", "match_name_1"),
        ]
        chains = (("VH", "1"),)

    for label, idx in chains:
        for i in (1, 2, 3, 4):
            spec.append((f"{label} FW{i} Nucleotide", f"fr{i}_{idx}"))
        for i in (1, 2, 3):
            spec.append((f"{label} CDR{i} Nucleotide", f"cdr{i}_{idx}"))
        for i in (1, 2, 3, 4):
            spec.append((f"{label} FW{i} Amino Acid", f"fr{i}_aa_{idx}"))
        for i in (1, 2, 3):
            spec.append((f"{label} CDR{i} Amino Acid", f"cdr{i}_aa_{idx}"))

    spec.append(("Count", "count"))
    if not single_chain:
        spec += [("Votes VL", "votes_1"), ("Votes VH", "votes_2")]
    else:
        spec += [("Votes VH", "votes_1")]

    spec += [("ROI Cluster", "roi_cluster"), ("HCDR3 Cluster", "hcdr3_cluster")]

    lchains = (("L", "1"), ("H", "2")) if not single_chain else (("H", "1"),)
    for label, idx in lchains:
        for i in (1, 2, 3):
            spec.append((f"Liability {label}CDR{i}", f"liability_string_cdr{i}_aa_{idx}"))
    spec.append(("Liability Quantification", "liability_quant_lcdr1_3_hcdr1_3"))

    for name in (
        "isoelectric_point", "high_viscosity_index", "charge_symmetric_parameter",
        "cdr3_aa_2_charge", "cdr3_aa_2_hydropathy", "cdr3_aa_2_length",
    ):
        spec.append((f"Biophysical: {name}", name))

    spec.append(("Confidence Score", "confidence_score"))
    return spec


def build_output_table(annotated: pd.DataFrame, single_chain: bool = False) -> pd.DataFrame:
    """Build the supplementary-format output table from an annotated consensus set.

    ``annotated`` is the re-annotated consensus (one row per consensus), and must
    carry a ``count`` (or ``cluster_size``) abundance column. Rows are ranked by
    abundance; ROI and HCDR3 cluster labels and the confidence score are derived.
    """
    df = annotated.copy()
    if "count" not in df.columns:
        df["count"] = df.get("cluster_size", 1)
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(1).astype(int)

    df = df.sort_values("count", ascending=False, kind="stable").reset_index(drop=True)
    df["cluster_rank"] = df.index + 1

    # ROI cluster: keep the clustering label if present, else the rank.
    df["roi_cluster"] = df["cluster_numeric"] if "cluster_numeric" in df.columns else df["cluster_rank"]
    # HCDR3 cluster: one label per distinct heavy-chain CDR3 amino-acid sequence.
    h3 = df["HCDR3_AA"] if "HCDR3_AA" in df.columns else df.get("cdr3_aa_2", df.get("cdr3_aa_1", ""))
    df["hcdr3_cluster"] = pd.Series(h3).astype(str).map(
        {v: i + 1 for i, v in enumerate(pd.Series(h3).astype(str).drop_duplicates())}
    )

    df["confidence_score"] = _confidence_scores(df["cluster_rank"], df["count"])

    out = pd.DataFrame()
    for display, source in _column_spec(single_chain):
        out[display] = df[source] if source in df.columns else ""
    return out
