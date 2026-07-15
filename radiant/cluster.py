# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""Cluster annotated reads by combined-CDR nucleotide distance.

Greedy single-linkage: the most abundant unassigned sequence seeds a cluster,
and anything within LD of that seed joins it (Levenshtein over the combined-CDR
nucleotide region, LD<=3 for scFv and <=2 for VHH). Clusters below min_size are
dropped.
"""

from __future__ import annotations

from dataclasses import dataclass

import editdistance
import pandas as pd

DEFAULT_ROI = "Merged_CDRs_NUC"
DEFAULT_CUTOFF = 3
DEFAULT_MIN_SIZE = 5
DEFAULT_METRIC = "levenshtein"
METRICS = ("levenshtein", "hamming")


def _hamming(a: str, b: str) -> int:
    """Hamming distance (substitutions only). Sequences of unequal length are
    treated as far apart: a length difference cannot be bridged by substitutions,
    so the distance is reported as the longer length, which exceeds any sane
    cutoff and keeps them in separate clusters."""
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(1 for x, y in zip(a, b) if x != y)


def _distance(a: str, b: str, metric: str) -> int:
    """Distance between two region strings under the selected metric."""
    if metric == "hamming":
        return _hamming(a, b)
    return editdistance.eval(a, b)  # levenshtein (default)


def cluster_by_roi(
    df: pd.DataFrame,
    roi_col: str = DEFAULT_ROI,
    cutoff: int = DEFAULT_CUTOFF,
    count_col: str = "count",
    metric: str = DEFAULT_METRIC,
) -> pd.DataFrame:
    """Assign a ``cluster_numeric`` label to each row.

    Rows with a missing region are dropped. The seed of each cluster is the
    most-abundant unassigned row; membership is decided by distance to that
    seed only (not transitively), matching the reference implementation.

    Parameters
    ----------
    df:
        Annotated reads; must contain ``roi_col`` and ``count_col``.
    roi_col:
        Region used for clustering (default combined-CDR nucleotides).
    cutoff:
        Maximum edit distance to the seed for cluster membership.
    count_col:
        Abundance column used to order seeds.

    Returns
    -------
    pandas.DataFrame
        Copy of the input (region-complete rows only) with a ``cluster_numeric``
        column, ordered by descending abundance.
    """
    work = df.dropna(subset=[roi_col]).copy()
    work = work[work[roi_col].astype(str).str.len() > 0]
    work = work.sort_values(count_col, ascending=False, kind="stable").reset_index(drop=True)

    rois = work[roi_col].astype(str).tolist()
    labels: list[int] = [-1] * len(work)
    next_label = 0
    for i in range(len(work)):
        if labels[i] != -1:
            continue
        labels[i] = next_label
        seed = rois[i]
        for j in range(i + 1, len(work)):
            if labels[j] != -1:
                continue
            if _distance(rois[j], seed, metric) <= cutoff:
                labels[j] = next_label
        next_label += 1

    work["cluster_numeric"] = labels
    return work


def cluster_by_scaffold_and_roi(
    df: pd.DataFrame,
    scaffold_col: str,
    roi_col: str = DEFAULT_ROI,
    cutoff: int = DEFAULT_CUTOFF,
    count_col: str = "count",
    metric: str = DEFAULT_METRIC,
) -> pd.DataFrame:
    """Cluster within each scaffold-identity group, then combine.

    Reads are first grouped by scaffold identity (``scaffold_col``) and clustered
    independently within each group, matching the manuscript's "clustered by
    scaffold identity and Levenshtein or hamming distance similarity to the
    seed". Cluster labels are made globally unique across scaffolds.
    """
    parts: list[pd.DataFrame] = []
    offset = 0
    for _, group in df.groupby(scaffold_col, sort=True):
        clustered = cluster_by_roi(
            group, roi_col=roi_col, cutoff=cutoff, count_col=count_col, metric=metric
        )
        if clustered.empty:
            continue
        clustered = clustered.copy()
        clustered["cluster_numeric"] = clustered["cluster_numeric"] + offset
        offset = int(clustered["cluster_numeric"].max()) + 1
        parts.append(clustered)
    if not parts:
        return df.iloc[0:0].assign(cluster_numeric=pd.Series(dtype=int))
    return pd.concat(parts, ignore_index=True)


@dataclass(frozen=True)
class ClusterView:
    """A retained cluster: its label, members (abundance-sorted), and read size.

    ``size`` is the number of reads (sum of the ``count`` column), not the number
    of distinct sequences, so the minimum-size threshold is applied on read depth
    exactly as in the reference pipeline (where the annotated table is per-read).
    """

    cluster_numeric: int
    rank: int
    members: pd.DataFrame
    size: int

    @property
    def seed(self) -> pd.Series:
        """The most-abundant member (the reference for consensus)."""
        return self.members.iloc[0]


def filter_and_rank_clusters(
    clustered: pd.DataFrame,
    min_size: int = DEFAULT_MIN_SIZE,
    count_col: str = "count",
) -> list[ClusterView]:
    """Keep clusters with at least ``min_size`` reads, ranked by read depth.

    Cluster size is the sum of ``count`` (reads) across members; members are
    ordered by descending abundance so the first row is the seed.
    """
    read_sizes = clustered.groupby("cluster_numeric")[count_col].sum()
    kept = [c for c, n in read_sizes.items() if n >= min_size]
    ordered = sorted(kept, key=lambda c: (-int(read_sizes[c]), int(c)))

    views: list[ClusterView] = []
    for rank, cluster_numeric in enumerate(ordered, start=1):
        members = clustered[clustered["cluster_numeric"] == cluster_numeric]
        members = members.sort_values(count_col, ascending=False, kind="stable")
        views.append(
            ClusterView(
                int(cluster_numeric),
                rank,
                members.reset_index(drop=True),
                int(read_sizes[cluster_numeric]),
            )
        )
    return views
