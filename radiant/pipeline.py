# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""End-to-end RADiANT pipeline: trim, filter, annotate, cluster, consensus.

This orchestrates the open stages together with the neutral annotator and
returns one row per retained cluster: the representative read, the polished
consensus, the cluster size, and the annotated combined-CDR region. It is the
same sequence of steps used for the manuscript results, and it is
library-agnostic: the scaffold database can be a proprietary or a germline
family database.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field

import pandas as pd

from . import annotator
from .annotator import DEFAULT_MIN_VOTES
from .cluster import (
    DEFAULT_CUTOFF,
    DEFAULT_METRIC,
    DEFAULT_MIN_SIZE,
    DEFAULT_ROI,
    cluster_by_scaffold_and_roi,
    filter_and_rank_clusters,
)
from .consensus import DEFAULT_RACON_ERROR, build_cluster_consensus
from .output_table import build_output_table
from .picking import is_functional, pick_leads, select_functional
from .preprocess import (
    CUTADAPT_ERROR,
    MEAN_Q,
    MIN_LEN_SCFV,
    MIN_LEN_VHH,
    PRETRIM_MIN_LEN,
    filter_length_quality,
    trim_adapters,
)


def _sample_name(path: str) -> str:
    """Derive a sample name from a FASTQ path (strips .gz, .fastq, .fq)."""
    base = os.path.basename(path)
    for ext in (".gz", ".fastq", ".fq"):
        if base.endswith(ext):
            base = base[: -len(ext)]
    return base


@dataclass
class PipelineConfig:
    """Parameters for a run. Defaults match the manuscript settings."""

    scaffold_db: str
    single_chain: bool = True
    roi: str = DEFAULT_ROI
    cutoff: int = DEFAULT_CUTOFF
    metric: str = DEFAULT_METRIC
    min_size: int = DEFAULT_MIN_SIZE
    # Length gate; when None it is set by format (scFv 650 bp, VHH 300 bp).
    min_length: int | None = None
    mean_q: int = MEAN_Q
    cutadapt_error: int = CUTADAPT_ERROR
    racon_error: float = DEFAULT_RACON_ERROR
    min_votes: int = DEFAULT_MIN_VOTES
    # Column carrying scaffold identity used to group reads before clustering.
    scaffold_col: str = "match_name_1_2"
    cpus: int = 1
    trim: bool = True
    adapters: dict[str, list[str]] = field(default_factory=dict)
    pick: bool = False
    top_n: int | None = 10
    output_table: bool = False

    @property
    def effective_min_length(self) -> int:
        if self.min_length is not None:
            return self.min_length
        return MIN_LEN_VHH if self.single_chain else MIN_LEN_SCFV

    @property
    def effective_scaffold_col(self) -> str:
        # Single-domain annotation carries scaffold identity in match_name_1.
        if self.scaffold_col != "match_name_1_2":
            return self.scaffold_col
        return "match_name_1" if self.single_chain else "match_name_1_2"


@dataclass
class PipelineResult:
    """Outputs of a run."""

    picks: pd.DataFrame
    annotated: pd.DataFrame
    stage_counts: dict[str, int]
    leads: pd.DataFrame | None = None
    output_table: pd.DataFrame | None = None


def _reannotate_consensus(
    picks: pd.DataFrame, workdir: str, config: "PipelineConfig"
) -> pd.DataFrame:
    """Re-annotate the consensus sequences and carry the cluster size across.

    Mirrors the manuscript's "a final annotation is performed on the consensus
    sequences" step. Returns the annotated consensus (one row per consensus),
    from which leads and the operator-facing output table are derived.
    """
    valid = picks[picks["consensus_seq"].notna() & (picks["consensus_seq"].astype(str).str.len() > 0)]
    if valid.empty:
        return pd.DataFrame()

    cons_dir = os.path.join(workdir, "consensus_annotate_in")
    cons_out = os.path.join(workdir, "consensus_annotate_out")
    os.makedirs(cons_dir, exist_ok=True)
    fasta = os.path.join(cons_dir, "consensus.fasta")
    size_by_seq: dict[str, int] = {}
    with open(fasta, "w") as fh:
        for _, row in valid.iterrows():
            seq = str(row["consensus_seq"]).upper()
            fh.write(f">{row['consensus_id']}\n{seq}\n")
            size_by_seq[seq] = int(row["cluster_size"])

    annotated = annotator.annotate(
        input_folder=cons_dir,
        output_folder=cons_out,
        scaffold_db=config.scaffold_db,
        single_chain=config.single_chain,
        cpus=config.cpus,
        min_votes=config.min_votes,
    )
    if annotated.empty:
        return pd.DataFrame()
    annotated["cluster_size"] = annotated["read"].astype(str).str.upper().map(size_by_seq)
    annotated = annotated.dropna(subset=["cluster_size"])
    annotated["cluster_size"] = annotated["cluster_size"].astype(int)
    annotated["count"] = annotated["cluster_size"]
    return annotated


def _write_read_accounting(workdir: str, stage_counts: dict[str, int]) -> str:
    """Write ``read_accounting.csv``: reads surviving each pipeline stage.

    Tabulates the flow raw -> quality-filtered -> annotated -> functional
    (no stop-codon) -> clustered (>= minimum cluster size), with a convenience
    percent-nonfunctional column. Emitted on every run so the read accounting is
    a reproducible pipeline output rather than only printed to the console.
    """
    ann = stage_counts.get("annotated_reads", 0)
    func = stage_counts.get("not_stop_codon_reads", 0)
    row = {
        "raw_reads": stage_counts.get("raw_reads", 0),
        "quality_filtered_reads": stage_counts.get("filtered_reads", 0),
        "annotated_reads": ann,
        "functional_reads": func,
        "pct_nonfunctional": round(100.0 * (ann - func) / ann, 2) if ann else 0.0,
        "clusters_retained": stage_counts.get("clusters_retained", 0),
        "clustered_reads": stage_counts.get("clustered_reads", 0),
        "output_rows": stage_counts.get(
            "output_table_rows", stage_counts.get("leads", 0)
        ),
    }
    path = os.path.join(workdir, "read_accounting.csv")
    pd.DataFrame([row]).to_csv(path, index=False)
    return path


def run(fastqs: list[str], workdir: str, config: PipelineConfig) -> PipelineResult:
    """Run the full pipeline over one or more FASTQ files.

    Parameters
    ----------
    fastqs:
        Input FASTQ paths (one sample; multiple files are concatenated).
    workdir:
        Working directory for intermediates and per-cluster consensus files.
    config:
        Run parameters.
    """
    os.makedirs(workdir, exist_ok=True)
    ann_in = os.path.join(workdir, "annotate_in")
    ann_out = os.path.join(workdir, "annotate_out")
    os.makedirs(ann_in, exist_ok=True)

    stage_counts: dict[str, int] = {}
    prepared: list[str] = []
    for i, fq in enumerate(fastqs):
        current = fq
        if config.trim:
            sample = _sample_name(fq)
            adapters = config.adapters.get(sample)
            if adapters:
                trimmed = os.path.join(workdir, f"{sample}.trimmed.fastq")
                current = trim_adapters(
                    current, trimmed, adapters, PRETRIM_MIN_LEN,
                    config.cpus, config.cutadapt_error,
                )
        filtered = os.path.join(ann_in, f"sample_{i}.filtered.fastq")
        kept, total = filter_length_quality(
            current,
            filtered,
            config.effective_min_length,
            config.mean_q,
        )
        stage_counts["raw_reads"] = stage_counts.get("raw_reads", 0) + total
        stage_counts["filtered_reads"] = stage_counts.get("filtered_reads", 0) + kept
        prepared.append(filtered)

    annotated = annotator.annotate(
        input_folder=ann_in,
        output_folder=ann_out,
        scaffold_db=config.scaffold_db,
        single_chain=config.single_chain,
        cpus=config.cpus,
        min_votes=config.min_votes,
    )
    stage_counts["annotated_reads"] = len(annotated)

    if "clone_id" not in annotated.columns or annotated["clone_id"].isna().all() or (
        annotated["clone_id"].astype(str).str.len() == 0
    ).all():
        annotated = annotated.reset_index(drop=True)
        annotated["clone_id"] = [f"read_{i:06d}" for i in range(len(annotated))]

    # Remove reads carrying a stop-codon before clustering (Methods: "remove those
    # with stop-codons not resulting from frameshifts").
    issue = annotated.get("sequence_issue", pd.Series("", index=annotated.index)).fillna("")
    to_cluster = annotated[~issue.astype(str).str.contains("Stop codon")]
    stage_counts["not_stop_codon_reads"] = len(to_cluster)

    clustered = cluster_by_scaffold_and_roi(
        to_cluster,
        scaffold_col=config.effective_scaffold_col,
        roi_col=config.roi,
        cutoff=config.cutoff,
        metric=config.metric,
    )
    clusters = filter_and_rank_clusters(clustered, min_size=config.min_size)
    stage_counts["clusters_retained"] = len(clusters)
    stage_counts["clustered_reads"] = sum(c.size for c in clusters)

    rows = []
    for view in clusters:
        # Reference the most-abundant functional read as the seed (manuscript),
        # while keeping every read in the cluster for the racon pileup.
        members = view.members
        func_mask = members.apply(lambda r: is_functional(r, config.single_chain), axis=1)
        ordered = pd.concat([members[func_mask], members[~func_mask]]).reset_index(drop=True)
        cons = build_cluster_consensus(
            ordered,
            workdir=os.path.join(workdir, "consensus", f"cluster_{view.rank:04d}"),
            racon_error=config.racon_error,
        )
        seed = ordered.iloc[0]
        rows.append(
            {
                "cluster_rank": view.rank,
                "cluster_size": view.size,
                "consensus_id": cons.consensus_id,
                "representative_read": cons.representative_read,
                "consensus_seq": cons.consensus_seq,
                config.roi: seed.get(config.roi, ""),
                "Merged_CDRs_AA": seed.get("Merged_CDRs_AA", ""),
                "match_name": seed.get("match_name_1", seed.get("match_name_2", "")),
            }
        )
    picks = pd.DataFrame(rows)

    leads = None
    output_table = None
    if (config.pick or config.output_table) and not picks.empty:
        reannotated = _reannotate_consensus(picks, workdir, config)
        if not reannotated.empty:
            if config.pick:
                leads = pick_leads(reannotated, config.single_chain, top_n=config.top_n)
                stage_counts["leads"] = len(leads)
            if config.output_table:
                functional = select_functional(reannotated, config.single_chain)
                output_table = build_output_table(functional, config.single_chain)
                stage_counts["output_table_rows"] = len(output_table)

    # Read accounting: tabulate how many reads survive each stage. Written as a
    # first-class pipeline output for every run.
    _write_read_accounting(workdir, stage_counts)

    return PipelineResult(
        picks=picks, annotated=annotated, stage_counts=stage_counts,
        leads=leads, output_table=output_table,
    )
