# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""Per-cluster consensus by minimap2 alignment and racon polishing.

For each cluster, member reads are aligned to the most-abundant read with
``minimap2 -a -N 0`` and the resulting primary alignments are polished with
``racon -e 0.7`` to produce a consensus sequence. The representative reported
for the cluster is the most-abundant real read (its identifier), not a
synthetic sequence; the polished consensus is retained separately as an
accuracy measure.
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
from dataclasses import dataclass

import pandas as pd
import pysam


def _resolve(tool: str, env_var: str) -> str:
    """Locate an external tool via env override, PATH, or common locations."""
    override = os.environ.get(env_var)
    if override and os.path.isfile(override) and os.access(override, os.X_OK):
        return override
    found = shutil.which(tool)
    if found:
        return found
    for candidate in (
        os.path.expanduser(f"~/.local/bin/{tool}"),
        f"/opt/homebrew/bin/{tool}",
        f"/usr/local/bin/{tool}",
    ):
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    raise FileNotFoundError(
        f"{tool} not found. Install it or set {env_var}. See env/setup_tools.sh."
    )


@dataclass(frozen=True)
class ConsensusResult:
    """Consensus for one cluster."""

    consensus_id: str
    consensus_seq: str | None
    representative_read: str
    n_reads: int


DEFAULT_RACON_ERROR = 0.7


def build_cluster_consensus(
    members: pd.DataFrame,
    workdir: str,
    read_col: str = "read",
    id_col: str = "clone_id",
    racon_error: float = DEFAULT_RACON_ERROR,
) -> ConsensusResult:
    """Align cluster members to the most-abundant read and polish with racon.

    Parameters
    ----------
    members:
        Cluster rows, ordered so the first row is the most-abundant read.
    workdir:
        Directory for per-cluster intermediate files.
    read_col, id_col:
        Column names for the nucleotide read and its identifier.

    Returns
    -------
    ConsensusResult
        The polished consensus (or ``None`` if racon produced no output), the
        representative read sequence, and the number of member reads.
    """
    minimap2 = _resolve("minimap2", "MINIMAP2_PATH")
    racon = _resolve("racon", "RACON_PATH")
    os.makedirs(workdir, exist_ok=True)

    seed = members.iloc[0]
    consensus_id = str(seed[id_col])
    reference_seq = str(seed[read_col])

    ref_path = os.path.join(workdir, "reference.fa")
    with open(ref_path, "w") as fh:
        fh.write(f">reference\n{reference_seq}\n")

    # Write one read per actual observation (expanded by ``count``) so racon
    # polishes over the true read depth, matching the reference pipeline, rather
    # than over one copy per distinct sequence.
    reads_path = os.path.join(workdir, "reads.fastq")
    with open(reads_path, "w") as fh:
        for i, (_, row) in enumerate(members.iterrows()):
            seq = str(row[read_col])
            qual = "I" * len(seq)
            copies = int(row["count"]) if "count" in row and pd.notna(row["count"]) else 1
            for j in range(max(copies, 1)):
                fh.write(f"@{row[id_col]}_{i}_{j}\n{seq}\n+\n{qual}\n")

    sam_gz = os.path.join(workdir, "aln.sam.gz")
    with gzip.open(sam_gz, "wb", compresslevel=3) as out:
        proc = subprocess.Popen(
            [minimap2, "-a", "-N", "0", ref_path, reads_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        assert proc.stdout is not None
        shutil.copyfileobj(proc.stdout, out)
        proc.stdout.close()
        proc.wait()

    primary_sam = os.path.join(workdir, "primary.sam")
    pysam.view("-h", "-F", "0x900", "-o", primary_sam, sam_gz, save_stdout=primary_sam)

    result = subprocess.run(
        [racon, reads_path, primary_sam, ref_path, "-e", str(racon_error)],
        capture_output=True,
        text=True,
        stdin=subprocess.DEVNULL,
    )
    consensus_seq: str | None
    lines = [ln for ln in result.stdout.splitlines() if ln.strip()]
    # racon emits FASTA: header line then sequence line(s). Take the sequence.
    if len(lines) >= 2 and lines[0].startswith(">"):
        consensus_seq = "".join(lines[1:])
    else:
        consensus_seq = None

    return ConsensusResult(
        consensus_id=consensus_id,
        consensus_seq=consensus_seq,
        representative_read=reference_seq,
        n_reads=len(members),
    )
