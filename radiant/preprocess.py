# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""Adapter trimming and length/quality filtering of ONT reads.

Trimming uses ``cutadapt`` with the linked, anchored adapters supplied per sample,
in demultiplexing mode: only reads that carry an adapter are kept (reads with no
adapter are discarded). Quality/length filtering uses ``fastq-filter`` to keep
reads at least ``min_length`` bp long with a mean quality of at least ``mean_q``.
"""

from __future__ import annotations

import os
import shutil
import subprocess

# Length gates by format (scFv 650 bp, VHH 300 bp), applied at the quality-filter step.
MIN_LEN_SCFV = 650
MIN_LEN_VHH = 300
# cutadapt minimum length at the adapter-trim step (before quality filtering).
PRETRIM_MIN_LEN = 300
# Quality gate: mean read quality >= MEAN_Q (fastq-filter -Q).
MEAN_Q = 30
CUTADAPT_ERROR = 4  # cutadapt -e (allowed errors)


def _count_reads(path: str) -> int:
    if not os.path.exists(path):
        return 0
    with open(path) as fh:
        return sum(1 for _ in fh) // 4


def trim_adapters(
    in_fastq: str,
    out_fastq: str,
    adapters: list[str],
    min_length: int = PRETRIM_MIN_LEN,
    cpus: int = 2,
    cutadapt_error: int = CUTADAPT_ERROR,
) -> str:
    """Trim linked adapters with cutadapt, keeping only adapter-matched reads.

    ``adapters`` is ``[forward, reverse]``. cutadapt is run in demultiplexing mode
    (``{name}`` in the output): reads matching the forward adapter go to ``-1``,
    the rest to ``-unknown``; the reverse adapter is then applied to ``-unknown``.
    The kept output is the union of the two adapter-matched (``-1``) files; reads
    with no adapter are discarded. This is the reference two-pass trimming.
    """
    cutadapt = shutil.which("cutadapt")
    if cutadapt is None:
        raise FileNotFoundError("cutadapt not found. See env/setup_tools.sh.")
    forward = adapters[0]
    reverse = adapters[1] if len(adapters) > 1 else adapters[0]
    wd = os.path.dirname(out_fastq) or "."

    fwd_tmpl = os.path.join(wd, "_pretrim_fwd-{name}.fastq")
    subprocess.run(
        [cutadapt, "--minimum-length", str(min_length), "-a", forward,
         "-e", str(cutadapt_error), "--cores", str(cpus), "-o", fwd_tmpl, in_fastq],
        check=True, capture_output=True, text=True,
    )
    fwd_matched = os.path.join(wd, "_pretrim_fwd-1.fastq")
    fwd_unknown = os.path.join(wd, "_pretrim_fwd-unknown.fastq")

    rev_tmpl = os.path.join(wd, "_pretrim_rev-{name}.fastq")
    if os.path.exists(fwd_unknown):
        subprocess.run(
            [cutadapt, "--minimum-length", str(min_length), "-a", reverse,
             "-e", str(cutadapt_error), "--cores", str(cpus), "-o", rev_tmpl, fwd_unknown],
            check=True, capture_output=True, text=True,
        )
    rev_matched = os.path.join(wd, "_pretrim_rev-1.fastq")

    with open(out_fastq, "w") as out:
        for part in (fwd_matched, rev_matched):
            if os.path.exists(part):
                with open(part) as fh:
                    shutil.copyfileobj(fh, out)
    return out_fastq


def filter_length_quality(
    in_fastq: str,
    out_fastq: str,
    min_length: int = MIN_LEN_SCFV,
    mean_q: int = MEAN_Q,
) -> tuple[int, int]:
    """Keep reads at least ``min_length`` bp with mean quality >= ``mean_q``.

    Uses ``fastq-filter`` (``-l`` minimum length, ``-Q`` minimum mean quality),
    matching the reference pipeline.

    Returns
    -------
    (kept, total)
        Number of reads written and number read.
    """
    fastq_filter = shutil.which("fastq-filter")
    if fastq_filter is None:
        raise FileNotFoundError("fastq-filter not found. See env/setup_tools.sh.")
    total = _count_reads(in_fastq)
    subprocess.run(
        [fastq_filter, "-Q", str(mean_q), "-l", str(min_length), "-o", out_fastq, in_fastq],
        check=True, capture_output=True, text=True,
    )
    return _count_reads(out_fastq), total
