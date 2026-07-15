# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""radiant-annotator, run-only interface to the antibody annotator.

This package exposes a single ``annotate_folder`` entry point used by the RADiANT
pipeline (``radiant.annotator``). The annotation logic is provided as a
**compiled, run-only binary** (``bin/annotator``); no source is distributed. The
binary reads FASTA/FASTQ from an input folder and writes an annotation CSV with
framework/CDR fields, including the combined-CDR nucleotide column used for
clustering.

Binary resolution order:
  1. ``$RADIANT_ANNOTATOR_BIN`` if set and executable.
  2. ``annotator`` on ``$PATH``.
  3. the ``bin/annotator`` bundled inside this package.

Build the binary from source with ``packaging/radiant_annotator/build_annotator.sh``.
"""

from __future__ import annotations

import glob
import os
import shutil
import stat
import subprocess

__all__ = ["annotate_folder", "__version__"]
__version__ = "1.0.0"


def _resolve_binary() -> str:
    env = os.environ.get("RADIANT_ANNOTATOR_BIN")
    if env and os.path.exists(env) and os.access(env, os.X_OK):
        return env
    onpath = shutil.which("annotator")
    if onpath:
        return onpath
    bundled = os.path.join(os.path.dirname(__file__), "bin", "annotator")
    if os.path.exists(bundled):
        os.chmod(bundled, os.stat(bundled).st_mode | stat.S_IEXEC)
        return bundled
    raise RuntimeError(
        "annotator binary not found. Set RADIANT_ANNOTATOR_BIN, put 'annotator' "
        "on PATH, or build it with packaging/radiant_annotator/build_annotator.sh."
    )


def annotate_folder(
    input_folder: str,
    output_folder: str,
    scaffold_db: str,
    single_chain: bool = True,
    cpus: int = 1,
    min_votes: int = 150,
):
    """Annotate every read file in ``input_folder`` against ``scaffold_db``.

    Invokes the compiled annotator binary. ``min_votes`` sets the annotator
    confidence threshold. Returns ``(csv_path, output_folder)``.
    """
    os.makedirs(output_folder, exist_ok=True)
    binary = _resolve_binary()
    cmd = [
        binary,
        "--input", input_folder,
        "--output", output_folder,
        "--db", scaffold_db,
        "--cpus", str(cpus),
        "--min-votes", str(min_votes),
    ]
    if single_chain:
        cmd.append("--single-chain")
    subprocess.run(cmd, check=True)
    hits = glob.glob(os.path.join(output_folder, "*output.csv")) or glob.glob(
        os.path.join(output_folder, "*.csv")
    )
    if not hits:
        raise RuntimeError(f"annotator produced no CSV in {output_folder}")
    return hits[0]
