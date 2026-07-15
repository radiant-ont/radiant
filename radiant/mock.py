# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""Simulate an Oxford Nanopore read set from germline-derived antibody clones.

Full-length heavy-chain sequences are assembled from germline frameworks and
OAS CDRs, then perturbed with an error model that mimics the R10.4.1 profile
(substitutions plus indels enriched in homopolymer runs). Because every read
derives from a known clone, the input clones form an exact truth set against
which the workflow's recovered consensus sequences can be scored.

The simulation is fully deterministic given a seed, so a run is reproducible.
"""

from __future__ import annotations

import random
from dataclasses import dataclass

from .codon import back_translate
from .germline import GermlineFrameworks
from .oas import OasClone

_BASES = ("A", "C", "G", "T")


@dataclass(frozen=True)
class SimClone:
    """A simulated clone: identity, true sequence, and intended read depth."""

    clone_id: str
    v_gene: str
    depth: int
    true_nt: str
    true_aa_cdrs: tuple[str, str, str]


@dataclass(frozen=True)
class ErrorModel:
    """Per-base error probabilities for the simulated ONT profile."""

    substitution: float = 0.008
    deletion: float = 0.003
    insertion: float = 0.002
    homopolymer_multiplier: float = 6.0
    homopolymer_min_run: int = 3


def assemble_heavy_chain(clone: OasClone, fw: GermlineFrameworks) -> str:
    """Assemble a full-length heavy-chain nucleotide sequence for a clone."""
    return "".join(
        [
            fw.fr1,
            back_translate(clone.cdr1_aa),
            fw.fr2,
            back_translate(clone.cdr2_aa),
            fw.fr3,
            back_translate(clone.cdr3_aa),
            fw.fr4,
        ]
    )


def select_clones(
    clones: list[OasClone],
    frameworks: dict[str, GermlineFrameworks],
    n_clones: int,
) -> list[OasClone]:
    """Pick distinct, buildable clones deterministically (most abundant first).

    Clones are de-duplicated at the CDR level and restricted to genes present
    in the framework reference, then ranked by redundancy so the selection is
    stable across runs.
    """
    seen: set[tuple[str, str, str]] = set()
    ranked = sorted(
        clones,
        key=lambda c: (-c.redundancy, c.v_gene, c.cdr3_aa),
    )
    chosen: list[OasClone] = []
    for clone in ranked:
        if clone.v_gene not in frameworks or clone.cdr_key in seen:
            continue
        seen.add(clone.cdr_key)
        chosen.append(clone)
        if len(chosen) >= n_clones:
            break
    return chosen


def _depth_for_rank(rank: int, n_clones: int) -> int:
    """Assign a read depth that spans the abundance spectrum.

    Most clones sit above the minimum cluster size so they are recoverable; a
    tail sits below it so filtering behaviour can also be demonstrated.
    """
    if rank < n_clones * 0.15:
        return 30
    if rank < n_clones * 0.6:
        return 12
    if rank < n_clones * 0.85:
        return 6
    return 3


def build_sim_clones(
    clones: list[OasClone],
    frameworks: dict[str, GermlineFrameworks],
    n_clones: int,
) -> list[SimClone]:
    """Turn selected OAS clones into simulated clones with true sequences."""
    selected = select_clones(clones, frameworks, n_clones)
    sim: list[SimClone] = []
    for rank, clone in enumerate(selected):
        fw = frameworks[clone.v_gene]
        sim.append(
            SimClone(
                clone_id=f"clone_{rank:03d}_{clone.v_gene}",
                v_gene=clone.v_gene,
                depth=_depth_for_rank(rank, len(selected)),
                true_nt=assemble_heavy_chain(clone, fw),
                true_aa_cdrs=(clone.cdr1_aa, clone.cdr2_aa, clone.cdr3_aa),
            )
        )
    return sim


def _homopolymer_run_length(seq: str, index: int) -> int:
    """Length of the homopolymer run ending at ``index``."""
    run = 1
    j = index - 1
    while j >= 0 and seq[j] == seq[index]:
        run += 1
        j -= 1
    return run


def apply_errors(
    seq: str,
    model: ErrorModel,
    rng: random.Random,
    base_quality: int = 40,
    error_quality: int = 10,
) -> tuple[str, list[int]]:
    """Perturb ``seq`` with the ONT-like error model.

    Returns the perturbed sequence together with a per-output-base Phred
    quality list. Quality is assigned as each base is emitted (low at
    substituted and inserted bases, high otherwise), so it stays aligned to the
    output even across indels rather than being inferred by position afterwards.
    """
    out: list[str] = []
    quals: list[int] = []
    for i, base in enumerate(seq):
        run = _homopolymer_run_length(seq, i)
        in_homopolymer = run >= model.homopolymer_min_run
        boost = model.homopolymer_multiplier if in_homopolymer else 1.0

        if rng.random() < model.deletion * boost:
            continue  # drop this base
        if rng.random() < model.substitution:
            out.append(rng.choice([b for b in _BASES if b != base]))
            quals.append(error_quality)
        else:
            out.append(base)
            quals.append(base_quality)
        if rng.random() < model.insertion * boost:
            out.append(base)  # duplicate (homopolymer-style insertion)
            quals.append(error_quality)
    return "".join(out), quals


def simulate_reads(
    sim_clones: list[SimClone],
    model: ErrorModel,
    seed: int,
    base_quality: int = 40,
    error_quality: int = 10,
):
    """Yield ``(read_id, sequence, quality)`` tuples for all simulated reads.

    Quality strings carry a high baseline with lower values at mismatched
    positions, so mean-quality filters behave sensibly while the read still
    reflects the introduced errors.
    """
    rng = random.Random(seed)
    for clone in sim_clones:
        for replicate in range(clone.depth):
            read, quals = apply_errors(
                clone.true_nt, model, rng, base_quality, error_quality
            )
            qual_str = "".join(chr(min(q, 93) + 33) for q in quals)
            read_id = f"{clone.clone_id}_read{replicate:03d}"
            yield read_id, read, qual_str
