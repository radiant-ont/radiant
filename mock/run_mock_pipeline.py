#!/usr/bin/env python3
# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""Run the full RADiANT pipeline on the germline mock set and score recovery.

The mock reads carry no adapters, so trimming is disabled; every other stage
(filter, annotate, cluster, consensus) runs exactly as for a real campaign.
Recovered consensus sequences are scored against the simulated truth set by
matching each recovered cluster to its nearest true clone at the amino-acid
combined-CDR level.

Example
-------
    RADIANT_ANNOTATOR_BACKEND=radiant_annotator_backend \
    RACON_PATH=$HOME/.local/bin/racon \
    python mock/run_mock_pipeline.py \
        --fastq data/generated/germline_mock.fastq \
        --scaffold-db data/generated/germline_scaffold_db.txt \
        --truth data/generated/germline_mock_truth.csv \
        --workdir data/runs/germline_mock
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

import editdistance

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from radiant.pipeline import PipelineConfig, run


def _translate(nt: str) -> str:
    from Bio.Seq import Seq

    trimmed = nt[: len(nt) - (len(nt) % 3)]
    return str(Seq(trimmed).translate())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--scaffold-db", required=True)
    parser.add_argument("--truth", required=True)
    parser.add_argument("--workdir", required=True)
    args = parser.parse_args()

    config = PipelineConfig(scaffold_db=args.scaffold_db, single_chain=True, trim=False)
    result = run([args.fastq], args.workdir, config)

    print("Stage counts:")
    for stage, n in result.stage_counts.items():
        print(f"  {stage}: {n:,}")

    with open(args.truth) as fh:
        truth = list(csv.DictReader(fh))
    truth_above = [t for t in truth if int(t["depth"]) >= 5]
    truth_cdr3 = {t["cdr3_aa"]: t for t in truth}

    picks = result.picks
    recovered = 0
    exact_consensus = 0
    for _, pick in picks.iterrows():
        merged_aa = str(pick.get("Merged_CDRs_AA", ""))
        # nearest true clone by amino-acid combined-CDR distance
        best = min(
            truth,
            key=lambda t: editdistance.eval(
                merged_aa, t["cdr1_aa"] + t["cdr2_aa"] + t["cdr3_aa"]
            ),
        )
        best_cdr = best["cdr1_aa"] + best["cdr2_aa"] + best["cdr3_aa"]
        dist = editdistance.eval(merged_aa, best_cdr)
        if dist <= 1:
            recovered += 1
        cons = pick.get("consensus_seq")
        if cons and _translate(str(cons)).find(best["cdr3_aa"]) >= 0:
            exact_consensus += 1

    print()
    print(f"Clusters retained: {len(picks)}")
    print(f"True clones at depth >= 5 (recoverable): {len(truth_above)}/{len(truth)}")
    print(f"Clusters matching a true clone (combined-CDR AA dist <= 1): {recovered}/{len(picks)}")
    print(f"Consensus sequences whose translation contains the true HCDR3: {exact_consensus}/{len(picks)}")


if __name__ == "__main__":
    main()
