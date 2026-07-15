#!/usr/bin/env python3
# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""Generate a mock ONT FASTQ (plus truth set) from germline-derived clones.

Example
-------
    python mock/simulate_ont.py \
        --oas /path/to/oas_human_db_subset_1M.csv \
        --germline-ref scaffolds/reference/germline_frameworks.txt \
        --clones 48 --seed 1 \
        --out-fastq data/generated/germline_mock.fastq \
        --out-truth data/generated/germline_mock_truth.csv
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from radiant.germline import load_germline_frameworks
from radiant.mock import ErrorModel, build_sim_clones, simulate_reads
from radiant.oas import load_oas_clones


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--oas", required=True)
    parser.add_argument("--germline-ref", required=True)
    parser.add_argument("--clones", type=int, default=48)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--out-fastq", required=True)
    parser.add_argument("--out-truth", required=True)
    args = parser.parse_args()

    clones = load_oas_clones(args.oas)
    frameworks = load_germline_frameworks(args.germline_ref)
    sim_clones = build_sim_clones(clones, frameworks, args.clones)

    n_reads = 0
    with open(args.out_fastq, "w") as fq:
        for read_id, seq, qual in simulate_reads(sim_clones, ErrorModel(), args.seed):
            fq.write(f"@{read_id}\n{seq}\n+\n{qual}\n")
            n_reads += 1

    with open(args.out_truth, "w", newline="") as tf:
        writer = csv.writer(tf)
        writer.writerow(
            ["clone_id", "v_gene", "depth", "true_len_nt", "cdr1_aa", "cdr2_aa", "cdr3_aa"]
        )
        for c in sim_clones:
            writer.writerow(
                [c.clone_id, c.v_gene, c.depth, len(c.true_nt), *c.true_aa_cdrs]
            )

    above = sum(1 for c in sim_clones if c.depth >= 5)
    print(f"Simulated {len(sim_clones)} clones, {n_reads:,} reads -> {args.out_fastq}")
    print(f"Truth set -> {args.out_truth}")
    print(f"Clones at depth >= 5 (recoverable): {above}/{len(sim_clones)}")


if __name__ == "__main__":
    main()
