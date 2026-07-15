#!/usr/bin/env python3
# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""Run the RADiANT pipeline on the manuscript ONT campaign samples.

Reads the per-sample linked adapters from an ``adapter_seqs`` CSV, then trims,
filters, annotates (against the scaffold database), clusters, and builds a
consensus for each sample. Writes a picks CSV per sample and prints per-stage
read counts and the top cluster picks.

Example
-------
    RADIANT_ANNOTATOR_BACKEND=radiant_annotator_backend \
    RACON_PATH=$HOME/.local/bin/racon MINIMAP2_PATH=$(which minimap2) \
    python campaigns/paper_ont/run_ont_only.py \
        --data-dir /path/to/reads \
        --adapter-csv /path/to/adapter_seqs.csv \
        --db /path/to/scaffold_db.txt \
        --samples target1 \
        --out data/runs/paper_ont
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from radiant.pipeline import PipelineConfig, run


def load_adapters(path: str) -> dict[str, list[str]]:
    """Return ``{sample: [forward_adapter, reverse_adapter]}`` from the CSV."""
    adapters: dict[str, list[str]] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            adapters[row["sample"]] = [row["for_adapter"], row["rev_adapter"]]
    return adapters


def find_fastq(data_dir: str, sample: str) -> str:
    for name in os.listdir(data_dir):
        if name.startswith(sample) and (name.endswith(".fastq.gz") or name.endswith(".fastq")):
            return os.path.join(data_dir, name)
    raise FileNotFoundError(f"No FASTQ for sample {sample!r} in {data_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--adapter-csv", required=True)
    parser.add_argument("--db", required=True)
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--single-chain", action="store_true", help="VHH/single-domain.")
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--pick", action="store_true",
                        help="Run final lead selection (functional + 100%% HCDR3 dedup + rank).")
    parser.add_argument("--top-n", type=int, default=None,
                        help="Keep only the top N leads (default: keep all ranked).")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    adapters = load_adapters(args.adapter_csv)
    os.makedirs(args.out, exist_ok=True)

    for sample in args.samples:
        fastq = find_fastq(args.data_dir, sample)
        print("=" * 70)
        print(f"Sample: {sample}")
        print(f"FASTQ: {fastq}")
        config = PipelineConfig(
            scaffold_db=args.db,
            single_chain=args.single_chain,
            trim=True,
            adapters=adapters,
            cpus=args.cpus,
            pick=args.pick,
            top_n=args.top_n,
        )
        result = run([fastq], os.path.join(args.out, sample), config)

        print("Per-stage read counts:")
        for stage, n in result.stage_counts.items():
            print(f"  {stage}: {n:,}")

        picks_path = os.path.join(args.out, f"{sample}_picks.csv")
        result.picks.to_csv(picks_path, index=False)
        print(f"Picks -> {picks_path}")

        if result.leads is not None:
            leads_path = os.path.join(args.out, f"{sample}_leads.csv")
            result.leads.to_csv(leads_path, index=False)
            print(f"Leads (functional, 100%% HCDR3 dedup, ranked): {len(result.leads)} -> {leads_path}")

        cols = [c for c in ("cluster_rank", "cluster_size", "match_name", "Merged_CDRs_AA")
                if c in result.picks.columns]
        if len(result.picks):
            print(f"Top {args.top} cluster picks:")
            print(result.picks[cols].head(args.top).to_string(index=False))
        print()


if __name__ == "__main__":
    main()
