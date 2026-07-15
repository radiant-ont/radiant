# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""Command-line entry point for the RADiANT ONT pipeline.

Every stage parameter is an explicit flag with the manuscript default, so a run
is fully specified by its command line. Running with different flags will give
different results; the exact manuscript settings are the defaults printed at the
top of every run and by --help.

Usage:
    python -m radiant --scaffold-db DB.txt --workdir OUT reads.fastq [reads2.fastq ...]
"""

from __future__ import annotations

import argparse
import json
import os

from .cluster import DEFAULT_CUTOFF, DEFAULT_METRIC, DEFAULT_MIN_SIZE, DEFAULT_ROI, METRICS
from .consensus import DEFAULT_RACON_ERROR
from .annotator import DEFAULT_MIN_VOTES
from .preprocess import CUTADAPT_ERROR, MEAN_Q, MIN_LEN_SCFV, MIN_LEN_VHH
from .pipeline import PipelineConfig, run


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="radiant",
        description="RADiANT ONT pipeline: trim, filter, annotate, cluster, "
        "consensus, pick. Defaults reproduce the manuscript.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("fastqs", nargs="+", help="input FASTQ file(s) for one sample")
    p.add_argument("--scaffold-db", required=True, help="scaffold/germline database")
    p.add_argument("--workdir", required=True, help="output/working directory")

    fmt = p.add_argument_group("library format")
    fmt.add_argument("--format", choices=("scfv", "vhh"), default="scfv",
                     help="chain format (scfv = paired VL+VH, vhh = single domain)")

    filt = p.add_argument_group("trim / filter")
    filt.add_argument("--no-trim", action="store_true", help="skip adapter trimming")
    filt.add_argument("--cutadapt-error", type=int, default=CUTADAPT_ERROR,
                      help="cutadapt allowed errors (-e)")
    filt.add_argument("--min-length", type=int, default=None,
                      help=f"minimum read length (default {MIN_LEN_SCFV} scfv / {MIN_LEN_VHH} vhh)")
    filt.add_argument("--mean-quality", type=int, default=MEAN_Q,
                      help="minimum mean read quality (fastq-filter -Q)")

    ann = p.add_argument_group("annotate")
    ann.add_argument("--min-votes", type=int, default=DEFAULT_MIN_VOTES,
                     help="annotator confidence threshold")

    clu = p.add_argument_group("cluster")
    clu.add_argument("--region", default=DEFAULT_ROI,
                     help="column clustered on (e.g. Merged_CDRs_NUC, Merged_CDRs_AA, "
                     "sequence_1_2, sequence_aa_1_2)")
    clu.add_argument("--distance-metric", choices=METRICS, default=DEFAULT_METRIC,
                     help="sequence distance metric for clustering")
    clu.add_argument("--distance", type=int, default=DEFAULT_CUTOFF,
                     help="maximum distance to the cluster seed")
    clu.add_argument("--min-cluster-size", type=int, default=DEFAULT_MIN_SIZE,
                     help="drop clusters with fewer than this many reads")
    clu.add_argument("--scaffold-col", default="match_name_1_2",
                     help="column grouping reads by scaffold before clustering")

    con = p.add_argument_group("consensus")
    con.add_argument("--racon-error", type=float, default=DEFAULT_RACON_ERROR,
                     help="racon error threshold (-e)")

    out = p.add_argument_group("output")
    out.add_argument("--pick", action="store_true",
                     help="re-annotate consensus and select top leads per HCDR3")
    out.add_argument("--top-n", type=int, default=10, help="leads to keep when --pick")
    out.add_argument("--output-table", action="store_true",
                     help="write the supplementary-format output table (one row per "
                     "cluster: sequences, frameworks/CDRs, votes, liabilities, "
                     "biophysical properties, confidence score)")
    out.add_argument("--cpus", type=int, default=1, help="worker processes")
    return p


def config_from_args(args: argparse.Namespace) -> PipelineConfig:
    return PipelineConfig(
        scaffold_db=args.scaffold_db,
        single_chain=(args.format == "vhh"),
        roi=args.region,
        cutoff=args.distance,
        metric=args.distance_metric,
        min_size=args.min_cluster_size,
        min_length=args.min_length,
        mean_q=args.mean_quality,
        cutadapt_error=args.cutadapt_error,
        racon_error=args.racon_error,
        min_votes=args.min_votes,
        scaffold_col=args.scaffold_col,
        cpus=args.cpus,
        trim=not args.no_trim,
        pick=args.pick,
        top_n=args.top_n,
        output_table=args.output_table,
    )


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    config = config_from_args(args)

    resolved = {
        "format": args.format,
        "trim": not args.no_trim,
        "cutadapt_error": config.cutadapt_error,
        "min_length": config.effective_min_length,
        "mean_quality": config.mean_q,
        "min_votes": config.min_votes,
        "region": config.roi,
        "distance_metric": config.metric,
        "distance": config.cutoff,
        "min_cluster_size": config.min_size,
        "scaffold_col": config.effective_scaffold_col,
        "racon_error": config.racon_error,
    }
    print("== RADiANT parameters (as run) ==")
    for k, v in resolved.items():
        print(f"  {k}: {v}")
    print()

    os.makedirs(args.workdir, exist_ok=True)
    result = run(args.fastqs, args.workdir, config)

    picks_path = os.path.join(args.workdir, "picks.csv")
    result.picks.to_csv(picks_path, index=False)
    with open(os.path.join(args.workdir, "run_parameters.json"), "w") as fh:
        json.dump(resolved, fh, indent=2)
    if result.leads is not None:
        result.leads.to_csv(os.path.join(args.workdir, "leads.csv"), index=False)
    if result.output_table is not None:
        table_path = os.path.join(args.workdir, "output_table.csv")
        result.output_table.to_csv(table_path, index=False)
        print(f"output table written to {table_path}")

    print("== stage counts ==")
    for k, v in result.stage_counts.items():
        print(f"  {k}: {v}")
    print(f"\nclusters written to {picks_path}")
    print(f"parameters recorded in {os.path.join(args.workdir, 'run_parameters.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
