#!/usr/bin/env python3
# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""Compare the ONT pipeline output for Target 1 against its reference.

The reference (data/Sanger/target1_reference.csv) is the annotated comparison
set; the validation shows the ONT cluster picks match it, reproducing the paper's
ONT-side numbers (9/10 shared, 7/10 concordant picks).

Method: combined-CDR nucleotide clustering, Levenshtein LD<=3, minimum 5 reads;
top pick = most-abundant full-length AA per HCDR3, functional-only.
"""

from __future__ import annotations

import glob
import json
import os

import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = os.path.join(REPO, "data", "runs", "paper_ont")
SANGER = os.path.join(REPO, "data", "Sanger")
MIN_READS = int(os.environ.get("RADIANT_MIN_READS", "5"))
AA = "sequence_aa_1_2"

# One reproduction target (Target 1) plus the germline mock benchmark.
SAMPLE = "target1"
PAPER_SHARED, PAPER_PICK = "9/10", "7/10"
METHOD_LABEL = "combined-CDR nucleotide, Levenshtein LD<=3, >= 5 reads"


def functional(d: pd.DataFrame) -> pd.DataFrame:
    return d[
        (d["functional_1"].astype(str) == "functional")
        & (d["functional_2"].astype(str) == "functional")
        & (d.get("sequence_issue", pd.Series("", index=d.index)).fillna("") == "")
    ]


def cluster_sizes(picks_csv: str) -> dict[str, int]:
    picks = pd.read_csv(picks_csv)
    size: dict[str, int] = {}
    for _, r in picks.iterrows():
        for col in ("consensus_seq", "representative_read"):
            if col in picks.columns and isinstance(r.get(col), str):
                size[str(r[col]).upper()] = int(r["cluster_size"])
    return size


def main() -> None:
    run = os.path.join(RUNS, SAMPLE)

    # ONT side: the pipeline's per-cluster consensus (combined-CDR nucleotide
    # clustering, LD<=3), filtered to clusters of at least MIN_READS members.
    cons_csv = glob.glob(os.path.join(run, "consensus_annotate_out", "*output.csv"))[0]
    cons = functional(pd.read_csv(cons_csv, low_memory=False)).copy()
    size = cluster_sizes(os.path.join(RUNS, f"{SAMPLE}_picks.csv"))
    cons["sz"] = cons["read"].astype(str).str.upper().map(size).fillna(1).astype(int)
    cons = cons[cons["sz"] >= MIN_READS]
    ont_top = (
        cons.sort_values("sz", ascending=False, kind="stable")
        .groupby("cdr3_aa_2", sort=False)[AA].first().to_dict()
    )
    ont_h3 = set(cons["cdr3_aa_2"].astype(str))

    # Reference side (not re-annotated here). Top pick = most-abundant full-length
    # per HCDR3, ties broken by earliest read order.
    ref = pd.read_csv(os.path.join(SANGER, "target1_reference.csv"))
    ab = (
        ref.groupby(["cdr3_aa_2", AA], sort=False)
        .agg(n=("read_order", "size"), first_order=("read_order", "min"))
        .reset_index()
    )
    s_pick = (
        ab.sort_values(["n", "first_order"], ascending=[False, True], kind="stable")
        .groupby("cdr3_aa_2", sort=False)[AA].first().to_dict()
    )
    clusters = list(s_pick)

    shared = sum(1 for h in clusters if h in ont_h3)
    pick = sum(1 for h in clusters if ont_top.get(h) == s_pick[h])
    n = len(clusters)

    print("\nRADiANT ONT pipeline validation, Target 1")
    print(f"(clustering: {METHOD_LABEL})\n")
    print(f"  HCDR3 clusters (reference):          {n}")
    print(f"  Shared with ONT:                     {shared}/{n}   (paper: {PAPER_SHARED})")
    print(f"  Concordant top cluster picks:        {pick}/{n}   (paper: {PAPER_PICK})")
    ok = (f"{shared}/{n}" == PAPER_SHARED) and (f"{pick}/{n}" == PAPER_PICK)
    print(f"\n  {'MATCH, reproduces the paper ONT numbers.' if ok else 'Deviation from the paper numbers.'}")

    with open(os.path.join(REPO, "validation", "results.json"), "w") as fh:
        json.dump({"target": "Target 1", "shared": f"{shared}/{n}", "pick": f"{pick}/{n}",
                   "paper_shared": PAPER_SHARED, "paper_pick": PAPER_PICK, "match": ok}, fh, indent=2)


if __name__ == "__main__":
    main()
