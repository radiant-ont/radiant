#!/usr/bin/env python3
# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""Build a germline-family scaffold database from an OAS subset.

Example
-------
    python scaffolds/build_scaffold_db.py \
        --oas /path/to/oas_human_db_subset_1M.csv \
        --germline-ref scaffolds/reference/germline_frameworks.txt \
        --out data/generated/germline_scaffold_db.txt
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from radiant.germline import load_germline_frameworks
from radiant.oas import load_oas_clones
from radiant.scaffold_db import build_germline_scaffold_db, write_scaffold_db


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--oas", required=True, help="OAS-style CSV subset.")
    parser.add_argument(
        "--germline-ref",
        required=True,
        help="Space-separated germline framework reference.",
    )
    parser.add_argument("--out", required=True, help="Output scaffold DB path.")
    args = parser.parse_args()

    clones = load_oas_clones(args.oas)
    frameworks = load_germline_frameworks(args.germline_ref)
    entries = build_germline_scaffold_db(clones, frameworks)
    write_scaffold_db(entries, args.out)

    genes = sorted({e[0].removeprefix("germline_") for e in entries})
    print(f"Loaded {len(clones):,} OAS clones")
    print(f"Loaded {len(frameworks):,} germline framework genes")
    print(f"Wrote {len(entries)} scaffold entries to {args.out}")
    print(f"Genes: {', '.join(genes)}")


if __name__ == "__main__":
    main()
