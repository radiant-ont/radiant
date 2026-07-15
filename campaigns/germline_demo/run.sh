#!/usr/bin/env bash
# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

# Germline-family demonstration: build a non-proprietary scaffold database from
# a public OAS subset, simulate a matched ONT read set, and run the identical
# RADiANT pipeline against it. Recovered consensus sequences are scored against
# the simulated truth set. This shows the workflow is not specific to any one
# antibody library.
#
# Prerequisites: run env/setup_tools.sh, and provide an annotator backend
# (the distributed component, or a development backend via
# RADIANT_ANNOTATOR_BACKEND).
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO}"

OAS="${OAS:-oas_human_db_subset.csv}"
REF="scaffolds/reference/germline_frameworks.txt"
GEN="data/generated"
mkdir -p "${GEN}" data/runs

# 1. Germline-family scaffold database (one entry per V gene).
python scaffolds/build_scaffold_db.py \
    --oas "${OAS}" --germline-ref "${REF}" \
    --out "${GEN}/germline_scaffold_db.txt"

# 2. Matched mock ONT read set with a known truth set (seed fixed for reproducibility).
python mock/simulate_ont.py \
    --oas "${OAS}" --germline-ref "${REF}" \
    --clones 48 --seed 1 \
    --out-fastq "${GEN}/germline_mock.fastq" \
    --out-truth "${GEN}/germline_mock_truth.csv"

# 3. Run the full pipeline and score recovery against the truth set.
python mock/run_mock_pipeline.py \
    --fastq "${GEN}/germline_mock.fastq" \
    --scaffold-db "${GEN}/germline_scaffold_db.txt" \
    --truth "${GEN}/germline_mock_truth.csv" \
    --workdir data/runs/germline_mock
