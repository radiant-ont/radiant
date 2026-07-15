# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""Unit tests for the open, annotator-independent parts of the workflow."""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from radiant.codon import back_translate
from radiant.germline import GermlineFrameworks, load_germline_frameworks
from radiant.mock import (
    ErrorModel,
    apply_errors,
    assemble_heavy_chain,
    build_sim_clones,
    simulate_reads,
)
from radiant.oas import OasClone
from radiant.scaffold_db import build_germline_scaffold_db


def test_back_translate_is_deterministic_and_correct_length():
    assert back_translate("AR") == "GCCCGG"
    assert len(back_translate("ARDYW")) == 15
    assert back_translate("ar") == back_translate("AR")  # case-insensitive


def _toy_frameworks():
    return {
        "IGHV3-23": GermlineFrameworks(
            name="IGHV3-23",
            fr1="GAGGTGCAGCTG",
            fr2="ATGAGCTGGGTC",
            fr3="TACTACGCAGAC",
            fr4="TGGGGCCAAGGG",
        )
    }


def _toy_clone(cdr3="ARDYW", red=5):
    return OasClone("IGHV3-23", "GFTFSSYA", "ISGSGGST", cdr3, red)


def test_germline_loader_strips_gaps(tmp_path):
    ref = tmp_path / "ref.txt"
    ref.write_text(
        "# comment\n"
        "human-IGHV3-23 gag-gtg CAG atg-agc TGG tac-tac nnnGCG tgg-ggc\n"
    )
    fw = load_germline_frameworks(str(ref))
    assert "IGHV3-23" in fw
    assert fw["IGHV3-23"].fr1 == "GAGGTG"  # dash removed, upper-cased
    assert "N" not in fw["IGHV3-23"].fr3


def test_scaffold_db_entry_shape():
    entries = build_germline_scaffold_db([_toy_clone()], _toy_frameworks())
    assert len(entries) == 1
    name, fields = entries[0]
    assert name == "germline_IGHV3-23"
    assert len(fields) == 7  # FR1 CDR1 FR2 CDR2 FR3 CDR3 FR4
    assert fields[5] == back_translate("ARDYW")


def test_assemble_heavy_chain_concatenates_all_regions():
    seq = assemble_heavy_chain(_toy_clone(), _toy_frameworks()["IGHV3-23"])
    assert seq.startswith("GAGGTGCAGCTG")
    assert back_translate("ARDYW") in seq


def test_apply_errors_no_error_model_is_identity():
    import random

    quiet = ErrorModel(substitution=0.0, deletion=0.0, insertion=0.0)
    seq = "ACGTACGTACGT"
    out, quals = apply_errors(seq, quiet, random.Random(0))
    assert out == seq
    assert len(quals) == len(seq)


def test_picking_functional_hcdr3_dedup_rank():
    import pandas as pd

    from radiant.picking import pick_leads

    df = pd.DataFrame(
        [
            # two rows share an HCDR3 -> collapse to the larger; keep abundance order
            {"cluster_size": 100, "functional_1": "functional", "sequence_issue": "",
             "Merged_CDRs_AA": "AAAABBBBCARWY", "HCDR3_AA": "CARWY"},
            {"cluster_size": 30, "functional_1": "functional", "sequence_issue": "",
             "Merged_CDRs_AA": "AAAABBBBCARWY", "HCDR3_AA": "CARWY"},
            {"cluster_size": 60, "functional_1": "functional", "sequence_issue": "",
             "Merged_CDRs_AA": "AAAABBBBCARGG", "HCDR3_AA": "CARGG"},
            # non-functional (stop) is dropped
            {"cluster_size": 90, "functional_1": "functional", "sequence_issue": "",
             "Merged_CDRs_AA": "AAAABBBBCA*GG", "HCDR3_AA": "CA*GG"},
            # flagged non-functional is dropped
            {"cluster_size": 80, "functional_1": "non-functional", "sequence_issue": "",
             "Merged_CDRs_AA": "AAAABBBBCARHH", "HCDR3_AA": "CARHH"},
        ]
    )
    leads = pick_leads(df, single_chain=True, top_n=None)
    assert list(leads["HCDR3_AA"]) == ["CARWY", "CARGG"]  # dedup + abundance rank
    assert list(leads["pick_rank"]) == [1, 2]
    assert abs(leads["rel_freq"].sum() - 1.0) < 1e-9


def test_picking_top_n():
    import pandas as pd

    from radiant.picking import pick_leads

    df = pd.DataFrame(
        [
            {"cluster_size": s, "functional_1": "functional", "sequence_issue": "",
             "Merged_CDRs_AA": f"AAAABBBBCAR{aa}", "HCDR3_AA": f"CAR{aa}"}
            for s, aa in [(50, "AA"), (40, "CC"), (30, "DD")]
        ]
    )
    leads = pick_leads(df, single_chain=True, top_n=2)
    assert len(leads) == 2
    assert list(leads["HCDR3_AA"]) == ["CARAA", "CARCC"]


def test_simulation_is_reproducible():
    clones = [_toy_clone(cdr3=c, red=r) for c, r in [("ARDYW", 9), ("ARGGW", 3)]]
    sim = build_sim_clones(clones, _toy_frameworks(), n_clones=2)
    reads_a = list(simulate_reads(sim, ErrorModel(), seed=7))
    reads_b = list(simulate_reads(sim, ErrorModel(), seed=7))
    assert reads_a == reads_b
    assert len(reads_a) == sum(c.depth for c in sim)
