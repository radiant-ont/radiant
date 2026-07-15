# RADiANT

Recovery of antibody lead sequences from Oxford Nanopore (ONT) reads of
polyclonal display outputs.

RADiANT takes raw ONT FASTQ from an antibody display campaign and returns ranked,
annotated lead clones. It is a general analysis pipeline: any group can run their
own ONT reads through it and get the per-cluster output table. It is
library-agnostic, working against either designed-scaffold databases or germline
(natural-repertoire) family databases. The proprietary annotator is a
separately-distributed run-only component; the rest of the pipeline is open.

## Quickstart: reproduce the germline example

Run the bundled germline example end to end on any machine — no proprietary data
required. It exercises every stage (filter, annotate, cluster, consensus, output)
and reproduces the annotated output.

```bash
# 1. Get the code
git clone https://github.com/radiant-ont/radiant ~/radiant_reader
cd ~/radiant_reader

# 2. External tools (cutadapt, fastq-filter, minimap2, racon) + Python deps
bash env/setup_tools.sh
pip install -e .

# 3. The run-only annotator binary for your platform (see "The annotator")
gh release download v0.1.0 --repo radiant-ont/radiant --pattern '*macosx_11_0_arm64*'  # macOS arm64
# gh release download v0.1.0 --repo radiant-ont/radiant --pattern '*manylinux*'         # Linux x86_64
pip install ./radiant_annotator-1.0.0-*.whl

# 4a. Score recovery against the known truth set
python mock/run_mock_pipeline.py \
    --fastq data/generated/germline_mock.fastq \
    --scaffold-db data/generated/germline_scaffold_db.txt \
    --truth data/generated/germline_mock_truth.csv \
    --workdir runs/germline_mock

# 4b. ...or write the supplementary-format output table (one row per cluster,
#     with per-region frameworks/CDRs, liabilities, biophysical, confidence)
python -m radiant --scaffold-db data/generated/germline_scaffold_db.txt \
    --workdir runs/germline_table --format vhh --no-trim --output-table --pick \
    data/generated/germline_mock.fastq
```

Expected, deterministic: 585 reads -> 270 functional -> 15 clusters; 13/15 match a
simulated clone and 14/15 carry the correct HCDR3. Step 4b writes
`runs/germline_table/output_table.csv`. The germline mock carries no adapters, so
trimming is skipped (`--no-trim`); on real ONT data you also supply the adapter
sequences to trim.

## What it does

Raw ONT FASTQ in, ranked annotated clones out:

1. Trim adapters (cutadapt, linked anchored adapters).
2. Length/quality filter.
3. Annotate reads against a scaffold database (the `annotator` component).
4. Cluster reads by the combined-CDR nucleotide sequence.
5. Build a consensus per cluster (minimap2 + racon).
6. Re-annotate the consensus, pick the top clone per HCDR3, and write the output table.

## Default parameters

These are the code defaults, and the values reported in the paper.

| Stage | Parameter | Value |
|-------|-----------|-------|
| Trim | cutadapt, linked anchored adapters (fwd + rev), keep adapter-matched | `-e 4` |
| Filter | minimum length | 650 bp (scFv), 300 bp (VHH) |
| Filter | quality | mean read quality >= 30 (`fastq-filter -Q 30`) |
| Annotate | remove reads with a stop-codon | before clustering |
| Cluster | region | combined CDRs, nucleotide |
| Cluster | distance | Levenshtein LD <= 3 (scFv), <= 2 (VHH) |
| Cluster | linkage | single-linkage to the most-abundant seed |
| Cluster | minimum size | 5 reads |
| Consensus | aligner | minimap2 v2.28 |
| Consensus | polisher | racon v1.4.3, max error rate 0.7 |
| Compare | HCDR3 clusters | 100% HCDR3 amino-acid identity |
| Compare | top pick | most-abundant full-length AA per HCDR3 |

## External tools

cutadapt, plus the pinned minimap2 and racon versions:

```bash
bash env/setup_tools.sh
export MINIMAP2_PATH="$HOME/.local/bin/minimap2-2.28"
export RACON_PATH="$HOME/.local/bin/racon-1.4.3"
```

## The annotator

The annotation step is a separate, run-only compiled component
(`radiant-annotator`) that assigns framework/CDR fields against a scaffold
database. It is distributed as a wheel that bundles the compiled binary. Download
the wheel for your platform from the repository's Releases, then install it (the
wheel must be downloaded first — `pip install <name>.whl` reads a local file, not a
remote asset):

```bash
# macOS (Apple silicon / arm64)
gh release download v0.1.0 --repo radiant-ont/radiant --pattern '*macosx_11_0_arm64*'
pip install ./radiant_annotator-1.0.0-py3-none-macosx_11_0_arm64.whl

# Linux (x86_64)
gh release download v0.1.0 --repo radiant-ont/radiant --pattern '*manylinux*'
pip install ./radiant_annotator-1.0.0-py3-none-manylinux_2_34_x86_64.whl
```

Without the GitHub CLI, install directly from the asset URL, e.g.
`pip install "https://github.com/radiant-ont/radiant/releases/download/v0.1.0/radiant_annotator-1.0.0-py3-none-macosx_11_0_arm64.whl"`.

The wheel is platform-specific (it contains a native binary): a macOS arm64 build
is provided, with a Linux x86_64 build alongside it. To point at a binary you have
elsewhere, set `RADIANT_ANNOTATOR_BIN=/path/to/annotator`. The annotator is only
needed for the annotation stage; the trimming, clustering, and consensus code in
this repository is open. See `packaging/radiant_annotator/`.

## Running the pipeline

Every stage parameter is an explicit command-line flag with the manuscript value
as its default, so a run is fully specified by its command line:

```bash
python -m radiant --scaffold-db <DB> --workdir <OUTDIR> <reads.fastq> [...]
python -m radiant --help          # every parameter, with defaults
```

Each run prints the parameters it used and writes them to
`<OUTDIR>/run_parameters.json`, and writes `<OUTDIR>/read_accounting.csv`
tabulating how many reads survive each stage (raw, quality-filtered, annotated,
functional, clustered). The distance metric (`--distance-metric
levenshtein|hamming`), clustering region (`--region`), distance (`--distance`),
minimum cluster size (`--min-cluster-size`), quality gates, annotator confidence
(`--min-votes`), and racon error (`--racon-error`) are all flags; `--help` lists
them with defaults.

Add `--output-table` to also write `<OUTDIR>/output_table.csv`, the
supplementary-format table with one row per cluster: full-length and per-chain
(VL/VH) amino-acid and nucleotide sequences, per-region frameworks and CDRs,
germline match names, counts and votes, ROI and HCDR3 cluster labels,
liabilities, biophysical properties, and a power-law confidence score.

## Validation

Two checks, covering a designed-scaffold library and natural germline families:

```bash
# designed scaffold (Target 1)
python validation/compare.py     # 9/10 shared, 7/10 picks

# natural germline families (mock dataset, no proprietary data)
python mock/run_mock_pipeline.py --fastq data/generated/germline_mock.fastq \
    --scaffold-db data/generated/germline_scaffold_db.txt \
    --truth data/generated/germline_mock_truth.csv --workdir data/runs/germline_mock
```

Outputs are deterministic (identical on rerun). From 585 simulated reads the mock
retains 15 clusters (270 functional reads); 13/15 match a simulated clone by
combined-CDR amino acid and 14/15 carry the correct HCDR3 (41 of 48 simulated
clones reach the depth >= 5 needed to form a cluster).

## Bring your own annotation database

The pipeline is library-agnostic: annotation is driven entirely by the
`--scaffold-db` you supply. Point it at a designed-scaffold database (as the
manuscript did) or at a germline (natural-repertoire) family database — the same
plain space-separated format (`name FR1 CDR1 FR2 CDR2 FR3 CDR3 FR4`, nucleotide)
serves both. A ready-to-use germline database is included so the full workflow can
be run and its output reproduced without any additional data.

## Germline / natural repertoires

For natural repertoires, use a germline-family scaffold database. A ready-to-use
one is included (`data/generated/germline_scaffold_db.txt`), so ONT reads can be
run against germline scaffolds directly:

```bash
python -m radiant --scaffold-db data/generated/germline_scaffold_db.txt \
    --workdir <OUTDIR> --output-table <reads.fastq>
```

To build your own from an OAS-style table of germline calls and CDRs, an example
input (`data/generated/oas_human_db_subset.csv`, columns `v_call, cdr1_aa,
cdr2_aa, cdr3_aa`) and the human germline frameworks
(`scaffolds/reference/germline_frameworks.txt`) are included:

```bash
python scaffolds/build_scaffold_db.py \
    --oas data/generated/oas_human_db_subset.csv \
    --germline-ref scaffolds/reference/germline_frameworks.txt \
    --out data/generated/germline_scaffold_db.txt
```

The scaffold database is a plain space-separated table
(`name FR1 CDR1 FR2 CDR2 FR3 CDR3 FR4`, nucleotide), so any germline or designed
set can be supplied in that format.

## Layout

| Path | Contents |
|------|----------|
| `radiant/` | pipeline stages: preprocess, cluster, consensus, picking, annotator interface |
| `validation/` | `compare.py` (ONT vs the Sanger reference) |
| `mock/` | germline read simulation + mock benchmark |
| `scaffolds/` | scaffold-database construction (OAS to germline families) |
| `campaigns/` | pinned run configs |
| `packaging/` | the run-only annotator wheel |
| `env/` | tool provisioning |
| `tests/` | unit tests |

## Tests

```bash
python -m pytest tests/ -q
```

## License

Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved. This
software is proprietary and is provided solely to evaluate and reproduce the
results in the associated publication. You may download and run it for that
purpose; redistribution, commercial use, derivative works, and reverse
engineering of the compiled annotator are not permitted. See [LICENSE](LICENSE)
for the full terms.
