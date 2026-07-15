# radiant-annotator

The **annotator** used by the RADiANT ONT pipeline (`radiant.annotator`),
distributed as a **run-only compiled binary**, no source is included. It assigns
framework and CDR fields to each read against a scaffold database, including the
combined-CDR nucleotide region used for clustering.

## Install

```
pip install radiant_annotator-1.0.0-*.whl
```

The wheel bundles the compiled `annotator` binary. To point at a binary
elsewhere, set `RADIANT_ANNOTATOR_BIN=/path/to/annotator`.

## CLI contract

```
annotator --input <folder> --output <folder> --db <scaffold_db> [--cpus N] [--single-chain]
```

Writes an annotation CSV (`<output>/annotation_output.csv`).

## Availability

The compiled annotator binary is provided by the maintainers (contact the
corresponding author). It is intentionally run-only; the pipeline that consumes
its output (this repository) is fully open.
