# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""Neutral interface to the antibody annotator.

The annotator assigns framework and CDR fields to each read against a scaffold
database. Its implementation is distributed as a separate compiled component
and is intentionally not part of this open package. This module locates a
backend and exposes a single ``annotate`` entry point returning a table of
annotated reads (including the combined-CDR nucleotide region used downstream).

Backend resolution order:

1. ``radiant_annotator`` if importable (the distributed compiled component).
2. A module named by the ``RADIANT_ANNOTATOR_BACKEND`` environment variable
   (used for development against a local build).

A backend must provide a callable::

    annotate_folder(input_folder, output_folder, scaffold_db, single_chain,
                    cpus[, min_votes]) -> str  # path to an annotation CSV

``min_votes`` is optional: it is passed only to backends whose signature accepts
it (or accepts ``**kwargs``), so older backends keep working unchanged.
"""

from __future__ import annotations

import importlib
import inspect
import os

import pandas as pd

DEFAULT_MIN_VOTES = 150


def _load_backend():
    try:
        return importlib.import_module("radiant_annotator")
    except ImportError:
        pass
    name = os.environ.get("RADIANT_ANNOTATOR_BACKEND")
    if name:
        return importlib.import_module(name)
    raise RuntimeError(
        "No annotator backend available. Install the compiled annotator "
        "(radiant-annotator) or set RADIANT_ANNOTATOR_BACKEND to a "
        "development backend module."
    )


def annotate(
    input_folder: str,
    output_folder: str,
    scaffold_db: str,
    single_chain: bool = True,
    cpus: int = 1,
    min_votes: int = DEFAULT_MIN_VOTES,
) -> pd.DataFrame:
    """Annotate every read file in ``input_folder`` against ``scaffold_db``.

    ``min_votes`` sets the annotator confidence threshold; it is forwarded only
    if the resolved backend accepts it (otherwise the backend default is used).

    Returns
    -------
    pandas.DataFrame
        Annotated reads. Includes at least ``read``, ``count`` and the
        combined-CDR nucleotide column used for clustering.
    """
    backend = _load_backend()
    kwargs = dict(
        input_folder=input_folder,
        output_folder=output_folder,
        scaffold_db=scaffold_db,
        single_chain=single_chain,
        cpus=cpus,
    )
    try:
        params = inspect.signature(backend.annotate_folder).parameters
        accepts = any(p.kind == p.VAR_KEYWORD for p in params.values()) or "min_votes" in params
    except (TypeError, ValueError):
        accepts = True
    if accepts:
        kwargs["min_votes"] = min_votes
    result = backend.annotate_folder(**kwargs)
    # Backends may return the CSV path, or a (csv_path, output_folder) tuple.
    csv_path = result[0] if isinstance(result, tuple) else result
    return pd.read_csv(csv_path)
