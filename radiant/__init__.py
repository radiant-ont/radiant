# Copyright (c) 2026 Specifica, an IQVIA business. All rights reserved.
# Licensed for reproduction use only; redistribution prohibited. See LICENSE.

"""RADiANT: reproducible Oxford Nanopore antibody-lead recovery.

This package accompanies the RADiANT manuscript. It provides the open,
library-agnostic parts of the workflow (scaffold-database construction,
mock read simulation, trimming, filtering, clustering, and consensus
building) together with a thin, neutral ``annotator`` interface. The
antibody annotator itself is distributed separately as a compiled
component so that the workflow can be run without exposing that source.

The clustering and consensus steps operate on read sequences alone and are
independent of library origin: the workflow accepts any scaffold database
supplied as space-separated framework and CDR fields, whether proprietary
scaffolds or germline families.
"""

__version__ = "0.1.0"
