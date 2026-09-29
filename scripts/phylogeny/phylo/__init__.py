"""phylo — phylogenetic inference on single-cell WGS panels.

Scanpy-style API mirroring `cellspec`:
    - pp : QC + per-worm somatic panel construction
    - tl : G0 diagnostics, tree builders, DNA-internal validation
    - pl : plotting
    - utils : config + run management

Intended to fold into `cellspec.tl.phylo`.
"""

from . import pl, pp, tl, utils

__all__ = ["pl", "pp", "tl", "utils"]
__version__ = "0.0.1"
