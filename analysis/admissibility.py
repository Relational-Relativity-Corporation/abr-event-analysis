"""
analysis/admissibility.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

Declared edge-image admissibility check.

Tests whether the declared edge-image is admissible under the V6
directional admissibility condition.

An admissible declared edge-image is ASYMMETRIC: no declared edge vector
has its negation as another declared edge vector.

A symmetric declared edge-image signals inadmissible structure — a reverse
edge was declared without satisfying the distinctness axiom and independent
provenance requirement.

This is an admissibility check, not an analysis function.
Returns (is_admissible: bool, witness: Optional[list])
  True  = admissible (asymmetric)
  False = inadmissible (symmetric — declaration error)

Declared projection:
  Preserves: admissibility classification and first inadmissibility witness.
  Discards:  row ordering; linear combinations not present as edges.
"""

import numpy as np
from typing import Optional, List, Tuple

SVD_TOLERANCE = 1e-10


def declared_edge_image_admissibility_check(
    delta_field: np.ndarray,
    tol: float = SVD_TOLERANCE,
) -> Tuple[bool, Optional[List[float]]]:
    """
    Check whether the declared edge-image is admissible (asymmetric).

    Returns:
        (True, None)           — admissible
        (False, witness_row)   — inadmissible; witness is the edge vector
                                 whose negation was found among declared edges
    """
    if delta_field.size == 0:
        return True, None

    rows = delta_field.tolist()

    for row in rows:
        neg = [-v for v in row]
        found = any(
            all(abs(a - b) < tol for a, b in zip(neg, other))
            for other in rows
        )
        if found:
            return False, row

    return True, None
