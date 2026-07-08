"""
analysis/admissibility.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

Declared edge-image admissibility check.
Returns (is_admissible: bool, witness).
True = admissible (asymmetric). False = inadmissible (declaration error).
"""

import numpy as np

SVD_TOLERANCE = 1e-10


def declared_edge_image_admissibility_check(delta_field, tol=SVD_TOLERANCE):
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
