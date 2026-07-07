"""
analysis/rank.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

rank(Im Delta) and rank(Im Sigma) — declared C projections.

Preserves: dimension of the span of declared edge vectors.
Discards:  individual edge values and their ordering.

Wording: rank(Im Delta) is the dimension of the span — not the count of
linearly independent directed differences. Those are equivalent only when
every edge vector is linearly independent of the others.

SVD_TOLERANCE declared before execution. Not adjusted post-hoc.
"""

import numpy as np

SVD_TOLERANCE = 1e-10


def im_delta_rank(delta_field: np.ndarray):
    """
    Compute rank(Im Delta) — dimension of span of declared edge vectors.

    Declared projection:
      Preserves: dimension of span and singular value spectrum.
      Discards:  individual edge values and ordering.

    Returns:
        rank: int
        singular_values: numpy array
    """
    if delta_field.size == 0:
        return 0, np.array([])
    U, S, Vt = np.linalg.svd(delta_field, full_matrices=False)
    rank = int(np.sum(S > SVD_TOLERANCE))
    return rank, S


def im_sigma_rank(sigma_field: np.ndarray):
    """
    Compute rank(Im Sigma) — dimension of span of Sigma output vectors.

    Declared projection:
      Preserves: dimension of span and singular value spectrum.
      Discards:  individual edge values and ordering.
    """
    return im_delta_rank(sigma_field)
