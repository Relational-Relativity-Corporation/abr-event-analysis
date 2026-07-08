"""
analysis/rank.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

rank(Im Delta) and rank(Im Sigma) -- declared C projections.
SVD_TOLERANCE declared before execution. Not adjusted post-hoc.
"""

import numpy as np

SVD_TOLERANCE = 1e-10


def im_delta_rank(delta_field):
    if delta_field.size == 0:
        return 0, np.array([])
    U, S, Vt = np.linalg.svd(delta_field, full_matrices=False)
    rank = int(np.sum(S > SVD_TOLERANCE))
    return rank, S


def im_sigma_rank(sigma_field):
    return im_delta_rank(sigma_field)
