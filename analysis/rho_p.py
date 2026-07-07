"""
analysis/rho_p.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

rho_P ratio: rank(Im Sigma) / propagation_capacity.

Primary Region: rho_P << 1.
Atomic/Molecular transition: rho_P ~= 1 — OPEN CONDITION.
The transition threshold is not yet derived as a computable criterion.
rho_P is computable. The threshold remains interpretive.
"""

import numpy as np
from typing import List, Tuple
from .rank import im_sigma_rank


def propagation_capacity(edges: List[Tuple[int, int]]) -> int:
    """Count of edges with at least one forward neighbor (adj+ not empty)."""
    targets = {t for _, t in edges}
    sources = {s for s, _ in edges}
    # Edge e has adj+ if target of e is a source of some other edge
    count = 0
    for _, t in edges:
        if t in sources:
            count += 1
    return count


def rho_p_ratio(
    sigma_field: np.ndarray,
    edges: List[Tuple[int, int]],
) -> float:
    """
    rho_P = rank(Im Sigma) / propagation_capacity.

    Declared projection:
      Preserves: ratio of Sigma output span dimension to propagation capacity.
      Discards:  individual edge values.

    OPEN CONDITION: rho_P ~= 1 transition threshold not yet derived.
    """
    rank_sigma, _ = im_sigma_rank(sigma_field)
    c_x = propagation_capacity(edges)
    if c_x == 0:
        return 0.0
    return rank_sigma / c_x
