"""
analysis/expression.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

Expression condition and failure mode detection.

Two legitimate failure modes for admissible declared structures:
  FM1 — Differentiation Collapse: rank(Im Delta) = 0
  FM2 — Relational Isolation: adj+(e) = empty and adj-(e) = empty for all e

Circulation Cancellation is an inadmissibility signal, not a legitimate
failure mode. It requires symmetric adjacency, which requires inadmissible
reverse edges under the V6 directional admissibility condition.
"""

import numpy as np
from enum import Enum
from typing import List, Tuple

from .rank import im_delta_rank, im_sigma_rank, SVD_TOLERANCE
from ..operators.primary import antisymmetric_term, build_adjacency


class FailureMode(Enum):
    NONE = "none"
    DIFFERENTIATION_COLLAPSE = "differentiation_collapse"
    RELATIONAL_ISOLATION = "relational_isolation"
    CIRCULATION_CANCELLATION = "circulation_cancellation"  # inadmissibility signal


def detect_failure_mode(
    delta_field: np.ndarray,
    sigma_field: np.ndarray,
    edges: List[Tuple[int, int]],
    rho_base: float = 0.2,
) -> FailureMode:
    """Detect which failure mode is present, if any."""
    # FM1
    rank_delta, _ = im_delta_rank(delta_field)
    if rank_delta == 0:
        return FailureMode.DIFFERENTIATION_COLLAPSE

    # FM2
    adj_plus, adj_minus = build_adjacency(edges)
    all_isolated = all(
        len(adj_plus[e]) == 0 and len(adj_minus[e]) == 0
        for e in range(len(edges))
    )
    if all_isolated:
        return FailureMode.RELATIONAL_ISOLATION

    # Circulation cancellation (inadmissibility signal)
    asym = antisymmetric_term(delta_field, edges, rho_base)
    asym_all_zero = np.all(np.abs(asym) <= SVD_TOLERANCE)
    rank_sigma, _ = im_sigma_rank(sigma_field)
    if asym_all_zero and rank_sigma > 0:
        return FailureMode.CIRCULATION_CANCELLATION

    return FailureMode.NONE


def expression_condition(
    delta_field: np.ndarray,
    sigma_field: np.ndarray,
    edges: List[Tuple[int, int]],
    rho_base: float = 0.2,
) -> Tuple[bool, int, bool]:
    """
    Expression is admissible when:
      1. rank(Im Sigma) > 0
      2. Antisymmetric term non-zero on at least one edge

    Returns: (expressed, rank_sigma, has_nonzero_asym)
    """
    rank_sigma, _ = im_sigma_rank(sigma_field)
    asym = antisymmetric_term(delta_field, edges, rho_base)
    has_nonzero_asym = np.any(np.abs(asym) > SVD_TOLERANCE)
    expressed = rank_sigma > 0 and has_nonzero_asym
    return expressed, rank_sigma, has_nonzero_asym
