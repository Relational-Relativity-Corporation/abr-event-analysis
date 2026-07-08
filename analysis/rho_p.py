"""
analysis/rho_p.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

rho_P ratio: rank(Im Sigma) / propagation_capacity.
Primary Region: rho_P << 1.
Atomic/Molecular transition: rho_P ~= 1 -- OPEN CONDITION.
"""

import numpy as np
from typing import List, Tuple
from analysis.rank import im_sigma_rank


def propagation_capacity(edges):
    sources = {s for s, _ in edges}
    return sum(1 for _, t in edges if t in sources)


def rho_p_ratio(sigma_field, edges):
    rank_sigma, _ = im_sigma_rank(sigma_field)
    c_x = propagation_capacity(edges)
    if c_x == 0:
        return 0.0
    return rank_sigma / c_x
