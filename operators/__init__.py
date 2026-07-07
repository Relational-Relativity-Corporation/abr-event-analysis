"""
operators/__init__.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

V6 primary operator kernel: Delta and Sigma.
Exports: operator_delta, operator_sigma, antisymmetric_term,
         compute_rho, operator_e_primary.
"""
from .primary import (
    operator_delta,
    operator_sigma,
    antisymmetric_term,
    compute_rho,
    operator_e_primary,
)
