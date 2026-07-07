"""
analysis/__init__.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

Analysis functions: rank(Im Delta), declared edge-image admissibility,
Sigma antisymmetric term summary, expression condition, rho_P ratio.

All functions are declared C projections of operator outputs.
What each projection preserves and discards is stated in its docstring.
"""
from .rank import im_delta_rank, im_sigma_rank, SVD_TOLERANCE
from .admissibility import declared_edge_image_admissibility_check
from .expression import expression_condition, detect_failure_mode, FailureMode
from .rho_p import rho_p_ratio
