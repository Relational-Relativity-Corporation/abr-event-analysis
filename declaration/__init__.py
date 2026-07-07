"""
declaration/__init__.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

M declaration package.
Exports: OBSERVABLE_COMPONENTS, load_events, declare_relations, check_comparison_target.
"""
from .observable import OBSERVABLE_COMPONENTS, ObservableVector
from .relations import declare_relations
from .targets import check_comparison_target, DECLARED_TARGETS
from .admissibility import AdmissibilityError
