"""
declaration/admissibility.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

AdmissibilityError: raised when a quantity cannot be traced to an observable
through M. Declaration is not a substitute for traceability.
"""


class AdmissibilityError(Exception):
    """
    Raised when a node, edge, comparison target, or quantity cannot be
    traced to an observable through M.

    Declaration is never a substitute for traceability to an observable
    through M. A quantity that computes but cannot be so traced is
    fabricated within D.
    """
    pass
