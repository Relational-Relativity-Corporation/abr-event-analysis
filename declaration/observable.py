"""
declaration/observable.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

M declaration: the observable vector for a single D0 -> K- pi+ decay event.

Every component must be directly traceable to a detector measurement.
No theoretical quantity. No statistical summary. No sub-observable construct.

Declared M:
    x[v] = (vertex_x, vertex_y, vertex_z,
             kaon_px, kaon_py, kaon_pz,
             pion_px, pion_py, pion_pz,
             decay_length, invariant_mass, decay_time)

12 dimensions. One observable vector per declared event.
"""

import numpy as np
from dataclasses import dataclass
from .admissibility import AdmissibilityError

# Observable vector components — declared before any operator acts.
# Each component must be directly measurable from detector data.
# No quantity derived from a probability distribution.
# No quantity derived from theoretical scaffolding.
OBSERVABLE_COMPONENTS = [
    "vertex_x",       # decay vertex x-coordinate (mm) — from track intersection
    "vertex_y",       # decay vertex y-coordinate (mm)
    "vertex_z",       # decay vertex z-coordinate (mm)
    "kaon_px",        # kaon track x-momentum (MeV/c) — from track curvature
    "kaon_py",        # kaon track y-momentum (MeV/c)
    "kaon_pz",        # kaon track z-momentum (MeV/c)
    "pion_px",        # pion track x-momentum (MeV/c)
    "pion_py",        # pion track y-momentum (MeV/c)
    "pion_pz",        # pion track z-momentum (MeV/c)
    "decay_length",   # distance from production to decay vertex (mm)
    "invariant_mass", # reconstructed D0 mass (MeV/c^2) — from track momenta
    "decay_time",     # decay time (ps) — from decay_length / momentum
]

N_COMPONENTS = len(OBSERVABLE_COMPONENTS)  # 12

# Inadmissible quantities — not permitted as components of x[v].
# Each encodes theoretical or statistical structure, not a direct measurement.
INADMISSIBLE_COMPONENTS = [
    "probability",          # not an observable
    "amplitude",            # QM amplitude — not directly measurable
    "wave_function",        # QM wave function — not directly measurable
    "quark_content",        # sub-observable — below declared boundary
    "color_charge",         # sub-observable — not measurable
    "ckm_element",          # theoretical — quark mixing matrix
    "branching_ratio",      # statistical — ratio over ensemble
    "cross_section",        # statistical — ratio over ensemble
    "decay_rate",           # statistical — inverse lifetime (mean over ensemble)
    "distribution_mean",    # statistical summary — not individual event
    "distribution_width",   # statistical summary — not individual event
]


@dataclass
class ObservableVector:
    """
    The declared observable vector for one D0 -> K- pi+ decay event.

    All values must be finite (in D).
    All values must be traceable to a detector measurement through M.
    Declaration is not a substitute for traceability.
    """
    event_id: int
    data: np.ndarray  # shape (N_COMPONENTS,)

    def __post_init__(self):
        if len(self.data) != N_COMPONENTS:
            raise AdmissibilityError(
                f"Observable vector must have {N_COMPONENTS} components. "
                f"Got {len(self.data)}. Every component of x[v] must be "
                f"declared and traceable to a detector measurement through M."
            )
        if not np.all(np.isfinite(self.data)):
            raise AdmissibilityError(
                "Observable vector contains non-finite values. "
                "All quantities must be in D := { x | |x[i]| < inf }."
            )

    def __getitem__(self, component: str) -> float:
        if component not in OBSERVABLE_COMPONENTS:
            raise AdmissibilityError(
                f"Component '{component}' is not a declared observable. "
                f"Declared components: {OBSERVABLE_COMPONENTS}"
            )
        return float(self.data[OBSERVABLE_COMPONENTS.index(component)])
