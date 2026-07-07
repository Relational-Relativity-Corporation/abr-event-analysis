"""
declaration/targets.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

Declared comparison targets — all PDG-traceable.
All targets declared BEFORE any operator runs.
A target declared after results are observed is inadmissible.

Per Section 9.1 of primary_operators_delta_sigma_v6.md:
Every comparison target must be traceable to an observable through M
by the same standard as operator inputs.
"""

from .admissibility import AdmissibilityError


def check_comparison_target(value, description: str, pdg_reference: str):
    """
    Admissibility check for comparison targets.
    Raises AdmissibilityError if pdg_reference is empty.
    """
    if not pdg_reference or not pdg_reference.strip():
        raise AdmissibilityError(
            f"Comparison target '{description}' has no declared PDG reference. "
            f"Every comparison target must be traceable to a PDG detector "
            f"measurement through M. Declare the reference before execution."
        )
    return value


# ── Declared admissible comparison targets ────────────────────────────────
# All declared here, before any operator runs.
# Source cited explicitly for each.

DECLARED_TARGETS = {

    "D0_MASS_MEV": check_comparison_target(
        value=1864.84,
        description="D0 meson mass (MeV/c^2)",
        pdg_reference="PDG 2022, Particle Listings: D0, Table of masses"
    ),

    "D0_LIFETIME_PS": check_comparison_target(
        value=0.4101,
        description="D0 meson mean lifetime (ps)",
        pdg_reference="PDG 2022, Particle Listings: D0, mean lifetime"
    ),

    "KAON_MASS_MEV": check_comparison_target(
        value=493.677,
        description="K- meson mass (MeV/c^2)",
        pdg_reference="PDG 2022, Particle Listings: K+/K-, mass"
    ),

    "PION_MASS_MEV": check_comparison_target(
        value=139.570,
        description="pi+ meson mass (MeV/c^2)",
        pdg_reference="PDG 2022, Particle Listings: pi+/pi-, mass"
    ),

    "EPSILON_CP": check_comparison_target(
        value=2.228e-3,
        description="CP violation parameter |epsilon| for neutral kaon system",
        pdg_reference="PDG 2022, CP Violation in Kaons, Table 1"
    ),

    "MAGNETIC_MOMENT_RATIO": check_comparison_target(
        value=-1.4599,
        description="Proton-to-neutron magnetic moment ratio mu_p / mu_n",
        pdg_reference=(
            "PDG 2022, Particle Listings: "
            "proton mu_p = 2.79285, neutron mu_n = -1.91304 nuclear magnetons"
        )
    ),
}
