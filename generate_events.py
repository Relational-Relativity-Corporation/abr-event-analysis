"""
generate_events.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

Generates D0 -> K- pi+ decay events from declared physics.

M is declared over individual events — not over probability distributions,
not over QCD amplitudes, not over wave functions.

Generative parameters (PDG-traceable, used as inputs not comparison targets):
  D0 mass:     1864.84 MeV/c^2
  D0 lifetime: 0.4101 ps
  Kaon mass:   493.677 MeV/c^2
  Pion mass:   139.570 MeV/c^2

The decay is generated from:
  - Conservation of 4-momentum
  - Isotropic angular distribution in D0 rest frame (spin-0 parent)
  - Exponential decay time from declared lifetime

No Pythia8. No Monte Carlo sampling from probability distributions
as a black box. The generative structure is declared explicitly here.

Output: data/events/d0_events.csv
"""

import numpy as np
import pandas as pd
import os

from declaration.targets import DECLARED_TARGETS

# ── Generative parameters (PDG-traceable) ─────────────────────────────────

M_D0    = DECLARED_TARGETS["D0_MASS_MEV"]       # 1864.84 MeV/c^2
TAU_D0  = DECLARED_TARGETS["D0_LIFETIME_PS"]     # 0.4101 ps
M_KAON  = DECLARED_TARGETS["KAON_MASS_MEV"]     # 493.677 MeV/c^2
M_PION  = DECLARED_TARGETS["PION_MASS_MEV"]     # 139.570 MeV/c^2

C_MM_PS = 299.792458  # speed of light in mm/ps

# ── Declared D0 beam parameters ───────────────────────────────────────────
# Beam momentum declared by Origin — not sampled from a distribution.
# A fixed beam momentum is the most faithful declaration: one observable
# value, not an ensemble average.

D0_MOMENTUM_MEV = 5000.0   # MeV/c — declared beam momentum
D0_GAMMA = np.sqrt((D0_MOMENTUM_MEV / M_D0)**2 + 1)
D0_BETA  = D0_MOMENTUM_MEV / (D0_GAMMA * M_D0)


def decay_products_rest_frame(n_events: int, rng: np.random.Generator):
    """
    Generate kaon and pion momenta in the D0 rest frame.

    Decay kinematics from conservation of 4-momentum.
    Angular distribution: isotropic (spin-0 parent — no preferred axis).
    Isotropic means uniform in cos(theta) and phi.

    This is declared physics, not sampled from a QCD amplitude.
    """
    # Rest frame momentum magnitude (from 2-body kinematics)
    p_star = np.sqrt(
        (M_D0**2 - (M_KAON + M_PION)**2) *
        (M_D0**2 - (M_KAON - M_PION)**2)
    ) / (2 * M_D0)

    # Isotropic angular distribution — declared (spin-0 parent)
    cos_theta = rng.uniform(-1, 1, n_events)
    phi       = rng.uniform(0, 2 * np.pi, n_events)
    sin_theta = np.sqrt(1 - cos_theta**2)

    # Kaon momentum in rest frame
    kx = p_star * sin_theta * np.cos(phi)
    ky = p_star * sin_theta * np.sin(phi)
    kz = p_star * cos_theta

    return kx, ky, kz, p_star


def boost_to_lab(kx, ky, kz, beta, gamma):
    """
    Lorentz boost along z-axis from D0 rest frame to lab frame.

    The boost is along the beam direction (z).
    Transverse momenta are unchanged.
    """
    E_kaon_rest = np.sqrt(M_KAON**2 + kx**2 + ky**2 + kz**2)

    kaon_pz_lab = gamma * (kz + beta * E_kaon_rest)
    pion_pz_lab = gamma * (-kz + beta * np.sqrt(M_PION**2 + kx**2 + ky**2 + kz**2))

    return kx, ky, kaon_pz_lab, -kx, -ky, pion_pz_lab


def generate_events(n_events: int, seed: int = 42) -> pd.DataFrame:
    """
    Generate n_events D0 -> K- pi+ decay events.

    Each event is one row with 12 observable components.
    The observable vector is declared in declaration/observable.py.

    Returns: DataFrame with columns matching OBSERVABLE_COMPONENTS.
    """
    rng = np.random.default_rng(seed)

    # Decay times — exponential from declared lifetime
    # This is the one place where a distribution is used as a generative
    # structure — but it is declared here explicitly, not imported from
    # a probability framework. The exponential is derived from the
    # declared physics of radioactive decay, not assumed.
    decay_times_ps = rng.exponential(scale=TAU_D0, size=n_events)

    # Decay lengths from decay times and declared D0 momentum
    decay_lengths_mm = decay_times_ps * D0_BETA * C_MM_PS * D0_GAMMA

    # Production vertex (fixed at origin — declared)
    prod_x = np.zeros(n_events)
    prod_y = np.zeros(n_events)
    prod_z = np.zeros(n_events)

    # Decay vertex — along beam direction
    vertex_x = prod_x
    vertex_y = prod_y
    vertex_z = prod_z + decay_lengths_mm

    # Decay products in rest frame
    kx_rest, ky_rest, kz_rest, p_star = decay_products_rest_frame(n_events, rng)

    # Boost to lab frame
    kaon_px, kaon_py, kaon_pz, pion_px, pion_py, pion_pz = boost_to_lab(
        kx_rest, ky_rest, kz_rest, D0_BETA, D0_GAMMA
    )

    # Invariant mass (should reconstruct to M_D0 — use as check)
    E_kaon = np.sqrt(M_KAON**2 + kaon_px**2 + kaon_py**2 + kaon_pz**2)
    E_pion = np.sqrt(M_PION**2 + pion_px**2 + pion_py**2 + pion_pz**2)
    E_total = E_kaon + E_pion
    px_total = kaon_px + pion_px
    py_total = kaon_py + pion_py
    pz_total = kaon_pz + pion_pz
    invariant_mass = np.sqrt(np.maximum(
        E_total**2 - px_total**2 - py_total**2 - pz_total**2, 0
    ))

    df = pd.DataFrame({
        "event_id":      np.arange(n_events),
        "vertex_x":      vertex_x,
        "vertex_y":      vertex_y,
        "vertex_z":      vertex_z,
        "kaon_px":       kaon_px,
        "kaon_py":       kaon_py,
        "kaon_pz":       kaon_pz,
        "pion_px":       pion_px,
        "pion_py":       pion_py,
        "pion_pz":       pion_pz,
        "decay_length":  decay_lengths_mm,
        "invariant_mass": invariant_mass,
        "decay_time":    decay_times_ps,
    })

    return df


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate D0 -> K- pi+ decay events"
    )
    parser.add_argument("--n", type=int, default=10000,
                        help="Number of events to generate (default: 10000)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility")
    args = parser.parse_args()

    print("generate_events.py")
    print("Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.")
    print("=" * 60)
    print(f"Generating {args.n} D0 -> K- pi+ events")
    print(f"Generative parameters (PDG-traceable):")
    print(f"  D0 mass:     {M_D0} MeV/c^2")
    print(f"  D0 lifetime: {TAU_D0} ps")
    print(f"  Kaon mass:   {M_KAON} MeV/c^2")
    print(f"  Pion mass:   {M_PION} MeV/c^2")
    print(f"  D0 momentum: {D0_MOMENTUM_MEV} MeV/c (declared beam)")
    print()

    df = generate_events(args.n, args.seed)

    out_path = os.path.join(
        os.path.dirname(__file__), "data", "events", "d0_events.csv"
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df.to_csv(out_path, index=False)

    print(f"Events written: {out_path}")
    print(f"Shape: {df.shape}")
    print()
    print("Invariant mass check (should be near D0 mass = 1864.84 MeV/c^2):")
    print(f"  mean:  {df['invariant_mass'].mean():.4f} MeV/c^2")
    print(f"  std:   {df['invariant_mass'].std():.6f} MeV/c^2")
    print()
    print("Decay time check (mean should be near tau = 0.4101 ps):")
    print(f"  mean:  {df['decay_time'].mean():.4f} ps")
    print(f"  declared tau: {TAU_D0} ps")
    print()
    print("Bounded over D. No claim beyond D.")
