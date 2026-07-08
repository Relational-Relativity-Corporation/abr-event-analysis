"""
analyze_spin1_recovery.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

Question 1 -- Recovery: Does Delta/Sigma recover the QM spin-1 prediction?
Fixed: both kaon and pion use the same full 3D boost function.
Pion rest-frame momentum is (-kx, -ky, -kz) by conservation (2-body decay).
"""

import numpy as np
import pandas as pd
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from declaration.observable import OBSERVABLE_COMPONENTS
from declaration.relations import declare_relations, RelationProvenance
from declaration.targets import DECLARED_TARGETS
from operators.primary import operator_delta, antisymmetric_term
from analysis.rank import im_delta_rank, SVD_TOLERANCE

EVENTS_PATH_SPIN0 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "events", "d0_events.csv")
EVENTS_PATH_SPIN1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "events", "dstar_events.csv")
REPORT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "spin1_recovery_report.md")

M_D0    = DECLARED_TARGETS["D0_MASS_MEV"]
M_KAON  = DECLARED_TARGETS["KAON_MASS_MEV"]
M_PION  = DECLARED_TARGETS["PION_MASS_MEV"]
M_DSTAR = 2010.26
M_PISOFT = 139.570
ALPHA_QM = 1.0
QM_COS2_MEAN = 0.4  # <cos^2> for 1 + cos^2(theta) distribution


def boost_to_lab(px_r, py_r, pz_r, mass, par_px, par_py, par_pz, par_mass):
    """Full 3D boost from parent rest frame to lab frame."""
    par_p = np.sqrt(par_px**2 + par_py**2 + par_pz**2)
    par_E = np.sqrt(par_mass**2 + par_p**2)
    beta = par_p / par_E
    gamma = par_E / par_mass
    ux, uy, uz = par_px/par_p, par_py/par_p, par_pz/par_p
    E_r = np.sqrt(mass**2 + px_r**2 + py_r**2 + pz_r**2)
    p_par = px_r*ux + py_r*uy + pz_r*uz
    p_par_lab = gamma * (p_par + beta * E_r)
    px_perp = px_r - p_par*ux
    py_perp = py_r - p_par*uy
    pz_perp = pz_r - p_par*uz
    return px_perp + p_par_lab*ux, py_perp + p_par_lab*uy, pz_perp + p_par_lab*uz


def boost_to_rest(px_l, py_l, pz_l, mass, par_px, par_py, par_pz, par_mass):
    """Full 3D inverse boost from lab to parent rest frame."""
    par_p = np.sqrt(par_px**2 + par_py**2 + par_pz**2)
    par_E = np.sqrt(par_mass**2 + par_p**2)
    beta = par_p / par_E
    gamma = par_E / par_mass
    ux, uy, uz = par_px/par_p, par_py/par_p, par_pz/par_p
    E_l = np.sqrt(mass**2 + px_l**2 + py_l**2 + pz_l**2)
    p_par = px_l*ux + py_l*uy + pz_l*uz
    p_par_r = gamma * (p_par - beta * E_l)
    px_perp = px_l - p_par*ux
    py_perp = py_l - p_par*uy
    pz_perp = pz_l - p_par*uz
    return px_perp + p_par_r*ux, py_perp + p_par_r*uy, pz_perp + p_par_r*uz


def sample_cos_theta_spin1(n, rng, alpha=1.0):
    """Sample cos(theta) from 1 + alpha*cos^2(theta) via rejection sampling."""
    f_max = 1.0 + abs(alpha)
    samples = []
    while len(samples) < n:
        ct = rng.uniform(-1, 1, n * 3)
        u = rng.uniform(0, f_max, len(ct))
        samples.extend(ct[u < (1.0 + alpha * ct**2)].tolist())
    return np.array(samples[:n])


def generate_dstar_events(n_events, seed=42):
    """
    Generate D*+ -> D0 pi+ -> (K- pi+) pi+ events.
    Both kaon and pion use the same full 3D boost from D0 rest frame to lab.
    """
    rng = np.random.default_rng(seed)
    C_MM_PS = 299.792458
    TAU_D0 = DECLARED_TARGETS["D0_LIFETIME_PS"]

    # D* -> D0 rest-frame momentum
    p_dstar_rf = np.sqrt(
        (M_DSTAR**2 - (M_D0 + M_PISOFT)**2) *
        (M_DSTAR**2 - (M_D0 - M_PISOFT)**2)
    ) / (2 * M_DSTAR)

    # D* moves along z in lab
    dstar_p_lab = 5000.0
    dstar_px_lab = np.zeros(n_events)
    dstar_py_lab = np.zeros(n_events)
    dstar_pz_lab = np.full(n_events, dstar_p_lab)

    # Sample D0 direction in D* rest frame from 1 + cos^2(theta)
    cos_t = sample_cos_theta_spin1(n_events, rng, alpha=ALPHA_QM)
    phi = rng.uniform(0, 2 * np.pi, n_events)
    sin_t = np.sqrt(1 - cos_t**2)

    d0_px_rf = p_dstar_rf * sin_t * np.cos(phi)
    d0_py_rf = p_dstar_rf * sin_t * np.sin(phi)
    d0_pz_rf = p_dstar_rf * cos_t

    # Boost D0 to lab
    d0_px, d0_py, d0_pz = boost_to_lab(
        d0_px_rf, d0_py_rf, d0_pz_rf, M_D0,
        dstar_px_lab, dstar_py_lab, dstar_pz_lab, M_DSTAR
    )

    # D0 -> K- pi+ in D0 rest frame (isotropic)
    p_star = np.sqrt(
        (M_D0**2 - (M_KAON + M_PION)**2) *
        (M_D0**2 - (M_KAON - M_PION)**2)
    ) / (2 * M_D0)

    ct2 = rng.uniform(-1, 1, n_events)
    phi2 = rng.uniform(0, 2 * np.pi, n_events)
    st2 = np.sqrt(1 - ct2**2)

    kx_rf = p_star * st2 * np.cos(phi2)
    ky_rf = p_star * st2 * np.sin(phi2)
    kz_rf = p_star * ct2

    # Boost kaon to lab
    kaon_px, kaon_py, kaon_pz = boost_to_lab(
        kx_rf, ky_rf, kz_rf, M_KAON,
        d0_px, d0_py, d0_pz, M_D0
    )

    # Boost pion to lab -- equal and opposite in D0 rest frame
    pion_px, pion_py, pion_pz = boost_to_lab(
        -kx_rf, -ky_rf, -kz_rf, M_PION,
        d0_px, d0_py, d0_pz, M_D0
    )

    # Invariant mass check
    E_k = np.sqrt(M_KAON**2 + kaon_px**2 + kaon_py**2 + kaon_pz**2)
    E_p = np.sqrt(M_PION**2 + pion_px**2 + pion_py**2 + pion_pz**2)
    inv_mass = np.sqrt(np.maximum(
        (E_k+E_p)**2 - (kaon_px+pion_px)**2
        - (kaon_py+pion_py)**2 - (kaon_pz+pion_pz)**2, 0
    ))

    # Decay length and vertex
    d0_p_mag = np.sqrt(d0_px**2 + d0_py**2 + d0_pz**2)
    d0_E = np.sqrt(M_D0**2 + d0_p_mag**2)
    d0_beta_v = d0_p_mag / d0_E
    d0_gamma_v = d0_E / M_D0

    decay_times = rng.exponential(scale=TAU_D0, size=n_events)
    decay_lengths = decay_times * d0_beta_v * C_MM_PS * d0_gamma_v

    return pd.DataFrame({
        "event_id":        np.arange(n_events),
        "vertex_x":        decay_lengths * d0_px/d0_p_mag,
        "vertex_y":        decay_lengths * d0_py/d0_p_mag,
        "vertex_z":        decay_lengths * d0_pz/d0_p_mag,
        "kaon_px":         kaon_px, "kaon_py": kaon_py, "kaon_pz": kaon_pz,
        "pion_px":         pion_px, "pion_py": pion_py, "pion_pz": pion_pz,
        "decay_length":    decay_lengths,
        "invariant_mass":  inv_mass,
        "decay_time":      decay_times,
        "cos_theta_dstar": cos_t,
    })


def load_field(path):
    df = pd.read_csv(path)
    x = {int(row["event_id"]): np.array([row[c] for c in OBSERVABLE_COMPONENTS])
         for _, row in df.iterrows()}
    return df, x


def delta_rank(x, edges, comps):
    idx = [OBSERVABLE_COMPONENTS.index(c) for c in comps]
    field = operator_delta({e: v[idx] for e, v in x.items()}, edges)
    r, sv = im_delta_rank(field)
    return r, field, sv


def make_edges(x):
    return declare_relations(
        list(x.keys()), RelationProvenance.TEMPORAL, "declared production order"
    )


def run():
    print("\nanalyze_spin1_recovery.py (fixed)")
    print("Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.")
    print("=" * 60)
    print()

    print("Generating D*+ -> D0 pi+ events...")
    df1 = generate_dstar_events(10000)
    df1.to_csv(EVENTS_PATH_SPIN1, index=False)

    imean = df1["invariant_mass"].mean()
    istd  = df1["invariant_mass"].std()
    print(f"  Invariant mass: mean={imean:.4f}  std={istd:.6f}  "
          f"(should be {M_D0:.2f}, std~0)")

    if abs(imean - M_D0) > 1.0:
        print("ERROR: invariant mass wrong. Stopping.")
        return
    print("  Invariant mass check: PASS")
    print()

    df0, x0 = load_field(EVENTS_PATH_SPIN0)
    _,   x1 = load_field(EVENTS_PATH_SPIN1)
    e0, e1 = make_edges(x0), make_edges(x1)

    # ── Step 1: rank ─────────────────────────────────────────────────────

    print("=" * 60)
    print("STEP 1: rank(Im Delta)")
    r0, _, _ = delta_rank(x0, e0, OBSERVABLE_COMPONENTS)
    r1, _, _ = delta_rank(x1, e1, OBSERVABLE_COMPONENTS)
    r0k, f0k, _ = delta_rank(x0, e0, ["kaon_px","kaon_py","kaon_pz"])
    r1k, f1k, _ = delta_rank(x1, e1, ["kaon_px","kaon_py","kaon_pz"])
    print(f"  Full 12-component: spin-0={r0}  spin-1={r1}")
    print(f"  Kaon momentum only: spin-0={r0k}  spin-1={r1k}")
    print()

    # ── Step 2: angular recovery ──────────────────────────────────────────

    print("=" * 60)
    print("STEP 2: Angular distribution recovery")

    def reco_cos2(df, x_field):
        cts = []
        for _, row in df.iterrows():
            kpx,kpy,kpz = row["kaon_px"],row["kaon_py"],row["kaon_pz"]
            ppx,ppy,ppz = row["pion_px"],row["pion_py"],row["pion_pz"]
            d0px,d0py,d0pz = kpx+ppx, kpy+ppy, kpz+ppz
            kx_r,ky_r,kz_r = boost_to_rest(
                kpx,kpy,kpz, M_KAON, d0px,d0py,d0pz, M_D0)
            d0p = np.sqrt(d0px**2+d0py**2+d0pz**2)
            km  = np.sqrt(kx_r**2+ky_r**2+kz_r**2)
            ct  = (kx_r*d0px/d0p + ky_r*d0py/d0p + kz_r*d0pz/d0p) / km
            cts.append(ct)
        return np.array(cts)

    ct0 = reco_cos2(df0, x0)
    ct1 = reco_cos2(df1, x1)
    truth_cos2 = (df1["cos_theta_dstar"].values**2).mean()

    print(f"  Generator truth <cos^2>: {truth_cos2:.4f}  (QM: {QM_COS2_MEAN:.4f})")
    print(f"  Spin-0 reco <cos^2>: {(ct0**2).mean():.4f}  (isotropic: 0.3333)")
    print(f"  Spin-1 reco <cos^2>: {(ct1**2).mean():.4f}  (QM: {QM_COS2_MEAN:.4f})")
    print()

    reco_c2 = (ct1**2).mean()

    # ── Step 3: Sigma distinction ─────────────────────────────────────────

    print("=" * 60)
    print("STEP 3: Sigma antisymmetric term -- spin-0 vs spin-1")

    x_a0 = {eid: np.array([ct]) for eid,ct in zip(x0.keys(), ct0)}
    x_a1 = {eid: np.array([ct]) for eid,ct in zip(x1.keys(), ct1)}

    fa0 = operator_delta(x_a0, e0)
    fa1 = operator_delta(x_a1, e1)

    as0 = np.abs(antisymmetric_term(fa0, e0, 0.2)).flatten()
    as1 = np.abs(antisymmetric_term(fa1, e1, 0.2)).flatten()

    nz0 = as0[as0 > SVD_TOLERANCE]
    nz1 = as1[as1 > SVD_TOLERANCE]

    print(f"  Spin-0: mean={nz0.mean():.6f}  std={nz0.std():.6f}")
    print(f"  Spin-1: mean={nz1.mean():.6f}  std={nz1.std():.6f}")
    ratio = nz1.mean() / nz0.mean()
    print(f"  Ratio spin-1/spin-0: {ratio:.4f}")
    print()

    # ── Step 4: checks ────────────────────────────────────────────────────

    print("=" * 60)
    print("STEP 4: Recovery assessment")

    err1  = abs(reco_c2 - QM_COS2_MEAN)
    c1    = err1 < 0.02
    c2    = abs(ratio - 1.0) > 0.05
    c3    = r0 != r1

    print(f"  Check 1 <cos^2> matches QM: reco={reco_c2:.4f} "
          f"QM={QM_COS2_MEAN:.4f} err={err1:.4f} "
          f"-> {'PASS' if c1 else 'FAIL'}")
    print(f"  Check 2 Sigma distinguishes: ratio={ratio:.4f} "
          f"dev={abs(ratio-1.0):.4f} "
          f"-> {'PASS' if c2 else 'FAIL'}")
    print(f"  Check 3 rank differs: {r0} vs {r1} "
          f"-> {'PASS' if c3 else 'SAME'}")
    print()

    overall = "SUPPORTED" if (c1 and c2) else "PARTIAL" if c1 else "NOT SUPPORTED"
    print(f"OVERALL: Question 1 (Recovery) -- {overall}")
    print()
    if overall == "SUPPORTED":
        print("  The operator framework recovers the QM spin-1 prediction.")
        print("  Ready for Question 2: does the framework detect structure")
        print("  BEYOND what the QM prediction accounts for?")
    elif overall == "PARTIAL":
        print("  <cos^2> matches QM but Sigma does not distinguish.")
        print("  Delta rank may carry the distinction.")
    else:
        print("  Recovery not established. Inspect declaration.")

    os.makedirs(os.path.dirname(REPORT_PATH), exist_ok=True)
    with open(REPORT_PATH, "w", encoding="utf-8") as f:
        f.write("# Spin-1 Recovery Analysis (fixed)\n\n")
        f.write("**Metatron Dynamics, Inc.** V6. Bounded over D. No claim beyond D.\n\n")
        f.write(f"Invariant mass mean: {imean:.4f} MeV/c^2\n\n")
        f.write(f"Spin-0 rank: {r0}  Spin-1 rank: {r1}\n\n")
        f.write(f"Generator truth <cos^2>: {truth_cos2:.4f} (QM: {QM_COS2_MEAN})\n\n")
        f.write(f"Spin-0 reco <cos^2>: {(ct0**2).mean():.4f}\n\n")
        f.write(f"Spin-1 reco <cos^2>: {reco_c2:.4f} (QM: {QM_COS2_MEAN})\n\n")
        f.write(f"Sigma ratio: {ratio:.4f}\n\n")
        f.write(f"Check 1: {'PASS' if c1 else 'FAIL'}\n\n")
        f.write(f"Check 2: {'PASS' if c2 else 'FAIL'}\n\n")
        f.write(f"Check 3: {'PASS' if c3 else 'SAME'}\n\n")
        f.write(f"Overall: {overall}\n\n")

    print(f"Report written: {REPORT_PATH}")
    print("\nBounded over D. No claim beyond D.")


if __name__ == "__main__":
    run()
