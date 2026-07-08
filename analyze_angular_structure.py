"""
analyze_angular_structure.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

Compares the antisymmetric term of Sigma over the kaon momentum
angular components against the QM isotropic prediction.

QM prediction (Born rule, spin-0 parent):
  The angular distribution of decay products is isotropic in the
  D0 rest frame. |psi|^2 is uniform in cos(theta) and phi.
  This predicts no preferred direction -- no angular asymmetry.

Operator question:
  Does Sigma detect asymmetric relational structure in the kaon
  momentum field that the isotropic prediction would not produce?

Method:
  1. Declare M over kaon momentum components only (px, py, pz).
  2. Compute Delta and Sigma over the declared temporal edges.
  3. Compute the antisymmetric term of Sigma at each edge.
  4. Compare the magnitude distribution against what isotropic
     sampling would produce (declared C projection).
  5. Compute the angular field in the D0 rest frame directly
     and apply the same analysis.

All comparison targets declared before execution.
No result adjusted post-hoc.

Declared projection of antisymmetric term:
  Preserves: magnitude of antisymmetric circulation per edge,
             distribution of magnitudes across declared edges.
  Discards:  individual edge component values, edge ordering.
"""

import numpy as np
import pandas as pd
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

from declaration.observable import OBSERVABLE_COMPONENTS
from declaration.relations import declare_relations, RelationProvenance
from declaration.targets import DECLARED_TARGETS
from operators.primary import operator_delta, operator_sigma, antisymmetric_term
from analysis.rank import im_delta_rank, im_sigma_rank, SVD_TOLERANCE

EVENTS_PATH = os.path.join(
    os.path.dirname(__file__), "data", "events", "d0_events.csv"
)
REPORT_PATH = os.path.join(
    os.path.dirname(__file__), "results", "angular_structure_report.md"
)

# Declared generative parameters
M_D0   = DECLARED_TARGETS["D0_MASS_MEV"]
M_KAON = DECLARED_TARGETS["KAON_MASS_MEV"]
M_PION = DECLARED_TARGETS["PION_MASS_MEV"]
D0_MOMENTUM_MEV = 5000.0
D0_GAMMA = np.sqrt((D0_MOMENTUM_MEV / M_D0)**2 + 1)
D0_BETA  = D0_MOMENTUM_MEV / (D0_GAMMA * M_D0)

# QM prediction: isotropic angular distribution in D0 rest frame.
# For a spin-0 parent decaying to two pseudoscalars:
#   d(Gamma)/d(cos theta) = constant
#   d(Gamma)/d(phi) = constant
# This is the Born rule prediction for this decay.
# It is a declared C projection: it discards the individual event
# angular values and retains only the distribution shape.
QM_PREDICTION = "isotropic: uniform in cos(theta) and phi in D0 rest frame"


def boost_to_rest_frame(kaon_px, kaon_py, kaon_pz, beta, gamma):
    """
    Inverse Lorentz boost: lab frame -> D0 rest frame.
    Boost is along z-axis (beam direction).
    """
    E_kaon_lab = np.sqrt(M_KAON**2 + kaon_px**2 + kaon_py**2 + kaon_pz**2)
    kaon_pz_rest = gamma * (kaon_pz - beta * E_kaon_lab)
    # Transverse components unchanged
    return kaon_px, kaon_py, kaon_pz_rest


def compute_angular_components(kaon_px, kaon_py, kaon_pz_rest):
    """
    Compute cos(theta) and phi in D0 rest frame.

    theta: polar angle from beam axis (z)
    phi: azimuthal angle in transverse plane

    Declared projection:
      Preserves: angular direction of kaon in D0 rest frame.
      Discards: momentum magnitude (p_star is fixed by kinematics).
    """
    p_mag = np.sqrt(kaon_px**2 + kaon_py**2 + kaon_pz_rest**2)
    cos_theta = kaon_pz_rest / p_mag
    phi = np.arctan2(kaon_py, kaon_px)
    return cos_theta, phi


def rank_of_components(x_sub, edges):
    field = operator_delta(x_sub, edges)
    rank, sv = im_delta_rank(field)
    return rank, field, sv


def antisymmetric_analysis(delta_field, edges, label):
    """
    Compute antisymmetric term of Sigma and report its distribution.

    Declared projection:
      Preserves: magnitude distribution of antisymmetric term.
      Discards:  individual component values and edge ordering.
    """
    asym = antisymmetric_term(delta_field, edges, rho_base=0.2)
    # Magnitude per edge (norm across components)
    magnitudes = np.linalg.norm(asym, axis=1)
    nonzero = magnitudes[magnitudes > SVD_TOLERANCE]

    print(f"  {label}:")
    print(f"    Edges with nonzero antisymmetric term: "
          f"{len(nonzero)} / {len(magnitudes)}")
    if len(nonzero) > 0:
        print(f"    Magnitude -- mean:   {nonzero.mean():.6f}")
        print(f"    Magnitude -- std:    {nonzero.std():.6f}")
        print(f"    Magnitude -- max:    {nonzero.max():.6f}")
        print(f"    Magnitude -- median: {np.median(nonzero):.6f}")
    return magnitudes, nonzero


def run():
    print("\nanalyze_angular_structure.py")
    print("Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.")
    print("=" * 60)
    print()
    print(f"QM prediction: {QM_PREDICTION}")
    print()

    if not os.path.exists(EVENTS_PATH):
        print(f"Event data not found: {EVENTS_PATH}")
        sys.exit(1)

    df = pd.read_csv(EVENTS_PATH)
    n_events = len(df)
    print(f"Events loaded: {n_events}")

    # Build full observable field
    x_full = {}
    for _, row in df.iterrows():
        eid = int(row["event_id"])
        vec = np.array([row[c] for c in OBSERVABLE_COMPONENTS], dtype=float)
        x_full[eid] = vec

    event_ids = list(x_full.keys())
    edges = declare_relations(
        event_ids=event_ids,
        relation_type=RelationProvenance.TEMPORAL,
        direction_basis="declared production order"
    )
    print(f"Declared edges: {len(edges)}")
    print()

    # ── Step 1: Lab-frame kaon momentum ──────────────────────────────────

    print("=" * 60)
    print("STEP 1: Kaon momentum in lab frame (px, py, pz)")
    print()

    kaon_idx = [OBSERVABLE_COMPONENTS.index(c)
                for c in ["kaon_px", "kaon_py", "kaon_pz"]]
    x_kaon_lab = {eid: vec[kaon_idx] for eid, vec in x_full.items()}
    rank_lab, field_lab, sv_lab = rank_of_components(x_kaon_lab, edges)

    print(f"  rank(Im Delta): {rank_lab}")
    print(f"  Singular values: {sv_lab[:6].round(4)}")
    print()
    mags_lab, nonzero_lab = antisymmetric_analysis(field_lab, edges,
                                                    "Sigma antisymmetric term")
    print()

    # ── Step 2: Angular components in D0 rest frame ───────────────────────

    print("=" * 60)
    print("STEP 2: Angular components in D0 rest frame (cos_theta, phi)")
    print()
    print(f"  QM prediction: uniform distribution in cos(theta) and phi")
    print(f"  If isotropic: no preferred direction, Sigma should detect")
    print(f"  only the relational contrast between adjacent events,")
    print(f"  which should be symmetric on average.")
    print()

    # Compute rest-frame angular components for each event
    x_angular = {}
    cos_thetas = []
    phis = []

    for _, row in df.iterrows():
        eid = int(row["event_id"])
        kpx = row["kaon_px"]
        kpy = row["kaon_py"]
        kpz = row["kaon_pz"]

        # Boost to rest frame
        kpx_r, kpy_r, kpz_r = boost_to_rest_frame(kpx, kpy, kpz,
                                                    D0_BETA, D0_GAMMA)
        cos_t, phi = compute_angular_components(kpx_r, kpy_r, kpz_r)

        x_angular[eid] = np.array([cos_t, phi])
        cos_thetas.append(cos_t)
        phis.append(phi)

    cos_thetas = np.array(cos_thetas)
    phis = np.array(phis)

    print(f"  cos(theta) -- mean: {cos_thetas.mean():.4f}  "
          f"std: {cos_thetas.std():.4f}  "
          f"(isotropic predicts mean=0, std=1/sqrt(3)={1/np.sqrt(3):.4f})")
    print(f"  phi        -- mean: {phis.mean():.4f}  "
          f"std: {phis.std():.4f}  "
          f"(isotropic predicts mean=0, std=pi/sqrt(3)={np.pi/np.sqrt(3):.4f})")
    print()

    rank_ang, field_ang, sv_ang = rank_of_components(x_angular, edges)
    print(f"  rank(Im Delta) over angular components: {rank_ang}")
    print(f"  Singular values: {sv_ang[:4].round(4)}")
    print()
    mags_ang, nonzero_ang = antisymmetric_analysis(field_ang, edges,
                                                    "Sigma antisymmetric term")
    print()

    # ── Step 3: Isotropy test ─────────────────────────────────────────────

    print("=" * 60)
    print("STEP 3: Isotropy test via Delta rank over angular bins")
    print()
    print("  If the distribution is isotropic, dividing the angular")
    print("  space into bins and computing rank(Im Delta) over the")
    print("  bin field should produce rank <= 2 (the angular space")
    print("  has 2 independent directions: cos_theta and phi).")
    print()

    # Bin the angular distribution into a 10x10 grid
    N_BINS = 10
    cos_bins = np.linspace(-1, 1, N_BINS + 1)
    phi_bins = np.linspace(-np.pi, np.pi, N_BINS + 1)

    # Count events per bin -- this is a C projection
    # Preserves: bin counts (frequency structure)
    # Discards: individual event angular values
    bin_counts = np.zeros((N_BINS, N_BINS))
    for ct, ph in zip(cos_thetas, phis):
        i = min(int((ct + 1) / 2 * N_BINS), N_BINS - 1)
        j = min(int((ph + np.pi) / (2 * np.pi) * N_BINS), N_BINS - 1)
        bin_counts[i, j] += 1

    expected_per_bin = n_events / (N_BINS * N_BINS)
    deviation = bin_counts - expected_per_bin
    max_dev = np.abs(deviation).max()
    mean_dev = np.abs(deviation).mean()
    chi2 = np.sum((bin_counts - expected_per_bin)**2 / expected_per_bin)

    print(f"  Angular bin grid: {N_BINS}x{N_BINS}")
    print(f"  Expected per bin: {expected_per_bin:.1f}")
    print(f"  Max deviation from expected: {max_dev:.1f}")
    print(f"  Mean |deviation|: {mean_dev:.1f}")
    print(f"  Chi-squared (isotropy test): {chi2:.2f}")
    print(f"  Chi-squared / dof: {chi2 / (N_BINS*N_BINS - 1):.3f}")
    print(f"  (Near 1.0 = consistent with isotropy)")
    print()

    # ── Step 4: Sigma output on isotropically-shuffled field ─────────────

    print("=" * 60)
    print("STEP 4: Sigma comparison -- declared vs shuffled angular field")
    print()
    print("  Shuffle: randomly permute the angular components across events.")
    print("  A shuffled field is isotropic by construction.")
    print("  Compare Sigma antisymmetric term magnitudes.")
    print("  If declared field differs from shuffled: non-isotropic structure")
    print("  is present that the QM prediction does not account for.")
    print()

    rng = np.random.default_rng(seed=42)
    perm = rng.permutation(n_events)

    cos_shuffled = cos_thetas[perm]
    phi_shuffled = phis[perm]

    x_shuffled = {}
    for i, eid in enumerate(event_ids):
        x_shuffled[eid] = np.array([cos_shuffled[i], phi_shuffled[i]])

    _, field_shuffled, _ = rank_of_components(x_shuffled, edges)
    mags_shuf, nonzero_shuf = antisymmetric_analysis(
        field_shuffled, edges, "Shuffled (isotropic) antisymmetric term"
    )

    print()
    print("  Declared angular field:")
    if len(nonzero_ang) > 0:
        print(f"    Mean antisymmetric magnitude: {nonzero_ang.mean():.6f}")
        print(f"    Std:                          {nonzero_ang.std():.6f}")
    print()
    print("  Shuffled (isotropic by construction) field:")
    if len(nonzero_shuf) > 0:
        print(f"    Mean antisymmetric magnitude: {nonzero_shuf.mean():.6f}")
        print(f"    Std:                          {nonzero_shuf.std():.6f}")
    print()

    if len(nonzero_ang) > 0 and len(nonzero_shuf) > 0:
        ratio = nonzero_ang.mean() / nonzero_shuf.mean()
        print(f"  Ratio (declared / shuffled): {ratio:.4f}")
        if abs(ratio - 1.0) < 0.05:
            finding = "CONSISTENT WITH ISOTROPY"
            detail = ("Sigma detects no more antisymmetric structure in the "
                      "declared angular field than in the shuffled field.")
        elif ratio > 1.0:
            finding = "EXCEEDS ISOTROPIC BASELINE"
            detail = ("Sigma detects more antisymmetric structure in the "
                      "declared field than isotropy predicts. Additional "
                      "relational structure present beyond QM Born rule.")
        else:
            finding = "BELOW ISOTROPIC BASELINE"
            detail = ("Sigma detects less antisymmetric structure than "
                      "the shuffled field. Inspect declaration.")
        print()
        print(f"  Finding: {finding}")
        print(f"  {detail}")

    # ── Summary ───────────────────────────────────────────────────────────

    print()
    print("=" * 60)
    print("SUMMARY")
    print()
    print(f"  rank(Im Delta) -- kaon lab momentum:    {rank_lab}")
    print(f"  rank(Im Delta) -- angular (rest frame): {rank_ang}")
    print(f"  Chi-squared / dof (isotropy):           {chi2/(N_BINS*N_BINS-1):.3f}")
    if len(nonzero_ang) > 0 and len(nonzero_shuf) > 0:
        print(f"  Sigma asym ratio (declared/shuffled):   {ratio:.4f}")
    print()
    print("  Statistical reduction discards:")
    print("    -- The full 3D momentum vector of the decay products")
    print("    -- The angular distribution in the D0 rest frame")
    print("    -- The event-to-event variation in all of the above")
    print()
    print("  Operators act on: the relational field over individual")
    print("  events before any of this reduction occurs.")

    # Write report
    os.makedirs(os.path.dirname(REPORT_PATH), exist_ok=True)
    with open(REPORT_PATH, "w", encoding="utf-8") as f:
        f.write("# Angular Structure Analysis\n\n")
        f.write("**Metatron Dynamics, Inc.** V6. "
                "Bounded over D. No claim beyond D.\n\n")
        f.write("## QM prediction\n\n")
        f.write(f"{QM_PREDICTION}\n\n")
        f.write("## Results\n\n")
        f.write(f"n_events: {n_events}\n\n")
        f.write(f"rank(Im Delta) kaon lab momentum: {rank_lab}\n\n")
        f.write(f"rank(Im Delta) angular rest frame: {rank_ang}\n\n")
        f.write(f"cos(theta) mean: {cos_thetas.mean():.4f} "
                f"(isotropic: 0)\n\n")
        f.write(f"phi mean: {phis.mean():.4f} (isotropic: 0)\n\n")
        f.write(f"Chi-squared / dof: {chi2/(N_BINS*N_BINS-1):.3f} "
                f"(isotropic: ~1.0)\n\n")
        if len(nonzero_ang) > 0 and len(nonzero_shuf) > 0:
            f.write(f"Sigma asym ratio (declared/shuffled): {ratio:.4f}\n\n")
            f.write(f"Finding: {finding}\n\n")
            f.write(f"{detail}\n\n")
        f.write("## Declared projection\n\n")
        f.write("Preserves: rank structure, antisymmetric term magnitude "
                "distribution, chi-squared vs isotropy.\n\n")
        f.write("Discards: individual edge values and ordering.\n\n")

    print()
    print(f"Report written: {REPORT_PATH}")
    print("\nBounded over D. No claim beyond D.")


if __name__ == "__main__":
    run()
