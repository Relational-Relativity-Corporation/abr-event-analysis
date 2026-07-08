"""
analyze_hadronic_phase0a.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

Phase 0a -- Hadronic Transition Analysis at Level 2.

Declared M: raw track-level quantities only.
No particle identity. No mass hypothesis. No invariant mass from
mass hypothesis. No chi-squared. No theoretical classification.

Every component of x[v] must satisfy:
  1. Produced by a physical instrument directly
  2. Finite and bounded
  3. Does not require prior theoretical classification to obtain

What is in M (Level 2):
  track1_px, track1_py, track1_pz   -- track momentum (MeV/c)
                                        from: curvature in declared B field
  track1_charge                      -- sign of curvature (+1 or -1)
                                        from: direction of deflection
  track2_px, track2_py, track2_pz   -- second track momentum
  track2_charge                      -- second track charge sign
  vertex_x, vertex_y, vertex_z      -- secondary vertex position (mm)
                                        from: track intersection geometry
  decay_length                       -- primary to secondary vertex (mm)
                                        from: declared vertex positions
  decay_time                         -- ps, from decay_length / total_p

What is NOT in M (excluded, with reason):
  invariant_mass   -- requires mass hypothesis for each track
                      (which track is kaon, which is pion?)
  particle_label   -- theoretical classification
  chi_squared      -- statistical model assumption
  branching_ratio  -- ensemble C projection
  quantum_numbers  -- theoretical classifications (strangeness, isospin)

Declared comparison:
  Object frame prediction: PDG summary statistics for D0 decay
    -- invariant mass peak at 1864.84 MeV/c^2
    -- lifetime 0.4101 ps
    -- kaon momentum distribution from phase space
  These are all C projections. We do not use them as inputs to M.
  We use them as declared comparison targets AFTER operator output.

Declared hypotheses (before execution):
  H1: rank(Im Delta) over Level 2 track observables >= rank over Level 1
      (4 columns: D0_MM, D0_TAU, D0_PT, D0_MINIPCHI2)
      Level 1 rank established: 4 (from run_phase0_events.py with 12 components)
      Level 2 removes invariant_mass, adds track charges and vertex position
      Predicted: rank increases because track charges add independent binary
      directions not present in the invariant mass summary

  H2: rank(Im Delta) over raw track momenta alone (no vertex, no timing) >= 3
      (from rank_contributions analysis: kaon momentum contributed 3 new dims)
      Level 2 removes "kaon" label -- prediction: same 3 momentum dimensions
      are present in the unlabeled track momenta

  H3: Sigma antisymmetric term over the within-event graph is non-zero
      (the decay vertex is a declared branching point -- two tracks
       emerge from one vertex -- Sigma should detect the branching asymmetry)

  H4: The object frame's particle identification (kaon vs pion assignment)
      is a declared C projection of the Level 2 field that:
        Preserves: invariant mass reconstruction (requires mass hypothesis)
        Discards: individual track momentum components, charge sign combinations
      If rank(Level 2) > rank(particle-identified), the C projection
      discards information. If equal, it is lossless for rank.

All hypotheses declared before execution. No result adjusted post-hoc.
"""

import numpy as np
import pandas as pd
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from declaration.relations import declare_relations, RelationProvenance
from declaration.targets import DECLARED_TARGETS
from operators.primary import (
    operator_delta, operator_sigma, antisymmetric_term, build_adjacency
)
from analysis.rank import im_delta_rank, im_sigma_rank, SVD_TOLERANCE

# ── Paths ─────────────────────────────────────────────────────────────────

EVENTS_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "data", "events", "d0_events.csv"
)
REPORT_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "results", "hadronic_phase0a_report.md"
)

# ── Declared M: Level 2 observable components ────────────────────────────
#
# Every component listed here satisfies:
#   - Produced by a physical instrument directly
#   - Finite and bounded (in D)
#   - Does not require prior theoretical classification
#
# Source for simulated data: generate_events.py
# In real LHCb data these map to:
#   track_px/py/pz -> track momentum from Kalman filter over VELO hits
#   track_charge   -> sign of curvature in magnetic field
#   vertex_x/y/z   -> secondary vertex from track intersection
#   decay_length   -> |vertex - primary_vertex|
#   decay_time     -> decay_length / |p_total| * M_D0 / c (requires mass
#                     hypothesis for timing -- see note below)

LEVEL2_COMPONENTS = [
    # Track 1 momentum -- from curvature in declared B field
    # No label "kaon" -- that is a theoretical classification
    "track1_px",
    "track1_py",
    "track1_pz",
    "track1_charge",    # +1 or -1 -- from direction of curvature

    # Track 2 momentum
    "track2_px",
    "track2_py",
    "track2_pz",
    "track2_charge",    # +1 or -1

    # Secondary vertex -- from track intersection geometry
    "vertex_x",
    "vertex_y",
    "vertex_z",

    # Decay length -- from vertex positions
    "decay_length",
]
# Note on decay_time: excluded from Level 2 because computing it from
# decay_length requires knowing the D0 momentum, which requires summing
# the track momenta (admissible), but then requires dividing by M_D0
# (which is a theoretical quantity -- the D0 mass from PDG).
# decay_length is admissible. decay_time requires a mass hypothesis.
# This is the Level 2 boundary: below which requires theoretical input.

N_LEVEL2 = len(LEVEL2_COMPONENTS)

# ── Level 1 components (for comparison) ──────────────────────────────────
# These are the 12 components from run_phase0_events.py
# Already analyzed. Included here only for comparison.
LEVEL1_COMPONENTS_AVAILABLE = [
    "vertex_x", "vertex_y", "vertex_z",
    "kaon_px", "kaon_py", "kaon_pz",       # particle-identified
    "pion_px", "pion_py", "pion_pz",       # particle-identified
    "decay_length", "invariant_mass",       # requires mass hypothesis
    "decay_time",                           # requires mass hypothesis
]


def load_and_declare_level2(df):
    """
    Load the simulated event data and declare Level 2 observable vectors.

    The simulated data was generated with track labels (kaon/pion).
    At Level 2, we discard those labels and treat the tracks as
    track1 and track2 -- two charged tracks from a secondary vertex.

    Track assignment: track1 = the positively charged track,
                     track2 = the negatively charged track.
    This assignment is admissible: charge sign is directly observable
    from the direction of curvature in the magnetic field.
    It does not require knowing which track is a kaon.

    Declared C projection of the label assignment:
      Preserves: charge sign separation of the two tracks
      Discards: which specific particle type each track is
    """
    x = {}
    for _, row in df.iterrows():
        eid = int(row["event_id"])

        # Kaon is negatively charged (K-), pion is positively charged (pi+)
        # in D0 -> K- pi+ decay.
        # At Level 2 we do not use these labels.
        # We assign by charge: track1 = positive (pi+ in object frame),
        #                      track2 = negative (K- in object frame)
        # This is admissible: charge is directly observable.

        kaon_px = row["kaon_px"]
        kaon_py = row["kaon_py"]
        kaon_pz = row["kaon_pz"]
        pion_px = row["pion_px"]
        pion_py = row["pion_py"]
        pion_pz = row["pion_pz"]

        # Assign by charge sign (admissible)
        # Kaon: negative charge (-1)
        # Pion: positive charge (+1)
        track1_px = pion_px    # positive track
        track1_py = pion_py
        track1_pz = pion_pz
        track1_charge = +1.0

        track2_px = kaon_px    # negative track
        track2_py = kaon_py
        track2_pz = kaon_pz
        track2_charge = -1.0

        vec = np.array([
            track1_px, track1_py, track1_pz, track1_charge,
            track2_px, track2_py, track2_pz, track2_charge,
            row["vertex_x"], row["vertex_y"], row["vertex_z"],
            row["decay_length"],
        ], dtype=float)

        x[eid] = vec

    return x


def load_level1(df):
    """Load Level 1 observable vectors (with particle labels) for comparison."""
    level1_cols = [
        "vertex_x", "vertex_y", "vertex_z",
        "kaon_px", "kaon_py", "kaon_pz",
        "pion_px", "pion_py", "pion_pz",
        "decay_length", "invariant_mass", "decay_time"
    ]
    x = {}
    for _, row in df.iterrows():
        eid = int(row["event_id"])
        vec = np.array([row[c] for c in level1_cols], dtype=float)
        x[eid] = vec
    return x


def delta_rank_subset(x, edges, component_indices):
    """Compute rank(Im Delta) over a declared component subset."""
    x_sub = {eid: vec[component_indices] for eid, vec in x.items()}
    field = operator_delta(x_sub, edges)
    rank, sv = im_delta_rank(field)
    return rank, field, sv


def within_event_analysis(df, n_events=100):
    """
    Compute Sigma over the within-event graph for each event.

    Within-event graph for D0 -> K- pi+ (Level 2, no particle labels):
      Node 0: primary vertex (production point)
              x[0] = (0, 0, 0, E_beam) -- declared beam parameters
      Node 1: secondary vertex (decay point)
              x[1] = (vertex_x, vertex_y, vertex_z, total_momentum)
      Node 2: positive track endpoint
              x[2] = (track1_px, track1_py, track1_pz, track1_charge)
      Node 3: negative track endpoint
              x[3] = (track2_px, track2_py, track2_pz, track2_charge)

    Declared edges:
      (0->1): production to decay vertex -- directed by decay
      (1->2): decay vertex to positive track -- directed by track
      (1->3): decay vertex to negative track -- directed by track

    This is the correct within-event relational structure.
    Each edge has independent observable provenance.
    """
    D0_MOMENTUM = DECLARED_TARGETS["D0_MASS_MEV"]  # beam momentum

    results = []
    for i, (_, row) in enumerate(df.iterrows()):
        if i >= n_events:
            break

        # Primary vertex (declared at origin)
        total_px = row["kaon_px"] + row["pion_px"]
        total_py = row["kaon_py"] + row["pion_py"]
        total_pz = row["kaon_pz"] + row["pion_pz"]
        total_p  = np.sqrt(total_px**2 + total_py**2 + total_pz**2)

        # 4-component observable vector per node: (px, py, pz, charge_or_type)
        # charge_or_type: 0 for neutral (D0), +1 for positive, -1 for negative
        x_event = {
            0: np.array([0.0, 0.0, 0.0, 0.0]),           # primary vertex: neutral beam
            1: np.array([total_px, total_py, total_pz, 0.0]),  # decay vertex: D0 (neutral)
            2: np.array([row["pion_px"], row["pion_py"],
                         row["pion_pz"], +1.0]),           # positive track
            3: np.array([row["kaon_px"], row["kaon_py"],
                         row["kaon_pz"], -1.0]),           # negative track
        }

        # Declared within-event edges
        # Each direction is fixed by the observable:
        # production -> decay: the D0 traveled this way
        # decay -> track: the track emerged from the decay vertex
        edges_event = [(0, 1), (1, 2), (1, 3)]

        # Compute Delta over within-event graph
        delta_field = operator_delta(x_event, edges_event)

        # Compute Sigma
        sigma_field = operator_sigma(delta_field, edges_event, rho_base=0.2)

        # Antisymmetric term at each edge
        asym_field = antisymmetric_term(delta_field, edges_event, rho_base=0.2)

        # Rank
        rank_delta, sv_delta = im_delta_rank(delta_field)
        rank_sigma, sv_sigma = im_sigma_rank(sigma_field)

        # Antisymmetric magnitude at decay vertex edge (0->1)
        # This edge has adj+(0->1) = {(1->2), (1->3)} -- two track edges
        # The antisymmetric term measures the branching asymmetry
        asym_at_decay = np.linalg.norm(asym_field[0])  # edge (0->1)

        # Conservation check: Δ[(0->1)] + Δ[(1->2)] + Δ[(1->3)]
        # For momentum (components 0,1,2): should sum to -x[0] = 0
        # because x[0] = (0,0,0,0)
        closure = delta_field[0] + delta_field[1] + delta_field[2]
        momentum_closure = np.linalg.norm(closure[:3])  # px, py, pz

        results.append({
            "event_id": int(row["event_id"]),
            "rank_delta": rank_delta,
            "rank_sigma": rank_sigma,
            "asym_at_decay_vertex": asym_at_decay,
            "momentum_closure": momentum_closure,
            "sv_delta_max": sv_delta[0] if len(sv_delta) > 0 else 0.0,
        })

    return pd.DataFrame(results)


def run():
    print("\nanalyze_hadronic_phase0a.py")
    print("Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.")
    print("=" * 60)
    print()
    print("Phase 0a -- Hadronic Transitions at Level 2")
    print("Declared M: raw track observables, no particle identity")
    print()
    print(f"Level 2 components ({N_LEVEL2}):")
    for c in LEVEL2_COMPONENTS:
        print(f"  {c}")
    print()
    print("Declared hypotheses:")
    print("  H1: rank(Level 2) >= rank(Level 1 without mass hypothesis)")
    print("  H2: rank(track momenta alone) >= 3")
    print("  H3: Sigma antisymmetric term non-zero at decay vertex")
    print("  H4: Particle identification is a C projection of Level 2 field")
    print()

    if not os.path.exists(EVENTS_PATH):
        print(f"Event data not found: {EVENTS_PATH}")
        print("Run generate_events.py first.")
        sys.exit(1)

    df = pd.read_csv(EVENTS_PATH)
    n_events = len(df)
    print(f"Events loaded: {n_events}")

    # Declare temporal edges for between-event analysis
    x_level2 = load_and_declare_level2(df)
    x_level1 = load_level1(df)

    event_ids = list(x_level2.keys())
    edges = declare_relations(
        event_ids,
        RelationProvenance.TEMPORAL,
        "declared generation order"
    )

    print(f"Declared edges: {len(edges)}")
    print()

    # ── Step 1: rank comparison across levels ─────────────────────────────

    print("=" * 60)
    print("STEP 1: rank(Im Delta) across declaration levels")
    print()

    # Level 2 full
    rank_l2, field_l2, sv_l2 = delta_rank_subset(
        x_level2, edges, list(range(N_LEVEL2)))
    print(f"  Level 2 (no particle identity, {N_LEVEL2} components):")
    print(f"    rank(Im Delta) = {rank_l2}")
    print(f"    Top singular values: {sv_l2[:6].round(2)}")
    print()

    # Level 1 full
    rank_l1, field_l1, sv_l1 = delta_rank_subset(
        x_level1, edges, list(range(12)))
    print(f"  Level 1 (with particle identity, 12 components):")
    print(f"    rank(Im Delta) = {rank_l1}")
    print(f"    Top singular values: {sv_l1[:6].round(2)}")
    print()

    print(f"  Information gap (Level 2 - Level 1): {rank_l2 - rank_l1}")
    print()

    # ── Step 2: component contribution analysis ───────────────────────────

    print("=" * 60)
    print("STEP 2: Component contribution to rank at Level 2")
    print()

    # Track momentum groups
    track1_idx = [0, 1, 2]   # track1_px, track1_py, track1_pz
    track2_idx = [4, 5, 6]   # track2_px, track2_py, track2_pz
    charge_idx = [3, 7]      # track1_charge, track2_charge
    vertex_idx = [8, 9, 10]  # vertex_x, vertex_y, vertex_z
    length_idx = [11]        # decay_length

    groups = {
        "track1 momentum (px,py,pz)":  track1_idx,
        "track2 momentum (px,py,pz)":  track2_idx,
        "charges (+1/-1)":             charge_idx,
        "vertex position (x,y,z)":     vertex_idx,
        "decay_length":                length_idx,
        "all track momenta":           track1_idx + track2_idx,
        "momenta + charges":           track1_idx + track2_idx + charge_idx,
        "momenta + vertex":            track1_idx + track2_idx + vertex_idx,
        "full Level 2":                list(range(N_LEVEL2)),
    }

    for name, idx in groups.items():
        r, _, _ = delta_rank_subset(x_level2, edges, idx)
        print(f"  {name:<40} rank = {r}")

    print()

    # ── Step 3: within-event graph analysis ──────────────────────────────

    print("=" * 60)
    print("STEP 3: Within-event graph analysis (Sigma at decay vertex)")
    print()
    print("  Within-event graph per event:")
    print("    Node 0: primary vertex  (production point)")
    print("    Node 1: secondary vertex (decay point)")
    print("    Node 2: positive track endpoint")
    print("    Node 3: negative track endpoint")
    print("    Edges: (0->1), (1->2), (1->3)")
    print("    No particle labels. Direction from observable.")
    print()

    within_results = within_event_analysis(df, n_events=min(1000, n_events))

    print(f"  Events analyzed: {len(within_results)}")
    print()
    print(f"  rank(Im Delta) per event:")
    rank_counts = within_results["rank_delta"].value_counts().sort_index()
    for r, count in rank_counts.items():
        print(f"    rank = {r}: {count} events ({100*count/len(within_results):.1f}%)")
    print()
    print(f"  Sigma antisymmetric term at decay vertex:")
    print(f"    Mean:   {within_results['asym_at_decay_vertex'].mean():.4f}")
    print(f"    Std:    {within_results['asym_at_decay_vertex'].std():.4f}")
    print(f"    Nonzero (> {SVD_TOLERANCE}): "
          f"{(within_results['asym_at_decay_vertex'] > SVD_TOLERANCE).sum()} "
          f"/ {len(within_results)}")
    print()
    print(f"  Momentum closure (should be ~0 for conservation):")
    print(f"    Mean |closure|: {within_results['momentum_closure'].mean():.6f} MeV/c")
    print(f"    Max  |closure|: {within_results['momentum_closure'].max():.6f} MeV/c")
    print()

    # ── Step 4: Particle identification as C projection ───────────────────

    print("=" * 60)
    print("STEP 4: Particle identification as declared C projection")
    print()
    print("  Object frame applies particle identification to Level 2 field:")
    print("  kaon = negative track (track2 in Level 2)")
    print("  pion = positive track (track1 in Level 2)")
    print()
    print("  This is a C projection:")
    print("    Preserves: track identity assignment, enables mass hypothesis")
    print("    Discards:  charge sign as independent observable dimension")
    print("               (charge is absorbed into the identity label)")
    print()

    # Compare rank with and without charge components
    rank_no_charge, _, _ = delta_rank_subset(
        x_level2, edges,
        track1_idx + track2_idx + vertex_idx + length_idx
    )
    rank_with_charge, _, _ = delta_rank_subset(
        x_level2, edges,
        track1_idx + track2_idx + charge_idx + vertex_idx + length_idx
    )

    print(f"  rank without charge components: {rank_no_charge}")
    print(f"  rank with charge components:    {rank_with_charge}")
    print(f"  Charge contribution to rank:    {rank_with_charge - rank_no_charge}")
    print()

    if rank_with_charge > rank_no_charge:
        print("  FINDING: Charge components contribute independent rank.")
        print("  Particle identification (which absorbs charge into label)")
        print("  discards these dimensions from the declared field.")
        print("  The C projection is LOSSY -- it discards information.")
    else:
        print("  FINDING: Charge components are redundant given momenta.")
        print("  Particle identification is a LOSSLESS C projection")
        print("  for rank purposes (though it still discards individual values).")
    print()

    # ── Step 5: Hypothesis assessment ────────────────────────────────────

    print("=" * 60)
    print("STEP 5: Hypothesis assessment")
    print()

    rank_track_momenta, _, _ = delta_rank_subset(
        x_level2, edges, track1_idx + track2_idx)

    h1 = rank_l2 >= rank_l1
    h2 = rank_track_momenta >= 3
    h3 = (within_results["asym_at_decay_vertex"] > SVD_TOLERANCE).all()
    h4_lossy = rank_with_charge > rank_no_charge

    print(f"  H1: rank(Level 2) >= rank(Level 1)")
    print(f"      Level 2: {rank_l2}, Level 1: {rank_l1}")
    print(f"      Result: {'SUPPORTED' if h1 else 'NOT SUPPORTED'}")
    print()
    print(f"  H2: rank(track momenta alone) >= 3")
    print(f"      rank(track momenta): {rank_track_momenta}")
    print(f"      Result: {'SUPPORTED' if h2 else 'NOT SUPPORTED'}")
    print()
    print(f"  H3: Sigma antisymmetric term non-zero at decay vertex")
    print(f"      All events nonzero: {h3}")
    print(f"      Result: {'SUPPORTED' if h3 else 'NOT SUPPORTED'}")
    print()
    print(f"  H4: Particle identification is a lossy C projection")
    print(f"      Charge adds rank: {h4_lossy}")
    print(f"      Result: {'SUPPORTED (lossy)' if h4_lossy else 'NOT SUPPORTED (lossless)'}")
    print()

    # ── Step 6: Object frame comparison ──────────────────────────────────

    print("=" * 60)
    print("STEP 6: Object frame vs relational frame -- information content")
    print()
    print("  Object frame (PDG summary -- Level 0):")
    print(f"    rank(Im Delta) = 3  [established in abr-primary-operators]")
    print(f"    Components: mass, charge, lifetime")
    print(f"    C projections applied: all individual event values discarded")
    print()
    print("  Object frame (Level 1 -- particle-identified per-event):")
    print(f"    rank(Im Delta) = {rank_l1}")
    print(f"    Components: vertex, labeled momenta, mass, time")
    print(f"    C projections applied: mass hypothesis for invariant_mass")
    print()
    print("  Relational frame (Level 2 -- no particle identity):")
    print(f"    rank(Im Delta) = {rank_l2}")
    print(f"    Components: track momenta, charges, vertex, decay length")
    print(f"    C projections applied: charge-sign ordering of tracks")
    print()
    print("  Information progression:")
    print(f"    PDG summary -> Level 1: +{rank_l1 - 3} dimensions recovered")
    print(f"    Level 1 -> Level 2:     +{rank_l2 - rank_l1} dimensions recovered")
    print(f"    Total above PDG:         +{rank_l2 - 3} dimensions")
    print()
    print("  Each step toward raw observables recovers relational structure")
    print("  that the object frame's C projections discard.")

    # ── Write report ──────────────────────────────────────────────────────

    os.makedirs(os.path.dirname(REPORT_PATH), exist_ok=True)
    with open(REPORT_PATH, "w", encoding="utf-8") as f:
        f.write("# Hadronic Phase 0a Report\n\n")
        f.write("**Metatron Dynamics, Inc.** V6. "
                "Bounded over D. No claim beyond D.\n\n")
        f.write("## Declared M\n\n")
        f.write("Level 2: raw track observables, no particle identity.\n\n")
        f.write(f"Components: {LEVEL2_COMPONENTS}\n\n")
        f.write("## Results\n\n")
        f.write(f"PDG summary rank (Level 0): 3\n\n")
        f.write(f"Level 1 rank (particle-identified): {rank_l1}\n\n")
        f.write(f"Level 2 rank (no particle identity): {rank_l2}\n\n")
        f.write(f"Track momenta rank: {rank_track_momenta}\n\n")
        f.write(f"Charge contribution: {rank_with_charge - rank_no_charge}\n\n")
        f.write("## Within-event analysis\n\n")
        f.write(f"Events analyzed: {len(within_results)}\n\n")
        f.write(f"Sigma asym at decay vertex -- "
                f"mean: {within_results['asym_at_decay_vertex'].mean():.4f}, "
                f"std: {within_results['asym_at_decay_vertex'].std():.4f}\n\n")
        f.write(f"Momentum closure -- "
                f"mean: {within_results['momentum_closure'].mean():.6f} MeV/c\n\n")
        f.write("## Hypotheses\n\n")
        f.write(f"H1 (Level 2 >= Level 1): "
                f"{'SUPPORTED' if h1 else 'NOT SUPPORTED'}\n\n")
        f.write(f"H2 (track momenta rank >= 3): "
                f"{'SUPPORTED' if h2 else 'NOT SUPPORTED'}\n\n")
        f.write(f"H3 (Sigma nonzero at decay vertex): "
                f"{'SUPPORTED' if h3 else 'NOT SUPPORTED'}\n\n")
        f.write(f"H4 (particle ID is lossy C projection): "
                f"{'SUPPORTED' if h4_lossy else 'NOT SUPPORTED'}\n\n")
        f.write("## Declared projection\n\n")
        f.write("Preserves: rank per level, Sigma magnitudes, "
                "momentum closure.\n\n")
        f.write("Discards: individual edge values and ordering.\n\n")
        f.write("## Object frame comparison\n\n")
        f.write(f"PDG (Level 0) -> Level 1: +{rank_l1 - 3} dimensions\n\n")
        f.write(f"Level 1 -> Level 2: +{rank_l2 - rank_l1} dimensions\n\n")
        f.write(f"Total above PDG: +{rank_l2 - 3} dimensions\n\n")

    print()
    print(f"Report written: {REPORT_PATH}")
    print("\nBounded over D. No claim beyond D.")


if __name__ == "__main__":
    run()
