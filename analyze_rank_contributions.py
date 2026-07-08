"""
analyze_rank_contributions.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

Identifies which observable components contribute the 3 additional
independent dimensions found in Phase 0 per-event analysis.

Method: incremental rank analysis.
Start with the PDG-equivalent component subset.
Add components one at a time.
Record when rank increases.

A rank increase when component C is added means C contributes an
independent direction of relational contrast not present in the
previously declared components.

This is a declared C projection of the Delta output field.
Preserves: which components contribute independent rank dimensions.
Discards: the specific values of those dimensions.

All declared before execution. No result adjusted post-hoc.
"""

import numpy as np
import pandas as pd
import os
import sys

# Add repo root to path
sys.path.insert(0, os.path.dirname(__file__))

from declaration.observable import OBSERVABLE_COMPONENTS, N_COMPONENTS
from declaration.relations import declare_relations, RelationProvenance
from operators.primary import operator_delta
from analysis.rank import im_delta_rank, SVD_TOLERANCE

EVENTS_PATH = os.path.join(
    os.path.dirname(__file__), "data", "events", "d0_events.csv"
)
REPORT_PATH = os.path.join(
    os.path.dirname(__file__), "results", "rank_contributions_report.md"
)

# PDG-equivalent components:
# The PDG summary for D0 preserves mass, lifetime, and charge.
# In our observable vector, the closest equivalents are:
#   invariant_mass  -- corresponds to PDG mass measurement
#   decay_time      -- corresponds to PDG lifetime measurement
#   decay_length    -- related to lifetime via momentum
# These 3 are what conventional analysis retains.
PDG_EQUIVALENT = ["invariant_mass", "decay_time", "decay_length"]

# Remaining components -- what conventional analysis discards
DISCARDED_BY_STATISTICS = [c for c in OBSERVABLE_COMPONENTS
                            if c not in PDG_EQUIVALENT]


def rank_of_subset(x, edges, component_subset):
    """
    Compute rank(Im Delta) restricted to the declared component subset.

    Declared projection:
      Preserves: dimension of span for the declared subset.
      Discards:  components not in the subset.
    """
    # Get indices of declared subset
    indices = [OBSERVABLE_COMPONENTS.index(c) for c in component_subset]

    # Build restricted observable field
    x_sub = {eid: vec[indices] for eid, vec in x.items()}

    # Apply Delta
    field = operator_delta(x_sub, edges)
    rank, sv = im_delta_rank(field)
    return rank, sv


def run():
    print("\nanalyze_rank_contributions.py")
    print("Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.")
    print("=" * 60)
    print()

    if not os.path.exists(EVENTS_PATH):
        print(f"Event data not found: {EVENTS_PATH}")
        print("Run generate_events.py first.")
        sys.exit(1)

    df = pd.read_csv(EVENTS_PATH)
    n_events = len(df)
    print(f"Events loaded: {n_events}")

    # Build observable field
    x = {}
    for _, row in df.iterrows():
        event_id = int(row["event_id"])
        vec = np.array([row[c] for c in OBSERVABLE_COMPONENTS], dtype=float)
        x[event_id] = vec

    # Declare relations
    event_ids = list(x.keys())
    edges = declare_relations(
        event_ids=event_ids,
        relation_type=RelationProvenance.TEMPORAL,
        direction_basis="declared production order"
    )

    print(f"Declared edges: {len(edges)}")
    print()

    # ── Step 1: PDG-equivalent baseline ──────────────────────────────────

    print("=" * 60)
    print("STEP 1: PDG-equivalent components only")
    print(f"Components: {PDG_EQUIVALENT}")
    rank_pdg, _ = rank_of_subset(x, edges, PDG_EQUIVALENT)
    print(f"rank(Im Delta): {rank_pdg}")
    print(f"(Established PDG summary rank: 3)")
    print()

    # ── Step 2: Full 12-component field ──────────────────────────────────

    print("=" * 60)
    print("STEP 2: Full 12-component observable field")
    print(f"Components: {OBSERVABLE_COMPONENTS}")
    rank_full, _ = rank_of_subset(x, edges, OBSERVABLE_COMPONENTS)
    print(f"rank(Im Delta): {rank_full}")
    print(f"Information gap: {rank_full - rank_pdg} dimensions")
    print()

    # ── Step 3: Incremental rank analysis ────────────────────────────────

    print("=" * 60)
    print("STEP 3: Incremental rank analysis")
    print("Adding discarded components one at a time.")
    print("A rank increase signals an independent relational direction.")
    print()

    current_components = list(PDG_EQUIVALENT)
    current_rank = rank_pdg
    contributions = []

    print(f"  Starting rank: {current_rank}")
    print(f"  Starting components: {current_components}")
    print()

    for component in DISCARDED_BY_STATISTICS:
        test_components = current_components + [component]
        test_rank, _ = rank_of_subset(x, edges, test_components)
        delta_rank = test_rank - current_rank

        status = "NEW DIMENSION" if delta_rank > 0 else "redundant"
        print(f"  + {component:<20} rank: {current_rank} -> {test_rank}  [{status}]")

        contributions.append({
            "component": component,
            "rank_before": current_rank,
            "rank_after": test_rank,
            "rank_increase": delta_rank,
            "contributes_new_dimension": delta_rank > 0,
        })

        if delta_rank > 0:
            current_components.append(component)
            current_rank = test_rank

    print()
    new_dims = [c for c in contributions if c["contributes_new_dimension"]]
    redundant = [c for c in contributions if not c["contributes_new_dimension"]]

    print("=" * 60)
    print("RESULT: Components contributing new independent dimensions:")
    for c in new_dims:
        print(f"  {c['component']}")
    print()
    print("Components redundant given prior declared components:")
    for c in redundant:
        print(f"  {c['component']}")
    print()

    # ── Step 4: Group analysis ────────────────────────────────────────────

    print("=" * 60)
    print("STEP 4: Group rank analysis")
    print()

    groups = {
        "decay_vertex_position": ["vertex_x", "vertex_y", "vertex_z"],
        "kaon_momentum":         ["kaon_px", "kaon_py", "kaon_pz"],
        "pion_momentum":         ["pion_px", "pion_py", "pion_pz"],
        "pdg_equivalent":        PDG_EQUIVALENT,
    }

    for group_name, group_components in groups.items():
        rank_group, _ = rank_of_subset(x, edges, group_components)
        print(f"  {group_name:<30} rank: {rank_group}")

    print()

    # Vertex + PDG
    combined = groups["decay_vertex_position"] + PDG_EQUIVALENT
    rank_combined, _ = rank_of_subset(x, edges, combined)
    print(f"  vertex + pdg_equivalent            rank: {rank_combined}")

    # Momenta + PDG
    combined2 = groups["kaon_momentum"] + groups["pion_momentum"] + PDG_EQUIVALENT
    rank_combined2, _ = rank_of_subset(x, edges, combined2)
    print(f"  momenta + pdg_equivalent           rank: {rank_combined2}")

    print()
    print("=" * 60)
    print("INTERPRETATION (declared projection):")
    print(f"  PDG summary preserves: {rank_pdg} independent dimension(s)")
    print(f"  Full per-event field:  {rank_full} independent dimensions")
    print(f"  Additional dimensions: {rank_full - rank_pdg}")
    print()
    print("  What statistical reduction discards:")
    for c in new_dims:
        print(f"    {c['component']} -- opens 1 new independent direction")
    print()

    # ── Write report ──────────────────────────────────────────────────────

    os.makedirs(os.path.dirname(REPORT_PATH), exist_ok=True)
    with open(REPORT_PATH, "w", encoding="utf-8") as f:
        f.write("# Rank Contribution Analysis\n\n")
        f.write("**Metatron Dynamics, Inc.** V6. Bounded over D. No claim beyond D.\n\n")
        f.write("## Declaration\n\n")
        f.write(f"n_events: {n_events}\n\n")
        f.write(f"PDG-equivalent components: {PDG_EQUIVALENT}\n\n")
        f.write(f"Full observable vector: {OBSERVABLE_COMPONENTS}\n\n")
        f.write("## Results\n\n")
        f.write(f"PDG-equivalent rank: {rank_pdg}\n\n")
        f.write(f"Full per-event rank: {rank_full}\n\n")
        f.write(f"Information gap: {rank_full - rank_pdg} dimensions\n\n")
        f.write("## Components contributing new independent dimensions\n\n")
        for c in new_dims:
            f.write(f"- {c['component']}: rank {c['rank_before']} -> {c['rank_after']}\n")
        f.write("\n")
        f.write("## Components redundant given prior declared components\n\n")
        for c in redundant:
            f.write(f"- {c['component']}\n")
        f.write("\n")
        f.write("## Group rank analysis\n\n")
        for group_name, group_components in groups.items():
            rank_g, _ = rank_of_subset(x, edges, group_components)
            f.write(f"- {group_name}: rank {rank_g}\n")
        f.write("\n")
        f.write("## Declared projection\n\n")
        f.write("Preserves: which components contribute independent rank dimensions.\n\n")
        f.write("Discards: the specific values of those dimensions and their ordering.\n\n")

    print(f"Report written: {REPORT_PATH}")
    print("\nBounded over D. No claim beyond D.")


if __name__ == "__main__":
    run()
