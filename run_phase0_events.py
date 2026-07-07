"""
run_phase0_events.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

Phase 0: run Delta and Sigma over individual D0 -> K- pi+ events.

Declared hypothesis:
  rank(Im Delta) over individual events > rank(Im Delta) over PDG summary.
  PDG summary rank = 3 (established in abr-primary-operators Phase 0).
  Per-event rank measures what statistical reduction discards.

Admissibility checks run before any operator acts.
All comparison targets declared before execution.
No result adjusted post-hoc.
"""

import numpy as np
import pandas as pd
import os
import sys

from declaration.observable import OBSERVABLE_COMPONENTS, N_COMPONENTS
from declaration.relations import declare_relations, RelationProvenance
from declaration.targets import DECLARED_TARGETS
from operators import operator_e_primary
from analysis import (
    im_delta_rank, im_sigma_rank,
    declared_edge_image_admissibility_check,
    expression_condition, detect_failure_mode,
    rho_p_ratio,
)

EVENTS_PATH = os.path.join(
    os.path.dirname(__file__), "data", "events", "d0_events.csv"
)
REPORT_PATH = os.path.join(
    os.path.dirname(__file__), "results", "phase0_events_report.md"
)

# Declared before execution
PDG_RANK = 3  # rank(Im Delta) over PDG summary statistics — Phase 0 result


def run():
    print("\nrun_phase0_events.py")
    print("Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.")
    print("=" * 60)
    print(f"Observable components ({N_COMPONENTS}): {OBSERVABLE_COMPONENTS}")
    print(f"PDG summary rank (established): {PDG_RANK}")
    print()

    if not os.path.exists(EVENTS_PATH):
        print(f"Event data not found: {EVENTS_PATH}")
        print("Run generate_events.py first.")
        sys.exit(1)

    df = pd.read_csv(EVENTS_PATH)
    n_events = len(df)
    print(f"Events loaded: {n_events}")

    # Build observable field: event_id -> observable vector
    x = {}
    for _, row in df.iterrows():
        event_id = int(row["event_id"])
        vec = np.array([row[c] for c in OBSERVABLE_COMPONENTS], dtype=float)
        x[event_id] = vec

    # Declare relations — temporal, strictly directed
    event_ids = list(x.keys())
    edges = declare_relations(
        event_ids=event_ids,
        relation_type=RelationProvenance.TEMPORAL,
        direction_basis=(
            "Event production order — earlier event index precedes later "
            "event index in declared generation sequence. Direction is "
            "determined by the declared generative order, not by a choice."
        )
    )

    print(f"Declared edges: {len(edges)}")
    print(f"Relation type: temporal (strictly directed)")
    print()

    # Run operators
    print("Running Delta and Sigma...")
    delta_field, sigma_field = operator_e_primary(x, edges)
    print(f"Delta field shape: {delta_field.shape}")

    # Analysis
    rank_delta, sv_delta = im_delta_rank(delta_field)
    rank_sigma, sv_sigma = im_sigma_rank(sigma_field)

    admissible, witness = declared_edge_image_admissibility_check(delta_field)
    expressed, _, has_asym = expression_condition(
        delta_field, sigma_field, edges
    )
    mode = detect_failure_mode(delta_field, sigma_field, edges)
    rho_p = rho_p_ratio(sigma_field, edges)

    # Results
    print(f"\nrank(Im Delta): {rank_delta}")
    print(f"rank(Im Sigma): {rank_sigma}")
    print(f"Declared edge-image admissible: {admissible}")
    print(f"Expression condition: {expressed}")
    print(f"Failure mode: {mode.value}")
    print(f"rho_P: {rho_p:.4f}")
    print()

    # Hypothesis
    gap = rank_delta - PDG_RANK
    print(f"PDG summary rank:    {PDG_RANK}")
    print(f"Per-event rank:      {rank_delta}")
    print(f"Information gap:     {gap} dimensions")
    print()
    if gap > 0:
        print(f"SUPPORTED: per-event rank ({rank_delta}) > PDG rank ({PDG_RANK})")
        print(f"Statistical reduction discarded {gap} independent dimensions.")
    elif gap == 0:
        print(f"NOT SUPPORTED: per-event rank = PDG rank = {PDG_RANK}")
        print(f"No additional relational information recovered from individual events.")
    else:
        print(f"ANOMALOUS: per-event rank ({rank_delta}) < PDG rank ({PDG_RANK})")
        print(f"Inspect declaration.")

    # Write report
    os.makedirs(os.path.dirname(REPORT_PATH), exist_ok=True)
    with open(REPORT_PATH, "w", encoding="utf-8") as f:
        f.write("# Phase 0 — Per-Event Analysis Report\n\n")
        f.write("**Metatron Dynamics, Inc.** V6. Bounded over D. No claim beyond D.\n\n")
        f.write("## Declaration\n\n")
        f.write(f"Observable vector: {OBSERVABLE_COMPONENTS}\n\n")
        f.write(f"n_events: {n_events}\n\n")
        f.write(f"Relation type: temporal, strictly directed\n\n")
        f.write("## Operator outputs\n\n")
        f.write(f"rank(Im Delta): {rank_delta}\n\n")
        f.write(f"rank(Im Sigma): {rank_sigma}\n\n")
        f.write(f"Declared edge-image admissible: {admissible}\n\n")
        f.write(f"Expression: {expressed}\n\n")
        f.write(f"Failure mode: {mode.value}\n\n")
        f.write(f"rho_P: {rho_p:.4f}\n\n")
        f.write("## Hypothesis\n\n")
        f.write(f"PDG summary rank (established): {PDG_RANK}\n\n")
        f.write(f"Per-event rank: {rank_delta}\n\n")
        f.write(f"Information gap: {gap} dimensions\n\n")
        f.write(f"Result: {'SUPPORTED' if gap > 0 else 'NOT SUPPORTED' if gap == 0 else 'ANOMALOUS'}\n\n")
        f.write("## Comparison targets (declared before execution)\n\n")
        for name, val in DECLARED_TARGETS.items():
            f.write(f"- {name}: {val}\n")
        f.write("\n")

    print(f"Report written: {REPORT_PATH}")
    print("\nBounded over D. No claim beyond D.")


if __name__ == "__main__":
    run()
