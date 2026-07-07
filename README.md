# abr-event-analysis

**Metatron Dynamics, Inc.** V6. Bounded over D. No claim beyond D.

## Purpose

Applies the V6 primary operator kernel (Δ and Σ) to individual observable
events — before any statistical reduction — and compares the operator output
against predictions from quantum mechanics, quantum chromodynamics, and
conventional probability-based analysis.

M is declared over individual detector events, not over statistical summaries.
Every C projection is declared explicitly with stated preserved and discarded
invariants. No probability distribution is used as an input to M.

## Domain

D0 → K⁻π⁺ decay. Individual simulated events with full observable vectors:

    x[v] = (vertex_x, vertex_y, vertex_z,
             kaon_px, kaon_py, kaon_pz,
             pion_px, pion_py, pion_pz,
             decay_length, invariant_mass, decay_time)

12 dimensions per event. No statistical reduction before Δ acts.

## Declared comparison targets

All comparison targets declared before execution, PDG-traceable:

- CP violation parameter |ε| = 2.228e-3 (PDG 2022)
- Proton-neutron magnetic moment ratio = -1.4599 (PDG 2022)
- D0 lifetime τ = 0.4101 ps (PDG 2022, used as generative parameter)
- D0 mass m = 1864.84 MeV/c² (PDG 2022, used as generative parameter)

## Structure

    declaration/        M declaration: observable vector, relations, provenance
    operators/          Δ and Σ operator implementations (Python, V6)
    analysis/           rank(Im Δ), declared edge-image admissibility, Σ output
    comparison/         QM, QCD, and probability comparison — all declared before execution
    data/pdg/           PDG source data (CSV) — summary statistics
    data/events/        Generated per-event data (local only, not in git)
    results/            Generated analysis reports (local only, not in git)
    tests/              Test suite — written before execution
    docs/               Kernel documents (V6)
    kernel/             operators_primary.rs reference copy

## Admissibility

Every quantity must be traceable to an observable through M.
Every C projection must state what it preserves and discards.
No probability distribution is admissible as an input to M.
No statistical summary is admissible as a substitute for individual event data.
Declaration is never a substitute for traceability.

## Kernel

V6 primary kernel: E_primary(x) = Σ(Δ(x))
Grounding: operators_primary.rs, primary_operators_delta_sigma_v6.md

## Origin

Robin Bruce Macomber
Metatron Dynamics, Inc.
Triad role assignment declared per session.
