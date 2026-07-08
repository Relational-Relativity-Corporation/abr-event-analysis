# The Mathematical Interior of the Primary Region

**Metatron Dynamics, Inc.** V6. Bounded over D. No claim beyond D.

**Grounding:** All formalism in this document is derived from the V6 primary
kernel (`operators_primary.rs`), the V5 ABR kernel definitions in `operators_v5.rs`,
and `operators_notation_and_constraint_v6.md`. No new operators are introduced.
No new primitives are introduced. Every quantity is traceable to the V6 declarations
or to a declared projection of them.

---

## Purpose

The Primary Region vocabulary (`canonical_primary_region_vocabulary.md`) and the
relational regions framework (`relational_regions_formal_boundaries.md`) define
the Primary Region descriptively: the space of minimally admissible relational
configurations, where persistence is not assumed and expression is contingent.

This document makes that description precise in operator terms. It derives, from
the V6 operator definitions, the following:

1. The **marginal admissibility condition** — when Δ is at its minimum non-trivial
   output
2. The **ρ_P ratio** — defined as a computable operator invariant
3. The **expression failure modes** — the operator conditions under which expression
   does not appear
4. The **relational direction unification** — V6 finding that spatial and
   persistence domains are subject to the same admissibility condition on direction
5. The **boundary cancellation asymmetry** — V5 finding with consequences for
   Primary Region persistence structure, now understood as the general case of
   directional admissibility

No temporal language is used. All conditions are structural, not causal.

---

## 1. Marginal Admissibility

### 1.1 The Image of Δ

Let M be declared and let x = M(o) be the observable field over declared relations
in D. The operator Δ produces, for each declared edge e = (s, t):

    Δ(x)[e] = x[s] − x[t]

The **declared edge-image of Δ** is the set of vectors {Δ(x)[e]} over all declared
edges. Im(Δ) is the span of {Δ(x)[e]} over the declared EdgeField, as a vector
space over ℝ. The **rank of Im(Δ)** — written rank(Im Δ) — is the dimension of
that span.

*Wording:* rank(Im Δ) is the dimension of the span — not the count of linearly
independent directed differences. Those are equivalent only when every edge vector
is linearly independent of the others. The correct statement in all cases is:
dimension of the span.

This is a declared projection of Δ's output — it discards the individual values
of Δ(x)[e] and retains the dimension of their span. What it preserves and what
it discards are stated here.

Under the V6 directional admissibility condition, the declared edge-image of any
admissible declared structure is asymmetric: no declared edge vector has its
negation as another declared edge vector. A symmetric declared edge-image signals
inadmissible structure — a reverse edge declared without satisfying the distinctness
axiom and independent provenance requirement.

### 1.2 Marginal Admissibility Condition

**Marginal admissibility** is the condition:

    rank(Im Δ) = 1

The span of the declared edge vectors is one-dimensional. The observable field is
not uniform (rank zero is excluded) but produces only a single direction of
relational contrast across the entire declared structure.

Below this — rank(Im Δ) = 0 — Δ collapses: all declared differences are zero,
and Σ(Δ(x)) = 0 identically. Expression is undefined over the observable. This
is **differentiation collapse**.

At rank 1, Σ acts on a one-dimensional image. Expression is contingent on whether
Σ produces a non-zero antisymmetric component from that single direction. This
contingency — expression may appear or fail — is the operator grounding of the
Primary Region vocabulary statement that expression is not presumed.

### 1.3 Primary Region Boundary

The Primary Region is the interval:

    1 ≤ rank(Im Δ) < threshold(ρ_P → 1)

where the lower bound is differentiation collapse and the upper bound is the
Atomic/Molecular transition defined in Section 2.4. Both bounds are operator-
derived quantities, not assumed scales.

### 1.4 Cold-Start Grounding

The V5 implementation provides a direct instantiation of marginal admissibility
in the persistence domain. On cold start, A_persistence = E_current (E_prior =
zero EdgeField). The entire current relational state is declared new — no prior
distinguishable difference exists. This is the persistence-domain analogue of
rank(Im Δ) = 1: exactly one relational step of distinguishable structure is
present, with no accumulation from prior cycles.

Note: cold start does not prove rank(Im Δ) = 1. It is the persistence analogue
of the marginal condition — both share the structural property of minimal prior
distinguishable structure, not identical operator definitions.

Implementation grounding: `v5_cold_start_a_persistence_equals_current`.

---

## 2. The ρ_P Ratio as an Operator Invariant

### 2.1 Descriptive Origin

The relational regions document defines:

    ρ_P = D_r / C_x  ≪  1

where D_r is described as "distinguishable relational reconfigurations" and C_x
as "relational complexity." Both are now given operator definitions computable
from V6 outputs without introducing new primitives.

### 2.2 D_r — Dimension of Σ Output Span

    D_r := rank(Im Σ)

Im Σ is the span of Σ outputs over the declared relational structure, as a vector
space over ℝ. D_r is the dimension of that span — how many linearly independent
directions the circulation operator produces over the declared relations.

Σ is defined as:

    Σ(g)[e] = g[e] + ρ[e] · ( Σ_{adj⁺(e)} g − Σ_{adj⁻(e)} g )

The output has two components: the pass-through term g[e] and the antisymmetric
term ρ[e] · (adj⁺ sum − adj⁻ sum). rank(Im Σ) counts the dimension of the
combined output space. rank(Im Σ) > 0 does not by itself guarantee a non-zero
antisymmetric component — the pass-through term alone can produce non-zero Σ
output while the antisymmetric term is zero. This distinction is relevant to the
expression threshold in Section 4.

Under the V6 directional admissibility condition, admissible declared structures
are asymmetric. For any admissible structure with non-trivial adjacency (not FM2),
the antisymmetric term will be non-zero, and D_r > 0 follows from the asymmetry
of the declared structure rather than from cancellation.

### 2.3 C_x — Propagation Capacity

    C_x := |{ e ∈ declared edges : adj⁺(e) ≠ ∅ }|

C_x is the count of declared edges with at least one declared forward neighbor —
edges that have somewhere to propagate their value in the Σ adjacency sum. This
is directly computable from DeclaredRelations: for each edge e, adj⁺(e) is the
list of edges whose source node equals the target node of e.

C_x measures the propagation capacity of the declared relational structure. A
large C_x relative to D_r means the declared structure offers many propagation
opportunities but produces few distinguishable circulation outputs — the regime
of the Primary Region.

Implementation grounding: `v5_b_terminal_accumulates_nothing`.

### 2.4 The Ratio

    ρ_P = rank(Im Σ) / |{ e : adj⁺(e) ≠ ∅ }|

In the Primary Region, rank(Im Σ) is small (often 1, at marginal admissibility),
while C_x is large relative to the configuration. ρ_P ≪ 1 is a derivable
consequence of the operator structure over declared relations.

The Atomic/Molecular transition — ρ_P ≈ 1 — occurs when rank(Im Σ) grows
commensurate with C_x. At this point B becomes admissible, Σ transitions to R
receiving B(Δ(x)), and the full V5 ABR kernel activates. The exact operator
condition at which this transition occurs is an open condition (Section 6, item 1).

### 2.5 Constraint

Both D_r and C_x are derived from declared relations with provenance traceable to
an observable through M. A ρ_P computed from undeclared or assumed relations is
fabricated within D.

---

## 3. Expression Failure Modes

### 3.1 Overview

The Primary Region vocabulary identifies boundary failure as legitimate
observable data — not error but result. Two structural failure modes are legitimate
outcomes of admissible declared structures. A third condition — Circulation
Cancellation — is now characterized as an inadmissibility signal rather than a
legitimate failure mode.

**V6 update:** Under the directional admissibility condition, Failure Mode 3
(Circulation Cancellation) requires a symmetrically balanced adjacency structure.
Such a structure requires inadmissible reverse edges — edges declared without
satisfying the distinctness axiom and independent provenance requirement. FM3 is
therefore retained as a named condition for detection purposes, but it signals a
declaration error rather than a legitimate operator outcome of an admissible
declared structure.

### 3.2 Failure Mode 1 — Differentiation Collapse

**Operator condition:**

    rank(Im Δ) = 0

All declared edges produce zero directed difference. The observable field is
uniform across all declared relations at the declared resolution.

**Consequence:** Σ(Δ(x)) = 0 identically. ρ = 0 everywhere. Expression is
undefined — the composition has no non-trivial input.

**Vocabulary mapping:** *differentiation collapse* — Δ fails.

**Distinguishing note:** Differentiation collapse is not the same as the observable
being zero. A uniform non-zero field collapses Δ equally: Δ(x)[e] = c − c = 0
for any constant c. Collapse is a property of the relational structure of the
observable across declared edges, not of its magnitude.

**Implementation grounding:** `v5_stable_field_a_persistence_zero`.

### 3.3 Failure Mode 2 — Relational Isolation

**Operator condition:**

    adj⁺(e) = ∅  and  adj⁻(e) = ∅  for all declared edges e

Every edge is an island — no declared predecessor, no declared successor. The
declared relational structure has no interior.

**Consequence:** The antisymmetric term of Σ at every edge reduces to:

    ρ[e] · ( Σ_{adj⁺(e)} g − Σ_{adj⁻(e)} g ) = ρ[e] · (0 − 0) = 0

Σ reduces to the identity on Δ(x). No local antisymmetry is detectable.

**Why the condition requires both adj⁺ and adj⁻ empty:** When adj⁺(e) = ∅ and
adj⁻(e) ≠ ∅, the antisymmetric term becomes ρ[e] · (0 − Σ_{adj⁻} g), which is
generally non-zero. The predecessor accumulation still produces an asymmetric
contribution. Full isolation requires both conditions.

**Vocabulary mapping:** *circulation insufficiency* — Σ has no relational
structure to circulate.

**Implementation grounding:** `v5_b_terminal_accumulates_nothing`;
`pred_nonempty_succ_empty_produces_nonzero_antisymmetric_term` (V6 — closes
the open condition from V3 §6 item 1: pred ≠ ∅, succ = ∅ produces non-zero
antisymmetric term, confirming full isolation requires both conditions).

### 3.4 Circulation Cancellation — Inadmissibility Condition

**Condition:**

    ρ[e] · ( Σ_{adj⁺(e)} Δ(x) − Σ_{adj⁻(e)} Δ(x) ) = 0  for all e
    rank(Im Δ) > 0

The antisymmetric term of Σ is zero for all declared edges while Δ(x) is
non-trivially non-zero.

**V6 characterization:** This condition requires a symmetrically balanced
adjacency structure — one where forward and backward neighbor sums cancel at
every edge. Such a structure requires declared reverse edges. Under the V6
directional admissibility condition, reverse edges are admissible only when each
direction has independent observable provenance and satisfies the distinctness
axiom — in which case they are distinct relations, not symmetric pairs, and
cancellation is not structurally guaranteed.

When Circulation Cancellation is detected in an operator output, it signals that
the declared structure contains inadmissible reverse edges — edges whose negation
was declared without independent provenance. It is a declaration error to be
resolved, not a legitimate operator outcome.

**Retained as a named condition** for detection purposes. The implementation
(`detect_failure_mode()` in `operators_primary.rs`) detects it and returns
`FailureMode::CirculationCancellation` — indicating inadmissible structure.

**Vocabulary mapping:** *coherence failure* — retained as vocabulary term, now
understood as a signal of declaration error rather than a legitimate boundary
failure of an admissible declared structure.

**Relationship to the spatial boundary note:** The V5 boundary cancellation
theorem establishes that B cancels at the boundary of an open declared relational
structure, holding by the Lemma A(t,s) = −A(s,t) for node-valued spatial loci.
That theorem is not Circulation Cancellation — it is a structural property of
open boundaries. Circulation Cancellation is an interior condition. The V6
directional admissibility condition clarifies why the Lemma holds for spatial
node-valued loci (both directions share node values) but not for persistence
edge-valued loci (each direction is independently computed). See Section 5.

---

## 4. Expression Threshold

**Definition:** Expression is non-trivial antisymmetric circulation output of Σ
— not merely non-zero Σ output.

Expression is admissible when the following two conditions hold jointly:

1. rank(Im Σ) > 0 — the Σ output span is non-trivial
2. The antisymmetric term of Σ is non-zero on at least one declared edge

**The conditions are non-redundant.** Condition 1 is satisfied by the pass-through
term g[e] alone. Condition 2 is the additional requirement distinguishing
expression from mere pass-through.

Failure Mode 1 (differentiation collapse) is excluded by Condition 1.
Failure Mode 2 (relational isolation) is excluded by Condition 2.
Circulation Cancellation (inadmissibility signal) is excluded by Condition 2.

**V6 note:** Under the directional admissibility condition, Condition 2 is always
satisfied for admissible declared structures that are not in FM1 or FM2. An
admissible declared structure is asymmetric; the antisymmetric term is non-zero
wherever adjacency is declared. The two conditions remain stated separately to
preserve the distinction between presence (Condition 1) and asymmetric
organization (Condition 2).

---

## 5. Relational Direction Unification

**V6 finding:** The spatial and persistence domains are subject to the same
admissibility condition on relational direction.

Every declared relation has exactly one admissible direction — the direction
traceable to an observable through M. The reverse direction requires independent
observable provenance. If it does not have it, it is inadmissible.

The V5 finding that A_persistence[(s,t)] ≠ −A_persistence[(t,s)] — and that the
boundary cancellation theorem does not transfer to the persistence domain — is not
a special property of the persistence domain. It is the general property of all
declared relations becoming visible where directionality is explicit.

The Lemma A(t,s) = −A(s,t) holds for spatial node-valued loci because both
directions happen to share the same node values: A(s,t) = x[s] − x[t] and
A(t,s) = x[t] − x[s] = −A(s,t). This algebraic relationship is a consequence
of the shared node values, not a property of observables in general. For
persistence edge-valued loci, E_current[(s,t)] and E_current[(t,s)] are
independently computed by the R operator — no algebraic relationship exists,
and the Lemma does not hold.

The V6 directional admissibility condition generalizes this: in both domains,
the reverse direction is inadmissible unless it has independent observable
provenance. The spatial Lemma does not make the reverse direction admissible —
it is an algebraic identity between two inadmissible declarations that happen
to cancel, not a property of admissible structure.

Implementation grounding: `v5_boundary_cancellation_does_not_transfer_to_persistence`;
`evolution_direction_is_current_minus_prior`;
`bidirectional_pair_is_inadmissible` (all in V6 `operators_primary.rs`).

---

## 6. Correspondence Table

| Vocabulary Term | Operator Condition | V6 Grounding |
|---|---|---|
| differentiable | rank(Im Δ) ≥ 1 | `v5_changing_field_a_persistence_nonzero` |
| admissible | rank(Im Δ) ≥ 1 ∧ at least one non-isolated edge ∧ declared edge-image asymmetric | Open DAG and directed chain structures in V6 tests |
| marginal persistence | rank(Im Δ) = 1 | Cold-start analogy; `v5_cold_start_a_persistence_equals_current` |
| expression threshold | rank(Im Σ) > 0 ∧ antisymmetric term non-zero on ≥ 1 edge | Section 4; conditions proved non-redundant |
| differentiation collapse | rank(Im Δ) = 0 | `v5_stable_field_a_persistence_zero` |
| circulation insufficiency | adj⁺(e) = ∅ ∧ adj⁻(e) = ∅ for all e (Failure Mode 2) | `v5_b_terminal_accumulates_nothing`; `pred_nonempty_succ_empty_produces_nonzero_antisymmetric_term` |
| coherence failure | antisymmetric term = 0 for all e; rank(Im Δ) > 0 — inadmissibility signal in V6 | `circulation_cancellation_is_inadmissibility_signal` |
| expression dropout | FM1 or FM2 active | Implied by Sections 3.2–3.3 |
| ρ_P ≪ 1 | rank(Im Σ) / \|{e : adj⁺(e) ≠ ∅}\| ≪ 1 | Computable from DeclaredRelations and Σ output |
| ρ_P ≈ 1 | rank(Im Σ) / \|{e : adj⁺(e) ≠ ∅}\| ≈ 1 | Transition condition; threshold value open (Section 7, item 1) |

---

## 7. Open Conditions

1. **Threshold characterization of ρ_P → 1.** The exact operator condition at
   which the Primary Region transitions to the Atomic/Molecular Region is not yet
   stated as a computable criterion. ρ_P is computable; the transition threshold
   is not yet derived. No specific numerical value of ρ_P should be treated as a
   declared transition threshold.

2. **Expression dropout distinction.** The vocabulary distinguishes expression
   dropout from coherence failure. Under V6, coherence failure is recharacterized
   as an inadmissibility signal rather than a legitimate failure mode. Whether
   expression dropout has a residual distinct operator definition — beyond the two
   legitimate failure modes (FM1, FM2) — is open.

3. **Independence of Failure Modes 1 and 2.** Both produce zero antisymmetric Σ
   output by structurally different conditions. Whether they span distinct
   equivalence classes of relational failure under a declared projection is not
   yet demonstrated.

4. **Stage 2 of the hadronic test.** The derivation that the magnetic moment ratio
   is a function of Δ(x) over the declared hadronic observable field has not been
   performed. See `primary_operators_delta_sigma_v6.md` §10.

5. **Persistence domain ρ_P.** The ρ_P ratio is defined for the primary kernel
   (Σ over spatial relations). Whether an analogous ratio over the persistence
   phase (Phase 2 of the ABR kernel) is admissible — and what it would characterize
   — is open. The relational direction unification in Section 5 suggests the
   persistence-domain ratio would not reduce to the spatial ratio by simple
   substitution.

6. **Span symmetry vs declared edge-image symmetry.** The current implementation
   checks declared edge-image admissibility (edge-set symmetry). Whether span
   symmetry — whether −v ∈ span(declared edge vectors) — provides additional
   analytical power is not yet determined.

---

## 8. Structural Observation

The document revolves around three operator-derived quantities:

    rank(Im Δ)           — dimension of the span of declared edge vectors
    rank(Im Σ)           — dimension of the span of circulation outputs
    |{e : adj⁺(e) ≠ ∅}| — propagation capacity

Together these three quantities characterize most of the Primary Region structure
defined in this document. A more compact formulation — possibly a single invariant
derived from these three — may be available without introducing new primitives.
This is a structural observation, not a claim or a proposed formulation.

---

## 9. Change Log

**V3 → V6**

| Section | Change | Reason |
|---|---|---|
| Header | Version V3 → V6. Grounding updated to V6 primary kernel. Role references removed from status and closing. | V6 release convention. |
| 1.1 | rank(Im A) → rank(Im Δ). Wording corrected: "dimension of the span" replaces "count of linearly independent directed differences." Note on declared edge-image asymmetry under V6 admissibility added. | V6 primary kernel uses Δ; rank wording corrected — dimension of span is the correct statement in all cases. |
| 1.2 | Marginal admissibility restated in terms of Δ and Σ. | Consistent with V6 primary kernel. |
| 2.2 | D_r redefined as rank(Im Σ) rather than rank(Im R). Wording corrected: "dimension of the span" replaces "linearly independent outputs." V6 note on admissible structures added. | Primary kernel uses Σ; rank wording corrected. |
| 2.3 | C_x redefined using adj⁺(e) rather than succ(e) — consistent with primary kernel adjacency notation. | Notation alignment with V6. |
| 3.1 | Overview updated: FM3 recharacterized as inadmissibility signal, not a legitimate failure mode of an admissible declared structure. | Position B: FM3 requires inadmissible reverse edges. |
| 3.3 | FM2 grounding updated: pred gap (V3 open condition 1) closed by V6 test `pred_nonempty_succ_empty_produces_nonzero_antisymmetric_term`. | Open condition closed in V6 implementation. |
| 3.4 | FM3 renamed "Circulation Cancellation — Inadmissibility Condition." V6 characterization added: signals declaration error, not legitimate operator outcome. Retained as named detection condition. | Position B and V6 directional admissibility. |
| 4 | Expression threshold: V6 note added — Condition 2 always satisfied for admissible structures not in FM1/FM2. | Consequence of Position B. |
| 5 (new) | Relational Direction Unification: V6 finding that spatial and persistence domains are subject to the same admissibility condition. Algebraic explanation of why the Lemma holds spatially but not persistently. | V6 Position B; generalizes V5 boundary cancellation finding. |
| 6 | Correspondence table updated: FM3 → inadmissibility signal. FM2 grounding updated. Admissibility condition updated to include asymmetric declared edge-image. | V6 changes throughout. |
| 7 | Open conditions renumbered. V3 item 1 (pred gap) closed. Items reordered. New item 6 (span vs edge-set symmetry). | Pred gap closed; new open condition from V6. |
| 8 | Structural observation: role reference removed. Wording updated to match V6 invariant names. | Role references removed; notation consistency. |
| 9 | Change log updated V3 → V6. Role references removed throughout. | V6 release. |

---

## Closing Statement

> The Primary Region is the interval 1 ≤ rank(Im Δ) < threshold(ρ_P → 1),
> where rank(Im Δ) is the dimension of the span of declared edge vectors from
> the V6 Δ operator over the declared EdgeField; ρ_P = rank(Im Σ) / |{e : adj⁺(e) ≠ ∅}|
> is derived from V6 Σ operator outputs; and expression is admissible when
> rank(Im Σ) > 0 with a non-zero antisymmetric term — two conditions proved
> non-redundant. Two structural failure modes — differentiation collapse and
> relational isolation — are legitimate outcomes of admissible declared structures.
> Circulation cancellation is an inadmissibility signal rather than a legitimate
> failure mode. The spatial and persistence domains are unified under one
> admissibility condition on relational direction. Six open conditions are declared.
> Three operator-derived invariants characterize the full Primary Region structure.
> Bounded over D. No claim beyond D.

---

*Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.*
