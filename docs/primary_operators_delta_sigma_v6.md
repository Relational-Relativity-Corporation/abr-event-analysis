# The Primary Operators: Δ and Σ

**Metatron Dynamics, Inc.** V6. Bounded over D. No claim beyond D.

**Status:** V6 release.
**Changes from prior version:** See Section 13.
Relational evolution direction stated as a constraint (Section 3.2). Declared
edge-image distinguished from span (Section 3). rank(Im Δ) wording corrected
throughout. FM3 condition aligned to formal document. ρ_P ≈ 1 marked as open.
**Change log:** See Section 13.

**Grounding:** All formalism in this document is derived from the V6 primary
kernel (`operators_primary.rs`), the V5 ABR kernel definitions in `operators_v5.rs`
and `operators_notation_and_constraint_v6.md`, and the Primary Region formalism
(`primary_region_formal_interior_v6.md`). No new primitives are
introduced beyond the renaming of A to Δ and the identification of Σ as R applied
directly to Δ(x) without prior accumulation. Every quantity is traceable to a
declared observable through M.

---

## 1. Declaration

Two questions are irreducible at the smallest observable scale:

1. **Does anything distinguishable exist across the declared relations?**
2. **Is what exists symmetrically or asymmetrically organized?**

No smaller question can be asked of a declared relational structure without
presupposing an answer to one of these. Every operator in every scale of the
framework — including the full V6 kernel chain — is built on top of these two questions.

The operators that ask them are **Δ** and **Σ**.

---

## 2. Domain and Measurement

**D** := { x ∈ ℝⁿ | n < ∞, |x[i]| < ∞ ∀ i }

All quantifiers bounded over D. No claim beyond D.

**M : O → D** declared by Origin before any operator acts. The operators act on
M(o) only — on what M produces from the observable, and nothing else. There is
nothing to evaluate before M is declared.

A statement is admissible within D if and only if every quantity in it can be
traced to an observable through M. Declaration is the act of stating that
traceability and taking responsibility for it. Declaration is never a substitute
for traceability. A quantity that is declared but not traceable to an observable
through M is not inadmissible — it is not about anything in D.

**Declared relational structure:** a finite set of nodes V, each a declared
observable locus, and a finite set of directed edges E ⊆ V × V, each a declared
relation traceable to an observable through M. No assumption of connectivity,
interior, or path structure beyond what is declared.

---

## 3. Operator Δ — Directed Difference

    Δ(x)[e] = x[s] − x[t]    for each declared edge e = (s, t)

Δ takes the directed difference of the observable field across each declared
relation. It produces an edge-valued field over E.

**Constraint:** take the directed difference and nothing else. No relation and no
direction may be added that the declaration did not trace to an observable through
M. No transformation that alters pairwise differences may precede Δ unless
declared within M. Admissible before Δ: uniform shift (x + c), declared unit
scale (x / s). All else requires declaration with preserved and discarded
invariants stated.

**Relationship to V5 A:** Δ is A under a declared renaming. The operator formula
is identical. The renaming reflects the role of Δ as the irreducible primitive
of the primary operator set — making explicit what was implicit in V5: that the
directed difference is the most primitive admissible operation on any declared
relational structure.

**What Δ produces:**

The declared edge-image of Δ — the set of vectors {Δ(x)[e] : e ∈ E} over the
declared EdgeField. This is a declared projection of Δ's output: it discards the
individual edge indices and retains the vectors as a set. What it preserves and
what it discards are stated here.

Two quantities are immediately computable from the declared edge-image:

**rank(Im Δ)** — the dimension of the span of the declared edge vectors. This
is the dimension of the vector space spanned by {Δ(x)[e] : e ∈ E} over ℝ.

*Wording note:* rank(Im Δ) is the dimension of the
span — not "the count of linearly independent directed differences." Those are
equivalent only when every edge vector is linearly independent of the others.
The correct statement in all cases is: dimension of the span.

- rank(Im Δ) = 0: the observable is undifferentiated across all declared
  relations. No relational contrast exists at the declared resolution.
- rank(Im Δ) = 1: the span is one-dimensional. Marginal admissibility.
- rank(Im Δ) > 1: the span has multiple independent directions of contrast.

**Declared edge-image asymmetry** — whether the declared edge-image is closed
under negation.

The declared edge-image is **symmetric** when for every declared edge vector
v = Δ(x)[e], the vector −v also appears as a declared edge vector (i.e. every
declared edge (s,t) has a declared reverse edge (t,s)).
The declared edge-image is **asymmetric** when any declared edge vector v has
no matching −v among the declared edge vectors.

*Scope clarification:* This is edge-set symmetry —
symmetry over the explicit set of declared edge vectors. It is not span symmetry.
−v could belong to the span of declared edge vectors as a linear combination
without appearing as an explicit declared edge vector. This function tests the
former, not the latter. The implementation (`declared_edge_image_admissibility_check()` in
`operators_primary.rs`) computes exactly this: edge-set symmetry.

Symmetry of the declared edge-image means every declared directed relation has
a declared reverse — every contrast is matched by an explicitly declared reversal.
Asymmetry means at least one declared relation has no declared reverse.

**Distinction from operator antisymmetry (Σ):** Declared edge-image asymmetry
is a property of the explicit set of declared edge vectors. Operator antisymmetry,
measured by Σ, is a local property at each individual declared edge — whether the
adjacency-weighted circulation at that edge is non-trivially non-zero. These are
related but not identical:

- Declared edge-image asymmetry is a necessary condition for non-trivial Σ output
  — if the declared edge-image is globally symmetric, every directed difference
  is matched by its explicit reversal, and Σ's antisymmetric term will cancel at
  every edge.
- Declared edge-image asymmetry is not sufficient for non-trivial Σ output — a
  globally asymmetric declared edge-image can still produce zero local antisymmetry
  at every edge if the adjacency structure is symmetrically balanced (Failure
  Mode 3).

Declared edge-image admissibility is an admissibility check, not an analysis
function. A symmetric declared edge-image is an inadmissibility signal — it means
a reverse edge was declared without satisfying the distinctness axiom and independent
provenance requirement. Under the declared admissibility conditions, the edge-image
of any admissible declared structure is asymmetric. Σ then determines whether the
declared adjacency structure makes that asymmetry locally detectable.

### 3.1 Ring Topology — Derived Inadmissibility

The inadmissibility of ring topology is not an asserted prohibition. It is derived
from the admissibility condition alone.

**The derivation:**

A ring on nodes {0, 1, ..., n−1} requires a directed edge (n−1, 0) — the closing
edge that connects the last node back to the first.

Every declared edge must be traceable to an observable through M. The closing edge
(n−1, 0) declares that node n−1 continues into node 0 — that the state at node 0
is relationally dependent on the state at node n−1.

Following the declared chain: node 0 depends on node 1, node 1 on node 2, and so
on through node n−1, which depends on node 0 through the closing edge. This is a
circular dependency: x[0] is partly constituted by a chain that includes x[0]
as its own input.

A quantity that is partly constituted by itself is not traceable to an observable
through M. An observable is something that can be measured — it has a definite
value independent of the relational structure declared over it. A self-referential
quantity requires solving a fixed-point equation to have a definite value. That
fixed point is not an observable. It is a mathematical construction imposed on the
field.

The closing edge (n−1, 0) therefore cannot be traced to an observable through M.
It is inadmissible. And since the ring requires that edge, the ring is inadmissible.

**The shorter form:** A ring declares that a node is its own relational predecessor
through a chain. Within the declared admissibility conditions, no quantity has
itself as a relational predecessor — such a quantity would require solving a
fixed-point equation to have a definite value, which is not traceable to an
observable through M. The closing edge has no admissible provenance. The ring
fails the admissibility condition at the closing edge, and the closing edge is
the ring.

**Consequence for implementation:** `has_undirected_cycle()` in
`operators_primary.rs` detects any undirected cycle — including diamond DAGs which
are admissible when declared with provenance. A ring is detected by this function,
but the function does not distinguish rings from diamonds. The caller is responsible
for determining admissibility: a ring closing edge has no admissible provenance; a
diamond DAG has declared provenance at each edge. The derivation above is the test.

### 3.2 Relational Direction — Spatial and Persistence Unified

Every declared relation has exactly one admissible direction — the direction
traceable to an observable through M. This holds without distinction for spatial
relations, persistence relations, and hadronic transitions.

**Distinctness axiom (declared condition of this framework):**

If both (s,t) and (t,s) are independently observed and declared, they are treated
as distinct relations with distinct provenance — not as reverse observations of
the same relation. This is a declared axiom of the framework, not a derived
universal truth. Whether it holds for a specific observable is a question for
Origin at declaration time. If Origin determines that (s,t) and (t,s) are reverse
observations of the same underlying relation rather than independent relations,
then declaring both would be inadmissible regardless of independent provenance.

**Consequence for the spatial domain:**

Within the declared admissibility conditions, a reverse edge (t,s) is admissible
only when the direction t→s is independently traced to an observable through M
and treated as a distinct relation under the distinctness axiom. If the reverse
edge is derived by negating (s,t) rather than independently observed, it is
inadmissible. A symmetric declared edge-image therefore signals inadmissible
structure — not a legitimate declared state — unless each direction has
independent observable provenance and is treated as a distinct relation.

**Consequence for the persistence domain:**

A persistence relation connects E_prior to E_current across one declared relational
step. The direction is fixed:

    A_persistence[e] = E_current[e] − E_prior[e]

This is not a convention. E_prior is the state before the declared relational step.
E_current is the state after it. The ordering is determined by the observable.

The reverse — E_prior[e] − E_current[e] — would declare evolution running backward
through the declared relational step. That is not traceable to an observable through
M. A persistence relation declared in the reverse direction is fabricated within D.

**Unification:**

The spatial and persistence domains are subject to the same admissibility condition
on relational direction. The V5 finding that A_persistence[(s,t)] ≠ −A_persistence[(t,s)]
is not a special property of the persistence domain — it is the general property
of all declared relations becoming visible where directionality is explicit. The
boundary cancellation theorem does not transfer to the persistence domain because
the Lemma A(t,s) = −A(s,t) holds only for spatial node-valued loci where both
directions share the same node values. It does not hold where directions are
independently edge-valued.

V5 finding confirmed in `operators_v5.rs`:
`v5_boundary_cancellation_does_not_transfer_to_persistence`.

At the Primary Region, persistence is not yet active. This constraint is stated
here as a forward declaration for any M that introduces relational evolution
across declared event steps.

---

## 4. Operator Σ — Local Antisymmetry

Σ is the V5 R operator applied directly to Δ(x), without prior accumulation
through B. The operator formula is preserved exactly from V5. What changes is
the input: Σ receives Δ(x) where R receives B(A(x)). Σ is not a new operator —
it is R at the Primary Region limit, identified separately to mark the distinct
input condition.

**Architectural note:** B is absent from the primary composition — it is not
instantiated as an identity operator. The correct statement is: Σ is evaluated
directly on Δ(x) without an intermediate accumulation step. This is an
architectural choice that reflects the Primary Region condition: accumulation
along paths is not admissible where path structure is not yet confirmed.

With g = Δ(x) as input and ρ derived from Δ(x):

    Σ(g)[e] = g[e] + ρ[e] · ( Σ_{e' ∈ adj⁺(e)} g[e'] − Σ_{e' ∈ adj⁻(e)} g[e'] )

where:

**adj⁺(e)** — the set of declared edges whose source node equals the target node
of e. Edges that continue forward from e in the declared structure.

**adj⁻(e)** — the set of declared edges whose target node equals the source node
of e. Edges that arrive at the source of e in the declared structure.

**ρ[e]** — local contrast weight derived from Δ(x) at the source node s of e:

    ρ[e] = ρ_base · m[s] / (1 + m[s])

where m[s] = max{ |Δ(x)[e']| : e' incident to s } — the largest directed
difference at the source node. ρ is derived per edge from Δ(x); no aggregation
beyond the node. ρ ∈ [0, ρ_base) for all declared edges.

**Constraint:** couple only by the asymmetry present in the immediately declared
neighborhood of e. Do not assume adjacency that was not declared. Do not assume
the two directions are equal — symmetry holds only if it was declared. No
relation and no direction may be added that the declaration did not trace to an
observable through M.

**What Σ produces:**

A signed edge-valued field that measures, at each declared relation, whether the
immediately adjacent declared differences are asymmetrically distributed around e.

Σ(g)[e] = g[e] when adj⁺(e) and adj⁻(e) are both empty, or when the
antisymmetric term is exactly zero — Σ passes through Δ(x) unchanged.

Σ(g)[e] ≠ g[e] when the adjacent declared differences are asymmetrically
distributed around e — Σ amplifies or reduces the directed difference at e
according to its local relational neighborhood.

**Why Σ differs from V5 R in practice:** In the full V5 chain, R receives
B(A(x)) — values accumulated along paths through the declared interior. Σ
receives Δ(x) directly. The adjacency sums in Σ are over immediate neighbors
only — edges sharing a node with e — not over accumulated path values. This
makes Σ well-defined on any declared relational structure, including those with
no interior: single edges, isolated pairs, and minimally connected configurations
where B would have nothing to accumulate.

---

## 5. The Primary Kernel

    E_primary(x) = Σ(Δ(x))

This is the complete primary kernel. Two operators. One composition. No path
structure assumed. No interior declared beyond immediate adjacency.

**Expression** under the primary kernel is admissible when:

1. rank(Im Σ) > 0 — the dimension of the span of Σ output is non-trivial
2. The antisymmetric term of Σ is non-zero on at least one declared edge —
   what is present is asymmetrically organized in its local neighborhood

These two conditions are non-redundant. Condition 1 is satisfied by the
pass-through term g[e] alone, regardless of the antisymmetric term. Condition 2
is the additional requirement that distinguishes expression from mere presence.

---

## 6. Failure Modes

Three structural failure modes are derivable from the primary operators alone.

**Failure Mode 1 — Differentiation Collapse**

    rank(Im Δ) = 0

The observable is undifferentiated across all declared relations. Δ(x)[e] = 0
for all e. Σ(Δ(x)) = 0 identically. E_primary is undefined — the composition
has no non-trivial input.

Vocabulary term: *differentiation collapse.*

Distinguishing note: collapse is not the same as the observable being zero.
A uniform non-zero field collapses Δ equally: Δ(x)[e] = c − c = 0 for any
constant c. Collapse is a property of the relational structure of the observable
across declared edges, not of its magnitude.

**Failure Mode 2 — Relational Isolation**

    adj⁺(e) = ∅  and  adj⁻(e) = ∅  for all declared edges e

Every declared edge is isolated — no declared predecessor, no declared successor.
The antisymmetric term of Σ at every edge:

    ρ[e] · ( Σ_{adj⁺} g − Σ_{adj⁻} g ) = ρ[e] · (0 − 0) = 0

Σ reduces to the identity on Δ(x). No local antisymmetry is detectable because
no adjacency structure is declared through which asymmetry could manifest.

Vocabulary term: *circulation insufficiency.*

**Failure Mode 3 — Circulation Cancellation**

    ρ[e] · ( Σ_{adj⁺(e)} Δ(x) − Σ_{adj⁻(e)} Δ(x) ) = 0  for all e
    rank(Im Δ) > 0

The antisymmetric term of Σ is zero for all declared edges while the dimension
of the span of Δ output is non-trivially non-zero. The declared adjacency
structure is symmetrically balanced — every forward neighbor's directed difference
exactly cancels the backward neighbor's. Contrast is present but symmetrically
organized — no directional asymmetry is locally detectable.

*Condition alignment:* The primary condition is
rank(Im Δ) > 0 — contrast is present. rank(Im Σ) > 0 follows as a consequence:
when the antisymmetric term is zero, Σ reduces to the identity on Δ(x), so the
dimension of the span of Σ output equals the dimension of the span of Δ output.
Both are equivalent under FM3 conditions, but rank(Im Δ) > 0 is the stated
primary condition.

Vocabulary term: *coherence failure.*

Note: the declared edge-image of Δ may be globally asymmetric while Σ produces
zero antisymmetric output at every edge — this is exactly the distinction between
declared edge-image asymmetry and operator antisymmetry established in Section 3.
Failure Mode 3 is the case where global asymmetry exists in the declared
edge-image but the local adjacency structure cancels it edge by edge.

---

## 7. The Primary Region Characterized

The Primary Region is the space of declared relational configurations where:

    1 ≤ rank(Im Δ) < threshold(ρ_P → 1)

and E_primary = Σ(Δ(x)) is the admissible kernel.

The ρ_P ratio in primary operator terms:

    ρ_P = rank(Im Σ) / |{e ∈ E : adj⁺(e) ≠ ∅}|

where rank(Im Σ) is the dimension of the span of Σ outputs over ℝ, and the
denominator counts declared edges with at least one forward neighbor — the
propagation capacity of the immediate adjacency structure.

In the Primary Region: rank(Im Σ) is small (often 1, at marginal admissibility)
while propagation capacity may be larger. ρ_P ≪ 1.

The Atomic/Molecular transition — ρ_P ≈ 1 — occurs when the declared relational
structure has confirmed interior extent: enough consistently non-zero rank(Im Δ)
across multiple declared steps that accumulation along paths adds distinguishable
structure. At that point B becomes admissible, Σ receives B(Δ(x)) rather than
Δ(x) directly, and the full V5 kernel activates.

**Open condition:** The exact operator condition at
which ρ_P ≈ 1 constitutes the Atomic/Molecular transition is not yet derived as
a computable criterion. ρ_P is computable from declared operator outputs. The
threshold at which it signals a transition remains interpretive until formally
derived. No specific numerical value of ρ_P should be treated as a declared
transition threshold. See Open Condition 2 in Section 11.

---

## 8. Relationship to V5 and V6

**Σ is R.** The antisymmetric circulation formula is identical across both scales.
The difference is the input. At the Primary Region, R receives Δ(x) directly —
named Σ to mark this input condition. At the Atomic/Molecular scale and above, R
receives B(Δ(x)) — the accumulated path values. The operator formula does not
change. The input changes when B activates.

**Δ is A.** The directed difference formula is identical across all scales. Δ is
A renamed to mark its role as the irreducible primitive.

**B activates when persistence is confirmed.** B encodes the admissibility of
path accumulation — the claim that relational structure persists along declared
paths long enough to accumulate meaningfully. This claim is not admissible at the
Primary Region. It becomes admissible when ρ_P approaches 1 across confirmed
relational steps.

**B is absent, not identity.** At the Primary Region, B is absent from the
composition — it is not instantiated as an identity operator. The correct
statement: Σ is evaluated directly on Δ(x) without an intermediate accumulation
step. Saying "B = identity" is mathematically stronger than the implementation
establishes and is not the correct characterization.

**V6 primary kernel.** `operators_primary.rs` (V6) implements Δ and Σ as a
separately declared kernel with its own type system, failure mode detection,
expression condition, and test suite. The relationship to V5 is stated rather
than assumed. V5 (`operators_v5.rs`) is unchanged — V6 adds the primary kernel
alongside it.

**The operator hierarchy, stated without ambiguity:**

At the Primary Region:

    Input: x = M(o)
    Step 1: g = Δ(x)           [directed difference — A under renaming]
    Step 2: E_primary = Σ(g)   [antisymmetric circulation applied directly to Δ(x)]

At the Atomic/Molecular scale and above (B active):

    Input: x = M(o)
    Step 1: g₁ = Δ(x) = A(x)       [same operator; V5 name A]
    Step 2: g₂ = B(g₁)              [accumulation along declared paths]
    Step 3: E = R(g₂, ρ(g₁))        [same antisymmetric circulation formula; V5 name R]

Σ and R are the same operator. The name Σ marks the condition: R applied without
prior B accumulation. No new operator is introduced.

| Concept | Primary Region | Atomic/Molecular+ |
|---|---|---|
| Directed difference | Δ(x)[e] = x[s] − x[t] | A(x)[e] = x[s] − x[t] (identical) |
| Accumulation | Absent (not an identity operator) | B(g)[e] = g[e] + Σ_{succ} g[f] |
| Antisymmetric circulation | Σ(g)[e] = R formula applied to Δ(x) | R(g)[e] = same formula applied to B(A(x)) |
| Kernel | E_primary = Σ(Δ(x)) | E = R(B(A(x)), ρ(A(x))) |
| ρ source | Derived from Δ(x) at each node | Derived from A(x) at each node (identical) |

---

## 9. Sub-Observable Structure and Predictive Power

A quantity below the declared observable boundary is admissible in the formal
language if and only if it encodes Δ or Σ structure — directed difference
structure or local antisymmetry structure — that produces traceable consequences
at the declared observable level through M.

**Admissible below the observable boundary:** quantities that encode Δ or Σ
structure producing testable consequences at the observable level. Quarks are
admissible under this criterion — they encode the constraint and circulation
structure of hadronic observables, producing correct predictions at the hadronic
boundary without themselves being observable through any declared M.

**Inadmissible below the observable boundary:** quantities that neither correspond
to observables nor produce traceable consequences at the declared observable level.
These add structure to the formal language that has no operator-level consequence —
they are fabricated within D even if they compute.

The test for admissibility of sub-observable structure:

> Does removing this quantity from the formal description change the predictions
> at the declared observable level?

If not, the quantity is redundant — the relational structure it encodes is already
present in the observable field and does not require sub-observable entities to be
stated.

This criterion does not declare quarks unnecessary. It declares the condition
under which they would be unnecessary: if Δ and Σ applied to the hadronic
observable field recover the same predictions that QCD derives from quark
sub-structure. Whether that condition is met is exactly what the hadronic test
in Section 10 investigates.

### 9.1 Admissibility of Comparison Targets

The admissibility condition on M applies equally to comparison targets — the
quantities against which operator outputs are evaluated in any declared test.

**Standing constraint:** A quantity used to evaluate the output of Δ or Σ must
be traceable to an observable through M by the same standard as operator inputs.
A comparison target that is not itself traceable to a detector measurement or
declared observable through M is inadmissible as a test target, regardless of
its predictive success in other frameworks.

This constraint has three immediate consequences:

**Consequence 1 — Theoretical quantities are not admissible targets.**
Group-theoretic numbers (SU(3) multiplet dimensions), matrix elements derived
from quark assignments (CKM entries), and symmetry-group representations are
not observables through M. They may appear in the formal language as admissible
sub-observable encodings under the criterion of Section 9 — but they may not
serve as the comparison target of a declared test.

**Consequence 2 — Admissible comparison targets are measured quantities or
quantities computable from declared observables.**
Examples: directly measured quantities traceable to PDG detector data (particle
masses, decay branching ratios, CP violation parameter ε, magnetic moment ratios);
structural conditions on operator output computable from declared observables;
pass/fail conditions stated entirely in terms of operator outputs.

**Consequence 3 — The comparison target must be declared before the test runs.**
A comparison target stated after operator outputs are observed is a post-hoc
reinterpretation. The admissibility condition requires that the target be
traceable to M before execution.

*Constraint summary: what enters the operators must be traceable to an observable
through M. What the operator output is compared against must be traceable to an
observable through M by the same standard. Declaration is never a substitute for
traceability. This holds for inputs and targets equally.*

---

## 10. The Hadronic Test

The primary operators enable a specific, falsifiable test against quantum
chromodynamics at the hadronic observable boundary.

**The claim under test:** The hadronic mass spectrum, decay structure, and
conserved quantum number organization — which QCD derives from quark sub-structure
— are recoverable from E_primary = Σ(Δ(x)) applied to the declared hadronic
observable field, without invoking quark-level quantities.

This is a declared research program, not an established result. Three stages are
required. Only Stage 1 is complete as of this document.

**Declared M for the hadronic test:**

M maps from hadronic scattering and decay experiments to the declared observable
field. Each node v ∈ V is a hadronic configuration characterized by its directly
measurable conserved quantities:

- Mass-energy (MeV/c²) — from Particle Data Group tables
- Electric charge (units of e)
- Isospin and its third component
- Strangeness
- Baryon number

Each directed edge e = (s, t) ∈ E is a declared hadronic transition — a decay or
scattering process — traceable to detector events through M.

Every quantity in this M declaration is traceable to a detector event. No
quark-level quantity appears.

**Stage 1 — Δ declared over hadronic observables (complete):**

    Δ(x)[e] = x[s] − x[t]

where x[v] is the vector of conserved quantum numbers at hadronic configuration v.
Δ(x)[e] is the change in conserved quantum numbers across the declared transition e.

The admissible output of Stage 1: rank(Im Δ), the dimension of the span of the
declared edge vectors, and declared edge-image asymmetry. These are computable
directly from the declared observable field and constitute the declared basis for
Stage 2.

Correspondence between these operator-derived invariants and existing theoretical
classifications — including SU(3) multiplet structure — is a secondary
interpretive step, not a declared test target. Per Section 9.1, SU(3) multiplet
dimensions are group-theoretic quantities not directly traceable to detector
measurements through M.

**Stage 2 — Derivation that target observables are functions of Δ(x)
(not yet complete; declared as next step):**

Before any numerical prediction can be made, it must be demonstrated that the
target observable — for example, the magnetic moment ratio — is a function of the
Δ(x) vector over the declared hadronic edges.

For the proton-neutron magnetic moment ratio:

- x[p] = (938.3 MeV, +1, I=½, I₃=+½, S=0, B=1)
- x[n] = (939.6 MeV, 0, I=½, I₃=−½, S=0, B=1)
- Δ(x)[e] = x[p] − x[n] = (−1.3 MeV, +1, 0, +1, 0, 0)

over the declared isospin transition edge e = (p, n).

The question Stage 2 must answer: is the magnetic moment ratio derivable as a
function of the charge and isospin components of this Δ vector, without
constituent quark spin assignments? That derivation has not yet been performed.

**Stage 3 — Numerical comparison to measured values (contingent on Stage 2):**

The measured value of the proton-to-neutron magnetic moment ratio is −1.4599 —
directly traceable to PDG measurements through M. This is an admissible
comparison target under Section 9.1.

If Stage 2 produces a derivation that the ratio is a function of Δ(x), Stage 3
compares the primary operator prediction to this measured value. If Stage 2 fails,
that is a reportable falsification: sub-observable quark structure encodes
genuinely necessary relational information that Δ and Σ cannot recover from the
hadronic level alone.

Both outcomes carry equal weight as declared results.

**Falsification condition for the full hadronic test:**

If Δ and Σ over the declared hadronic observable field do not recover the directly
measured target quantities — if Stage 2 fails, or if Σ does not recover the
measured CP violation parameter ε from the declared kaon transition subgraph —
then quark sub-structure encodes genuinely necessary relational information that
the primary operators cannot recover from hadronic observables alone.

Note: correspondence with or failure to reproduce SU(3) multiplet dimensions is
not a falsification condition. SU(3) dimensions are not admissible test targets
under Section 9.1.

---

## 11. Open Conditions

1. **Persistence domain extension (elevated).** The V5 finding that boundary
   cancellation does not transfer to the persistence domain — because
   A_persistence[(s,t)] and A_persistence[(t,s)] are independently edge-valued —
   suggests that Σ over persistence loci may have fundamentally different algebraic
   behavior than Σ over spatial loci. If so, the primary operators as stated here
   apply to the spatial domain only, and their extension to Phase 2 (persistence)
   requires a separate derivation. This is the highest-priority open condition.

2. **Threshold characterization of ρ_P → 1 (open condition, not a computable
   criterion).** The exact operator condition at which the Primary Region
   transitions to the Atomic/Molecular Region — and B becomes admissible — is
   not yet stated as a computable criterion. ρ_P is computable; the transition
   threshold is not yet derived. No specific value of ρ_P should be treated as
   a declared threshold until this derivation is complete.

3. **Independence of Failure Modes 2 and 3.** Both produce zero antisymmetric
   Σ output by structurally different conditions. Whether they span distinct
   equivalence classes of relational failure is not yet demonstrated.

4. **Stage 2 of the hadronic test.** The derivation that the magnetic moment
   ratio is a function of Δ(x) over the declared hadronic observable field has
   not been performed. This is the next declared step of the research program.

5. **B activation condition.** The precise ρ_P threshold at which accumulation
   along declared paths adds distinguishable structure — and B therefore becomes
   admissible — is stated descriptively but not yet derived as a computable
   operator criterion.

6. **Span symmetry vs edge-set symmetry.** The current implementation computes
   declared edge-image asymmetry (edge-set symmetry). Whether span symmetry —
   whether −v ∈ span(declared edge vectors) — provides additional analytical
   power for the hadronic test or the per-event extension is not yet determined.
   If it does, a separate linear-algebraic test would be required.

---

## 12. Canonical Statements

> **Δ is the operator that asks whether anything distinguishable exists across
> the declared relations. Σ is the operator that asks whether what exists is
> symmetrically or asymmetrically organized in its local neighborhood. Σ is R
> applied directly to Δ(x) without prior accumulation — B is absent from the
> composition, not an identity operator. Together they constitute the irreducible
> primary kernel E_primary = Σ(Δ(x)), applicable at the smallest observable
> scale. B activates when persistence is confirmed, interposing accumulation
> between Δ and Σ, recovering the full V5 kernel. Bounded over D.
> No claim beyond D.**

> **A quantity below the declared observable boundary is admissible in the formal
> language if and only if it encodes Δ or Σ structure producing traceable
> consequences at the declared observable level through M. Whether quark
> sub-structure satisfies this criterion for hadronic observables is the subject
> of the declared hadronic test. No conclusion is reached before the test is run.**

> **The complete admissibility structure of the framework:**
> - **Inputs** require observable traceability through M.
> - **Operators** require declared behavior with stated constraints.
> - **Outputs** require declared projections with stated preserved and discarded
>   invariants.
> - **Evaluation targets** require observable traceability through M by the same
>   standard as inputs.
>
> *This structure is closed. No asymmetry remains between what enters the
> operators, what the operators do, what they produce, and what their output is
> compared against. All four are subject to the same admissibility condition.
> Bounded over D. No claim beyond D.*

> **Ring topology is inadmissible not by assertion but by derivation. The closing
> edge of a ring declares a node as its own relational predecessor through a chain.
> Within the declared admissibility conditions, no quantity has itself as a
> relational predecessor — such a quantity requires solving a fixed-point equation
> to have a definite value, which is not traceable to an observable through M.
> The closing edge is therefore inadmissible. The ring fails the admissibility
> condition at the closing edge, and the closing edge is the ring.**

---

## 13. Change Log

**V5 → V6**

| Section | Change | Reason |
|---|---|---|
| 3 | rank(Im Δ) wording corrected: "dimension of the span of declared edge vectors" replaces "count of linearly independent directed differences" throughout. Wording note added. | The two are equivalent only when every edge vector is linearly independent — the correct statement in all cases is dimension of the span. |
| 3 | Asymmetry definition narrowed from span symmetry to declared edge-image symmetry. Scope clarification added distinguishing edge-set symmetry from span symmetry. Function renamed `declared_edge_image_admissibility_check()` in implementation — recast as admissibility check rather than analysis function. Open Condition 6 added for span symmetry question. | The implementation checks whether −v appears as an explicit declared edge vector — not whether −v is in the span. The specification is narrowed to match what is actually computed. |
| 3.1 (new) | Ring topology inadmissibility derived from admissibility condition. Closes the logical gap between assertion and derivation. Canonical statement added in Section 12. | The closing edge of a ring declares a node as its own relational predecessor through a chain — not traceable to an observable through M. The prohibition follows from the admissibility condition; it need not be asserted separately. |
| 3.2 (new) | Relational evolution direction stated as a constraint: E_current − E_prior is the only admissible direction. Reverse direction is not traceable to an observable through M. Forward reference to V5 boundary cancellation finding. | The persistence direction is determined by the observable, not by convention. This was implicit in V5; it is now explicit. |
| 4 | "B = identity" replaced with accurate architectural statement: B is absent from the composition, not an identity operator. | B is absent — not instantiated as an identity operator. The distinction matters: one is an architectural choice, the other is a mathematical property. |
| 6 FM3 | Condition aligned to formal document: rank(Im Δ) > 0 stated as primary condition. rank(Im Σ) > 0 identified as consequence. | rank(Im Δ) > 0 is the stated primary condition in the formal document. rank(Im Σ) > 0 follows as a consequence when the antisymmetric term is zero. |
| 7 | ρ_P ≈ 1 transition explicitly marked as open condition. No specific numerical value to be treated as declared threshold. | The transition threshold is not yet derived as a computable criterion. ρ_P is computable; the threshold remains interpretive until formally derived. |
| 8 | "B = identity" corrected to "B is absent from the composition." Table updated. V6 primary kernel reference added. | Consistency with Section 4 correction. |
| 11 | Open Condition 2 expanded with explicit statement that ρ_P is computable but the transition threshold is not yet derived. Open Condition 6 added (span vs edge-set symmetry). | Matches corrections in Sections 7 and 3 respectively. |
| 3.1 | Ring prohibition canonical statement scoped to framework: "within the declared admissibility conditions" replaces "no observable has itself as a relational predecessor." | Universal claim scoped to framework — the derivation holds within the admissibility conditions; it is not a claim about nature. |
| 3.2 | Section renamed "Relational Direction — Spatial and Persistence Unified." Distinctness axiom stated explicitly as a declared condition of the framework. Spatial bidirectionality updated: admissible only with independent provenance and distinctness axiom satisfied. Spatial and persistence domains unified under one admissibility condition. | Verifier warning: implicit axiom made explicit; universal claims scoped to framework. |
| 3 | `declared_edge_image_asymmetry()` renamed `declared_edge_image_admissibility_check()`. Recast as admissibility check: symmetric edge-image is inadmissibility signal, not analysis output. | Consistent with Position B: symmetry signals declaration error under the admissibility conditions. |
| 12 | Ring canonical statement scoped: "within the declared admissibility conditions" replaces "no observable." | Scope discipline on universal claims. |
| 13 | Role references (Generator, Verifier, Origin as process role) removed from document. Operator description documents state what operators are and what constraints hold — not who performs what function in a workflow. | Role separation belongs in `role_separation_and_operator_application_v6.md`. |
| 12 | Canonical statement updated: "B is absent from the composition, not an identity operator." Ring prohibition derivation added as fourth canonical statement. | All above |

---

*Metatron Dynamics, Inc. — V6. Grounded in V6 primary kernel.*
*Bounded over D. No claim beyond D.*
