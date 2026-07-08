# The Operators — Plain Language

**Metatron Dynamics, Inc.** V6. Bounded over D. No claim beyond D.

---

A statement is admissible when every quantity in it can be traced to an observable through a declared measurement mapping M.
When it cannot, the statement is inadmissible — regardless of whether it computes, sounds reasonable, or produces a result.

Declarability is necessary but not sufficient. A quantity that is declared but not traceable to an observable through M is not inadmissible — it is simply not about anything in D. It exists only within the formal language.

Declaration is never a substitute for traceability to an observable through M. Declaration is the act of stating traceability to an observable through M and taking responsibility for it. A relation, quantity, or structure that cannot be traced to an observable through M is fabricated within D — regardless of who declared it, regardless of whether it computes, and regardless of whether it sounds reasonable. The admissibility condition is not relaxed by the act of declaration. It is what declaration is required to satisfy.

The operators below make those constraints explicit.
They do not add structure. They require that structure already present be declared before it is used.

The reader remains Origin.
This document states constraints.
Verification remains available to you directly.

---

## Before anything is measured

A domain must be declared before any operator acts.
A measurement mapping M must be declared before any operator acts.
The operators act on M(o) only — on what M produces from the observable, and nothing else.

*Constraint: there is nothing to evaluate before something is declared.*

---

## Every relation has one direction

Every declared relation has exactly one admissible direction — the direction traceable to an observable through M.

The direction of a relation is not a choice. It is determined by what was observed. If you observed something moving from s to t, the direction is s → t. Declaring the reverse — t → s — requires a separate observable that supports that direction independently. If you do not have one, the reverse direction is not traceable to an observable through M. It is inadmissible.

This applies to every kind of declared relation: spatial relations between measurement loci, transitions between particle configurations, evolution from one state to the next. In every case, the direction is determined by the observable, not by the declaring party's preference.

**What this means for pairs of relations:**

If both (s,t) and (t,s) are declared, each direction must be independently traceable to an observable through M. When both directions have independent provenance, they are treated as distinct relations — not as the same relation observed from two sides. If one direction is simply the mathematical negation of the other with no independent observable supporting it, it is inadmissible.

Within the declared admissibility conditions of this framework, a declared structure where every edge has a matching reverse edge signals a declaration error — not a legitimate relational configuration. The operators will detect this and report it.

**Why rings are inadmissible:**

A ring declares that a node is its own relational predecessor through a chain: 0 → 1 → 2 → ... → n−1 → 0. Following this chain, the state at node 0 depends on a chain that includes node 0 as its own input. Within the declared admissibility conditions, no quantity has itself as a relational predecessor — such a quantity would need to solve an equation involving itself to have a definite value, and that solution is not an observable. The closing edge of a ring therefore cannot be traced to an observable through M. It is inadmissible. The ring fails at the closing edge, and the closing edge is the ring.

*Constraint: every declared relation must have a single admissible direction independently traceable to an observable through M. The reverse direction requires independent provenance. Without it, it is inadmissible.*

---

## Δ and Σ — the primary operators

At the smallest scale of declared relational structure — where there is no confirmed path interior and no accumulated relational history — two questions are irreducible:

1. Does anything distinguishable exist across the declared relations?
2. Is what exists symmetrically or asymmetrically organized in its immediate neighborhood?

**Δ — directed difference** asks the first question.

Δ takes the directed difference of the observable field across each declared relation: Δ(x)[e] = x[s] − x[t] for each declared edge e = (s, t). It produces one value per declared relation — the contrast across that relation in the direction it was declared.

*Constraint: take the directed difference and nothing else. No relation and no direction may be added that the declaration did not trace to an observable through M.*

**Σ — local antisymmetry** asks the second question.

Σ takes the Δ output and asks: at each declared relation, is the immediately adjacent contrast distributed asymmetrically around it? It adds the differences from edges that continue forward and subtracts the differences from edges that arrived from behind, scaled by local contrast. It acts on immediate neighbors only — no path accumulation, no assumed interior.

*Constraint: couple only by the asymmetry present in the immediately declared neighborhood. Do not assume adjacency that was not declared.*

**Together** they constitute the primary kernel: E_primary = Σ(Δ(x)). Two operators. One composition. No path structure assumed.

Under the directional admissibility condition, every admissible declared structure is asymmetric. This means Σ will always detect non-zero local antisymmetry wherever adjacency is declared and the observable is not uniform. Symmetry — where every declared relation has a matching reverse — is a signal that the declaration contains inadmissible structure, not a legitimate state the operators can act on.

---

## A — directed difference (ABR kernel)

A measures directed difference across declared relations.

In the spatial domain, a relation is a directed edge (s, t) whose existence is traceable to an observable through M. A takes the difference x[s] − x[t] and nothing else.

In the relational-evolutionary domain (V6), a relation connects two complete kernel output states across one declared relational step. A takes the difference E_current[e] − E_prior[e] — the same directed-difference formula, applied to edge-valued loci rather than node-valued loci. The direction is fixed: E_current − E_prior. E_prior is the state before the declared relational step; E_current is the state after it. The reverse direction is not traceable to an observable through M.

*Constraint: no relation and no direction may be added that the declaration did not trace to an observable through M. This holds for both spatial and persistence loci.*

---

## B — accumulation along declared continuation

B accumulates along declared continuation and nowhere else.

At each edge, B adds the values of edges that continue forward from it.
A terminal edge accumulates nothing.
No boundary is closed to supply continuation that was not declared.

B is absent from the primary kernel — it is not an identity operator that happens to do nothing. It is simply not invoked. B activates when persistence is confirmed: when enough consistently non-zero relational contrast has been established across declared steps that accumulation along paths adds distinguishable structure.

*Constraint: accumulation follows declared structure. It does not supply structure.*

---

## R — coupling through observed asymmetry (ABR kernel)

R couples declared relations through observed asymmetry.

At each edge, R adds the difference between what continues forward and what arrives from behind, scaled by local contrast. Where declared relation families couple, the asymmetry of each contributes to the other according to the declared coupling.

Under the directional admissibility condition, every admissible declared structure is asymmetric — each declared relation has a single admissible direction, and the reverse requires independent provenance. Symmetry is therefore not a legitimate declared state within the admissibility conditions of this framework. Do not assume the two directions are equal: without independent provenance for each direction, only one direction is admissible.

*Constraint: couple by the asymmetry present. Every declared relation has one admissible direction. The reverse requires independent provenance.*

---

## ρ — local contrast

ρ scales coupling according to local contrast.

At each node, ρ is derived from the largest gradient present at that node in the operator's output.
ρ does not aggregate beyond the node.

*Constraint: coupling strength is derived locally. It is not assigned globally.*

---

## C — declared projection

C reports a declared projection and states what was discarded.

Any reduction of the field — to one value per node, to a bound, to a variance — is a projection.
A projection is admissible when what it preserves and what it discards are stated.

*Constraint: no silent reduction.*

---

## Before Δ / Before A

Differences may not be altered before measurement unless the alteration is declared within M.

Admissible before Δ or A: uniform shift, declared unit scale.
Everything else requires declaration with preserved and discarded invariants stated.

*Constraint: do not change what you are about to measure without saying so.*

---

## What this produces

**At the Primary Region** — where persistence is not yet confirmed and path accumulation has not been established — the primary kernel applies:

Δ, then Σ: the operators ask whether anything distinguishable exists, and whether what exists is asymmetrically organized in its immediate neighborhood. The result reflects the declared relational structure of the observable at its most minimal. No interior is assumed. No history is carried.

**When persistence is confirmed** — the ABR kernel applies in two phases:

Phase 1 (spatial): A → B → R over declared spatial relations produces the spatial relational field.

Phase 2 (persistence): A → B → R over declared persistence relations — each connecting the prior complete kernel output to the current one across one relational step — produces the relational-evolution field over one declared relational step.

Both phases use the same operator formulas. What changes between phases is what the operators act over, not how they act. The direction of every relation — spatial or persistence — is determined by the observable and fixed by the declaration.

The result is not an interpretation.
It is a function of the declaration and the observable.
Change the declaration and the result changes.
The same declaration on the same observable produces the same result.

---

## What remains with the reader

The operators do not interpret findings.
The operators do not assign meaning to departures.
The operators do not determine what constitutes a significant result.

Those determinations remain with whoever declared the domain.

*Constraint is not authority.*
*The reader remains Origin.*

---

*Bounded over D. No claim beyond D.*
*Metatron Dynamics, Inc. V6.*

---

**V5 → V6 changes:** New section: "Every relation has one direction" — states the directional admissibility condition and distinctness axiom in plain language; derives ring inadmissibility from first principles. New section: "Δ and Σ — the primary operators" — plain language description of the primary kernel, including the consequence of Position B for Σ. A section: persistence direction stated as fixed (E_current − E_prior); reverse direction named as inadmissible. B section: absence from primary kernel distinguished from identity operator. R section: symmetry recharacterized — within the admissibility conditions, every admissible declared structure is asymmetric; symmetry signals inadmissible structure. Before Δ / Before A: heading updated to cover both operators. What this produces: primary kernel output added alongside ABR kernel output, with condition for each. No universal claims about observables — all claims scoped to "within the declared admissibility conditions of this framework."
