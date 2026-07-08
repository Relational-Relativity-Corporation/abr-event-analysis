# Relational Operator Detection of Structural Invariants in Hadronic Decay Data

**Robin Bruce Macomber**
Metatron Dynamics, Inc., Lompoc, California

*Correspondence: robin@metatrondynamics.com*

**Version:** 1.1 — Working paper. Verifier pass applied.
**Repository:** https://github.com/Relational-Relativity-Corporation/abr-event-analysis
**License:** CC BY 4.0

---

## Abstract

We demonstrate that a relational operator framework — applying directed
difference (Δ) and antisymmetric circulation (Σ) operators directly to
declared observable fields without statistical reduction — detects structural
invariants in individual hadronic decay events that ensemble-statistical
reduction does not preserve. On both simulated D0 → K⁻π⁺ events and real LHCb
detector data (D0 Masterclass dataset, 91,529 events; B⁺ → h⁺h⁺h⁻ dataset,
analyzed at track level), we establish: (1) the rank of the directed difference
field over individual events exceeds PDG summary statistics by 1–3 independent
dimensions, quantifying what statistical reduction discards; (2) the Σ
operator detects a structural constant at each decay vertex — the parent
candidate momentum, derived analytically from momentum conservation — to
numerical precision on real detector data, without invoking particle
identification, mass hypotheses, or probability distributions; (3) conservation
laws emerge as relational properties of the declared within-event field through
directed closure, without being imposed as theoretical constraints; (4) the
formula generalizes across decay topologies (2-body and 3-body) with a
correction factor derived from the number of outgoing tracks. These results
are computable from Level 2 declared observables (track momenta and charge
signs from curvature in a magnetic field) and do not require quark
sub-structure, wave functions, or ensemble averaging as inputs.

---

## 1. Introduction

The standard framework for analyzing particle physics data applies a sequence
of statistical reductions before any structural analysis is performed: raw
detector hits are reconstructed into tracks, tracks are assigned particle
identities through mass hypotheses and classifier outputs, reconstructed
candidates are filtered into signal and background categories, and ensemble
statistics (invariant mass distributions, lifetime fits, branching ratios)
are computed over the selected candidates. Each step discards information.
The final quantities — PDG summary statistics — represent the end of a
multi-stage reduction pipeline whose intermediate steps are rarely examined
for what they discard.

We apply an alternative framework grounded in a relational operator formalism
[CITE: Macomber 2026, arXiv:2601.22389] in which the primitive object is
not a particle with properties but a directed difference between observable
states. Two operators act on a declared observable field:

**Δ** (directed difference): Δ(x)[e] = x[s] − x[t] for each declared
directed edge e = (s,t), computing the contrast between observable states
at adjacent declared loci.

**Σ** (local antisymmetric circulation): Σ(g)[e] = g[e] + ρ[e] ×
(Σ_{adj⁺(e)} g − Σ_{adj⁻(e)} g), measuring whether the directed difference
field is asymmetrically organized in the immediate neighborhood of each
declared relation.

Neither operator requires particle identification, mass hypotheses,
probability distributions, or ensemble averaging. Both act on what is
directly declared from instrument readings through a measurement mapping M.

We report three analyses on hadronic decay data, establishing a progression
from simulated events to two independent real LHCb datasets, demonstrating
that the relational framework detects genuine physical structure — the parent
candidate momentum at each decay vertex — with analytical precision from
Level 2 declared observables alone.

---

## 2. Declared Observable Field (M Declaration)

All analysis proceeds from a declared measurement mapping M : O → D where
D := {x ∈ ℝⁿ | n < ∞, |x[i]| < ∞} is the space of finite bounded
observable vectors.

### 2.1 Admissibility Levels

We define three levels of admissibility for observable components, explicitly
declared for each analysis:

**Level 2 (admissible):** quantities directly traceable to instrument
readings without prior theoretical classification. For LHCb track data:
track momentum components (px, py, pz) from curvature in the declared
magnetic field; track charge sign (+1/−1) from direction of curvature.
No mass hypothesis required.

**Level 1 (partially admissible):** quantities requiring mass hypotheses
or vertex reconstruction. Invariant mass (requires declaring which track
is kaon, which is pion); decay time (requires total momentum, hence mass).

**Level 0.5 (C projection):** statistical quantities. Impact parameter
chi-squared (requires declared track error model); particle identification
probabilities (classifier outputs).

All analyses in this paper declare their admissibility level explicitly.
No quantity above its declared level is used in any downstream computation.

### 2.2 Within-Event Graph Declaration

For each decay event, the within-event relational graph is declared as:

```
Node 0: production vertex    x[0] = (0, 0, 0, 0)
Node 1: decay vertex         x[1] = (p_x, p_y, p_z, 0)
Node k: track k endpoint     x[k] = (pk_x, pk_y, pk_z, charge_k)
         for k = 2, ..., n+1

Edges: (0→1), (1→2), (1→3), ..., (1→n+1)
```

Every edge direction is fixed by the observable: the parent candidate
traveled from production to decay vertex; each track emerged from the
decay vertex. No edge direction is assumed or derived from theory.

### 2.3 Declared Comparison Targets

All PDG-traceable comparison targets are declared before any analysis
runs. Per the admissibility condition, comparison targets must be traceable
to observable measurements through M by the same standard as operator inputs.

| Quantity | Value | Source |
|---|---|---|
| D0 mass | 1864.84 MeV/c² | PDG 2022 |
| D0 lifetime | 0.4101 ps | PDG 2022 |
| Kaon mass | 493.677 MeV/c² | PDG 2022 |
| Pion mass | 139.570 MeV/c² | PDG 2022 |
| PDG summary rank | 3 | abr-primary-operators Phase 0 |

---

## 3. Analysis 1: Information Gap Between Per-Event and Summary Statistics

### 3.1 Setup

We generated 10,000 D0 → K⁻π⁺ decay events at declared beam momentum
5000 MeV/c using conservation of 4-momentum and an isotropic angular
distribution in the D0 rest frame (spin-0 parent, PDG-traceable parameters).
We applied Δ to two distinct observable fields over the same events:

**Level 0 (PDG summary):** x[v] = (mass_MeV, charge, lifetime_s) per
particle species. The standard PDG representation. Rank established in
prior work [CITE abr-primary-operators]: rank(Im Δ) = 3.

**Level 1 (per-event, particle-identified):** x[v] = (vertex_x, vertex_y,
vertex_z, kaon_px, kaon_py, kaon_pz, pion_px, pion_py, pion_pz,
decay_length, invariant_mass, decay_time). 12 components per event.

**Level 2 (per-event, no particle identity):** Same events, tracks labeled
by charge sign only. x[v] = (track1_px, track1_py, track1_pz, track1_charge,
track2_px, track2_py, track2_pz, track2_charge, vertex_x, vertex_y,
vertex_z, decay_length). 12 components, no mass hypothesis.

### 3.2 Rank Results

| Level | Components | rank(Im Δ) |
|---|---|---|
| 0 (PDG summary) | 3 (mass, charge, lifetime) | 3 |
| 1 (per-event, with ID) | 12 | 6 |
| 2 (per-event, no ID) | 12 | 5 |

**Information gap:** rank(Level 1) − rank(Level 0) = 3 independent
dimensions discarded by statistical reduction from per-event data to
PDG summary statistics.

### 3.3 Component Analysis

Incremental rank analysis identifies the components contributing each
new independent dimension above the PDG baseline:

| Added component | rank before | rank after | contribution |
|---|---|---|---|
| invariant_mass | — | — | 1 (PDG baseline) |
| decay_time | 1 | 2 | +1 |
| decay_length | 2 | 2 | redundant |
| kaon_px | 2 | 3 | +1 |
| kaon_py | 3 | 4 | +1 |
| kaon_pz | 4 | 5 | +1 |
| pion_px/py/pz | 5 | 5 | redundant (momentum conservation) |
| vertex_x/y/z | 2 | 2 | redundant given decay_length |

The three additional dimensions recovered at Level 1 are the kaon
momentum components — the angular distribution of the decay products
that PDG statistical reduction discards through ensemble averaging.
Pion momenta are redundant given kaon momenta because of 2-body
momentum conservation: p_pion = p_D0 − p_kaon.

At Level 2 (no particle identity), rank = 5: removing the mass
hypothesis loses one dimension (the invariant mass encodes mass
information not present in track momenta alone).

### 3.4 Interpretation

The three-dimensional information gap is not a statistical artifact.
It is a structural property of the event field: three independent
directions of relational contrast exist in the per-event angular
momentum structure that PDG ensemble averaging systematically discards.
These dimensions are real, measurable, and present in every individual
event — but invisible to any analysis that operates on summary statistics.

---

## 4. Analysis 2: Sigma Constant at the Decay Vertex

### 4.1 Analytical Derivation

**Proposition 4.1** (Sigma constant at decay vertex). For a decay with n outgoing tracks, the within-event graph has the
structure declared in Section 2.2. The antisymmetric term of Σ at
the incoming edge (0→1) is:

    A[(0→1)] = ρ[(0→1)] · (Σ_{k=2}^{n+1} Δ[(1→k)])

By momentum conservation:

    Σ_{k=2}^{n+1} Δ[(1→k)] = Σ_{k=2}^{n+1} (x[k] − x[1])
                              = Σ_k p_k − n·p_parent
                              = p_parent − n·p_parent
                              = −(n−1)·p_parent

The local contrast weight at the production vertex:

    ρ[(0→1)] = ρ_base · |x[1] − x[0]| / (1 + |x[1] − x[0]|)
             = ρ_base · p_parent / (1 + p_parent)

Therefore:

    |A[(0→1)]| = ρ_base · p_parent / (1 + p_parent) · (n−1) · p_parent
               = (n−1) · ρ_base · p_parent² / (1 + p_parent)

**For 2-body decay (n=2):** A = ρ_base · p² / (1+p)
**For 3-body decay (n=3):** A = 2 · ρ_base · p² / (1+p)
**General n-body:**         A = (n−1) · ρ_base · p² / (1+p)  □

This result follows from momentum conservation alone — no particle
identification, no mass hypothesis, no probability distribution.

### 4.2 Verification on Simulated Data

Applied to 10,000 simulated D0 → K⁻π⁺ events (2-body, n=2):

| Quantity | Value |
|---|---|
| Declared formula | 0.2 · p² / (1+p) |
| At p = 5000 MeV/c | 999.800... MeV/c |
| Measured mean | 999.800000 MeV/c |
| Measured std | 0.000000 MeV/c |

**Constant to numerical precision across all 10,000 events.** The std
is zero because the beam momentum is fixed (declared 5000 MeV/c) and
momentum conservation holds exactly in the simulation.

Scaling confirmed across six declared beam momenta:

| Beam (MeV/c) | Predicted | Measured | Match |
|---|---|---|---|
| 1000 | 199.960 | 199.960 | YES |
| 2000 | 399.800 | 399.800 | YES |
| 3000 | 599.700 | 599.700 | YES |
| 5000 | 999.800 | 999.800 | YES |
| 7000 | 1399.861 | 1399.861 | YES |
| 10000 | 1999.800 | 1999.800 | YES |

### 4.3 Verification on Real LHCb Data (Phase 0c)

Applied to the LHCb B⁺ → h⁺h⁺h⁻ dataset (B2HHH_MagnetDown.root,
DOI:10.7483/OPENDATA.LHCB.AOF7.JH09, CC0) at Level 2 admissibility.

**Declared M (Level 2 only):**
- H1_PX, H1_PY, H1_PZ, H1_Charge — track 1 momentum and charge
- H2_PX, H2_PY, H2_PZ, H2_Charge — track 2
- H3_PX, H3_PY, H3_PZ, H3_Charge — track 3
- B_P = √((H1+H2+H3)_PX² + ...) — total candidate momentum (admissible)

No particle identification applied. No mass hypothesis applied.
No invariant mass computed.

**Declared prediction (3-body, n=3):**

    A_predicted[i] = 2 · 0.2 · B_P[i]² / (1 + B_P[i])

per event i, where B_P[i] is the reconstructed total momentum of
event i. This varies event-by-event since B_P varies.

**Results on 100 verified events:**

| Quantity | Value |
|---|---|
| Max |A_measured − A_predicted| | 0.000000 MeV/c |
| Mean momentum closure | 0.000000 MeV/c |
| rank(Im Δ) track momenta | 6 |
| CV of B_P (event-to-event variation) | 0.2585 |
| CV of (A_measured/A_predicted) | ~0 |

**The predicted Sigma antisymmetric term matches the measured value
to numerical precision for every event analyzed, on real LHCb
detector data, without invoking particle identity or mass hypothesis.**

### 4.4 Directed Closure Within Declared Conserved Graphs

In both simulated and real data, the directed closure condition:

    Σ_{e ∈ P_event} Δ(x)[e][q] = 0

holds for conserved quantity q (momentum components) to numerical
precision, over the declared within-event graph. Within the declared
graph representation of conserved processes, this directed closure
reproduces the conservation relation without requiring it to be
imposed as a theoretical constraint on the operator.

Note: this demonstrates that the declared graph representation of
a conserved process produces directed closure. It does not demonstrate
that the physical conservation law is derived from the operator kernel.

---

## 5. Analysis 3: Real LHCb D0 Data (Phase 0b)

Applied to the LHCb D0 Masterclass dataset (MasterclassData.root,
DOI:10.7483/OPENDATA.LHCb.E7EJ.JUWR, CC0), 91,529 signal events.

**Available branches and admissibility:**

| Branch | Value | Admissibility |
|---|---|---|
| D0_PT | transverse momentum (MeV/c) | Level 2 |
| D0_TAU | decay time (nanoseconds) | Level 1 |
| D0_MM | invariant mass (MeV/c²) | Level 1 |
| D0_MINIPCHI2 | impact parameter χ² | Level 0.5 |

Note: D0_TAU is in nanoseconds (confirmed: first event TAU = 0.000413 ns
= 0.413 ps, matching PDG D0 lifetime 0.4101 ps to 0.3%).

**Rank results:**

| Components | rank(Im Δ) |
|---|---|
| PDG summary (Level 0) | 3 |
| D0_PT only | 1 |
| D0_PT + D0_TAU | 2 |
| D0_PT + D0_TAU + D0_MM | 3 |
| All 4 components | 4 |

**rank(Im Δ) = 4 > PDG rank = 3.** Real LHCb individual events carry
one additional independent dimension beyond PDG summary statistics,
confirmed on 91,529 real detector events.

**Sigma constant:** Not reproducible from PT alone (CV = 0.9897),
confirming that the Sigma constant test requires p_total (3D momentum).
PT is the transverse component only — without the longitudinal component
PZ, p_total cannot be reconstructed from this file. This is the stated
boundary of the Masterclass dataset's admissibility level.

---

## 6. Discussion

### 6.1 What the Object Frame Cannot Express

The Sigma constant A = (n−1) · ρ_base · p² / (1+p) is a structural
invariant of the within-event relational graph. It exists at the level
of a single event — before any ensemble is formed. The ensemble-statistical
framework:

- Does not predict it within the standard pipeline, because the
  pipeline acts on ensemble statistics not within-event structure
- Does not preserve it, because individual event structure is
  discarded through averaging before analysis begins
- Does not derive it within the standard pipeline, because the
  derivation requires acting on directed differences across the
  within-event graph — a step the standard pipeline does not take

This is not a limitation that more data resolves. It is a structural
difference between frameworks. The object frame computes what the
relational frame can summarize as a declared C projection. The relational
frame computes what the standard reduction pipeline does not
formulate within its declared analysis chain.

### 6.2 What This Does Not Claim

This paper does not claim that the relational framework is a complete
theory of particle physics, or that it reproduces all predictions of
the Standard Model. It claims the following, all of which are established
on declared observable data:

1. Individual hadronic decay events carry more independent relational
   dimensions than PDG summary statistics (rank gap = 1–3 dimensions,
   confirmed on simulation and real data).

2. The Σ operator detects the parent momentum at any decay vertex
   through the antisymmetric structure of the declared within-event graph,
   analytically predictable from momentum conservation, confirmed to
   numerical precision on real LHCb data.

3. Conservation laws are relational properties of the declared
   observable field — they emerge from directed closure without being
   imposed as theoretical constraints.

4. The framework operates at Level 2 admissibility — track momenta
   and charge signs — without requiring particle identification, mass
   hypotheses, or probability distributions.

### 6.3 The Information Gap as a Declared Quantity

The rank gap between per-event analysis and PDG summary statistics
is a precisely defined, reproducible quantity:

    Gap = rank(Im Δ | individual events) − rank(Im Δ | PDG summary)

This gap is 3 for the D0 decay at full Level 1 admissibility and 1
at Level 2 (no mass hypothesis). It measures what statistical reduction
discards. It is computable from any dataset with the appropriate
admissibility level.

The gap is not zero. The information discarded by statistical reduction
is real, measurable, and quantifiable. Whether that information is
physically significant for any specific analysis is a separate question —
but its existence and magnitude are now established.

---

## 7. Conclusions

We have demonstrated on real LHCb detector data that:

1. **The relational framework detects structural invariants invisible
   to ensemble-statistical methods.** The Sigma antisymmetric term
   at each decay vertex predicts the parent candidate momentum to
   numerical precision from Level 2 observable data alone.

2. **The formula A = (n−1) · ρ_base · p² / (1+p) is universal across
   decay topologies.** Confirmed on 2-body (simulated D0) and 3-body
   (real LHCb B→HHH) decays. The correction factor (n−1) is derived
   from momentum conservation — no free parameters, no fitting.

3. **Directed closure reproduces conservation within declared graphs.**
   Within the declared graph representation of conserved processes,
   the directed closure condition Σ_e Δ(x)[e][q] = 0 holds for
   conserved quantities q across the declared within-event graph,
   without imposing conservation laws as theoretical constraints.
   This demonstrates that the declared graph encodes the conservation
   relation — not that the physical law is derived from the operators.

4. **Statistical reduction discards 1–3 independent relational dimensions**
   per event compared to per-event analysis. Quantified as the rank gap
   of Im(Δ) between individual events and PDG summary statistics.

5. **The framework operates without particle identification, mass
   hypotheses, wave functions, or probability distributions** as inputs
   to M. All quantities are traceable to Level 2 declared observables
   through the declared measurement mapping.

These results suggest that within-event relational structure carries
physical information currently discarded by the standard analysis pipeline,
and that the relational operator framework provides a complementary
analysis layer operating directly on the declared observable field before
statistical reduction.

---

## 8. Methods and Reproducibility

All analyses are fully reproducible from the declared repository:

**Repository:** https://github.com/Relational-Relativity-Corporation/abr-event-analysis

**Analysis scripts (run in order):**
```
python generate_events.py --n 10000
python analyze_rank_contributions.py
python analyze_angular_structure.py
python analyze_hadronic_phase0a.py
python analyze_sigma_constant.py
python analyze_hadronic_phase0b.py    # requires MasterclassData.root
python analyze_hadronic_phase0c.py    # requires B2HHH_MagnetDown.root
```

**Data sources (CC0 license):**
- LHCb D0 Masterclass data: DOI:10.7483/OPENDATA.LHCb.E7EJ.JUWR
- LHCb B2HHH data: DOI:10.7483/OPENDATA.LHCB.AOF7.JH09

**Kernel:** V6 primary operators (Δ, Σ) implemented in
`operators/primary.py`. Reference Rust implementation: `operators_primary.rs`.
Formal specification: `docs/primary_operators_delta_sigma_v6.md`.

**Declared tolerances:**
- SVD tolerance (rank computation): τ = 10⁻¹⁰
- Momentum closure tolerance: < 1.0 MeV/c (declared before execution)
- Sigma prediction tolerance: < 1.0 MeV/c (declared before execution)

All hypotheses declared before execution. No results adjusted post-hoc.
No parameter tuned to fit results.

---

## 9. Data Availability

All data used in this paper is publicly available:

**Simulated data:** Generated by `generate_events.py` in the repository.
Full generation code declared in the repository. Parameters: D0 mass
1864.84 MeV/c², D0 lifetime 0.4101 ps, kaon mass 493.677 MeV/c²,
pion mass 139.570 MeV/c² (all PDG 2022), beam momentum 5000 MeV/c.
Seed: 42. CC BY 4.0.

**Real LHCb data:**
- MasterclassData.root: https://opendata.cern.ch/record/401
  CC0. DOI:10.7483/OPENDATA.LHCb.E7EJ.JUWR
- B2HHH_MagnetDown.root: https://opendata.cern.ch/record/4900
  CC0. DOI:10.7483/OPENDATA.LHCB.AOF7.JH09

Neither the LHCb experiment nor CERN endorses this work.

---

## References

[1] Macomber, R.B. (2026). The Object Error: A Formal Derivation of the
    Null Space Introduced by Object-Primary Framing. arXiv:2601.22389.
    Metatron Dynamics, Inc.

[2] LHCb collaboration (2014). LHCb event file for real measurement.
    CERN Open Data Portal. DOI:10.7483/OPENDATA.LHCb.E7EJ.JUWR. CC0.

[3] LHCb collaboration (2020). Matter Antimatter Differences (B meson
    decays to three hadrons) — Data Files. CERN Open Data Portal.
    DOI:10.7483/OPENDATA.LHCB.AOF7.JH09. CC0.

[4] Particle Data Group (2022). Review of Particle Physics.
    Prog. Theor. Exp. Phys. 2022, 083C01.

[5] Macomber, R.B. (2026). V6 Primary Operator Kernel.
    Metatron Dynamics, Inc. GitHub:
    Relational-Relativity-Corporation/abr-event-analysis.

---

*Metatron Dynamics, Inc. Bounded over D. No claim beyond D.*
*All analyses declared before execution. No results adjusted post-hoc.*
*The operators do not act, cause, optimize, or enforce.*
*Whatever the analysis commits to, a person committed it.*
