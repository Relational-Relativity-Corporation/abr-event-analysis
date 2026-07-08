# abr-event-analysis

**Metatron Dynamics, Inc.** V6. Bounded over D. No claim beyond D.

Relational operator analysis of hadronic decay data.
Applies the V6 primary operator kernel (Δ, Σ) directly to declared
observable fields before statistical reduction.

---

## Findings

Three analyses on hadronic decay data establish the following results,
all declared before execution and confirmed on real LHCb detector data:

**1. Information gap: per-event rank exceeds PDG summary statistics**

The rank of the directed difference field Im(Δ) over individual decay
events exceeds the rank over PDG summary statistics by 1–3 independent
dimensions. Statistical reduction discards real relational structure.

```
Level 0  (PDG summary, 3 components):              rank = 3
Level 1  (per-event, particle-identified, 12 comp): rank = 6
Level 2  (per-event, no particle ID, 12 comp):      rank = 5
Level 1.5 (real LHCb D0 data, 4 components):        rank = 4
```

**2. Sigma constant at decay vertices (confirmed on real data)**

The Σ operator detects the parent candidate momentum at any decay vertex
through the antisymmetric structure of the declared within-event graph.
The formula is derived analytically from momentum conservation:

```
A = (n-1) · rho_base · p_parent² / (1 + p_parent)
```

where n is the number of outgoing tracks and p_parent is the total
candidate momentum. Confirmed to numerical precision:

- Simulated D0 → K⁻π⁺ (2-body): std = 0.000000 MeV/c across 10,000 events
- Real LHCb B → h⁺h⁺h⁻ (3-body): max error = 0.000000 MeV/c on real detector data

No particle identification. No mass hypothesis. No probability distributions.

**3. Directed closure reproduces conservation**

Within the declared graph representation of conserved processes, the
directed closure condition Σ_e Δ(x)[e][q] = 0 holds for conserved
quantities q to numerical precision, without imposing conservation laws
as theoretical constraints.

**4. General formula across decay topologies**

The Sigma formula generalizes across decay topologies with a declared
correction factor (n-1) derived from momentum conservation:

```
2-body (D0 → Kπ):   A = 1 · rho_base · p² / (1+p)
3-body (B → HHH):   A = 2 · rho_base · p² / (1+p)
n-body:              A = (n-1) · rho_base · p² / (1+p)
```

---

## Data

All data used in this analysis is publicly available under CC0 license.

| File | Source | DOI | Size |
|---|---|---|---|
| `data/MasterclassData.root` | CERN Open Data record 401 | 10.7483/OPENDATA.LHCb.E7EJ.JUWR | ~1.2 MB |
| `data/B2HHH_MagnetDown.root` | CERN Open Data record 4900 | 10.7483/OPENDATA.LHCB.AOF7.JH09 | ~660 MB |

Download:
```powershell
# D0 Masterclass data (~1.2 MB)
Invoke-WebRequest -Uri "https://opendata.cern.ch/eos/opendata/lhcb/LHCb-Masterclass/D0lifetime/data_mypy.root" -OutFile "data\MasterclassData.root"

# B->HHH data (~660 MB)
Invoke-WebRequest -Uri "https://opendata.cern.ch/eos/opendata/lhcb/AntimatterMatters2017/data/B2HHH_MagnetDown.root" -OutFile "data\B2HHH_MagnetDown.root"
```

Neither the LHCb experiment nor CERN endorses this work.

---

## Reproduce

```powershell
# Install dependencies
pip install -r requirements.txt

# Run full analysis chain
.\reproduce.ps1

# Simulation only (no data files required)
.\reproduce.ps1 -SkipRealData
```

Results are written to `results/`.

---

## Analysis Chain

| Script | Description | Data required |
|---|---|---|
| `generate_events.py` | Generate 10,000 simulated D0 events | None |
| `analyze_rank_contributions.py` | Information gap: rank per component | Simulated |
| `analyze_angular_structure.py` | Isotropy confirmation vs QM prediction | Simulated |
| `analyze_hadronic_phase0a.py` | Level 2 analysis, no particle identity | Simulated |
| `analyze_sigma_constant.py` | Sigma constant derivation and scaling | Simulated |
| `analyze_hadronic_phase0b.py` | Real LHCb D0 data, 91,529 events | MasterclassData.root |
| `analyze_hadronic_phase0c.py` | Real LHCb B→HHH, Level 2, Sigma constant | B2HHH_MagnetDown.root |

---

## Repository Structure

```
abr-event-analysis/
  declaration/          M declaration package
    observable.py       Observable vector components (Level 2 admissible)
    relations.py        Directed relation declaration with provenance
    targets.py          PDG comparison targets (declared before execution)
    admissibility.py    AdmissibilityError

  operators/            V6 primary kernel (Python)
    primary.py          Delta, Sigma, rho, antisymmetric_term

  analysis/             Declared C projections of operator output
    rank.py             rank(Im Delta), rank(Im Sigma)
    admissibility.py    Declared edge-image admissibility check
    expression.py       Expression condition, failure mode detection
    rho_p.py            rho_P ratio

  data/                 Data directory (large files not in git)
    events/             Generated event CSV files
    MasterclassData.root   Real LHCb D0 data (download separately)
    B2HHH_MagnetDown.root  Real LHCb B->HHH data (download separately)

  results/              Generated analysis reports (not in git)

  docs/                 V6 kernel documents
    operators_notation_and_constraint_v6.md
    primary_operators_delta_sigma_v6.md
    primary_region_formal_interior_v6.md
    role_separation_and_operator_application_v6.md
    abr_operators_plain_v6.md

  kernel/               Reference Rust implementation
    operators_primary.rs   V6 primary kernel, 36/36 tests passing

  paper_relational_hadronic_v1.md   Research paper (working draft)
  reproduce.ps1                      Full analysis chain (one command)
  generate_events.py
  analyze_*.py
```

---

## Kernel

**V6 primary kernel:** E_primary(x) = Σ(Δ(x))

The operators act on declared observable fields. No path accumulation.
No persistence. No statistical reduction before the operators act.

**Formal specification:** `docs/primary_operators_delta_sigma_v6.md`

**Reference implementation:** `kernel/operators_primary.rs`
(Rust, 36/36 tests passing, includes all failure mode detection and
directional admissibility checks)

**Python implementation:** `operators/primary.py`

---

## Admissibility

Every quantity must satisfy:

1. Produced by a physical instrument directly
2. Finite and bounded (in D)
3. Does not require prior theoretical classification to obtain

**What is in M (Level 2):**
- Track momentum components (px, py, pz) — from curvature in declared B field
- Track charge sign (+1 or -1) — from direction of curvature

**What is not in M:**
- Particle identity (kaon, pion, D0) — theoretical classification
- Invariant mass — requires mass hypothesis
- Probability distributions — ensemble C projections
- Wave functions — not observable

Declaration is never a substitute for traceability to an observable
through M. No claim beyond D.

---

## Declared Constraints

Every admissibility condition, relational direction constraint,
and operator constraint is formally stated in the V6 kernel documents
under `docs/`. Key constraints:

- Every declared relation has exactly one admissible direction
- Symmetric declared edge-images signal inadmissible structure
- Ring topology is inadmissible (derived, not asserted)
- Correlation is admissible only as a declared C projection with
  stated preserved and discarded invariants
- No observable exists outside D

---

## Citation

```
Macomber, R.B. (2026). Relational Operator Detection of Structural
Invariants in Hadronic Decay Data. Metatron Dynamics, Inc.
Working paper. GitHub: Relational-Relativity-Corporation/abr-event-analysis
```

Related work:
```
Macomber, R.B. (2026). The Object Error: A Formal Derivation of the
Null Space Introduced by Object-Primary Framing. arXiv:2601.22389.
```

---

## License

Code: MIT
Documents and paper: CC BY 4.0
Data: CC0 (CERN Open Data)

---

*Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.*
*All analyses declared before execution. No results adjusted post-hoc.*
