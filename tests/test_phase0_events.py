"""
tests/test_phase0_events.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

Tests for Phase 0 per-event analysis.
Written before execution — per the admissibility condition on declared tests.
"""

import numpy as np
import pytest
from declaration.observable import ObservableVector, OBSERVABLE_COMPONENTS, N_COMPONENTS
from declaration.relations import declare_relations, RelationProvenance
from declaration.admissibility import AdmissibilityError
from operators import operator_delta, operator_e_primary
from analysis import (
    im_delta_rank, declared_edge_image_admissibility_check,
    detect_failure_mode, FailureMode,
)


def make_event_field(n=5):
    """Synthetic event field for testing — strictly directed."""
    x = {}
    for i in range(n):
        vec = np.array([
            float(i),          # vertex_x
            float(i) * 0.1,    # vertex_y
            float(i) * 10.0,   # vertex_z
            100.0 + i,         # kaon_px
            50.0 - i,          # kaon_py
            200.0 + i * 2,     # kaon_pz
            -100.0 - i,        # pion_px
            -50.0 + i,         # pion_py
            -200.0 - i * 2,    # pion_pz
            float(i) * 5.0,    # decay_length
            1864.84,           # invariant_mass
            float(i) * 0.05,   # decay_time
        ])
        x[i] = vec
    return x


def make_edges(n=5):
    event_ids = list(range(n))
    return declare_relations(
        event_ids=event_ids,
        relation_type=RelationProvenance.TEMPORAL,
        direction_basis="declared production order"
    )


class TestObservableVector:

    def test_correct_dimension(self):
        x = make_event_field()
        v = ObservableVector(0, x[0])
        assert len(v.data) == N_COMPONENTS

    def test_rejects_wrong_dimension(self):
        with pytest.raises(AdmissibilityError):
            ObservableVector(0, np.zeros(5))

    def test_rejects_non_finite(self):
        vec = np.ones(N_COMPONENTS)
        vec[0] = np.inf
        with pytest.raises(AdmissibilityError):
            ObservableVector(0, vec)


class TestDeclareRelations:

    def test_temporal_strictly_directed(self):
        edges = make_edges(5)
        assert len(edges) == 4
        # All edges go forward
        for s, t in edges:
            assert t == s + 1

    def test_no_reverse_edges(self):
        edges = make_edges(5)
        edge_set = set(edges)
        for s, t in edges:
            assert (t, s) not in edge_set, \
                "Reverse edge found — inadmissible under directional admissibility"

    def test_rejects_empty_direction_basis(self):
        with pytest.raises(AdmissibilityError):
            declare_relations([0, 1, 2], RelationProvenance.TEMPORAL, "")

    def test_rejects_single_event(self):
        with pytest.raises(AdmissibilityError):
            declare_relations([0], RelationProvenance.TEMPORAL, "order")


class TestDelta:

    def test_directed_difference_correct(self):
        x = make_event_field()
        edges = make_edges()
        field = operator_delta(x, edges)
        for e, (s, t) in enumerate(edges):
            expected = x[s] - x[t]
            np.testing.assert_allclose(field[e], expected)

    def test_shape(self):
        x = make_event_field(10)
        edges = make_edges(10)
        field = operator_delta(x, edges)
        assert field.shape == (9, N_COMPONENTS)

    def test_all_finite(self):
        x = make_event_field()
        edges = make_edges()
        field = operator_delta(x, edges)
        assert np.all(np.isfinite(field))


class TestAdmissibility:

    def test_directed_chain_is_admissible(self):
        x = make_event_field()
        edges = make_edges()
        field = operator_delta(x, edges)
        admissible, witness = declared_edge_image_admissibility_check(field)
        assert admissible, "Directed chain must be admissible (asymmetric)"
        assert witness is None

    def test_per_event_rank_exceeds_summary(self):
        """Per-event rank should exceed PDG summary rank of 3."""
        x = make_event_field(20)
        edges = make_edges(20)
        field = operator_delta(x, edges)
        rank, _ = im_delta_rank(field)
        PDG_RANK = 3
        assert rank >= PDG_RANK, \
            f"Per-event rank ({rank}) should be >= PDG summary rank ({PDG_RANK})"


class TestFailureModes:

    def test_fm1_uniform_field(self):
        x = {i: np.ones(N_COMPONENTS) for i in range(5)}
        edges = make_edges()
        d, s = operator_e_primary(x, edges)
        mode = detect_failure_mode(d, s, edges)
        assert mode == FailureMode.DIFFERENTIATION_COLLAPSE

    def test_no_failure_gradient_field(self):
        x = make_event_field()
        edges = make_edges()
        d, s = operator_e_primary(x, edges)
        mode = detect_failure_mode(d, s, edges)
        assert mode != FailureMode.DIFFERENTIATION_COLLAPSE
