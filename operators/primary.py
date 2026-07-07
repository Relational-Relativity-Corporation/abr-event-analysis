"""
operators/primary.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

V6 primary operator kernel: Delta (Δ) and Sigma (Σ).

Grounding: operators_primary.rs, primary_operators_delta_sigma_v6.md.

Δ(x)[e] = x[s] - x[t]  — directed difference
Σ(g)[e] = g[e] + rho[e] * (forward_sum - backward_sum)  — local antisymmetry

No B accumulation. No path structure assumed.
Every declared relation has exactly one admissible direction.
Symmetric declared edge-images signal inadmissible structure.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional


SVD_TOLERANCE = 1e-10  # declared before execution, not adjusted post-hoc


def operator_delta(
    x: Dict[int, np.ndarray],
    edges: List[Tuple[int, int]],
) -> np.ndarray:
    """
    Delta(x)[e] = x[s] - x[t] for each declared edge e = (s, t).

    Args:
        x: dict mapping event_id -> observable vector (numpy array)
        edges: list of (source, target) pairs — declared directed relations

    Returns:
        delta_field: numpy array, shape (n_edges, n_components)

    Constraint: directed difference only. No relation and no direction may
    be added that the declaration did not trace to an observable through M.
    Admissible before Delta: uniform shift, declared unit scale.
    """
    if len(edges) == 0:
        return np.empty((0, 0))
    return np.array([x[s] - x[t] for s, t in edges], dtype=float)


def build_adjacency(
    edges: List[Tuple[int, int]],
) -> Tuple[List[List[int]], List[List[int]]]:
    """
    Build adj+ and adj- for Sigma.

    adj+[e]: indices of edges whose source = target of e
    adj-[e]: indices of edges whose target = source of e

    Constraint: adjacency is declared, not assumed. Only edges present
    in the declared edge list appear here.
    """
    n_edges = len(edges)
    adj_plus  = [[] for _ in range(n_edges)]
    adj_minus = [[] for _ in range(n_edges)]

    for e, (_, et) in enumerate(edges):
        for f, (fs, _) in enumerate(edges):
            if fs == et:
                adj_plus[e].append(f)

    for e, (es, _) in enumerate(edges):
        for f, (_, ft) in enumerate(edges):
            if ft == es:
                adj_minus[e].append(f)

    return adj_plus, adj_minus


def compute_rho(
    delta_field: np.ndarray,
    edges: List[Tuple[int, int]],
    rho_base: float = 0.2,
) -> np.ndarray:
    """
    rho[e] = rho_base * m[s] / (1 + m[s])
    m[s] = max |delta_field[e']| for e' incident to source node of e.

    Derived per edge from delta_field. No aggregation beyond the node.
    rho in [0, rho_base) for all declared edges.
    """
    if len(edges) == 0:
        return np.array([])

    n_edges = len(edges)
    node_incident: Dict[int, List[int]] = {}
    for e_idx, (s, t) in enumerate(edges):
        node_incident.setdefault(s, []).append(e_idx)
        node_incident.setdefault(t, []).append(e_idx)

    rho = np.zeros(n_edges)
    for e_idx, (s, _) in enumerate(edges):
        incident = node_incident.get(s, [])
        if incident:
            m_s = max(np.max(np.abs(delta_field[i])) for i in incident)
        else:
            m_s = 0.0
        rho[e_idx] = rho_base * m_s / (1.0 + m_s)

    return rho


def operator_sigma(
    delta_field: np.ndarray,
    edges: List[Tuple[int, int]],
    rho_base: float = 0.2,
) -> np.ndarray:
    """
    Sigma(g)[e] = g[e] + rho[e] * (adj+ sum - adj- sum)

    g = delta_field — output of Delta applied to x.
    Sigma receives Delta(x) directly. No B accumulation.
    Immediate adjacency only.

    Constraint: immediate adjacency only. No path accumulation.
    Direction determined by the observable — not assumed.
    """
    if len(edges) == 0:
        return np.empty((0, 0))

    adj_plus, adj_minus = build_adjacency(edges)
    rho = compute_rho(delta_field, edges, rho_base)
    result = delta_field.copy()

    for e_idx in range(len(edges)):
        n_comp = delta_field.shape[1]
        forward = (
            np.sum([delta_field[f] for f in adj_plus[e_idx]], axis=0)
            if adj_plus[e_idx] else np.zeros(n_comp)
        )
        backward = (
            np.sum([delta_field[b] for b in adj_minus[e_idx]], axis=0)
            if adj_minus[e_idx] else np.zeros(n_comp)
        )
        result[e_idx] += rho[e_idx] * (forward - backward)

    return result


def antisymmetric_term(
    delta_field: np.ndarray,
    edges: List[Tuple[int, int]],
    rho_base: float = 0.2,
) -> np.ndarray:
    """
    Extracts the antisymmetric term of Sigma separately from pass-through.

    antisymmetric_term[e] = rho[e] * (adj+ sum - adj- sum)

    Declared projection:
      Preserves: circulation contribution per declared edge.
      Discards:  pass-through term g[e].
    """
    if len(edges) == 0:
        return np.empty((0, 0))

    adj_plus, adj_minus = build_adjacency(edges)
    rho = compute_rho(delta_field, edges, rho_base)
    asym = np.zeros_like(delta_field)

    for e_idx in range(len(edges)):
        n_comp = delta_field.shape[1]
        forward = (
            np.sum([delta_field[f] for f in adj_plus[e_idx]], axis=0)
            if adj_plus[e_idx] else np.zeros(n_comp)
        )
        backward = (
            np.sum([delta_field[b] for b in adj_minus[e_idx]], axis=0)
            if adj_minus[e_idx] else np.zeros(n_comp)
        )
        asym[e_idx] = rho[e_idx] * (forward - backward)

    return asym


def operator_e_primary(
    x: Dict[int, np.ndarray],
    edges: List[Tuple[int, int]],
    rho_base: float = 0.2,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    E_primary(x) = Sigma(Delta(x))

    Returns (delta_field, sigma_field).
    Both preserved for analysis — provenance maintained.
    """
    delta = operator_delta(x, edges)
    sigma = operator_sigma(delta, edges, rho_base)
    return delta, sigma
