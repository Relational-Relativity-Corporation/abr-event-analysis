"""
declaration/relations.py
Metatron Dynamics, Inc. V6. Bounded over D. No claim beyond D.

Declares the relational structure over the event field.

Each relation is a directed edge (s, t) where s and t are event indices.
Direction must be traceable to an observable through M.
No relation is fabricated. No direction is assumed.

Declared relation types for the event field:
  - temporal: event s precedes event t in declared production order
  - spatial: events s and t share a declared spatial proximity
  - kinematic: events s and t share a declared kinematic relation

Origin declares which relation type is admissible for each analysis.
"""

from dataclasses import dataclass
from enum import Enum
from typing import List, Tuple
from .admissibility import AdmissibilityError


class RelationProvenance(Enum):
    """
    Every declared relation must carry one of these provenances.
    Provenance is not annotation — it is the admissibility condition.
    """
    TEMPORAL   = "temporal"    # direction: earlier event -> later event
    SPATIAL    = "spatial"     # direction: declared from spatial proximity observable
    KINEMATIC  = "kinematic"   # direction: declared from kinematic observable


@dataclass
class DeclaredRelation:
    source: int              # event index
    target: int              # event index
    provenance: RelationProvenance
    direction_basis: str     # what observable determines the direction


def declare_relations(
    event_ids: List[int],
    relation_type: RelationProvenance,
    direction_basis: str,
) -> List[Tuple[int, int]]:
    """
    Declare a set of relations over the event field.

    For temporal relations: each event is connected to the next in
    declared production order. Direction is s -> t where s precedes t.

    Args:
        event_ids: ordered list of event indices
        relation_type: the provenance of all declared relations
        direction_basis: the observable that determines direction

    Returns:
        list of (source, target) tuples — the declared edge set

    Raises:
        AdmissibilityError if relation_type is not declared or
        direction_basis is empty.
    """
    if not direction_basis.strip():
        raise AdmissibilityError(
            "direction_basis must state what observable determines the "
            "direction of each declared relation. Empty string is inadmissible "
            "— declaration is not a substitute for traceability."
        )

    if len(event_ids) < 2:
        raise AdmissibilityError(
            "At least two events required to declare a relation. "
            "A single event has no declared relation to act over."
        )

    if relation_type == RelationProvenance.TEMPORAL:
        # Strictly directed: earlier -> later.
        # No reverse edge — the reverse direction would require declaring
        # that a later event precedes an earlier one, which is not traceable
        # to an observable through M under the directional admissibility condition.
        edges = [(event_ids[i], event_ids[i+1])
                 for i in range(len(event_ids) - 1)]
    else:
        raise AdmissibilityError(
            f"Relation type '{relation_type}' not yet implemented. "
            f"Origin must declare the specific observable that determines "
            f"direction before this relation type is used."
        )

    return edges
