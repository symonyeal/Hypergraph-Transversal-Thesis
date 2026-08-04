"""SageMath-backed implementations.

Imported only when :func:`fka.backends.has_sage` is true. Nothing else in the
package imports this module at load time, so ``import fka`` stays fast and
works on a plain CPython install.

Every function here has a pure-Python counterpart that is used otherwise; these
exist to (a) give GAP's authoritative group names and (b) serve as an
independent implementation to cross-check the built-in search against. The test
suite runs the cross-checks only when Sage is importable and skips them
otherwise.
"""

from __future__ import annotations

from typing import Any, Sequence

from ..groups import GroupIdentification, PermutationGroup, cycle_notation
from ..hypergraph import Hypergraph, verts

__all__ = [
    "name_group",
    "sage_automorphism_group",
    "to_sage_matrix",
    "to_incidence_structure",
    "sage_graph_classes",
]


def _require_sage():
    try:
        import sage.all as sage  # noqa: F401
    except ImportError as exc:  # pragma: no cover - only on non-Sage installs
        raise RuntimeError(
            "the SageMath backend was called but sage.all is not importable; "
            "run under `sage -python`"
        ) from exc
    return sage


def to_sage_matrix(H: Hypergraph):
    """The incidence matrix of ``H`` as a Sage integer matrix."""
    sage = _require_sage()
    return sage.Matrix(sage.ZZ, H.to_incidence())


def to_incidence_structure(H: Hypergraph, *, include_isolated: bool = False):
    """``H`` as a Sage :class:`IncidenceStructure`.

    Points are the *active* vertices by default. Passing
    ``include_isolated=True`` keeps the full ground set, which lets isolated
    vertices permute among themselves -- see :mod:`fka.automorphism` for why
    that is not the default.
    """
    sage = _require_sage()
    from sage.combinat.designs.incidence_structures import IncidenceStructure

    target, mapping = (H, list(range(H.n))) if include_isolated else H.compact()
    points = list(range(target.n))
    blocks = [list(verts(e)) for e in target.edges]
    if not blocks:
        # IncidenceStructure rejects an empty block list on some Sage versions.
        raise ValueError("hypergraph has no edges; Aut is the full symmetric group")
    return IncidenceStructure(points, blocks), mapping


def name_group(group: PermutationGroup) -> GroupIdentification:
    """Name a permutation group using GAP's ``StructureDescription``.

    Falls back to the pure-Python catalogue if GAP raises -- some structure
    descriptions are genuinely expensive, and a missing name should not abort
    an experiment that already has the order and generators.
    """
    sage = _require_sage()
    from ..groups import identify

    if group.order == 1:
        return GroupIdentification("C1", 1, True, "sage-trivial")
    try:
        gens = [cycle_notation(g, one_indexed=True) for g in (group.generators or [])]
        gens = [g for g in gens if g != "()"]
        if not gens:
            return GroupIdentification("C1", 1, True, "sage-trivial")
        G = sage.PermutationGroup([sage.PermutationGroupElement(g) for g in gens])
    except Exception:
        return identify(group)

    sage_order = int(G.order())
    if sage_order != group.order:
        # Disagreement means the generating set was mis-transcribed; that is a
        # correctness bug, not a reason to silently relabel with the fallback.
        raise ValueError(
            f"Sage read the generators as a group of order {sage_order} "
            f"but the search found {group.order} automorphisms"
        )
    try:
        name = str(G.structure_description())
        return GroupIdentification(name, group.order, bool(G.is_abelian()), "sage-gap")
    except Exception:
        return identify(group)


def _sage_permutation_label(perm, mapping: Sequence[int]) -> str:
    """Translate Sage cycles on compact 0-based points to thesis labels.

    ``IncidenceStructure`` preserves the point labels supplied to it.  Here
    those labels are ``0..n-1``; subtracting one (as the legacy adapter did)
    rotated every generator and mapped point 0 to the last vertex.
    """
    cycles = [
        "(" + ",".join(str(mapping[int(i)] + 1) for i in cyc) + ")"
        for cyc in perm.cycle_tuples()
        if len(cyc) > 1
    ]
    return "".join(cycles) or "()"


def sage_automorphism_group(H: Hypergraph, *, include_isolated: bool = False) -> dict[str, Any]:
    """``Aut(H)`` via Sage's ``IncidenceStructure``, for cross-checking.

    Returns order, structure name and generators in cycle notation on the
    caller's 1-indexed vertex numbering.
    """
    _require_sage()
    structure, mapping = to_incidence_structure(H, include_isolated=include_isolated)
    A = structure.automorphism_group()
    gens = [_sage_permutation_label(perm, mapping) for perm in A.gens()]
    return {
        "order": int(A.order()),
        "name": str(A.structure_description()) if int(A.order()) > 1 else "C1",
        "generators": gens,
        "backend": "sage",
    }


def sage_graph_classes(edges: Sequence[tuple[int, int]], n: int) -> dict[str, bool]:
    """Graph-class recognition using Sage's polynomial-time algorithms.

    The pure-Python versions in :mod:`fka.graphs` enumerate induced subgraphs
    and are capped at 20 vertices; this has no such limit.
    """
    sage = _require_sage()
    G = sage.Graph([list(range(n)), list(edges)], format="vertices_and_edges")
    return {
        "cograph": bool(G.is_cograph()),
        "chordal": bool(G.is_chordal()),
        "bipartite": bool(G.is_bipartite()),
        "split": bool(G.is_split()),
        "weakly_chordal": bool(G.is_weakly_chordal()),
        "perfect": bool(G.is_perfect()),
        "chordal_bipartite": bool(G.is_chordal_bipartite()),
        "triangle_free": bool(G.is_triangle_free()),
    }
