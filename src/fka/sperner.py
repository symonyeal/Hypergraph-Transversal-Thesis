"""Sperner reduction -- Step 1 of FK-A (thesis Algorithm 1, p.13).

A hypergraph is *Sperner* when its edge set is an antichain. FK-A reduces both
inputs to Sperner form at the top of every recursive call, because contraction
(``H_0``) can create edges that contain one another even when the parent was
Sperner (thesis Section 2.2.2, p.15).

The thesis' implementation splits this into two filters (Section 5.1.1, p.46):

* ``Unique(f, g)`` -- duplicate edge removal;
* ``remove_superset(np_arr)`` -- drop any edge that contains another.

Both are kept here as separate steps so the visualiser can report exactly which
edges each node discarded, which is the point of the post-decomposition filter
described in the thesis.

Two implementations of superset removal are provided and are required by the
test suite to agree on every instance:

``superset_reduce_pairwise``
    The direct definition: compare every ordered pair. This is the reference.

``superset_reduce_sorted``
    The optimisation described in the legacy README -- process edges in size
    order so that each edge is only tested against strictly smaller candidates.
    Same output, fewer comparisons.
"""

from __future__ import annotations

from dataclasses import dataclass

from .hypergraph import Hypergraph, popcount

__all__ = [
    "SpernerReduction",
    "dedupe",
    "superset_reduce_pairwise",
    "superset_reduce_sorted",
    "sperner_reduce",
]


@dataclass(frozen=True, slots=True)
class SpernerReduction:
    """Result of reducing a hypergraph to Sperner form.

    ``removed_duplicate`` and ``removed_superset`` are the discarded edges, kept
    for the recursion-tree visualisation ("Non Sperner Edges" in the thesis'
    node dumps). ``removed`` is their concatenation.
    """

    reduced: Hypergraph
    removed_duplicate: tuple[int, ...]
    removed_superset: tuple[int, ...]

    @property
    def removed(self) -> tuple[int, ...]:
        return self.removed_duplicate + self.removed_superset

    @property
    def changed(self) -> bool:
        return bool(self.removed)


def dedupe(edges: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Split ``edges`` into (unique edges, duplicates dropped).

    Order-stable so the dropped list is reproducible.
    """
    seen: set[int] = set()
    kept: list[int] = []
    dropped: list[int] = []
    for e in edges:
        if e in seen:
            dropped.append(e)
        else:
            seen.add(e)
            kept.append(e)
    return tuple(kept), tuple(dropped)


def superset_reduce_pairwise(edges: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Reference superset removal: drop ``a`` whenever some distinct ``b`` has ``b subseteq a``.

    Assumes ``edges`` is already deduplicated -- otherwise two equal edges would
    each be a (non-strict) superset of the other and *both* would be dropped.
    ``sperner_reduce`` guarantees this ordering.
    """
    kept: list[int] = []
    dropped: list[int] = []
    for i, a in enumerate(edges):
        if any((b & a) == b for j, b in enumerate(edges) if i != j):
            dropped.append(a)
        else:
            kept.append(a)
    return tuple(kept), tuple(dropped)


def superset_reduce_sorted(edges: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Superset removal in ascending size order.

    An edge can only contain a *strictly smaller* one (equal-size containment
    would mean equality, excluded by deduplication), so processing smallest
    first means each edge is tested only against the minimal edges already
    kept. That is the sorting optimisation the legacy README describes, with
    the ordering reversed: keeping the survivors as the comparison set is what
    makes the pass single, rather than sorting descending and scanning columns.

    Returns the same *set* as ``superset_reduce_pairwise``; the kept edges come
    back in the caller's original order so the two are directly comparable.
    """
    order = sorted(range(len(edges)), key=lambda i: (popcount(edges[i]), edges[i]))
    minimal: list[int] = []
    dropped_idx: set[int] = set()
    for i in order:
        a = edges[i]
        if any((b & a) == b for b in minimal):
            dropped_idx.add(i)
        else:
            minimal.append(a)
    kept = tuple(e for i, e in enumerate(edges) if i not in dropped_idx)
    dropped = tuple(e for i, e in enumerate(edges) if i in dropped_idx)
    return kept, dropped


def sperner_reduce(H: Hypergraph, *, method: str = "sorted") -> SpernerReduction:
    """Reduce ``H`` to its minimal edges, reporting what was removed.

    ``method`` selects the superset-removal implementation: ``"sorted"``
    (default, the optimised pass) or ``"pairwise"`` (the reference). They
    produce identical output; the switch exists so the test suite can assert
    that, and so a future optimisation can be benchmarked against the
    definition rather than against itself.
    """
    unique, dup = dedupe(H.edges)
    if method == "sorted":
        kept, sup = superset_reduce_sorted(unique)
    elif method == "pairwise":
        kept, sup = superset_reduce_pairwise(unique)
    else:
        raise ValueError(f"unknown superset removal method {method!r}")
    return SpernerReduction(
        reduced=Hypergraph(H.n, tuple(sorted(kept))),
        removed_duplicate=dup,
        removed_superset=sup,
    )
