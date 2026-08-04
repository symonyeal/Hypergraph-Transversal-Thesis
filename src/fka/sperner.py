"""Sperner reduction -- FK-A Step 1 (thesis Algorithm 1, p.13).

A hypergraph is *Sperner* when its edges form an antichain. FK-A reduces both
sides at the top of every call, because contraction can create containments even
when the parent had none (thesis §2.2.2, p.15). FK-B calls the same operation
``Irredundant``.

Split into the thesis' two filters (§5.1.1, p.46) -- ``Unique(f, g)`` and
``remove_superset(np_arr)`` -- so the visualiser can report exactly what each
node discarded.

Two superset removals are provided and the tests require them to agree
everywhere: ``remove_superset_pairwise`` is the definition, and
``remove_superset_sorted`` is the size-ordered pass the legacy README describes.
"""

from __future__ import annotations

from dataclasses import dataclass

from .hypergraph import Hypergraph, popcount

__all__ = [
    "SpernerReduction",
    "unique",
    "remove_superset_pairwise",
    "remove_superset_sorted",
    "sperner_reduce",
]


@dataclass(frozen=True, slots=True)
class SpernerReduction:
    """What a reduction kept and what it threw away.

    The discarded edges are the "Non Sperner Edges" of the thesis' node dumps.
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


def unique(edges: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Split into (kept, duplicates dropped), order-stable (``Unique(f, g)``)."""
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


def remove_superset_pairwise(
    edges: tuple[int, ...],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """The definition: drop ``a`` when some distinct ``b`` has ``b subseteq a``.

    Assumes ``edges`` is already deduplicated -- two equal edges are each a
    non-strict superset of the other, so both would be dropped.
    :func:`sperner_reduce` guarantees the ordering.
    """
    kept: list[int] = []
    dropped: list[int] = []
    for i, a in enumerate(edges):
        if any((b & a) == b for j, b in enumerate(edges) if i != j):
            dropped.append(a)
        else:
            kept.append(a)
    return tuple(kept), tuple(dropped)


def remove_superset_sorted(
    edges: tuple[int, ...],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Superset removal in ascending size order.

    An edge can only contain a *strictly* smaller one -- equal-size containment
    is equality, excluded by deduplication -- so smallest-first means each edge
    is tested only against the minimal edges already kept. That is the legacy
    README's sorting optimisation with the order reversed: keeping the survivors
    as the comparison set is what makes it a single pass.

    Same *set* as :func:`remove_superset_pairwise`; kept edges come back in the
    caller's order so the two are directly comparable.
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
    """Reduce ``H`` to its minimal edges, reporting what went.

    ``method`` picks the superset removal: ``"sorted"`` (the optimised pass) or
    ``"pairwise"`` (the definition). Identical output; the switch exists so the
    tests can assert that, and so a future optimisation is benchmarked against
    the definition rather than against itself.
    """
    kept, dup = unique(H.edges)
    if method == "sorted":
        kept, sup = remove_superset_sorted(kept)
    elif method == "pairwise":
        kept, sup = remove_superset_pairwise(kept)
    else:
        raise ValueError(f"unknown superset removal method {method!r}")
    return SpernerReduction(
        reduced=Hypergraph(H.n, tuple(sorted(kept))),
        removed_duplicate=dup,
        removed_superset=sup,
    )
