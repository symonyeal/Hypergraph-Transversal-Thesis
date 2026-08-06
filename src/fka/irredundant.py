"""``Irredundant`` (Irredundant.m) -- FK-A Step 1, thesis Algorithm 1 p.13.

A monotone Boolean function is irredundant when its terms form an antichain
(hypergraph reading: Sperner). Both algorithms reduce both sides at the top of
every call, because ``phi_x_0`` can create containments even when the parent had
none (thesis §2.2.2, p.15).

Split into the thesis' two filters (§5.1.1, p.46) so a node can report exactly
what it discarded.
"""

from __future__ import annotations

from dataclasses import dataclass

from .hypergraph import Hypergraph, popcount

__all__ = ["Reduction", "unique", "remove_superset", "irredundant"]


@dataclass(frozen=True, slots=True)
class Reduction:
    """What a reduction kept and what it threw away."""

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
    """Split into (kept, duplicates dropped), order-stable."""
    seen: set[int] = set()
    kept: list[int] = []
    dropped: list[int] = []
    for t in edges:
        if t in seen:
            dropped.append(t)
        else:
            seen.add(t)
            kept.append(t)
    return tuple(kept), tuple(dropped)


def remove_superset(
    edges: tuple[int, ...],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Drop ``a`` when some distinct ``b`` has ``b subseteq a``.

    Ascending size order, as Irredundant.m sorts: a term can only contain a
    *strictly* smaller one, so each is tested against the minimal terms already
    kept rather than against all the others. Assumes ``edges`` is deduplicated,
    which :func:`irredundant` guarantees -- two equal terms each contain the
    other, so both would go.
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
    kept = tuple(t for i, t in enumerate(edges) if i not in dropped_idx)
    dropped = tuple(t for i, t in enumerate(edges) if i in dropped_idx)
    return kept, dropped


def irredundant(MBF: Hypergraph) -> Reduction:
    """Reduce ``MBF`` to its minimal terms, reporting what went."""
    kept, dup = unique(MBF.edges)
    kept, sup = remove_superset(kept)
    return Reduction(
        reduced=Hypergraph(MBF.n, tuple(sorted(kept))),
        removed_duplicate=dup,
        removed_superset=sup,
    )
