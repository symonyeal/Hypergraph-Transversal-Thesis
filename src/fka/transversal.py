"""``Tr(H)`` by Berge's method -- the independent oracle.

Both algorithms *decide* whether ``G = Tr(H)`` without ever building ``Tr(H)``.
Testing them needs something that does build it, by a different route, or a
shared misreading of the conventions would go undetected.

Berge's sequential method: start from ``Tr({}) = {{}}`` and fold in one edge at
a time, extending every partial transversal that does not yet meet it and
re-minimising. Exponential in the worst case -- the output can be -- but exact,
and fast enough for the sizes here.

The degenerate conventions fall out of the fold rather than being special-cased,
and agree with the monotone-Boolean reading used throughout FK96:
``Tr({}) = {{}}`` and ``Tr({{}}) = {}``.
"""

from __future__ import annotations

from typing import Iterable, Optional

from .hypergraph import Hypergraph, popcount, verts
from .sperner import sperner_reduce

__all__ = [
    "transversal",
    "minimal_elements",
    "is_transversal",
    "is_dual_oracle",
    "find_certificate",
]


def minimal_elements(sets: Iterable[int]) -> tuple[int, ...]:
    """Keep only the inclusion-minimal bitsets, ascending."""
    items = sorted(set(sets), key=lambda s: (popcount(s), s))
    kept: list[int] = []
    for s in items:
        if not any((k & s) == k for k in kept):
            kept.append(s)
    return tuple(sorted(kept))


def transversal(H: Hypergraph) -> Hypergraph:
    """``Tr(H)``: the minimal transversals (minimal hitting sets).

    Edges are folded in shortest-first. A short edge multiplies the working set
    by its size, so paying that while the set is still small is markedly cheaper.
    """
    current: set[int] = {0}
    for e in sorted(H.edges, key=popcount):
        if e == 0:
            # Nothing can hit the empty edge: Tr collapses.
            return Hypergraph.empty(H.n)
        nxt: set[int] = set()
        vs = verts(e)
        for t in current:
            if t & e:
                nxt.add(t)
            else:
                for v in vs:
                    nxt.add(t | (1 << v))
        current = set(minimal_elements(nxt))
    return Hypergraph(H.n, minimal_elements(current))


def is_transversal(H: Hypergraph, S: int) -> bool:
    """True iff ``S`` meets every edge of ``H``."""
    return all(S & e for e in H.edges)


def is_dual_oracle(G: Hypergraph, H: Hypergraph) -> bool:
    """Decide ``G = Tr(H)`` by construction, independently of FK-A and FK-B.

    Both sides are Sperner-reduced first, as the algorithms do at their root, so
    the comparison is between antichains.
    """
    return transversal(sperner_reduce(H).reduced).edges == sperner_reduce(G).reduced.edges


def find_certificate(G: Hypergraph, H: Hypergraph) -> Optional[int]:
    """Search for an ``S`` satisfying thesis equation 2.1, or ``None``.

    ``S`` must meet every edge of ``G`` and contain no edge of ``H``. The second
    half rearranges: ``S`` contains no edge of ``H`` exactly when ``V \\ S``
    meets every edge of ``H``. So the search is over transversals of ``G`` whose
    complement is a transversal of ``H``.

    Only the *minimal* transversals of ``G`` need testing: if a larger ``S``
    works, any minimal ``T`` inside it works too, since ``V \\ T`` contains
    ``V \\ S``. That bounds the search by ``|Tr(G)|`` rather than ``2^n``.

    ``None`` means no certificate of this form exists, which for a degenerate
    instance such as ``G = {{}}`` is not the same as the pair being dual. Use
    :func:`is_dual_oracle` to decide duality.
    """
    full = (1 << max(G.n, H.n)) - 1
    for T in transversal(sperner_reduce(G).reduced).edges:
        if is_transversal(H, full & ~T):
            return T
    return None
