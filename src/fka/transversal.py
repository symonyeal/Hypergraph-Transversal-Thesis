"""Ground-truth transversal computation.

FK-A is a *decision* procedure: it answers "is ``G = Tr(H)``?" without ever
building ``Tr(H)``. To test it we need an independent oracle that does build it,
by a completely different route -- otherwise a shared misunderstanding of the
conventions would go undetected.

``transversal`` uses Berge's sequential method: start from ``Tr(empty) = {{}}``
and fold in one hyperedge at a time, extending every partial transversal that
does not yet meet the new edge and re-minimising. This is exponential in the
worst case (the output itself can be), but it is exact and fast enough for the
instance sizes in the thesis.

Degenerate conventions fall out of the fold rather than being special-cased,
and they agree with the monotone-Boolean reading used in FK96:

* ``Tr({})   = {{}}``  -- the constant 0 dualises to the constant 1;
* ``Tr({{}}) = {}``    -- and back again.
"""

from __future__ import annotations

from typing import Iterable, Optional

from .hypergraph import Hypergraph, bitset_to_vertices, popcount
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
    """``Tr(H)``: the hypergraph of minimal transversals (minimal hitting sets).

    Edges of ``H`` are folded in shortest-first, which keeps the intermediate
    families small -- a short edge multiplies the working set by its size, so
    paying that cost while the set is still small is markedly cheaper.
    """
    current: set[int] = {0}
    for e in sorted(H.edges, key=popcount):
        if e == 0:
            # The empty edge cannot be hit by anything: Tr collapses to nothing.
            return Hypergraph.empty(H.n)
        nxt: set[int] = set()
        vertices = bitset_to_vertices(e)
        for t in current:
            if t & e:
                nxt.add(t)
            else:
                for v in vertices:
                    nxt.add(t | (1 << v))
        current = set(minimal_elements(nxt))
    return Hypergraph(H.n, minimal_elements(current))


def is_transversal(H: Hypergraph, S: int) -> bool:
    """True iff the vertex set ``S`` meets every edge of ``H``."""
    return all(S & e for e in H.edges)


def is_dual_oracle(G: Hypergraph, H: Hypergraph) -> bool:
    """Decide ``G = Tr(H)`` by construction, independently of FK-A.

    Both sides are Sperner-reduced first, matching what FK-A does at its root,
    so the comparison is between antichains.
    """
    Gr = sperner_reduce(G).reduced
    Hr = sperner_reduce(H).reduced
    return transversal(Hr).edges == Gr.edges


def find_certificate(G: Hypergraph, H: Hypergraph) -> Optional[int]:
    """Search for a set ``S`` satisfying thesis equation 2.1, or ``None``.

    ``S`` must meet every edge of ``G`` and contain no edge of ``H``. The second
    half rearranges: ``S`` contains no edge of ``H`` exactly when ``V \\ S``
    meets every edge of ``H``. So the search is over transversals of ``G`` whose
    complement is a transversal of ``H``.

    It suffices to test the *minimal* transversals of ``G``: if a larger
    transversal ``S`` works, any minimal ``T`` inside it works too, since
    ``V \\ T`` contains ``V \\ S``. That bounds the search by ``|Tr(G)|`` rather
    than ``2^n``.

    Returning ``None`` means no certificate of this form exists -- which for
    degenerate instances such as ``G = {{}}`` is not the same as the pair being
    dual. Use :func:`is_dual_oracle` to decide duality.
    """
    full = (1 << max(G.n, H.n)) - 1
    for T in transversal(sperner_reduce(G).reduced).edges:
        if is_transversal(H, full & ~T):
            return T
    return None
