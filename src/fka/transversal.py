"""``Tr`` by Berge's method (berge.m) -- the independent oracle.

Both algorithms *decide* whether ``cnf = Tr(dnf)`` without ever building
``Tr(dnf)``. Testing them needs something that does build it, by a different
route, or a shared misreading of the conventions would go undetected.

Start from ``Tr({}) = {{}}`` and fold in one term at a time, extending every
partial transversal that does not yet meet it and re-minimising. Exponential in
the worst case -- the output can be -- but exact.

The degenerate conventions fall out of the fold rather than being special-cased:
``Tr({}) = {{}}`` and ``Tr({{}}) = {}``, as in FK96.
"""

from __future__ import annotations

from typing import Iterable, Optional

from .hypergraph import Hypergraph, popcount, verts
from .irredundant import irredundant

__all__ = [
    "transversal",
    "minimal_elements",
    "is_transversal",
    "is_dual_oracle",
    "find_CA",
]


def minimal_elements(sets: Iterable[int]) -> tuple[int, ...]:
    """Keep only the inclusion-minimal bitsets, ascending."""
    items = sorted(set(sets), key=lambda s: (popcount(s), s))
    kept: list[int] = []
    for s in items:
        if not any((k & s) == k for k in kept):
            kept.append(s)
    return tuple(sorted(kept))


def transversal(MBF: Hypergraph) -> Hypergraph:
    """``Tr(MBF)``: the minimal transversals, i.e. the other normal form.

    Terms are folded in shortest-first. A short term multiplies the working set
    by its size, so paying that while the set is still small is markedly cheaper.
    """
    current: set[int] = {0}
    for t in sorted(MBF.edges, key=popcount):
        if t == 0:
            # Nothing can hit the empty term: Tr collapses.
            return Hypergraph.empty(MBF.n)
        nxt: set[int] = set()
        vs = verts(t)
        for s in current:
            if s & t:
                nxt.add(s)
            else:
                for v in vs:
                    nxt.add(s | (1 << v))
        current = set(minimal_elements(nxt))
    return Hypergraph(MBF.n, minimal_elements(current))


def is_transversal(MBF: Hypergraph, S: int) -> bool:
    """True iff ``S`` meets every term of ``MBF``."""
    return all(S & t for t in MBF.edges)


def is_dual_oracle(cnf: Hypergraph, dnf: Hypergraph) -> bool:
    """Decide ``cnf = Tr(dnf)`` by construction, independently of FK-A and FK-B.

    Both sides are reduced first, as the algorithms do at their root, so the
    comparison is between antichains.
    """
    return transversal(irredundant(dnf).reduced).edges == irredundant(cnf).reduced.edges


def find_CA(cnf: Hypergraph, dnf: Hypergraph) -> Optional[int]:
    """Search for a conflicting assignment satisfying eq. 2.1, or ``None``.

    ``S`` must meet every clause and contain no monomial. The second half
    rearranges: ``S`` contains no monomial exactly when ``V \\ S`` meets them
    all. So the search is over ``Tr(cnf)`` whose complement is a transversal of
    ``dnf``.

    Only the *minimal* transversals need testing: if a larger ``S`` works, any
    minimal ``T`` inside it works too, since ``V \\ T`` contains ``V \\ S``. That
    bounds the search by ``|Tr(cnf)|`` rather than ``2^n``.

    ``None`` means no assignment of this form exists, which for a degenerate
    instance such as ``cnf = {{}}`` is not the same as the pair being dual; use
    :func:`is_dual_oracle` for that.
    """
    full = (1 << max(cnf.n, dnf.n)) - 1
    for T in transversal(irredundant(cnf).reduced).edges:
        if is_transversal(dnf, full & ~T):
            return T
    return None
