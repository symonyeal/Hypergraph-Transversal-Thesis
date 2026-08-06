"""Monotone composition ``T[B_1..B_t]``, and Berge's identities for it.

Thesis Section 4.1 gives the two operators used to build the alternating trees:
``A v B`` (term-set union) and ``A ^ B`` (all ``s u t``), with
``Tr(A v B) = Tr(A) ^ Tr(B)`` and ``Tr(A ^ B) = Tr(A) v Tr(B)`` when the
operands are variable-disjoint [Ber87].

Both are the extreme cases of one operation. Substituting a hypergraph ``B_i``
for variable ``i`` of a template ``T`` -- monotone Boolean composition -- gives
``v`` when ``T`` is a set of singletons and ``^`` when ``T`` is a single edge,
and the two identities become the single law

    Tr(T[B_1..B_t]) = Tr(T)[Tr(B_1)..Tr(B_t)].

Writing the operators this way is what makes the whole GK/TK/Hertel/Mark zoo one
object: a family is a template plus a seed, and its transversals never have to
be computed, only composed. :mod:`fka.transversal` stays the independent oracle
the tests check this against.

Every function assumes variable-disjoint blocks, which is what makes the
identities exact rather than "up to Irredundant".
"""

from __future__ import annotations

from typing import Sequence

from .hypergraph import Hypergraph, verts

__all__ = ["shift", "substitute", "disjoint_vee", "disjoint_wedge", "single"]


def shift(MBF: Hypergraph, offset: int, n_total: int) -> Hypergraph:
    """Re-place ``MBF`` on variables ``offset .. offset + MBF.n - 1`` of ``n_total``."""
    return Hypergraph.from_bitsets(n_total, (t << offset for t in MBF.edges))


def substitute(T: Hypergraph, blocks: Sequence[Hypergraph]) -> Hypergraph:
    """``T[B_1,...,B_t]``: variable ``i`` of ``T`` replaced by ``blocks[i]``.

    ``E = { union_{i in e} f_i : e in E(T), f_i in E(B_i) }`` on the disjoint
    union of the blocks' variable sets. Irredundant if ``T`` and every block is.
    """
    if len(blocks) != T.n:
        raise ValueError(f"T has {T.n} vertices but {len(blocks)} blocks given")
    offsets, total = [], 0
    for B in blocks:
        offsets.append(total)
        total += B.n
    placed = [shift(B, off, total) for B, off in zip(blocks, offsets)]
    out: list[int] = []
    for e in T.edges:
        acc = [0]
        for i in verts(e):
            acc = [a | f for a in acc for f in placed[i].edges]
        out.extend(acc)
    return Hypergraph.from_bitsets(total, out)


def disjoint_vee(A: Hypergraph, B: Hypergraph) -> Hypergraph:
    """``A v B`` on disjoint supports (thesis 4.1).

    Distinct from :meth:`Hypergraph.vee`, which unions term sets on one shared
    ground set; here the two are shifted apart first.
    """
    return substitute(Hypergraph.from_sets(2, [[0], [1]]), [A, B])


def disjoint_wedge(A: Hypergraph, B: Hypergraph) -> Hypergraph:
    """``A ^ B`` on disjoint supports (thesis 4.1)."""
    return substitute(Hypergraph.from_sets(2, [[0, 1]]), [A, B])


def single() -> Hypergraph:
    """``{{v}}``: one variable, one term. Self-dual, and the trivial seed."""
    return Hypergraph.from_sets(1, [[0]])
