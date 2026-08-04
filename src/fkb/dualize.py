"""Dualisation by repeated equivalence testing.

Port of ``FK_Dualization.m`` with ``Maximum_False_Point.m``,
``Minimality_Check.m`` and ``Initialize_CNF.m``.

This is what the FK-B reference actually uses FK-B *for*. :func:`fkb.fk_b` only
decides whether a given ``C`` equals a given ``D``; here that decision is turned
into a generator of ``Tr(D)`` by the standard self-reduction:

1. start from a small guess at ``C``;
2. ask FK-B whether ``C = D``. No conflicting assignment means ``C`` is already
   ``Tr(D)`` and the loop is done;
3. otherwise the assignment ``S`` has ``C(S) = 1`` and ``D(S) = 0``. Grow it to
   a *maximal* false point of ``D`` -- as many variables true as possible while
   ``D`` stays false -- and add its complement to ``C`` as a new clause;
4. repeat.

The complement of a maximal false point is a minimal transversal of ``D``, so
every pass adds a genuine member of ``Tr(D)`` and none is ever added twice.
That bounds the loop by ``|Tr(D)|`` passes.

Compared with :func:`fka.transversal.transversal`, which folds monomials in one
at a time by Berge's method, this pays one equivalence test per output edge.
Berge is far faster on the instances in this project and remains the oracle; the
value of this route is that it is *incremental* -- it produces correct partial
output and can be stopped early -- and that it exercises FK-B on the sequence of
subproblems the reference benchmarks time.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from fka.hypergraph import Hypergraph, verts
from fka.sperner import sperner_reduce

from .algorithm import cnf_value, dnf_value, fk_b

__all__ = ["Dualisation", "dualise", "maximum_false_point", "initial_cnf"]


def _greedy_cover(dnf: Hypergraph, allowed: int) -> Optional[int]:
    """A transversal of ``dnf`` drawn only from the variables in ``allowed``.

    Repeatedly takes the allowed variable lying in the most surviving monomials
    and drops the monomials it meets -- the greedy set cover that
    ``Maximum_False_Point.m`` and ``Initialize_CNF.m`` both run, differing only
    in the tie-break. ``None`` when some monomial has no allowed variable at
    all, so no transversal of this shape exists.
    """
    left = list(dnf.edges)
    picked = 0
    while left:
        counts: dict[int, int] = {}
        for m in left:
            for v in verts(m & allowed):
                counts[v] = counts.get(v, 0) + 1
        if not counts:
            return None
        best = max(counts, key=lambda v: (counts[v], -v))
        picked |= 1 << best
        left = [m for m in left if not m & (1 << best)]
    return picked


def _minimise(cover: int, dnf: Hypergraph) -> int:
    """Drop variables ``cover`` does not need (Minimality_Check.m).

    Highest index first, matching the MATLAB's descending loop, so the result is
    reproducible rather than merely minimal.
    """
    for v in reversed(verts(cover)):
        without = cover & ~(1 << v)
        if all(m & without for m in dnf.edges):
            cover = without
    return cover


def maximum_false_point(dnf: Hypergraph, S: int) -> int:
    """Grow ``S`` to a maximal set on which ``dnf`` is still false.

    ``Maximum_False_Point.m``: find a minimal set of variables whose being false
    already falsifies every monomial, and set everything else true. The
    complement of the result is therefore a *minimal* transversal of ``dnf``,
    which is what makes it a legitimate new clause.

    ``S`` fixes which variables stay true: the cover is drawn from outside it,
    so nothing ``S`` already sets is turned off. ``S`` must satisfy
    ``dnf(S) = 0``, which is exactly what FK-B hands back, and that is also what
    guarantees the cover exists -- no monomial can be a subset of ``S``, so
    every one of them still has a variable available.
    """
    if dnf_value(dnf, S):
        raise ValueError("the seed must be a false point of the DNF")
    ground = dnf.support() | S
    off = _greedy_cover(dnf, allowed=ground & ~S)
    if off is None:  # pragma: no cover -- excluded by dnf(S) = 0
        raise AssertionError("a false point always leaves every monomial coverable")
    return ground & ~_minimise(off, dnf)


def initial_cnf(dnf: Hypergraph, n_clause: int = 1) -> Hypergraph:
    """A starting ``C`` of a few real members of ``Tr(dnf)`` (Initialize_CNF.m).

    The MATLAB draws them at random. Here the greedy cover is deterministic, so
    successive clauses are forced apart by blocking the variables already used,
    and the search stops as soon as no further transversal avoids them --
    ``n_clause`` is an upper bound, not a promise. A deterministic start makes a
    dualisation run reproducible, which the random one is not.
    """
    clauses: list[int] = []
    blocked = 0
    while len(clauses) < n_clause:
        cover = _greedy_cover(dnf, allowed=dnf.support() & ~blocked)
        if cover is None:
            break
        clause = _minimise(cover, dnf)
        if clause in clauses:
            break
        clauses.append(clause)
        blocked |= clause
    return sperner_reduce(Hypergraph.from_bitsets(dnf.n, clauses)).reduced


@dataclass(slots=True)
class Dualisation:
    """The outcome of a dualisation run, and its trace.

    ``sizes`` is ``|C|`` at the start of each pass -- the curve
    ``Computing_CPUtimes.m`` plots to compare FK, MFK and IMFK. ``nodes`` is the
    FK-B recursion-tree size of each equivalence test, which is where the work
    actually goes.
    """

    cnf: Hypergraph
    complete: bool
    passes: int = 0
    sizes: list[int] = field(default_factory=list)
    nodes: list[int] = field(default_factory=list)

    @property
    def total_nodes(self) -> int:
        return sum(self.nodes)


def dualise(
    dnf: Hypergraph,
    *,
    variant: str = "faithful",
    split_rule: str = "most_frequent",
    max_passes: Optional[int] = None,
) -> Dualisation:
    """Build ``Tr(dnf)`` by repeated FK-B equivalence tests.

    ``max_passes`` stops early and reports ``complete=False``; the partial
    ``cnf`` is still a subset of ``Tr(dnf)``, never a wrong answer.
    """
    dnf = sperner_reduce(dnf).reduced
    if not dnf.edges:
        # Tr of the constant 0 is the constant 1: one empty clause.
        return Dualisation(cnf=Hypergraph(dnf.n, (0,)), complete=True)
    if 0 in dnf.edges:
        return Dualisation(cnf=Hypergraph.empty(dnf.n), complete=True)

    cnf = initial_cnf(dnf)
    run = Dualisation(cnf=cnf, complete=False)
    while max_passes is None or run.passes < max_passes:
        run.sizes.append(len(cnf))
        tree = fk_b(cnf, dnf, variant=variant, split_rule=split_rule)
        run.nodes.append(len(tree))
        run.passes += 1
        if tree.dual:
            run.cnf, run.complete = cnf, True
            return run

        S = tree.verdict.certificate
        if S is None or not (cnf_value(cnf, S) and not dnf_value(dnf, S)):
            # Only condition 1 yields the other polarity, and it cannot fire
            # here: every clause of C is a transversal of D by construction.
            raise AssertionError(
                f"dualisation expected a C-true/D-false assignment, got {S!r}"
            )
        clause = dnf.support() & ~maximum_false_point(dnf, S)
        if cnf.contains_edge(clause):
            raise AssertionError(
                "dualisation repeated a clause; the maximal false point is not maximal"
            )
        cnf = Hypergraph.from_bitsets(cnf.n, cnf.edges + (clause,))
        run.cnf = cnf
    return run
