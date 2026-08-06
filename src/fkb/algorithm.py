"""The Fredman-Khachiyan algorithm FK-B.

A port of the MATLAB reference in ``FK- B/FK-master/src`` (``FK_B.m``,
``MFK_B.m``, ``Check_Conditions.m``, ``Easy_case.m``, ``Choose_SplitVar.m``,
``mu_*.m``), following the listing in *How to apply SAT-solving for the
equivalence test of monotone normal forms*, p.5. Line-number comments below cite
that listing. The restrictions ``phi_x_0``, ``phi_x_1``, ``A_c_x`` and ``A_m_x``
are shared with FK-A and live in :mod:`fka.hypergraph`.

FK-B decides the same question as FK-A on the same pair: ``C = D`` exactly when
``cnf = Tr(dnf)``. What it returns when they differ is a conflicting assignment
``S`` with ``C(S) != D(S)``.

The split is where the two part company. FK-A always recurses twice. FK-B first
asks whether ``x`` is at most mu-frequent on one side; if it is, it recurses
once on that side's restriction and then once per clause (or per monomial) of
the other side, forcing that whole term. Those are the ``c1..ck`` and ``m1..mk``
branches, and they are why an FK-B node can have many children. When neither
side qualifies it falls back to the plain two-way split, which is FK-A's step in
all but name.

``variant``
    ``"faithful"`` is ``FK_B.m``: stop at the first conflicting assignment.
    ``"multiple"`` is ``MFK_B.m``: return every one the deciding node can see.
    The control flow, and so the tree, is identical; only the size of the answer
    changes.

``split_rule``
    ``"most_frequent"`` is ``Choose_SplitVar.m``'s ``'mostFreq'``. ``"degree_sum"``
    is FK-A's rule, so both algorithms can be run on identical split choices.
    Both break ties toward the lowest index.

Three deviations from the MATLAB, all corrections:

* an empty term is kept, not deleted -- see :func:`fka.hypergraph.phi_x_0`;
* ``A_c_x``/``A_m_x`` return ordinary term sets where the MATLAB returns the
  scalar ``false`` (and one path calls an undefined ``stop``);
* "no conflicting assignment" is ``None``, not ``[]``. The MATLAB uses ``[]``
  for both, so a genuine empty-set conflict is indistinguishable from
  equivalence.
"""

from __future__ import annotations

import math
from typing import Iterable, Optional

from fka.algorithm import eps, eps_v
from fka.hypergraph import (
    A_c_x,
    A_m_x,
    Hypergraph,
    phi_x_0,
    phi_x_1,
    popcount,
    verts,
)
from fka.irredundant import irredundant
from fka.transversal import transversal
from fka.tree import Node, Tree, Verdict

__all__ = [
    "fk_b",
    "is_dual",
    "mu",
    "mu_frequent",
    "choose_split_var",
    "check_conditions",
    "easy_case",
    "cnf_value",
    "dnf_value",
    "is_conflict",
    "VARIANTS",
    "SPLIT_RULES",
    "BRANCHES",
]

VARIANTS = ("faithful", "multiple")
SPLIT_RULES = ("most_frequent", "degree_sum")
#: The three shapes a split can take, recorded on each internal node.
BRANCHES = ("mu_D", "mu_C", "split")


# ----------------------------------------------------------------------
# evaluating the two normal forms
# ----------------------------------------------------------------------
def cnf_value(cnf: Hypergraph, S: int) -> bool:
    """``C(S)``: true when every clause meets ``S``."""
    return all(c & S for c in cnf.edges)


def dnf_value(dnf: Hypergraph, S: int) -> bool:
    """``D(S)``: true when some monomial is contained in ``S``."""
    return any((m & S) == m for m in dnf.edges)


def is_conflict(cnf: Hypergraph, dnf: Hypergraph, S: int) -> bool:
    """Whether ``C(S) != D(S)``.

    Once condition 1 holds, ``D <= C`` pointwise and the only reachable polarity
    is ``C(S) = 1, D(S) = 0`` -- which is thesis equation 2.1, so FK-A and FK-B
    agree on what a counterexample looks like. Condition 1 failing is the one
    case that yields the other polarity.
    """
    return cnf_value(cnf, S) != dnf_value(dnf, S)


# ----------------------------------------------------------------------
# split-variable selection
# ----------------------------------------------------------------------
def mu(n: int) -> float:
    """``mu(n) = log n / log log n`` (mu_function.m).

    ``log log n`` is undefined at or below ``e``. FK-B only evaluates this with
    ``nC * nD >= 9``, so the guard keeps direct callers honest rather than
    papering over a live case; infinity makes no variable mu-infrequent, which
    is the safe answer.
    """
    if n < 3:
        return math.inf
    return math.log(n) / math.log(math.log(n))


def mu_frequent(x: int, A: Hypergraph, B: Hypergraph) -> bool:
    """``|{t in A : x in t}| / |A| <= 1 / mu(|A| * |B|)`` (mu_frequent_in_A.m).

    Read as "``x`` is at most mu-frequent in ``A``": rare enough on that side
    that FK-B can afford one subproblem per term of the other side.
    """
    if not len(A):
        return False
    return A.degree(x) / len(A) <= 1.0 / mu(len(A) * len(B))


def choose_split_var(
    cnf: Hypergraph, dnf: Hypergraph, rule: str = "most_frequent"
) -> Optional[int]:
    """The split variable, or ``None`` when no variable has positive degree.

    Only the joint support is considered; ``phi_x_0`` leaves variables behind in
    no term at all, and splitting on one of those would recurse forever.
    """
    if rule not in SPLIT_RULES:
        raise ValueError(f"unknown split rule {rule!r}; expected one of {SPLIT_RULES}")

    live = verts(cnf.support() | dnf.support())
    if not live:
        return None

    if rule == "degree_sum":
        return max(live, key=lambda v: (cnf.degree(v) + dnf.degree(v), -v))

    # 'mostFreq': walk down both degree orders together and take the first
    # variable that has appeared in the top i of each.
    by_c = sorted(live, key=lambda v: (-cnf.degree(v), v))
    by_d = sorted(live, key=lambda v: (-dnf.degree(v), v))
    for i in range(1, len(live) + 1):
        common = set(by_c[:i]) & set(by_d[:i])
        if common:
            return min(common)
    return None


# ----------------------------------------------------------------------
# conditions and easy cases
# ----------------------------------------------------------------------
def _label(t: int) -> str:
    return "{" + ",".join(str(v + 1) for v in verts(t)) + "}"


def _conflict(
    reason: str,
    detail: str,
    found: Iterable[int],
    variant: str,
    witness_edges: tuple[int, ...] = (),
) -> Verdict:
    """A not-equivalent verdict carrying one or all of the assignments found."""
    seen = tuple(dict.fromkeys(found))
    return Verdict(
        dual=False,
        reason=reason,
        detail=detail,
        CA=seen[0],
        witness_edges=witness_edges,
        CAs=seen if variant == "multiple" else (),
    )


def _constant(MBF: Hypergraph, *, clauses: bool) -> Optional[int]:
    """``0`` or ``1`` when ``MBF`` is a constant function, else ``None``.

    An empty CNF is the constant 1 and a CNF holding the empty clause is the
    constant 0; a DNF reads the other way round.
    """
    if not MBF.edges:
        return 1 if clauses else 0
    if 0 in MBF.edges:
        return 0 if clauses else 1
    return None


def _one_per_clause(cnf: Hypergraph) -> int:
    """The lowest variable of each clause -- Check_Conditions.m's transversal."""
    S = 0
    for c in cnf.edges:
        S |= c & -c
    return S


def _decide_constant(
    cnf: Hypergraph,
    dnf: Hypergraph,
    c: Optional[int],
    d: Optional[int],
    V: int,
    variant: str,
) -> Verdict:
    """Settle an instance in which at least one side is a constant function."""
    if c is not None and d is not None:
        if c == d:
            return Verdict(
                dual=True, reason="constant", detail=f"both sides are the constant {c}"
            )
        # C = 1, D = 0 conflicts at the empty set; C = 0, D = 1 at all of V.
        S = 0 if c == 1 else V
    elif c == 1:
        S = 0  # D has no empty monomial, so D is false on the empty set
    elif c == 0:
        S = V  # D has a monomial and it lies inside V, so D is true on V
    elif d == 0:
        S = _one_per_clause(cnf)  # meets every clause, contains no monomial
    else:
        S = 0  # d == 1: C has a non-empty clause, so C is false on the empty set

    sides = []
    if c is not None:
        sides.append(f"C is the constant {c}")
    if d is not None:
        sides.append(f"D is the constant {d}")
    return _conflict("constant", "; ".join(sides), (S,), variant)


def check_conditions(
    cnf: Hypergraph, dnf: Hypergraph, *, variant: str = "faithful"
) -> Optional[Verdict]:
    """FK-B's three necessary conditions (Check_Conditions.m).

    1. every clause meets every monomial;
    2. both sides mention exactly the same variables;
    3. ``max |m| <= nC`` and ``max |c| <= nD``.

    ``None`` when all hold, otherwise a decided verdict. Unlike the MATLAB, that
    verdict can be ``dual=True``: a degenerate side settles the instance outright
    rather than falling through to a recursion with nothing left to split.
    ``multiple`` collects every assignment a failed condition exposes, matching
    Multiple_Check_Conditions.m. Both sides are assumed irredundant.
    """
    V = cnf.support() | dnf.support()

    c_const = _constant(cnf, clauses=True)
    d_const = _constant(dnf, clauses=False)
    if c_const is not None or d_const is not None:
        return _decide_constant(cnf, dnf, c_const, d_const, V, variant)

    # 1: a clause disjoint from a monomial. Setting that monomial true makes D
    #    true while the clause it misses keeps C false -- the one place FK-B
    #    reports a conflict of the opposite polarity to equation 2.1.
    missed = [m for m in dnf.edges if any(not (c & m) for c in cnf.edges)]
    if missed:
        m = missed[0]
        c = next(c for c in cnf.edges if not (c & m))
        return _conflict(
            "cond_1_disjoint",
            f"condition 1: clause {_label(c)} and monomial {_label(m)} are disjoint",
            missed,
            variant,
            witness_edges=(c, m),
        )

    # 2: unequal supports. Whichever side has the extra variables, the term
    #    holding one of them yields an assignment directly.
    if cnf.support() != dnf.support():
        extra_d = dnf.support() & ~cnf.support()
        if extra_d:
            # m without the variables C never mentions still meets every clause
            # (they cannot have used them), and is a proper subset of m, so it
            # contains no monomial.
            found = [m & ~extra_d for m in dnf.edges if m & extra_d]
            extra, where = extra_d, "D"
        else:
            # Everything C mentions except c, plus the extras. Every other
            # clause sticks out of c, c itself holds an extra, and every
            # monomial keeps a variable of c that is excluded.
            extra_c = cnf.support() & ~dnf.support()
            found = [(cnf.support() & ~c) | extra_c for c in cnf.edges if c & extra_c]
            extra, where = extra_c, "C"
        return _conflict(
            "cond_2_support",
            f"condition 2: {_label(extra)} appears only in {where}",
            found,
            variant,
        )

    # 3: a term longer than the other side's term count. Picking one variable of
    #    it per term of the other side lands inside it, so the bound is what
    #    makes the choice a *proper* subset -- which is what forces the conflict.
    if max(dnf.edge_sizes()) > len(cnf):
        m = max(dnf.edges, key=popcount)
        S = 0
        for c in cnf.edges:
            hit = c & m
            S |= hit & -hit
        return _conflict(
            "cond_3_size",
            f"condition 3: monomial {_label(m)} has {popcount(m)} variables "
            f"but nC = {len(cnf)}",
            (S,),
            variant,
            witness_edges=(m,),
        )
    if max(cnf.edge_sizes()) > len(dnf):
        c = max(cnf.edges, key=popcount)
        T = 0
        for m in dnf.edges:
            hit = m & c
            T |= hit & -hit
        return _conflict(
            "cond_3_size",
            f"condition 3: clause {_label(c)} has {popcount(c)} variables "
            f"but nD = {len(dnf)}",
            (V & ~T,),
            variant,
            witness_edges=(c,),
        )
    return None


def easy_case(
    cnf: Hypergraph, dnf: Hypergraph, *, variant: str = "faithful"
) -> Verdict:
    """Decide ``min(nC, nD) <= 2`` outright (Easy_case.m).

    Reached only after :func:`check_conditions` passes, so ``D <= C`` pointwise
    and the only possible conflict has ``C(S) = 1, D(S) = 0``. One exists iff the
    two sides differ, and it can be read off the small side:

    * ``nD <= 2``: for each ``t`` in ``Tr(D)`` try ``S = V \\ t``. It contains no
      monomial because ``t`` meets them all, so only ``C(S)`` needs testing.
    * ``nC <= 2``: for each ``s`` in ``Tr(C)`` try ``S = s``. It meets every
      clause by construction, so only ``D(S)`` needs testing.

    Both are complete: any ``S`` with ``D(S) = 0`` sits inside some ``V \\ t``,
    any ``S`` with ``C(S) = 1`` contains some ``s``, and monotonicity carries the
    conflict to the extreme candidate. A two-term side has at most ``|t1|*|t2|``
    minimal transversals, so this is the cheap decision the name promises rather
    than the ``2^n`` ``dec2bin`` search CA_CNF2.m runs.
    """
    V = cnf.support() | dnf.support()
    found: list[int] = []
    if len(dnf) <= 2:
        for t in transversal(dnf).edges:
            S = V & ~t
            if cnf_value(cnf, S):
                found.append(S)
                if variant != "multiple":
                    break
    else:
        for s in transversal(cnf).edges:
            if not dnf_value(dnf, s):
                found.append(s)
                if variant != "multiple":
                    break

    small = min(len(cnf), len(dnf))
    if not found:
        return Verdict(
            dual=True,
            reason="easy_case",
            detail=f"min(nC, nD) = {small} <= 2 and the two sides agree",
        )
    return _conflict(
        "easy_case",
        f"min(nC, nD) = {small} <= 2 and the two sides differ",
        found,
        variant,
    )


# ----------------------------------------------------------------------
# the algorithm
# ----------------------------------------------------------------------
def _set_false(cnf: Hypergraph, dnf: Hypergraph, s: int) -> tuple[Hypergraph, Hypergraph]:
    """Both sides with the variables of ``s`` forced FALSE (A_c_x.m)."""
    return A_c_x(cnf, s, clauses=True), A_c_x(dnf, s, clauses=False)


def _set_true(cnf: Hypergraph, dnf: Hypergraph, s: int) -> tuple[Hypergraph, Hypergraph]:
    """Both sides with the variables of ``s`` forced TRUE (A_m_x.m)."""
    return A_m_x(cnf, s, clauses=True), A_m_x(dnf, s, clauses=False)


def _lift(child: Verdict, add: int, label: str) -> Verdict:
    """Carry a subproblem's assignment up, restoring what it forced.

    ``add`` is the set the parent held TRUE while the subproblem ran: the split
    variable for ``x=1`` and for the per-clause branch, the whole monomial for
    the per-monomial branch, nothing for ``x=0``. Variables the subproblem never
    saw are absent from its assignment, so the union is well defined.
    """
    return Verdict(
        dual=False,
        reason="child_failed",
        detail=f"subproblem {label} returned a conflicting assignment",
        CA=None if child.CA is None else child.CA | add,
        witness_edges=child.witness_edges,
        CAs=tuple(S | add for S in child.CAs),
    )


def fk_b(
    cnf: Hypergraph,
    dnf: Hypergraph,
    *,
    variant: str = "faithful",
    split_rule: str = "most_frequent",
    redundancy_check: bool = True,
    instance: str = "",
    validate: bool = False,
    max_nodes: Optional[int] = None,
) -> Tree:
    """Decide whether the CNF ``cnf`` and the DNF ``dnf`` are the same function.

    Equivalently -- and this is how the result is reported, so FK-A and FK-B are
    directly comparable -- whether ``cnf = Tr(dnf)``. ``tree.dual`` is ``True``
    when no conflicting assignment exists; otherwise ``tree.verdict.CA`` is one.

    ``redundancy_check=False`` skips the per-node ``Irredundant`` that FK-B's
    conditions    assume; it exists to measure what that reduction costs and buys, not as a
    supported way to get an answer. ``validate=True`` re-checks every assignment
    against ``C(S) != D(S)`` at the node that produced it. ``max_nodes`` aborts
    once the tree passes that size.
    """
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant {variant!r}; expected one of {VARIANTS}")
    if split_rule not in SPLIT_RULES:
        raise ValueError(
            f"unknown split rule {split_rule!r}; expected one of {SPLIT_RULES}"
        )
    if cnf.n != dnf.n:
        raise ValueError(
            f"C is on {cnf.n} variables but D is on {dnf.n}; "
            "both sides must share a ground set"
        )

    tree = Tree(
        instance=instance, algorithm="FK-B", variant=variant, split_rule=split_rule
    )
    counter = {"next_id": 1}

    def recurse(
        cnf_in: Hypergraph,
        dnf_in: Hypergraph,
        parent_id: Optional[int],
        path: tuple[str, ...],
    ) -> Verdict:
        if max_nodes is not None and len(tree) >= max_nodes:
            raise RuntimeError(
                f"recursion tree exceeded max_nodes={max_nodes} on instance "
                f"{instance!r}; raise the limit or use a smaller instance"
            )

        node_id = counter["next_id"]
        counter["next_id"] += 1

        if redundancy_check:
            rc, rd = irredundant(cnf_in), irredundant(dnf_in)
            C, D, cut_C, cut_D = rc.reduced, rd.reduced, rc.removed, rd.removed
        else:
            C, D, cut_C, cut_D = cnf_in, dnf_in, (), ()

        node = Node(
            node_id=node_id,
            parent_id=parent_id,
            path=path,
            cnf_in=cnf_in,
            dnf_in=dnf_in,
            cnf=C,
            dnf=D,
            cut_C=cut_C,
            cut_D=cut_D,
            epsilon=eps(C, D),
        )
        tree.add(node)

        def settle(verdict: Verdict) -> Verdict:
            if validate and not verdict.dual:
                seen = verdict.CAs or (() if verdict.CA is None else (verdict.CA,))
                for S in seen:
                    if not is_conflict(C, D, S):
                        raise AssertionError(
                            f"node {node_id}: {_label(S)} is not a conflicting "
                            f"assignment for {C.label()} / {D.label()}"
                        )
            node.verdict = verdict
            return verdict

        decided = check_conditions(C, D, variant=variant)
        if decided is not None:
            return settle(decided)

        if min(len(C), len(D)) <= 2:
            return settle(easy_case(C, D, variant=variant))

        x = choose_split_var(C, D, split_rule)
        if x is None:
            # Unreachable: condition 2 passed, so both supports are equal and
            # non-empty, and every variable in them has positive degree.
            raise AssertionError(
                f"node {node_id}: no variable of positive degree in "
                f"{C.label()} / {D.label()}"
            )
        node.split_var = x
        node.split_var_frequency = eps_v(C, D, x)
        bit = 1 << x

        C_x_0, C_x_1 = phi_x_0(C, x), phi_x_1(C, x)
        D_x_0, D_x_1 = phi_x_0(D, x), phi_x_1(D, x)

        if mu_frequent(x, D, C):
            node.branch = "mu_D"
            # line 8: the x = 0 restriction.
            v = recurse(*_set_false(C, D, bit), node_id, path + ("x=0",))
            if not v.dual:
                return settle(_lift(v, 0, "x=0"))
            # lines 10-12: one subproblem per clause of C_x_0, that clause forced
            # FALSE. D_x_1 can be dropped from them because the x = 0 branch just
            # proved it equal to a CNF this clause falsifies.
            for i, c in enumerate(C_x_0.edges, start=1):
                v = recurse(*_set_false(C_x_1, D_x_0, c), node_id, path + (f"c{i}",))
                if not v.dual:
                    return settle(_lift(v, bit, f"c{i}"))

        elif mu_frequent(x, C, D):
            node.branch = "mu_C"
            # line 14: the x = 1 restriction.
            v = recurse(*_set_true(C, D, bit), node_id, path + ("x=1",))
            if not v.dual:
                return settle(_lift(v, bit, "x=1"))
            # lines 16-17: the mirror image, one subproblem per monomial of
            # D_x_0 with that monomial forced TRUE.
            for i, m in enumerate(D_x_0.edges, start=1):
                v = recurse(*_set_true(C_x_0, D_x_1, m), node_id, path + (f"m{i}",))
                if not v.dual:
                    return settle(_lift(v, m, f"m{i}"))

        else:
            node.branch = "split"
            # lines 19-23: neither side is sparse enough, so split both ways.
            v = recurse(*_set_false(C, D, bit), node_id, path + ("x=0",))
            if not v.dual:
                return settle(_lift(v, 0, "x=0"))
            v = recurse(*_set_true(C, D, bit), node_id, path + ("x=1",))
            if not v.dual:
                return settle(_lift(v, bit, "x=1"))

        return settle(
            Verdict(
                dual=True,
                reason="equivalent",
                detail=f"every subproblem of the split on v{x + 1} agreed",
            )
        )

    verdict = recurse(cnf, dnf, None, ())
    tree.dual = verdict.dual
    tree.verdict = verdict
    return tree


def is_dual(
    cnf: Hypergraph,
    dnf: Hypergraph,
    *,
    variant: str = "faithful",
    split_rule: str = "most_frequent",
) -> bool:
    """``True`` iff ``cnf = Tr(dnf)``. Convenience wrapper that discards the tree."""
    return bool(fk_b(cnf, dnf, variant=variant, split_rule=split_rule).dual)
