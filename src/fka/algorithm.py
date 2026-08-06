"""The Fredman-Khachiyan algorithm FK-A.

Thesis Algorithm 1 (p.13), step for step:

1. ``Irredundant`` both sides;
2. check conditions (i)-(iv); any failure means the pair is not dual;
3. resolve trivial instances;
4. otherwise pick a frequent variable ``x`` and recurse on
   ``L: (C_x_1, D_x_0 v D_x_1)`` and ``R: (D_x_1, C_x_0 v C_x_1)``, TRUE only if
   both are.

Two axes are parameterised so the variants can be compared on one instance:

``variant``
    ``"faithful"`` stops at Step 3's ``nC * nD <= 1``, as FK96 states it.
    ``"modified"`` adds the ``D_{1,k}`` check of §5.1.1 (p.47): with ``nC = 1``
    on a term of size ``k``, duality is decided in ``O(k)`` by testing whether
    the other side is precisely that term's ``k`` singletons -- saving the
    ``k-1`` further recursions the faithful version needs.

``split_rule``
    ``"degree_sum"`` takes ``argmax_x (deg_C(x) + deg_D(x))``; it produced the
    published trees, so it is the default. ``"frequency"`` takes
    ``argmax_x eps_v``, the quantity Lemma 2.3.2 (p.16) bounds below by
    ``1/log n``. Both break ties to the lowest index.

Neither switch affects the answer -- every frequent variable gives a valid
decomposition. They change the shape and size of the tree, which is the object
under study.
"""

from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
from typing import Optional

from .hypergraph import Hypergraph, phi_x_0, phi_x_1, popcount, verts
from .irredundant import irredundant
from .transversal import find_CA
from .tree import Node, Tree, Verdict

__all__ = [
    "fk_a",
    "is_dual",
    "eps",
    "eps_v",
    "eps_exact",
    "choose_split_var",
    "check_conditions",
    "verify_CA",
    "VARIANTS",
    "SPLIT_RULES",
]

VARIANTS = ("faithful", "modified")
SPLIT_RULES = ("degree_sum", "frequency")


# ----------------------------------------------------------------------
# frequency and split-variable selection
# ----------------------------------------------------------------------
def eps_v(cnf: Hypergraph, dnf: Hypergraph, x: int) -> float:
    """``eps(x) = max(deg_C(x)/nC, deg_D(x)/nD)`` (thesis §2.2.2, p.14).

    An empty side contributes 0 rather than dividing by zero.
    """
    fc = cnf.degree(x) / len(cnf) if len(cnf) else 0.0
    fd = dnf.degree(x) / len(dnf) if len(dnf) else 0.0
    return max(fc, fd)


def eps(cnf: Hypergraph, dnf: Hypergraph) -> float:
    """``eps(C, D)``: the frequency of the most frequent variable."""
    n = max(cnf.n, dnf.n)
    if n == 0:
        return 0.0
    return max(eps_v(cnf, dnf, x) for x in range(n))


def eps_exact(cnf: Hypergraph, dnf: Hypergraph) -> Fraction:
    """:func:`eps` as an exact fraction.

    The thesis quotes values like ``5/11`` and ``3/7`` (§5.2), so the instance
    baselines and the Chapter 4 bound store the fraction; the algorithm itself
    only ever compares, so it uses the float.
    """
    best = Fraction(0)
    for x in range(max(cnf.n, dnf.n)):
        if len(cnf):
            best = max(best, Fraction(cnf.degree(x), len(cnf)))
        if len(dnf):
            best = max(best, Fraction(dnf.degree(x), len(dnf)))
    return best


def choose_split_var(
    cnf: Hypergraph, dnf: Hypergraph, rule: str = "degree_sum"
) -> Optional[int]:
    """The split variable, or ``None`` if every variable has degree 0.

    Ties go to the lowest index, which is what reproduces the published trees.
    """
    if rule not in SPLIT_RULES:
        raise ValueError(f"unknown split rule {rule!r}; expected one of {SPLIT_RULES}")
    n = max(cnf.n, dnf.n)
    best: Optional[int] = None
    best_key = -1.0
    for x in range(n):
        dc, dd = cnf.degree(x), dnf.degree(x)
        if dc + dd == 0:
            continue
        key = float(dc + dd) if rule == "degree_sum" else eps_v(cnf, dnf, x)
        if key > best_key:
            best_key, best = key, x
    return best


# ----------------------------------------------------------------------
# conditions
# ----------------------------------------------------------------------
def check_conditions(cnf: Hypergraph, dnf: Hypergraph) -> Optional[Verdict]:
    """Algorithm 1 Step 2, conditions (i)-(iv).

    ``None`` when all hold, else a not-dual verdict naming the one that failed.
    Both sides are assumed irredundant. FK-B numbers its own conditions 1-3 in a
    different order and drops (iv); see :func:`fkb.algorithm.check_conditions`.
    """
    # (i) equal supports (Prop. 2.2.2).
    if cnf.support() != dnf.support():
        diff = cnf.support() ^ dnf.support()
        return Verdict(
            dual=False,
            reason="cond_i_support",
            detail=(
                "condition (i): the ground sets differ; variables "
                f"{[v + 1 for v in verts(diff)]} appear on one side only"
            ),
        )

    # (ii) each side's largest term is bounded by the other's term count
    #      (Prop. 2.2.3).
    if cnf.edges and len(dnf) < max(cnf.edge_sizes()):
        big = max(cnf.edges, key=popcount)
        return Verdict(
            dual=False,
            reason="cond_ii_size",
            detail=(
                f"condition (ii): C has a term of size {popcount(big)} "
                f"but nD = {len(dnf)}"
            ),
            witness_edges=(big,),
        )
    if dnf.edges and len(cnf) < max(dnf.edge_sizes()):
        big = max(dnf.edges, key=popcount)
        return Verdict(
            dual=False,
            reason="cond_ii_size",
            detail=(
                f"condition (ii): D has a term of size {popcount(big)} "
                f"but nC = {len(cnf)}"
            ),
            witness_edges=(big,),
        )

    # (iii) every cross pair intersects (Prop. 2.2.1) -- the witness Algorithm 1
    #       returns at line 6.
    for a in cnf.edges:
        for b in dnf.edges:
            if not (a & b):
                return Verdict(
                    dual=False,
                    reason="cond_iii_disjoint",
                    detail=(
                        f"condition (iii): clause {[v + 1 for v in verts(a)]} and "
                        f"monomial {[v + 1 for v in verts(b)]} are disjoint"
                    ),
                    witness_edges=(a, b),
                )

    # (iv) the sum of 2^-|t| over both sides is at least 1 (Prop. 2.2.4), in
    #      exact integer arithmetic. Floats here can land just under 1.0 through
    #      rounding and reject a genuine dual pair; scaling by 2^L for the
    #      largest term makes every summand an integer.
    all_terms = cnf.edges + dnf.edges
    scale = max((popcount(t) for t in all_terms), default=0)
    if sum(1 << (scale - popcount(t)) for t in all_terms) < (1 << scale):
        return Verdict(
            dual=False,
            reason="cond_iv_frequency",
            detail="condition (iv): sum of 2^-|t| over both sides is below 1",
        )
    return None


# ----------------------------------------------------------------------
# conflicting assignments
# ----------------------------------------------------------------------
def verify_CA(cnf: Hypergraph, dnf: Hypergraph, S: int) -> bool:
    """Check thesis equation 2.1 (p.11): ``S`` meets every clause and contains
    no monomial, which certifies that the pair is not dual."""
    if any(not (S & t) for t in cnf.edges):
        return False
    if any((t & S) == t for t in dnf.edges):
        return False
    return True


def _trivial_CA(cnf: Hypergraph, dnf: Hypergraph) -> Optional[int]:
    """A conflicting assignment when ``nC = 1`` or ``nD = 1``, else ``None``.

    ``nC = 1``: the single clause ``c`` is dual to its ``|c|`` singletons. If
    not, some ``u`` in ``c`` has ``{u}`` not in ``D``, and ``S = {u}`` meets
    ``c`` while containing no monomial -- every monomial meets ``c`` by
    condition (iii), so none is empty and none equals ``{u}``.

    ``nD = 1`` mirrors it: some ``w`` in the single monomial has ``{w}`` not in
    ``C``, and ``S = V \\ {w}`` meets every clause while excluding it.
    """
    if len(cnf) == 1:
        for u in verts(cnf.edges[0]):
            if verify_CA(cnf, dnf, 1 << u):
                return 1 << u
    if len(dnf) == 1:
        full = (1 << max(cnf.n, dnf.n)) - 1
        for w in verts(dnf.edges[0]):
            if verify_CA(cnf, dnf, full & ~(1 << w)):
                return full & ~(1 << w)
    if len(cnf) == 0 and verify_CA(cnf, dnf, 0):
        # C is the constant 1; the empty set vacuously meets its every clause.
        return 0
    if len(dnf) == 0:
        full = (1 << max(cnf.n, dnf.n)) - 1
        if verify_CA(cnf, dnf, full):
            return full
    return None


# ----------------------------------------------------------------------
# trivial cases
# ----------------------------------------------------------------------
def _is_singletons(MBF: Hypergraph, t: int) -> bool:
    """True iff ``MBF`` is exactly ``{{v} : v in t}``, i.e. ``MBF = Tr({t})``."""
    want = {1 << v for v in verts(t)}
    return len(MBF) == len(want) and set(MBF.edges) == want


def _trivial_case(
    cnf: Hypergraph, dnf: Hypergraph, variant: str
) -> Optional[Verdict]:
    """Step 3, and under ``variant="modified"`` Step 3'.

    The degenerate cases follow the monotone-Boolean reading of FK96: no terms
    is the constant 0 read as a DNF, ``{{}}`` is the constant 1, and the two are
    dual. See ``docs/dictionary.md``.
    """
    if len(cnf) * len(dnf) <= 1:
        if len(cnf) == 0:
            dual = len(dnf) == 1 and dnf.edges[0] == 0
        elif len(dnf) == 0:
            dual = len(cnf) == 1 and cnf.edges[0] == 0
        else:  # nC = nD = 1
            c, m = cnf.edges[0], dnf.edges[0]
            dual = popcount(c) == 1 and c == m
        return Verdict(
            dual=dual,
            reason="trivial",
            detail=f"nC * nD = {len(cnf) * len(dnf)} <= 1",
            CA=None if dual else _trivial_CA(cnf, dnf),
        )

    if variant != "modified":
        return None

    # D_{1,k} / D_{k,1}: thesis §5.1.1 (p.47).
    for one, other, name, other_name in ((cnf, dnf, "C", "D"), (dnf, cnf, "D", "C")):
        if len(one) == 1:
            t = one.edges[0]
            dual = _is_singletons(other, t)
            return Verdict(
                dual=dual,
                reason="single_edge",
                detail=(
                    f"n{name} = 1 with a term of size {popcount(t)}; "
                    f"{other_name} {'is' if dual else 'is not'} its "
                    f"{popcount(t)} singletons"
                ),
                CA=None if dual else _trivial_CA(cnf, dnf),
            )
    return None


# ----------------------------------------------------------------------
# the algorithm
# ----------------------------------------------------------------------
def fk_a(
    cnf: Hypergraph,
    dnf: Hypergraph,
    *,
    variant: str = "modified",
    split_rule: str = "degree_sum",
    instance: str = "",
    validate: bool = False,
    max_nodes: Optional[int] = None,
) -> Tree:
    """Decide whether ``cnf = Tr(dnf)``, returning the full recursion tree.

    ``validate=True`` re-checks every conflicting assignment against equation
    2.1 at the node that built it; the tests turn it on. ``max_nodes`` aborts
    once the tree passes that size -- the tree is the thing that can blow up.
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
        instance=instance, algorithm="FK-A", variant=variant, split_rule=split_rule
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

        # Step 1: Irredundant, recording what it removed.
        rc, rd = irredundant(cnf_in), irredundant(dnf_in)
        C, D = rc.reduced, rd.reduced

        node = Node(
            node_id=node_id,
            parent_id=parent_id,
            path=path,
            cnf_in=cnf_in,
            dnf_in=dnf_in,
            cnf=C,
            dnf=D,
            cut_C=rc.removed,
            cut_D=rd.removed,
            epsilon=eps(C, D),
        )
        tree.add(node)

        def settle(verdict: Verdict) -> Verdict:
            if validate and verdict.CA is not None:
                if not verify_CA(C, D, verdict.CA):
                    raise AssertionError(
                        f"node {node_id}: {[v + 1 for v in verts(verdict.CA)]} "
                        f"fails equation 2.1 for {C.label()} / {D.label()}"
                    )
            node.verdict = verdict
            return verdict

        # Step 2: conditions.
        failed = check_conditions(C, D)
        if failed is not None:
            # Algorithm 1 reports a structural witness here. Equation 2.1's set
            # assignment is more useful in a research artifact -- it can be
            # checked independently at every ancestor after the lifts -- but not
            # every degenerate pair admits one in this orientation, so the
            # structural verdict stands when the exact search returns None.
            if failed.CA is None:
                found = find_CA(C, D)
                if found is not None:
                    failed = replace(failed, CA=found)
            return settle(failed)

        # Step 3: trivial instances, and the modified single-term rule.
        trivial = _trivial_case(C, D, variant)
        if trivial is not None:
            return settle(trivial)

        # Step 4: recursive decomposition.
        x = choose_split_var(C, D, split_rule)
        if x is None:
            # Unreachable for irredundant inputs past Step 2: every term is
            # non-empty there, so some variable has positive degree.
            raise AssertionError(
                f"node {node_id}: no variable of positive degree in "
                f"{C.label()} / {D.label()}"
            )
        node.split_var = x
        node.split_var_frequency = eps_v(C, D, x)

        C_x_0, C_x_1 = phi_x_0(C, x), phi_x_1(C, x)
        D_x_0, D_x_1 = phi_x_0(D, x), phi_x_1(D, x)

        left = recurse(C_x_1, D_x_0.vee(D_x_1), node_id, path + ("L",))
        if not left.dual:
            # The child's S misses x, so S + {x} still meets every clause --
            # those through x are met by x -- and still contains no monomial.
            return settle(
                Verdict(
                    dual=False,
                    reason="child_failed",
                    detail=f"left subproblem (C_x_1, D_x_0 v D_x_1) on v{x + 1} failed",
                    CA=None if left.CA is None else left.CA | (1 << x),
                    witness_edges=left.witness_edges,
                )
            )

        right = recurse(D_x_1, C_x_0.vee(C_x_1), node_id, path + ("R",))
        if not right.dual:
            # The right subproblem has the sides swapped, so its assignment is
            # complemented within V \ {x} (thesis p.11: S works for (C,D) iff
            # V \ S works for (D,C)).
            CA = None
            if right.CA is not None:
                CA = (((1 << C.n) - 1) & ~(1 << x)) & ~right.CA
            return settle(
                Verdict(
                    dual=False,
                    reason="child_failed",
                    detail=f"right subproblem (D_x_1, C_x_0 v C_x_1) on v{x + 1} failed",
                    CA=CA,
                    witness_edges=right.witness_edges,
                )
            )

        return settle(
            Verdict(dual=True, reason="dual", detail="both subproblems returned TRUE")
        )

    verdict = recurse(cnf, dnf, None, ())
    tree.dual = verdict.dual
    tree.verdict = verdict
    return tree


def is_dual(
    cnf: Hypergraph,
    dnf: Hypergraph,
    *,
    variant: str = "modified",
    split_rule: str = "degree_sum",
) -> bool:
    """``True`` iff ``cnf = Tr(dnf)``. Discards the tree."""
    return bool(fk_a(cnf, dnf, variant=variant, split_rule=split_rule).dual)
