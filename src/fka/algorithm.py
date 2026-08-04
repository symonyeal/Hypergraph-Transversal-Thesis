"""The Fredman-Khachiyan algorithm FK-A.

Thesis Algorithm 1 (p.13), step for step:

1. Sperner-reduce ``G`` and ``H``.
2. Check conditions (i)-(iv); any failure means ``G != Tr(H)``.
3. Resolve trivial instances.
4. Otherwise pick a frequent vertex ``v`` and recurse on
   ``L: (G_1, H_0 or H_1)`` and ``R: (H_1, G_0 or G_1)``, TRUE only if both are.

Two axes are parameterised rather than hard-coded, so the variants can be
compared on one instance:

``variant``
    ``"faithful"`` stops at Step 3's ``|G|*|H| <= 1``, exactly as FK96 states
    it. ``"modified"`` adds the ``D_{1,k}`` check of §5.1.1 (p.47): with
    ``|G| = 1`` on an edge of size ``k``, duality is decided in ``O(k)`` by
    testing whether ``H`` is precisely that edge's ``k`` singletons -- saving
    the ``k-1`` further recursions the faithful version needs.

``pivot_rule``
    ``"degree_sum"`` takes ``argmax_v (deg_G(v) + deg_H(v))``. That is what the
    legacy implementation did, so it is what produced the published trees, so it
    is the default. ``"frequency"`` takes ``argmax_v eps_v``, the quantity
    Lemma 2.3.2 (p.16) actually bounds below by ``1/log n``. Both break ties to
    the lowest index.

Neither switch affects the answer -- every frequent vertex gives a valid
decomposition. They change the shape and size of the tree, which is the object
under study.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Optional

from .hypergraph import Hypergraph, popcount, verts
from .sperner import sperner_reduce
from .transversal import find_certificate
from .tree import RecursionNode, RecursionTree, Verdict

__all__ = [
    "fk_a",
    "is_dual",
    "eps",
    "eps_v",
    "choose_pivot",
    "check_conditions",
    "verify_certificate",
    "VARIANTS",
    "PIVOT_RULES",
]

VARIANTS = ("faithful", "modified")
PIVOT_RULES = ("degree_sum", "frequency")


# ----------------------------------------------------------------------
# frequency and pivot selection
# ----------------------------------------------------------------------
def eps_v(G: Hypergraph, H: Hypergraph, v: int) -> float:
    """``eps(v) = max(deg_G(v)/|G|, deg_H(v)/|H|)`` (thesis §2.2.2, p.14).

    An empty side contributes 0 rather than dividing by zero.
    """
    fg = G.degree(v) / len(G) if len(G) else 0.0
    fh = H.degree(v) / len(H) if len(H) else 0.0
    return max(fg, fh)


def eps(G: Hypergraph, H: Hypergraph) -> float:
    """``eps(G, H)``: the frequency of the most frequent vertex.

    What the thesis reports per instance -- 3/7 for Fano, 1/2 for 6-A..6-D,
    5/11 for 7-Ver, 2/5 for 8-Ver (§5.2).
    """
    n = max(G.n, H.n)
    if n == 0:
        return 0.0
    return max(eps_v(G, H, v) for v in range(n))


def choose_pivot(
    G: Hypergraph, H: Hypergraph, rule: str = "degree_sum"
) -> Optional[int]:
    """The split vertex, or ``None`` if every vertex has degree 0.

    Ties go to the lowest index, matching the legacy ``np.argmax`` so the
    published trees reproduce exactly.
    """
    n = max(G.n, H.n)
    best: Optional[int] = None
    best_key = -1.0
    for v in range(n):
        dg, dh = G.degree(v), H.degree(v)
        if dg + dh == 0:
            continue
        if rule == "degree_sum":
            key = float(dg + dh)
        elif rule == "frequency":
            key = eps_v(G, H, v)
        else:
            raise ValueError(
                f"unknown pivot rule {rule!r}; expected one of {PIVOT_RULES}"
            )
        if key > best_key:
            best_key, best = key, v
    return best


# ----------------------------------------------------------------------
# conditions
# ----------------------------------------------------------------------
def check_conditions(G: Hypergraph, H: Hypergraph) -> Optional[Verdict]:
    """Algorithm 1 Step 2, conditions (i)-(iv).

    ``None`` when all hold, else a not-dual :class:`Verdict` naming the one that
    failed. Both sides are assumed Sperner-reduced.

    FK-B numbers its own conditions 1-3 in a different order and drops (iv); see
    :func:`fkb.algorithm.check_conditions`.
    """
    # (i) V(G) = V(H): the two sides must span the same vertices (Prop. 2.2.2).
    if G.support() != H.support():
        diff = G.support() ^ H.support()
        return Verdict(
            dual=False,
            reason="cond_i_support",
            detail=(
                "condition (i): the ground sets differ; vertices "
                f"{[v + 1 for v in verts(diff)]} appear on one side only"
            ),
        )

    # (ii) each side's largest edge is bounded by the other's edge count
    #      (Prop. 2.2.3).
    if G.edges and len(H) < max(G.edge_sizes()):
        big = max(G.edges, key=popcount)
        return Verdict(
            dual=False,
            reason="cond_ii_size",
            detail=(
                f"condition (ii): G has an edge of size {popcount(big)} "
                f"but |H| = {len(H)}"
            ),
            witness_edges=(big,),
        )
    if H.edges and len(G) < max(H.edge_sizes()):
        big = max(H.edges, key=popcount)
        return Verdict(
            dual=False,
            reason="cond_ii_size",
            detail=(
                f"condition (ii): H has an edge of size {popcount(big)} "
                f"but |G| = {len(G)}"
            ),
            witness_edges=(big,),
        )

    # (iii) every cross pair intersects (Prop. 2.2.1). This is the witness the
    #       thesis returns at Algorithm 1 line 6.
    for a in G.edges:
        for b in H.edges:
            if not (a & b):
                return Verdict(
                    dual=False,
                    reason="cond_iii_disjoint",
                    detail=(
                        f"condition (iii): edges {[v + 1 for v in verts(a)]} in G "
                        f"and {[v + 1 for v in verts(b)]} in H are disjoint"
                    ),
                    witness_edges=(a, b),
                )

    # (iv) the sum of 2^-|e| over both sides is at least 1 (Prop. 2.2.4), in
    #      exact integer arithmetic. Comparing floats here can land just under
    #      1.0 through rounding and reject a genuine transversal pair; scaling
    #      by 2^L for the largest edge L makes every term an integer.
    all_edges = G.edges + H.edges
    scale = max((popcount(e) for e in all_edges), default=0)
    if sum(1 << (scale - popcount(e)) for e in all_edges) < (1 << scale):
        return Verdict(
            dual=False,
            reason="cond_iv_frequency",
            detail="condition (iv): sum of 2^-|e| over both sides is below 1",
        )
    return None


# ----------------------------------------------------------------------
# certificates
# ----------------------------------------------------------------------
def verify_certificate(G: Hypergraph, H: Hypergraph, S: int) -> bool:
    """Check thesis equation 2.1 (p.11): ``S`` meets every edge of ``G`` and
    contains none of ``H``, which certifies ``G != Tr(H)``."""
    if any(not (S & e) for e in G.edges):
        return False
    if any((e & S) == e for e in H.edges):
        return False
    return True


def _trivial_certificate(G: Hypergraph, H: Hypergraph) -> Optional[int]:
    """A certificate for a non-dual instance with ``|G| = 1`` or ``|H| = 1``.

    ``|G| = 1``: the single edge ``e`` is dual to its ``|e|`` singletons. If not,
    some ``u`` in ``e`` has ``{u}`` not in ``H``, and ``S = {u}`` meets ``e``
    while containing no edge of ``H`` -- every edge of ``H`` meets ``e`` by
    condition (iii), so none is empty and none equals ``{u}``.

    ``|H| = 1`` mirrors it: some ``w`` in the single edge has ``{w}`` not in
    ``G``, and ``S = V \\ {w}`` meets every edge of ``G`` while excluding it.

    ``None`` when no such construction applies, rather than guessing.
    """
    if len(G) == 1:
        for u in verts(G.edges[0]):
            if verify_certificate(G, H, 1 << u):
                return 1 << u
    if len(H) == 1:
        full = (1 << max(G.n, H.n)) - 1
        for w in verts(H.edges[0]):
            if verify_certificate(G, H, full & ~(1 << w)):
                return full & ~(1 << w)
    if len(G) == 0 and verify_certificate(G, H, 0):
        # G is the constant 0; the empty set vacuously meets its every edge.
        return 0
    if len(H) == 0:
        full = (1 << max(G.n, H.n)) - 1
        if verify_certificate(G, H, full):
            return full
    return None


# ----------------------------------------------------------------------
# trivial cases
# ----------------------------------------------------------------------
def _is_singletons(H: Hypergraph, e: int) -> bool:
    """True iff ``H`` is exactly ``{{v} : v in e}``, i.e. ``H = Tr({e})``."""
    want = {1 << v for v in verts(e)}
    return len(H) == len(want) and set(H.edges) == want


def _trivial_case(G: Hypergraph, H: Hypergraph, variant: str) -> Optional[Verdict]:
    """Step 3, and under ``variant="modified"`` Step 3'.

    A :class:`Verdict` when the instance is decided here, else ``None``.

    The degenerate hypergraphs follow the monotone-Boolean reading of FK96: the
    edgeless hypergraph is the constant 0, ``{{}}`` is the constant 1, and the
    two are dual.
    """
    if len(G) * len(H) <= 1:
        if len(G) == 0:
            dual = len(H) == 1 and H.edges[0] == 0
        elif len(H) == 0:
            dual = len(G) == 1 and G.edges[0] == 0
        else:  # |G| = |H| = 1
            e, f = G.edges[0], H.edges[0]
            dual = popcount(e) == 1 and e == f
        return Verdict(
            dual=dual,
            reason="trivial",
            detail=f"|G|*|H| = {len(G) * len(H)} <= 1",
            certificate=None if dual else _trivial_certificate(G, H),
        )

    if variant != "modified":
        return None

    # D_{1,k} / D_{k,1}: thesis §5.1.1 (p.47).
    for one, other, name, other_name in ((G, H, "G", "H"), (H, G, "H", "G")):
        if len(one) == 1:
            e = one.edges[0]
            dual = _is_singletons(other, e)
            return Verdict(
                dual=dual,
                reason="single_edge",
                detail=(
                    f"|{name}| = 1 with an edge of size {popcount(e)}; "
                    f"{other_name} {'is' if dual else 'is not'} its "
                    f"{popcount(e)} singletons"
                ),
                certificate=None if dual else _trivial_certificate(G, H),
            )
    return None


# ----------------------------------------------------------------------
# the algorithm
# ----------------------------------------------------------------------
def fk_a(
    G: Hypergraph,
    H: Hypergraph,
    *,
    variant: str = "modified",
    pivot_rule: str = "degree_sum",
    instance: str = "",
    validate: bool = False,
    max_nodes: Optional[int] = None,
) -> RecursionTree:
    """Decide whether ``G = Tr(H)``, returning the full recursion tree.

    ``validate=True`` re-checks every certificate against equation 2.1 at the
    node that built it. Off by default because it is only useful when changing
    the algorithm; the tests turn it on. ``max_nodes`` aborts once the tree
    passes that size -- the tree is the thing that can blow up.
    """
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant {variant!r}; expected one of {VARIANTS}")
    if pivot_rule not in PIVOT_RULES:
        raise ValueError(
            f"unknown pivot rule {pivot_rule!r}; expected one of {PIVOT_RULES}"
        )
    if G.n != H.n:
        raise ValueError(
            f"G is on {G.n} vertices but H is on {H.n}; "
            "both sides must share a ground set"
        )

    tree = RecursionTree(
        instance=instance, algorithm="FK-A", variant=variant, pivot_rule=pivot_rule
    )
    counter = {"next_id": 1}

    def recurse(
        G_in: Hypergraph,
        H_in: Hypergraph,
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

        # Step 1: Sperner reduction, recording what it removed.
        rg, rh = sperner_reduce(G_in), sperner_reduce(H_in)
        G, H = rg.reduced, rh.reduced

        node = RecursionNode(
            node_id=node_id,
            parent_id=parent_id,
            path=path,
            G_in=G_in,
            H_in=H_in,
            G=G,
            H=H,
            removed_G=rg.removed,
            removed_H=rh.removed,
            epsilon=eps(G, H),
        )
        tree.add(node)

        def settle(verdict: Verdict) -> Verdict:
            if validate and verdict.certificate is not None:
                if not verify_certificate(G, H, verdict.certificate):
                    raise AssertionError(
                        f"node {node_id}: certificate "
                        f"{[v + 1 for v in verts(verdict.certificate)]} "
                        f"fails equation 2.1 for {G.label()} / {H.label()}"
                    )
            node.verdict = verdict
            return verdict

        # Step 2: conditions.
        failed = check_conditions(G, H)
        if failed is not None:
            # Algorithm 1 reports a structural witness here. For a research
            # artifact, equation 2.1's set certificate is more useful: it can be
            # checked independently at every ancestor after the lifts. Not every
            # degenerate pair admits one in this orientation, so keep the
            # structural verdict when the exact search returns None.
            if failed.certificate is None:
                found = find_certificate(G, H)
                if found is not None:
                    failed = replace(failed, certificate=found)
            return settle(failed)

        # Step 3: trivial instances, and the modified single-edge rule.
        trivial = _trivial_case(G, H, variant)
        if trivial is not None:
            return settle(trivial)

        # Step 4: recursive decomposition.
        pivot = choose_pivot(G, H, pivot_rule)
        if pivot is None:
            # Unreachable for Sperner-reduced inputs past Step 2: every edge is
            # non-empty there, so some vertex has positive degree.
            raise AssertionError(
                f"node {node_id}: no vertex of positive degree in "
                f"{G.label()} / {H.label()}"
            )
        node.pivot = pivot
        node.pivot_frequency = eps_v(G, H, pivot)

        G0, G1 = G.contraction(pivot), G.deletion(pivot)
        H0, H1 = H.contraction(pivot), H.deletion(pivot)

        # L: (G_1, H_0 or H_1)
        left = recurse(G1, H0.union(H1), node_id, path + ("L",))
        if not left.dual:
            # S_L misses the pivot, so S_L + {v} still meets every edge of G --
            # edges through v are hit by v -- and still contains no edge of H.
            return settle(
                Verdict(
                    dual=False,
                    reason="child_failed",
                    detail=f"left subproblem (G1, H0 or H1) on pivot v{pivot + 1} failed",
                    certificate=(
                        None
                        if left.certificate is None
                        else left.certificate | (1 << pivot)
                    ),
                    witness_edges=left.witness_edges,
                )
            )

        # R: (H_1, G_0 or G_1)
        right = recurse(H1, G0.union(G1), node_id, path + ("R",))
        if not right.dual:
            # The right subproblem has the sides swapped, so its certificate is
            # complemented within V \ {v} (thesis p.11: S works for (G,H) iff
            # V \ S works for (H,G)).
            cert = None
            if right.certificate is not None:
                cert = (((1 << G.n) - 1) & ~(1 << pivot)) & ~right.certificate
            return settle(
                Verdict(
                    dual=False,
                    reason="child_failed",
                    detail=f"right subproblem (H1, G0 or G1) on pivot v{pivot + 1} failed",
                    certificate=cert,
                    witness_edges=right.witness_edges,
                )
            )

        return settle(
            Verdict(dual=True, reason="dual", detail="both subproblems returned TRUE")
        )

    verdict = recurse(G, H, None, ())
    tree.dual = verdict.dual
    tree.verdict = verdict
    return tree


def is_dual(
    G: Hypergraph,
    H: Hypergraph,
    *,
    variant: str = "modified",
    pivot_rule: str = "degree_sum",
) -> bool:
    """``True`` iff ``G = Tr(H)``. Discards the tree."""
    return bool(fk_a(G, H, variant=variant, pivot_rule=pivot_rule).dual)
