"""``Aut(MBF)``: the variable permutations mapping the term set onto itself.

Thesis Definition 3.1.3 (p.18). Annotating every recursion node with it is the
central experiment of Chapter 3.

Isolated variables are excluded by default. ``phi_x_0`` drops variables from
terms, so subproblems deep in a recursion routinely have variables in no term;
those are freely permutable and would multiply ``|Aut(MBF)|`` by ``k!``, drowning
the structure being measured. The thesis is explicit that "isolated vertices are
not permitted" (p.53). ``include_isolated=True`` keeps them, to reconcile a
result computed under the other convention.

:func:`term_set_automorphisms` is the one search. Given one term set it is
``Aut(MBF)``; given two it is the pair group of a recursion node, where a clause
may not be mapped onto a monomial.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations
from math import factorial
from typing import Any, Optional, Sequence

from .groups import (
    MAX_GROUP_ORDER,
    Perm,
    PermutationGroup,
    cycle_notation,
    identify,
)
from .hypergraph import Hypergraph, popcount, verts

__all__ = [
    "AutomorphismResult",
    "automorphism_group",
    "hypergraph_automorphisms",
    "term_set_automorphisms",
]

#: Rounds of colour refinement. Each round costs one pass over the terms and can
#: only split classes further, so it converges; three is where it stops paying
#: on everything in the instance library.
REFINEMENT_ROUNDS = 3


def _colours(families: Sequence[Hypergraph], n: int) -> list[tuple]:
    """A permutation-invariant colour per variable.

    Seeded with, per family, the variable's degree and the sizes of the terms
    holding it, then refined by the multiset of its co-members' colours --
    1-dimensional Weisfeiler-Leman. Refinement never separates two variables
    that some automorphism identifies, so it only prunes branches that could not
    have succeeded.

    Family membership is part of the colour, which is what keeps the two sides
    of a pair apart.
    """
    seed: list[list] = [[] for _ in range(n)]
    incident: list[list[tuple[int, int]]] = [[] for _ in range(n)]
    for f, MBF in enumerate(families):
        sizes: list[list[int]] = [[] for _ in range(n)]
        for t in MBF.edges:
            k = popcount(t)
            for v in verts(t):
                sizes[v].append(k)
                incident[v].append((f, t))
        for v in range(n):
            seed[v].append((len(sizes[v]), tuple(sorted(sizes[v]))))

    current: list[tuple] = [tuple(c) for c in seed]
    for _ in range(REFINEMENT_ROUNDS):
        fresh = []
        for v in range(n):
            sig = sorted(
                (f, tuple(sorted(current[u] for u in verts(t) if u != v)))
                for f, t in incident[v]
            )
            fresh.append((current[v], tuple(sig)))
        # Compress to small integers keyed by first appearance, so the tuples do
        # not grow without bound across rounds.
        codes: dict[tuple, int] = {}
        compressed: list[tuple] = [(codes.setdefault(c, len(codes)),) for c in fresh]
        if compressed == current:
            break
        current = compressed
    return current


def term_set_automorphisms(
    families: Sequence[Hypergraph],
    *,
    cap: int = MAX_GROUP_ORDER,
    partial: bool = False,
) -> tuple[list[Perm], bool]:
    """Permutations mapping each family's term set onto itself, and whether the
    enumeration stopped early.

    One family is ``Aut(MBF)``. Two is the pair group of a recursion node, where
    the sides need not be a dual pair and ``Aut(C) = Aut(D)`` (Lemma 3.3.1) no
    longer applies.

    Backtracking over the colour classes, with partial-term feasibility: a map
    is abandoned as soon as some term's assigned image cannot be extended to a
    term of the same size in the same family.

    ``partial=True`` stops at ``cap`` and returns what it has, which is a
    subgroup; ``partial=False`` raises instead. A subgroup still *proves*
    primitivity -- every block of a group is a block of each of its subgroups --
    so the cap can only under-report, never invent symmetry.
    """
    if not families:
        raise ValueError("need at least one term set")
    n = families[0].n
    if any(f.n != n for f in families):
        raise ValueError(
            f"families are on {[f.n for f in families]} variables; they must share one"
        )
    if n == 0:
        return [()], False

    sets = [set(f.edges) for f in families]
    if not any(sets):
        # No terms at all: every permutation qualifies. Enumerating n! is
        # pointless and explosive; say so rather than pretend to compute it.
        if factorial(n) > cap:
            if partial:
                found = [tuple(p) for _, p in zip(range(cap), permutations(range(n)))]
                return found, True
            raise ValueError(
                f"no terms, so Aut is the full symmetric group S{n} of order "
                f"{factorial(n)}, above the cap of {cap}"
            )
        return [tuple(p) for p in permutations(range(n))], False

    by_size: list[dict[int, list[int]]] = []
    for edges in sets:
        sized: dict[int, list[int]] = {}
        for t in edges:
            sized.setdefault(popcount(t), []).append(t)
        by_size.append(sized)

    incident: list[list[tuple[int, int]]] = [[] for _ in range(n)]
    for f, edges in enumerate(sets):
        for t in edges:
            for v in verts(t):
                incident[v].append((f, t))

    colours = _colours(families, n)
    classes: dict[tuple, list[int]] = {}
    for v in range(n):
        classes.setdefault(colours[v], []).append(v)

    # Assign the most constrained variables first: smallest colour class, then
    # highest degree.
    order = sorted(
        range(n), key=lambda v: (len(classes[colours[v]]), -len(incident[v]), v)
    )

    img = [-1] * n
    used = 0
    found: list[Perm] = []
    capped = False

    def feasible(v: int) -> bool:
        """Check every term through ``v`` after ``v`` has just been assigned."""
        for f, t in incident[v]:
            image = 0
            complete = True
            for u in verts(t):
                if img[u] < 0:
                    complete = False
                else:
                    image |= 1 << img[u]
            if complete:
                if image not in sets[f]:
                    return False
            else:
                candidates = by_size[f].get(popcount(t))
                if not candidates:
                    return False
                if not any((image & c) == image for c in candidates):
                    return False
        return True

    def backtrack(depth: int) -> None:
        nonlocal used, capped
        if len(found) >= cap:
            if not partial:
                raise ValueError(
                    f"automorphism group exceeded {cap} elements; "
                    "compact the input or raise the cap"
                )
            capped = True
            return
        if depth == n:
            found.append(tuple(img))
            return
        v = order[depth]
        for target in classes[colours[v]]:
            if used >> target & 1:
                continue
            img[v] = target
            used |= 1 << target
            if feasible(v):
                backtrack(depth + 1)
            used &= ~(1 << target)
            img[v] = -1
            if capped:
                return

    backtrack(0)
    return sorted(found), capped


def hypergraph_automorphisms(
    MBF: Hypergraph, *, cap: int = MAX_GROUP_ORDER
) -> list[Perm]:
    """``Aut(MBF)`` as permutations of ``0 .. MBF.n-1``. Raises above ``cap``."""
    return term_set_automorphisms([MBF], cap=cap)[0]


# ----------------------------------------------------------------------
# public result type
# ----------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class AutomorphismResult:
    """``Aut(MBF)`` with everything the recursion-tree annotation needs.

    ``generators`` and ``orbits`` are expressed in the *original* 1-indexed
    variable numbering of the function that was passed in, not the compacted
    numbering used internally, so they can be read against the instance.
    """

    order: int
    name: str
    generators: tuple[str, ...]
    orbits: tuple[tuple[int, ...], ...]
    active_vertices: tuple[int, ...]
    isolated_vertices: tuple[int, ...]
    abelian: bool
    backend: str
    naming_method: str

    @property
    def is_trivial(self) -> bool:
        return self.order == 1

    def to_json(self) -> dict[str, Any]:
        return {
            "order": self.order,
            "name": self.name,
            "generators": list(self.generators),
            "orbits": [list(o) for o in self.orbits],
            "active_vertices": [v + 1 for v in self.active_vertices],
            "isolated_vertices": [v + 1 for v in self.isolated_vertices],
            "abelian": self.abelian,
            "backend": self.backend,
            "naming_method": self.naming_method,
        }


def automorphism_group(
    MBF: Hypergraph,
    *,
    include_isolated: bool = False,
    backend: Optional[str] = None,
    cap: int = MAX_GROUP_ORDER,
) -> AutomorphismResult:
    """Compute and name ``Aut(MBF)``.

    ``backend`` selects the implementation: ``"python"`` for the built-in
    search, ``"sage"`` to delegate group naming to GAP through SageMath, or
    ``None`` (default) to use SageMath when it is importable and fall back
    otherwise. The group *order*, *generators* and *orbits* agree between
    backends; only the structure *name* can differ, and only in that the Python
    catalogue may return an ``order N`` placeholder where GAP gives a full
    decomposition.
    """
    from .backends import resolve_backend

    chosen = resolve_backend(backend)

    if include_isolated:
        target, mapping = MBF, list(range(MBF.n))
    else:
        target, mapping = MBF.compact()

    isolated = MBF.isolated_vertices()

    group = PermutationGroup(
        degree=target.n, elements=hypergraph_automorphisms(target, cap=cap)
    )
    group.generators = group.minimal_generators()

    if chosen == "sage":
        from .backends import sage_backend

        naming = sage_backend.name_group(group)
    else:
        naming = identify(group)

    # Translate generators and orbits back to the caller's variable numbering.
    def to_original(p: Perm) -> Perm:
        out = list(range(MBF.n))
        for i, j in enumerate(p):
            out[mapping[i]] = mapping[j]
        return tuple(out)

    return AutomorphismResult(
        order=group.order,
        name=naming.name,
        generators=tuple(cycle_notation(to_original(g)) for g in group.generators),
        orbits=tuple(
            tuple(sorted(mapping[i] + 1 for i in orb)) for orb in group.orbits()
        ),
        active_vertices=tuple(mapping),
        isolated_vertices=isolated,
        abelian=naming.abelian,
        backend=chosen,
        naming_method=naming.method,
    )
