"""``Aut(H)``: the vertex permutations mapping the edge set onto itself.

Thesis Definition 3.1.3 (p.18). Annotating every recursion node with it is the
central experiment of Chapter 3.

Isolated vertices are excluded by default. Contraction removes vertices from
edges, so subhypergraphs deep in a recursion routinely have vertices in no edge;
those are freely permutable and would multiply ``|Aut(H)|`` by ``k!``, drowning
the structure being measured. The thesis is explicit that "isolated vertices are
not permitted" (p.53).

The two legacy notebooks disagreed here and would report different orders for
the same node: the property-tests notebook restricted to active vertices
(``include_isolated=False``), while the Sage notebook appended a singleton block
per zero column, inflating the group. ``include_isolated=True`` reproduces the
second convention for reconciling an old result. It is not the default.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from typing import Any, Optional, Sequence

from .groups import (
    MAX_GROUP_ORDER,
    GroupIdentification,
    Perm,
    PermutationGroup,
    cycle_notation,
    identify,
)
from .hypergraph import Hypergraph, verts, popcount

__all__ = ["AutomorphismResult", "automorphism_group", "hypergraph_automorphisms"]


# ----------------------------------------------------------------------
# vertex colouring (search pruning)
# ----------------------------------------------------------------------
def _initial_colours(H: Hypergraph) -> list[tuple]:
    """A permutation-invariant signature per vertex: degree and incident sizes."""
    sizes: list[list[int]] = [[] for _ in range(H.n)]
    for e in H.edges:
        k = popcount(e)
        for v in verts(e):
            sizes[v].append(k)
    return [(len(sizes[v]), tuple(sorted(sizes[v]))) for v in range(H.n)]


def _refine_colours(H: Hypergraph, colours: list[tuple], rounds: int = 3) -> list[tuple]:
    """Iteratively refine by the multiset of neighbour colours within each edge.

    Standard 1-dimensional Weisfeiler-Leman style refinement. It never
    separates vertices that are actually equivalent, so it only ever prunes
    search branches that could not succeed.
    """
    incident: list[list[int]] = [[] for _ in range(H.n)]
    for e in H.edges:
        for v in verts(e):
            incident[v].append(e)
    current = list(colours)
    for _ in range(rounds):
        new: list[tuple] = []
        for v in range(H.n):
            sig = []
            for e in incident[v]:
                sig.append(tuple(sorted(current[u] for u in verts(e) if u != v)))
            new.append((current[v], tuple(sorted(sig))))
        # Compress to small integers keyed by first appearance, so the tuples
        # do not grow without bound across rounds.
        codes: dict[tuple, int] = {}
        compressed = []
        for c in new:
            compressed.append((codes.setdefault(c, len(codes)),))
        if compressed == current:
            break
        current = compressed
    return current


# ----------------------------------------------------------------------
# the search
# ----------------------------------------------------------------------
def hypergraph_automorphisms(
    H: Hypergraph, *, cap: int = MAX_GROUP_ORDER
) -> list[Perm]:
    """Every automorphism of ``H``, as permutations of ``0 .. H.n-1``.

    Backtracking with colour-class pruning and partial-edge feasibility: a
    partial map is abandoned as soon as some edge's assigned image cannot be
    extended to any edge of the same size. Exhaustive and exact.

    Raises :class:`ValueError` if the group exceeds ``cap`` elements -- which in
    practice means the input has a large free part (many isolated or
    interchangeable vertices) and should be compacted first.
    """
    n = H.n
    if n == 0:
        return [()]
    edge_set = set(H.edges)
    if not edge_set:
        # No edges: every permutation is an automorphism. Enumerating n! is
        # pointless and explosive; say so rather than pretend to compute it.
        from math import factorial

        if factorial(n) > cap:
            raise ValueError(
                f"hypergraph has no edges, so Aut is the full symmetric group S{n} "
                f"of order {factorial(n)}, above the cap of {cap}"
            )
        return [tuple(p) for p in _all_permutations(n)]

    by_size: dict[int, list[int]] = {}
    for e in edge_set:
        by_size.setdefault(popcount(e), []).append(e)

    colours = _refine_colours(H, _initial_colours(H))
    classes: dict[tuple, list[int]] = {}
    for v in range(n):
        classes.setdefault(colours[v], []).append(v)

    incident: list[list[int]] = [[] for _ in range(n)]
    for e in edge_set:
        for v in verts(e):
            incident[v].append(e)

    # Assign the most constrained vertices first: smallest colour class, then
    # highest degree.
    order = sorted(
        range(n), key=lambda v: (len(classes[colours[v]]), -len(incident[v]), v)
    )

    img = [-1] * n
    used = 0
    found: list[Perm] = []

    def feasible(v: int) -> bool:
        """Check every edge through ``v`` after ``v`` has just been assigned."""
        for e in incident[v]:
            partial = 0
            complete = True
            for u in verts(e):
                if img[u] < 0:
                    complete = False
                else:
                    partial |= 1 << img[u]
            if complete:
                if partial not in edge_set:
                    return False
            else:
                candidates = by_size.get(popcount(e))
                if not candidates:
                    return False
                if not any((partial & f) == partial for f in candidates):
                    return False
        return True

    def backtrack(depth: int) -> None:
        nonlocal used
        if len(found) > cap:
            raise ValueError(
                f"automorphism group exceeded {cap} elements; "
                "compact the hypergraph or raise the cap"
            )
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

    backtrack(0)
    return sorted(found)


def _all_permutations(n: int) -> list[list[int]]:
    from itertools import permutations

    return [list(p) for p in permutations(range(n))]


# ----------------------------------------------------------------------
# public result type
# ----------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class AutomorphismResult:
    """``Aut(H)`` with everything the recursion-tree annotation needs.

    ``generators`` and ``orbits`` are expressed in the *original* 1-indexed
    vertex numbering of the hypergraph that was passed in, not the compacted
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
    H: Hypergraph,
    *,
    include_isolated: bool = False,
    backend: Optional[str] = None,
    cap: int = MAX_GROUP_ORDER,
) -> AutomorphismResult:
    """Compute and name ``Aut(H)``.

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
        target, mapping = H, list(range(H.n))
    else:
        target, mapping = H.compact()

    isolated = H.isolated_vertices()

    perms = hypergraph_automorphisms(target, cap=cap)
    group = PermutationGroup(degree=target.n, elements=perms)
    group.generators = group.minimal_generators()

    if chosen == "sage":
        from .backends import sage_backend

        naming = sage_backend.name_group(group)
    else:
        naming = identify(group)

    # Translate generators and orbits back to the caller's vertex numbering.
    def to_original(p: Perm) -> Perm:
        out = list(range(H.n))
        for i, j in enumerate(p):
            out[mapping[i]] = mapping[j]
        return tuple(out)

    gens = tuple(cycle_notation(to_original(g)) for g in group.generators)
    orbits = tuple(
        tuple(sorted(mapping[i] + 1 for i in orb)) for orb in group.orbits()
    )

    return AutomorphismResult(
        order=group.order,
        name=naming.name,
        generators=gens,
        orbits=orbits,
        active_vertices=tuple(mapping),
        isolated_vertices=isolated,
        abelian=naming.abelian,
        backend=chosen,
        naming_method=naming.method,
    )
