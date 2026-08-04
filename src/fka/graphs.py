"""Primal graphs and the graph classes the thesis tracks.

Chapter 5 annotates recursion-tree nodes with the class of the input
hypergraph's primal graph -- cograph, chordal, split, bipartite, chordal
bipartite, perfect -- because those classes are subclasses of
distance-hereditary graphs and therefore carry twin-vertex structure
(thesis Section 5.2, p.48).

Everything here is exact and runs in plain Python. Several tests enumerate
induced subgraphs, which is exponential in the vertex count; that is fine for
the instances in the thesis (at most 14 vertices) but is guarded by
:data:`MAX_ENUMERATION_VERTICES` so a larger instance fails loudly rather than
hanging. SageMath has polynomial-time versions of some of these; see
:mod:`fka.backends.sage_backend`.

Vertices are 0-indexed integers, matching :mod:`fka.hypergraph`.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Iterable, Iterator, Optional

import networkx as nx

from .hypergraph import Hypergraph, bitset_to_vertices

__all__ = ["Graph", "primal_graph", "GraphClasses", "classify_graph"]

#: Above this many vertices, the induced-subgraph enumerations below are
#: refused rather than attempted.
MAX_ENUMERATION_VERTICES = 20


@dataclass(frozen=True, slots=True)
class Graph:
    """A simple undirected graph on ``0 .. n-1``, adjacency stored as bitsets."""

    n: int
    adj: tuple[int, ...]

    @classmethod
    def from_edges(cls, n: int, edges: Iterable[tuple[int, int]]) -> "Graph":
        adj = [0] * n
        for u, v in edges:
            if u == v:
                continue
            adj[u] |= 1 << v
            adj[v] |= 1 << u
        return cls(n, tuple(adj))

    def edges(self) -> list[tuple[int, int]]:
        return [
            (u, v)
            for u in range(self.n)
            for v in bitset_to_vertices(self.adj[u])
            if u < v
        ]

    @property
    def size(self) -> int:
        return len(self.edges())

    def has_edge(self, u: int, v: int) -> bool:
        return bool(self.adj[u] >> v & 1)

    def degree(self, v: int) -> int:
        return self.adj[v].bit_count()

    def complement(self) -> "Graph":
        full = (1 << self.n) - 1
        return Graph(
            self.n,
            tuple(full & ~self.adj[v] & ~(1 << v) for v in range(self.n)),
        )

    def to_networkx(self) -> nx.Graph:
        g = nx.Graph()
        g.add_nodes_from(range(self.n))
        g.add_edges_from(self.edges())
        return g

    def induced(self, vertices: Iterable[int]) -> "Graph":
        """Induced subgraph, re-indexed onto ``vertices`` in the given order."""
        vs = list(vertices)
        idx = {v: i for i, v in enumerate(vs)}
        adj = [0] * len(vs)
        for i, u in enumerate(vs):
            for v in bitset_to_vertices(self.adj[u]):
                if v in idx:
                    adj[i] |= 1 << idx[v]
        return Graph(len(vs), tuple(adj))

    def maximal_cliques(self) -> list[tuple[int, ...]]:
        return [tuple(sorted(c)) for c in nx.find_cliques(self.to_networkx())]

    def is_connected(self) -> bool:
        if self.n == 0:
            return True
        seen = 1
        frontier = [0]
        while frontier:
            v = frontier.pop()
            for u in bitset_to_vertices(self.adj[v] & ~seen):
                seen |= 1 << u
                frontier.append(u)
        return seen == (1 << self.n) - 1


def primal_graph(H: Hypergraph) -> Graph:
    """The primal (2-section) graph of ``H``.

    Two vertices are adjacent iff some hyperedge contains both. Isolated
    vertices of ``H`` stay as isolated vertices here, so the graph is always on
    the same ground set as ``H``.
    """
    adj = [0] * H.n
    for e in H.edges:
        vs = bitset_to_vertices(e)
        for u in vs:
            adj[u] |= e & ~(1 << u)
    return Graph(H.n, tuple(adj))


# ----------------------------------------------------------------------
# induced-subgraph tests
# ----------------------------------------------------------------------
def _check_size(g: Graph, what: str) -> None:
    if g.n > MAX_ENUMERATION_VERTICES:
        raise ValueError(
            f"{what} enumerates induced subgraphs and is capped at "
            f"{MAX_ENUMERATION_VERTICES} vertices; this graph has {g.n}. "
            "Use the SageMath backend or raise MAX_ENUMERATION_VERTICES."
        )


def find_induced_p4(g: Graph) -> Optional[tuple[int, ...]]:
    """Return the vertices of an induced ``P4``, or ``None``.

    ``P4`` is the path on 4 vertices; a graph is a cograph exactly when it has
    none (Corneil-Lerchs-Stewart). Checking all 4-subsets is ``O(n^4)``, which
    beats writing a cotree decomposition for graphs this size.
    """
    _check_size(g, "cograph testing")
    for quad in combinations(range(g.n), 4):
        sub = g.induced(quad)
        degs = sorted(sub.degree(i) for i in range(4))
        # A P4 is the unique 4-vertex graph with 3 edges and degrees 1,1,2,2.
        if sub.size == 3 and degs == [1, 1, 2, 2] and sub.is_connected():
            return quad
    return None


def is_cograph(g: Graph) -> bool:
    """``P4``-free."""
    return find_induced_p4(g) is None


def _induced_cycles_at_least(g: Graph, length: int) -> Iterator[tuple[int, ...]]:
    """Yield vertex sets inducing a chordless cycle of size >= ``length``."""
    _check_size(g, "induced cycle enumeration")
    for k in range(length, g.n + 1):
        for subset in combinations(range(g.n), k):
            sub = g.induced(subset)
            # A chordless k-cycle is exactly: k edges, every degree 2, connected.
            if sub.size != k:
                continue
            if any(sub.degree(i) != 2 for i in range(k)):
                continue
            if sub.is_connected():
                yield subset


def is_chordal(g: Graph) -> bool:
    return nx.is_chordal(g.to_networkx())


def is_bipartite(g: Graph) -> bool:
    return nx.is_bipartite(g.to_networkx())


def is_split(g: Graph) -> bool:
    """Split iff both the graph and its complement are chordal (Foldes-Hammer)."""
    return is_chordal(g) and is_chordal(g.complement())


def is_weakly_chordal(g: Graph) -> bool:
    """No induced cycle of length >= 5 in the graph or its complement."""
    if next(_induced_cycles_at_least(g, 5), None) is not None:
        return False
    return next(_induced_cycles_at_least(g.complement(), 5), None) is None


def is_perfect(g: Graph) -> bool:
    """Strong Perfect Graph Theorem: no odd hole and no odd antihole."""
    for cyc in _induced_cycles_at_least(g, 5):
        if len(cyc) % 2 == 1:
            return False
    for cyc in _induced_cycles_at_least(g.complement(), 5):
        if len(cyc) % 2 == 1:
            return False
    return True


def is_chordal_bipartite(g: Graph) -> bool:
    """Bipartite, with every cycle of length >= 6 chorded.

    Equivalently: bipartite with no induced cycle of length >= 6. (Bipartite
    graphs have no odd cycles, so length 5 cannot occur.)
    """
    if not is_bipartite(g):
        return False
    return next(_induced_cycles_at_least(g, 6), None) is None


def is_triangle_free(g: Graph) -> bool:
    return all(
        not (g.adj[u] & g.adj[v])
        for u, v in g.edges()
    )


@dataclass(frozen=True, slots=True)
class GraphClasses:
    """Which classes a primal graph belongs to."""

    n: int
    m: int
    cograph: bool
    chordal: bool
    bipartite: bool
    split: bool
    weakly_chordal: bool
    perfect: bool
    chordal_bipartite: bool
    triangle_free: bool

    def names(self) -> list[str]:
        """Class names in the order the thesis lists them on tree nodes."""
        out = []
        for label, flag in (
            ("Cograph", self.cograph),
            ("Chordal", self.chordal),
            ("Split", self.split),
            ("Bipartite", self.bipartite),
            ("Chordal bipartite", self.chordal_bipartite),
            ("Weakly chordal", self.weakly_chordal),
            ("Perfect", self.perfect),
            ("Triangle-free", self.triangle_free),
        ):
            if flag:
                out.append(label)
        return out

    def to_json(self) -> dict:
        return {
            "n": self.n,
            "m": self.m,
            "cograph": self.cograph,
            "chordal": self.chordal,
            "bipartite": self.bipartite,
            "split": self.split,
            "weakly_chordal": self.weakly_chordal,
            "perfect": self.perfect,
            "chordal_bipartite": self.chordal_bipartite,
            "triangle_free": self.triangle_free,
            "names": self.names(),
        }


def classify_graph(g: Graph) -> GraphClasses:
    """Run every class test on ``g``."""
    return GraphClasses(
        n=g.n,
        m=g.size,
        cograph=is_cograph(g),
        chordal=is_chordal(g),
        bipartite=is_bipartite(g),
        split=is_split(g),
        weakly_chordal=is_weakly_chordal(g),
        perfect=is_perfect(g),
        chordal_bipartite=is_chordal_bipartite(g),
        triangle_free=is_triangle_free(g),
    )
