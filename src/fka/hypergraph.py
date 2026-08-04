"""The hypergraph type, stored as bitsets.

Hyperedge ``e`` is the integer whose bit ``i`` is set iff vertex ``i`` is in
``e``. Vertices are 0-indexed here; the thesis numbers from 1, so every
human-facing path converts with ``+1``.

Bitsets rather than the legacy incidence rows because subset, superset and
disjointness become single machine words -- ``a & b == a`` is ``a subseteq b``
exactly -- the FK-A split is ``e & ~bit``, and edges become hashable, so
deduplication is a ``set`` pass. Incidence matrices stay the interchange format
for SageMath and MATLAB, so ``from_incidence`` / ``to_incidence`` round-trip.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Iterator, Sequence

__all__ = ["Hypergraph", "verts", "bits", "popcount"]


def popcount(x: int) -> int:
    """Number of vertices in the hyperedge ``x``."""
    return x.bit_count()


def bits(vertices: Iterable[int]) -> int:
    """Pack 0-indexed vertex ids into a bitset."""
    e = 0
    for v in vertices:
        e |= 1 << v
    return e


def verts(e: int) -> tuple[int, ...]:
    """Unpack a bitset into ascending 0-indexed vertex ids."""
    out = []
    while e:
        low = e & -e
        out.append(low.bit_length() - 1)
        e ^= low
    return tuple(out)


@dataclass(frozen=True, slots=True)
class Hypergraph:
    """An immutable hypergraph on vertices ``0 .. n-1``.

    ``edges`` is deduplicated and sorted, so equal edge sets on equal ground
    sets compare and hash equal. The experiments need that: the thesis observes
    FK-A regenerating identical subproblems (p.53), and detecting one is an
    equality test on this type.

    ``n`` is carried explicitly rather than inferred. Contraction removes
    vertices from edges, so a subhypergraph deep in the recursion can have
    vertices in no edge at all, and ``n`` keeps every node's coordinates aligned
    with the root.
    """

    n: int
    edges: tuple[int, ...]

    def __post_init__(self) -> None:
        if self.n < 0:
            raise ValueError(f"vertex count must be non-negative, got {self.n}")
        mask = (1 << self.n) - 1
        for e in self.edges:
            if e & ~mask:
                raise ValueError(
                    f"hyperedge {bin(e)} references a vertex outside 0..{self.n - 1}"
                )
        # Direct construction is public and useful for the degenerate ``{{}}``,
        # so enforce the canonical invariant here too, not only in from_bitsets.
        canonical = tuple(sorted(set(self.edges), key=verts))
        if canonical != self.edges:
            object.__setattr__(self, "edges", canonical)

    # ------------------------------------------------------------------
    # construction
    # ------------------------------------------------------------------
    @classmethod
    def from_bitsets(cls, n: int, edges: Iterable[int]) -> "Hypergraph":
        """Build from raw bitsets, deduplicating and sorting lexicographically.

        Sorted by the unpacked vertex tuple, not numerically: ``{4,5}`` is a
        smaller integer than ``{1,5,6}`` but comes later in the vertex-list
        notation the thesis and every export use.
        """
        return cls(n, tuple(sorted(set(edges), key=verts)))

    @classmethod
    def from_sets(
        cls,
        n: int,
        edges: Iterable[Sequence[int]],
        *,
        one_indexed: bool = False,
    ) -> "Hypergraph":
        """Build from vertex lists. ``one_indexed`` accepts the thesis' ``v1..vn``."""
        shift = 1 if one_indexed else 0
        packed = []
        for edge in edges:
            for v in edge:
                if not 0 <= v - shift < n:
                    raise ValueError(
                        f"vertex {v} out of range for a hypergraph on {n} vertices "
                        f"({'1' if one_indexed else '0'}-indexed input)"
                    )
            packed.append(bits(v - shift for v in edge))
        return cls.from_bitsets(n, packed)

    @classmethod
    def from_incidence(cls, matrix: Sequence[Sequence[int]]) -> "Hypergraph":
        """Build from a 0/1 incidence matrix, one row per hyperedge.

        Takes anything row-iterable, including numpy arrays and Sage matrices.
        """
        rows = [list(row) for row in matrix]
        if not rows:
            return cls(0, ())
        widths = {len(r) for r in rows}
        if len(widths) != 1:
            raise ValueError(f"incidence matrix has ragged rows: widths {sorted(widths)}")
        n = widths.pop()
        packed = []
        for i, row in enumerate(rows):
            for j, val in enumerate(row):
                if int(val) not in (0, 1):
                    raise ValueError(
                        f"incidence matrix entry ({i},{j}) is {val!r}; expected 0 or 1"
                    )
            packed.append(bits(j for j, val in enumerate(row) if int(val)))
        return cls.from_bitsets(n, packed)

    @classmethod
    def empty(cls, n: int) -> "Hypergraph":
        """No edges -- FK-A's ``|G| = 0``, and the constant 0 read as a DNF."""
        return cls(n, ())

    # ------------------------------------------------------------------
    # export
    # ------------------------------------------------------------------
    def to_incidence(self) -> list[list[int]]:
        """A 0/1 incidence matrix, as a list of rows."""
        return [[(e >> j) & 1 for j in range(self.n)] for e in self.edges]

    def to_sets(self, *, one_indexed: bool = False) -> list[list[int]]:
        """Vertex lists, ascending within each edge."""
        shift = 1 if one_indexed else 0
        return [[v + shift for v in verts(e)] for e in self.edges]

    def label(self, *, one_indexed: bool = True) -> str:
        """Compact human form, e.g. ``{{1,3},{2,4}}``."""
        if not self.edges:
            return "{}"
        return (
            "{"
            + ",".join(
                "{" + ",".join(str(v) for v in edge) + "}"
                for edge in self.to_sets(one_indexed=one_indexed)
            )
            + "}"
        )

    # ------------------------------------------------------------------
    # basic structure
    # ------------------------------------------------------------------
    def __len__(self) -> int:
        return len(self.edges)

    def __iter__(self) -> Iterator[int]:
        return iter(self.edges)

    def __bool__(self) -> bool:
        # __len__ would give this anyway; stating it keeps the algorithms'
        # empty-instance checks explicit.
        return bool(self.edges)

    @property
    def mask(self) -> int:
        """Bitset of the whole ground set ``0 .. n-1``."""
        return (1 << self.n) - 1

    def support(self) -> int:
        """Vertices in at least one edge -- the ``union of e_i`` of Prop. 2.2.2.

        A transversal pair must have equal supports.
        """
        u = 0
        for e in self.edges:
            u |= e
        return u

    def degree(self, v: int) -> int:
        """Number of edges containing vertex ``v``."""
        bit = 1 << v
        return sum(1 for e in self.edges if e & bit)

    def degrees(self) -> list[int]:
        """Vertex degrees: the column sums of the incidence matrix."""
        out = [0] * self.n
        for e in self.edges:
            for v in verts(e):
                out[v] += 1
        return out

    def edge_sizes(self) -> list[int]:
        """Edge cardinalities: the row sums of the incidence matrix."""
        return [popcount(e) for e in self.edges]

    def isolated_vertices(self) -> tuple[int, ...]:
        """Vertices in no edge.

        Not allowed at the root, but contraction creates them mid-recursion
        (thesis p.53), so nodes report them rather than silently renumbering.
        """
        sup = self.support()
        return tuple(v for v in range(self.n) if not (sup >> v) & 1)

    # ------------------------------------------------------------------
    # FK-A primitives
    # ------------------------------------------------------------------
    def contraction(self, v: int) -> "Hypergraph":
        """``H/v``: edges containing ``v``, with ``v`` deleted (Def. 2.2.5).

        This is ``H_0``. It need not be Sperner.
        """
        bit = 1 << v
        return Hypergraph.from_bitsets(
            self.n, (e & ~bit for e in self.edges if e & bit)
        )

    def deletion(self, v: int) -> "Hypergraph":
        """``H \\ v``: edges avoiding ``v`` (Def. 2.2.6).

        This is ``H_1``. A subfamily of a Sperner family, so it stays Sperner.
        """
        bit = 1 << v
        return Hypergraph.from_bitsets(
            self.n, (e for e in self.edges if not e & bit)
        )

    def union(self, other: "Hypergraph") -> "Hypergraph":
        """Edge-set union -- the ``H_0 or H_1`` of the FK-A subproblems.

        Deliberately does *not* Sperner-reduce: FK-A reduces at the top of the
        next call, and keeping the two apart is what lets the visualiser report
        which edges each node's reduction removed.
        """
        if self.n != other.n:
            raise ValueError(
                f"cannot union hypergraphs on {self.n} and {other.n} vertices"
            )
        return Hypergraph.from_bitsets(self.n, self.edges + other.edges)

    def restricted_to(self, vertices: Sequence[int]) -> tuple["Hypergraph", list[int]]:
        """Re-index onto ``vertices``, returning the map new -> old.

        For handing a subhypergraph to a backend (SageMath, GAP) that would
        otherwise treat isolated vertices as real points and inflate the group.
        """
        remap = {old: new for new, old in enumerate(vertices)}
        packed = [
            bits(remap[v] for v in verts(e) if v in remap) for e in self.edges
        ]
        return Hypergraph.from_bitsets(len(vertices), packed), list(vertices)

    def compact(self) -> tuple["Hypergraph", list[int]]:
        """``restricted_to`` the support: drop isolated vertices."""
        return self.restricted_to(verts(self.support()))

    # ------------------------------------------------------------------
    # predicates
    # ------------------------------------------------------------------
    def is_sperner(self) -> bool:
        """True iff the edge set is an antichain."""
        for i, a in enumerate(self.edges):
            for j, b in enumerate(self.edges):
                if i != j and (a & b) == a:
                    return False
        return True

    def contains_edge(self, e: int) -> bool:
        return e in self.edges

    def __repr__(self) -> str:
        return f"Hypergraph(n={self.n}, edges={self.to_sets(one_indexed=True)})"
