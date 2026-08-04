"""Hypergraph representation for the FK-A study.

A hypergraph on ``n`` vertices is stored as a tuple of *bitsets*: hyperedge
``e`` is the Python integer whose bit ``i`` is set iff vertex ``i`` belongs to
``e``.  Vertices are 0-indexed internally; the thesis numbers them from 1, so
every human-facing path (labels, JSON, reports) converts with ``+1``.

Why bitsets rather than the legacy numpy incidence matrix:

* subset / superset / disjointness are single machine words -- ``a & b == a``
  tests ``a subseteq b`` exactly, with no float or dtype surprises;
* the FK-A split is ``e & ~bit`` and a partition on ``e & bit``, so the whole
  recursion is integer arithmetic;
* hyperedges become hashable, so duplicate removal is a ``dict``/``set`` pass
  rather than a sort.

Incidence matrices are still the interchange format (they are what the thesis,
SageMath and the MATLAB FK-B code all speak), so ``from_incidence`` /
``to_incidence`` are lossless round trips.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Iterator, Sequence

__all__ = ["Hypergraph", "bitset_to_vertices", "vertices_to_bitset", "popcount"]


def popcount(x: int) -> int:
    """Number of vertices in the hyperedge bitset ``x``."""
    return x.bit_count()


def vertices_to_bitset(vertices: Iterable[int]) -> int:
    """Pack 0-indexed vertex ids into a bitset."""
    e = 0
    for v in vertices:
        e |= 1 << v
    return e


def bitset_to_vertices(e: int) -> tuple[int, ...]:
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

    ``edges`` is deduplicated and sorted, so two hypergraphs with the same edge
    set on the same ground set compare and hash equal.  This matters in the
    experiments: the thesis observes repeated subproblems in the FK-A recursion
    tree (Section 5.2.2, thesis p.53), and detecting them is an equality test on
    this type.

    The vertex count ``n`` is carried explicitly rather than inferred from the
    edges.  FK-A contracts vertices out of hyperedges, so a subhypergraph deep
    in the recursion can have vertices that appear in no edge; ``n`` keeps the
    coordinates of every node aligned with the root instance.
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
        # Direct construction is intentionally public (and useful for the
        # degenerate ``{{}}`` family), so enforce the same canonical invariant
        # as ``from_bitsets`` here as well.
        canonical = tuple(sorted(set(self.edges), key=bitset_to_vertices))
        if canonical != self.edges:
            object.__setattr__(self, "edges", canonical)

    # ------------------------------------------------------------------
    # construction
    # ------------------------------------------------------------------
    @classmethod
    def from_bitsets(cls, n: int, edges: Iterable[int]) -> "Hypergraph":
        """Build from raw bitsets, deduplicating and sorting lexicographically.

        Numeric bitset order is an implementation detail: for example ``{4,5}``
        has a smaller integer than ``{1,5,6}``, despite coming later in the
        vertex-list notation used by the thesis and every exported format.
        Canonicalising by the unpacked vertex tuple keeps displays, JSON and
        incidence round trips stable and human-readable.
        """
        return cls(n, tuple(sorted(set(edges), key=bitset_to_vertices)))

    @classmethod
    def from_sets(
        cls,
        n: int,
        edges: Iterable[Sequence[int]],
        *,
        one_indexed: bool = False,
    ) -> "Hypergraph":
        """Build from vertex lists.

        ``one_indexed=True`` accepts the thesis' ``v1..vn`` numbering directly.
        """
        shift = 1 if one_indexed else 0
        packed = []
        for edge in edges:
            bits = 0
            for v in edge:
                idx = v - shift
                if not 0 <= idx < n:
                    raise ValueError(
                        f"vertex {v} out of range for a hypergraph on {n} vertices "
                        f"({'1' if one_indexed else '0'}-indexed input)"
                    )
                bits |= 1 << idx
            packed.append(bits)
        return cls.from_bitsets(n, packed)

    @classmethod
    def from_incidence(cls, matrix: Sequence[Sequence[int]]) -> "Hypergraph":
        """Build from a 0/1 incidence matrix: one row per hyperedge.

        Accepts anything row-iterable, including numpy arrays and Sage matrices,
        so long as entries compare equal to 0 or 1.
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
            bits = 0
            for j, val in enumerate(row):
                v = int(val)
                if v not in (0, 1):
                    raise ValueError(
                        f"incidence matrix entry ({i},{j}) is {val!r}; expected 0 or 1"
                    )
                if v:
                    bits |= 1 << j
            packed.append(bits)
        return cls.from_bitsets(n, packed)

    @classmethod
    def empty(cls, n: int) -> "Hypergraph":
        """The hypergraph with no edges -- FK-A's ``|G| = 0`` case."""
        return cls(n, ())

    # ------------------------------------------------------------------
    # export
    # ------------------------------------------------------------------
    def to_incidence(self) -> list[list[int]]:
        """Render as a 0/1 incidence matrix (list of rows)."""
        return [[(e >> j) & 1 for j in range(self.n)] for e in self.edges]

    def to_sets(self, *, one_indexed: bool = False) -> list[list[int]]:
        """Render as vertex lists, ascending within each edge."""
        shift = 1 if one_indexed else 0
        return [[v + shift for v in bitset_to_vertices(e)] for e in self.edges]

    def label(self, *, one_indexed: bool = True) -> str:
        """Compact human-readable form, e.g. ``{{1,3},{2,4}}``."""
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
        # An edgeless hypergraph is falsy; without this, `if not H` would use
        # __len__ anyway, but stating it keeps the empty-instance checks in the
        # algorithm explicit.
        return bool(self.edges)

    @property
    def mask(self) -> int:
        """Bitset of the whole ground set ``0 .. n-1``."""
        return (1 << self.n) - 1

    def support(self) -> int:
        """Bitset of vertices that appear in at least one hyperedge.

        This is the ``union of e_i`` of thesis Proposition 2.2.2 (p.12): a
        transversal pair must have equal supports.
        """
        u = 0
        for e in self.edges:
            u |= e
        return u

    def degree(self, v: int) -> int:
        """Number of hyperedges containing vertex ``v`` (0-indexed)."""
        bit = 1 << v
        return sum(1 for e in self.edges if e & bit)

    def degrees(self) -> list[int]:
        """Vertex degrees, i.e. the column sums of the incidence matrix."""
        out = [0] * self.n
        for e in self.edges:
            for v in bitset_to_vertices(e):
                out[v] += 1
        return out

    def edge_sizes(self) -> list[int]:
        """Hyperedge cardinalities, i.e. the row sums of the incidence matrix."""
        return [popcount(e) for e in self.edges]

    def isolated_vertices(self) -> tuple[int, ...]:
        """Vertices in no hyperedge.

        FK-A does not permit isolated vertices at the root, but contraction
        creates them mid-recursion (thesis p.53), so nodes report them rather
        than silently renumbering.
        """
        sup = self.support()
        return tuple(v for v in range(self.n) if not (sup >> v) & 1)

    # ------------------------------------------------------------------
    # FK-A primitives
    # ------------------------------------------------------------------
    def contraction(self, v: int) -> "Hypergraph":
        """``H/v``: edges containing ``v``, with ``v`` deleted (thesis Def. 2.2.5).

        This is ``H_0`` in the FK-A decomposition. It need not be Sperner.
        """
        bit = 1 << v
        return Hypergraph.from_bitsets(
            self.n, (e & ~bit for e in self.edges if e & bit)
        )

    def deletion(self, v: int) -> "Hypergraph":
        """``H \\ v``: edges avoiding ``v`` (thesis Def. 2.2.6).

        This is ``H_1``. Being a subfamily of a Sperner family, it stays Sperner.
        """
        bit = 1 << v
        return Hypergraph.from_bitsets(
            self.n, (e for e in self.edges if not e & bit)
        )

    def union(self, other: "Hypergraph") -> "Hypergraph":
        """Edge-set union -- the ``H_0 or H_1`` of the FK-A subproblems.

        Note this is a plain union: it does *not* Sperner-reduce. FK-A applies
        the reduction at the top of the next recursive call, and keeping the two
        steps apart is what lets the visualiser report which edges the reduction
        removed at each node.
        """
        if self.n != other.n:
            raise ValueError(
                f"cannot union hypergraphs on {self.n} and {other.n} vertices"
            )
        return Hypergraph.from_bitsets(self.n, self.edges + other.edges)

    def restricted_to(self, vertices: Sequence[int]) -> tuple["Hypergraph", list[int]]:
        """Drop unused coordinates.

        Returns the hypergraph re-indexed onto ``vertices`` together with the
        map new-index -> old-index. Used when handing a subhypergraph to a
        backend (SageMath, GAP) that would otherwise treat isolated vertices as
        real points and inflate the automorphism group.
        """
        remap = {old: new for new, old in enumerate(vertices)}
        packed = []
        for e in self.edges:
            bits = 0
            for v in bitset_to_vertices(e):
                if v in remap:
                    bits |= 1 << remap[v]
            packed.append(bits)
        return Hypergraph.from_bitsets(len(vertices), packed), list(vertices)

    def compact(self) -> tuple["Hypergraph", list[int]]:
        """``restricted_to`` the support: drop isolated vertices."""
        return self.restricted_to(bitset_to_vertices(self.support()))

    # ------------------------------------------------------------------
    # predicates
    # ------------------------------------------------------------------
    def is_sperner(self) -> bool:
        """True iff no hyperedge contains another (an antichain)."""
        for i, a in enumerate(self.edges):
            for j, b in enumerate(self.edges):
                if i != j and (a & b) == a:
                    return False
        return True

    def contains_edge(self, e: int) -> bool:
        return e in self.edges

    def __repr__(self) -> str:
        return f"Hypergraph(n={self.n}, edges={self.to_sets(one_indexed=True)})"
