"""The monotone Boolean function, stored as bitsets.

A term is the integer whose bit ``i`` is set iff variable ``i`` is in it.
Variables are 0-indexed; the thesis numbers from 1, so human-facing paths
convert with ``+1``.

Bitsets make subset, superset and disjointness single machine words --
``a & b == a`` is ``a subseteq b`` -- and terms hashable, so deduplication is a
``set`` pass. Incidence matrices remain the interchange format for SageMath and
MATLAB.

``phi_x_0``, ``phi_x_1``, ``A_c_x`` and ``A_m_x`` are the restrictions both
algorithms split on; see ``docs/dictionary.md`` for the hypergraph reading.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Iterator, Sequence

__all__ = [
    "Hypergraph",
    "verts",
    "bits",
    "popcount",
    "phi_x_0",
    "phi_x_1",
    "A_c_x",
    "A_m_x",
]


def popcount(x: int) -> int:
    """``|t|``, the number of variables in term ``x``."""
    return x.bit_count()


def bits(variables: Iterable[int]) -> int:
    """Pack 0-indexed variable ids into a bitset."""
    t = 0
    for v in variables:
        t |= 1 << v
    return t


def verts(t: int) -> tuple[int, ...]:
    """Unpack a bitset into ascending 0-indexed variable ids."""
    out = []
    while t:
        low = t & -t
        out.append(low.bit_length() - 1)
        t ^= low
    return tuple(out)


@dataclass(frozen=True, slots=True)
class Hypergraph:
    """A monotone Boolean function on variables ``0 .. n-1``, as its term set.

    ``edges`` is deduplicated and sorted, so equal term sets on equal ground
    sets compare and hash equal -- which is how repeated subproblems are
    detected (thesis p.53).

    ``n`` is carried explicitly rather than inferred: ``phi_x_0`` drops
    variables from terms, so a subproblem deep in the recursion can have
    variables in no term at all, and ``n`` keeps every node's coordinates
    aligned with the root.
    """

    n: int
    edges: tuple[int, ...]

    def __post_init__(self) -> None:
        if self.n < 0:
            raise ValueError(f"variable count must be non-negative, got {self.n}")
        mask = (1 << self.n) - 1
        for t in self.edges:
            if t & ~mask:
                raise ValueError(
                    f"term {bin(t)} references a variable outside 0..{self.n - 1}"
                )
        # Direct construction is public and useful for the degenerate ``{{}}``,
        # so the canonical invariant is enforced here too, not only in
        # from_bitsets.
        canonical = tuple(sorted(set(self.edges), key=verts))
        if canonical != self.edges:
            object.__setattr__(self, "edges", canonical)

    # ------------------------------------------------------------------
    # construction
    # ------------------------------------------------------------------
    @classmethod
    def from_bitsets(cls, n: int, edges: Iterable[int]) -> "Hypergraph":
        """Build from raw bitsets, deduplicating and sorting lexicographically.

        Sorted by the unpacked variable tuple, not numerically: ``{4,5}`` is a
        smaller integer than ``{1,5,6}`` but comes later in the notation the
        thesis and every export use.
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
        """Build from variable lists. ``one_indexed`` accepts the thesis' ``v1..vn``."""
        shift = 1 if one_indexed else 0
        packed = []
        for term in edges:
            for v in term:
                if not 0 <= v - shift < n:
                    raise ValueError(
                        f"variable {v} out of range for {n} variables "
                        f"({'1' if one_indexed else '0'}-indexed input)"
                    )
            packed.append(bits(v - shift for v in term))
        return cls.from_bitsets(n, packed)

    @classmethod
    def from_incidence(cls, matrix: Sequence[Sequence[int]]) -> "Hypergraph":
        """Build from a 0/1 incidence matrix, one row per term.

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
        """No terms: the constant 1 read as a CNF, the constant 0 read as a DNF."""
        return cls(n, ())

    # ------------------------------------------------------------------
    # export
    # ------------------------------------------------------------------
    def to_incidence(self) -> list[list[int]]:
        """A 0/1 incidence matrix, as a list of rows."""
        return [[(t >> j) & 1 for j in range(self.n)] for t in self.edges]

    def to_sets(self, *, one_indexed: bool = False) -> list[list[int]]:
        """Variable lists, ascending within each term."""
        shift = 1 if one_indexed else 0
        return [[v + shift for v in verts(t)] for t in self.edges]

    def label(self, *, one_indexed: bool = True) -> str:
        """Compact human form, e.g. ``{{1,3},{2,4}}``."""
        if not self.edges:
            return "{}"
        return (
            "{"
            + ",".join(
                "{" + ",".join(str(v) for v in term) + "}"
                for term in self.to_sets(one_indexed=one_indexed)
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
        return bool(self.edges)

    @property
    def mask(self) -> int:
        """Bitset of the whole ground set ``0 .. n-1``."""
        return (1 << self.n) - 1

    def support(self) -> int:
        """Variables in at least one term. A dual pair must have equal supports."""
        u = 0
        for t in self.edges:
            u |= t
        return u

    def degree(self, v: int) -> int:
        """Number of terms containing variable ``v``."""
        bit = 1 << v
        return sum(1 for t in self.edges if t & bit)

    def degrees(self) -> list[int]:
        """Column sums of the incidence matrix."""
        out = [0] * self.n
        for t in self.edges:
            for v in verts(t):
                out[v] += 1
        return out

    def edge_sizes(self) -> list[int]:
        """Row sums of the incidence matrix."""
        return [popcount(t) for t in self.edges]

    def isolated_vertices(self) -> tuple[int, ...]:
        """Variables in no term.

        Not allowed at the root, but ``phi_x_0`` creates them mid-recursion
        (thesis p.53), so nodes report them rather than silently renumbering.
        """
        sup = self.support()
        return tuple(v for v in range(self.n) if not (sup >> v) & 1)

    def vee(self, other: "Hypergraph") -> "Hypergraph":
        """``f v g``: term-set union, FK-A's ``[phi_x_0; phi_x_1]``.

        Deliberately does *not* reduce: both algorithms call
        :func:`fka.irredundant.irredundant` at the top of the next call, and
        keeping the two apart is what lets a node report what it discarded.
        """
        if self.n != other.n:
            raise ValueError(
                f"cannot combine functions on {self.n} and {other.n} variables"
            )
        return Hypergraph.from_bitsets(self.n, self.edges + other.edges)

    def restricted_to(self, variables: Sequence[int]) -> tuple["Hypergraph", list[int]]:
        """Re-index onto ``variables``, returning the map new -> old.

        For handing a subproblem to a backend (SageMath, GAP) that would
        otherwise treat isolated variables as real points and inflate the group.
        """
        remap = {old: new for new, old in enumerate(variables)}
        packed = [
            bits(remap[v] for v in verts(t) if v in remap) for t in self.edges
        ]
        return Hypergraph.from_bitsets(len(variables), packed), list(variables)

    def compact(self) -> tuple["Hypergraph", list[int]]:
        """``restricted_to`` the support: drop isolated variables."""
        return self.restricted_to(verts(self.support()))

    # ------------------------------------------------------------------
    # predicates
    # ------------------------------------------------------------------
    def is_sperner(self) -> bool:
        """True iff the term set is an antichain, i.e. irredundant."""
        for i, a in enumerate(self.edges):
            for j, b in enumerate(self.edges):
                if i != j and (a & b) == a:
                    return False
        return True

    def contains_edge(self, t: int) -> bool:
        return t in self.edges

    def __repr__(self) -> str:
        return f"Hypergraph(n={self.n}, edges={self.to_sets(one_indexed=True)})"


# ----------------------------------------------------------------------
# restrictions (phi_x_0.m, phi_x_1.m, A_c_x.m, A_m_x.m)
# ----------------------------------------------------------------------
def phi_x_1(MBF: Hypergraph, x: int) -> Hypergraph:
    """``phi_{x,1}``: the terms of ``MBF`` that do not contain ``x``."""
    bit = 1 << x
    return Hypergraph.from_bitsets(MBF.n, (t for t in MBF.edges if not t & bit))


def phi_x_0(MBF: Hypergraph, x: int) -> Hypergraph:
    """``phi_{x,0}``: the terms of ``MBF`` containing ``x``, with ``x`` deleted.

    A term ``{x}`` becomes the empty term and is *kept*. ``phi_x_0.m`` drops it,
    which loses the fact that the side became a constant -- reachable whenever a
    singleton term is split on.
    """
    bit = 1 << x
    return Hypergraph.from_bitsets(MBF.n, (t & ~bit for t in MBF.edges if t & bit))


def A_c_x(MBF: Hypergraph, c: int, *, clauses: bool) -> Hypergraph:
    """Set the variables of ``c`` FALSE (A_c_x.m).

    A CNF loses them from every clause; a DNF loses the monomials meeting them.
    ``clauses`` says which ``MBF`` is. The MATLAB returns the scalar ``false``
    when a side collapses; here the constants are ordinary term sets and
    :func:`fkb.algorithm.check_conditions` reads them off.
    """
    if clauses:
        return Hypergraph.from_bitsets(MBF.n, (t & ~c for t in MBF.edges))
    return Hypergraph.from_bitsets(MBF.n, (t for t in MBF.edges if not t & c))


def A_m_x(MBF: Hypergraph, m: int, *, clauses: bool) -> Hypergraph:
    """Set the variables of ``m`` TRUE (A_m_x.m). The mirror of :func:`A_c_x`."""
    if clauses:
        return Hypergraph.from_bitsets(MBF.n, (t for t in MBF.edges if not t & m))
    return Hypergraph.from_bitsets(MBF.n, (t & ~m for t in MBF.edges))
