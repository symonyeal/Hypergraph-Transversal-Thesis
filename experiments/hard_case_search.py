"""Search the vertex-transitive self-transversal hypergraphs for hard cases.

    sage -python experiments/hard_case_search.py            degrees 3..13
    sage -python experiments/hard_case_search.py 3 15        a chosen range

Writes ``results/hard-case-search.json``. Requires SageMath: the enumeration is
over SageMath's transitive-group database, and the identifications
(automorphism group, primitivity, canonical labelling) are Sage's.

Why this class, and why it is a search rather than a brute force
----------------------------------------------------------------
The thesis' two hard classes are the Fano plane (Section 5.2, self-transversal,
``Aut = PSL(3,2)``) and the alternating trees (Section 4.2, ``Tr(f_k) = g_k``,
frequency exactly at the Fredman-Khachiyan floor). Both are *self-transversal*
or a tight dual pair, which is the condition that leaves FK no excuse: the
answer is "dual", so every branch has to be closed, and ``|Tr(H)| = |H|`` means
no part of the tree can be attributed to an exponential output.

A hypergraph is self-transversal exactly when it is a **maximal intersecting
family**: for every ``S``, exactly one of ``S``, ``V\\S`` contains an edge.
Three facts turn that class into a finite search.

**L1.** Fix an edge ``e`` of minimum size ``s``. Every other edge meets ``e``, so
with ``D`` the maximum degree, ``m <= 1 + s(D-1)``, hence

    eps(H) = D/m > 1/s.

**L2.** FK96 guarantees ``eps >= 1/log2(|G|+|H|)`` on a dual pair, and thesis
Section 4.3 calls a family hard when ``eps`` sits on that floor. Writing

    T(H) = eps * log2(2m)   >= 1,

L1 gives ``T > log2(2m)/s``: a tight instance needs edges of size about the
log of their own number. The Fano plane has ``s = 3`` against ``log2(14) =
3.81``, so ``T = 1.632`` -- the tightest instance in the library.

**L3.** If ``H`` is ``s``-uniform and vertex-transitive then ``D = ms/n``, and
L1 collapses to ``n <= s^2 - 1``. So for each ground set only a few edge sizes
are possible at all, and every candidate is a union of orbits of a transitive
group acting on subsets.

The search enumerates exactly that: every transitive group of each degree,
every orbit of ``s``-subsets allowed by L3, and every intersecting *pair* of
orbits (the non-uniform candidates). Candidates are filtered by two cheap
necessary conditions before the exact test, which decides ``Tr(H) = H`` over
all ``2^n`` subsets at once rather than by building ``Tr(H)``.

What it finds
-------------
See ``docs/hard-cases.md``. Up to 15 points the class splits in two: the
instances whose automorphism group is *imprimitive* are compositions
(``MAJ3[MAJ3]``, ``MAJ5[MAJ3]``, ``MAJ3[MAJ5]``), hard for the same reason the
alternating trees are -- FK does not notice that they decompose. Among the
*primitive* ones, the Fano plane ``S(2,3,7)`` and the Witt design
``S(4,5,11) = W11`` are the two non-trivial families, and W11 is the harder.
"""

from __future__ import annotations

import json
import math
import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

try:
    import sage.all  # noqa: F401
    from sage.all import TransitiveGroups
    from sage.combinat.designs.incidence_structures import IncidenceStructure
except ImportError:  # pragma: no cover - the program is Sage-only by design
    raise SystemExit(
        "hard_case_search.py needs SageMath: it enumerates SageMath's "
        "transitive-group database. Run it under `sage -python`, or see "
        "docs/sagemath.md."
    )

import numpy as np  # noqa: E402

from fka.algorithm import fk_a  # noqa: E402
from fka.hypergraph import Hypergraph, bits, verts  # noqa: E402
from fka.properties import analyse  # noqa: E402
from fka.transversal import is_dual_oracle  # noqa: E402
from fkb.algorithm import fk_b  # noqa: E402

RESULTS = Path(__file__).resolve().parents[1] / "results"

#: Above this many edges the two algorithms are not run and the edge list is
#: not recorded: the point of the search is the class, and the majority families
#: (``C(2k+1, k+1)`` edges) are known to be loose -- ``eps`` is about 1/2, so
#: they are nowhere near the FK floor whatever their tree looks like, and their
#: blocks follow from ``n``.
RUN_BELOW = 400


# ----------------------------------------------------------------------
# group orbits
# ----------------------------------------------------------------------
def permute(e: int, p: list[int]) -> int:
    """Image of the bitset ``e`` under the 0-based permutation ``p``."""
    out = 0
    while e:
        low = e & -e
        out |= 1 << p[low.bit_length() - 1]
        e ^= low
    return out


def orbits_on_subsets(gens: list[list[int]], n: int, k: int) -> list[frozenset[int]]:
    """Every orbit of ``<gens>`` on the ``k``-subsets of ``0..n-1``."""
    seen: set[int] = set()
    out: list[frozenset[int]] = []
    for c in combinations(range(n), k):
        e = bits(c)
        if e in seen:
            continue
        orbit = {e}
        frontier = [e]
        while frontier:
            x = frontier.pop()
            for p in gens:
                y = permute(x, p)
                if y not in orbit:
                    orbit.add(y)
                    frontier.append(y)
        seen |= orbit
        out.append(frozenset(orbit))
    return out


# ----------------------------------------------------------------------
# the three tests, cheapest first
# ----------------------------------------------------------------------
def intersecting(edges) -> bool:
    es = list(edges)
    return all(a & b for i, a in enumerate(es) for b in es[i + 1:])


def cross_intersecting(a, b) -> bool:
    return all(x & y for x in a for y in b)


def edges_are_minimal_transversals(edges) -> bool:
    """Necessary for ``Tr(H) = H``: every edge must stay minimal as a
    transversal, so for each edge ``e`` and each ``x`` in ``e`` some edge must
    meet ``e`` in ``x`` alone. Rejects almost every intersecting orbit."""
    es = list(edges)
    for e in es:
        x = e
        while x:
            low = x & -x
            x ^= low
            rest = e & ~low
            if not any((f & rest) == 0 and (f & low) for f in es):
                return False
    return True


def is_self_transversal(H: Hypergraph) -> bool:
    """Exact. ``Tr(H) = H`` iff no ``S`` has neither ``S`` nor ``V\\S``
    containing an edge -- that ``S`` is precisely a thesis equation 2.1
    certificate against ``(H, H)``. Tested over all ``2^n`` subsets at once;
    complementation is index reversal, since ``V\\S`` is ``FULL - S``."""
    if H.n > 22:
        return is_dual_oracle(H, H)
    masks = np.arange(1 << H.n, dtype=np.int64)
    holds = np.zeros(1 << H.n, dtype=bool)
    for e in H.edges:
        holds |= (masks & e) == e
    return not bool((~holds & ~holds[::-1]).any())


# ----------------------------------------------------------------------
def structure(H: Hypergraph) -> IncidenceStructure:
    target, _ = H.compact()
    return IncidenceStructure(
        list(range(target.n)), [list(verts(e)) for e in target.edges]
    )


def canonical_key(H: Hypergraph):
    """Isomorphism-invariant key, from Sage's canonical labelling."""
    I = structure(H)
    C = I.relabel(I.canonical_label(), inplace=False)
    return (H.n, tuple(sorted(tuple(sorted(b)) for b in C.blocks())))


def describe(H: Hypergraph) -> dict:
    """Everything Sage can say about one hit, plus both algorithms' trees."""
    I = structure(H)
    A = I.automorphism_group()
    m = len(H)
    eps = max(H.degrees()) / m
    row = {
        "n": H.n,
        "m": m,
        "edge_sizes": sorted(set(H.edge_sizes())),
        "degrees": sorted(set(H.degrees())),
        "epsilon": eps,
        "T": eps * math.log2(2 * m),
        "aut_order": int(A.order()),
        "aut_name": str(A.structure_description()),
        "primitive": bool(A.is_primitive()),
        "t_design": list(I.is_t_design(return_parameters=True)[1]),
        "edges": H.to_sets(one_indexed=True),
    }
    props = analyse(H, graph_classes=False)
    row["conformal"] = props.conformal
    row["read_once"] = props.read_once
    if m >= RUN_BELOW:
        # Every family this large in the searched range is a majority function,
        # whose blocks are reconstructible from n alone. Recording 6435 blocks
        # of size 8 costs 235 KB and says nothing the parameters do not, so the
        # construction is recorded instead. Cf. the fk-b.multiple withdrawal in
        # _archive/20260804-redundant-fkb-multiple-artifacts/.
        half = (H.n + 1) // 2
        is_majority = row["edge_sizes"] == [half] and m == math.comb(H.n, half)
        row["edges"] = None
        row["edges_omitted"] = (
            f"majority({H.n}): every {half}-subset of {H.n} points"
            if is_majority
            else f"{m} edges, above the {RUN_BELOW}-edge recording threshold"
        )
    if m < RUN_BELOW:
        a = fk_a(H, H, variant="modified")
        b = fk_b(H, H, variant="faithful")
        if not (a.verdict.dual and b.verdict.dual):
            raise AssertionError("an algorithm disagreed with the exact test")
        row["fka_nodes"] = len(a)
        row["fkb_nodes"] = len(b)
        row["fka_per_edge"] = len(a) / (2 * m)
        row["fkb_per_edge"] = len(b) / (2 * m)
    return row


def search(n_from: int, n_to: int, pair_max: int = 12) -> list[dict]:
    found: dict = {}
    for n in range(n_from, n_to + 1):
        groups = list(TransitiveGroups(n))
        print(f"degree {n}: {len(groups)} transitive groups", flush=True)
        for T in groups:
            gens = [[int(g(i + 1)) - 1 for i in range(n)] for g in T.gens()]
            if not gens:
                gens = [list(range(n))]
            label = f"{n}T{T.transitive_number()}"
            intersecting_orbits = []
            for k in range(2, (n + 3) // 2 + 1):
                if n > k * k - 1:          # L3
                    continue
                for orbit in orbits_on_subsets(gens, n, k):
                    if not intersecting(orbit):
                        continue
                    intersecting_orbits.append(orbit)
                    if not edges_are_minimal_transversals(orbit):
                        continue
                    H = Hypergraph.from_bitsets(n, orbit)
                    if is_self_transversal(H):
                        found.setdefault(canonical_key(H), (H, set()))[1].add(label)
            if n > pair_max:
                continue
            for a, b in combinations(intersecting_orbits, 2):
                if not cross_intersecting(a, b):
                    continue
                edges = list(a | b)
                if any(x != y and (x & y) == x for x in edges for y in edges):
                    continue               # a self-transversal family is Sperner
                if not edges_are_minimal_transversals(edges):
                    continue
                H = Hypergraph.from_bitsets(n, edges)
                if is_self_transversal(H):
                    found.setdefault(canonical_key(H), (H, set()))[1].add(label)

    rows = []
    for H, labels in found.values():
        row = describe(H)
        row["groups"] = sorted(labels)
        rows.append(row)
    return sorted(rows, key=lambda r: (-r.get("fka_per_edge", 0), r["n"]))


def main(argv: list[str]) -> int:
    n_from = int(argv[1]) if len(argv) > 1 else 3
    n_to = int(argv[2]) if len(argv) > 2 else 13
    rows = search(n_from, n_to)

    print(f"\n{len(rows)} vertex-transitive self-transversal hypergraphs, "
          f"degrees {n_from}..{n_to}\n")
    head = (f"{'n':>3}{'m':>6}{'sizes':>9}{'eps':>8}{'T':>7}{'FK-A':>7}{'FK-B':>7}"
            f"{'A/2m':>7}{'B/2m':>7}{'prim':>6}{'|Aut|':>10}  Aut")
    print(head)
    print("-" * len(head))
    for r in rows:
        a, b = r.get("fka_nodes"), r.get("fkb_nodes")
        print(f"{r['n']:>3}{r['m']:>6}{str(r['edge_sizes']):>9}{r['epsilon']:>8.4f}"
              f"{r['T']:>7.3f}{str(a or '-'):>7}{str(b or '-'):>7}"
              f"{r.get('fka_per_edge', float('nan')):>7.2f}"
              f"{r.get('fkb_per_edge', float('nan')):>7.2f}"
              f"{str(r['primitive']):>6}{r['aut_order']:>10}  {r['aut_name'][:34]}")

    RESULTS.mkdir(parents=True, exist_ok=True)
    out = RESULTS / "hard-case-search.json"
    out.write_text(
        json.dumps({"degrees": [n_from, n_to], "families": rows}, indent=1) + "\n",
        encoding="utf-8",
    )
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
