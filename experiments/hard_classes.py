"""Chapter 4's hard classes: the frequency threshold, and symmetry through the
whole recursion.

    python experiments/hard_classes.py   ->  results/hard-classes.{json,md}

Two measurements over one instance set:

1. Placement.  ``c = eps * log2 N`` against FK96's window ``[1, 2]``
   (:mod:`fka.bound`), beside FK-A and FK-B tree sizes.  ``exp`` is
   ``log(FK-A nodes)/log N``, the empirical degree of the tree in the instance
   size -- what a superpolynomial lower bound is a claim about.

2. Symmetry, exhaustively.  Every node of both trees, deduplicated on the
   compacted ``(cnf, dnf)``, classified as read-once / transitive / primitive.
   Read-once nodes are not enumerated: their group is the cotree's, so any
   transitive subgroup preserves its block system (thesis 4.4.15).  For the rest
   the enumeration is capped -- a partial subgroup with no nontrivial block
   still proves primitivity, since blocks of a group are blocks of every
   subgroup, so the cap can only under-report.
"""

from __future__ import annotations

import itertools
import json
import math
import time
from pathlib import Path

from fka.algorithm import eps_exact, fk_a
from fka.hypergraph import Hypergraph, verts
from fka.instances import list_ids, load, matching, threshold, witt_11
from fka.graphs import primal_graph
from fka.properties import is_conformal, is_read_once
from fka.selfdual import fano, is_primitive, is_transitive, pair_automorphisms
from fka.substitution import disjoint_vee, disjoint_wedge, single, substitute
from fka.bound import shortest_edge_ratio, threshold_ratio
from fka.transversal import transversal
from fkb.algorithm import fk_b

RESULTS = Path(__file__).resolve().parents[1] / "results"
#: Automorphisms enumerated per node before the order is left unreported.
CAP = 120
#: Above this many hyperedges the incidence graph is too big to enumerate in
#: reasonable time, so the node is recorded as undetermined rather than guessed.
MAX_EDGES = 100

X = Hypergraph.from_sets(4, [[0, 1], [2, 3]])
Y = Hypergraph.from_sets(4, [[0, 2], [0, 3], [1, 2], [1, 3]])


def word(w):
    atoms = {"X": X, "Y": Y}
    MBF = atoms[w[-1]]
    for ch in reversed(w[:-1]):
        MBF = substitute(atoms[ch], [MBF] * 4)
    return MBF


def swap(w):
    return w.translate(str.maketrans("XY", "YX"))


def gk(k):
    f = Hypergraph.from_sets(2, [[0], [1]])
    for _ in range(k - 1):
        f = disjoint_vee(disjoint_wedge(f, f), disjoint_wedge(f, f))
    return f


def tk(k):
    S = single()
    for _ in range(k):
        S = disjoint_wedge(disjoint_vee(S, S), disjoint_vee(S, S))
    return S


def _has_induced_p4(g) -> bool:
    """Direct induced-P4 scan, O(n^4).

    ``fka.graphs.is_cograph`` refuses graphs over 20 vertices because its
    recogniser enumerates induced subgraphs, and the deeper alternating-tree
    members are larger than that. Only P4-freeness is needed here, and checking
    the 4-subsets one at a time is cheap enough at these sizes.
    """
    n = g.n
    for quad in itertools.combinations(range(n), 4):
        deg = [sum(1 for w in quad if w != v and g.adj[v] >> w & 1) for v in quad]
        if sorted(deg) == [1, 1, 2, 2]:
            ends = [v for v, d in zip(quad, deg) if d == 1]
            if not (g.adj[ends[0]] >> ends[1]) & 1:
                return True
    return False


def read_once(MBF: Hypergraph) -> bool:
    """``is_read_once``, falling back to the direct P4 scan past its cap."""
    if not len(MBF):
        return True
    try:
        return is_read_once(MBF)
    except ValueError:
        g = primal_graph(MBF)
        return is_conformal(MBF, g) and not _has_induced_p4(g)


def classify(cnf: Hypergraph, dnf: Hypergraph):
    """(order or None if capped, transitive, primitive, read-once) on the support."""
    keep = verts(cnf.support() | dnf.support())
    if len(keep) < 2:
        return 1, False, False, True, False
    C, _ = cnf.restricted_to(keep)
    D, _ = dnf.restricted_to(keep)
    if read_once(C):
        return None, True, False, True, False   # 4.4.15: imprimitive or symmetric
    if len(C) + len(D) > MAX_EDGES:
        return None, False, False, False, True  # undetermined
    perms, capped = pair_automorphisms(C, D, cap=CAP)
    tr = is_transitive(perms, len(keep))
    return ((None if capped else len(perms)), tr,
            tr and is_primitive(perms, len(keep)), False, False)


def walk(tree, cache):
    by_depth = {}
    for node in tree.nodes.values():
        key = (node.cnf.edges, node.dnf.edges)
        if key not in cache:
            cache[key] = classify(node.cnf, node.dnf)
        _, tr, prim, ro, unknown = cache[key]
        rec = by_depth.setdefault(
            node.depth,
            {"nodes": 0, "primitive": 0, "read_once": 0, "undetermined": 0})
        rec["nodes"] += 1
        rec["primitive"] += int(prim)
        rec["read_once"] += int(ro)
        rec["undetermined"] += int(unknown)
    return by_depth


def repeat_rate(tree):
    seen, rep = set(), 0
    for node in tree.nodes.values():
        k = (node.cnf.edges, node.dnf.edges)
        if k in seen:
            rep += 1
        else:
            seen.add(k)
    return rep / max(len(tree), 1)


def cases():
    for i in list_ids():
        inst = load(i)
        yield i, inst.cnf, inst.dnf, "library"
    W = witt_11()
    F = fano()
    yield "fano-plane", F, F, "self-dual"
    yield "w11", W, W, "self-dual"
    for k in (2, 3):
        yield f"gk-f{k}", gk(k), transversal(gk(k)), "read-once"
    for k in (1, 2):
        yield f"tk-s{k}", tk(k), transversal(tk(k)), "read-once"
    for w in ("XY", "YX"):
        yield f"word-{w.lower()}", word(w), word(swap(w)), "read-once"
    yield "matching-12", matching(12), transversal(matching(12)), "output"
    yield "threshold-10", threshold(10), transversal(threshold(10)), "shape"


def main():
    rows = []
    for name, cnf, dnf, cls in cases():
        t0 = time.time()
        ta, tb = fk_a(cnf, dnf, variant="modified"), fk_b(cnf, dnf, variant="faithful")
        cache = {}
        da, db = walk(ta, cache), walk(tb, cache)
        prim = [d for d, r in list(da.items()) + list(db.items()) if r["primitive"]]
        N = len(cnf) + len(dnf)
        top = ta.nodes[min(ta.nodes)]
        root = cache[(top.cnf.edges, top.dnf.edges)]
        rows.append({
            "instance": name, "class": cls, "n_vertices": max(cnf.n, dnf.n), "N": N,
            "eps": str(eps_exact(cnf, dnf)), "c": round(threshold_ratio(cnf, dnf), 4),
            "s_over_log2N": round(shortest_edge_ratio(cnf, dnf), 4),
            "fka_nodes": len(ta), "fka_depth": ta.depth,
            "fkb_nodes": len(tb), "fkb_depth": tb.depth,
            "ratio": round(len(ta) / len(tb), 3),
            "exponent": round(math.log(len(ta)) / math.log(N), 3),
            "fka_repeat_rate": round(repeat_rate(ta), 4),
            "distinct_pairs": len(cache),
            "root_primitive": root[2], "root_read_once": root[3],
            "last_primitive_depth": max(prim) if prim else -1,
            "primitive_nodes_fka": sum(r["primitive"] for r in da.values()),
            "primitive_nodes_fkb": sum(r["primitive"] for r in db.values()),
            "read_once_nodes_fka": sum(r["read_once"] for r in da.values()),
            "undetermined_nodes": (sum(r["undetermined"] for r in da.values())
                                   + sum(r["undetermined"] for r in db.values())),
            "fka_by_depth": {str(d): da[d] for d in sorted(da)},
            "fkb_by_depth": {str(d): db[d] for d in sorted(db)},
            "seconds": round(time.time() - t0, 1),
        })
        print(f"  {name:<16} c={rows[-1]['c']:>6.2f} FK-A={len(ta):>7} "
              f"FK-B={len(tb):>6} {rows[-1]['ratio']:>7.2f}x  "
              f"prim<=d{rows[-1]['last_primitive_depth']:<3} "
              f"[{rows[-1]['seconds']:.0f}s]", flush=True)

    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "hard-classes.json").write_text(
        json.dumps(rows, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# Hard classes: placement and symmetry depth",
        "",
        "Regenerate with `python experiments/hard_classes.py`.",
        "",
        "`c = eps * log2 N` (FK96 window is `[1, 2]`); `exp = log(FK-A nodes)/log N`;",
        "`prim depth` is the deepest recursion level still carrying a primitive,",
        "non-symmetric automorphism group.",
        "",
        "| instance | class | \\|V\\| | N | eps | c | s/log2N | FK-A | FK-B | A/B | exp | repeat | prim depth | RO root | undet |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :-: | ---: |",
    ]
    for r in rows:
        lines.append(
            f"| `{r['instance']}` | {r['class']} | {r['n_vertices']} | {r['N']} | "
            f"{r['eps']} | {r['c']:.2f} | {r['s_over_log2N']:.2f} | {r['fka_nodes']} | "
            f"{r['fkb_nodes']} | {r['ratio']:.2f}x | {r['exponent']:.2f} | "
            f"{r['fka_repeat_rate']:.1%} | {r['last_primitive_depth']} | "
            f"{'yes' if r['root_read_once'] else 'no'} | "
            f"{r['undetermined_nodes']} |")
    (RESULTS / "hard-classes.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nwrote {RESULTS/'hard-classes.json'} and {RESULTS/'hard-classes.md'}")


if __name__ == "__main__":
    main()
