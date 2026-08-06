"""Self-transversal hypergraphs beyond the Fano plane, by automorphism orbits.

    python experiments/selfdual_orbits.py   ->  results/selfdual-orbits.md

`MBF = Tr(MBF)` iff its term set is intersecting and not 2-colourable, decided by one
`2^n` scan (:func:`fka.selfdual.is_self_dual`). Candidates are unions of orbits
of a transitive group, so the group sits inside `Aut` by construction and
regularity is automatic; `Aut` is then computed exactly and tested for
primitivity, which certifies non-read-once by itself.

Pruning is from Prop. 2.2.4 specialised to `MBF = Tr(MBF)`: an r-uniform self-dual
family needs `|E| >= 2^(r-1)`.

The one find worth the name is `S(4,5,11)` on 11 points -- 66 blocks of size 5,
`|Aut| = 7920 = |M_11|`, primitive. On 9 points every self-dual hypergraph found
is imprimitive, the largest being `S_3 wr S_3` of order 1296, which is
majority-of-majorities and decomposes.
"""

from __future__ import annotations

import itertools
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fka.hypergraph import Hypergraph  # noqa: E402
from fka.properties import is_read_once  # noqa: E402
from fka.selfdual import (automorphisms, is_intersecting, is_primitive,  # noqa: E402
                          is_self_dual, is_transitive, subset_orbits)

RESULTS = Path(__file__).resolve().parents[1] / "results"

#: Transitive groups to build orbits under, as point permutations.
GROUPS = {
    9: {"Z9": [tuple((v + 1) % 9 for v in range(9))],
        "Z3xZ3": [tuple((v // 3) * 3 + (v + 1) % 3 for v in range(9)),
                  tuple((v + 3) % 9 for v in range(9))],
        "AGL(1,9)": [tuple((v + 1) % 9 for v in range(9)),
                     tuple((2 * v) % 9 for v in range(9))]},
    11: {"Z11": [tuple((v + 1) % 11 for v in range(11))],
         "AGL(1,11)": [tuple((v + 1) % 11 for v in range(11)),
                       tuple((3 * v) % 11 for v in range(11))]},
    13: {"Z13": [tuple((v + 1) % 13 for v in range(13))],
         "F39": [tuple((v + 1) % 13 for v in range(13)),
                 tuple((3 * v) % 13 for v in range(13))]},
}


def main():
    out, hits = [], {}
    log = lambda s: (print(s, flush=True), out.append(s))  # noqa: E731
    log("# Self-transversal hypergraphs by orbit search")
    log("")
    log("Regenerate with `python experiments/selfdual_orbits.py`.")
    log("")
    log("Self-dual = intersecting and not 2-colourable, decided by the 2^n scan.")
    log("")
    t0 = time.time()
    for n, groups in GROUPS.items():
        for gname, gens in groups.items():
            for k in range(3, n // 2 + 2):
                orbs = subset_orbits(gens, n, k)
                need = 2 ** (k - 1)
                for count in (1, 2, 3):
                    for combo in itertools.combinations(orbs, count):
                        E = [e for o in combo for e in o]
                        if len(E) < need:
                            continue
                        MBF = Hypergraph.from_bitsets(n, E)
                        if not is_intersecting(MBF) or not is_self_dual(MBF):
                            continue
                        hits.setdefault((n, tuple(sorted(E))), (gname, k))
            log(f"- scanned n = {n}, {gname}")

    log("")
    log("| n | \\|E\\| | size | degrees | \\|Aut\\| | transitive | primitive | read-once | found under |")
    log("| ---: | ---: | ---: | ---: | ---: | :-: | :-: | :-: | --- |")
    best = None
    for (n, E), (gname, k) in sorted(hits.items(), key=lambda kv: (kv[0][0], len(kv[0][1]))):
        MBF = Hypergraph.from_bitsets(n, list(E))
        perms, capped = automorphisms(MBF, cap=20000)
        tr = is_transitive(perms, n)
        pr = tr and is_primitive(perms, n)
        order = f"> {len(perms)}" if capped else str(len(perms))
        log(f"| {n} | {len(MBF)} | {sorted(set(MBF.edge_sizes()))} | "
            f"{sorted(set(MBF.degrees()))} | {order} | {'y' if tr else 'n'} | "
            f"{'y' if pr else 'n'} | {'y' if is_read_once(MBF) else 'n'} | {gname} |")
        if pr and not capped and (best is None or len(perms) > best[1]):
            best = (MBF, len(perms), n)

    log("")
    if best is not None:
        MBF, order, n = best
        log(f"Largest primitive group found: order {order} on {n} points.")
        log("")
        log(f"    {MBF.label()}")
    log("")
    log(f"_{time.time() - t0:.0f}s_")
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "selfdual-orbits.md").write_text("\n".join(out) + "\n", encoding="utf-8")
    print(f"\nwrote {RESULTS / 'selfdual-orbits.md'}")


if __name__ == "__main__":
    main()
