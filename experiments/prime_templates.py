"""Is any template at the FK-A frequency threshold *not* read-once?

    python experiments/prime_templates.py   ->  results/prime-templates.md

A template ``T`` generates ``Gamma_0 = {{v}}``, ``Gamma_{j+1} = T[Gamma_j,...]``.
From the recursions in :mod:`fka.bound`, the threshold ratio converges only
when ``eps_T * r_T = 1``; since ``eps_T >= r_T / t`` (average degree), that forces
``t = r^2`` with ``T`` r-uniform *and regular*, and the same of ``Tr(T)``. Then

    c_inf = max(log2 |E(T)|, log2 |E(Tr T)|) / (r - 1).

Counting also shows every block of ``T`` meets every block of ``Tr(T)`` in
exactly one point, so ``Tr(T)`` is the set of *perfect selectors*: r points
meeting every block once. That makes the whole search Berge-free -- selectors
are read off directly, and duality is decided by the ``2^t`` subset scan.

Normalisation. Any valid ``T`` has a perfect selector, so relabel it to
``{0,..,r-1}``; then every block carries exactly one of those points plus an
``(r-1)``-set from the remaining ``r^2 - r``, which is the parametrisation
enumerated here.

Result: for ``r = 2`` and ``r = 3`` every convergent template is read-once, so
no substitution-closed family at the threshold escapes the alternating trees.
``r >= 4`` is open.
"""

from __future__ import annotations

import itertools
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fka.graphs import is_cograph, primal_graph  # noqa: E402
from fka.hypergraph import Hypergraph, bits  # noqa: E402
from fka.properties import is_conformal, is_read_once  # noqa: E402
from fka.bound import log2_int  # noqa: E402

RESULTS = Path(__file__).resolve().parents[1] / "results"


def selectors(blocks, t, r):
    """r-subsets meeting every block in exactly one point."""
    return [bits(s) for s in itertools.combinations(range(t), r)
            if all((bits(s) & b).bit_count() == 1 for b in blocks)]


def is_dual(E, F, t):
    """No S with S containing no E-block and its complement no F-block."""
    full = (1 << t) - 1
    return not any(
        not any((b & S) == b for b in E) and not any((b & (full ^ S)) == b for b in F)
        for S in range(1 << t))


def search(r, dmax, log):
    t = r * r
    rest = list(range(r, t))
    parts = [tuple(c) for c in itertools.combinations(rest, r - 1)]
    by_deg = defaultdict(list)
    log(f"\n## r = {r}, t = {t}\n")
    log(f"c_inf <= 2 needs |E| <= {4 ** (r - 1)}, and |E| = {r}d.\n")
    log("| d | \\|E\\| | groups tried | sel(E) nonempty | closed + regular | DUAL | not read-once |")
    log("| ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    found = []
    for d in range(1, dmax + 1):
        groups = list(itertools.combinations(parts, d))
        by_deg.clear()
        for g in groups:
            dd = [0] * t
            for p in g:
                for v in p:
                    dd[v] += 1
            by_deg[tuple(dd[r:])].append(g)
        tried = nsel = nclosed = ndual = nprime = 0
        target = (d,) * (t - r)
        for combo in itertools.product(groups, repeat=r - 1):
            deg = [0] * (t - r)
            ok = True
            for g in combo:
                for p in g:
                    for v in p:
                        deg[v - r] += 1
                        if deg[v - r] > d:
                            ok = False
            if not ok:
                continue
            need = tuple(d - x for x in deg)
            for last in by_deg.get(need, ()):
                tried += 1
                E = [bits((i,) + p) for i, g in enumerate(list(combo) + [last])
                     for p in g]
                if len(set(E)) != r * d:
                    continue
                F = selectors(E, t, r)
                if not F:
                    continue
                nsel += 1
                degF = [sum(1 for b in F if b >> v & 1) for v in range(t)]
                if len(set(degF)) != 1 or set(selectors(F, t, r)) != set(E):
                    continue
                nclosed += 1
                if not is_dual(E, F, t):
                    continue
                ndual += 1
                MBF = Hypergraph.from_bitsets(t, E)
                g = primal_graph(MBF)
                ro = is_read_once(MBF, g)
                nprime += int(not ro)
                found.append((d, MBF, Hypergraph.from_bitsets(t, F), ro,
                              is_conformal(MBF, g), is_cograph(g)))
        log(f"| {d} | {r*d} | {tried} | {nsel} | {nclosed} | {ndual} | {nprime} |")
    return found


def main():
    out = []
    log = lambda s: (print(s, flush=True), out.append(s))  # noqa: E731
    log("# Convergent templates at the FK-A frequency threshold")
    log("")
    log("Regenerate with `python experiments/prime_templates.py`.")
    log("")
    log("A template is r-uniform and regular on t = r^2 vertices with an")
    log("r-uniform regular dual; c_inf = max(log2|E|, log2|Tr E|)/(r-1).")
    t0 = time.time()
    every = []
    for r, dmax in ((2, 3), (3, 5)):
        every += [(r, *f) for f in search(r, dmax, log)]
    log("\n## Templates found\n")
    log("| r | \\|E\\| | \\|Tr(E)\\| | c_inf | read-once | E |")
    log("| ---: | ---: | ---: | ---: | :-: | --- |")
    seen = set()
    for r, d, MBF, F, ro, conf, cog in every:
        key = (r, len(MBF), len(F), MBF.label())
        if key in seen:
            continue
        seen.add(key)
        c = max(log2_int(len(MBF)), log2_int(len(F))) / (r - 1)
        log(f"| {r} | {len(MBF)} | {len(F)} | {c:.4f} | "
            f"{'yes' if ro else 'NO'} | `{MBF.label()}` |")
    log("")
    if any(not f[4] for f in every):
        log("**A non-read-once template exists.**")
    else:
        log("Every convergent template for r = 2 and r = 3 is read-once: the")
        log("substitution route to the threshold cannot leave the alternating")
        log("trees. r >= 4 is not settled here.")
    log(f"\n_{time.time() - t0:.0f}s_")
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "prime-templates.md").write_text("\n".join(out) + "\n", encoding="utf-8")
    print(f"\nwrote {RESULTS / 'prime-templates.md'}")


if __name__ == "__main__":
    main()
