"""FK-A against FK-B: tree size, and what happens to symmetry inside each tree.

    python experiments/compare_algorithms.py

Writes ``results/fk-a-vs-fk-b.json`` and ``results/fk-a-vs-fk-b.md``.

Chapter 5 of the thesis studies FK-A by annotating each recursion-tree node with
the automorphism group of its two subhypergraphs. The same annotation applies to
an FK-B tree without modification -- both are :class:`fka.tree.RecursionTree` --
so the two decompositions can be compared node for node on the same instance,
which is what this program does.

Three things are collected per instance and per algorithm:

``shape``
    node count, depth, leaf count, and for FK-B the split-branch mix.
``groups``
    the multiset of ``Aut(G)`` names over all nodes, the number of nodes whose
    two sides disagree, and how far the group falls from the root's.
``symmetry decay``
    the root group order against the mean and minimum over the tree. The thesis'
    reading of a hard instance is one whose symmetry survives decomposition:
    the subproblems stay as symmetric as the parent, so none of them is easier
    than what it came from.

Nothing here decides duality. The oracle does that, and both algorithms are
checked against it first; if either disagrees the program stops rather than
reporting statistics about a wrong answer.
"""

from __future__ import annotations

import json
import statistics
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fka.algorithm import fk_a  # noqa: E402
from fka.analysis import AnalysisOptions, annotate  # noqa: E402
from fka.instances import load_all, self_dual_fano  # noqa: E402
from fka.transversal import is_dual_oracle  # noqa: E402
from fka.tree import RecursionTree  # noqa: E402
from fkb.algorithm import fk_b  # noqa: E402

#: Group analysis is exponential in the worst case, and the point here is the
#: shape of the tree, not a record-breaking group computation.
OPTIONS = AnalysisOptions(graph_classes=False, skip_group_above=40)

#: Above this many nodes a tree is reported by shape alone. Annotation memoises
#: on the (G, H) pair, but a tree with tens of thousands of nodes still holds
#: thousands of distinct subproblems, and at up to a second each that is hours
#: for a number the smaller instances already establish. Which instances were
#: annotated is recorded per row, so a shape-only row is never mistaken for a
#: symmetry measurement of zero.
ANNOTATE_BELOW = 1000

RESULTS = Path(__file__).resolve().parents[1] / "results"

#: The variant each algorithm's committed artifacts were written at, so an
#: already-annotated tree can be reused instead of re-derived. Annotating
#: ``sdfp-sd-2``'s FK-A tree alone costs several minutes, and ``fka run --all``
#: has already paid it.
ARTIFACT = {"FK-A": "fk-a.modified", "FK-B": "fk-b.faithful"}


def profile(tree: RecursionTree, annotated: bool) -> dict:
    """Reduce a tree to the numbers being compared.

    ``annotated=False`` reports shape only, and says so, rather than reporting
    an absent group distribution as an empty one.
    """
    orders: list[int] = []
    names: Counter[str] = Counter()
    disagree = 0
    analysed = 0
    for node in tree:
        a = node.analysis or {}
        aut_g = a.get("aut_G")
        if aut_g:
            orders.append(aut_g["order"])
            names[aut_g["name"]] += 1
            analysed += 1
        if a.get("aut_agree") is False:
            disagree += 1

    root = (tree.root.analysis or {}).get("aut_G") or {}
    root_order = root.get("order")
    out = {
        "nodes": len(tree),
        "leaves": len(tree.leaves),
        "depth": tree.depth,
        "dual": tree.dual,
        "annotated": annotated,
        "root_group": root.get("name"),
        "root_order": root_order,
        "nodes_analysed": analysed,
        "aut_disagreements": disagree,
        "group_names": dict(names.most_common()),
        "trivial_group_nodes": names.get("1", 0) + names.get("C1", 0),
    }
    if orders:
        out["order_mean"] = round(statistics.fmean(orders), 3)
        out["order_min"] = min(orders)
        out["order_max"] = max(orders)
        # How much of the root's symmetry the average subproblem still carries.
        if root_order:
            out["symmetry_retained"] = round(statistics.fmean(orders) / root_order, 4)
    branches = tree.summary().get("branches")
    if branches:
        out["branches"] = branches
    return out


def run(name: str, G, H, instance_id: str) -> dict:
    truth = is_dual_oracle(G, H)
    row: dict = {"instance": instance_id, "n": G.n, "edges_G": len(G), "edges_H": len(H),
                 "oracle_dual": truth, "algorithms": {}}
    for label, build in (
        ("FK-A", lambda: fk_a(G, H, variant="modified", instance=instance_id)),
        ("FK-B", lambda: fk_b(G, H, variant="faithful", instance=instance_id)),
    ):
        # Time the algorithm itself, always: that is a number being compared.
        t0 = time.perf_counter()
        tree = build()
        elapsed = time.perf_counter() - t0
        if tree.dual != truth:
            raise SystemExit(
                f"{label} disagrees with the oracle on {instance_id}: "
                f"{tree.dual} vs {truth}"
            )

        # Reuse the committed artifact's annotation when it is for this exact
        # tree, and fall back to computing it only for instances that have none.
        annotated = False
        artifact = RESULTS / f"{instance_id}.{ARTIFACT[label]}.json"
        if artifact.exists():
            stored = RecursionTree.load(artifact)
            if len(stored) == len(tree) and any(n.analysis for n in stored):
                tree, annotated = stored, True
        if not annotated and len(tree) < ANNOTATE_BELOW:
            annotate(tree, OPTIONS)
            annotated = True

        row["algorithms"][label] = profile(tree, annotated) | {
            "seconds": round(elapsed, 4)
        }
    a, b = row["algorithms"]["FK-A"], row["algorithms"]["FK-B"]
    row["node_ratio"] = round(a["nodes"] / b["nodes"], 3)
    print(
        f"{name:16} n={G.n:3} |G|={len(G):4} |H|={len(H):4} "
        f"fk-a={a['nodes']:6} fk-b={b['nodes']:5} ratio={row['node_ratio']:6.2f}x "
        f"root={a['root_group'] or '-'}"
        + ("" if a["annotated"] and b["annotated"] else "  (shape only)")
    )
    return row


def main() -> int:
    rows = []
    for inst in load_all():
        rows.append(run(inst.id, inst.G, inst.H, inst.id))

    # SDFP(3) is not in the instance library -- 365 edges makes every baseline
    # refresh and test run expensive -- but it is where the gap is widest, so it
    # is included in the comparison. Both its trees are past ANNOTATE_BELOW, so
    # it contributes node counts only.
    E = self_dual_fano(3)
    rows.append(run("sdfp-sd-3", E, E, "sdfp-sd-3"))

    out = RESULTS
    out.mkdir(parents=True, exist_ok=True)
    (out / "fk-a-vs-fk-b.json").write_text(
        json.dumps({"instances": rows}, indent=2), encoding="utf-8"
    )
    (out / "fk-a-vs-fk-b.md").write_text(render(rows), encoding="utf-8")
    print(f"\nwrote {out / 'fk-a-vs-fk-b.json'}")
    print(f"wrote {out / 'fk-a-vs-fk-b.md'}")
    return 0


def render(rows: list[dict]) -> str:
    lines = [
        "# FK-A against FK-B",
        "",
        "Regenerate with `python experiments/compare_algorithms.py`. Every row is",
        "checked against the brute-force transversal oracle before it is reported.",
        "",
        "`retained` is the mean `|Aut(G)|` over the tree divided by the root's: how",
        "much of the instance's symmetry the average subproblem still carries.",
        "",
        "| instance | n | \\|G\\| | \\|H\\| | dual | root Aut | FK-A nodes | FK-B nodes | ratio | FK-A retained | FK-B retained |",
        "| --- | ---: | ---: | ---: | :-: | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in rows:
        a, b = r["algorithms"]["FK-A"], r["algorithms"]["FK-B"]
        note = "" if a["annotated"] and b["annotated"] else " ‡"
        lines.append(
            f"| `{r['instance']}`{note} | {r['n']} | {r['edges_G']} | {r['edges_H']} | "
            f"{'yes' if r['oracle_dual'] else 'no'} | "
            f"{a['root_group'] or '—'} ({a['root_order'] or '—'}) | "
            f"{a['nodes']} | {b['nodes']} | {r['node_ratio']}x | "
            f"{a.get('symmetry_retained', '—')} | {b.get('symmetry_retained', '—')} |"
        )
    lines += ["", "‡ too large to annotate; node counts only. See `ANNOTATE_BELOW`."]

    lines += ["", "## Automorphism groups seen inside each tree", "",
              "Group of the `G` side, counted over every node whose group was",
              "computed. A node above `skip_group_above` edges is not counted.", "",
              "| instance | FK-A groups | FK-B groups |", "| --- | --- | --- |"]
    for r in rows:
        a, b = r["algorithms"]["FK-A"], r["algorithms"]["FK-B"]

        def fmt(p):
            if not p["annotated"]:
                return "not annotated"
            items = list(p["group_names"].items())[:4]
            return ", ".join(f"`{k}`×{v}" for k, v in items) or "—"

        lines.append(f"| `{r['instance']}` | {fmt(a)} | {fmt(b)} |")

    lines += ["", "## FK-B split branches", "",
              "`mu_D` and `mu_C` are the per-clause and per-monomial branches, taken",
              "when the split variable is at most mu-frequent on one side. `split` is",
              "the plain two-way split, which is FK-A's step.", "",
              "| instance | branches |", "| --- | --- |"]
    for r in rows:
        b = r["algorithms"]["FK-B"]
        got = b.get("branches") or {}
        lines.append(
            f"| `{r['instance']}` | "
            + (", ".join(f"`{k}`×{v}" for k, v in sorted(got.items())) or "— (no split)")
            + " |"
        )
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
