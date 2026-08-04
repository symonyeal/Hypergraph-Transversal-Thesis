"""Command-line interface: ``python -m fkb <command>``.

    python -m fkb run fano                  FK-B + annotation + HTML report
    python -m fkb run --all --all-variants  every instance, every variant
    python -m fkb verify                    check FK-B against the oracle
    python -m fkb compare --all             FK-A vs FK-B vs the oracle
    python -m fkb benchmark                 both algorithms on the scaling families
    python -m fkb dualize fano              generate Tr(H) by repeated FK-B

The instance library, the ``env``, ``list`` and ``show`` commands, and the
baseline refresh are not duplicated here -- they belong to the instances, not
to either algorithm. Use ``python -m fka`` for those.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Optional, Sequence

from fka.algorithm import fk_a
from fka.analysis import AnalysisOptions, annotate
from fka.dot import write_dot
from fka.instances import (
    Instance,
    list_ids,
    load,
    load_all,
    matching,
    self_dual_fano,
    threshold,
)
from fka.report import write_report
from fka.transversal import is_dual_oracle, transversal

from . import __version__
from .algorithm import SPLIT_RULES, VARIANTS, fk_b
from .dualize import dualise


def _results_dir(arg: Optional[str]) -> Path:
    if arg:
        return Path(arg)
    return Path(__file__).resolve().parents[2] / "results"


def _targets(args: argparse.Namespace) -> Optional[list[Instance]]:
    if args.all:
        return load_all()
    if args.instance:
        return [load(i) for i in args.instance]
    print("give an instance id, or --all. Available: " + ", ".join(list_ids()))
    return None


# ----------------------------------------------------------------------
def _run_one(inst: Instance, args: argparse.Namespace, outdir: Path) -> dict:
    tree = fk_b(
        inst.G,
        inst.H,
        variant=args.variant,
        split_rule=args.split_rule,
        instance=inst.id,
        validate=args.validate,
        max_nodes=args.max_nodes,
    )
    if not args.no_annotate:
        annotate(
            tree,
            AnalysisOptions(graph_classes=not args.no_graph_classes, backend=args.backend),
        )

    stem = f"{inst.id}.fk-b.{args.variant}"
    outputs = [
        tree.save(outdir / f"{stem}.json"),
        write_report(tree, outdir / f"{stem}.html", title=inst.name, instance=inst),
    ]
    if args.dot:
        outputs.append(write_dot(tree, outdir / f"{stem}.dot", view=args.dot_view))

    summary = tree.summary()
    flag = ""
    baseline = (inst.expected.get("fkb", {}) or {}).get(args.variant, {})
    if baseline and args.split_rule == "most_frequent":
        if baseline.get("nodes") != summary["nodes"]:
            flag = f"  !! baseline says {baseline['nodes']} nodes"
    expected_dual = inst.expected.get("dual")
    if expected_dual is not None and expected_dual != tree.dual:
        flag += f"  !! baseline says dual={expected_dual}"

    branches = summary.get("branches", {})
    print(
        f"{inst.id:16} dual={str(tree.dual):5} nodes={summary['nodes']:4} "
        f"depth={summary['depth']:2} branches={branches or '-'}{flag}"
    )
    for p in outputs:
        print(f"                 -> {p}")
    return summary


def cmd_run(args: argparse.Namespace) -> int:
    targets = _targets(args)
    if targets is None:
        return 2
    outdir = _results_dir(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    variants = VARIANTS if args.all_variants else (args.variant,)
    for inst in targets:
        for variant in variants:
            args.variant = variant
            _run_one(inst, args, outdir)
    return 0


def cmd_verify(args: argparse.Namespace) -> int:
    """Check FK-B against the brute-force oracle on every instance and setting."""
    failures = 0
    for inst in load_all():
        truth = is_dual_oracle(inst.G, inst.H)
        for variant in VARIANTS:
            for rule in SPLIT_RULES:
                tree = fk_b(
                    inst.G,
                    inst.H,
                    variant=variant,
                    split_rule=rule,
                    instance=inst.id,
                    validate=True,
                )
                ok = tree.dual == truth
                failures += not ok
                print(
                    f"{'ok  ' if ok else 'FAIL'} {inst.id:16} {variant:9} {rule:14} "
                    f"fk-b={str(tree.dual):5} oracle={str(truth):5} nodes={len(tree)}"
                )
    print(f"\n{'all agree with the oracle' if not failures else f'{failures} FAILURES'}")
    return 1 if failures else 0


def cmd_compare(args: argparse.Namespace) -> int:
    """Run FK-A, FK-B and the oracle on the same instances and tabulate them.

    Node counts are not comparable as a quality measure -- the two algorithms
    generate different subproblems, and FK-B's per-clause branches are cheaper
    per node than FK-A's. They are reported because a change in either column
    is the thing worth noticing.
    """
    targets = _targets(args)
    if targets is None:
        return 2

    head = (
        "instance",
        "oracle",
        "fk-a",
        "nodes",
        "depth",
        "fk-b",
        "nodes",
        "depth",
        "branches",
        "agree",
    )
    rows = []
    failures = 0
    for inst in targets:
        truth = is_dual_oracle(inst.G, inst.H)
        a = fk_a(inst.G, inst.H, variant=args.fka_variant, instance=inst.id, validate=True)
        b = fk_b(
            inst.G,
            inst.H,
            variant=args.variant,
            split_rule=args.split_rule,
            instance=inst.id,
            validate=True,
        )
        agree = a.dual == truth == b.dual
        failures += not agree
        sb = b.summary()
        rows.append(
            (
                inst.id,
                str(truth),
                str(a.dual),
                str(len(a)),
                str(a.depth),
                str(b.dual),
                str(len(b)),
                str(b.depth),
                ",".join(f"{k}={v}" for k, v in sorted(sb.get("branches", {}).items())) or "-",
                "yes" if agree else "NO",
            )
        )

    widths = [max(len(r[i]) for r in (*rows, head)) for i in range(len(head))]
    line = "  ".join(h.ljust(w) for h, w in zip(head, widths))
    print(line)
    print("-" * len(line))
    for r in rows:
        print("  ".join(c.ljust(w) for c, w in zip(r, widths)))
    print(
        f"\n{len(rows)} instances; "
        + ("FK-A, FK-B and the oracle agree on all of them" if not failures
           else f"{failures} DISAGREEMENTS")
    )
    return 1 if failures else 0


#: The scaling families, smallest first. ``self_dual_fano`` is the shape that
#: separates the two algorithms most sharply; the others are the standard
#: dualisation benchmarks from ``StandardProblems.m``.
FAMILIES = {
    "sdfp": ("SDFP(k), self-dualised Fano planes", self_dual_fano, (1, 2, 3)),
    "matching": ("M(v), disjoint pairs", matching, (4, 6, 8, 10, 12)),
    "threshold": ("TH(v), all pairs with an even endpoint", threshold, (4, 6, 8, 10)),
}


def cmd_benchmark(args: argparse.Namespace) -> int:
    """Run both algorithms on the self-transversal scaling families.

    Each family member is a genuine transversal pair by construction, so the
    answer is always TRUE and no branch can be cut short by an early
    counterexample. That is what makes the node counts comparable: both
    algorithms are being asked to do the whole job.
    """
    names = args.family or list(FAMILIES)
    head = ("family", "k/v", "n", "|G|", "|H|", "fk-a nodes", "sec",
            "fk-b nodes", "sec", "ratio")
    rows = []

    def run(algorithm, G, H):
        t0 = time.perf_counter()
        try:
            tree = algorithm(G, H)
        except RuntimeError:  # the max_nodes guard
            return None, time.perf_counter() - t0
        return tree, time.perf_counter() - t0

    for name in names:
        label, build, sizes = FAMILIES[name]
        for size in sizes:
            E = build(size)
            # SDFP is self-transversal, so it pairs with itself; the other
            # families are DNFs and pair with their computed dual.
            G, H = (E, E) if name == "sdfp" else (transversal(E), E)

            a, a_s = run(
                lambda g, h: fk_a(g, h, variant=args.fka_variant,
                                  max_nodes=args.max_nodes), G, H)
            b, b_s = run(
                lambda g, h: fk_b(g, h, variant=args.variant,
                                  split_rule=args.split_rule,
                                  max_nodes=args.max_nodes), G, H)
            if (a is not None and not a.dual) or (b is not None and not b.dual):
                print(f"!! {name}({size}) is not a transversal pair; aborting")
                return 1

            cap = f">{args.max_nodes}"
            rows.append((
                name, str(size), str(G.n), str(len(G)), str(len(H)),
                cap if a is None else str(len(a)), f"{a_s:.2f}",
                cap if b is None else str(len(b)), f"{b_s:.2f}",
                "-" if a is None or b is None else f"{len(a) / len(b):.1f}x",
            ))

    widths = [max(len(r[i]) for r in (*rows, head)) for i in range(len(head))]
    line = "  ".join(h.ljust(w) for h, w in zip(head, widths))
    print(line)
    print("-" * len(line))
    for r in rows:
        print("  ".join(c.ljust(w) for c, w in zip(r, widths)))
    print("\nratio is FK-A nodes / FK-B nodes; every pair here answers TRUE, so")
    print("neither algorithm can stop early on a counterexample.")
    for name in names:
        print(f"  {name:10} {FAMILIES[name][0]}")
    return 0


def cmd_dualize(args: argparse.Namespace) -> int:
    """Generate ``Tr(H)`` by repeated FK-B equivalence tests, and check it."""
    targets = _targets(args)
    if targets is None:
        return 2
    failures = 0
    for inst in targets:
        H = inst.H if args.side == "H" else inst.G
        t0 = time.perf_counter()
        run = dualise(H, variant=args.variant, split_rule=args.split_rule,
                      max_passes=args.max_passes)
        elapsed = time.perf_counter() - t0
        truth = transversal(H)
        ok = run.complete and run.cnf.edges == truth.edges
        failures += not ok
        print(
            f"{'ok  ' if ok else 'FAIL'} {inst.id:16} side={args.side} "
            f"|Tr|={len(truth):4} produced={len(run.cnf):4} passes={run.passes:4} "
            f"fk-b nodes={run.total_nodes:6} {elapsed:6.2f}s"
        )
    print(f"\n{'all match the Berge oracle' if not failures else f'{failures} FAILURES'}")
    return 1 if failures else 0


# ----------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fkb",
        description="Fredman-Khachiyan algorithm B: recursion trees and conflicting assignments.",
    )
    p.add_argument("--version", action="version", version=f"fkb {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    def add_algorithm_options(s: argparse.ArgumentParser) -> None:
        s.add_argument("instance", nargs="*", help="instance ids (or use --all)")
        s.add_argument("--all", action="store_true", help="use every instance")
        s.add_argument("--variant", choices=VARIANTS, default="faithful")
        s.add_argument("--split-rule", choices=SPLIT_RULES, default="most_frequent")

    s = sub.add_parser("run", help="run FK-B and write a report")
    add_algorithm_options(s)
    s.add_argument(
        "--all-variants",
        action="store_true",
        help=f"write one set of artifacts per variant ({', '.join(VARIANTS)})",
    )
    s.add_argument("--backend", choices=("python", "sage"), default=None)
    s.add_argument("--out", default=None, help="output directory (default: results/)")
    s.add_argument("--no-annotate", action="store_true", help="skip automorphism analysis")
    s.add_argument(
        "--no-graph-classes", action="store_true", help="skip primal graph classification"
    )
    s.add_argument("--dot", action="store_true", help="also write Graphviz DOT")
    s.add_argument(
        "--dot-view",
        choices=("structure", "automorphism", "properties"),
        default="structure",
    )
    s.add_argument(
        "--validate", action="store_true", help="verify assignments as they are built"
    )
    s.add_argument("--max-nodes", type=int, default=None, help="abort above this many nodes")
    s.set_defaults(func=cmd_run)

    sub.add_parser(
        "verify", help="cross-check FK-B against the brute-force oracle"
    ).set_defaults(func=cmd_verify)

    s = sub.add_parser("compare", help="run FK-A, FK-B and the oracle side by side")
    add_algorithm_options(s)
    s.add_argument(
        "--fka-variant",
        choices=("faithful", "modified"),
        default="modified",
        help="which FK-A variant to put in the comparison",
    )
    s.set_defaults(func=cmd_compare)

    s = sub.add_parser("benchmark", help="both algorithms on the scaling families")
    s.add_argument("family", nargs="*", choices=list(FAMILIES) or None, default=None)
    s.add_argument("--variant", choices=VARIANTS, default="faithful")
    s.add_argument("--split-rule", choices=SPLIT_RULES, default="most_frequent")
    s.add_argument("--fka-variant", choices=("faithful", "modified"), default="modified")
    s.add_argument(
        "--max-nodes",
        type=int,
        default=200_000,
        help="report a size as capped rather than running past this many nodes",
    )
    s.set_defaults(func=cmd_benchmark)

    s = sub.add_parser("dualize", help="generate Tr(H) by repeated FK-B tests")
    add_algorithm_options(s)
    s.add_argument("--side", choices=("G", "H"), default="H")
    s.add_argument("--max-passes", type=int, default=None)
    s.set_defaults(func=cmd_dualize)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
