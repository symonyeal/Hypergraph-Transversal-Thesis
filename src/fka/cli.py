"""Command-line interface: ``python -m fka <command>``.

The whole experiment workflow lives here so that a result can be reproduced
from a shell line recorded in the research log, rather than from remembering
which notebook cells to run in which order.

    python -m fka env                       what this interpreter can do
    python -m fka list                      the instance library
    python -m fka show fano                 one instance, with its baseline
    python -m fka run fano                  FK-A + annotation + HTML report
    python -m fka run --all --all-variants  every instance, every variant
    python -m fka verify                    check FK-A against the oracle
    python -m fka refresh                   recompute instance baselines

Everything works identically under ``sage -python -m fka ...``; see
:mod:`fka.backends`.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional, Sequence

from . import __version__
from .analysis import AnalysisOptions, annotate
from .algorithm import SPLIT_RULES, VARIANTS, fk_a
from .backends import describe
from .dot import write_dot
from .instances import Instance, list_ids, load, load_all, refresh_expected
from .report import write_report
from .transversal import is_dual_oracle


def _results_dir(arg: Optional[str]) -> Path:
    if arg:
        return Path(arg)
    return Path(__file__).resolve().parents[2] / "results"


# ----------------------------------------------------------------------
def cmd_env(args: argparse.Namespace) -> int:
    info = describe()
    width = max(len(k) for k in info)
    for k, v in info.items():
        print(f"{k:<{width}}  {v}")
    if not info["sage_available"]:
        print(
            "\nSageMath is not importable here, so group naming uses the built-in "
            "catalogue and graph-class tests use the capped pure-Python versions. "
            "See docs/sagemath.md."
        )
    return 0


def cmd_list(args: argparse.Namespace) -> int:
    rows = []
    for inst in load_all():
        e = inst.expected
        rows.append(
            (
                inst.id,
                inst.name,
                inst.provenance,
                str(e.get("dual", "?")),
                str(e.get("epsilon", "?")),
                f"{e.get('nC','?')}/{e.get('nD','?')}",
                str((e.get("fka", {}).get("modified", {}) or {}).get("nodes", "?")),
                str((e.get("fkb", {}).get("faithful", {}) or {}).get("nodes", "?")),
            )
        )
    head = (
        "id", "name", "provenance", "dual", "eps", "nC/nD", "fk-a nodes", "fk-b nodes"
    )
    widths = [max(len(str(r[i])) for r in (*rows, head)) for i in range(len(head))]
    line = "  ".join(h.ljust(w) for h, w in zip(head, widths))
    print(line)
    print("-" * len(line))
    for r in rows:
        print("  ".join(str(c).ljust(w) for c, w in zip(r, widths)))
    print(f"\n{len(rows)} instances in {load_all()[0].path.parent if rows else '(none)'}")
    return 0


def cmd_show(args: argparse.Namespace) -> int:
    inst = load(args.instance)
    print(f"{inst.name}  [{inst.id}]")
    print(f"  family      {inst.family}")
    print(f"  source      {inst.source}")
    print(f"  provenance  {inst.provenance}")
    print(f"  variables   {inst.n_vertices}")
    print(f"  C ({len(inst.cnf_edges)} clauses)   {inst.cnf.label()}")
    print(f"  D ({len(inst.dnf_edges)} monomials) {inst.dnf.label()}")
    if inst.expected:
        print("  expected    " + json.dumps(inst.expected))
    if args.notes and inst.notes:
        print("\n  " + inst.notes.replace(". ", ".\n  "))
    return 0


def _run_one(inst: Instance, args: argparse.Namespace, outdir: Path) -> dict:
    tree = fk_a(
        inst.cnf,
        inst.dnf,
        variant=args.variant,
        split_rule=args.split_rule,
        instance=inst.id,
        validate=args.validate,
        max_nodes=args.max_nodes,
    )
    if not args.no_annotate:
        annotate(
            tree,
            AnalysisOptions(
                graph_classes=not args.no_graph_classes,
                backend=args.backend,
            ),
        )

    stem = f"{inst.id}.fk-a.{args.variant}"
    json_path = tree.save(outdir / f"{stem}.json")
    html_path = write_report(tree, outdir / f"{stem}.html", title=inst.name, instance=inst)
    outputs = [json_path, html_path]
    if args.dot:
        outputs.append(write_dot(tree, outdir / f"{stem}.dot", view=args.dot_view))

    summary = tree.summary()
    baseline = (inst.expected.get("fka", {}) or {}).get(args.variant, {})
    flag = ""
    if baseline and args.split_rule == "degree_sum":
        if baseline.get("nodes") != summary["nodes"]:
            flag = f"  !! baseline says {baseline['nodes']} nodes"
    expected_dual = inst.expected.get("dual")
    if expected_dual is not None and expected_dual != tree.dual:
        flag += f"  !! baseline says dual={expected_dual}"

    print(
        f"{inst.id:16} dual={str(tree.dual):5} nodes={summary['nodes']:4} "
        f"depth={summary['depth']:2}{flag}"
    )
    for p in outputs:
        print(f"                 -> {p}")
    return summary


def cmd_run(args: argparse.Namespace) -> int:
    outdir = _results_dir(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    if args.all:
        targets = load_all()
    elif args.instance:
        targets = [load(i) for i in args.instance]
    else:
        print("give an instance id, or --all. Available: " + ", ".join(list_ids()))
        return 2
    variants = VARIANTS if args.all_variants else (args.variant,)
    for inst in targets:
        for variant in variants:
            args.variant = variant
            _run_one(inst, args, outdir)
    return 0


def cmd_verify(args: argparse.Namespace) -> int:
    """Check FK-A against the brute-force oracle on every instance and variant."""
    failures = 0
    for inst in load_all():
        truth = is_dual_oracle(inst.cnf, inst.dnf)
        for variant in VARIANTS:
            for rule in SPLIT_RULES:
                tree = fk_a(
                    inst.cnf,
                    inst.dnf,
                    variant=variant,
                    split_rule=rule,
                    instance=inst.id,
                    validate=True,
                )
                ok = tree.dual == truth
                if not ok:
                    failures += 1
                print(
                    f"{'ok  ' if ok else 'FAIL'} {inst.id:16} {variant:9} {rule:10} "
                    f"fk-a={str(tree.dual):5} oracle={str(truth):5} nodes={len(tree)}"
                )
    print(f"\n{'all agree with the oracle' if not failures else f'{failures} FAILURES'}")
    return 1 if failures else 0


def cmd_refresh(args: argparse.Namespace) -> int:
    for inst in load_all():
        before = dict(inst.expected)
        refresh_expected(inst)
        inst.save()
        changed = "changed" if before != inst.expected else "unchanged"
        print(f"{inst.id:16} {changed}")
        if changed == "changed" and before:
            print(f"                 was {json.dumps(before)}")
            print(f"                 now {json.dumps(inst.expected)}")
    return 0


# ----------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fka",
            description=(
            "Fredman-Khachiyan algorithm A: recursion trees and symmetry."
        ),
    )
    p.add_argument("--version", action="version", version=f"fka {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    sub.add_parser("env", help="report interpreter and backend capabilities").set_defaults(
        func=cmd_env
    )
    sub.add_parser("list", help="list the instance library").set_defaults(func=cmd_list)

    s = sub.add_parser("show", help="print one instance and its baseline")
    s.add_argument("instance")
    s.add_argument(
        "--notes",
        action="store_true",
        help="also print research notes and interpretation attached to the instance",
    )
    s.set_defaults(func=cmd_show)

    s = sub.add_parser("run", help="run FK-A and write a report")
    s.add_argument("instance", nargs="*", help="instance ids (or use --all)")
    s.add_argument("--all", action="store_true", help="run every instance")
    s.add_argument("--variant", choices=VARIANTS, default="modified")
    s.add_argument(
        "--all-variants",
        action="store_true",
        help=f"write one set of artifacts per variant ({', '.join(VARIANTS)})",
    )
    s.add_argument("--split-rule", choices=SPLIT_RULES, default="degree_sum")
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
    s.add_argument("--validate", action="store_true", help="verify certificates as they are built")
    s.add_argument("--max-nodes", type=int, default=None, help="abort above this many nodes")
    s.set_defaults(func=cmd_run)

    sub.add_parser(
        "verify", help="cross-check FK-A against the brute-force oracle"
    ).set_defaults(func=cmd_verify)
    sub.add_parser(
        "refresh", help="recompute the expected block of every instance file"
    ).set_defaults(func=cmd_refresh)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
