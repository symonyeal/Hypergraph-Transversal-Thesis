"""Graphviz DOT export for recursion trees.

The HTML report in :mod:`fka.report` is the primary output and needs nothing
installed. DOT is kept for the cases it is genuinely better at: dropping a
vector figure into the LaTeX source of a paper, and feeding the tree to other
Graphviz-aware tooling.

Writing the DOT text is pure Python. Rendering it to PNG/PDF/SVG needs the
``dot`` binary, which this module does not require and does not shell out to --
:func:`to_dot` returns the source and it is the caller's business what to do
with it::

    dot -Tpdf results/fano.dot -o results/fano.pdf
"""

from __future__ import annotations

from pathlib import Path

from .tree import RecursionTree

__all__ = ["to_dot", "write_dot"]

_VIEWS = ("structure", "automorphism", "properties")


def _escape(s: str) -> str:
    return s.replace("\\", "\\\\").replace('"', '\\"')


def _node_label(tree: RecursionTree, node_id: int, view: str) -> str:
    node = tree[node_id]
    a = node.analysis or {}
    lines = [f"Node {node.node_id}"]
    if view == "automorphism":
        lines.append(a.get("aut_label") or "-")
        aut = a.get("aut_G")
        if aut:
            lines.append(f"order {aut['order']}")
    elif view == "properties":
        props = a.get("props_G")
        if props:
            lines.append(" / ".join(props.get("labels", [])) or "-")
            classes = (props.get("graph_classes") or {}).get("names", [])
            if classes:
                lines.append(", ".join(classes[:2]))
    else:
        lines.append(f"|G|={len(node.G)}  |H|={len(node.H)}")
        if node.pivot is not None:
            lines.append(f"split {node.pivot_label()}")
        elif node.verdict is not None:
            lines.append(node.verdict.reason)
    return "\\n".join(_escape(x) for x in lines)


def to_dot(tree: RecursionTree, *, view: str = "structure", rankdir: str = "TB") -> str:
    """Render ``tree`` as Graphviz DOT source.

    ``view`` selects the node labels, matching the HTML report's three views.
    Left (``L``) edges are solid and right (``R``) edges dashed, as in the
    report and the thesis figures.
    """
    if view not in _VIEWS:
        raise ValueError(f"unknown view {view!r}; expected one of {_VIEWS}")

    out = [
        "digraph FKA {",
        f"  rankdir={rankdir};",
        '  graph [fontname="Helvetica", labelloc="t", '
        f'label="{_escape(tree.instance or "FK-A")}  ({view} view)"];',
        '  node [shape=box, style="rounded,filled", fillcolor="#ffffff", '
        'fontname="Helvetica", fontsize=10];',
        '  edge [fontname="Helvetica", fontsize=9, color="#888888"];',
    ]
    for node in tree:
        colour = "#ffffff"
        if node.verdict is not None and node.is_leaf:
            colour = "#eaf3ef" if node.verdict.dual else "#f7e8e4"
        out.append(
            f'  n{node.node_id} [label="{_node_label(tree, node.node_id, view)}", '
            f'fillcolor="{colour}"];'
        )
    for node in tree:
        for child_id in node.children:
            child = tree[child_id]
            style = "solid" if child.branch == "L" else "dashed"
            out.append(
                f'  n{node.node_id} -> n{child_id} '
                f'[label="{child.branch}", style={style}];'
            )
    out.append("}")
    return "\n".join(out) + "\n"


def write_dot(
    tree: RecursionTree, path: str | Path, *, view: str = "structure"
) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(to_dot(tree, view=view), encoding="utf-8")
    return target
