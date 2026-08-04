"""Graphviz DOT export.

:mod:`fka.report` is the primary output and needs nothing installed. DOT is kept
for what it is genuinely better at: a vector figure for the LaTeX source of a
paper, and feeding the tree to other Graphviz-aware tooling.

Writing the text is pure Python; rendering it needs the ``dot`` binary, which
this module neither requires nor shells out to. :func:`to_dot` returns the
source and what happens next is the caller's business::

    dot -Tpdf results/fano.fk-a.modified.dot -o fano.pdf
"""

from __future__ import annotations

from pathlib import Path

from .tree import RecursionTree

__all__ = ["to_dot", "write_dot"]

_VIEWS = ("structure", "automorphism", "properties")


def _escape(s: str) -> str:
    return s.replace("\\", "\\\\").replace('"', '\\"')


def _sides(tree: RecursionTree) -> tuple[str, str]:
    """The two side names this algorithm uses, e.g. ``("G", "H")``."""
    from .report import ALGORITHM_STYLE

    style = ALGORITHM_STYLE.get(tree.algorithm or "FK-A", ALGORITHM_STYLE["FK-A"])
    a, b = style["sides"].split(",")
    return a, b


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
        left, right = _sides(tree)
        lines.append(f"|{left}|={len(node.G)}  |{right}|={len(node.H)}")
        if node.pivot is not None:
            lines.append(f"split {node.pivot_label()}")
        elif node.verdict is not None:
            lines.append(node.verdict.reason)
    return "\\n".join(_escape(x) for x in lines)


def to_dot(tree: RecursionTree, *, view: str = "structure", rankdir: str = "TB") -> str:
    """Render ``tree`` as Graphviz DOT source.

    ``view`` selects the node labels, matching the HTML report's three views.
    Edge styles come from :meth:`RecursionNode.branch_style`, so FK-A's ``L``
    is solid and ``R`` dashed as in the thesis figures, and FK-B's extra
    per-clause and per-monomial subproblems are dotted.
    """
    if view not in _VIEWS:
        raise ValueError(f"unknown view {view!r}; expected one of {_VIEWS}")

    algorithm = tree.algorithm or "FK-A"
    out = [
        f"digraph {algorithm.replace('-', '')} {{",
        f"  rankdir={rankdir};",
        '  graph [fontname="Helvetica", labelloc="t", '
        f'label="{_escape(tree.instance or algorithm)}  ({algorithm}, {view} view)"];',
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
            out.append(
                f'  n{node.node_id} -> n{child_id} '
                f'[label="{_escape(child.branch)}", style={child.branch_style()}];'
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
