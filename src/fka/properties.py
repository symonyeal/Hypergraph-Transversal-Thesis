"""Structural properties tracked across the recursion.

The annotations on the thesis' automorphism trees: conformality,
alpha-acyclicity, read-once, and the primal graph's class (:mod:`fka.graphs`).

"Conformal" and "normal" name the same predicate -- conformal is every
**maximal** clique of the primal graph lying in a term, normal is every clique
doing so, and since every clique sits inside a maximal one they cannot disagree.
Only :func:`is_conformal` exists; the read-once literature (Gurvich;
Golumbic-Mintz-Rotics) states the criterion with the other word.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

from .graphs import Graph, GraphClasses, classify_graph, is_cograph, primal_graph
from .hypergraph import Hypergraph, bits, verts

__all__ = [
    "is_conformal",
    "nonconformal_clique",
    "is_alpha_acyclic",
    "is_read_once",
    "Properties",
    "analyse",
]


def nonconformal_clique(MBF: Hypergraph, g: Optional[Graph] = None) -> Optional[tuple[int, ...]]:
    """A maximal clique of the primal graph inside no term, or ``None``.

    Returning the offending clique rather than a bare boolean is what makes the
    result usable: the thesis reports *which* clique breaks conformality when
    discussing the Fano plane and the asymmetric instances (Section 5.2.2).
    """
    g = primal_graph(MBF) if g is None else g
    for clique in g.maximal_cliques():
        c = bits(clique)
        if not any((c & e) == c for e in MBF.edges):
            return tuple(v + 1 for v in clique)
    return None


def is_conformal(MBF: Hypergraph, g: Optional[Graph] = None) -> bool:
    """Every maximal clique of the primal graph is contained in a term."""
    return nonconformal_clique(MBF, g) is None


def is_alpha_acyclic(MBF: Hypergraph) -> bool:
    """Decide alpha-acyclicity by GYO reduction.

    Repeatedly apply, in any order:

    1. delete a variable that lies in exactly one term ("ear removal");
    2. delete a term contained in another (and any empty term).

    ``MBF`` is alpha-acyclic exactly when this strips the hypergraph to nothing.
    The reduction is confluent, so the order of the steps does not matter.
    """
    edges = list(dict.fromkeys(MBF.edges))
    changed = True
    while changed and edges:
        changed = False

        # (2) drop empty edges and edges contained in another.
        keep: list[int] = []
        for i, e in enumerate(edges):
            if e == 0:
                changed = True
                continue
            if any((e & f) == e for j, f in enumerate(edges) if i != j and e != f):
                changed = True
                continue
            keep.append(e)
        # Duplicates: keep one copy of each.
        deduped = list(dict.fromkeys(keep))
        if len(deduped) != len(keep):
            changed = True
        edges = deduped
        if not edges:
            break

        # (1) delete vertices lying in exactly one edge.
        counts: dict[int, int] = {}
        for e in edges:
            for v in verts(e):
                counts[v] = counts.get(v, 0) + 1
        lonely = 0
        for v, c in counts.items():
            if c == 1:
                lonely |= 1 << v
        if lonely:
            edges = [e & ~lonely for e in edges]
            changed = True

    return not edges


def is_read_once(MBF: Hypergraph, g: Optional[Graph] = None) -> bool:
    """Decide whether ``MBF`` is a read-once hypergraph.

    Gurvich's criterion, as used by Golumbic-Mintz-Rotics [GMR06] and cited in
    thesis Section 4.4: a positive DNF is read-once iff it is normal (here:
    conformal) and its primal graph is P4-free (a cograph).

    An edgeless hypergraph is read-once vacuously.
    """
    if not MBF.edges:
        return True
    g = primal_graph(MBF) if g is None else g
    return is_conformal(MBF, g) and is_cograph(g)


@dataclass(frozen=True, slots=True)
class Properties:
    """Everything the annotation records about one side of a node."""

    n_vertices: int
    n_edges: int
    edge_sizes: tuple[int, ...]
    degrees: tuple[int, ...]
    irredundant: bool
    conformal: bool
    nonconformal_clique: Optional[tuple[int, ...]]
    alpha_acyclic: bool
    read_once: bool
    graph_classes: Optional[GraphClasses]
    isolated_vertices: tuple[int, ...]

    def labels(self) -> list[str]:
        """Short tags, in the style of the thesis' node annotations."""
        out = []
        if self.read_once:
            out.append("ROH")
        if self.alpha_acyclic:
            out.append("alpha-acyclic")
        if self.conformal:
            out.append("Conformal")
        else:
            out.append("Not conformal")
        return out

    def to_json(self) -> dict[str, Any]:
        return {
            "n_vertices": self.n_vertices,
            "n_edges": self.n_edges,
            "edge_sizes": list(self.edge_sizes),
            "degrees": list(self.degrees),
            "irredundant": self.irredundant,
            "conformal": self.conformal,
            "nonconformal_clique": (
                None if self.nonconformal_clique is None else list(self.nonconformal_clique)
            ),
            "alpha_acyclic": self.alpha_acyclic,
            "read_once": self.read_once,
            "graph_classes": (
                None if self.graph_classes is None else self.graph_classes.to_json()
            ),
            "isolated_vertices": [v + 1 for v in self.isolated_vertices],
            "labels": self.labels(),
        }


def analyse(MBF: Hypergraph, *, graph_classes: bool = True) -> Properties:
    """Run every property test on ``MBF``.

    ``graph_classes=False`` skips the primal-graph classification, which is the
    expensive part (it enumerates induced subgraphs). The hypergraph-level
    properties are cheap and always computed.
    """
    g = primal_graph(MBF)
    bad_clique = nonconformal_clique(MBF, g)
    conformal = bad_clique is None
    return Properties(
        n_vertices=MBF.n,
        n_edges=len(MBF),
        edge_sizes=tuple(MBF.edge_sizes()),
        degrees=tuple(MBF.degrees()),
        irredundant=MBF.is_sperner(),
        conformal=conformal,
        nonconformal_clique=bad_clique,
        alpha_acyclic=is_alpha_acyclic(MBF),
        read_once=(conformal and is_cograph(g)) if MBF.edges else True,
        graph_classes=classify_graph(g) if graph_classes else None,
        isolated_vertices=MBF.isolated_vertices(),
    )
