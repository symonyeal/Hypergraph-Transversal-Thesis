"""Hypergraph properties tracked across the recursion.

The annotations on the thesis' automorphism trees: conformality,
alpha-acyclicity, read-once, and the primal graph's class (:mod:`fka.graphs`).

"Conformal" and "normal" are the same predicate, though the legacy notebook
computed both -- *conformal* as every **maximal** clique of the primal graph
lying in a hyperedge, *normal* as every clique doing so. Every clique sits
inside a maximal one and containment is transitive, so the two can never
disagree, and the "normal" pass enumerated an exponentially larger set for
nothing. Only :func:`is_conformal` is computed; :func:`is_normal` is a
documented alias, because the read-once literature (Gurvich;
Golumbic-Mintz-Rotics) states the criterion with that word.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

from .graphs import Graph, GraphClasses, classify_graph, is_cograph, primal_graph
from .hypergraph import Hypergraph, bits, verts

__all__ = [
    "is_conformal",
    "is_normal",
    "nonconformal_clique",
    "is_alpha_acyclic",
    "is_read_once",
    "HypergraphProperties",
    "analyse",
]


def nonconformal_clique(H: Hypergraph, g: Optional[Graph] = None) -> Optional[tuple[int, ...]]:
    """A maximal clique of the primal graph inside no hyperedge, or ``None``.

    Returning the offending clique rather than a bare boolean is what makes the
    result usable: the thesis reports *which* clique breaks conformality when
    discussing the Fano plane and the asymmetric instances (Section 5.2.2).
    """
    g = primal_graph(H) if g is None else g
    for clique in g.maximal_cliques():
        c = bits(clique)
        if not any((c & e) == c for e in H.edges):
            return tuple(v + 1 for v in clique)
    return None


def is_conformal(H: Hypergraph, g: Optional[Graph] = None) -> bool:
    """Every maximal clique of the primal graph is contained in a hyperedge."""
    return nonconformal_clique(H, g) is None


#: Alias -- see the module docstring for why these coincide.
is_normal = is_conformal


def is_alpha_acyclic(H: Hypergraph) -> bool:
    """Decide alpha-acyclicity by GYO reduction.

    Repeatedly apply, in any order:

    1. delete a vertex that lies in exactly one hyperedge ("ear removal");
    2. delete a hyperedge contained in another (and any empty hyperedge).

    ``H`` is alpha-acyclic exactly when this strips the hypergraph to nothing.
    The reduction is confluent, so the order of the steps does not matter.
    """
    edges = list(dict.fromkeys(H.edges))
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


def is_read_once(H: Hypergraph, g: Optional[Graph] = None) -> bool:
    """Decide whether ``H`` is a read-once hypergraph.

    Gurvich's criterion, as used by Golumbic-Mintz-Rotics [GMR06] and cited in
    thesis Section 4.4: a positive DNF is read-once iff it is normal (here:
    conformal) and its primal graph is P4-free (a cograph).

    An edgeless hypergraph is read-once vacuously.
    """
    if not H.edges:
        return True
    g = primal_graph(H) if g is None else g
    return is_conformal(H, g) and is_cograph(g)


@dataclass(frozen=True, slots=True)
class HypergraphProperties:
    """Everything the recursion-tree annotation records about one hypergraph."""

    n_vertices: int
    n_edges: int
    edge_sizes: tuple[int, ...]
    degrees: tuple[int, ...]
    sperner: bool
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
            "sperner": self.sperner,
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


def analyse(H: Hypergraph, *, graph_classes: bool = True) -> HypergraphProperties:
    """Run every property test on ``H``.

    ``graph_classes=False`` skips the primal-graph classification, which is the
    expensive part (it enumerates induced subgraphs). The hypergraph-level
    properties are cheap and always computed.
    """
    g = primal_graph(H)
    bad_clique = nonconformal_clique(H, g)
    conformal = bad_clique is None
    return HypergraphProperties(
        n_vertices=H.n,
        n_edges=len(H),
        edge_sizes=tuple(H.edge_sizes()),
        degrees=tuple(H.degrees()),
        sperner=H.is_sperner(),
        conformal=conformal,
        nonconformal_clique=bad_clique,
        alpha_acyclic=is_alpha_acyclic(H),
        read_once=(conformal and is_cograph(g)) if H.edges else True,
        graph_classes=classify_graph(g) if graph_classes else None,
        isolated_vertices=H.isolated_vertices(),
    )
