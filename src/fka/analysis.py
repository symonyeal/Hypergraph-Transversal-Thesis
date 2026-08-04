"""Annotating a recursion tree -- the thesis' "automorphism tree".

Chapters 3 and 5 study what happens to symmetry as the instance is decomposed.
The artefact is a second tree of the same shape, each node carrying the
automorphism group and structural properties of its two hypergraphs rather than
the hypergraphs themselves (Figures 5.1-5.4).

A separate pass over a finished :class:`~fka.tree.RecursionTree`, not
interleaved with the recursion: the algorithms stay pure decision procedures, a
tree can be re-annotated with a different backend without re-running one, and
the expensive part can be skipped or capped on its own. It applies unchanged to
FK-B trees, which is what makes the two comparable.

The thesis records one group per node, observing that "the input hypergraphs G
and H (if they are transversals of each other) have the same automorphism
group" (p.51). Both sides are computed anyway and
:attr:`NodeAnalysis.aut_agree` records whether they matched -- measured over
3455 nodes of both algorithms, they always do; see ``docs/fk-a-vs-fk-b.md`` §7.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

from .automorphism import AutomorphismResult, automorphism_group
from .hypergraph import Hypergraph
from .properties import HypergraphProperties, analyse
from .tree import RecursionNode, RecursionTree

__all__ = ["NodeAnalysis", "annotate_node", "annotate", "AnalysisOptions"]


@dataclass(slots=True)
class AnalysisOptions:
    """Knobs for the annotation pass.

    ``max_vertices_for_graph_classes`` guards the primal-graph classification,
    which enumerates induced subgraphs (see :mod:`fka.graphs`). Nodes above the
    limit still get hypergraph-level properties and automorphism groups.

    ``skip_group_above`` skips the automorphism group when a node has more
    hyperedges than this. The search is fast on the thesis instances but is
    exponential in the worst case, and one runaway node should not stall a
    whole experiment.

    Both defaults are set from measurement, not from taste. On the benchmark
    instance ``sdfp-sd-2`` -- 16 vertices, 64 edges, ``|Aut| = 56448`` -- the
    group search takes about 19 seconds and the graph classification about 4.6,
    *per distinct node*, against roughly a millisecond on every thesis instance.
    A tree with a thousand nodes therefore stops being annotatable somewhere
    between those two regimes, so the caps sit below the benchmark and above
    everything in the published work: the largest thesis instance is
    ``trivial-aut-1`` at 35 edges and 8 vertices, and it is comfortably inside.
    Raise them deliberately, per run, when a specific group is the question.
    """

    graph_classes: bool = True
    max_vertices_for_graph_classes: int = 14
    skip_group_above: int = 40
    include_isolated: bool = False
    backend: Optional[str] = None


@dataclass(slots=True)
class NodeAnalysis:
    """Automorphism and property data for one recursion-tree node."""

    aut_G: Optional[AutomorphismResult]
    aut_H: Optional[AutomorphismResult]
    props_G: HypergraphProperties
    props_H: HypergraphProperties

    @property
    def aut_agree(self) -> Optional[bool]:
        """Whether ``Aut(G)`` and ``Aut(H)`` have the same order and name.

        ``None`` when either group was skipped.
        """
        if self.aut_G is None or self.aut_H is None:
            return None
        return (
            self.aut_G.order == self.aut_H.order
            and self.aut_G.name == self.aut_H.name
        )

    def label(self) -> str:
        """The one-line group label used on tree nodes.

        Matches the thesis' convention of showing a single group per node when
        the two sides agree, and flagging the disagreement when they do not.
        """
        if self.aut_G is None or self.aut_H is None:
            return ""
        if self.aut_agree:
            return self.aut_G.name
        return f"G: {self.aut_G.name} / H: {self.aut_H.name}"

    def to_json(self) -> dict[str, Any]:
        return {
            "aut_G": None if self.aut_G is None else self.aut_G.to_json(),
            "aut_H": None if self.aut_H is None else self.aut_H.to_json(),
            "aut_agree": self.aut_agree,
            "aut_label": self.label(),
            "props_G": self.props_G.to_json(),
            "props_H": self.props_H.to_json(),
        }


def _group_for(
    H: Hypergraph, opts: AnalysisOptions
) -> Optional[AutomorphismResult]:
    if len(H) > opts.skip_group_above:
        return None
    if not H.edges:
        # Every permutation of the (empty) support is an automorphism; on the
        # compacted hypergraph that is the trivial group, which is the honest
        # answer and avoids the symmetric-group blow-up.
        return None
    try:
        return automorphism_group(
            H, include_isolated=opts.include_isolated, backend=opts.backend
        )
    except ValueError:
        # Group above the enumeration cap -- recorded as unknown rather than
        # aborting the whole annotation pass.
        return None


def annotate_node(node: RecursionNode, opts: Optional[AnalysisOptions] = None) -> NodeAnalysis:
    """Compute and attach the analysis for a single node."""
    opts = opts or AnalysisOptions()
    do_classes = opts.graph_classes and node.G.n <= opts.max_vertices_for_graph_classes
    result = NodeAnalysis(
        aut_G=_group_for(node.G, opts),
        aut_H=_group_for(node.H, opts),
        props_G=analyse(node.G, graph_classes=do_classes),
        props_H=analyse(node.H, graph_classes=do_classes),
    )
    node.analysis = result.to_json()
    return result


def annotate(
    tree: RecursionTree, opts: Optional[AnalysisOptions] = None
) -> RecursionTree:
    """Annotate every node of ``tree`` in place, and return it.

    Nodes are memoised on their ``(G, H)`` pair: the thesis observes that FK-A
    regenerates identical subproblems -- "The resulting contraction
    subhypergraph at Node 8 is identical to that at Node 17. This repetition of
    subproblems contributes to the enlargement of the recursion tree" (p.53) --
    so on symmetric instances this avoids recomputing the same automorphism
    group many times over.
    """
    opts = opts or AnalysisOptions()
    cache: dict[tuple, dict[str, Any]] = {}
    for node in tree:
        key = (node.G.n, node.G.edges, node.H.edges)
        cached = cache.get(key)
        if cached is not None:
            node.analysis = cached
            continue
        annotate_node(node, opts)
        cache[key] = node.analysis
    return tree
