"""Annotating a recursion tree -- the thesis' "automorphism tree".

Chapters 3 and 5 study what happens to symmetry as the instance is decomposed.
The artefact is a second tree of the same shape, each node carrying the
automorphism group and structural properties of its two sides rather than the
sides themselves (Figures 5.1-5.4).

A separate pass over a finished :class:`~fka.tree.Tree`, not interleaved with
the recursion: the algorithms stay pure decision procedures, a tree can be
re-annotated with a different backend without re-running one, and the expensive
part can be capped on its own. It applies unchanged to FK-B trees, which is what
makes the two comparable.

The thesis records one group per node, observing that the two sides of a dual
pair "have the same automorphism group" (p.51). Both are computed anyway and
:attr:`NodeAnalysis.aut_agree` records whether they matched.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

from .automorphism import AutomorphismResult, automorphism_group
from .hypergraph import Hypergraph
from .properties import Properties, analyse
from .tree import Node, Tree

__all__ = ["NodeAnalysis", "annotate_node", "annotate", "AnalysisOptions"]


@dataclass(slots=True)
class AnalysisOptions:
    """Knobs for the annotation pass.

    Both caps are set from measurement. On ``sdfp-sd-2`` -- 16 variables, 64
    terms, ``|Aut| = 56448`` -- the group search takes about 19 seconds and the
    graph classification about 4.6, *per distinct node*, against roughly a
    millisecond on every thesis instance. They sit below that benchmark and above
    everything in the published work, whose largest instance is ``trivial-aut-1``
    at 35 terms. Raise them deliberately, per run, when a specific group is the
    question.
    """

    graph_classes: bool = True
    max_vertices_for_graph_classes: int = 14
    skip_group_above: int = 40
    include_isolated: bool = False
    backend: Optional[str] = None


@dataclass(slots=True)
class NodeAnalysis:
    """Automorphism and property data for one node."""

    aut_C: Optional[AutomorphismResult]
    aut_D: Optional[AutomorphismResult]
    props_C: Properties
    props_D: Properties

    @property
    def aut_agree(self) -> Optional[bool]:
        """Whether the two sides' groups have the same order and name.

        ``None`` when either was skipped.
        """
        if self.aut_C is None or self.aut_D is None:
            return None
        return (
            self.aut_C.order == self.aut_D.order
            and self.aut_C.name == self.aut_D.name
        )

    def label(self) -> str:
        """The one-line group label used on tree nodes.

        The thesis shows a single group per node when the two sides agree, and
        flags the disagreement when they do not.
        """
        if self.aut_C is None or self.aut_D is None:
            return ""
        if self.aut_agree:
            return self.aut_C.name
        return f"C: {self.aut_C.name} / D: {self.aut_D.name}"

    def to_json(self) -> dict[str, Any]:
        return {
            "aut_C": None if self.aut_C is None else self.aut_C.to_json(),
            "aut_D": None if self.aut_D is None else self.aut_D.to_json(),
            "aut_agree": self.aut_agree,
            "aut_label": self.label(),
            "props_C": self.props_C.to_json(),
            "props_D": self.props_D.to_json(),
        }


def _group_for(
    MBF: Hypergraph, opts: AnalysisOptions
) -> Optional[AutomorphismResult]:
    if len(MBF) > opts.skip_group_above:
        return None
    if not MBF.edges:
        # Every permutation of the (empty) support is an automorphism; on the
        # compacted function that is the trivial group, which is the honest
        # answer and avoids the symmetric-group blow-up.
        return None
    try:
        return automorphism_group(
            MBF, include_isolated=opts.include_isolated, backend=opts.backend
        )
    except ValueError:
        # Group above the enumeration cap -- recorded as unknown rather than
        # aborting the whole pass.
        return None


def annotate_node(node: Node, opts: Optional[AnalysisOptions] = None) -> NodeAnalysis:
    """Compute and attach the analysis for a single node."""
    opts = opts or AnalysisOptions()
    do_classes = (
        opts.graph_classes and node.cnf.n <= opts.max_vertices_for_graph_classes
    )
    result = NodeAnalysis(
        aut_C=_group_for(node.cnf, opts),
        aut_D=_group_for(node.dnf, opts),
        props_C=analyse(node.cnf, graph_classes=do_classes),
        props_D=analyse(node.dnf, graph_classes=do_classes),
    )
    node.analysis = result.to_json()
    return result


def annotate(tree: Tree, opts: Optional[AnalysisOptions] = None) -> Tree:
    """Annotate every node of ``tree`` in place, and return it.

    Nodes are memoised on their ``(cnf, dnf)`` pair: FK-A regenerates identical
    subproblems -- "This repetition of subproblems contributes to the enlargement
    of the recursion tree" (p.53) -- so on symmetric instances this avoids
    recomputing the same group many times over.
    """
    opts = opts or AnalysisOptions()
    cache: dict[tuple, dict[str, Any]] = {}
    for node in tree:
        key = (node.cnf.n, node.cnf.edges, node.dnf.edges)
        cached = cache.get(key)
        if cached is not None:
            node.analysis = cached
            continue
        annotate_node(node, opts)
        cache[key] = node.analysis
    return tree
