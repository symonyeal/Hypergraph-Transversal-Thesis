"""The FK-A recursion tree.

Every recursive call of FK-A becomes one :class:`RecursionNode`. The thesis
treats this tree as the primary object of study -- "The performance of FK-A is
determined by the number of subproblems it generates, and this is reflected in
the number of nodes of the recursion tree" (thesis p.12) -- and Chapter 5
annotates each node with the automorphism group and hypergraph properties of
its inputs.

The legacy implementation tracked this in module-level lists and dicts keyed by
a global ``x`` that was never assigned, so every node overwrote the same entry.
Here the tree is an explicit, serialisable structure: nodes carry their own
state, ``to_json``/``from_json`` round-trip losslessly, and analysis passes
(automorphism groups, primal-graph classes) attach to nodes afterwards rather
than being interleaved with the recursion.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator, Optional

from .hypergraph import Hypergraph, bitset_to_vertices

__all__ = ["RecursionNode", "RecursionTree", "Verdict"]


@dataclass(frozen=True, slots=True)
class Verdict:
    """Why a node concluded what it concluded.

    ``reason`` is a stable machine-readable code so experiment output can be
    grouped without string-matching prose:

    ``"dual"``               subtree confirmed duality
    ``"trivial"``            resolved by the ``|G|*|H| <= 1`` check (Alg. 1 Step 3)
    ``"single_edge"``        resolved by the modified ``D_{1,k}`` check (thesis 5.1.1)
    ``"cond_i_support"``     precondition (i) failed: ground sets differ
    ``"cond_ii_size"``       precondition (ii) failed: an edge is larger than the other side
    ``"cond_iii_disjoint"``  precondition (iii) failed: a disjoint edge pair exists
    ``"cond_iv_frequency"``  precondition (iv) failed: the ``2^-|e|`` sum is below 1
    ``"child_failed"``       a subproblem returned not-dual

    ``certificate`` is a bitset ``S`` satisfying thesis equation 2.1 -- it hits
    every edge of ``G`` and contains no edge of ``H`` -- when one could be
    constructed soundly. Precondition failures retain the thesis' structural
    witness (``witness_edges``) and also carry a set certificate when the exact
    search finds one. Never treat ``None`` as "no certificate exists"; it means
    none was built in this orientation.
    """

    dual: bool
    reason: str
    detail: str = ""
    certificate: Optional[int] = None
    witness_edges: tuple[int, ...] = ()

    def to_json(self) -> dict[str, Any]:
        return {
            "dual": self.dual,
            "reason": self.reason,
            "detail": self.detail,
            "certificate": (
                None
                if self.certificate is None
                else [v + 1 for v in bitset_to_vertices(self.certificate)]
            ),
            "witness_edges": [
                [v + 1 for v in bitset_to_vertices(e)] for e in self.witness_edges
            ],
        }


@dataclass(slots=True)
class RecursionNode:
    """One call of FK-A.

    ``G_in``/``H_in`` are the inputs as handed down by the parent;
    ``G``/``H`` are the same pair after Sperner reduction. The difference is
    exactly what the thesis' node dumps label "Non Sperner Edges".
    """

    node_id: int
    parent_id: Optional[int]
    path: tuple[str, ...]
    G_in: Hypergraph
    H_in: Hypergraph
    G: Hypergraph
    H: Hypergraph
    removed_G: tuple[int, ...] = ()
    removed_H: tuple[int, ...] = ()
    pivot: Optional[int] = None
    pivot_frequency: Optional[float] = None
    epsilon: Optional[float] = None
    children: list[int] = field(default_factory=list)
    verdict: Optional[Verdict] = None
    # Filled in by analysis passes (automorphism groups, property checks).
    analysis: dict[str, Any] = field(default_factory=dict)

    @property
    def depth(self) -> int:
        return len(self.path)

    @property
    def is_leaf(self) -> bool:
        return not self.children

    @property
    def branch(self) -> str:
        """``"L"``, ``"R"``, or ``"root"``."""
        return self.path[-1] if self.path else "root"

    def pivot_label(self) -> str:
        """The split vertex in the thesis' 1-indexed ``v_i`` notation."""
        return "" if self.pivot is None else f"v{self.pivot + 1}"

    def to_json(self) -> dict[str, Any]:
        return {
            "node_id": self.node_id,
            "parent_id": self.parent_id,
            "path": list(self.path),
            "depth": self.depth,
            "G": self.G.to_sets(one_indexed=True),
            "H": self.H.to_sets(one_indexed=True),
            "G_in": self.G_in.to_sets(one_indexed=True),
            "H_in": self.H_in.to_sets(one_indexed=True),
            "n": self.G.n,
            "removed_G": [
                [v + 1 for v in bitset_to_vertices(e)] for e in self.removed_G
            ],
            "removed_H": [
                [v + 1 for v in bitset_to_vertices(e)] for e in self.removed_H
            ],
            "pivot": None if self.pivot is None else self.pivot + 1,
            "pivot_frequency": self.pivot_frequency,
            "epsilon": self.epsilon,
            "children": list(self.children),
            "verdict": None if self.verdict is None else self.verdict.to_json(),
            "analysis": self.analysis,
        }

    @classmethod
    def from_json(cls, d: dict[str, Any]) -> "RecursionNode":
        n = d["n"]

        def hg(key: str) -> Hypergraph:
            return Hypergraph.from_sets(n, d[key], one_indexed=True)

        def bits(key: str) -> tuple[int, ...]:
            return tuple(
                Hypergraph.from_sets(n, [e], one_indexed=True).edges[0] if e else 0
                for e in d[key]
            )

        verdict = None
        if d.get("verdict"):
            v = d["verdict"]
            cert = v.get("certificate")
            verdict = Verdict(
                dual=v["dual"],
                reason=v["reason"],
                detail=v.get("detail", ""),
                certificate=(
                    None
                    if cert is None
                    else sum(1 << (x - 1) for x in cert)
                ),
                witness_edges=tuple(
                    sum(1 << (x - 1) for x in e) for e in v.get("witness_edges", [])
                ),
            )
        return cls(
            node_id=d["node_id"],
            parent_id=d["parent_id"],
            path=tuple(d["path"]),
            G_in=hg("G_in"),
            H_in=hg("H_in"),
            G=hg("G"),
            H=hg("H"),
            removed_G=bits("removed_G"),
            removed_H=bits("removed_H"),
            pivot=None if d["pivot"] is None else d["pivot"] - 1,
            pivot_frequency=d.get("pivot_frequency"),
            epsilon=d.get("epsilon"),
            children=list(d.get("children", [])),
            verdict=verdict,
            analysis=d.get("analysis", {}),
        )


@dataclass(slots=True)
class RecursionTree:
    """The full trace of one FK-A run."""

    nodes: dict[int, RecursionNode] = field(default_factory=dict)
    root_id: Optional[int] = None
    instance: str = ""
    variant: str = ""
    pivot_rule: str = ""
    dual: Optional[bool] = None
    verdict: Optional[Verdict] = None

    # ------------------------------------------------------------------
    def add(self, node: RecursionNode) -> RecursionNode:
        self.nodes[node.node_id] = node
        if node.parent_id is None:
            self.root_id = node.node_id
        else:
            self.nodes[node.parent_id].children.append(node.node_id)
        return node

    def __len__(self) -> int:
        return len(self.nodes)

    def __getitem__(self, node_id: int) -> RecursionNode:
        return self.nodes[node_id]

    def __iter__(self) -> Iterator[RecursionNode]:
        """Iterate in node-id order, which is the order the thesis numbers them."""
        for k in sorted(self.nodes):
            yield self.nodes[k]

    @property
    def root(self) -> RecursionNode:
        if self.root_id is None:
            raise ValueError("tree has no root")
        return self.nodes[self.root_id]

    @property
    def depth(self) -> int:
        """Maximum node depth; 0 for a single-node tree."""
        return max((node.depth for node in self.nodes.values()), default=0)

    @property
    def leaves(self) -> list[RecursionNode]:
        return [n for n in self if n.is_leaf]

    def children_of(self, node_id: int) -> list[RecursionNode]:
        return [self.nodes[c] for c in self.nodes[node_id].children]

    def preorder(self, start: Optional[int] = None) -> Iterator[RecursionNode]:
        """Depth-first, left subproblem before right."""
        start = self.root_id if start is None else start
        if start is None:
            return
        stack = [start]
        while stack:
            node = self.nodes[stack.pop()]
            yield node
            stack.extend(reversed(node.children))

    def summary(self) -> dict[str, Any]:
        """Aggregate statistics -- the numbers the experiments compare."""
        by_reason: dict[str, int] = {}
        for node in self:
            if node.verdict is not None:
                by_reason[node.verdict.reason] = by_reason.get(node.verdict.reason, 0) + 1
        return {
            "instance": self.instance,
            "variant": self.variant,
            "pivot_rule": self.pivot_rule,
            "dual": self.dual,
            "nodes": len(self.nodes),
            "leaves": len(self.leaves),
            "depth": self.depth,
            "root_epsilon": self.root.epsilon if self.nodes else None,
            "verdicts": by_reason,
        }

    # ------------------------------------------------------------------
    def to_json(self) -> dict[str, Any]:
        return {
            "instance": self.instance,
            "variant": self.variant,
            "pivot_rule": self.pivot_rule,
            "dual": self.dual,
            "verdict": None if self.verdict is None else self.verdict.to_json(),
            "root_id": self.root_id,
            "summary": self.summary(),
            "nodes": [n.to_json() for n in self],
        }

    def save(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_json(), indent=2), encoding="utf-8")
        return path

    @classmethod
    def from_json(cls, d: dict[str, Any]) -> "RecursionTree":
        tree = cls(
            instance=d.get("instance", ""),
            variant=d.get("variant", ""),
            pivot_rule=d.get("pivot_rule", ""),
            dual=d.get("dual"),
            root_id=d.get("root_id"),
        )
        for nd in d["nodes"]:
            node = RecursionNode.from_json(nd)
            tree.nodes[node.node_id] = node
        return tree

    @classmethod
    def load(cls, path: str | Path) -> "RecursionTree":
        return cls.from_json(json.loads(Path(path).read_text(encoding="utf-8")))
