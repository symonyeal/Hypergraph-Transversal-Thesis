"""The recursion tree, shared by FK-A and FK-B.

One :class:`Node` per recursive call. The thesis treats this tree as the primary
object of study -- "The performance of FK-A is determined by the number of
subproblems it generates, and this is reflected in the number of nodes of the
recursion tree" (p.12) -- and Chapter 5 annotates each node with the
automorphism group and properties of its inputs.

Algorithm-neutral: FK-A builds binary ``L``/``R`` trees, FK-B trees whose nodes
may have many children, and :attr:`Tree.algorithm` says which.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator, Optional

from .hypergraph import Hypergraph, verts

__all__ = ["Node", "Tree", "Verdict"]


@dataclass(frozen=True, slots=True)
class Verdict:
    """Why a node concluded what it concluded.

    ``reason`` is a stable code, so experiment output groups without
    string-matching prose. FK-A emits ``dual``, ``trivial``, ``single_edge``,
    ``cond_i_support`` / ``cond_ii_size`` / ``cond_iii_disjoint`` /
    ``cond_iv_frequency`` and ``child_failed``; FK-B numbers its own conditions
    1-3, so the two sets never collide.

    ``CA`` is a conflicting assignment: a bitset with ``C(S) != D(S)``. ``None``
    never means "none exists", only that none was built in this orientation --
    a failed condition then keeps the thesis' structural ``witness_edges``.
    ``CAs`` holds all of them under FK-B's ``multiple`` variant (MFK_B.m), and
    ``CAs[0] == CA`` when set.
    """

    dual: bool
    reason: str
    detail: str = ""
    CA: Optional[int] = None
    witness_edges: tuple[int, ...] = ()
    CAs: tuple[int, ...] = ()

    def to_json(self) -> dict[str, Any]:
        out = {
            "dual": self.dual,
            "reason": self.reason,
            "detail": self.detail,
            "CA": None if self.CA is None else [v + 1 for v in verts(self.CA)],
            "witness_edges": [
                [v + 1 for v in verts(t)] for t in self.witness_edges
            ],
        }
        if self.CAs:
            out["CAs"] = [[v + 1 for v in verts(s)] for s in self.CAs]
        return out


@dataclass(slots=True)
class Node:
    """One recursive call.

    ``cnf_in``/``dnf_in`` are the inputs as handed down by the parent;
    ``cnf``/``dnf`` are the same pair after ``Irredundant``. The difference is
    what the thesis' node dumps label "Non Sperner Edges".
    """

    node_id: int
    parent_id: Optional[int]
    path: tuple[str, ...]
    cnf_in: Hypergraph
    dnf_in: Hypergraph
    cnf: Hypergraph
    dnf: Hypergraph
    cut_C: tuple[int, ...] = ()
    cut_D: tuple[int, ...] = ()
    split_var: Optional[int] = None
    split_var_frequency: Optional[float] = None
    epsilon: Optional[float] = None
    #: Which split shape this node took. FK-A has only one and leaves it empty;
    #: FK-B records ``"mu_D"``, ``"mu_C"`` or ``"split"``, which is what decides
    #: how many children it has.
    branch: str = ""
    children: list[int] = field(default_factory=list)
    verdict: Optional[Verdict] = None
    #: Filled in by :mod:`fka.analysis`.
    analysis: dict[str, Any] = field(default_factory=dict)

    @property
    def depth(self) -> int:
        return len(self.path)

    @property
    def is_leaf(self) -> bool:
        return not self.children

    @property
    def step(self) -> str:
        """The step that produced this node, or ``"root"``.

        FK-A uses ``L``/``R``; FK-B uses ``x=0``, ``x=1``, ``c1``.. and ``m1``..
        Display only -- renderers needing a fixed vocabulary use
        :meth:`step_style`.
        """
        return self.path[-1] if self.path else "root"

    def step_style(self) -> str:
        """``solid``, ``dashed`` or ``dotted`` for this node's incoming edge.

        The first subproblem of a split is solid and the second dashed in both
        algorithms; FK-B's per-clause and per-monomial subproblems, which have
        no FK-A counterpart, are dotted.
        """
        s = self.step
        if s in ("L", "x=0"):
            return "solid"
        if s in ("R", "x=1"):
            return "dashed"
        return "dotted"

    def split_var_label(self) -> str:
        """The split variable in the thesis' 1-indexed ``v_i`` notation."""
        return "" if self.split_var is None else f"v{self.split_var + 1}"

    def to_json(self) -> dict[str, Any]:
        return {
            "node_id": self.node_id,
            "parent_id": self.parent_id,
            "path": list(self.path),
            "depth": self.depth,
            "cnf": self.cnf.to_sets(one_indexed=True),
            "dnf": self.dnf.to_sets(one_indexed=True),
            "cnf_in": self.cnf_in.to_sets(one_indexed=True),
            "dnf_in": self.dnf_in.to_sets(one_indexed=True),
            "n": self.cnf.n,
            "cut_C": [[v + 1 for v in verts(t)] for t in self.cut_C],
            "cut_D": [[v + 1 for v in verts(t)] for t in self.cut_D],
            "split_var": None if self.split_var is None else self.split_var + 1,
            "split_var_frequency": self.split_var_frequency,
            "epsilon": self.epsilon,
            "branch": self.branch,
            "children": list(self.children),
            "verdict": None if self.verdict is None else self.verdict.to_json(),
            "analysis": self.analysis,
        }

    @classmethod
    def from_json(cls, d: dict[str, Any]) -> "Node":
        n = d["n"]

        def mbf(key: str) -> Hypergraph:
            return Hypergraph.from_sets(n, d[key], one_indexed=True)

        def terms(key: str) -> tuple[int, ...]:
            return tuple(
                sum(1 << (v - 1) for v in term) for term in d[key]
            )

        verdict = None
        if d.get("verdict"):
            v = d["verdict"]
            CA = v.get("CA")
            verdict = Verdict(
                dual=v["dual"],
                reason=v["reason"],
                detail=v.get("detail", ""),
                CA=None if CA is None else sum(1 << (x - 1) for x in CA),
                witness_edges=tuple(
                    sum(1 << (x - 1) for x in t) for t in v.get("witness_edges", [])
                ),
                CAs=tuple(
                    sum(1 << (x - 1) for x in s) for s in v.get("CAs", [])
                ),
            )
        return cls(
            node_id=d["node_id"],
            parent_id=d["parent_id"],
            path=tuple(d["path"]),
            cnf_in=mbf("cnf_in"),
            dnf_in=mbf("dnf_in"),
            cnf=mbf("cnf"),
            dnf=mbf("dnf"),
            cut_C=terms("cut_C"),
            cut_D=terms("cut_D"),
            split_var=None if d["split_var"] is None else d["split_var"] - 1,
            split_var_frequency=d.get("split_var_frequency"),
            epsilon=d.get("epsilon"),
            branch=d.get("branch", ""),
            children=list(d.get("children", [])),
            verdict=verdict,
            analysis=d.get("analysis", {}),
        )


@dataclass(slots=True)
class Tree:
    """The full trace of one FK-A or FK-B run."""

    nodes: dict[int, Node] = field(default_factory=dict)
    root_id: Optional[int] = None
    instance: str = ""
    algorithm: str = "FK-A"
    variant: str = ""
    split_rule: str = ""
    dual: Optional[bool] = None
    verdict: Optional[Verdict] = None

    # ------------------------------------------------------------------
    def add(self, node: Node) -> Node:
        self.nodes[node.node_id] = node
        if node.parent_id is None:
            self.root_id = node.node_id
        else:
            self.nodes[node.parent_id].children.append(node.node_id)
        return node

    def __len__(self) -> int:
        return len(self.nodes)

    def __getitem__(self, node_id: int) -> Node:
        return self.nodes[node_id]

    def __iter__(self) -> Iterator[Node]:
        """Iterate in node-id order, which is the order the thesis numbers them."""
        for k in sorted(self.nodes):
            yield self.nodes[k]

    @property
    def root(self) -> Node:
        if self.root_id is None:
            raise ValueError("tree has no root")
        return self.nodes[self.root_id]

    @property
    def depth(self) -> int:
        return max((node.depth for node in self.nodes.values()), default=0)

    @property
    def leaves(self) -> list[Node]:
        return [n for n in self if n.is_leaf]

    def children_of(self, node_id: int) -> list[Node]:
        return [self.nodes[c] for c in self.nodes[node_id].children]

    def preorder(self, start: Optional[int] = None) -> Iterator[Node]:
        """Depth-first, first subproblem before second."""
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
        by_branch: dict[str, int] = {}
        for node in self:
            if node.verdict is not None:
                by_reason[node.verdict.reason] = by_reason.get(node.verdict.reason, 0) + 1
            if node.branch:
                by_branch[node.branch] = by_branch.get(node.branch, 0) + 1
        out = {
            "instance": self.instance,
            "algorithm": self.algorithm,
            "variant": self.variant,
            "split_rule": self.split_rule,
            "dual": self.dual,
            "nodes": len(self.nodes),
            "leaves": len(self.leaves),
            "depth": self.depth,
            "root_epsilon": self.root.epsilon if self.nodes else None,
            "verdicts": by_reason,
        }
        # Only FK-B has more than one branch shape, so this stays out of FK-A's
        # summaries rather than sitting there permanently empty.
        if by_branch:
            out["branches"] = by_branch
        return out

    # ------------------------------------------------------------------
    def to_json(self) -> dict[str, Any]:
        return {
            "instance": self.instance,
            "algorithm": self.algorithm,
            "variant": self.variant,
            "split_rule": self.split_rule,
            "dual": self.dual,
            "verdict": None if self.verdict is None else self.verdict.to_json(),
            "root_id": self.root_id,
            "summary": self.summary(),
            "nodes": [n.to_json() for n in self],
        }

    def save(self, path: str | Path, *, indent: Optional[int] = None) -> Path:
        """Write the tree as JSON.

        Compact by default. A tree is regenerated wholesale, never hand-edited,
        and pretty-printing one costs about 2.7x the bytes -- 7.4 MB against
        2.7 MB for ``sdfp-sd-2`` under FK-A. Pass ``indent=2`` when a specific
        small tree is going to be read by eye.
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        text = json.dumps(
            self.to_json(),
            indent=indent,
            separators=None if indent else (",", ":"),
        )
        path.write_text(text, encoding="utf-8")
        return path

    @classmethod
    def from_json(cls, d: dict[str, Any]) -> "Tree":
        tree = cls(
            instance=d.get("instance", ""),
            algorithm=d.get("algorithm", "FK-A"),
            variant=d.get("variant", ""),
            split_rule=d.get("split_rule", ""),
            dual=d.get("dual"),
            root_id=d.get("root_id"),
        )
        for nd in d["nodes"]:
            node = Node.from_json(nd)
            tree.nodes[node.node_id] = node
        return tree

    @classmethod
    def load(cls, path: str | Path) -> "Tree":
        return cls.from_json(json.loads(Path(path).read_text(encoding="utf-8")))
