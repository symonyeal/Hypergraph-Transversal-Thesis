"""The problem-instance library, and the generators that build families of it.

One JSON file per instance under ``data/instances/``. Data rather than literals
in notebook cells: defined once, diffable, carrying its provenance, and readable
by SageMath or MATLAB without importing this package. The format is documented
in ``data/README.md``.

The ``expected`` block is not hand-written. :func:`refresh_expected` recomputes
it -- duality from the oracle, ``epsilon`` as an exact fraction, node counts for
both algorithms and every variant -- and the tests assert the committed values
still hold. So it is a regression baseline, not a claim: a moved node count
shows up as a diff.

Provenance is ``verbatim`` (transcribed as published, even where demonstrably
wrong), ``corrected`` (a stated, minimal repair, with the evidence in ``notes``)
or ``derived`` (built here, from a generator or a run log). Nothing is silently
rewritten: a withdrawn instance goes to ``_archive/`` and is still readable
through :func:`load_archived`.
"""

from __future__ import annotations

import itertools
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Sequence

from .hypergraph import Hypergraph, bits
from .selfdual import FANO_LINES

__all__ = [
    "Instance",
    "instance_dir",
    "data_dir",
    "load_data",
    "archive_root",
    "load",
    "load_all",
    "load_archived",
    "list_ids",
    "sdfp",
    "self_dual_fano",
    "self_dualise",
    "witt_11",
    "matching",
    "threshold",
    "refresh_expected",
]


def instance_dir() -> Path:
    """Location of the instance JSON files.

    Resolved relative to the installed package so the library works from any
    working directory: ``src/fka/instances.py`` -> ``<repo>/data/instances``.
    """
    return Path(__file__).resolve().parents[2] / "data" / "instances"


@dataclass(slots=True)
class Instance:
    """One named ``(cnf, dnf)`` problem instance."""

    id: str
    name: str
    n_vertices: int
    cnf_edges: list[list[int]]
    dnf_edges: list[list[int]]
    family: str = ""
    source: str = ""
    provenance: str = "verbatim"
    notes: str = ""
    expected: dict[str, Any] = field(default_factory=dict)
    path: Optional[Path] = None

    @property
    def cnf(self) -> Hypergraph:
        return Hypergraph.from_sets(self.n_vertices, self.cnf_edges, one_indexed=True)

    @property
    def dnf(self) -> Hypergraph:
        return Hypergraph.from_sets(self.n_vertices, self.dnf_edges, one_indexed=True)

    def to_json(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "family": self.family,
            "source": self.source,
            "provenance": self.provenance,
            "n_vertices": self.n_vertices,
            "cnf": self.cnf_edges,
            "dnf": self.dnf_edges,
            "notes": self.notes,
            "expected": self.expected,
        }

    @classmethod
    def from_json(cls, d: dict[str, Any], path: Optional[Path] = None) -> "Instance":
        return cls(
            id=d["id"],
            name=d.get("name", d["id"]),
            n_vertices=d["n_vertices"],
            cnf_edges=[list(t) for t in d["cnf"]],
            dnf_edges=[list(t) for t in d["dnf"]],
            family=d.get("family", ""),
            source=d.get("source", ""),
            provenance=d.get("provenance", "verbatim"),
            notes=d.get("notes", ""),
            expected=d.get("expected", {}),
            path=path,
        )

    def save(self, path: Optional[Path] = None) -> Path:
        target = Path(path or self.path or (instance_dir() / f"{self.id}.json"))
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(self.to_json(), indent=2) + "\n", encoding="utf-8")
        self.path = target
        return target


def data_dir() -> Path:
    """The one folder all research inputs live under.

    ``instances/`` the ``(cnf, dnf)`` library, ``reference/`` the MATLAB FK-B
    authors' own recorded test vectors, ``baselines/`` the node counts and split
    sequences printed in the thesis. Superseded inputs move to ``_archive/`` and
    are read back through :func:`load_archived`.
    """
    return Path(__file__).resolve().parents[2] / "data"


def load_data(*parts: str) -> Any:
    """Read one JSON file from :func:`data_dir`, e.g. ``load_data("baselines",
    "published-trees.json")``."""
    path = data_dir().joinpath(*parts)
    if not path.exists():
        raise FileNotFoundError(f"no data file at {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def archive_root() -> Path:
    """Where withdrawn instances live, one dated subdirectory per withdrawal."""
    return Path(__file__).resolve().parents[2] / "_archive"


def load_archived(instance_id: str) -> Instance:
    """Load an instance that has been withdrawn from the live library.

    Withdrawn instances are kept for traceability -- principally the two thesis
    instances that are not transversal pairs as printed -- and are deliberately
    absent from :func:`list_ids`, so no experiment or baseline runs against one
    by accident. Reading one is an explicit act, and the withdrawal's own
    ``README.md`` states what it was and why.
    """
    matches = sorted(archive_root().glob(f"*/{instance_id}.json"))
    if not matches:
        available = ", ".join(
            sorted(p.stem for p in archive_root().glob("*/*.json"))
        ) or "(none)"
        raise FileNotFoundError(
            f"no archived instance {instance_id!r} under {archive_root()}. "
            f"Available: {available}"
        )
    path = matches[0]
    return Instance.from_json(json.loads(path.read_text(encoding="utf-8")), path)


def list_ids() -> list[str]:
    return sorted(p.stem for p in instance_dir().glob("*.json"))


def load(instance_id: str) -> Instance:
    path = instance_dir() / f"{instance_id}.json"
    if not path.exists():
        available = ", ".join(list_ids()) or "(none)"
        raise FileNotFoundError(
            f"no instance {instance_id!r} in {instance_dir()}. Available: {available}"
        )
    return Instance.from_json(json.loads(path.read_text(encoding="utf-8")), path)


def load_all() -> list[Instance]:
    return [load(i) for i in list_ids()]


# ----------------------------------------------------------------------
# generated families
# ----------------------------------------------------------------------
def sdfp(k: int = 1) -> Hypergraph:
    """``k`` disjoint copies of the Fano lines (thesis Def. 5.2.1, p.50).

    Only ``k = 1`` is self-dual; for ``k >= 2`` the copies are disjoint and
    ``Tr`` is the much larger family the thesis gives a size formula for. This
    is *not* :func:`self_dual_fano`, which dualises this construction.
    """
    if k < 1:
        raise ValueError(f"k must be at least 1, got {k}")
    n = 7 * k
    edges = [
        [v + 7 * copy for v in line] for copy in range(k) for line in FANO_LINES
    ]
    return Hypergraph.from_sets(n, edges, one_indexed=True)


def self_dual_fano(k: int = 2) -> Hypergraph:
    """``SDFP(k)``: the self-dualisation of ``k`` disjoint Fano planes.

    The reference's ``SDFP16_CNF_DNF.mat`` (``k = 2``) and ``SDFP23`` (``k = 3``).
    On ``7k + 2`` variables, with ``a``, ``b`` the two extra: ``{a,b}``, each of
    the ``7k`` lines plus ``b``, and each of the ``7^k`` cross-copy line choices
    plus ``a``. Self-dual, so ``(E, E)`` is a dual pair with no free parameters
    -- the benchmark shape, where the answer is yes and every branch must be
    explored to establish it.

    Variables use this project's :data:`FANO_LINES` labelling, not the reference
    snapshot's. The two are isomorphic, so node counts can differ from the MATLAB
    run by however much the numbering moves the split choice.
    """
    return self_dualise([list(line) for line in FANO_LINES], 7, k)


def self_dualise(base_edges: Sequence[Sequence[int]], n_base: int, k: int = 1) -> Hypergraph:
    """Self-dualise ``k`` disjoint copies of a 1-indexed base family.

    The construction :func:`self_dual_fano` names, with the base family a
    parameter: on ``n_base*k + 2`` vertices, writing ``a`` and ``b`` for the two
    extra ones, the edges are ``{a,b}``; every base edge of every copy plus
    ``b``; and every choice of one base edge per copy plus ``a``. The result is
    self-transversal whatever the base is, which is what makes it the standard
    way to turn a family into a dualisation benchmark: the answer is always
    "dual", so no branch can be closed early by a counterexample.

    ``|E| = |base|^k + k|base| + 1``, so ``k`` is limited by the base: two
    copies of the Fano plane's 7 lines give 64 edges, two of the 66 blocks of
    :func:`witt_11` give 4489.
    """
    if k < 1:
        raise ValueError(f"k must be at least 1, got {k}")
    n = n_base * k + 2
    a, b = n - 1, n  # 1-indexed
    edges = [[a, b]]
    for copy in range(k):
        for e in base_edges:
            edges.append(sorted(v + n_base * copy for v in e) + [b])
    for combo in itertools.product(base_edges, repeat=k):
        picked = sorted(v + n_base * c for c, e in enumerate(combo) for v in e)
        edges.append(picked + [a])
    return Hypergraph.from_sets(n, edges, one_indexed=True)


#: ``g(x) = -1 + x^2 - x^3 + x^4 + x^5``, one of the two degree-5 factors of
#: ``x^11 - 1`` over ``GF(3)``. The cyclic code it generates is the perfect
#: ``[11, 6, 5]`` ternary Golay code.
GOLAY3_GENERATOR: tuple[int, ...] = (2, 0, 1, 2, 1, 1)


def witt_11() -> Hypergraph:
    """``W11``: the 66 blocks of the Steiner system ``S(4, 5, 11)``.

    Built as the supports of the 132 minimum-weight words of the perfect
    ``[11, 6, 5]`` ternary Golay code, which is the standard construction of
    the design: every 4-subset of the 11 points lies in exactly one block,
    ``Aut = M11`` of order 7920, every vertex has degree 30.

    Like the Fano plane -- which is the same construction on the perfect
    ``[7, 4, 3]`` Hamming code -- it is *self-transversal*: ``Tr(W11) = W11``.
    Unlike the Fano plane it is hard for FK-B as well as FK-A; see
    ``docs/hard-cases.md``. The third perfect code, the binary Golay
    ``[23, 12, 7]``, gives ``S(4, 7, 23)``, which is intersecting but **not**
    self-transversal, so the family stops here.
    """
    n = 11
    degree = len(GOLAY3_GENERATOR) - 1
    blocks = []
    for message in itertools.product(range(3), repeat=n - degree):
        word = [0] * n
        for i, coefficient in enumerate(message):
            if coefficient:
                for j, g in enumerate(GOLAY3_GENERATOR):
                    word[(i + j) % n] = (word[(i + j) % n] + coefficient * g) % 3
        support = [i for i, x in enumerate(word) if x]
        if len(support) == 5:
            blocks.append(bits(support))
    H = Hypergraph.from_bitsets(n, blocks)
    if len(H) != 66:
        raise AssertionError(f"S(4,5,11) must have 66 blocks, built {len(H)}")
    return H


def matching(v: int) -> Hypergraph:
    """``M(v)``: the perfect matching on ``v`` vertices (StandardProblems.m).

    ``v/2`` disjoint pairs ``{1,2}, {3,4}, ...``. Its transversal is the
    ``2^(v/2)`` sets that pick one endpoint per pair, so it is the standard
    example of output-exponential dualisation.
    """
    if v < 2 or v % 2:
        raise ValueError(f"v must be even and at least 2, got {v}")
    return Hypergraph.from_sets(
        v, [[i, i + 1] for i in range(1, v, 2)], one_indexed=True
    )


def threshold(v: int) -> Hypergraph:
    """``TH(v)``: the threshold graph on ``v`` vertices (StandardProblems.m).

    All pairs ``{i, j}`` with ``i < j`` and ``j`` even.
    """
    if v < 2 or v % 2:
        raise ValueError(f"v must be even and at least 2, got {v}")
    return Hypergraph.from_sets(
        v,
        [[i, j] for j in range(2, v + 1, 2) for i in range(1, j)],
        one_indexed=True,
    )


# ----------------------------------------------------------------------
# expected-value maintenance
# ----------------------------------------------------------------------
def refresh_expected(inst: Instance, *, node_counts: bool = True) -> dict[str, Any]:
    """Recompute the ``expected`` block for ``inst``.

    Duality comes from the brute-force oracle, never from either algorithm --
    the whole point of the baseline is to check them against something
    independent. Node counts are recorded per algorithm and per variant, at
    each one's default split rule.
    """
    from .algorithm import eps_exact, fk_a
    from .transversal import is_dual_oracle

    # ``fkb`` ships alongside this package and imports the model from it. The
    # import is function-local so the two never form an import cycle.
    from fkb.algorithm import fk_b

    cnf, dnf = inst.cnf, inst.dnf
    expected: dict[str, Any] = {
        "dual": is_dual_oracle(cnf, dnf),
        "epsilon": str(eps_exact(cnf, dnf)),
        "nC": len(cnf),
        "nD": len(dnf),
    }
    if node_counts:
        expected["fka"] = {
            variant: _tree_shape(fk_a(cnf, dnf, variant=variant, instance=inst.id))
            for variant in ("faithful", "modified")
        }
        expected["fkb"] = {
            variant: _tree_shape(fk_b(cnf, dnf, variant=variant, instance=inst.id))
            for variant in ("faithful", "multiple")
        }
    inst.expected = expected
    return expected


def _tree_shape(tree: Any) -> dict[str, int]:
    return {"nodes": len(tree), "depth": tree.depth}
