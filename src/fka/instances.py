"""The problem-instance library.

Each instance lives in one JSON file under ``data/instances/``. Keeping them as
data rather than as literals pasted into notebook cells means an instance is
defined once, is diffable, carries its provenance, and can be read by SageMath,
MATLAB or anything else without importing this package.

File format::

    {
      "id":          "f2g2",
      "name":        "f2 / g2",
      "family":      "read-once (GK-family)",
      "source":      "Islam (2024) MSc thesis, Section 5.2.1, p.48",
      "provenance":  "verbatim" | "corrected" | "derived",
      "n_vertices":  8,
      "G":           [[1, 3], [1, 4], ...],       1-indexed vertex lists
      "H":           [[1, 2, 7, 8], ...],
      "notes":       "free text",
      "expected":    { ... }                       machine-maintained, see below
    }

The ``expected`` block is *not* hand-written. ``fka instances refresh``
recomputes it from the instance itself -- duality against the brute-force
oracle, the frequency ``epsilon`` as an exact fraction, and the FK-A node counts
for each variant -- and the test suite then asserts the committed values still
hold. That makes it a regression baseline rather than a claim: if a change to
the algorithm moves a node count, the diff says so.

Provenance
----------
``verbatim``
    Transcribed exactly as published, even where that is demonstrably wrong.
``corrected``
    A published instance with a stated, minimal repair. Two exist, both because
    the pair as printed is not a transversal pair at all; each carries the
    correction and the evidence for it in ``notes``. The verbatim form is kept
    alongside, so nothing is silently rewritten.
``derived``
    Constructed here (from a generator, or recovered from a run log).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterator, Optional

from .hypergraph import Hypergraph

__all__ = ["Instance", "instance_dir", "load", "load_all", "list_ids", "sdfp", "refresh_expected"]


def instance_dir() -> Path:
    """Location of the instance JSON files.

    Resolved relative to the installed package so the library works from any
    working directory: ``src/fka/instances.py`` -> ``<repo>/data/instances``.
    """
    return Path(__file__).resolve().parents[2] / "data" / "instances"


@dataclass(slots=True)
class Instance:
    """One named ``(G, H)`` problem instance."""

    id: str
    name: str
    n_vertices: int
    G_edges: list[list[int]]
    H_edges: list[list[int]]
    family: str = ""
    source: str = ""
    provenance: str = "verbatim"
    notes: str = ""
    expected: dict[str, Any] = field(default_factory=dict)
    path: Optional[Path] = None

    @property
    def G(self) -> Hypergraph:
        return Hypergraph.from_sets(self.n_vertices, self.G_edges, one_indexed=True)

    @property
    def H(self) -> Hypergraph:
        return Hypergraph.from_sets(self.n_vertices, self.H_edges, one_indexed=True)

    def to_json(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "family": self.family,
            "source": self.source,
            "provenance": self.provenance,
            "n_vertices": self.n_vertices,
            "G": self.G_edges,
            "H": self.H_edges,
            "notes": self.notes,
            "expected": self.expected,
        }

    @classmethod
    def from_json(cls, d: dict[str, Any], path: Optional[Path] = None) -> "Instance":
        return cls(
            id=d["id"],
            name=d.get("name", d["id"]),
            n_vertices=d["n_vertices"],
            G_edges=[list(e) for e in d["G"]],
            H_edges=[list(e) for e in d["H"]],
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
#: The seven lines of the Fano plane, 1-indexed (thesis Definition 5.2.1, p.50).
FANO_LINES: tuple[tuple[int, ...], ...] = (
    (1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 7),
    (1, 5, 6), (2, 6, 7), (1, 3, 7),
)


def sdfp(k: int = 1) -> Hypergraph:
    """The self-dual Fano-plane hypergraph ``SDFP`` on ``k`` disjoint copies.

    Thesis Definition 5.2.1 (p.50): ``E`` is ``k`` disjoint copies of the Fano
    plane's line set, giving ``7k`` edges on ``7k`` vertices.

    ``k = 1`` is the self-transversal base case shown in the thesis
    visualisations. For ``k >= 2`` the copies are disjoint, so the hypergraph is
    *not* self-transversal -- its transversal is the much larger family the
    thesis quotes a size formula for; use :func:`fka.transversal.transversal` to
    build it, and expect that to be expensive.
    """
    if k < 1:
        raise ValueError(f"k must be at least 1, got {k}")
    n = 7 * k
    edges = [
        [v + 7 * copy for v in line] for copy in range(k) for line in FANO_LINES
    ]
    return Hypergraph.from_sets(n, edges, one_indexed=True)


# ----------------------------------------------------------------------
# expected-value maintenance
# ----------------------------------------------------------------------
def exact_epsilon(G: Hypergraph, H: Hypergraph) -> Fraction:
    """``epsilon(G, H)`` as an exact fraction, for comparing with the thesis.

    :func:`fka.algorithm.epsilon` returns a float, which is what the algorithm
    needs; the thesis quotes values like ``5/11`` and ``3/7``, so the instance
    baselines store the fraction.
    """
    n = max(G.n, H.n)
    best = Fraction(0)
    for v in range(n):
        if len(G):
            best = max(best, Fraction(G.degree(v), len(G)))
        if len(H):
            best = max(best, Fraction(H.degree(v), len(H)))
    return best


def refresh_expected(inst: Instance, *, node_counts: bool = True) -> dict[str, Any]:
    """Recompute the ``expected`` block for ``inst``.

    Duality comes from the brute-force oracle, never from FK-A -- the whole
    point of the baseline is to check FK-A against something independent.
    """
    from .algorithm import fk_a
    from .transversal import is_dual_oracle

    G, H = inst.G, inst.H
    expected: dict[str, Any] = {
        "dual": is_dual_oracle(G, H),
        "epsilon": str(exact_epsilon(G, H)),
        "n_edges_G": len(G),
        "n_edges_H": len(H),
    }
    if node_counts:
        counts: dict[str, dict[str, int]] = {}
        for variant in ("faithful", "modified"):
            tree = fk_a(G, H, variant=variant, instance=inst.id)
            counts[variant] = {"nodes": len(tree), "depth": tree.depth}
        expected["fka"] = counts
    inst.expected = expected
    return expected
