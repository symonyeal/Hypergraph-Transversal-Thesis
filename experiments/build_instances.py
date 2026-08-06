"""Author the instance library under ``data/instances/``.

Run once to create the files, and again after adding an instance below. The
``expected`` block of each file is computed here, not typed: duality comes from
the brute-force oracle in :mod:`fka.transversal`, ``epsilon`` from the exact
fraction, and node counts from both FK-A variants.

    python experiments/build_instances.py

Adding an instance means adding one entry to ``DEFINITIONS`` and re-running.
Editing an existing instance's ``cnf``/``dnf`` is a research decision -- say why in
``notes`` and set ``provenance`` accordingly.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fka.instances import (  # noqa: E402
    Instance,
    archive_root,
    instance_dir,
    refresh_expected,
)
from fka.selfdual import FANO_LINES  # noqa: E402

THESIS = "Islam, 'Tracing the Effects of Symmetry in Hypergraphs On the Fredman-Khachiyan Algorithm', MSc thesis, SFU 2024"

DEFINITIONS: list[dict] = [
    {
        "id": "f2g2",
        "name": "f2 / g2",
        "family": "read-once, GK-family (Hertel's Mark's class)",
        "source": f"{THESIS}, Section 5.2.1, p.48",
        "provenance": "verbatim",
        "n_vertices": 8,
        "cnf": [[1, 3], [1, 4], [2, 3], [2, 4], [5, 7], [5, 8], [6, 7], [6, 8]],
        "dnf": [[1, 2, 7, 8], [3, 4, 5, 6], [1, 2, 5, 6], [3, 4, 7, 8]],
        "notes": (
            "The worked example of thesis Figure 5.1. Its primal graph is two "
            "disjoint C4s, so Aut is D4 wr S2 of order 128 -- the thesis writes "
            "this as (((C2 x C2 x C2 x C2) : C2) : C2) : C2 and notes the wreath "
            "form on p.53. Reproduces the published recursion tree exactly: "
            "25 nodes under the modified variant, with split vertices "
            "v1 v2 v5 v6 - - - v3 - v4 v5 v6 - - - - v3 - v4 v5 v6 - - - -."
        ),
    },
    {
        "id": "fano",
        "name": "Fano plane (SDFP, k=1)",
        "family": "self-transversal projective plane",
        "source": f"{THESIS}, Definition 5.2.1, p.50",
        "provenance": "verbatim",
        "n_vertices": 7,
        "cnf": [list(line) for line in FANO_LINES],
        "dnf": [list(line) for line in FANO_LINES],
        "notes": (
            "Self-transversal: Tr(Fano) = Fano. Aut is PSL(3,2) of order 168. "
            "The thesis uses it to show the frequency bound is tight at n = 7, "
            "with every vertex at epsilon = 3/7 (Observation 5.2.3, p.54). "
            "Not conformal; its primal graph is the complete graph K7. "
            "Reproduces thesis Figure 5.2: 35 nodes under the modified variant."
        ),
    },
    {
        "id": "k3",
        "name": "K3",
        "family": "self-transversal projective plane (smallest)",
        "source": f"{THESIS}, Observation 5.2.3, p.54",
        "provenance": "verbatim",
        "n_vertices": 3,
        "cnf": [[1, 2], [2, 3], [1, 3]],
        "dnf": [[1, 2], [2, 3], [1, 3]],
        "notes": (
            "The complete 2-uniform graph on 3 vertices, self-transversal, all "
            "frequencies 2/3. The smallest finite projective plane and the "
            "n = 3 tightness witness for the frequency bound."
        ),
    },
    {
        "id": "6a-verbatim",
        "name": "6-A (as published)",
        "family": "reinforcement-learning search results of [SS24]",
        "source": f"{THESIS}, Section 5.2.1, p.49",
        "provenance": "verbatim",
        "n_vertices": 6,
        "cnf": [[1, 5, 6], [4, 5, 6], [1, 2], [2, 3, 4]],
        "dnf": [[2, 4, 6], [2, 5], [1, 4], [1, 3, 5]],
        "notes": (
            "NOT a transversal pair as printed. Tr(C) has five edges -- it "
            "additionally contains {2,6} and {1,3,6} -- while D contains {2,4,6}, "
            "which is a proper superset of {2,6} and so cannot be a minimal "
            "transversal. Kept verbatim for traceability; see instance '6a' for "
            "the repair."
        ),
    },
    {
        "id": "6a",
        "name": "6-A (corrected)",
        "family": "reinforcement-learning search results of [SS24]",
        "source": f"{THESIS}, Section 5.2.1, p.49, with a corrected hyperedge",
        "provenance": "corrected",
        "n_vertices": 6,
        "cnf": [[1, 5, 6], [4, 5], [1, 2], [2, 3, 4]],
        "dnf": [[2, 4, 6], [2, 5], [1, 4], [1, 3, 5]],
        "notes": (
            "Repair of '6a-verbatim': C's edge {v4,v5,v6} should be {v4,v5}. "
            "Evidence: Tr(D) as published equals exactly this C, the repaired "
            "pair is mutually transversal in both directions, and it reproduces "
            "the thesis' own reported values for the instance -- epsilon = 1/2 "
            "(p.55) and Aut = C2 x C2 at the root (Figure 5.4a). The verbatim "
            "pair reproduces neither. Note the thesis also quotes generators "
            "<(v1 v4),(v5 v6)> on p.18; (v1 v4) is not an automorphism of either "
            "form of C, though the group is C2 x C2 as stated."
        ),
    },
    {
        "id": "6b",
        "name": "6-B",
        "family": "reinforcement-learning search results of [SS24]",
        "source": f"{THESIS}, Section 5.2.1, p.49",
        "provenance": "verbatim",
        "n_vertices": 6,
        "cnf": [
            [1, 3, 5], [1, 2, 6], [1, 2, 3], [1, 5, 6],
            [3, 4, 5], [2, 3, 4], [4, 5, 6], [2, 4, 6],
        ],
        "dnf": [[2, 5], [1, 4], [3, 6]],
        "notes": (
            "D is a perfect matching on the three twin pairs. Aut(C) = Aut(D) = "
            "C2 x S4 of order 48, matching thesis Figure 5.4b. D is read-once "
            "and alpha-acyclic; C is read-once but not alpha-acyclic (no vertex "
            "lies in a single edge, so GYO reduction stalls immediately)."
        ),
    },
    {
        "id": "6c",
        "name": "6-C",
        "family": "reinforcement-learning search results of [SS24]",
        "source": f"{THESIS}, Section 5.2.1, p.49",
        "provenance": "verbatim",
        "n_vertices": 6,
        "cnf": [
            [3, 4], [1, 5], [2, 3], [4, 5],
            [5, 6], [3, 6], [2, 5], [1, 3],
        ],
        "dnf": [[3, 5], [1, 2, 4, 6]],
        "notes": (
            "C is 2-uniform, so it is literally a graph. Aut = C2 x S4 of order "
            "48, matching thesis Figure 5.4c."
        ),
    },
    {
        "id": "6d",
        "name": "6-D",
        "family": "reinforcement-learning search results of [SS24]",
        "source": f"{THESIS}, Section 5.2.1, p.49",
        "provenance": "verbatim",
        "n_vertices": 6,
        "cnf": [[1, 2, 6], [2, 4], [3, 4, 5], [1, 3, 5, 6]],
        "dnf": [[2, 3], [2, 5], [1, 4], [4, 6]],
        "notes": "Aut(C) = D4 of order 8, matching thesis Figure 5.4d.",
    },
    {
        "id": "7ver",
        "name": "7-Ver",
        "family": "reinforcement-learning search results of [SS24]",
        "source": f"{THESIS}, Section 5.2.1, p.50",
        "provenance": "verbatim",
        "n_vertices": 7,
        "cnf": [
            [1, 2, 6], [1, 5, 7], [3, 4, 7], [2, 5, 7], [1, 5, 6], [1, 3, 7],
            [3, 4, 5], [2, 3, 4], [2, 4, 6], [1, 6, 7], [2, 5, 6],
        ],
        "dnf": [
            [3, 6, 7], [3, 5, 6], [1, 2, 3], [2, 5, 7],
            [4, 6, 7], [1, 4, 5], [1, 2, 4],
        ],
        "notes": (
            "epsilon = 5/11, suboptimal against the Fano plane's 3/7. The only "
            "instance from [SS24] that is not conformal (thesis p.55)."
        ),
    },
    {
        "id": "8ver-verbatim",
        "name": "8-Ver (as published)",
        "family": "reinforcement-learning search results of [SS24]",
        "source": f"{THESIS}, Section 5.2.1, p.50",
        "provenance": "verbatim",
        "n_vertices": 8,
        "cnf": [[2, 6, 7, 8], [3, 7, 8], [1, 4, 5], [2, 5, 6], [1, 3, 4]],
        "dnf": [
            [3, 4, 6], [2, 3, 5], [8, 1, 2], [8, 3, 5], [4, 5, 7], [2, 4, 7],
            [1, 5, 7], [1, 2, 7], [3, 5, 7], [1, 3, 7], [8, 1, 6], [8, 4, 6],
            [3, 5, 6], [8, 2, 4], [8, 1, 5], [2, 3, 4], [4, 6, 7], [1, 3, 6],
            [8, 4, 5],
        ],
        "notes": (
            "NOT a transversal pair as printed. D lists 19 edges; Tr(C) has 20. "
            "D omits {1,2,3} and {1,6,7}, and includes {1,3,7}, which is not a "
            "transversal of C at all -- it misses the edge {2,5,6}. Kept "
            "verbatim; see instance '8ver' for the repair."
        ),
    },
    {
        "id": "8ver",
        "name": "8-Ver (corrected)",
        "family": "reinforcement-learning search results of [SS24]",
        "source": f"{THESIS}, Section 5.2.1, p.50, with D recomputed as Tr(C)",
        "provenance": "corrected",
        "n_vertices": 8,
        "cnf": [[2, 6, 7, 8], [3, 7, 8], [1, 4, 5], [2, 5, 6], [1, 3, 4]],
        "dnf": [
            [1, 2, 3], [2, 3, 4], [2, 3, 5], [1, 3, 6], [3, 4, 6], [3, 5, 6],
            [1, 2, 7], [2, 4, 7], [1, 5, 7], [3, 5, 7], [4, 5, 7], [1, 6, 7],
            [4, 6, 7], [1, 2, 8], [2, 4, 8], [1, 5, 8], [3, 5, 8], [4, 5, 8],
            [1, 6, 8], [4, 6, 8],
        ],
        "notes": (
            "D replaced by Tr(C), taking C as published. Evidence the repair is "
            "the intended instance: it is mutually transversal in both "
            "directions, it agrees with 18 of the 19 published edges, and it "
            "reproduces the thesis' reported epsilon(C,D) = 2/5 (p.55), which "
            "the published D does not."
        ),
    },
    {
        "id": "trivial-aut-1",
        "name": "Trivial-Aut 1",
        "family": "asymmetric hypergraph on 8 vertices",
        "source": "_archive/20260804-legacy-code/'#Trivial Aut 1 fk-A Decomposition.txt' run log, 2024",
        "provenance": "derived",
        "n_vertices": 8,
        "cnf": None,  # filled in below from the archived log
        "dnf": None,
        "notes": (
            "Recovered from the archived FK-A decomposition log. Both sides have "
            "trivial automorphism group, matching the thesis' description of the "
            "Asym-1 / Asym-2 constructions (Section 5.2.1, p.50), whose edge "
            "lists were never printed. This is not necessarily the thesis' "
            "Asym-1: that figure splits on v2 at the root, this splits on v8. "
            "The archived log agrees with this implementation on all 41 nodes "
            "and all 41 split vertices under the modified variant, so it is the "
            "instance the log was produced from."
        ),
    },
]


def _self_dual_fano_spec(k: int) -> dict:
    """The ``SDFP(k)`` benchmark as a ``(E, E)`` instance.

    Self-transversal, so both sides are the same hypergraph and the pair has no
    free parameters. Generated from :func:`fka.instances.self_dual_fano` rather
    than copied from the FK-B snapshot's ``.mat`` files: the construction is
    stated in full there, the generated form reproduces the published shape
    exactly, and generating it keeps the repository free of third-party data
    whose redistribution terms are not established.
    """
    from fka.instances import self_dual_fano

    edges = self_dual_fano(k).to_sets(one_indexed=True)
    return {
        "id": f"sdfp-sd-{k}",
        "name": f"SDFP({k}) — self-dualised Fano plane",
        "family": "self-dual dualisation benchmark",
        "source": (
            "Self-dualisation of thesis Definition 5.2.1's k disjoint Fano "
            "planes. Matches the shape of SDFP16_CNF_DNF.mat (k=2) and "
            "SDFP23_CNF_DNF.mat (k=3) in the FK-B MATLAB reference"
        ),
        "provenance": "derived",
        "n_vertices": 7 * k + 2,
        "cnf": edges,
        "dnf": edges,
        "notes": (
            f"On 7k+2 = {7 * k + 2} vertices with k = {k}. Writing a and b for "
            "the two extra vertices, the edges are {a,b}; every Fano line plus "
            "b; and every choice of one line per copy plus a. That is "
            f"7^k + 7k + 1 = {7**k + 7 * k + 1} edges. Self-transversal: "
            "Tr(E) = E, verified against the Berge oracle, so (E, E) is a "
            "transversal pair by construction and the answer is always TRUE. "
            "This is the shape that separates the two algorithms most sharply "
            "-- see docs/fk-a-vs-fk-b.md. Vertices use this project's "
            "FANO_LINES labelling, which is isomorphic to but not identical "
            "with the MATLAB snapshot's, so node counts may differ from a "
            "MATLAB run by the amount vertex numbering moves the split choice."
        ),
    }


DEFINITIONS.extend(_self_dual_fano_spec(k) for k in (1, 2))


def _witt_11_specs() -> list[dict]:
    """``W11 = S(4,5,11)``, and its self-dualisation.

    Found by ``experiments/hard_case_search.py``: over every transitive group
    of degree <= 15 in SageMath's database, and every orbit of subsets under
    it, ``W11`` is the vertex-transitive self-transversal hypergraph whose
    recursion trees are largest per unit of input -- ahead of the Fano plane on
    both algorithms. See ``docs/hard-cases.md``.
    """
    from fka.instances import self_dualise, witt_11

    W = witt_11()
    blocks = W.to_sets(one_indexed=True)
    sd = self_dualise(blocks, 11, 1).to_sets(one_indexed=True)
    return [
        {
            "id": "w11",
            "name": "W11 — Witt design S(4,5,11)",
            "family": "self-transversal Steiner system",
            "source": (
                "The design of the perfect [11,6,5] ternary Golay code: the "
                "supports of its 132 minimum-weight words. Generated by "
                "fka.instances.witt_11"
            ),
            "provenance": "derived",
            "n_vertices": 11,
            "cnf": blocks,
            "dnf": blocks,
            "notes": (
                "Self-transversal: Tr(W11) = W11, verified against the Berge "
                "oracle. 66 blocks of size 5 on 11 points, every 4-subset in "
                "exactly one block, every vertex of degree 30, epsilon = 5/11. "
                "Aut = M11 of order 7920, vertex-transitive, and like the Fano "
                "plane it is neither conformal nor read-once. The Fano plane is "
                "the same construction on the perfect [7,4,3] Hamming code; the "
                "third perfect code, the binary Golay [23,12,7], gives S(4,7,23), "
                "which is intersecting but NOT self-transversal -- 16744 subsets "
                "S of its 23 points satisfy thesis equation 2.1 -- so the family "
                "stops at 11. This is the library's hardest instance per unit of "
                "input for FK-B and the only one that is hard for both "
                "algorithms at once; see docs/hard-cases.md."
            ),
        },
        {
            "id": "w11-sd-1",
            "name": "SD-W11(1) — self-dualised Witt design",
            "family": "self-dual dualisation benchmark",
            "source": "Self-dualisation of W11, the construction of SDFP(k) with the Fano plane replaced by S(4,5,11)",
            "provenance": "derived",
            "n_vertices": 13,
            "cnf": sd,
            "dnf": sd,
            "notes": (
                "SDFP(k)'s construction on a different base: {a,b}, every block "
                "plus b, and every block plus a, for 2*66 + 1 = 133 edges on 13 "
                "vertices. Self-transversal, verified against the Berge oracle. "
                "The point of comparison is sdfp-sd-1, the same construction on "
                "the Fano plane: this one keeps FK-A's node count per edge and "
                "roughly doubles FK-B's. k = 2 is not in the library because it "
                "would have 66^2 + 133 = 4489 edges."
            ),
        },
    ]


def _word_spec(w: str) -> dict:
    """``T_w`` for a word in the atoms X = (v1^v2)v(v3^v4) and Y = X^d."""
    from fka.hypergraph import Hypergraph
    from fka.substitution import substitute

    atoms = {"X": Hypergraph.from_sets(4, [[0, 1], [2, 3]]),
             "Y": Hypergraph.from_sets(4, [[0, 2], [0, 3], [1, 2], [1, 3]])}

    def build(word: str) -> Hypergraph:
        D = atoms[word[-1]]
        for ch in reversed(word[:-1]):
            D = substitute(atoms[ch], [D] * 4)
        return D

    swapped = w.translate(str.maketrans("XY", "YX"))
    C, D = build(w), build(swapped)
    return {
        "id": f"word-{w.lower()}",
        "name": f"Alternating tree of word {w}",
        "family": "read-once, word-indexed alternating tree",
        "source": "Derived here; generalises thesis 4.2's GK- and TK-families",
        "provenance": "derived",
        "n_vertices": C.n,
        "cnf": C.to_sets(one_indexed=True),
        "dnf": D.to_sets(one_indexed=True),
        "notes": (
            f"T_{w} = substitution of the atoms along the word {w}, with "
            f"Tr(T_{w}) = T_{swapped} by (S[T])^d = S^d[T^d] and X^d = Y. "
            "Writing a(w) = log2|E(T_w)|, composition gives "
            "a(w) + a(wbar) = 3(2^k - 1) for every word of length k, so the "
            "family's threshold ratio max(a, abar)/(2^k - 1) is at least 3/2, "
            "and a(w) - a(wbar) is odd so it is never met. The constant words "
            "XX..X and YY..Y are exactly the GK- and TK-families and are the "
            "worst-balanced ones, both at 2; this mixed word sits at 5/3, "
            "strictly past the tight families of Observation 2.3.3."
        ),
    }


DEFINITIONS.extend(_witt_11_specs())
DEFINITIONS.extend(_word_spec(w) for w in ("XY", "YX"))


def _load_trivial_aut() -> tuple[list[list[int]], list[list[int]]]:
    """Parse the two incidence matrices out of the archived decomposition log."""
    import ast
    import re

    path = (
        Path(__file__).resolve().parents[1]
        / "_archive"
        / "20260804-legacy-code"
        / "#Trivial Aut 1 fk-A Decomposition.txt"
    )
    if not path.exists():
        raise FileNotFoundError(f"the Trivial-Aut decomposition log is not at {path}")
    text = path.read_text(encoding="utf-8", errors="replace")
    from fka.hypergraph import Hypergraph

    def grab(name: str) -> list[list[int]]:
        m = re.search(rf"{name}=\s*(\[\[.*?\]\])", text, re.S)
        if not m:
            raise ValueError(f"{name} not found in {path}")
        D = Hypergraph.from_incidence(ast.literal_eval(m.group(1)))
        return D.to_sets(one_indexed=True)

    return grab("np_f"), grab("np_g")


#: Instances withdrawn from the live library on 2026-08-04 because they are not
#: transversal pairs as printed. They are still authored here so the archived
#: files stay regenerable, but they are written to the archive, not to
#: ``data/instances/``. See that directory's README, and the corrected ``6a``
#: and ``8ver`` above.
WITHDRAWN = {"6a-verbatim", "8ver-verbatim"}

ARCHIVE = "20260804-verbatim-nontransversal"


def main() -> int:
    target = instance_dir()
    target.mkdir(parents=True, exist_ok=True)
    archive = archive_root() / ARCHIVE
    archive.mkdir(parents=True, exist_ok=True)

    live = 0
    for spec in DEFINITIONS:
        spec = dict(spec)
        if spec["id"] == "trivial-aut-1":
            spec["cnf"], spec["dnf"] = _load_trivial_aut()
        inst = Instance(
            id=spec["id"],
            name=spec["name"],
            n_vertices=spec["n_vertices"],
            cnf_edges=spec["cnf"],
            dnf_edges=spec["dnf"],
            family=spec.get("family", ""),
            source=spec.get("source", ""),
            provenance=spec.get("provenance", "verbatim"),
            notes=spec.get("notes", ""),
        )
        refresh_expected(inst)
        withdrawn = inst.id in WITHDRAWN
        path = inst.save((archive if withdrawn else target) / f"{inst.id}.json")
        live += not withdrawn
        e = inst.expected
        print(
            f"{inst.id:16} dual={str(e['dual']):5} eps={e['epsilon']:>5} "
            f"nC={e['nC']:3} nD={e['nD']:3} "
            f"fk-a={e['fka']['modified']['nodes']:5} fk-b={e['fkb']['faithful']['nodes']:5} "
            f"-> {'_archive/' + ARCHIVE if withdrawn else 'data/instances'}/{path.name}"
        )

    print(f"\n{live} live instances in {target}")
    print(f"{len(WITHDRAWN)} withdrawn instances in {archive}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
