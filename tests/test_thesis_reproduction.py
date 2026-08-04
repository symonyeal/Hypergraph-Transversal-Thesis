"""Reproduction of the published recursion trees.

These are the tests that make the rewrite trustworthy. The thesis prints, for
several instances, an automorphism tree whose nodes are numbered 1..N and
labelled with a split vertex. Both the node count and the full split-vertex
sequence are checked here against the figures, and against the archived run log
for the instance the thesis did not print.

Any change to Sperner reduction, pivot selection, base cases, or recursion
order will move these numbers. That is the point: they are the contract with
the published results.

All of these use ``variant="modified"`` and ``pivot_rule="degree_sum"``, which
is what the original implementation did and therefore what produced the figures.
"""

from __future__ import annotations

import pytest

from fka.algorithm import fk_a
from fka.instances import load
from fkb.algorithm import fk_b

#: instance id -> (node count, split vertices in node order, thesis reference)
#: "-" marks a node that terminated without splitting.
PUBLISHED_TREES = {
    "f2g2": (
        25,
        "v1 v2 v5 v6 - - - v3 - v4 v5 v6 - - - - v3 - v4 v5 v6 - - - -",
        "Figure 5.1 (p.52)",
    ),
    "fano": (
        35,
        "v1 v2 v3 - v4 v6 - - - v3 v5 - v6 - - v6 - - "
        "v2 v3 - v4 v6 - - - v3 v5 - v6 - - v6 - -",
        "Figure 5.2 (p.54)",
    ),
    "6a": (11, None, "Figure 5.4a (p.57)"),
    "6b": (13, None, "Figure 5.4b (p.57)"),
    "6c": (5, None, "Figure 5.4c (p.57)"),
    "6d": (9, None, "Figure 5.4d (p.57)"),
    "7ver": (21, None, "Figure 5 (p.63)"),
    "8ver": (21, None, "Figure 6 (p.64)"),
}


@pytest.mark.parametrize("instance_id", sorted(PUBLISHED_TREES))
def test_node_count_matches_thesis_figure(instance_id):
    nodes, _, reference = PUBLISHED_TREES[instance_id]
    inst = load(instance_id)
    tree = fk_a(inst.G, inst.H, variant="modified", pivot_rule="degree_sum")
    assert len(tree) == nodes, (
        f"{instance_id}: thesis {reference} numbers nodes 1..{nodes}, got {len(tree)}"
    )
    assert tree.dual is True


@pytest.mark.parametrize(
    "instance_id",
    [k for k, v in PUBLISHED_TREES.items() if v[1] is not None],
)
def test_split_vertices_match_thesis_figure(instance_id):
    _, splits, reference = PUBLISHED_TREES[instance_id]
    inst = load(instance_id)
    tree = fk_a(inst.G, inst.H, variant="modified", pivot_rule="degree_sum")
    got = " ".join(node.pivot_label() or "-" for node in tree)
    assert got == splits, f"{instance_id} ({reference})"


def test_trivial_aut_matches_the_archived_run_log():
    """The 2024 decomposition log recorded 41 nodes and an MFV per node.

    Reproducing it exactly ties this implementation to the original one on an
    instance with 31 and 35 hyperedges -- far past the point where agreement
    could be coincidence.
    """
    import re
    from pathlib import Path

    log_path = (
        Path(__file__).resolve().parents[1]
        / "_archive"
        / "20260804-legacy-code"
        / "#Trivial Aut 1 fk-A Decomposition.txt"
    )
    assert log_path.exists(), (
        f"the archived decomposition log is missing from {log_path.parent}; "
        "it is the only evidence tying this implementation to the original"
    )
    text = log_path.read_text(encoding="utf-8", errors="replace")
    logged: dict[int, str] = {}
    for block in re.split(r"Node Number: ", text)[1:]:
        node_id = int(block.split()[0])
        m = re.search(r"MFV= x_(\d+)", block)
        logged[node_id] = f"v{m.group(1)}" if m else "-"

    inst = load("trivial-aut-1")
    tree = fk_a(inst.G, inst.H, variant="modified", pivot_rule="degree_sum")
    got = {node.node_id: (node.pivot_label() or "-") for node in tree}

    assert len(tree) == len(logged) == 41
    assert got == logged


def test_archived_verbatim_instances_are_not_transversal_pairs():
    """Two published instances are not dual as printed; both repairs are recorded.

    The verbatim forms were archived on 2026-08-04 out of ``data/instances/``,
    because every run against them answers "not dual" -- a fact about the
    transcription, not about either algorithm. They are re-derived here from the
    archive so the finding stays executable: if a future edit "fixes" them in
    place, or drops them, this fails and says why.
    """
    from fka.instances import load_archived
    from fka.transversal import is_dual_oracle

    for verbatim_id, corrected_id in (
        ("6a-verbatim", "6a"),
        ("8ver-verbatim", "8ver"),
    ):
        verbatim = load_archived(verbatim_id)
        corrected = load(corrected_id)
        assert (verbatim.path.parent / "README.md").exists(), (
            "an archive must state what was withdrawn and why"
        )

        # Re-derived, not read from the stored baseline.
        assert is_dual_oracle(verbatim.G, verbatim.H) is False, verbatim_id
        assert is_dual_oracle(corrected.G, corrected.H) is True, corrected_id
        assert corrected.provenance == "corrected"
        assert corrected.notes.strip(), "a correction must state its evidence"

        # Both algorithms must agree with the oracle on the archived form too.
        assert fk_a(verbatim.G, verbatim.H).dual is False
        assert fk_b(verbatim.G, verbatim.H).dual is False


def test_archived_verbatim_instances_are_out_of_the_live_library():
    from fka.instances import list_ids

    assert "6a-verbatim" not in list_ids()
    assert "8ver-verbatim" not in list_ids()


def test_corrected_instances_reproduce_the_thesis_reported_epsilon():
    """The repairs are justified by matching the thesis' own reported values.

    6-A is quoted at epsilon = 1/2 and 8-Ver at 2/5 (thesis p.55). The verbatim
    8-Ver gives 8/19; only the repair gives 2/5. That is the evidence for the
    correction, so it is recomputed here rather than read from a baseline.
    """
    from fka.instances import exact_epsilon, load_archived

    assert load("6a").expected["epsilon"] == "1/2"
    assert load("8ver").expected["epsilon"] == "2/5"

    archived = load_archived("8ver-verbatim")
    assert str(exact_epsilon(archived.G, archived.H)) == "8/19"


def test_committed_baselines_still_hold():
    """Every instance file's ``expected`` block must match a fresh computation.

    Regenerate with ``python -m fka refresh`` and review the diff when this
    fails -- a moved node count is a real change in the algorithm.
    """
    from fka.instances import load_all, refresh_expected

    for inst in load_all():
        committed = dict(inst.expected)
        recomputed = refresh_expected(inst)
        assert committed == recomputed, (
            f"{inst.id}: baseline drift.\n"
            f"  committed  {committed}\n"
            f"  recomputed {recomputed}"
        )
