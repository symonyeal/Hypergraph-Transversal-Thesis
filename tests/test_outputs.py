"""Tree serialisation, analysis annotation, HTML/DOT output, and the CLI."""

from __future__ import annotations

import json
import re

import pytest

from fka.algorithm import fk_a
from fka.analysis import AnalysisOptions, annotate
from fka.dot import to_dot
from fka.instances import instance_dir, list_ids, load, sdfp
from fka.report import layout, render_html, write_report
from fka.tree import Tree


@pytest.fixture(scope="module")
def small_tree():
    inst = load("6d")
    return fk_a(inst.cnf, inst.dnf, instance="6d")


# ----------------------------------------------------------------------
# tree model
# ----------------------------------------------------------------------
def test_tree_structure_invariants(small_tree):
    tree = small_tree
    assert tree.root.parent_id is None
    assert tree.root.depth == 0
    for node in tree:
        for child_id in node.children:
            child = tree[child_id]
            assert child.parent_id == node.node_id
            assert child.depth == node.depth + 1
            assert child.step in ("L", "R")
        # An internal node split; a leaf did not.
        assert bool(node.children) == (node.split_var is not None)


def test_node_ids_follow_preorder(small_tree):
    """Node numbering must be the recursion order, matching the thesis figures."""
    assert [n.node_id for n in small_tree.preorder()] == sorted(small_tree.nodes)


def test_tree_json_round_trip(small_tree, tmp_path):
    path = small_tree.save(tmp_path / "t.json")
    restored = Tree.load(path)
    assert len(restored) == len(small_tree)
    assert restored.dual == small_tree.dual
    for a, b in zip(small_tree, restored):
        assert a.node_id == b.node_id
        assert a.cnf.edges == b.cnf.edges
        assert a.dnf.edges == b.dnf.edges
        assert a.split_var == b.split_var
        assert a.path == b.path


def test_summary_counts_verdicts(small_tree):
    s = small_tree.summary()
    assert s["nodes"] == len(small_tree)
    assert s["leaves"] == len(small_tree.leaves)
    assert sum(s["verdicts"].values()) == len(small_tree)


# ----------------------------------------------------------------------
# analysis
# ----------------------------------------------------------------------
def test_annotation_fills_every_node(small_tree):
    tree = fk_a(load("6d").cnf, load("6d").dnf, instance="6d")
    annotate(tree, AnalysisOptions(graph_classes=False))
    for node in tree:
        assert "props_C" in node.analysis
        assert "props_D" in node.analysis
        assert "aut_label" in node.analysis


def test_annotation_records_root_group_of_the_fano_plane():
    tree = fk_a(load("fano").cnf, load("fano").dnf, instance="fano")
    annotate(tree, AnalysisOptions(graph_classes=False))
    assert tree.root.analysis["aut_C"]["name"] == "PSL(3,2)"
    assert tree.root.analysis["aut_agree"] is True


def test_annotation_is_memoised_on_repeated_subproblems():
    """Thesis p.53 notes FK-A regenerates identical subproblems.

    Annotation caches on the (cnf, MBF) pair, so identical nodes must come out with
    identical analysis dicts.
    """
    tree = fk_a(load("f2g2").cnf, load("f2g2").dnf, instance="f2g2")
    annotate(tree, AnalysisOptions(graph_classes=False))
    seen: dict[tuple, str] = {}
    for node in tree:
        key = (node.cnf.edges, node.dnf.edges)
        blob = json.dumps(node.analysis, sort_keys=True)
        if key in seen:
            assert seen[key] == blob
        seen[key] = blob
    assert len(seen) < len(tree), "expected at least one repeated subproblem"


# ----------------------------------------------------------------------
# report
# ----------------------------------------------------------------------
def test_layout_places_every_node_and_respects_depth(small_tree):
    pos = layout(small_tree)
    assert set(pos) == set(small_tree.nodes)
    for node in small_tree:
        _, y = pos[node.node_id]
        assert y == pytest.approx(node.depth * 132)
    # Parents sit between their children horizontally.
    for node in small_tree:
        if len(node.children) == 2:
            px = pos[node.node_id][0]
            left, right = (pos[c][0] for c in node.children)
            assert min(left, right) <= px <= max(left, right)


def test_report_is_self_contained(small_tree, tmp_path):
    path = write_report(small_tree, tmp_path / "r.html", instance=load("6d"))
    page = path.read_text(encoding="utf-8")
    assert page.startswith("<!doctype html>")
    # No request may leave the page: a strict-CSP or offline viewer must work.
    assert not re.findall(r'(?:src|href)\s*=\s*"(?!#)', page)
    assert "http://" not in page and "https://" not in page
    assert "</script>" in page and page.count("<script>") == 1


def test_report_embeds_every_node(small_tree):
    page = render_html(small_tree)
    payload = json.loads(re.search(r"const NODES = (\[.*?\]);", page, re.S).group(1))
    assert len(payload) == len(small_tree)
    assert {p["id"] for p in payload} == set(small_tree.nodes)


def test_report_supports_both_themes(small_tree):
    page = render_html(small_tree)
    assert "prefers-color-scheme: dark" in page
    assert 'data-theme="dark"' in page and 'data-theme="light"' in page


def test_report_opens_readably_and_nodes_are_keyboard_accessible(small_tree):
    page = render_html(small_tree)
    assert "if (NODES.length) select(NODES[0].id)" in page
    assert "function readable()" in page
    assert 'role: "button", tabindex: "0"' in page


def test_report_escapes_instance_names(small_tree):
    page = render_html(small_tree, title='<script>alert("x")</script>')
    assert "<script>alert" not in page
    assert "&lt;script&gt;alert" in page


# ----------------------------------------------------------------------
# dot
# ----------------------------------------------------------------------
def test_dot_export(small_tree):
    src = to_dot(small_tree)
    assert src.startswith("digraph FKA {")
    assert src.rstrip().endswith("}")
    assert src.count(" -> ") == len(small_tree) - 1
    assert "style=dashed" in src  # right branches


def test_dot_rejects_unknown_view(small_tree):
    with pytest.raises(ValueError, match="unknown view"):
        to_dot(small_tree, view="nope")


# ----------------------------------------------------------------------
# instances
# ----------------------------------------------------------------------
def test_every_instance_file_parses_and_is_well_formed():
    for instance_id in list_ids():
        inst = load(instance_id)
        assert inst.id == instance_id
        assert inst.provenance in ("verbatim", "corrected", "derived")
        assert inst.source, f"{instance_id} has no source"
        assert inst.cnf.n == inst.dnf.n == inst.n_vertices
        assert inst.expected, f"{instance_id} has no baseline"


def test_missing_instance_names_the_alternatives():
    with pytest.raises(FileNotFoundError, match="Available:"):
        load("does-not-exist")


def test_instance_dir_resolves_from_the_package():
    assert instance_dir().is_dir()
    assert (instance_dir() / "fano.json").exists()


def test_sdfp_generator():
    one = sdfp(1)
    assert one.n == 7 and len(one) == 7
    assert one.edges == load("fano").cnf.edges
    two = sdfp(2)
    assert two.n == 14 and len(two) == 14
    # The copies are vertex-disjoint.
    first, second = two.edges[:7], two.edges[7:]
    assert not (sum(first) & sum(second) or (max(first) & max(second)))
    with pytest.raises(ValueError, match="at least 1"):
        sdfp(0)


# ----------------------------------------------------------------------
# cli
# ----------------------------------------------------------------------
def test_cli_env_and_list(capsys):
    from fka.cli import main

    assert main(["env"]) == 0
    assert "active_backend" in capsys.readouterr().out
    assert main(["list"]) == 0
    assert "fano" in capsys.readouterr().out


def test_cli_show(capsys):
    from fka.cli import main

    assert main(["show", "fano"]) == 0
    out = capsys.readouterr().out
    assert "Fano" in out and "PSL" not in out  # show reports data, not analysis

    assert main(["show", "fano", "--notes"]) == 0
    assert "PSL" in capsys.readouterr().out


def test_cli_run_writes_outputs(tmp_path, capsys):
    """Artifacts are named ``<instance>.<algorithm>.<variant>.<ext>``.

    The algorithm is in the stem so FK-A and FK-B runs of the same instance and
    variant do not overwrite each other.
    """
    from fka.cli import main

    rc = main(["run", "k3", "--out", str(tmp_path), "--no-annotate", "--dot"])
    assert rc == 0
    for ext in ("html", "json", "dot"):
        assert (tmp_path / f"k3.fk-a.modified.{ext}").exists()


def test_cli_run_without_target_explains(capsys):
    from fka.cli import main

    assert main(["run"]) == 2
    assert "--all" in capsys.readouterr().out


def test_cli_requesting_sage_when_absent_is_an_error():
    from fka.backends import has_sage, resolve_backend

    if has_sage():
        pytest.skip("SageMath is available here")
    with pytest.raises(RuntimeError, match="not importable"):
        resolve_backend("sage")
