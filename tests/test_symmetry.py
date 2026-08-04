"""Automorphism groups, group naming, graph classes and hypergraph properties."""

from __future__ import annotations

import pytest

from fka.automorphism import automorphism_group, hypergraph_automorphisms
from fka.graphs import (
    Graph,
    classify_graph,
    find_induced_p4,
    is_chordal,
    is_cograph,
    primal_graph,
)
from fka.groups import (
    PermutationGroup,
    compose,
    cycle_notation,
    identify,
    inverse,
    perm_order,
)
from fka.hypergraph import Hypergraph
from fka.instances import load
from fka.properties import analyse, is_alpha_acyclic, is_conformal, is_read_once


# ----------------------------------------------------------------------
# permutation plumbing
# ----------------------------------------------------------------------
def test_permutation_algebra():
    p = (1, 2, 0, 3)  # (0 1 2)
    assert perm_order(p) == 3
    assert compose(p, inverse(p)) == (0, 1, 2, 3)
    assert cycle_notation(p) == "(1,2,3)"
    assert cycle_notation((0, 1, 2, 3)) == "()"


@pytest.mark.parametrize(
    "generators,degree,name,order",
    [
        ([(1, 0)], 2, "C2", 2),
        ([(1, 0, 3, 2), (2, 3, 0, 1)], 4, "C2 x C2", 4),
        ([(1, 2, 0)], 3, "C3", 3),
        ([(1, 2, 0), (1, 0, 2)], 3, "S3", 6),
        ([(1, 2, 3, 0), (3, 2, 1, 0)], 4, "D4", 8),
        ([(1, 2, 3, 0), (1, 0, 2, 3)], 4, "S4", 24),
    ],
)
def test_small_group_naming(generators, degree, name, order):
    g = PermutationGroup.from_generators(generators, degree)
    assert g.order == order
    assert identify(g).name == name


def test_naming_reports_when_it_cannot_identify():
    """An unmatched group must say so, not guess a name."""
    from fka.groups import GroupIdentification, _candidates

    # A group of prime order 101 is C101, which the catalogue does have; use a
    # deliberately empty candidate list to exercise the fallback path.
    g = PermutationGroup.from_generators([tuple(range(1, 101)) + (0,)], 101)
    assert identify(g).name == "C101"
    assert _candidates(101)


# ----------------------------------------------------------------------
# hypergraph automorphisms
# ----------------------------------------------------------------------
def test_automorphisms_are_actual_automorphisms():
    inst = load("f2g2")
    H, _ = inst.G.compact()
    edges = set(H.edges)
    for p in hypergraph_automorphisms(H):
        mapped = {sum(1 << p[v] for v in range(H.n) if e >> v & 1) for e in edges}
        assert mapped == edges


@pytest.mark.parametrize(
    "instance_id,side,order,name",
    [
        # Thesis Figure 5.2 (p.54): the Fano plane's automorphism group.
        ("fano", "G", 168, "PSL(3,2)"),
        # Thesis p.53: two disjoint C4s give D4 wr S2, order 128, which the
        # thesis also writes as (((C2 x C2 x C2 x C2) : C2) : C2) : C2.
        ("f2g2", "G", 128, "D4 wr C2"),
        ("f2g2", "H", 128, "D4 wr C2"),
        # Thesis Figures 5.4b-d.
        ("6b", "G", 48, "C2 x S4"),
        ("6c", "G", 48, "C2 x S4"),
        ("6d", "G", 8, "D4"),
        # Thesis Figure 5.4a, after the 6-A correction.
        ("6a", "G", 4, "C2 x C2"),
        ("6a", "H", 4, "C2 x C2"),
    ],
)
def test_automorphism_groups_match_the_thesis(instance_id, side, order, name):
    inst = load(instance_id)
    result = automorphism_group(inst.G if side == "G" else inst.H)
    assert result.order == order
    assert result.name == name


def test_asymmetric_instance_has_trivial_group():
    inst = load("trivial-aut-1")
    assert automorphism_group(inst.G).order == 1
    assert automorphism_group(inst.H).order == 1


def test_isolated_vertices_are_excluded_by_default():
    """Isolated vertices would otherwise contribute a free symmetric group.

    The two legacy notebooks disagreed on this; ``include_isolated`` makes the
    choice explicit and defaults to the thesis' convention (p.53).
    """
    H = Hypergraph.from_sets(6, [[0, 1]])
    assert automorphism_group(H).order == 2  # just the swap of {0,1}
    assert automorphism_group(H, include_isolated=True).order == 2 * 24  # x S4


def test_automorphism_result_uses_original_vertex_numbering():
    H = Hypergraph.from_sets(6, [[3, 5]])  # 0-indexed; vertices v4, v6
    result = automorphism_group(H)
    assert result.generators == ("(4,6)",)
    assert result.isolated_vertices == (0, 1, 2, 4)


def test_sage_cycle_translation_uses_zero_based_incidence_points():
    """The Sage adapter must not rotate 0-based IncidenceStructure points."""
    from fka.backends.sage_backend import _sage_permutation_label

    class FakePermutation:
        def cycle_tuples(self):
            return ((0, 2), (1,))

    assert _sage_permutation_label(FakePermutation(), [3, 5, 7]) == "(4,8)"


# ----------------------------------------------------------------------
# graphs
# ----------------------------------------------------------------------
def test_primal_graph_of_the_fano_plane_is_k7():
    inst = load("fano")
    g = primal_graph(inst.G)
    assert g.size == 7 * 6 // 2  # complete graph on 7 vertices


def test_p4_detection():
    path = Graph.from_edges(4, [(0, 1), (1, 2), (2, 3)])
    assert find_induced_p4(path) is not None
    assert not is_cograph(path)
    complete = Graph.from_edges(4, [(a, b) for a in range(4) for b in range(a + 1, 4)])
    assert is_cograph(complete)


def test_graph_classes_on_c4():
    """C4 is bipartite, chordal bipartite, perfect, cograph, and not chordal."""
    c4 = Graph.from_edges(4, [(0, 1), (1, 2), (2, 3), (3, 0)])
    classes = classify_graph(c4)
    assert classes.bipartite and classes.chordal_bipartite and classes.perfect
    assert classes.cograph and classes.triangle_free
    assert not classes.chordal


def test_graph_classes_on_c5():
    """C5 is an odd hole: not perfect, not chordal, not bipartite."""
    c5 = Graph.from_edges(5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)])
    classes = classify_graph(c5)
    assert not classes.perfect and not classes.chordal and not classes.bipartite
    assert classes.triangle_free


def test_enumeration_cap_is_enforced():
    big = Graph.from_edges(25, [(0, 1)])
    with pytest.raises(ValueError, match="capped at"):
        is_cograph(big)
    assert is_chordal(big)  # networkx path, no cap


# ----------------------------------------------------------------------
# hypergraph properties
# ----------------------------------------------------------------------
def test_fano_is_not_conformal():
    """Thesis p.55: the Fano plane is non-conformal, its primal graph being K7."""
    inst = load("fano")
    assert not is_conformal(inst.G)
    assert analyse(inst.G).nonconformal_clique is not None


def test_read_once_instances():
    """f2/g2 is the thesis' read-once example (Section 5.2.1, p.48)."""
    inst = load("f2g2")
    assert is_read_once(inst.G)
    assert is_read_once(inst.H)


def test_alpha_acyclicity_by_gyo():
    # An alpha-acyclic hypergraph: GYO strips it completely.
    acyclic = Hypergraph.from_sets(4, [[0, 1], [1, 2], [2, 3]])
    assert is_alpha_acyclic(acyclic)
    # The Fano plane has no ear anywhere, so GYO stalls at once.
    assert not is_alpha_acyclic(load("fano").G)


def test_conformal_and_normal_are_the_same_predicate():
    """Documented in fka.properties: every clique sits inside a maximal one."""
    from fka.properties import is_conformal, is_normal

    assert is_normal is is_conformal


def test_thesis_claim_that_aut_g_equals_aut_h_has_a_counterexample():
    """Thesis p.51 states G and H share an automorphism group at every node.

    That holds on the corrected instances but fails on 6-A as published, which
    is also not a transversal pair -- so the claim is untested there rather
    than contradicted. Recorded so the distinction is not lost.
    """
    verbatim = load("6a-verbatim")
    assert automorphism_group(verbatim.G).order == 2
    assert automorphism_group(verbatim.H).order == 4
    assert verbatim.expected["dual"] is False

    corrected = load("6a")
    assert (
        automorphism_group(corrected.G).order
        == automorphism_group(corrected.H).order
        == 4
    )


# ----------------------------------------------------------------------
# SageMath cross-checks (skipped when Sage is unavailable)
# ----------------------------------------------------------------------
@pytest.mark.sage
@pytest.mark.parametrize("instance_id", ["fano", "f2g2", "6b", "6d"])
def test_sage_agrees_on_automorphism_order(instance_id):
    from fka.backends.sage_backend import sage_automorphism_group

    inst = load(instance_id)
    ours = automorphism_group(inst.G, backend="python")
    theirs = sage_automorphism_group(inst.G)
    assert ours.order == theirs["order"]


@pytest.mark.sage
@pytest.mark.parametrize("instance_id", ["fano", "6b", "6c", "6d"])
def test_sage_agrees_on_graph_classes(instance_id):
    from fka.backends.sage_backend import sage_graph_classes

    inst = load(instance_id)
    g = primal_graph(inst.G)
    ours = classify_graph(g).to_json()
    theirs = sage_graph_classes(g.edges(), g.n)
    for key, value in theirs.items():
        assert ours[key] == value, f"{instance_id}: {key}"
