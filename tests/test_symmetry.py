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


def test_naming_reports_when_it_cannot_identify(monkeypatch):
    """An unmatched group must say so, not guess a name."""
    import fka.groups as groups

    g = PermutationGroup.from_generators([tuple(range(1, 101)) + (0,)], 101)
    assert identify(g).name == "C101"

    # Starve every route to a name: no analytic match, no candidates, no
    # special case. The result must be the placeholder, not a guess.
    monkeypatch.setattr(groups, "_named_cyclic", lambda _: None)
    monkeypatch.setattr(groups, "_named_dihedral", lambda _: None)
    monkeypatch.setattr(groups, "_candidates", lambda _order: ())
    monkeypatch.setattr(groups, "_SPECIAL", {})
    result = groups.identify(g)
    assert result.method == "unidentified"
    assert result.name == "order 101 (abelian, unidentified)"


def test_cyclic_and_dihedral_are_named_without_being_built():
    """The two infinite families are recognised from invariants alone.

    ``_candidates`` stops at :data:`fka.groups.MAX_ATOM_ORDER`, so a larger
    cyclic or dihedral group has no constructed candidate to match against and
    must be named by :func:`_named_cyclic` / :func:`_named_dihedral`.
    """
    from fka.groups import MAX_ATOM_ORDER, _candidates, _dihedral

    n = MAX_ATOM_ORDER + 2
    big_cyclic = PermutationGroup.from_generators([tuple(range(1, n)) + (0,)], n)
    assert identify(big_cyclic).name == f"C{n}"
    assert not any(name == f"C{n}" for name, _ in _candidates(n))

    m = MAX_ATOM_ORDER  # D_m has order 2m, past the atom bound
    big_dihedral = _dihedral(m)
    assert identify(big_dihedral).name == f"D{m}"
    assert not any(name == f"D{m}" for name, _ in _candidates(2 * m))


# ----------------------------------------------------------------------
# hypergraph automorphisms
# ----------------------------------------------------------------------
def test_automorphisms_are_actual_automorphisms():
    inst = load("f2g2")
    MBF, _ = inst.cnf.compact()
    edges = set(MBF.edges)
    for p in hypergraph_automorphisms(MBF):
        mapped = {sum(1 << p[v] for v in range(MBF.n) if e >> v & 1) for e in edges}
        assert mapped == edges


@pytest.mark.parametrize(
    "instance_id,side,order,name",
    [
        # Thesis Figure 5.2 (p.54): the Fano plane's automorphism group.
        ("fano", "cnf", 168, "PSL(3,2)"),
        # Thesis p.53: two disjoint C4s give D4 wr S2, order 128, which the
        # thesis also writes as (((C2 x C2 x C2 x C2) : C2) : C2) : C2.
        ("f2g2", "cnf", 128, "D4 wr C2"),
        ("f2g2", "dnf", 128, "D4 wr C2"),
        # Thesis Figures 5.4b-d.
        ("6b", "cnf", 48, "C2 x S4"),
        ("6c", "cnf", 48, "C2 x S4"),
        ("6d", "cnf", 8, "D4"),
        # Thesis Figure 5.4a, after the 6-A correction.
        ("6a", "cnf", 4, "C2 x C2"),
        ("6a", "dnf", 4, "C2 x C2"),
    ],
)
def test_automorphism_groups_match_the_thesis(instance_id, side, order, name):
    """The order is backend-independent; the *name* is the package catalogue's.

    ``backend="python"`` is explicit because GAP's ``StructureDescription`` is
    not canonical for these groups: handed the search's generators for
    ``f2g2`` it answers ``(((C4 x C4) : C2) : C2) : C2``, and handed an
    abstract ``D4 wr C2`` it answers ``(((C2 x C2 x C2 x C2) : C2) : C2) : C2``
    -- the thesis' own second spelling, and an isomorphic group either way.
    Under the default backend this assertion would therefore pass on CPython
    and fail under Sage. :func:`test_sage_group_names_are_isomorphic` checks
    what the Sage backend actually owes: a group of the right isomorphism type.
    """
    inst = load(instance_id)
    result = automorphism_group(getattr(inst, side), backend="python")
    assert result.order == order
    assert result.name == name


@pytest.mark.sage
@pytest.mark.parametrize(
    "instance_id,side,order",
    [("fano", "G", 168), ("f2g2", "G", 128), ("6b", "G", 48), ("w11", "G", 7920)],
)
def test_sage_group_names_are_isomorphic(instance_id, side, order):
    """Sage's name may spell the group differently; it must not describe a
    different group. Checked against the catalogue name through GAP."""
    from sage.all import libgap

    from fka.backends.sage_backend import sage_automorphism_group

    inst = load(instance_id)
    H = inst.G if side == "G" else inst.H
    py = automorphism_group(H, backend="python")
    sg = sage_automorphism_group(H)
    assert sg["order"] == py.order == order

    gens = [libgap.eval(g) for g in sg["generators"] if g != "()"]
    from_sage = libgap.Group(gens)
    from_catalogue = libgap.Group(
        [libgap.eval(g) for g in py.generators if g != "()"]
    )
    assert libgap.IsomorphismGroups(from_sage, from_catalogue) != libgap.eval("fail")


def test_asymmetric_instance_has_trivial_group():
    inst = load("trivial-aut-1")
    assert automorphism_group(inst.cnf).order == 1
    assert automorphism_group(inst.dnf).order == 1


def test_isolated_vertices_are_excluded_by_default():
    """Isolated vertices would otherwise contribute a free symmetric group.

    ``include_isolated`` makes the
    choice explicit and defaults to the thesis' convention (p.53).
    """
    MBF = Hypergraph.from_sets(6, [[0, 1]])
    assert automorphism_group(MBF).order == 2  # just the swap of {0,1}
    assert automorphism_group(MBF, include_isolated=True).order == 2 * 24  # x S4


def test_automorphism_result_uses_original_vertex_numbering():
    MBF = Hypergraph.from_sets(6, [[3, 5]])  # 0-indexed; vertices v4, v6
    result = automorphism_group(MBF)
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
    g = primal_graph(inst.cnf)
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
    assert not is_conformal(inst.cnf)
    assert analyse(inst.cnf).nonconformal_clique is not None


def test_read_once_instances():
    """f2/g2 is the thesis' read-once example (Section 5.2.1, p.48)."""
    inst = load("f2g2")
    assert is_read_once(inst.cnf)
    assert is_read_once(inst.dnf)


def test_alpha_acyclicity_by_gyo():
    # An alpha-acyclic hypergraph: GYO strips it completely.
    acyclic = Hypergraph.from_sets(4, [[0, 1], [1, 2], [2, 3]])
    assert is_alpha_acyclic(acyclic)
    # The Fano plane has no ear anywhere, so GYO stalls at once.
    assert not is_alpha_acyclic(load("fano").cnf)


def test_conformal_and_normal_are_the_same_predicate():
    """Documented in fka.properties: every clique sits inside a maximal one."""
    from fka.properties import is_conformal, is_conformal

    assert is_conformal is is_conformal


def test_thesis_claim_that_aut_g_equals_aut_h_needs_a_transversal_pair():
    """Thesis p.51 states cnf and MBF share an automorphism group at every node.

    The apparent counterexample is 6-A *as published*, where Aut(cnf) = C2 and
    Aut(MBF) = C2 x C2. That pair is not a transversal pair, which is exactly the
    hypothesis the claim carries, so it never contradicted it. Recorded so the
    distinction is not lost -- and so a future edit cannot quietly promote the
    archived form back into the library and call it a counterexample.
    """
    from fka.instances import load_archived
    from fka.transversal import is_dual_oracle

    verbatim = load_archived("6a-verbatim")
    assert automorphism_group(verbatim.cnf).order == 2
    assert automorphism_group(verbatim.dnf).order == 4
    assert is_dual_oracle(verbatim.cnf, verbatim.dnf) is False

    corrected = load("6a")
    assert is_dual_oracle(corrected.cnf, corrected.dnf) is True
    assert (
        automorphism_group(corrected.cnf).order
        == automorphism_group(corrected.dnf).order
        == 4
    )


@pytest.mark.slow
@pytest.mark.parametrize("algorithm", ["fk-a", "fk-b"])
def test_aut_g_equals_aut_h_at_every_node_of_both_algorithms(algorithm):
    """The thesis' p.51 claim, extended to every node and to FK-B.

    The thesis states it of the input pair. Measured over the whole recursion of
    both algorithms on every live instance, it holds at every node where both
    groups were computed -- 3455 of them as of 2026-08-04, with zero
    disagreements. See docs/fk-a-vs-fk-b.md section 7.

    FK-B is the part the thesis never examined: it decomposes by forcing whole
    clauses and monomials rather than single vertices, so there was no reason in
    advance for its subproblems to preserve the property.
    """
    from fka.algorithm import fk_a
    from fka.analysis import AnalysisOptions, annotate
    from fka.instances import load_all
    from fkb.algorithm import fk_b

    build = fk_a if algorithm == "fk-a" else fk_b
    opts = AnalysisOptions(graph_classes=False, skip_group_above=20)

    checked = 0
    for inst in load_all():
        tree = build(inst.cnf, inst.dnf, instance=inst.id)
        if len(tree) > 200:  # sdfp-sd-2 is covered by the committed report
            continue
        annotate(tree, opts)
        for node in tree:
            agree = node.analysis.get("aut_agree")
            if agree is None:
                continue  # one side's group was above the cap
            assert agree is True, (
                f"{algorithm} {inst.id} node {node.node_id}: "
                f"{node.analysis['aut_C']['name']} vs {node.analysis['aut_D']['name']} "
                f"for {node.cnf.label()} / {node.dnf.label()}"
            )
            checked += 1
    assert checked > 100, f"only {checked} nodes had both groups; not meaningful"


# ----------------------------------------------------------------------
# SageMath cross-checks (skipped when Sage is unavailable)
# ----------------------------------------------------------------------
@pytest.mark.sage
@pytest.mark.parametrize("instance_id", ["fano", "f2g2", "6b", "6d"])
def test_sage_agrees_on_automorphism_order(instance_id):
    from fka.backends.sage_backend import sage_automorphism_group

    inst = load(instance_id)
    ours = automorphism_group(inst.cnf, backend="python")
    theirs = sage_automorphism_group(inst.cnf)
    assert ours.order == theirs["order"]


@pytest.mark.sage
@pytest.mark.parametrize("instance_id", ["fano", "6b", "6c", "6d"])
def test_sage_agrees_on_graph_classes(instance_id):
    from fka.backends.sage_backend import sage_graph_classes

    inst = load(instance_id)
    g = primal_graph(inst.cnf)
    ours = classify_graph(g).to_json()
    theirs = sage_graph_classes(g.edges(), g.n)
    for key, value in theirs.items():
        assert ours[key] == value, f"{instance_id}: {key}"
