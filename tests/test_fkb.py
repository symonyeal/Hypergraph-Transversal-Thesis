"""FK-B: agreement with the oracle, with FK-A, and conformance to the MATLAB source.

The MATLAB reference cannot be run on this machine, so "faithful" is pinned the
only way it can be: against the same brute-force oracle FK-A is held to, against
FK-A itself, and against a direct reading of each helper's specification. Where
the port deliberately departs from ``FK- B/FK-master/src`` -- kept empty terms,
constants as hypergraphs, ``None`` rather than ``[]`` for "no conflict" -- there
is a test naming the case.
"""

from __future__ import annotations

import itertools

import pytest

from fka.algorithm import fk_a, verify_certificate
from fka.hypergraph import Hypergraph, bitset_to_vertices
from fka.instances import load, load_all
from fka.transversal import is_dual_oracle, transversal
from fkb.algorithm import (
    BRANCHES,
    SPLIT_RULES,
    VARIANTS,
    c_set_0,
    c_set_1,
    check_conditions,
    choose_split_var,
    cnf_value,
    d_set_0,
    d_set_1,
    dnf_value,
    easy_case,
    fk_b,
    is_conflict,
    is_dual,
    mu,
    mu_frequent,
    phi_x_0,
    phi_x_1,
)


def hg(n, edges):
    return Hypergraph.from_sets(n, edges, one_indexed=True)


def brute_conflict(cnf, dnf):
    """Every conflicting assignment, by enumeration. The specification, slowly."""
    n = max(cnf.n, dnf.n)
    return {S for S in range(1 << n) if cnf_value(cnf, S) != dnf_value(dnf, S)}


# ----------------------------------------------------------------------
# correctness against the independent oracle
# ----------------------------------------------------------------------
@pytest.mark.parametrize("instance_id", [i.id for i in load_all()])
@pytest.mark.parametrize("variant", VARIANTS)
@pytest.mark.parametrize("split_rule", SPLIT_RULES)
def test_agrees_with_oracle_on_library(instance_id, variant, split_rule):
    inst = load(instance_id)
    tree = fk_b(
        inst.G, inst.H, variant=variant, split_rule=split_rule, validate=True
    )
    assert tree.dual == is_dual_oracle(inst.G, inst.H)


@pytest.mark.parametrize("instance_id", [i.id for i in load_all()])
def test_agrees_with_fk_a_on_library(instance_id):
    """The two algorithms decide the same predicate, so they must never differ."""
    inst = load(instance_id)
    assert fk_b(inst.G, inst.H).dual == fk_a(inst.G, inst.H).dual


@pytest.mark.slow
def test_agrees_with_oracle_on_random_instances(rng):
    from conftest import random_sperner

    seeded_dual = 0
    for _ in range(300):
        n = rng.randint(2, 7)
        G = random_sperner(n, rng.randint(1, 5), rng)
        if not G.edges:
            continue
        if rng.random() < 0.5:
            H = transversal(G)
            seeded_dual += 1
        else:
            H = random_sperner(n, rng.randint(1, 5), rng)
        truth = is_dual_oracle(G, H)
        for variant in VARIANTS:
            for rule in SPLIT_RULES:
                tree = fk_b(G, H, variant=variant, split_rule=rule, validate=True)
                assert tree.dual == truth, (
                    f"n={n} C={G.label()} D={H.label()} "
                    f"variant={variant} rule={rule}"
                )
    assert seeded_dual > 50, "not enough genuinely equivalent cases were exercised"


@pytest.mark.slow
def test_root_assignment_conflicts_on_the_root_instance(rng):
    """The lifts must reach the root, not just be locally consistent.

    Each ``c``/``m`` branch hands its child an instance with whole terms forced,
    and the child's answer is only a counterexample after those variables are
    put back. This checks the composed lift end to end.
    """
    from conftest import random_sperner

    checked = 0
    for _ in range(300):
        n = rng.randint(2, 7)
        G = random_sperner(n, rng.randint(1, 5), rng)
        H = random_sperner(n, rng.randint(1, 5), rng)
        if not G.edges or not H.edges:
            continue
        tree = fk_b(G, H, variant="multiple")
        if tree.dual or tree.verdict is None:
            continue
        for S in tree.verdict.certificates or (tree.verdict.certificate,):
            assert is_conflict(G, H, S), f"C={G.label()} D={H.label()} S={S:b}"
            checked += 1
    assert checked > 40, "not enough assignments were produced to be meaningful"


def test_is_dual_convenience_wrapper():
    inst = load("fano")
    assert is_dual(inst.G, inst.H) is True


# ----------------------------------------------------------------------
# the two variants
# ----------------------------------------------------------------------
@pytest.mark.parametrize("instance_id", [i.id for i in load_all()])
def test_multiple_variant_leaves_the_tree_alone(instance_id):
    """``MFK_B.m`` differs from ``FK_B.m`` only in how much it reports.

    Its ``Multiple_*`` helpers return non-empty in exactly the same situations,
    so control flow -- and therefore the recursion tree -- is untouched.
    """
    inst = load(instance_id)
    one = fk_b(inst.G, inst.H, variant="faithful")
    many = fk_b(inst.G, inst.H, variant="multiple")
    assert len(one) == len(many)
    assert one.dual == many.dual
    assert [n.pivot for n in one] == [n.pivot for n in many]
    assert [n.path for n in one] == [n.path for n in many]


@pytest.mark.parametrize(
    "cnf,dnf,reason",
    [
        # Two monomials each miss the only clause.
        (hg(4, [[1, 2]]), hg(4, [[3], [4]]), "cond_1_disjoint"),
        # Two monomials each hold the variable C never mentions.
        (hg(3, [[1, 2]]), hg(3, [[1, 3], [2, 3]]), "cond_2_support"),
    ],
)
def test_multiple_variant_reports_every_assignment_it_finds(cnf, dnf, reason):
    """``MFK_B.m``'s point: one pass can expose several conflicts, not just one."""
    one = fk_b(cnf, dnf, variant="faithful", validate=True)
    many = fk_b(cnf, dnf, variant="multiple", validate=True)
    assert one.dual is many.dual is False
    assert one.verdict.reason == many.verdict.reason == reason

    found = many.verdict.certificates
    assert len(found) > 1
    assert found[0] == many.verdict.certificate == one.verdict.certificate
    assert one.verdict.certificates == (), "faithful must not collect the rest"
    for S in found:
        assert is_conflict(cnf, dnf, S)


# ----------------------------------------------------------------------
# split selection
# ----------------------------------------------------------------------
def test_mu_matches_its_definition():
    import math

    assert mu(100) == pytest.approx(math.log(100) / math.log(math.log(100)))
    # log log n is undefined at or below e; FK-B never reaches it, and the
    # guard must not make a variable look mu-infrequent.
    assert mu(2) == math.inf
    assert mu_frequent(0, hg(2, [[1], [2]]), hg(2, [[1, 2]])) is False


def test_mu_frequent_reads_the_named_side():
    """``mu_frequent(x, A, B)`` measures ``x`` in ``A``, with the bound from both."""
    C = hg(4, [[1, 2], [1, 3], [1, 4], [2, 3]])  # v1 in three of four clauses
    D = hg(4, [[1], [2], [3], [4]])  # every variable in one of four monomials
    threshold = 1 / mu(len(C) * len(D))
    assert mu_frequent(0, C, D) is (C.degree(0) / len(C) <= threshold)
    assert mu_frequent(3, C, D) is (C.degree(3) / len(C) <= threshold)
    # A variable in every term of A is never mu-infrequent there; one in a
    # single term of a large A is.
    assert mu_frequent(0, C, D) is False
    assert mu_frequent(0, D, C) is True


def test_most_frequent_rule_matches_choose_splitvar():
    """'mostFreq' walks both degree orders down together, lowest index winning.

    Fano is vertex-transitive: every variable ties at degree 3 on both sides, so
    the first pass already intersects and v1 must win.
    """
    inst = load("fano")
    assert choose_split_var(inst.G, inst.H, "most_frequent") == 0
    assert choose_split_var(inst.G, inst.H, "degree_sum") == 0


def test_split_var_ignores_variables_no_term_uses():
    """Contraction leaves variables in no term; splitting on one never terminates."""
    C = hg(5, [[1, 2], [2, 3]])
    D = hg(5, [[2], [1, 3]])
    assert choose_split_var(C, D) in (0, 1, 2)
    assert choose_split_var(Hypergraph.empty(4), Hypergraph.empty(4)) is None


def test_unknown_split_rule_is_rejected():
    with pytest.raises(ValueError, match="unknown split rule"):
        choose_split_var(hg(2, [[1]]), hg(2, [[1]]), "nope")


# ----------------------------------------------------------------------
# decomposition and restriction
# ----------------------------------------------------------------------
def test_phi_splits_a_family_on_a_variable():
    f = hg(4, [[1, 2], [1, 3], [2, 4]])
    assert phi_x_1(f, 0).to_sets(one_indexed=True) == [[2, 4]]
    assert phi_x_0(f, 0).to_sets(one_indexed=True) == [[2], [3]]


def test_phi_x_0_keeps_a_term_its_variable_emptied():
    """Deviation from ``phi_x_0.m``, which deletes it.

    A clause ``{x}`` reduced by ``x = 0`` is the empty clause, and that is
    precisely the fact that the CNF became the constant 0. Dropping it loses the
    only evidence of it, and singleton clauses are reachable: every monomial
    must then contain ``x``, which makes ``x`` the most frequent variable and so
    the one split on.
    """
    assert phi_x_0(hg(3, [[1], [2, 3]]), 0).edges == (0,)


def test_restrictions_agree_with_evaluation():
    """``c_set_*``/``d_set_*`` must be the restrictions their names claim."""
    C = hg(4, [[1, 2], [2, 3], [3, 4]])
    D = hg(4, [[2, 3], [1, 3], [2, 4]])
    for x in range(4):
        bit = 1 << x
        rest = [S for S in range(1 << 4) if not S & bit]
        for S in rest:
            assert cnf_value(c_set_0(C, bit), S) == cnf_value(C, S)
            assert dnf_value(d_set_0(D, bit), S) == dnf_value(D, S)
            assert cnf_value(c_set_1(C, bit), S) == cnf_value(C, S | bit)
            assert dnf_value(d_set_1(D, bit), S) == dnf_value(D, S | bit)


def test_restriction_pieces_reassemble_the_x_zero_form():
    C = hg(4, [[1, 2], [2, 3], [1, 4]])
    assert phi_x_0(C, 0).union(phi_x_1(C, 0)) == c_set_0(C, 1)


# ----------------------------------------------------------------------
# conditions
# ----------------------------------------------------------------------
def test_condition_1_disjoint_clause_and_monomial():
    C = hg(3, [[1]])
    D = hg(3, [[2]])
    v = check_conditions(C, D)
    assert v.reason == "cond_1_disjoint"
    assert is_conflict(C, D, v.certificate)
    # This is the one condition whose conflict runs the other way: the monomial
    # is satisfied while the clause it misses is not.
    assert dnf_value(D, v.certificate) and not cnf_value(C, v.certificate)


def test_condition_2_extra_variable_on_either_side():
    for C, D in (
        (hg(3, [[1, 2]]), hg(3, [[1, 3], [2]])),  # v3 only in D
        (hg(3, [[1, 2], [1, 3]]), hg(3, [[1]])),  # v2, v3 only in C
    ):
        v = check_conditions(C, D)
        assert v.reason == "cond_2_support"
        assert is_conflict(C, D, v.certificate)
        assert verify_certificate(C, D, v.certificate)


def test_condition_3_term_longer_than_the_other_side():
    for C, D in (
        (hg(3, [[1, 2], [1, 3]]), hg(3, [[1, 2, 3]])),  # monomial longer than |C|
        (hg(3, [[1, 2, 3]]), hg(3, [[1, 2], [1, 3]])),  # clause longer than |D|
    ):
        v = check_conditions(C, D)
        assert v.reason == "cond_3_size"
        assert is_conflict(C, D, v.certificate)
        # Everything but condition 1 yields equation 2.1's certificate, so FK-A
        # would accept the same witness.
        assert verify_certificate(C, D, v.certificate)


@pytest.mark.parametrize(
    "cnf,dnf,dual",
    [
        (Hypergraph.empty(3), hg(3, [[1], [2]]), False),  # C = 1, D is not
        (hg(3, [[1], [2]]), Hypergraph.empty(3), False),  # D = 0, C is not
        (Hypergraph(3, (0,)), hg(3, [[1], [2]]), False),  # C = 0, D is not
        (hg(3, [[1], [2]]), Hypergraph(3, (0,)), False),  # D = 1, C is not
        (Hypergraph.empty(3), Hypergraph(3, (0,)), True),  # both the constant 1
        (Hypergraph(3, (0,)), Hypergraph.empty(3), True),  # both the constant 0
    ],
)
def test_constant_sides_are_decided_not_recursed(cnf, dnf, dual):
    """Deviation from ``A_c_x.m``/``A_m_x.m``, which return a MATLAB scalar here.

    A collapsed side is an ordinary hypergraph -- ``{}`` or ``{{}}`` -- and
    reading which constant it is settles the instance immediately.
    """
    v = check_conditions(cnf, dnf)
    assert v is not None and v.reason == "constant"
    assert v.dual is dual
    assert fk_b(cnf, dnf, validate=True).dual is dual
    if not dual:
        assert is_conflict(cnf, dnf, v.certificate)


def test_empty_set_is_a_conflict_not_an_absent_one():
    """Deviation from the MATLAB, which uses ``[]`` for both readings.

    ``C = 1``, ``D = 0`` conflicts at the empty assignment, and reporting that
    as "no conflict" would call a non-equivalent pair equivalent.
    """
    v = check_conditions(Hypergraph.empty(3), Hypergraph.empty(3))
    assert v.dual is False
    assert v.certificate == 0


# ----------------------------------------------------------------------
# easy case
# ----------------------------------------------------------------------
@pytest.mark.parametrize("n", [3, 4])
def test_easy_case_matches_exhaustive_search(n):
    """The rewritten easy case must decide exactly what CA_CNF2.m searched for.

    The MATLAB enumerates subsets with ``dec2bin`` loops; this walks the small
    side's minimal transversals instead. Both are complete, so on every pair
    that reaches the easy case they must return a conflict in exactly the same
    situations -- checked here against all ``2^n`` assignments.
    """
    universe = [tuple(e) for k in range(1, n + 1) for e in itertools.combinations(range(1, n + 1), k)]
    checked = 0
    for c_edges in itertools.combinations(universe, 2):
        C = hg(n, [list(e) for e in c_edges])
        if not C.is_sperner():
            continue
        for size in (1, 2, 3):
            for d_edges in itertools.combinations(universe, size):
                D = hg(n, [list(e) for e in d_edges])
                if not D.is_sperner():
                    continue
                if check_conditions(C, D) is not None:
                    continue  # the easy case is only reached once these pass
                v = easy_case(C, D, variant="multiple")
                truth = brute_conflict(C, D)
                assert v.dual == (not truth), f"C={C.label()} D={D.label()}"
                for S in v.certificates:
                    assert S in truth
                checked += 1
    assert checked > 20, "the easy case was never reached"


# ----------------------------------------------------------------------
# the recursion tree
# ----------------------------------------------------------------------
def test_tree_is_labelled_fk_b():
    inst = load("fano")
    tree = fk_b(inst.G, inst.H, instance="fano")
    assert tree.algorithm == "FK-B"
    assert tree.summary()["algorithm"] == "FK-B"
    assert fk_a(inst.G, inst.H).algorithm == "FK-A"


def test_tree_structure_invariants():
    inst = load("8ver")
    tree = fk_b(inst.G, inst.H, instance="8ver")
    assert tree.root.parent_id is None
    for node in tree:
        for child_id in node.children:
            child = tree[child_id]
            assert child.parent_id == node.node_id
            assert child.depth == node.depth + 1
            assert child.branch in ("x=0", "x=1") or child.branch[0] in ("c", "m")
        assert bool(node.children) == (node.pivot is not None)
        assert bool(node.split_branch) == (node.pivot is not None)
    assert [n.node_id for n in tree.preorder()] == sorted(tree.nodes)


def test_every_split_branch_is_exercised_by_the_library():
    """All three FK-B branches must be live, or the mu test is dead code."""
    seen = set()
    for inst in load_all():
        for node in fk_b(inst.G, inst.H):
            if node.split_branch:
                seen.add(node.split_branch)
    assert seen == set(BRANCHES), f"unexercised branches: {set(BRANCHES) - seen}"


def test_per_term_branches_give_a_node_more_than_two_children():
    """The shape that distinguishes an FK-B tree from an FK-A one."""
    inst = load("6b")
    tree = fk_b(inst.G, inst.H)
    wide = [n for n in tree if len(n.children) > 2]
    assert wide, "expected at least one node with a subproblem per clause"
    assert wide[0].split_branch in ("mu_C", "mu_D")


def test_tree_json_round_trip_keeps_the_fk_b_fields(tmp_path):
    from fka.tree import RecursionTree

    inst = load("8ver")
    tree = fk_b(inst.G, inst.H, instance="8ver", variant="multiple")
    restored = RecursionTree.load(tree.save(tmp_path / "b.json"))
    assert restored.algorithm == "FK-B"
    for a, b in zip(tree, restored):
        assert a.path == b.path
        assert a.split_branch == b.split_branch
        if a.verdict is not None:
            assert a.verdict.certificates == b.verdict.certificates


def test_arguments_are_validated():
    inst = load("k3")
    with pytest.raises(ValueError, match="unknown variant"):
        fk_b(inst.G, inst.H, variant="nope")
    with pytest.raises(ValueError, match="unknown split rule"):
        fk_b(inst.G, inst.H, split_rule="nope")
    with pytest.raises(ValueError, match="both sides must share a ground set"):
        fk_b(Hypergraph.from_sets(3, [[0]]), Hypergraph.from_sets(4, [[0]]))


def test_max_nodes_guard():
    inst = load("trivial-aut-1")
    with pytest.raises(RuntimeError, match="max_nodes"):
        fk_b(inst.G, inst.H, max_nodes=5)


def test_irredundant_switch_is_recorded_on_the_node():
    inst = load("f2g2")
    reduced = fk_b(inst.G, inst.H)
    assert any(n.removed_G or n.removed_H for n in reduced), (
        "no node needed reduction, so the switch below proves nothing"
    )
    assert all(not (n.removed_G or n.removed_H) for n in fk_b(inst.G, inst.H, irredundant=False))


# ----------------------------------------------------------------------
# output
# ----------------------------------------------------------------------
def test_report_and_dot_are_labelled_fk_b():
    from fka.dot import to_dot
    from fka.report import render_html

    inst = load("6b")
    tree = fk_b(inst.G, inst.H, instance="6b")
    page = render_html(tree, instance=inst)
    assert "FK-B" in page
    assert '"a":"C"' in page and '"b":"D"' in page
    assert "conflicting assignment" in page

    src = to_dot(tree)
    assert src.startswith("digraph FKB {")
    assert "FK-B" in src
    assert "|C|=" in src
    assert "style=dotted" in src  # the per-clause branch


def test_cli_run_and_compare(tmp_path, capsys):
    from fkb.cli import main

    assert main(["run", "k3", "--out", str(tmp_path), "--no-annotate", "--dot"]) == 0
    for ext in ("html", "json", "dot"):
        assert (tmp_path / f"k3.fk-b.faithful.{ext}").exists()

    capsys.readouterr()
    assert main(["compare", "fano", "6b", "k3"]) == 0
    out = capsys.readouterr().out
    assert "fk-a" in out and "fk-b" in out
    assert "agree on all of them" in out


def test_cli_without_target_explains(capsys):
    from fkb.cli import main

    assert main(["run"]) == 2
    assert "--all" in capsys.readouterr().out
