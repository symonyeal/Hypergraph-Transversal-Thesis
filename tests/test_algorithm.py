"""FK-A: agreement with the oracle, and the invariants each step must hold."""

from __future__ import annotations

import pytest

from fka.algorithm import (
    SPLIT_RULES,
    VARIANTS,
    check_conditions,
    choose_split_var,
    eps,
    fk_a,
    is_dual,
    verify_CA,
)
from fka.hypergraph import Hypergraph
from fka.algorithm import eps_exact
from fka.instances import load, load_all
from fka.transversal import is_dual_oracle, transversal


# ----------------------------------------------------------------------
# the invariants each step must hold
# ----------------------------------------------------------------------
def test_recursion_actually_happens():
    """The recursion must actually reach Step 4 on a non-trivial instance.

    It ran ``if check_base(...): return True`` immediately followed by
    ``if not check_base(...): return False``, so every call returned at the
    base-case check and Step 4 was dead code.
    """
    inst = load("f2g2")
    tree = fk_a(inst.cnf, inst.dnf)
    assert len(tree) > 1
    assert tree.root.split_var is not None


def test_child_results_are_propagated():
    """A child's non-dual verdict must reach the root, not be discarded."""
    inst = load("f2g2")
    cnf, dnf = inst.cnf, inst.dnf
    broken = Hypergraph(dnf.n, dnf.edges[:-1])  # drop a term: no longer dual
    assert not is_dual_oracle(cnf, broken)
    assert fk_a(cnf, broken).dual is False


def test_preconditions_handle_multi_vertex_edges():
    """Condition (iii) must handle terms of more than one variable.

    ``any(f[i] & g_row for g_row in g)`` builds numpy arrays and takes their
    truth value, which raises ValueError for arrays of length > 1.
    """
    cnf = Hypergraph.from_sets(4, [[0, 1], [2, 3]])
    MBF = transversal(cnf)
    assert check_conditions(cnf, MBF) is None


def test_condition_iv_is_exact_not_floating_point():
    """A transversal pair must never be rejected by rounding in the 2^-|e| sum.

    ``Tr`` of a single k-edge is k singletons, so the sum is exactly
    ``2^-k + k * 2^-1``; for k = 1 that is exactly 1 and must pass.
    """
    cnf = Hypergraph.from_sets(1, [[0]])
    MBF = transversal(cnf)
    assert check_conditions(cnf, MBF) is None
    assert fk_a(cnf, MBF).dual is True


# ----------------------------------------------------------------------
# specification conformance
# ----------------------------------------------------------------------
def test_split_rules_break_ties_to_lowest_index():
    """Fano is vertex-transitive, so every vertex ties and v1 must win."""
    inst = load("fano")
    for rule in SPLIT_RULES:
        assert choose_split_var(inst.cnf, inst.dnf, rule) == 0


def test_epsilon_matches_the_thesis_values():
    """Thesis Section 5.2: Fano 3/7, 6-A..6-D 1/2, 7-Ver 5/11, 8-Ver 2/5."""
    expected = {
        "fano": "3/7",
        "k3": "2/3",
        "6a": "1/2",
        "6b": "1/2",
        "6c": "1/2",
        "6d": "1/2",
        "7ver": "5/11",
        "8ver": "2/5",
    }
    for instance_id, want in expected.items():
        inst = load(instance_id)
        assert str(eps_exact(inst.cnf, inst.dnf)) == want, instance_id
        assert eps(inst.cnf, inst.dnf) == pytest.approx(
            float(eps_exact(inst.cnf, inst.dnf))
        )


def test_modified_variant_never_grows_the_tree():
    """The D_{1,k} rule removes recursion steps; it cannot add any.

    Thesis Section 5.1.1 (p.47): each such node would otherwise need up to
    ``k-1`` further recursions to reach D_{1,1}.
    """
    for inst in load_all():
        faithful = fk_a(inst.cnf, inst.dnf, variant="faithful")
        modified = fk_a(inst.cnf, inst.dnf, variant="modified")
        assert len(modified) <= len(faithful), inst.id
        assert modified.dual == faithful.dual, inst.id


def test_variant_and_split_rule_are_validated():
    inst = load("k3")
    with pytest.raises(ValueError, match="unknown variant"):
        fk_a(inst.cnf, inst.dnf, variant="nope")
    with pytest.raises(ValueError, match="unknown split rule"):
        fk_a(inst.cnf, inst.dnf, split_rule="nope")


def test_mismatched_ground_sets_are_rejected():
    with pytest.raises(ValueError, match="both sides must share a ground set"):
        fk_a(Hypergraph.from_sets(3, [[0]]), Hypergraph.from_sets(4, [[0]]))


def test_max_nodes_guard():
    inst = load("trivial-aut-1")
    with pytest.raises(RuntimeError, match="max_nodes"):
        fk_a(inst.cnf, inst.dnf, max_nodes=5)


# ----------------------------------------------------------------------
# correctness against the independent oracle
# ----------------------------------------------------------------------
@pytest.mark.parametrize("instance_id", [i.id for i in load_all()])
@pytest.mark.parametrize("variant", VARIANTS)
@pytest.mark.parametrize("split_rule", SPLIT_RULES)
def test_agrees_with_oracle_on_library(instance_id, variant, split_rule):
    inst = load(instance_id)
    truth = is_dual_oracle(inst.cnf, inst.dnf)
    tree = fk_a(
        inst.cnf, inst.dnf, variant=variant, split_rule=split_rule, validate=True
    )
    assert tree.dual == truth


@pytest.mark.slow
def test_agrees_with_oracle_on_random_instances(rng):
    """The real correctness check: many random pairs, both variants, both rules.

    Most random pairs are *not* dual, so genuinely-dual cases are seeded by
    building ``MBF = Tr(cnf)`` half the time.
    """
    from conftest import random_sperner

    checked_dual = 0
    for _ in range(400):
        n = rng.randint(2, 7)
        cnf = random_sperner(n, rng.randint(1, 5), rng)
        if not cnf.edges:
            continue
        if rng.random() < 0.5:
            MBF = transversal(cnf)
            checked_dual += 1
        else:
            MBF = random_sperner(n, rng.randint(1, 5), rng)
        truth = is_dual_oracle(cnf, MBF)
        for variant in VARIANTS:
            for rule in SPLIT_RULES:
                tree = fk_a(cnf, MBF, variant=variant, split_rule=rule, validate=True)
                assert tree.dual == truth, (
                    f"n={n} cnf={cnf.label()} MBF={MBF.label()} "
                    f"variant={variant} rule={rule}"
                )
    assert checked_dual > 50, "not enough genuinely dual cases were exercised"


@pytest.mark.slow
def test_certificates_certify_the_root_instance(rng):
    """A non-dual answer's certificate must satisfy equation 2.1 at the *root*.

    The recursion lifts certificates up through pivots and through Sperner
    reduction; this checks the lift end to end rather than node-locally.
    """
    from conftest import random_sperner

    checked = 0
    for _ in range(400):
        n = rng.randint(2, 7)
        cnf = random_sperner(n, rng.randint(1, 5), rng)
        MBF = random_sperner(n, rng.randint(1, 5), rng)
        if not cnf.edges or not MBF.edges:
            continue
        tree = fk_a(cnf, MBF)
        if tree.dual or tree.verdict is None or tree.verdict.CA is None:
            continue
        assert verify_CA(cnf, MBF, tree.verdict.CA), (
            f"bad certificate for cnf={cnf.label()} MBF={MBF.label()}"
        )
        checked += 1
    assert checked > 40, "not enough certificates were produced to be meaningful"


def test_is_dual_convenience_wrapper():
    inst = load("fano")
    assert is_dual(inst.cnf, inst.dnf) is True
