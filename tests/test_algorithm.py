"""FK-A: agreement with the oracle, thesis reproduction, and legacy-bug regressions."""

from __future__ import annotations

import pytest

from fka.algorithm import (
    PIVOT_RULES,
    VARIANTS,
    check_preconditions,
    choose_pivot,
    epsilon,
    fk_a,
    is_dual,
    verify_certificate,
)
from fka.hypergraph import Hypergraph
from fka.instances import exact_epsilon, load, load_all
from fka.transversal import is_dual_oracle, transversal


# ----------------------------------------------------------------------
# regressions against the legacy implementation's defects
# ----------------------------------------------------------------------
def test_recursion_actually_happens():
    """Regression: the legacy ``fk`` could never recurse.

    It ran ``if check_base(...): return True`` immediately followed by
    ``if not check_base(...): return False``, so every call returned at the
    base-case check and Step 4 was dead code.
    """
    inst = load("f2g2")
    tree = fk_a(inst.G, inst.H)
    assert len(tree) > 1
    assert tree.root.pivot is not None


def test_child_results_are_propagated():
    """Regression: the legacy ``fk`` discarded both recursive return values.

    It called ``fk(...)`` twice without binding the results and then
    unconditionally returned ``True``, so a non-dual instance deep in the tree
    was reported dual.
    """
    inst = load("f2g2")
    G, H = inst.G, inst.H
    broken = Hypergraph(H.n, H.edges[:-1])  # drop an edge: no longer dual
    assert not is_dual_oracle(G, broken)
    assert fk_a(G, broken).dual is False


def test_preconditions_handle_multi_vertex_edges():
    """Regression: the legacy ``check_pre`` raised on condition (iii).

    ``any(f[i] & g_row for g_row in g)`` builds numpy arrays and takes their
    truth value, which raises ValueError for arrays of length > 1.
    """
    G = Hypergraph.from_sets(4, [[0, 1], [2, 3]])
    H = transversal(G)
    assert check_preconditions(G, H) is None


def test_condition_iv_is_exact_not_floating_point():
    """A transversal pair must never be rejected by rounding in the 2^-|e| sum.

    ``Tr`` of a single k-edge is k singletons, so the sum is exactly
    ``2^-k + k * 2^-1``; for k = 1 that is exactly 1 and must pass.
    """
    G = Hypergraph.from_sets(1, [[0]])
    H = transversal(G)
    assert check_preconditions(G, H) is None
    assert fk_a(G, H).dual is True


# ----------------------------------------------------------------------
# specification conformance
# ----------------------------------------------------------------------
def test_pivot_rules_break_ties_to_lowest_index():
    """Fano is vertex-transitive, so every vertex ties and v1 must win."""
    inst = load("fano")
    for rule in PIVOT_RULES:
        assert choose_pivot(inst.G, inst.H, rule) == 0


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
        assert str(exact_epsilon(inst.G, inst.H)) == want, instance_id
        assert epsilon(inst.G, inst.H) == pytest.approx(
            float(exact_epsilon(inst.G, inst.H))
        )


def test_modified_variant_never_grows_the_tree():
    """The D_{1,k} rule removes recursion steps; it cannot add any.

    Thesis Section 5.1.1 (p.47): each such node would otherwise need up to
    ``k-1`` further recursions to reach D_{1,1}.
    """
    for inst in load_all():
        faithful = fk_a(inst.G, inst.H, variant="faithful")
        modified = fk_a(inst.G, inst.H, variant="modified")
        assert len(modified) <= len(faithful), inst.id
        assert modified.dual == faithful.dual, inst.id


def test_variant_and_pivot_rule_are_validated():
    inst = load("k3")
    with pytest.raises(ValueError, match="unknown variant"):
        fk_a(inst.G, inst.H, variant="nope")
    with pytest.raises(ValueError, match="unknown pivot rule"):
        fk_a(inst.G, inst.H, pivot_rule="nope")


def test_mismatched_ground_sets_are_rejected():
    with pytest.raises(ValueError, match="both sides must share a ground set"):
        fk_a(Hypergraph.from_sets(3, [[0]]), Hypergraph.from_sets(4, [[0]]))


def test_max_nodes_guard():
    inst = load("trivial-aut-1")
    with pytest.raises(RuntimeError, match="max_nodes"):
        fk_a(inst.G, inst.H, max_nodes=5)


# ----------------------------------------------------------------------
# correctness against the independent oracle
# ----------------------------------------------------------------------
@pytest.mark.parametrize("instance_id", [i.id for i in load_all()])
@pytest.mark.parametrize("variant", VARIANTS)
@pytest.mark.parametrize("pivot_rule", PIVOT_RULES)
def test_agrees_with_oracle_on_library(instance_id, variant, pivot_rule):
    inst = load(instance_id)
    truth = is_dual_oracle(inst.G, inst.H)
    tree = fk_a(
        inst.G, inst.H, variant=variant, pivot_rule=pivot_rule, validate=True
    )
    assert tree.dual == truth


@pytest.mark.slow
def test_agrees_with_oracle_on_random_instances(rng):
    """The real correctness check: many random pairs, both variants, both rules.

    Most random pairs are *not* dual, so genuinely-dual cases are seeded by
    building ``H = Tr(G)`` half the time.
    """
    from conftest import random_sperner

    checked_dual = 0
    for _ in range(400):
        n = rng.randint(2, 7)
        G = random_sperner(n, rng.randint(1, 5), rng)
        if not G.edges:
            continue
        if rng.random() < 0.5:
            H = transversal(G)
            checked_dual += 1
        else:
            H = random_sperner(n, rng.randint(1, 5), rng)
        truth = is_dual_oracle(G, H)
        for variant in VARIANTS:
            for rule in PIVOT_RULES:
                tree = fk_a(G, H, variant=variant, pivot_rule=rule, validate=True)
                assert tree.dual == truth, (
                    f"n={n} G={G.label()} H={H.label()} "
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
        G = random_sperner(n, rng.randint(1, 5), rng)
        H = random_sperner(n, rng.randint(1, 5), rng)
        if not G.edges or not H.edges:
            continue
        tree = fk_a(G, H)
        if tree.dual or tree.verdict is None or tree.verdict.certificate is None:
            continue
        assert verify_certificate(G, H, tree.verdict.certificate), (
            f"bad certificate for G={G.label()} H={H.label()}"
        )
        checked += 1
    assert checked > 40, "not enough certificates were produced to be meaningful"


def test_is_dual_convenience_wrapper():
    inst = load("fano")
    assert is_dual(inst.G, inst.H) is True
