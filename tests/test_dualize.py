"""Dualisation by repeated FK-B equivalence tests (``FK_Dualization.m``)."""

from __future__ import annotations

import pytest

from fka.hypergraph import Hypergraph
from fka.instances import load, load_all, matching, self_dual_fano, threshold
from fka.transversal import transversal
from fkb.algorithm import dnf_value
from fkb.dualize import dualise, initial_cnf, maximum_false_point


@pytest.mark.parametrize("instance_id", [i.id for i in load_all()])
@pytest.mark.parametrize("side", ["cnf", "dnf"])
def test_reproduces_the_oracle_transversal(instance_id, side):
    """The whole point: FK-B as a generator must produce exactly ``Tr(D)``.

    Checked against Berge's method, which shares no code with this route.
    """
    D = getattr(load(instance_id), side)
    run = dualise(D)
    assert run.complete
    assert run.cnf.edges == transversal(D).edges


@pytest.mark.parametrize(
    "family",
    [matching(6), matching(8), threshold(6), threshold(8), self_dual_fano(1)],
    ids=["M(6)", "M(8)", "TH(6)", "TH(8)", "SDFP(1)"],
)
def test_reproduces_the_oracle_on_the_benchmark_families(family):
    run = dualise(family)
    assert run.complete
    assert run.cnf.edges == transversal(family).edges


def test_one_pass_per_output_edge():
    """Each pass adds a distinct member of ``Tr(D)``, so passes == |Tr(D)|.

    That is the bound the self-reduction promises, and it is what makes the
    loop terminate rather than merely appear to.
    """
    D = load("trivial-aut-1").dnf
    run = dualise(D)
    assert run.passes == len(transversal(D))
    assert run.sizes == list(range(1, run.passes + 1))


def test_maximum_false_point_is_maximal_and_false():
    """Its complement must be a *minimal* transversal, or a clause repeats."""
    D = load("fano").dnf
    S = 0  # the empty assignment is false for any monomial-bearing DNF
    mfp = maximum_false_point(D, S)
    assert not dnf_value(D, mfp)
    # Maximal: turning on any further variable makes the DNF true.
    for v in range(D.n):
        if not (mfp >> v) & 1:
            assert dnf_value(D, mfp | (1 << v))
    assert (D.support() & ~mfp) in transversal(D).edges


def test_maximum_false_point_keeps_the_seed_true():
    D = load("6b").dnf
    S = next(s for s in range(1 << D.n) if not dnf_value(D, s) and s)
    mfp = maximum_false_point(D, S)
    assert mfp & S == S, "the seed's variables must survive"


def test_maximum_false_point_rejects_a_true_seed():
    D = load("fano").dnf
    with pytest.raises(ValueError, match="false point"):
        maximum_false_point(D, D.edges[0])


def test_initial_cnf_is_made_of_real_transversals():
    D = load("7ver").dnf
    truth = transversal(D)
    start = initial_cnf(D, n_clause=3)
    assert start.edges, "the start must not be empty"
    for c in start.edges:
        assert c in truth.edges


def test_partial_run_is_a_subset_never_a_wrong_answer():
    D = load("trivial-aut-1").dnf
    truth = set(transversal(D).edges)
    run = dualise(D, max_passes=5)
    assert not run.complete
    assert run.passes == 5
    assert set(run.cnf.edges) <= truth


def test_degenerate_inputs():
    """``Tr({}) = {{}}`` and ``Tr({{}}) = {}``, as elsewhere in the project."""
    empty = dualise(Hypergraph.empty(3))
    assert empty.complete and empty.cnf.edges == (0,)
    one = dualise(Hypergraph(3, (0,)))
    assert one.complete and one.cnf.edges == ()


@pytest.mark.parametrize("variant", ["faithful", "multiple"])
def test_both_variants_dualise_identically(variant):
    D = load("f2g2").dnf
    run = dualise(D, variant=variant)
    assert run.complete
    assert run.cnf.edges == transversal(D).edges


@pytest.mark.slow
def test_reproduces_the_oracle_on_random_instances(rng):
    from conftest import random_sperner

    checked = 0
    for _ in range(120):
        n = rng.randint(2, 7)
        D = random_sperner(n, rng.randint(1, 5), rng)
        if not D.edges:
            continue
        run = dualise(D)
        assert run.complete
        assert run.cnf.edges == transversal(D).edges, D.label()
        checked += 1
    assert checked > 80
