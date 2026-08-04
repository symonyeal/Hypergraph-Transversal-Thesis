"""Dualisation by repeated FK-B equivalence tests (``FK_Dualization.m``)."""

from __future__ import annotations

import pytest

from fka.hypergraph import Hypergraph
from fka.instances import load, load_all, matching, self_dual_fano, threshold
from fka.transversal import transversal
from fkb.algorithm import dnf_value
from fkb.dualize import dualise, initial_cnf, maximum_false_point


@pytest.mark.parametrize("instance_id", [i.id for i in load_all()])
@pytest.mark.parametrize("side", ["G", "H"])
def test_reproduces_the_oracle_transversal(instance_id, side):
    """The whole point: FK-B as a generator must produce exactly ``Tr(H)``.

    Checked against Berge's method, which shares no code with this route.
    """
    H = getattr(load(instance_id), side)
    run = dualise(H)
    assert run.complete
    assert run.cnf.edges == transversal(H).edges


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
    H = load("trivial-aut-1").H
    run = dualise(H)
    assert run.passes == len(transversal(H))
    assert run.sizes == list(range(1, run.passes + 1))


def test_maximum_false_point_is_maximal_and_false():
    """Its complement must be a *minimal* transversal, or a clause repeats."""
    H = load("fano").H
    S = 0  # the empty assignment is false for any monomial-bearing DNF
    mfp = maximum_false_point(H, S)
    assert not dnf_value(H, mfp)
    # Maximal: turning on any further variable makes the DNF true.
    for v in range(H.n):
        if not (mfp >> v) & 1:
            assert dnf_value(H, mfp | (1 << v))
    assert (H.support() & ~mfp) in transversal(H).edges


def test_maximum_false_point_keeps_the_seed_true():
    H = load("6b").H
    S = next(s for s in range(1 << H.n) if not dnf_value(H, s) and s)
    mfp = maximum_false_point(H, S)
    assert mfp & S == S, "the seed's variables must survive"


def test_maximum_false_point_rejects_a_true_seed():
    H = load("fano").H
    with pytest.raises(ValueError, match="false point"):
        maximum_false_point(H, H.edges[0])


def test_initial_cnf_is_made_of_real_transversals():
    H = load("7ver").H
    truth = transversal(H)
    start = initial_cnf(H, n_clause=3)
    assert start.edges, "the start must not be empty"
    for c in start.edges:
        assert c in truth.edges


def test_partial_run_is_a_subset_never_a_wrong_answer():
    H = load("trivial-aut-1").H
    truth = set(transversal(H).edges)
    run = dualise(H, max_passes=5)
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
    H = load("f2g2").H
    run = dualise(H, variant=variant)
    assert run.complete
    assert run.cnf.edges == transversal(H).edges


@pytest.mark.slow
def test_reproduces_the_oracle_on_random_instances(rng):
    from conftest import random_sperner

    checked = 0
    for _ in range(120):
        n = rng.randint(2, 7)
        H = random_sperner(n, rng.randint(1, 5), rng)
        if not H.edges:
            continue
        run = dualise(H)
        assert run.complete
        assert run.cnf.edges == transversal(H).edges, H.label()
        checked += 1
    assert checked > 80
