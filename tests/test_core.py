"""The term-set model, ``Irredundant``, and the transversal oracle."""

from __future__ import annotations

import pytest

from fka.hypergraph import Hypergraph, phi_x_0, phi_x_1, popcount, verts
from fka.irredundant import irredundant, remove_superset
from fka.transversal import find_CA, is_dual_oracle, minimal_elements, transversal


# ----------------------------------------------------------------------
# representation
# ----------------------------------------------------------------------
def test_incidence_round_trip():
    matrix = [
        [1, 0, 1, 0],
        [0, 1, 0, 1],
        [1, 1, 0, 0],
    ]
    MBF = Hypergraph.from_incidence(matrix)
    assert MBF.n == 4
    assert sorted(MBF.to_incidence()) == sorted(matrix)


def test_sets_round_trip_one_indexed():
    MBF = Hypergraph.from_sets(6, [[1, 5, 6], [4, 5], [1, 2]], one_indexed=True)
    assert MBF.to_sets(one_indexed=True) == [[1, 2], [1, 5, 6], [4, 5]]


def test_edges_are_deduplicated_and_sorted():
    MBF = Hypergraph.from_sets(3, [[0, 1], [1, 0], [2]])
    assert len(MBF) == 2


def test_rejects_vertex_out_of_range():
    with pytest.raises(ValueError, match="out of range"):
        Hypergraph.from_sets(3, [[1, 9]], one_indexed=True)


def test_rejects_ragged_incidence():
    with pytest.raises(ValueError, match="ragged"):
        Hypergraph.from_incidence([[1, 0], [1, 0, 1]])


def test_phi_x_matches_thesis_example():
    """Thesis Example 2.2.7 (p.15).

    ``V = {v1..v4}`` with terms ``{{v1,v2},{v2,v3,v4},{v1,v3}}``: ``phi_x_0`` on
    v1 gives ``{{v2},{v3}}`` and ``phi_x_1`` gives ``{{v2,v3,v4}}``.
    """
    MBF = Hypergraph.from_sets(4, [[1, 2], [2, 3, 4], [1, 3]], one_indexed=True)
    assert phi_x_0(MBF, 0).to_sets(one_indexed=True) == [[2], [3]]
    assert phi_x_1(MBF, 0).to_sets(one_indexed=True) == [[2, 3, 4]]


def test_phi_x_1_stays_irredundant_phi_x_0_need_not():
    """Thesis Section 2.2.2 (p.15): MBF\\v stays Sperner, MBF/v need not."""
    MBF = Hypergraph.from_sets(4, [[1, 2], [2, 3, 4], [1, 3]], one_indexed=True)
    assert MBF.is_sperner()
    assert phi_x_1(MBF, 0).is_sperner()
    both = phi_x_0(MBF, 0).vee(phi_x_1(MBF, 0))
    assert not both.is_sperner()  # {v2} is inside {v2,v3,v4}


def test_isolated_vertices_and_compact():
    MBF = Hypergraph.from_sets(6, [[1, 3]], one_indexed=True)
    assert MBF.isolated_vertices() == (1, 3, 4, 5)
    small, mapping = MBF.compact()
    assert small.n == 2
    assert small.to_sets() == [[0, 1]]
    assert mapping == [0, 2]


def test_popcount_and_bitset_helpers():
    e = 0b101101
    assert popcount(e) == 4
    assert verts(e) == (0, 2, 3, 5)


# ----------------------------------------------------------------------
# Irredundant
# ----------------------------------------------------------------------
def test_superset_removal_drops_the_all_ones_row():
    """A term containing every variable is a superset of all the others."""
    MBF = Hypergraph.from_incidence([[1, 1, 1], [1, 0, 0], [0, 1, 0]])
    reduced = irredundant(MBF).reduced
    assert reduced.to_sets() == [[0], [1]]


def test_superset_removal_keeps_one_of_each_duplicate():
    """Duplicates must be removed *before* the superset pass.

    Two equal terms are each a (non-strict) superset of the other, so a
    superset pass run first would delete both.
    """
    MBF = Hypergraph.from_sets(3, [[0, 1], [0, 1], [2]])
    result = irredundant(MBF)
    assert result.reduced.to_sets() == [[0, 1], [2]]


def test_reduction_reports_what_it_removed():
    MBF = Hypergraph.from_sets(4, [[0], [0, 1], [0, 1, 2], [3]])
    result = irredundant(MBF)
    assert result.reduced.to_sets() == [[0], [3]]
    assert len(result.removed_superset) == 2
    assert result.changed


def test_chain_reduces_to_the_minimum():
    MBF = Hypergraph.from_sets(4, [[0, 1, 2, 3], [0, 1, 2], [0, 1], [0]])
    assert irredundant(MBF).reduced.to_sets() == [[0]]


@pytest.mark.slow
def test_superset_removal_matches_the_definition(rng):
    """The size-ordered pass must equal the pairwise definition on every input.

    The definition lives here rather than in ``src``: it is what the optimised
    pass is checked against, not a second supported implementation.
    """

    def pairwise(edges):
        return tuple(
            a for i, a in enumerate(edges)
            if not any((b & a) == b for j, b in enumerate(edges) if i != j)
        )

    for _ in range(300):
        n = rng.randint(1, 9)
        edges = tuple(
            sorted(
                {
                    sum(1 << v for v in rng.sample(range(n), rng.randint(1, n)))
                    for _ in range(rng.randint(1, 10))
                }
            )
        )
        kept, dropped = remove_superset(edges)
        assert set(kept) == set(pairwise(edges))
        assert set(kept) | set(dropped) == set(edges)


@pytest.mark.slow
def test_reduction_is_idempotent_and_sperner(rng):
    for _ in range(200):
        n = rng.randint(1, 8)
        MBF = Hypergraph.from_sets(
            n,
            [rng.sample(range(n), rng.randint(1, n)) for _ in range(rng.randint(1, 8))],
        )
        once = irredundant(MBF).reduced
        assert once.is_sperner()
        assert irredundant(once).reduced == once


# ----------------------------------------------------------------------
# the transversal oracle
# ----------------------------------------------------------------------
def test_minimal_elements():
    assert minimal_elements([0b111, 0b011, 0b100, 0b011]) == (0b011, 0b100)


def test_degenerate_dual_conventions():
    """The edgeless hypergraph is constant 0; ``{{}}`` is constant 1; they are dual."""
    empty = Hypergraph.empty(3)
    unit = Hypergraph(3, (0,))
    assert transversal(empty).edges == unit.edges
    assert transversal(unit).edges == empty.edges
    assert is_dual_oracle(unit, empty)
    assert is_dual_oracle(empty, unit)
    assert not is_dual_oracle(empty, empty)


def test_transversal_of_a_single_edge_is_its_singletons():
    """``Tr({e})`` is the family of singletons of ``e`` -- the D_{1,k} base case."""
    MBF = Hypergraph.from_sets(5, [[0, 2, 4]])
    assert transversal(MBF).to_sets() == [[0], [2], [4]]


def test_fano_is_self_transversal():
    from fka.instances import FANO_LINES

    fano = Hypergraph.from_sets(7, [list(l) for l in FANO_LINES], one_indexed=True)
    assert transversal(fano).edges == fano.edges


@pytest.mark.slow
def test_transversal_is_an_involution_on_sperner_hypergraphs(rng):
    """``Tr(Tr(MBF)) = MBF`` for every Sperner ``MBF`` -- the duality the problem is about."""
    from conftest import random_sperner

    for _ in range(120):
        n = rng.randint(1, 7)
        MBF = random_sperner(n, rng.randint(1, 6), rng)
        if not MBF.edges:
            continue
        assert transversal(transversal(MBF)).edges == MBF.edges


def test_find_CA_agrees_with_duality():
    cnf = Hypergraph.from_sets(4, [[0, 1], [2, 3]])
    MBF = transversal(cnf)
    assert find_CA(cnf, MBF) is None
    broken = Hypergraph(MBF.n, MBF.edges[:-1])
    S = find_CA(cnf, broken)
    assert S is not None
    assert all(S & e for e in cnf.edges)
    assert all((e & S) != e for e in broken.edges)
