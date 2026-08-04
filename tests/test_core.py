"""Hypergraph representation, Sperner reduction, and the transversal oracle."""

from __future__ import annotations

import pytest

from fka.hypergraph import Hypergraph, verts, popcount
from fka.sperner import (
    sperner_reduce,
    remove_superset_pairwise,
    remove_superset_sorted,
)
from fka.transversal import find_certificate, is_dual_oracle, minimal_elements, transversal


# ----------------------------------------------------------------------
# representation
# ----------------------------------------------------------------------
def test_incidence_round_trip():
    matrix = [
        [1, 0, 1, 0],
        [0, 1, 0, 1],
        [1, 1, 0, 0],
    ]
    H = Hypergraph.from_incidence(matrix)
    assert H.n == 4
    assert sorted(H.to_incidence()) == sorted(matrix)


def test_sets_round_trip_one_indexed():
    H = Hypergraph.from_sets(6, [[1, 5, 6], [4, 5], [1, 2]], one_indexed=True)
    assert H.to_sets(one_indexed=True) == [[1, 2], [1, 5, 6], [4, 5]]


def test_edges_are_deduplicated_and_sorted():
    H = Hypergraph.from_sets(3, [[0, 1], [1, 0], [2]])
    assert len(H) == 2


def test_rejects_vertex_out_of_range():
    with pytest.raises(ValueError, match="out of range"):
        Hypergraph.from_sets(3, [[1, 9]], one_indexed=True)


def test_rejects_ragged_incidence():
    with pytest.raises(ValueError, match="ragged"):
        Hypergraph.from_incidence([[1, 0], [1, 0, 1]])


def test_contraction_and_deletion_match_thesis_example():
    """Thesis Example 2.2.7 (p.15).

    ``V = {v1..v4}``, ``E = {{v1,v2},{v2,v3,v4},{v1,v3}}``, contracting v1 gives
    ``{{v2},{v2,v3,v4},{v3}}`` and deleting it gives ``{{v2,v3,v4}}``.
    """
    H = Hypergraph.from_sets(4, [[1, 2], [2, 3, 4], [1, 3]], one_indexed=True)
    assert H.contraction(0).to_sets(one_indexed=True) == [[2], [3]]
    assert H.deletion(0).to_sets(one_indexed=True) == [[2, 3, 4]]


def test_deletion_preserves_sperner_contraction_need_not():
    """Thesis Section 2.2.2 (p.15): H\\v stays Sperner, H/v need not."""
    H = Hypergraph.from_sets(4, [[1, 2], [2, 3, 4], [1, 3]], one_indexed=True)
    assert H.is_sperner()
    assert H.deletion(0).is_sperner()
    contracted = H.contraction(0).union(H.deletion(0))
    assert not contracted.is_sperner()  # {v2} is inside {v2,v3,v4}


def test_isolated_vertices_and_compact():
    H = Hypergraph.from_sets(6, [[1, 3]], one_indexed=True)
    assert H.isolated_vertices() == (1, 3, 4, 5)
    small, mapping = H.compact()
    assert small.n == 2
    assert small.to_sets() == [[0, 1]]
    assert mapping == [0, 2]


def test_popcount_and_bitset_helpers():
    e = 0b101101
    assert popcount(e) == 4
    assert verts(e) == (0, 2, 3, 5)


# ----------------------------------------------------------------------
# Sperner reduction
# ----------------------------------------------------------------------
def test_superset_removal_drops_the_all_ones_row():
    """Regression: the legacy ``remove_supersets`` skipped rows with no zeros.

    It looked only at the zero columns of each row, so a row of all ones had
    nothing to compare and survived -- despite being a superset of everything.
    """
    H = Hypergraph.from_incidence([[1, 1, 1], [1, 0, 0], [0, 1, 0]])
    reduced = sperner_reduce(H).reduced
    assert reduced.to_sets() == [[0], [1]]


def test_superset_removal_keeps_one_of_each_duplicate():
    """Duplicates must be removed *before* the superset pass.

    Two equal edges are each a (non-strict) superset of the other, so a
    superset pass run first would delete both.
    """
    H = Hypergraph.from_sets(3, [[0, 1], [0, 1], [2]])
    result = sperner_reduce(H)
    assert result.reduced.to_sets() == [[0, 1], [2]]


def test_reduction_reports_what_it_removed():
    H = Hypergraph.from_sets(4, [[0], [0, 1], [0, 1, 2], [3]])
    result = sperner_reduce(H)
    assert result.reduced.to_sets() == [[0], [3]]
    assert len(result.removed_superset) == 2
    assert result.changed


def test_chain_reduces_to_the_minimum():
    H = Hypergraph.from_sets(4, [[0, 1, 2, 3], [0, 1, 2], [0, 1], [0]])
    assert sperner_reduce(H).reduced.to_sets() == [[0]]


@pytest.mark.slow
def test_both_superset_implementations_agree(rng):
    """The optimised pass must equal the definition on every input."""
    from conftest import random_sperner  # noqa: F401  (import for symmetry)

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
        a_keep, a_drop = remove_superset_pairwise(edges)
        b_keep, b_drop = remove_superset_sorted(edges)
        assert set(a_keep) == set(b_keep)
        assert set(a_drop) == set(b_drop)


@pytest.mark.slow
def test_reduction_is_idempotent_and_sperner(rng):
    for _ in range(200):
        n = rng.randint(1, 8)
        H = Hypergraph.from_sets(
            n,
            [rng.sample(range(n), rng.randint(1, n)) for _ in range(rng.randint(1, 8))],
        )
        once = sperner_reduce(H).reduced
        assert once.is_sperner()
        assert sperner_reduce(once).reduced == once


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
    H = Hypergraph.from_sets(5, [[0, 2, 4]])
    assert transversal(H).to_sets() == [[0], [2], [4]]


def test_fano_is_self_transversal():
    from fka.instances import FANO_LINES

    fano = Hypergraph.from_sets(7, [list(l) for l in FANO_LINES], one_indexed=True)
    assert transversal(fano).edges == fano.edges


@pytest.mark.slow
def test_transversal_is_an_involution_on_sperner_hypergraphs(rng):
    """``Tr(Tr(H)) = H`` for every Sperner ``H`` -- the duality the problem is about."""
    from conftest import random_sperner

    for _ in range(120):
        n = rng.randint(1, 7)
        H = random_sperner(n, rng.randint(1, 6), rng)
        if not H.edges:
            continue
        assert transversal(transversal(H)).edges == H.edges


def test_find_certificate_agrees_with_duality():
    G = Hypergraph.from_sets(4, [[0, 1], [2, 3]])
    H = transversal(G)
    assert find_certificate(G, H) is None
    broken = Hypergraph(H.n, H.edges[:-1])
    S = find_certificate(G, broken)
    assert S is not None
    assert all(S & e for e in G.edges)
    assert all((e & S) != e for e in broken.edges)
