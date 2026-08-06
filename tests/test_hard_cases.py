"""W11 = S(4,5,11), and the boundary of the family it belongs to.

The Fano plane and W11 are the same construction on two different perfect
codes -- Hamming [7,4,3] and ternary Golay [11,6,5] -- and both are
self-transversal. The third perfect code, the binary Golay [23,12,7], gives
S(4,7,23), which is intersecting but is *not* self-transversal; that boundary
is asserted here with an explicit certificate, because it is the reason
``docs/hard-cases.md`` claims a family of two rather than a family of three.
"""

from __future__ import annotations

from itertools import combinations, product

import pytest

from fka.algorithm import fk_a
from fka.automorphism import automorphism_group
from fka.hypergraph import Hypergraph, bits, popcount, verts
from fka.instances import (
    FANO_LINES,
    exact_epsilon,
    load,
    self_dual_fano,
    self_dualise,
    witt_11,
)
from fka.properties import analyse
from fka.transversal import is_dual_oracle, transversal
from fkb.algorithm import fk_b


# ----------------------------------------------------------------------
# the design
# ----------------------------------------------------------------------
def test_witt_11_is_the_steiner_system_s4_5_11():
    W = witt_11()
    assert W.n == 11
    assert len(W) == 66
    assert set(W.edge_sizes()) == {5}
    assert set(W.degrees()) == {30}
    # every 4-subset of the 11 points lies in exactly one block
    for quad in combinations(range(11), 4):
        q = bits(quad)
        assert sum(1 for e in W.edges if e & q == q) == 1


def test_witt_11_is_self_transversal():
    W = witt_11()
    assert transversal(W).edges == W.edges
    assert is_dual_oracle(W, W)


def test_witt_11_blocks_pairwise_meet():
    """Intersecting is necessary for self-transversality, and here the blocks
    meet in 1, 2 or 3 points -- never 0, and never 4 (that would repeat a
    quadruple)."""
    W = witt_11()
    sizes = {popcount(a & b) for i, a in enumerate(W.edges) for b in W.edges[i + 1:]}
    assert sizes == {1, 2, 3}


def test_witt_11_frequency_and_symmetry():
    W = witt_11()
    assert str(exact_epsilon(W, W)) == "5/11"
    aut = automorphism_group(W)
    assert aut.order == 7920                     # M11
    assert len(aut.orbits) == 1                  # vertex-transitive
    props = analyse(W, graph_classes=False)
    assert not props.conformal                   # as with the Fano plane
    assert not props.read_once


@pytest.mark.sage
def test_witt_11_group_is_m11_by_sage():
    from fka.backends.sage_backend import sage_automorphism_group

    result = sage_automorphism_group(witt_11())
    assert result["order"] == 7920
    assert result["name"] == "M11"


# ----------------------------------------------------------------------
# both algorithms, against the oracle
# ----------------------------------------------------------------------
@pytest.mark.parametrize("instance_id", ["w11", "w11-sd-1"])
def test_new_instances_agree_with_the_oracle(instance_id):
    inst = load(instance_id)
    G, H = inst.G, inst.H
    assert is_dual_oracle(G, H)
    for variant in ("faithful", "modified"):
        assert fk_a(G, H, variant=variant).verdict.dual
    for variant in ("faithful", "multiple"):
        assert fk_b(G, H, variant=variant).verdict.dual


def test_w11_is_harder_than_fano_for_both_algorithms():
    """The point of the instance. Per unit of input -- nodes over |G| + |H| --
    W11 costs both algorithms more than the Fano plane does, and the margin is
    largest for FK-B, which the Fano plane barely troubles."""
    fano, w11 = load("fano"), load("w11")
    for algorithm in (fk_a, fk_b):
        f = len(algorithm(fano.G, fano.H)) / (len(fano.G_edges) + len(fano.H_edges))
        w = len(algorithm(w11.G, w11.H)) / (len(w11.G_edges) + len(w11.H_edges))
        assert w > f


# ----------------------------------------------------------------------
# where the family stops
# ----------------------------------------------------------------------
#: ``x^11 + x^10 + x^6 + x^5 + x^4 + x^2 + 1``, a factor of ``x^23 - 1`` over
#: GF(2); the cyclic code it generates is the perfect [23,12,7] binary Golay
#: code, whose minimum-weight supports are the 253 blocks of S(4,7,23).
GOLAY2_GENERATOR = (1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1)

#: A subset of the 23 points containing no block and whose complement contains
#: no block: a certificate against ``Tr(S(4,7,23)) = S(4,7,23)`` in the form of
#: thesis equation 2.1. Found by ``experiments/hard_case_search.py``'s exact
#: test; 16744 of the 1352078 11-subsets are certificates.
GOLAY2_CERTIFICATE = (1, 2, 3, 4, 5, 6, 7, 8, 10, 14, 21)


def _witt_23() -> Hypergraph:
    n = 23
    degree = len(GOLAY2_GENERATOR) - 1
    blocks = []
    for message in product(range(2), repeat=n - degree):
        word = [0] * n
        for i, coefficient in enumerate(message):
            if coefficient:
                for j, g in enumerate(GOLAY2_GENERATOR):
                    word[(i + j) % n] ^= g
        support = [i for i, x in enumerate(word) if x]
        if len(support) == 7:
            blocks.append(bits(support))
    return Hypergraph.from_bitsets(n, blocks)


def test_binary_golay_design_is_intersecting_but_not_self_transversal():
    W = _witt_23()
    assert len(W) == 253 and set(W.edge_sizes()) == {7} and set(W.degrees()) == {77}
    for quad in combinations(range(9), 4):       # a sample: S(4,7,23) is a 4-design
        q = bits(quad)
        assert sum(1 for e in W.edges if e & q == q) == 1
    assert all(a & b for i, a in enumerate(W.edges) for b in W.edges[i + 1:])

    S = bits(v - 1 for v in GOLAY2_CERTIFICATE)
    complement = ((1 << 23) - 1) & ~S
    # Equation 2.1: S contains no edge, and neither does its complement. That
    # alone settles Tr(W) != W -- if W were self-transversal, one side of every
    # partition would contain a block -- and it settles it without building
    # Tr(W), which for 253 blocks on 23 points is not a bounded computation.
    assert not any((e & S) == e for e in W.edges)
    assert not any((e & complement) == e for e in W.edges)


# ----------------------------------------------------------------------
# the generalised construction still reproduces the Fano benchmark
# ----------------------------------------------------------------------
@pytest.mark.parametrize("k", [1, 2, 3])
def test_self_dualise_reproduces_self_dual_fano(k):
    lines = [list(line) for line in FANO_LINES]
    assert self_dualise(lines, 7, k).edges == self_dual_fano(k).edges


def test_self_dualised_w11_is_self_transversal():
    W = witt_11()
    S = self_dualise(W.to_sets(one_indexed=True), 11, 1)
    assert S.n == 13
    assert len(S) == 2 * 66 + 1
    assert is_dual_oracle(S, S)
    assert sorted(set(popcount(e) for e in S.edges)) == [2, 6]
    assert [len(verts(e)) for e in S.edges].count(2) == 1
