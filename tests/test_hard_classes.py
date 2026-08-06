"""Chapter 4's hard classes, as executable claims.

Every family here is built by substitution and its transversal is obtained by
composing duals, so the oracle in :mod:`fka.transversal` stays independent and
is what these assertions are checked against.
"""

from __future__ import annotations

import itertools
import math
from fractions import Fraction

import pytest

from fka.hypergraph import Hypergraph, phi_x_0
from fka.instances import witt_11
from fka.properties import is_read_once
from fka.selfdual import (automorphisms, fano, is_intersecting, is_primitive,
                          is_self_dual, is_transitive, pair_automorphisms)
from fka.substitution import disjoint_vee, disjoint_wedge, single, substitute
from fka.algorithm import eps_exact
from fka.bound import iterate, shortest_edge_ratio, threshold_ratio
from fka.transversal import transversal

X = Hypergraph.from_sets(4, [[0, 1], [2, 3]])                    # (v1^v2)v(v3^v4)
Y = Hypergraph.from_sets(4, [[0, 2], [0, 3], [1, 2], [1, 3]])    # (v1vv2)^(v3vv4)


def word(w):
    """T_w = w_1[w_2[...w_k]] over the atoms X, Y."""
    atoms = {"X": X, "Y": Y}
    MBF = atoms[w[-1]]
    for ch in reversed(w[:-1]):
        MBF = substitute(atoms[ch], [MBF] * 4)
    return MBF


def gk(k):
    f = Hypergraph.from_sets(2, [[0], [1]])
    for _ in range(k - 1):
        f = disjoint_vee(disjoint_wedge(f, f), disjoint_wedge(f, f))
    return f


def tk(k):
    S = single()
    for _ in range(k):
        S = disjoint_wedge(disjoint_vee(S, S), disjoint_vee(S, S))
    return S


# ----------------------------------------------------------------------
# Berge's identities and the substitution law they are the extremes of (4.1)
# ----------------------------------------------------------------------
@pytest.mark.parametrize("A,B", [
    (Hypergraph.from_sets(3, [[0], [1, 2]]), Hypergraph.from_sets(2, [[0, 1]])),
    (Hypergraph.from_sets(4, [[0, 1], [2, 3]]), Hypergraph.from_sets(3, [[0], [1], [2]])),
    (Hypergraph.from_sets(3, [[0, 1], [1, 2], [0, 2]]), Hypergraph.from_sets(2, [[0], [1]])),
])
def test_berge_identities(A, B):
    assert transversal(disjoint_vee(A, B)) == disjoint_wedge(transversal(A), transversal(B))
    assert transversal(disjoint_wedge(A, B)) == disjoint_vee(transversal(A), transversal(B))
    assert len(transversal(disjoint_vee(A, B))) == len(transversal(A)) * len(transversal(B))
    assert len(transversal(disjoint_wedge(A, B))) == len(transversal(A)) + len(transversal(B))


def test_substitution_duality():
    """Tr(T[B_i]) = Tr(T)[Tr(B_i)] -- the one law 4.1's two identities collapse to."""
    T = Hypergraph.from_sets(3, [[0, 1], [1, 2]])
    blocks = [Hypergraph.from_sets(2, [[0], [1]]),
              Hypergraph.from_sets(2, [[0, 1]]),
              Hypergraph.from_sets(3, [[0, 1], [2]])]
    assert transversal(substitute(T, blocks)) == substitute(
        transversal(T), [transversal(b) for b in blocks])


# ----------------------------------------------------------------------
# the GK- and TK-families of 4.2, and their frequency bound in 4.3
# ----------------------------------------------------------------------
@pytest.mark.parametrize("k", [1, 2, 3])
def test_gk_family(k):
    f = gk(k)
    g = transversal(f)
    assert len(f) == 2 ** (2 ** k - 1)
    assert len(g) == 2 ** (2 ** k - 2)
    assert eps_exact(f, g) == Fraction(1, 2 ** (k - 1))


@pytest.mark.parametrize("i", [1, 2, 3])
def test_tk_family(i):
    S = tk(i)
    r = transversal(S)
    assert len(S) == 2 ** (2 * (2 ** i - 1))
    assert len(r) == 2 ** (2 ** i - 1)
    assert eps_exact(S, r) == Fraction(1, 2 ** i)


def test_threshold_ratio_rises_to_two():
    """Observation 2.3.3: both tight families have c -> 2 and s/log2 N -> 1/2."""
    cs, ss = [], []
    for k in (1, 2, 3):
        f = gk(k)
        g = transversal(f)
        cs.append(threshold_ratio(f, g))
        ss.append(shortest_edge_ratio(f, g))
    assert cs == sorted(cs) and cs[-1] < 2
    assert ss == sorted(ss, reverse=True) and ss[-1] > 0.5
    assert math.isclose(cs[-1], 1.8962, abs_tol=1e-3)


def test_iterate_matches_the_explicit_recursion():
    for T in (X, Y):
        D = transversal(T)
        closed = iterate(T, D, 3)
        built = []
        B, Bd = single(), single()
        for _ in range(3):
            B = substitute(T, [B] * T.n)
            Bd = substitute(D, [Bd] * D.n)
            built.append(threshold_ratio(B, Bd))
        assert [round(x, 9) for x in closed] == [round(x, 9) for x in built]


# ----------------------------------------------------------------------
# the word families: GK and TK are the two extreme words
# ----------------------------------------------------------------------
@pytest.mark.parametrize("w", ["X", "Y", "XX", "XY", "YX", "YY"])
def test_word_dual_is_the_swapped_word(w):
    swapped = w.translate(str.maketrans("XY", "YX"))
    assert transversal(word(w)) == word(swapped)


@pytest.mark.parametrize("w", ["X", "Y", "XX", "XY", "YX", "YY"])
def test_word_edge_count_invariant(w):
    """a(w) + a(wbar) = 3(2^k - 1), which forces c_inf >= 3/2."""
    swapped = w.translate(str.maketrans("XY", "YX"))
    a, b = math.log2(len(word(w))), math.log2(len(word(swapped)))
    assert math.isclose(a + b, 3 * (2 ** len(w) - 1))


def test_word_families_beat_gk_and_tk():
    """XX..X is GK and YY..Y is TK at c_inf = 2; the mixed word XY sits at 5/3."""
    for w, want in (("X", 2.0), ("Y", 2.0), ("XY", 5 / 3), ("YX", 5 / 3)):
        T = word(w)
        D = word(w.translate(str.maketrans("XY", "YX")))
        r = T.edge_sizes()[0]
        assert math.isclose(max(math.log2(len(T)), math.log2(len(D))) / (r - 1),
                            want, abs_tol=1e-9)
    assert iterate(word("XY"), word("YX"), 6)[-1] < iterate(X, Y, 6)[-1]


# ----------------------------------------------------------------------
# self-duality, and primitivity as the read-once certificate
# ----------------------------------------------------------------------
def test_fano_is_self_dual_and_primitive():
    F = fano()
    assert is_self_dual(F) and is_intersecting(F)
    perms, capped = automorphisms(F)
    assert not capped and len(perms) == 168
    assert is_transitive(perms, 7) and is_primitive(perms, 7)
    assert not is_read_once(F)


def test_w11_is_a_steiner_system():
    W = witt_11()
    assert W.n == 11 and len(W) == 66
    assert set(W.edge_sizes()) == {5} and set(W.degrees()) == {30}
    cover = {}
    for e in W.to_sets():
        for q in itertools.combinations(e, 4):
            cover[q] = cover.get(q, 0) + 1
    assert len(cover) == 330 and set(cover.values()) == {1}
    assert {(a & b).bit_count() for a, b in itertools.combinations(W.edges, 2)} == {1, 2, 3}


def test_w11_is_self_dual_and_not_read_once():
    W = witt_11()
    assert is_self_dual(W) and is_intersecting(W)
    assert not is_read_once(W)


@pytest.mark.slow
def test_witt_automorphism_group_is_m11():
    """``|Aut| = 7920 = |M_11|``, transitive and primitive, hence not read-once.

    Enumerated in full rather than capped. A capped subgroup proves primitivity
    only once it is known to be *transitive*, and transitivity is not something
    a prefix of the enumeration need have -- a depth-first search returns its
    first solutions sharing a fixed prefix.
    """
    perms, capped = automorphisms(witt_11(), cap=9000)
    assert not capped and len(perms) == 7920
    assert is_transitive(perms, 11)
    assert is_primitive(perms, 11)


@pytest.mark.parametrize("MBF,order", [
    (gk(2), 128),
    (tk(1), 8),
    pytest.param(tk(2), 32768, marks=pytest.mark.slow),
])
def test_read_once_groups_are_iterated_wreath_products(MBF, order):
    """4.4.15: the cotree's group, hence imprimitive -- the contrapositive of
    the primitivity certificate."""
    perms, capped = automorphisms(MBF, cap=order + 1)
    assert not capped and len(perms) == order
    assert is_transitive(perms, MBF.n)
    assert not is_primitive(perms, MBF.n)
    assert is_read_once(MBF)


def test_one_split_collapses_fano_but_not_the_witt_design():
    """Fano's link is a perfect matching -- read-once, imprimitive -- while
    S(4,5,11)'s link is S(3,4,10) and stays primitive."""
    for MBF, prim_after in ((fano(), False), (witt_11(), True)):
        link, _ = phi_x_0(MBF, 0).compact()
        perms, _ = automorphisms(link, cap=2000)
        assert is_transitive(perms, link.n)
        assert is_primitive(perms, link.n) is prim_after
        assert is_read_once(link) is (not prim_after)


def test_pair_automorphisms_agree_with_the_single_hypergraph_group():
    F = fano()
    a, _ = automorphisms(F)
    b, _ = pair_automorphisms(F, F)
    assert sorted(a) == sorted(b)
