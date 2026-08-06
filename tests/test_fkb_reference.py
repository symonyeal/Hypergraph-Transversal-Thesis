"""Conformance with the MATLAB FK-B reference, and its benchmark families.

The reference implementation cannot be executed here, but it ships its own
recorded expectations: ``Unit_Tests_Equivalency.m`` and
``Unit_tests_Check_Conditions.m`` each state an input and, in a trailing
comment, the answer the authors observed. Those comments are the only direct
evidence of what the original returns, so they are transcribed into
``data/reference/`` and read from there rather than embedded here.

Two things are checked per vector, and they are deliberately different claims:

* the *decision* must match. Equivalent or not is the algorithm's contract and
  admits exactly one right answer, so this is an equality test.
* the reference's own conflicting assignment must be *accepted* by this port.
  It need not be the one produced -- a conflicting assignment is a witness, not
  a canonical value, and several usually exist -- but if the reference's answer
  failed this port's test, one of the two would be wrong about the semantics.

The benchmark families in ``StandardProblems.m`` and the ``*_CNF_DNF.mat``
snapshot are covered by generators rather than by copying the data: see
:func:`fka.instances.self_dual_fano`.
"""

from __future__ import annotations

import pytest

from fka.hypergraph import Hypergraph
from fka.instances import load_data, matching, self_dual_fano, threshold
from fka.transversal import is_dual_oracle, transversal
from fkb.algorithm import check_conditions, fk_b, is_conflict


def M(rows):
    return Hypergraph.from_incidence(rows)


def bits(vector):
    return sum(1 << i for i, x in enumerate(vector) if x)


# ----------------------------------------------------------------------
# Unit_Tests_Equivalency.m
# ----------------------------------------------------------------------
#: Inputs and recorded answers, from data/reference/. ``conflicting_assignment``
#: is null where the reference returned ``[]`` (equivalent).
_EQUIV = load_data("reference", "fkb-equivalency-vectors.json")["vectors"]
EQUIVALENCY_VECTORS = [
    (v["name"], v["cnf"], v["dnf"], v["conflicting_assignment"]) for v in _EQUIV
]


@pytest.mark.parametrize(
    "name,cnf,dnf,reference_ca",
    EQUIVALENCY_VECTORS,
    ids=[v[0] for v in EQUIVALENCY_VECTORS],
)
def test_matches_the_reference_equivalency_vectors(name, cnf, dnf, reference_ca):
    C, D = M(cnf), M(dnf)
    tree = fk_b(C, D, validate=True)

    # The decision has one right answer, and the oracle is asked too so a
    # transcription slip in the vector cannot pass silently.
    assert tree.dual is (reference_ca is None)
    assert tree.dual == is_dual_oracle(C, D)

    if reference_ca is not None:
        assert is_conflict(C, D, bits(reference_ca)), (
            "the assignment the MATLAB recorded is not a conflict under this port"
        )
        assert is_conflict(C, D, tree.verdict.CA)


def test_reference_vectors_5_and_6_are_the_same_pair_transposed():
    """``Unit_Tests_Equivalency.m``'s last two cases swap the roles of C and D.

    Both are recorded as equivalent, which is the consistency check: a CNF and
    DNF representing the same function still do so when the sides are exchanged.

    The vector's CNF is *not* irredundant as written -- one clause contains
    another -- so ``Tr(D)`` equals its reduction rather than the matrix itself.
    FK-B reduces at every node, which is why it decides the pair either way
    round, and it is why duality has to be asserted against the reduced form.
    """
    from fka.irredundant import irredundant

    _, cnf, dnf, _ = EQUIVALENCY_VECTORS[4]
    C, D = M(cnf), M(dnf)
    assert fk_b(C, D).dual is True
    assert fk_b(D, C).dual is True

    reduced = irredundant(C).reduced
    assert reduced.edges != C.edges, "this vector is expected to carry a redundancy"
    assert transversal(D).edges == reduced.edges
    assert transversal(reduced).edges == D.edges


# ----------------------------------------------------------------------
# Unit_tests_Check_Conditions.m
# ----------------------------------------------------------------------
#: ``Unit_tests_Check_Conditions.m`` is a scratch script: it calls the function
#: and records what the author saw in a trailing comment, with no assertion. For
#: these the comment reads as a 1-indexed variable list and that list is a
#: genuine conflicting assignment; the file's last two carry a bare "% 1" that is
#: not one in either direction, and are exercised without a reference witness.
_COND = load_data("reference", "fkb-condition-vectors.json")
CONDITION_VECTORS = [
    (v["cnf"], v["dnf"], v["conflicting_assignment"]) for v in _COND["vectors"]
]
CONDITION_VECTORS_WITHOUT_A_USABLE_WITNESS = [
    (v["cnf"], v["dnf"]) for v in _COND["vectors_without_a_usable_witness"]
]


@pytest.mark.parametrize("cnf,dnf,reference_ca", CONDITION_VECTORS)
def test_matches_the_reference_check_conditions_vectors(cnf, dnf, reference_ca):
    """Every recorded vector fails a condition, and the recorded witness holds.

    Which condition fired is not asserted: the MATLAB records only the
    assignment, and this port checks the degenerate cases before the three
    numbered ones, so the attribution can legitimately differ.
    """
    C, D = M(cnf), M(dnf)
    verdict = check_conditions(C, D)
    assert verdict is not None, "the reference recorded a failure here"
    assert verdict.dual is False
    assert verdict.reason.startswith(("cond_", "constant"))

    S = bits([1 if (i + 1) in reference_ca else 0 for i in range(max(C.n, D.n))])
    assert is_conflict(C, D, S), (
        "the assignment the MATLAB recorded is not a conflict under this port"
    )
    assert is_conflict(C, D, verdict.CA)
    assert fk_b(C, D, validate=True).dual is False


@pytest.mark.parametrize("cnf,dnf", CONDITION_VECTORS_WITHOUT_A_USABLE_WITNESS)
def test_reference_vectors_whose_recorded_answer_does_not_parse(cnf, dnf):
    """Still a regression: this port must decide them and produce a real witness."""
    C, D = M(cnf), M(dnf)
    verdict = check_conditions(C, D)
    assert verdict is not None and verdict.dual is False
    assert is_conflict(C, D, verdict.CA)
    assert fk_b(C, D, validate=True).dual is is_dual_oracle(C, D) is False
    # The bare "% 1" the file records is not a conflict in either direction.
    assert not is_conflict(C, D, bits([1] + [0] * (max(C.n, D.n) - 1)))


# ----------------------------------------------------------------------
# StandardProblems.m and the .mat benchmark shapes
# ----------------------------------------------------------------------
def test_matching_family():
    """``M(v)``: v/2 disjoint pairs, so ``|Tr| = 2^(v/2)``."""
    for v in (2, 4, 6, 8):
        MBF = matching(v)
        assert len(MBF) == v // 2
        assert MBF.edge_sizes() == [2] * (v // 2)
        assert len(transversal(MBF)) == 2 ** (v // 2)
        assert fk_b(transversal(MBF), MBF, validate=True).dual is True
    for bad in (0, 3):
        with pytest.raises(ValueError, match="even and at least 2"):
            matching(bad)


def test_threshold_family():
    for v in (2, 4, 6):
        MBF = threshold(v)
        assert all(size == 2 for size in MBF.edge_sizes())
        assert fk_b(transversal(MBF), MBF, validate=True).dual is True
    with pytest.raises(ValueError, match="even and at least 2"):
        threshold(5)


@pytest.mark.parametrize(
    "k,n,edges,sizes",
    [
        # Shapes read off SDFP16_CNF_DNF.mat and SDFP23_CNF_DNF.mat: one pair,
        # 7k lines with b, and 7^k line-tuples with a.
        (1, 9, 15, {2: 1, 4: 14}),
        (2, 16, 64, {2: 1, 4: 14, 7: 49}),
        (3, 23, 365, {2: 1, 4: 21, 10: 343}),
    ],
)
def test_self_dual_fano_matches_the_benchmark_shape(k, n, edges, sizes):
    import collections

    E = self_dual_fano(k)
    assert E.n == n
    assert len(E) == edges == 7**k + 7 * k + 1
    assert dict(collections.Counter(E.edge_sizes())) == sizes


@pytest.mark.parametrize("k", [1, 2])
def test_self_dual_fano_is_self_transversal(k):
    """The property that makes it a benchmark, checked against Berge, not FK.

    ``k = 3`` is excluded here only for runtime; ``python -m fkb benchmark``
    covers it.
    """
    E = self_dual_fano(k)
    assert transversal(E).edges == E.edges
    assert fk_b(E, E, validate=True).dual is True


def test_self_dual_fano_is_not_the_thesis_sdfp():
    """``sdfp`` is the k disjoint copies; ``self_dual_fano`` self-dualises them."""
    from fka.instances import sdfp

    assert sdfp(2).n == 14 and len(sdfp(2)) == 14
    assert self_dual_fano(2).n == 16 and len(self_dual_fano(2)) == 64
    with pytest.raises(ValueError, match="at least 1"):
        self_dual_fano(0)
