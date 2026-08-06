"""The FK-A frequency bound as a measured quantity, and its limit on a family.

Thesis 2.3 bounds the most frequent variable of a dual pair from below,
``eps >= 1/log2 N`` with ``N = nC + nD`` (Lemma 2.3.2, [FK96]), and
Observation 2.3.3 says that is tight to a factor of 2. The single number that
places a pair against both is

    c(C, D) = eps(C, D) * log2 N,          1 <= c,  c ~ 2 for the tight families.

Small ``c`` means FK-A's pivot guarantee is at its weakest, which is what makes
the alternating trees of Chapter 4 hard. Section 4.3 computes ``c -> 2`` for the
GK- and TK-families; :func:`iterate` reproduces that from the template alone.

Note. Thesis Corollary 2.3.1 states ``s >= log n`` for the shortest term,
but its own derivation ``n * 2^-s >= 1`` gives ``s <= log2 n``, which is the
[FK96] statement and the direction used here.
"""

from __future__ import annotations

import math

from .algorithm import eps_exact
from .hypergraph import Hypergraph, verts

__all__ = [
    "threshold_ratio",
    "shortest_edge_ratio",
    "log2_int",
    "iterate",
]


def threshold_ratio(cnf: Hypergraph, dnf: Hypergraph) -> float:
    """``c = eps * log2 N``. See the module docstring."""
    return float(eps_exact(cnf, dnf)) * math.log2(len(cnf) + len(dnf))


def shortest_edge_ratio(cnf: Hypergraph, dnf: Hypergraph) -> float:
    """``s / log2 N``, the other half of Observation 2.3.3 (tight families: 1/2)."""
    s = min(min(cnf.edge_sizes()), min(dnf.edge_sizes()))
    return s / math.log2(len(cnf) + len(dnf))


def log2_int(x: int) -> float:
    """``log2`` of an arbitrarily large integer.

    Edge counts here are doubly exponential -- ``|E(f_k)| = 2^(2^k - 1)`` -- so
    ``float(x)`` overflows within a handful of levels.
    """
    bl = x.bit_length()
    if bl <= 53:
        return math.log2(x)
    sh = bl - 53
    return sh + math.log2(x >> sh)


def iterate(T: Hypergraph, D: Hypergraph, levels: int = 8) -> list[float]:
    """``c_j`` for the family ``Gamma_0 = {{v}}``, ``Gamma_{j+1} = T[Gamma_j,...]``.

    ``D`` is ``Tr(T)``. With every block equal to the previous level ``B``,

        m'          = sum_{e in E(T)} m^|e|
        deg'(v/B_i) = deg_B(v) * sum_{e ni i} m^(|e|-1)

    and the same on the dual side, by the substitution law of
    :mod:`fka.substitution`. Exact integer arithmetic throughout; only the
    reported ratio is a float. No hypergraph is built, so levels far past what
    is representable are still cheap.

    ``c_j`` converges iff ``T`` and ``Tr(T)`` are ``r``-uniform and regular on
    ``t = r^2`` vertices, and then the limit is
    ``max(log2|E(T)|, log2|E(Tr T)|) / (r - 1)``.
    """
    m = md = 1
    leps = lepsd = 0.0          # frequencies in log2, so 2^-hundreds is fine
    out = []
    for _ in range(levels):
        m_new = sum(m ** len(verts(e)) for e in T.edges)
        md_new = sum(md ** len(verts(e)) for e in D.edges)
        fac = max(sum(m ** (len(verts(e)) - 1) for e in T.edges if e >> i & 1)
                  for i in range(T.n))
        facd = max(sum(md ** (len(verts(e)) - 1) for e in D.edges if e >> i & 1)
                   for i in range(D.n))
        leps += log2_int(m) + log2_int(fac) - log2_int(m_new)
        lepsd += log2_int(md) + log2_int(facd) - log2_int(md_new)
        m, md = m_new, md_new
        out.append(2 ** max(leps, lepsd) * log2_int(m + md))
    return out
