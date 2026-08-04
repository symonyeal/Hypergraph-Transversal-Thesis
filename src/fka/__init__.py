"""``fka`` -- the Fredman-Khachiyan algorithm A, and hypergraph symmetry under it.

Research code for *Tracing the Effects of Symmetry in Hypergraphs On the
Fredman-Khachiyan Algorithm* (Saimon Y. Islam, MSc thesis, Simon Fraser
University, 2024) and its continuation.

Layout
------
``hypergraph``   the bitset-backed :class:`Hypergraph` type
``sperner``      Sperner reduction (FK-A Step 1)
``algorithm``    FK-A itself, faithful and modified variants
``tree``         the recursion tree, serialisable
``transversal``  brute-force ``Tr(H)``, the independent oracle
``automorphism`` ``Aut(H)``
``groups``       permutation groups and small-group naming without SageMath
``graphs``       primal graphs and graph classes
``properties``   conformality, alpha-acyclicity, read-once
``analysis``     annotating a recursion tree into the thesis' automorphism tree
``report``       self-contained interactive HTML output
``dot``          Graphviz export
``instances``    the problem-instance library in ``data/instances``
``backends``     SageMath when present, pure Python otherwise
``cli``          ``python -m fka``

Quick start
-----------
::

    from fka import Hypergraph, fk_a, load

    inst = load("fano")
    tree = fk_a(inst.G, inst.H)
    print(tree.summary())

Nothing here imports SageMath at module load, so ``import fka`` works on a
plain CPython install; see :mod:`fka.backends`.
"""

from __future__ import annotations

__version__ = "1.0.0"

from .algorithm import epsilon, fk_a, is_dual
from .hypergraph import Hypergraph
from .instances import Instance, list_ids, load, load_all
from .transversal import is_dual_oracle, transversal
from .tree import RecursionNode, RecursionTree, Verdict

__all__ = [
    "__version__",
    "Hypergraph",
    "Instance",
    "RecursionNode",
    "RecursionTree",
    "Verdict",
    "epsilon",
    "fk_a",
    "is_dual",
    "is_dual_oracle",
    "list_ids",
    "load",
    "load_all",
    "transversal",
]
