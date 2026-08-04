"""``fka`` -- FK-A, and the model both algorithms share.

Research code for *Tracing the Effects of Symmetry in Hypergraphs On the
Fredman-Khachiyan Algorithm* (Saimon Y. Islam, MSc thesis, Simon Fraser
University, 2024) and its continuation. FK-B is :mod:`fkb`, which builds on this
package; the split is layout, not dependence.

``hypergraph``   the bitset-backed :class:`Hypergraph`
``sperner``      Sperner reduction -- FK-A Step 1, FK-B's ``Irredundant``
``algorithm``    FK-A, faithful and modified
``tree``         the recursion tree, serialisable, algorithm-neutral
``transversal``  ``Tr(H)`` by Berge's method -- the independent oracle
``automorphism`` ``Aut(H)``
``groups``       permutation groups and small-group naming without SageMath
``graphs``       primal graphs and graph classes
``properties``   conformality, alpha-acyclicity, read-once
``analysis``     annotating a tree into the thesis' automorphism tree
``report``       self-contained interactive HTML
``dot``          Graphviz export
``instances``    the library in ``data/instances``, and its generators
``backends``     SageMath when present, pure Python otherwise
``cli``          ``python -m fka``

::

    from fka import Hypergraph, fk_a, load

    inst = load("fano")
    print(fk_a(inst.G, inst.H).summary())

Nothing imports SageMath at module load, so ``import fka`` works on plain
CPython; see :mod:`fka.backends`.
"""

from __future__ import annotations

__version__ = "1.1.0"

from .algorithm import eps, fk_a, is_dual
from .hypergraph import Hypergraph
from .instances import Instance, list_ids, load, load_all, load_archived
from .transversal import is_dual_oracle, transversal
from .tree import RecursionNode, RecursionTree, Verdict

__all__ = [
    "__version__",
    "Hypergraph",
    "Instance",
    "RecursionNode",
    "RecursionTree",
    "Verdict",
    "eps",
    "fk_a",
    "is_dual",
    "is_dual_oracle",
    "list_ids",
    "load",
    "load_all",
    "load_archived",
    "transversal",
]
