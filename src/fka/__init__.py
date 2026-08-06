"""``fka`` -- FK-A, and the model both algorithms share.

Research code for *Tracing the Effects of Symmetry in Hypergraphs On the
Fredman-Khachiyan Algorithm* (Saimon Y. Islam, MSc thesis, Simon Fraser
University, 2024) and its continuation. FK-B is :mod:`fkb`, which builds on this
package; the split is layout, not dependence.

The code is written in the monotone-Boolean register of the FK-B reference;
``docs/dictionary.md`` maps it to the thesis' hypergraph notation.

``hypergraph``   the bitset-backed :class:`Hypergraph` and the restrictions
``irredundant``  ``Irredundant`` -- FK-A Step 1
``algorithm``    FK-A, faithful and modified
``tree``         the recursion tree, serialisable, algorithm-neutral
``transversal``  ``Tr`` by Berge's method -- the independent oracle
``automorphism`` ``Aut``
``groups``       permutation groups and small-group naming without SageMath
``graphs``       primal graphs and graph classes
``properties``   conformality, alpha-acyclicity, read-once
``selfdual``     self-transversal families, primitivity, the Witt design
``substitution`` monotone composition and Berge's identities
``bound``        the Chapter 2 frequency bound as a measured quantity
``analysis``     annotating a tree into the thesis' automorphism tree
``report``       self-contained interactive HTML
``dot``          Graphviz export
``instances``    the library in ``data/instances``, and its generators
``backends``     SageMath when present, pure Python otherwise
``cli``          ``python -m fka``

::

    from fka import fk_a, load

    inst = load("fano")
    print(fk_a(inst.cnf, inst.dnf).summary())

Nothing imports SageMath at module load, so ``import fka`` works on plain
CPython; see :mod:`fka.backends`.
"""

from __future__ import annotations

__version__ = "2.0.0"

from .algorithm import eps, fk_a, is_dual
from .hypergraph import Hypergraph
from .instances import Instance, list_ids, load, load_all, load_archived
from .transversal import is_dual_oracle, transversal
from .tree import Node, Tree, Verdict

__all__ = [
    "__version__",
    "Hypergraph",
    "Instance",
    "Node",
    "Tree",
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
