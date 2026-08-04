"""``fkb`` -- the Fredman-Khachiyan algorithm B, on the FK-A project's model.

Companion to :mod:`fka`. FK-B decides the same question on the same instance
library and writes the same reports; what differs is the algorithm. See
:mod:`fkb.algorithm` for the port and for how it departs from the MATLAB
reference in ``FK- B/FK-master/src``.

The hypergraph model, recursion tree, instance library, analysis passes, and
report writers live in :mod:`fka` and are shared, not duplicated: they are not
specific to either algorithm, and one of the points of running both is that the
trees are directly comparable because they are the same type.

Quick start
-----------
::

    from fka import load
    from fkb import fk_b

    inst = load("fano")
    tree = fk_b(inst.G, inst.H)
    print(tree.summary())

``python -m fkb compare --all`` runs FK-A and FK-B against the brute-force
oracle on every instance and tabulates the three answers.
"""

from __future__ import annotations

__version__ = "1.0.0"

from .algorithm import (
    SPLIT_RULES,
    VARIANTS,
    check_conditions,
    choose_split_var,
    easy_case,
    fk_b,
    is_conflict,
    is_dual,
    mu,
    mu_frequent,
)

__all__ = [
    "__version__",
    "SPLIT_RULES",
    "VARIANTS",
    "check_conditions",
    "choose_split_var",
    "easy_case",
    "fk_b",
    "is_conflict",
    "is_dual",
    "mu",
    "mu_frequent",
]
