"""``fkb`` -- the Fredman-Khachiyan algorithm B, on the FK-A project's model.

Companion to :mod:`fka`. FK-B decides the same question on the same instance
library and writes the same reports; what differs is the algorithm. See
:mod:`fkb.algorithm` for the port and how it departs from the MATLAB reference
in ``FK- B/FK-master/src``.

The model, tree, instance library, analysis passes and report writers live in
:mod:`fka` and are shared: they are specific to neither algorithm, and one point
of running both is that the trees are the same type and so directly comparable.

::

    from fka import load
    from fkb import fk_b

    inst = load("fano")
    print(fk_b(inst.cnf, inst.dnf).summary())

``python -m fkb compare --all`` runs FK-A, FK-B and the brute-force oracle on
every instance and tabulates the three answers.
"""

from __future__ import annotations

__version__ = "2.0.0"

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
