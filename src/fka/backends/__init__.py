"""Backend selection: SageMath when available, plain Python otherwise.

The research code originally ran only inside a SageMath 10.3 Jupyter kernel, so
none of it could be executed, tested or version-controlled on a machine without
Sage -- which on Windows means WSL or a conda environment. Everything in this
package therefore has a pure-Python implementation that is exact and needs only
``networkx``, and Sage is an *optional accelerator and cross-check*
rather than a hard requirement.

What Sage adds when it is present:

* ``StructureDescription`` from GAP, which names any finite group, where the
  pure-Python catalogue in :mod:`fka.groups` covers a fixed list and otherwise
  reports ``order N (unidentified)``;
* ``IncidenceStructure.automorphism_group`` as an independent check on the
  built-in automorphism search;
* polynomial-time graph-class recognisers, where :mod:`fka.graphs` enumerates
  induced subgraphs and is capped at 20 vertices.

Running under Sage
------------------
No code change is needed -- the same package, the same entry points::

    sage -python -m fka run f2g2        # Sage's Python
    python -m fka run f2g2              # plain CPython

``fka.backends.describe()`` reports which backend is live, and every result
object records the backend that produced it so archived output stays
attributable.
"""

from __future__ import annotations

import importlib.util
from functools import lru_cache
from typing import Optional

__all__ = ["has_sage", "resolve_backend", "describe", "BACKENDS"]

BACKENDS = ("python", "sage")


@lru_cache(maxsize=1)
def has_sage() -> bool:
    """True iff ``sage.all`` can be imported in this interpreter.

    Checked via :func:`importlib.util.find_spec` rather than a bare import:
    importing ``sage.all`` costs several seconds, and this is called on paths
    that may not need Sage at all.
    """
    try:
        return importlib.util.find_spec("sage.all") is not None
    except (ImportError, ValueError):
        return False


def resolve_backend(name: Optional[str] = None) -> str:
    """Normalise a backend request to ``"python"`` or ``"sage"``.

    ``None`` means "Sage if available". Asking for ``"sage"`` when it is not
    importable is an error rather than a silent downgrade: a run that was meant
    to cross-check against Sage must not quietly report Python results.
    """
    if name is None:
        return "sage" if has_sage() else "python"
    if name not in BACKENDS:
        raise ValueError(f"unknown backend {name!r}; expected one of {BACKENDS}")
    if name == "sage" and not has_sage():
        raise RuntimeError(
            "backend 'sage' requested but sage.all is not importable in this "
            "interpreter. Run under `sage -python`, or omit --backend to use "
            "the pure-Python implementation."
        )
    return name


def describe() -> dict[str, object]:
    """A short report on what this interpreter can do."""
    import sys

    info: dict[str, object] = {
        "python": sys.version.split()[0],
        "executable": sys.executable,
        "sage_available": has_sage(),
        "active_backend": resolve_backend(None),
    }
    if has_sage():
        try:
            from sage.version import version as sage_version

            info["sage_version"] = sage_version
        except Exception:  # pragma: no cover - depends on the Sage build
            info["sage_version"] = "unknown"
    for mod in ("numpy", "networkx", "matplotlib"):
        spec = importlib.util.find_spec(mod)
        info[f"{mod}_available"] = spec is not None
    return info
