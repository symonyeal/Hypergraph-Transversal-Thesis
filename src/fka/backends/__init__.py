"""Backend selection: SageMath when available, plain Python otherwise.

Everything here has an exact pure-Python implementation needing only
``networkx``; Sage is an optional accelerator and cross-check, never a
requirement. On Windows a Sage install means WSL or conda, so requiring it would
put the whole project out of reach of the machine it is maintained on.

Sage adds GAP's ``StructureDescription``, which names any finite group where
:mod:`fka.groups`' catalogue covers a fixed list; ``IncidenceStructure
.automorphism_group`` as an independent check on the built-in search; and
polynomial-time graph-class recognisers where :mod:`fka.graphs` enumerates
induced subgraphs under a variable cap.

No code change is needed to run under it -- same package, same entry points::

    sage -python -m fka run f2g2        # Sage's Python
    python -m fkb run f2g2              # plain CPython

:func:`describe` reports which backend is live, and every result object records
the one that produced it, so archived output stays attributable.
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
