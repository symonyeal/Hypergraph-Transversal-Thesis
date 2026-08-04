"""Entry point for ``python -m fkb`` (and ``sage -python -m fkb``)."""

from __future__ import annotations

from .cli import main

if __name__ == "__main__":
    raise SystemExit(main())
