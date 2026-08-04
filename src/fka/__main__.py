"""Entry point for ``python -m fka`` (and ``sage -python -m fka``)."""

from __future__ import annotations

from .cli import main

if __name__ == "__main__":
    raise SystemExit(main())
