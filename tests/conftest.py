"""Shared fixtures and helpers."""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from fka.backends import has_sage
from fka.hypergraph import Hypergraph
from fka.sperner import sperner_reduce


def _usable(path: Path) -> bool:
    """Whether pytest can own ``path``: create it, list it, and write inside it."""
    try:
        path.mkdir(parents=True, exist_ok=True)
        list(path.iterdir())
        probe = path / ".probe"
        probe.write_text("", encoding="utf-8")
        probe.unlink()
    except OSError:
        return False
    return True


def pytest_configure(config):
    """Keep temporary test artifacts in the shared work folder on this PC.

    The project contract forbids session-scoped ``%TEMP%`` files.  On Binder,
    Linux, and CI the normal pytest temporary directory remains appropriate and
    avoids embedding a machine-specific Windows path in ``pyproject.toml``.

    Each candidate is probed rather than assumed. A basetemp pytest cannot list
    or delete makes every ``tmp_path`` test error out before it runs, and one
    wedged that way is not recoverable from in-process, so an unusable candidate
    is passed over instead of taking the whole suite down with it.
    """
    if config.option.basetemp is not None:
        return
    project = Path(__file__).resolve().parents[1]
    shared = project.parents[2] / "Claude Func Folder"
    if not shared.is_dir():
        return
    for name in ("pytest-fka", "pytest-fk"):
        if _usable(shared / name):
            config.option.basetemp = str(shared / name)
            return


def pytest_collection_modifyitems(config, items):
    """Skip ``@pytest.mark.sage`` tests when SageMath is not importable."""
    if has_sage():
        return
    skip = pytest.mark.skip(reason="SageMath not importable in this interpreter")
    for item in items:
        if "sage" in item.keywords:
            item.add_marker(skip)


def random_sperner(n: int, m: int, rng: random.Random) -> Hypergraph:
    """A random Sperner hypergraph on ``n`` vertices, aiming for ``m`` edges.

    Built by drawing random non-empty subsets and Sperner-reducing, so the
    result is an antichain but its size is only approximately ``m``.
    """
    edges = []
    for _ in range(m):
        size = rng.randint(1, max(1, n // 2))
        edges.append(rng.sample(range(n), size))
    H = Hypergraph.from_sets(n, edges)
    return sperner_reduce(H).reduced


@pytest.fixture
def rng() -> random.Random:
    """Deterministic RNG, so a failure is reproducible from the seed alone."""
    return random.Random(20260804)
