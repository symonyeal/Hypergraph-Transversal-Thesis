"""Shared fixtures and helpers."""

from __future__ import annotations

import random

import pytest

from fka.backends import has_sage
from fka.hypergraph import Hypergraph
from fka.sperner import sperner_reduce


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
