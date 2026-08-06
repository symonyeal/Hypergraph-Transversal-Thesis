"""Shared fixtures and helpers."""

from __future__ import annotations

import random

import pytest

from fka.backends import has_sage
from fka.hypergraph import Hypergraph
from fka.irredundant import irredundant


def pytest_collection_modifyitems(config, items):
    """Skip ``@pytest.mark.sage`` tests when SageMath is not importable."""
    if has_sage():
        return
    skip = pytest.mark.skip(reason="SageMath not importable in this interpreter")
    for item in items:
        if "sage" in item.keywords:
            item.add_marker(skip)


def random_sperner(n: int, m: int, rng: random.Random) -> Hypergraph:
    """A random irredundant function on ``n`` variables, aiming for ``m`` terms.

    Built by drawing random non-empty subsets and reducing, so the result is an
    antichain but its size is only approximately ``m``.
    """
    terms = []
    for _ in range(m):
        size = rng.randint(1, max(1, n // 2))
        terms.append(rng.sample(range(n), size))
    return irredundant(Hypergraph.from_sets(n, terms)).reduced


@pytest.fixture
def rng() -> random.Random:
    """Deterministic RNG, so a failure is reproducible from the seed alone."""
    return random.Random(20260804)
