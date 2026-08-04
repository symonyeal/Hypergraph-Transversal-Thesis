# Hypergraphs and FK-A

[![Tests](https://github.com/symonyeal/Hypergraph-Transversal-Thesis/actions/workflows/tests.yml/badge.svg)](https://github.com/symonyeal/Hypergraph-Transversal-Thesis/actions/workflows/tests.yml)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/symonyeal/Hypergraph-Transversal-Thesis/main?urlpath=lab)

This is the maintained research implementation for Saimon Y. Islam's work on
hypergraph dualization, symmetry, and the Fredman-Khachiyan Algorithm A (FK-A).
The mathematical reference is [Saimon Islam Manuscript MSc Thesis (1).pdf](Saimon%20Islam%20Manuscript%20MSc%20Thesis%20(1).pdf).
The definitive repository is
[symonyeal/Hypergraph-Transversal-Thesis](https://github.com/symonyeal/Hypergraph-Transversal-Thesis).

The project replaces the old notebooks, free-form text dumps, global state, and
per-frame Graphviz images with:

- an exact bitset hypergraph model;
- faithful and thesis-modified FK-A variants;
- an independent transversal oracle and checkable non-duality certificates;
- structured, provenance-bearing JSON instances;
- deterministic JSON, DOT, and self-contained interactive HTML results;
- a pytest suite tied to the published recursion trees; and
- an optional SageMath backend for independent group and graph checks.

## Start here

From this directory, the source-tree tests require no installation:

```text
C:\Python314\python.exe -m pytest
```

Install the command-line package into the machine's existing Python (no project
virtual environment is required):

```text
C:\Python314\python.exe -m pip install -e ".[dev]"
```

Then use the same interpreter consistently:

```text
C:\Python314\python.exe -m fka env
C:\Python314\python.exe -m fka list
C:\Python314\python.exe -m fka show fano
C:\Python314\python.exe -m fka run fano
C:\Python314\python.exe -m fka verify
```

`run` writes a portable HTML explorer and machine-readable JSON to `results/`.
Open the HTML directly; it has no network, font, JavaScript, or image dependency.

## Research contracts

- `variant="faithful"` implements the FK96 base case. `variant="modified"`
  adds the thesis Section 5.1.1 `D_{1,k}` / `D_{k,1}` rule.
- `pivot_rule="degree_sum"` reproduces the legacy implementation and published
  tree numbering. `pivot_rule="frequency"` implements the normalized frequency
  quantity bounded in the thesis.
- Changes to preprocessing, pivot choice, base cases, or recursion order must
  preserve the published node counts and split sequences in
  `tests/test_thesis_reproduction.py` unless the research interpretation itself
  is deliberately revised.
- The exact oracle is independent of FK-A. A passing FK-A run alone is not a
  ground-truth proof.
- Corrections never overwrite a printed instance. `6a-verbatim` and
  `8ver-verbatim` preserve the thesis data; `6a` and `8ver` record the corrected
  pairs and the evidence in their JSON notes.

## Filing system

| Path | Purpose |
| --- | --- |
| `src/fka/` | Maintained library and CLI |
| `data/instances/` | Canonical research inputs and expected baselines |
| `tests/` | Unit, randomized, thesis-reproduction, output, and Sage checks |
| `results/` | Regenerable JSON, HTML, and DOT artifacts |
| `experiments/` | Explicit migration or one-off research programs |
| `docs/` | Architecture, Sage setup, and assessment notes |
| `Reading/` | Local-only third-party reference library; PDFs are not redistributed |
| `FK- B/` | Local-only third-party MATLAB FK-B reference snapshot |
| `_archive/` | Preserved superseded material with an archive note |

Do not add new research logic to notebooks or `.txt` logs. Add a module or an
explicit program, store inputs as JSON, and make the claim executable as a test.
See [docs/architecture.md](docs/architecture.md) for the repeatable workflow.

## SageMath

The package runs completely under CPython. Sage is an optional independent
backend, not an import-time dependency. Use the Binder badge above to run Sage
10.9 in a browser without installing it locally. GitHub Actions runs the same
Sage-backed checks for every definitive revision. On this Windows machine,
local Sage would require WSL 2; see [docs/sagemath.md](docs/sagemath.md).

## Authorship and license

The legacy implementation credits Saimon Islam and Lacey Liang and declares
GPL-3.0 licensing. The archived original notice is retained in
`_archive/20260804-legacy-code/LEGACY README.md`. The maintained package carries
the same `GPL-3.0-or-later` project metadata.
