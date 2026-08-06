# SageMath compatibility

The package runs completely under CPython. Sage is an optional independent
backend, never an import-time dependency: nothing imports `sage.all` at module
load, `--backend sage` fails loudly when Sage is absent, and the default selects
Sage only when `sage.all` is genuinely importable.

## Current status

As assessed on 2026-08-05, on the Linux workstation this repository is now
maintained from:

- Sage 10.9 is installed and `sage.all` is importable, so the Sage-marked tests
  **execute** rather than skip;
- the complete suite passes under both interpreters — 361 passed / 13 skipped
  under plain CPython, 373 passed / 1 skipped under Sage's own Python;
- `experiments/hard_case_search.py`, which needs Sage's transitive-group
  database, runs degrees 3–15 in under two minutes.

The superseded assessment — Windows CPython at `C:\Python314`, no `sage`
executable, no WSL, every Sage test skipping — is preserved in
[`_archive/20260805-stale-platform-assumptions/`](../_archive/20260805-stale-platform-assumptions/)
together with the Windows-specific instructions that went with it.

## Installing

Sage's own installation guide recommends conda-forge on Linux:

```text
mamba create -n sage -c conda-forge sage python=3.12
```

Then install this package into **Sage's** interpreter, which is simply that
environment's Python:

```text
<env>/bin/python -m pip install -e ".[dev]"
```

Two things about the conda-forge build are worth knowing, because both cost
time to rediscover:

- **There is no `sage -python` subcommand.** The distribution's own docs use it
  and it works for a source or binary install, but under conda Sage is a normal
  Python package inside the environment. Use that environment's `python`
  directly; `import sage.all` works from it.
- **The GAP `design` package is not included**, so `designs.WittDesign(...)`
  raises `RuntimeError: Error loading GAP package design`. Build the designs
  this project needs from codes instead — `fka.instances.witt_11` derives
  `S(4,5,11)` from the ternary Golay code, and needs no Sage at all.

## Browser and continuous-integration execution

- The Binder badge in `README.md` builds `Dockerfile` from the official
  `ghcr.io/sagemath/sage-binder-env:10.9` image and opens JupyterLab in a
  temporary browser session. Binder is ephemeral: download anything generated
  before ending the session.
- `.github/workflows/tests.yml` runs the complete suite and independent oracle
  verification in the same Sage 10.9 image on every push and pull request. Do
  not claim Sage verification for a revision unless its Sage job passes.

## Verification commands

Substitute the Sage environment's interpreter for `python`:

```text
python -m pytest                       the whole suite, Sage tests included
python -m fka env                      what this interpreter can do
python -m fka verify                   FK-A against the oracle
python -m fka run fano --backend sage  one instance, Sage-backed annotation
python experiments/hard_case_search.py 3 15
```

The expected signal is that the tests marked `sage` execute rather than skip.
They cross-check automorphism orders and graph classes against the pure-Python
implementations, and assert that `Aut(w11)` is `M11`.

## Supported Sage interfaces

The adapter uses public Sage interfaces documented in the current reference:

- `IncidenceStructure(...).automorphism_group()` for the independent hypergraph
  automorphism calculation, plus `.canonical_label()`, `.is_t_design()` and
  `.is_isomorphic()` in the hard-case search;
- GAP-backed `structure_description()` for finite-group names; and
- `Graph.is_cograph()`, `is_chordal()`, `is_chordal_bipartite()`,
  `is_perfect()`, and the other public graph predicates.

The adapter preserves Sage's 0-based incidence points and maps them back to the
thesis' 1-based vertex labels. A pure-Python regression test protects that
translation even when Sage is not installed.

## What the two backends do and do not agree on

Group **order, generators and orbits** agree exactly, and the suite asserts it
for every instance in the library, on both sides.

Group **names** need not agree, and that is not a defect in either backend.
GAP's `StructureDescription` is not canonical for 2-groups: handed the search's
generators for `f2g2` it returns `(((C4 x C4) : C2) : C2) : C2`, and handed an
abstract `D4 wr C2` it returns `(((C2 x C2 x C2 x C2) : C2) : C2) : C2` — the
thesis' own second spelling of the same group. `IsomorphismGroups` confirms the
two are isomorphic.

So `tests/test_symmetry.py` pins the *name* assertion to `backend="python"`,
where the name comes from this package's catalogue and is a stable contract, and
checks the Sage backend against what it actually owes: a group of the right
isomorphism type. Before 2026-08-05 that test asserted the catalogue name
against whichever backend happened to be live, so it passed on CPython and
failed under Sage — including in the Sage CI job.
