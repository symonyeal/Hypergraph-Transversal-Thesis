# SageMath compatibility

## Current status on this machine

As assessed on 2026-08-04:

- Windows CPython is `C:\Python314\python.exe`;
- the complete pure-Python suite passes;
- neither a `sage` executable nor WSL is installed; and
- Sage-only tests therefore skip rather than silently using the Python backend.

The maintained code does not import Sage at module load. `--backend sage`
fails clearly when Sage is absent, while the default selects Sage only when
`sage.all` is genuinely importable.

## Browser and continuous-integration execution

The definitive GitHub repository has two no-install Sage paths:

- The Binder badge in `README.md` builds `Dockerfile` from the official
  `ghcr.io/sagemath/sage-binder-env:10.9` image and opens JupyterLab in a
  temporary browser session.
- `.github/workflows/tests.yml` runs the complete suite and independent oracle
  verification in the same Sage 10.9 image for pushes and pull requests.

Binder is for interactive, ephemeral work. Download any generated files before
ending the session. GitHub Actions is the durable cross-backend record: do not
claim Sage verification for a revision unless its Sage job passes.

## Windows installation

SageMath's official 10.8 installation guide recommends WSL 2 plus Ubuntu on
Windows. Installing WSL is an administrator-level operating-system change and
can require a reboot, so it is intentionally not performed implicitly by this
project.

1. In an elevated Windows terminal, install WSL and Ubuntu following the
   official guide: <https://doc.sagemath.org/html/en/installation/index.html#windows>.
2. Complete the first-launch Ubuntu username and password setup.
3. In Ubuntu, follow the guide's current Miniforge commands and create the Sage
   distribution environment.
4. Activate Sage, change to this project through its WSL path, and install this
   package into Sage's own interpreter:

```text
cd "/mnt/c/Users/sislam/Downloads/New Work Folder (Former Codex)/My Projects/F U N/Hypergraphs and FK-A"
sage -python -m pip install -e ".[dev]"
```

Sage necessarily supplies its own supported Python runtime; this is not a
second project virtual environment and should not be used for the ordinary
Windows CPython workflow.

## Verification commands

Run all tests, including the Sage marker:

```text
sage -python -m pytest
sage -python -m fka env
sage -python -m fka verify
sage -python -m fka run fano --backend sage
```

The expected signal is that the eight tests marked `sage` execute rather than
skip. They compare automorphism orders for `fano`, `f2g2`, `6b`, and `6d`, and
graph classes for `fano`, `6b`, `6c`, and `6d`.

## Supported Sage interfaces

The adapter uses public Sage interfaces documented in the current reference:

- `IncidenceStructure(...).automorphism_group()` for the independent
  hypergraph automorphism calculation;
- GAP-backed `structure_description()` for finite-group names; and
- `Graph.is_cograph()`, `is_chordal()`, `is_chordal_bipartite()`,
  `is_perfect()`, and the other public graph predicates.

The adapter preserves Sage's 0-based incidence points and maps them back to the
thesis' 1-based vertex labels. A pure-Python regression test protects that
translation even when Sage is not installed.
