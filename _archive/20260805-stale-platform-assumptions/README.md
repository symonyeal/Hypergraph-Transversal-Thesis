# Withdrawn 2026-08-05 — the Windows-only platform assumptions

## What was here

The repository was written on a Windows machine with no SageMath, and three
places recorded that machine rather than the project. Maintenance moved to a
Linux workstation with Sage 10.9 installed, at which point all three were
describing a machine that no longer exists.

| File here | Was | Now |
| --- | --- | --- |
| `sagemath.md` | `docs/sagemath.md` as of 2026-08-04 | rewritten for the Sage-present state |
| `conftest-before.py` | `tests/conftest.py` as of 2026-08-04 | `pytest_configure` and `_usable` removed |
| `architecture-excerpt-before.md` | `docs/architecture.md` as of 2026-08-04 | interpreter paths made platform-neutral |

## 1. The dated Sage assessment

`docs/sagemath.md` opened with a status block asserting that Windows CPython is
`C:\Python314\python.exe`, that "neither a `sage` executable nor WSL is
installed", and that "Sage-only tests therefore skip rather than silently using
the Python backend". It carried a whole section of WSL 2 installation
instructions.

None of that is true of the maintenance machine now: Sage 10.9 is installed
from conda-forge and the Sage-marked tests execute. The assessment is kept here
because it is *dated evidence* — it is what was true on 2026-08-04, and the
claim "the eight `sage` tests skip" is a measurement, not an error.

The rewritten doc records what replaced it, including the two conda-build traps
that cost real time: there is no `sage -python` subcommand, and the GAP `design`
package is absent, so `designs.WittDesign(11)` raises and `S(4,5,11)` has to be
built from the ternary Golay code instead.

## 2. The interpreter paths

`README.md` and `docs/architecture.md` spelled every command as
`C:\Python314\python.exe -m ...` — 14 occurrences in the README's "Start here"
block alone. That is a machine-specific path in a public repository's first
screen of instructions. All occurrences are now plain `python`, which is correct
on every platform including under Sage's own interpreter. Nothing else in those
commands changed.

## 3. The pytest basetemp hook

`tests/conftest.py` carried a `pytest_configure` that redirected pytest's
`basetemp` into `<repo>/../../../Claude Func Folder/pytest-fka`, with a
docstring stating that "the project contract forbids session-scoped `%TEMP%`
files".

That contract is not recorded anywhere in this repository — not in `README.md`'s
research contracts, not in `docs/`, not in `pyproject.toml`. It was an
instruction belonging to the session that wrote the code, not to the project. On
any machine without that folder the hook was already a no-op, which is why it
never showed up as a failure; on Linux it has been dead code since the move.

Removing it restores pytest's default `tmp_path` behaviour, which is what Binder
and the GitHub Actions job were using in any case. The full function is
preserved in `conftest-before.py` if the constraint turns out to be real, in
which case it belongs in `pyproject.toml` as a documented setting rather than in
a probing hook.

## Nothing is lost

The Sage assessment is a dated record and is reproduced above in full. The
removed hook is preserved verbatim. The interpreter paths are a mechanical
substitution, reversible with one `sed`.
