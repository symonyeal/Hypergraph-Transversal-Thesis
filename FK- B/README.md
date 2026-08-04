# FK-B reference material

`FK-master/` is a third-party MATLAB implementation of FK-B by Nafiseh
Sedaghat. Its redistribution provenance and license have not been established,
so the snapshot stays **local and git-ignored** — see `.gitignore`.

## What was taken from it, and how

`src/fkb/` is an independent Python implementation written from the algorithm as
published (*How to apply SAT-solving for the equivalence test of monotone normal
forms*, p.5) and from this snapshot's stated behaviour. No file here is copied,
translated line-for-line, or redistributed. What crossed over is:

| From | To | Form |
| --- | --- | --- |
| `FK_B.m`, `MFK_B.m` | `fkb/algorithm.py` | the algorithm, as the `faithful` and `multiple` variants |
| `Check_Conditions.m`, `Easy_case.m`, `CA_CNF2.m` | `fkb/algorithm.py` | the three conditions and the small-case decision, re-derived rather than transcribed |
| `FK_Dualization.m`, `Maximum_False_Point.m`, `Minimality_Check.m`, `Initialize_CNF.m` | `fkb/dualize.py` | dualization by repeated equivalence testing |
| `StandardProblems.m` | `fka.instances.matching`, `.threshold` | generators |
| `SDFP16/23_CNF_DNF.mat` | `fka.instances.self_dual_fano` | the *construction*, reproduced exactly; the data is not copied |
| `Unit_Tests_Equivalency.m`, `Unit_tests_Check_Conditions.m` | `tests/test_fkb_reference.py` | the recorded input/answer vectors, as conformance tests |

`Irredundant.m` needed nothing: it is Sperner reduction, which `fka.sperner`
already had. `berge.m` likewise duplicates `fka.transversal`.

## Deviations, all of them corrections

Enumerated with a test each in `fkb/algorithm.py`'s docstring: kept empty terms
where `phi_x_0.m` deletes them, constants as ordinary hypergraphs where
`A_c_x.m`/`A_m_x.m` return a MATLAB scalar, and `None` rather than `[]` for "no
conflicting assignment".

Two observations about the snapshot itself, recorded because they affect how its
outputs should be read:

- Seven of its twelve `*_CNF_DNF.mat` benchmark files store a `cnf` and `dnf`
  with different column counts, so those two matrices are not over the same
  variable set and are not a transversal pair. Run
  `python experiments/import_fkb_benchmarks.py` to reproduce that finding.
- The last two vectors in `Unit_tests_Check_Conditions.m` carry a bare `% 1`
  comment that is not a conflicting assignment in either direction. The file has
  no assertions, so these read as stale notes.

## Nothing depends on the snapshot

`src/`, the instance baselines, and both test suites run from a clean checkout
with `FK- B/FK-master/` absent. Only `experiments/import_fkb_benchmarks.py`
reads it, and it exits cleanly when it is not there.
