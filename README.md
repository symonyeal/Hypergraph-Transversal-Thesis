# Hypergraph dualization: FK-A and FK-B

[![Tests](https://github.com/symonyeal/Hypergraph-Transversal-Thesis/actions/workflows/tests.yml/badge.svg)](https://github.com/symonyeal/Hypergraph-Transversal-Thesis/actions/workflows/tests.yml)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/symonyeal/Hypergraph-Transversal-Thesis/main?urlpath=lab)

This is the maintained research implementation for Saimon Y. Islam's work on
hypergraph dualization, symmetry, and the Fredman-Khachiyan algorithms. The
mathematical reference is [Saimon Islam Manuscript MSc Thesis (1).pdf](Saimon%20Islam%20Manuscript%20MSc%20Thesis%20(1).pdf).
The definitive repository is
[symonyeal/Hypergraph-Transversal-Thesis](https://github.com/symonyeal/Hypergraph-Transversal-Thesis).

Both algorithms decide the same question against one shared model, one instance
library, and one visualiser, so their recursion trees are directly comparable.
They state it in different languages, and the identification is the point:

```
G = Tr(H)          hypergraph dualization        FK-A's statement
f^d = g            monotone Boolean duality
C ≡ D              a CNF and a DNF agree         FK-B's statement
```

A monotone Boolean function **is** a hypergraph and its dual **is** the
transversal hypergraph, so the code is written in **one** vocabulary — the
monotone-Boolean register of the FK-B reference. `cnf` is the `G` side, `dnf`
the `H` side, and `fk_a(cnf, dnf)` and `fk_b(cnf, dnf)` take the same arguments
and return the same answer. [docs/dictionary.md](docs/dictionary.md) is the full
correspondence, and the only place the hypergraph names are recorded.

| | package | source | what it returns |
| --- | --- | --- | --- |
| **FK-A** | `src/fka/` | thesis Algorithm 1, p.13 | a conflicting assignment `S` satisfying equation 2.1 |
| **FK-B** | `src/fkb/` | *How to apply SAT-solving for the equivalence test of monotone normal forms*, p.5 | a conflicting assignment `S` with `C(S) ≠ D(S)` |

What the project provides:

- an exact bitset term-set model;
- faithful and thesis-modified FK-A variants, and a faithful FK-B port with its
  single- and multiple-assignment variants;
- FK-B run as a *generator* of `Tr(D)`, not only as a decision procedure;
- an independent transversal oracle and checkable conflicting assignments;
- structured, provenance-bearing JSON instances;
- deterministic JSON, DOT, and self-contained interactive HTML results;
- a pytest suite tied to the published recursion trees and to the MATLAB
  reference's own recorded test vectors; and
- an optional SageMath backend for independent group and graph checks.

## Start here

From this directory, the source-tree tests require no installation:

```text
python -m pytest
```

Install the command-line package into whichever interpreter you use (no project
virtual environment is required):

```text
python -m pip install -e ".[dev]"
```

Then use the same interpreter consistently. Under SageMath, substitute that
environment's Python and the Sage-backed checks come with it — see
[docs/sagemath.md](docs/sagemath.md):

```text
python -m fka env            what this interpreter can do
python -m fka list           the instance library
python -m fka show fano      one instance, with its baseline
python -m fka run fano       FK-A + annotation + HTML report
python -m fka verify         FK-A against the oracle

python -m fkb run fano       FK-B, same reports
python -m fkb verify         FK-B against the oracle
python -m fkb compare --all  FK-A vs FK-B vs the oracle
python -m fkb benchmark      both, on the scaling families
python -m fkb dualize --all  Tr(D) by repeated FK-B tests
```

`run` writes a portable HTML explorer and machine-readable JSON to `results/`.
Open the HTML directly; it has no network, font, JavaScript, or image dependency.

## Research contracts

- **Duality is decided by the oracle, never by either algorithm.**
  `fka.transversal.is_dual_oracle` builds `Tr(dnf)` by Berge's method and
  compares. A passing FK-A or FK-B run alone is not a ground-truth proof.
- **The two algorithms must never disagree.** They decide the same predicate;
  `tests/test_fkb.py` asserts agreement with each other and with the oracle on
  the library, on randomised pairs, and under every variant and split rule.
- FK-A: `variant="faithful"` implements the FK96 base case; `variant="modified"`
  adds the thesis Section 5.1.1 `D_{1,k}` / `D_{k,1}` rule.
  `split_rule="degree_sum"` reproduces the published tree numbering;
  `split_rule="frequency"` implements the normalized frequency quantity bounded
  in the thesis.
- FK-B: `variant="faithful"` is `FK_B.m`; `variant="multiple"` is `MFK_B.m` and
  returns every conflicting assignment a node can see, leaving the tree
  unchanged. `split_rule="most_frequent"` is `Choose_SplitVar.m`;
  `split_rule="degree_sum"` matches FK-A so both can run on identical splits.
- Changes to preprocessing, split choice, base cases, or recursion order must
  preserve the published node counts and split sequences in
  `tests/test_thesis_reproduction.py` unless the research interpretation itself
  is deliberately revised.
- Corrections never overwrite a printed instance. `6a` and `8ver` record the
  corrected pairs and the evidence in their JSON notes; the as-printed forms,
  which are **not** transversal pairs, are preserved under
  [`_archive/20260804-verbatim-nontransversal/`](_archive/20260804-verbatim-nontransversal/)
  and are still re-derived by the test suite.
- Deviations from the MATLAB FK-B reference are corrections, are enumerated in
  `fkb.algorithm`'s docstring, and each has a test naming the case.

## Filing system

| Path | Purpose |
| --- | --- |
| `src/fka/` | FK-A, and the shared model: hypergraph, tree, analysis, report, instances |
| `src/fkb/` | FK-B, and dualization by repeated equivalence testing |
| `data/` | Every research input: `instances/`, the MATLAB `reference/` vectors, thesis `baselines/` |
| `tests/` | Unit, randomized, thesis-reproduction, MATLAB-conformance, output, and Sage checks |
| `results/` | Every run of both algorithms and every variant, plus the cross-algorithm report |
| `experiments/` | Explicit migration or one-off research programs |
| `docs/` | Architecture, the FK-A/FK-B findings, the hard-case search, Sage setup, and assessment notes |
| `Reading/` | Local-only third-party reference library; PDFs are not redistributed |
| `FK- B/` | Local-only third-party MATLAB FK-B reference snapshot |
| `Dedekind/` | Third-party MIT snapshot of Jäkel's `d(9)` computation — the size of the space both algorithms search |
| `_archive/` | Preserved superseded material, each with an archive note |

`src/fkb` imports the shared model from `src/fka`; the dependency is one-way and
is about layout, not about one algorithm needing the other. Do not add new
research logic to notebooks or `.txt` logs. Add a module or an explicit program,
store inputs as JSON, and make the claim executable as a test. See
[docs/architecture.md](docs/architecture.md) for the repeatable workflow and
[docs/fk-a-vs-fk-b.md](docs/fk-a-vs-fk-b.md) for what the comparison found.
[docs/hard-cases.md](docs/hard-cases.md) reports the search for new hard
instances and the one it found: the Witt design `S(4,5,11)`, in the library as
`w11`. [docs/hard-classes.md](docs/hard-classes.md) states what is settled about
the families that make FK-A's frequency bound tight, and which kind of symmetry
causes it.

## SageMath

The package runs completely under CPython. Sage is an optional independent
backend, not an import-time dependency. Use the Binder badge above to run Sage
10.9 in a browser without installing it locally, and GitHub Actions runs the
same Sage-backed checks for every definitive revision.

Sage 10.9 is installed on the maintenance workstation, so the Sage-marked tests
execute rather than skip: 373 passed / 1 skipped under Sage's Python against
361 passed / 13 skipped under plain CPython.
[docs/sagemath.md](docs/sagemath.md) has the conda recipe and the two traps in
that build (no `sage -python`, no GAP `design` package).

## Authorship and license

The legacy implementation credits Saimon Islam and Lacey Liang and declares
GPL-3.0 licensing. The archived original notice is retained in
`_archive/20260804-legacy-code/LEGACY README.md`. The maintained package carries
the same `GPL-3.0-or-later` project metadata.

FK-B is an independent Python implementation written from the algorithm as
published and from the MATLAB reference's stated behaviour. The benchmark
families it is measured on are generated by `fka.instances`, not copied, and
the MATLAB snapshot in `FK- B/` stays local and git-ignored because its
redistribution terms are not established.

`Dedekind/` is the one vendored third party, kept because it carries an explicit
MIT licence permitting it. It is not covered by this project's GPL, nothing in
`src/` links against it, and its copyright notice travels with it — see
[Dedekind/README.md](Dedekind/README.md).
