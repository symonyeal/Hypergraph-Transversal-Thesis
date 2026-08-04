# Architecture and research workflow

## Why the maintained model is different

The legacy code represented every hyperedge as a mutable NumPy row and mixed
preprocessing, recursion, analysis, plotting, and notebook state. The maintained
model stores an edge as an integer bitset and a hypergraph as an immutable,
canonical tuple. Intersection, subset, deletion, and contraction are exact
integer operations, while incidence matrices remain a lossless interchange
format for SageMath and MATLAB.

## Layout

Two algorithms, one model. `src/fka` holds FK-A *and* everything neither
algorithm owns; `src/fkb` holds FK-B and imports that model. The dependency is
one-way and is about where files sit, not about one algorithm needing the other.

| Module | Owns |
| --- | --- |
| `fka/hypergraph.py`, `fka/sperner.py` | exact data operations |
| `fka/tree.py` | the serializable `RecursionTree`, shared by both algorithms |
| `fka/transversal.py` | the independent Berge-method oracle |
| `fka/analysis.py` | annotating completed nodes, memoizing repeated subproblems |
| `fka/report.py`, `fka/dot.py` | rendering a stored tree without changing it |
| `fka/instances.py` | the instance library, its generators, and its baselines |
| `fka/algorithm.py` | FK-A |
| `fkb/algorithm.py` | FK-B |
| `fkb/dualize.py` | `Tr(H)` by repeated FK-B equivalence tests |

Neither algorithm plots, computes automorphism groups, or touches the filesystem
while recursing. That separation is what makes a result reproducible and
reviewable months later, and it is what lets one visualiser serve both.

Start with [glossary.md](glossary.md) if the two algorithms' vocabularies are
not yet interchangeable in your head. A monotone Boolean function is a
hypergraph and its dual is the transversal hypergraph; that identification is
what lets one model, one library and one visualiser serve both.

## Adding a research instance

1. Add one entry to `experiments/build_instances.py`. Use a short stable `id`,
   1-indexed edge lists, a precise source, and one of `verbatim`, `corrected`,
   or `derived` for provenance.
2. If a printed source needs repair, keep the as-printed form and state exactly
   what changed and why in the corrected file's `notes`.
3. Rerun `C:\Python314\python.exe experiments/build_instances.py`. It recomputes
   every `expected` block: duality from the oracle, `epsilon` as an exact
   fraction, and node counts for both algorithms and all four variants. Review
   every changed value; never accept the diff mechanically.
4. Run `C:\Python314\python.exe -m pytest`, then
   `C:\Python314\python.exe -m fka verify` and
   `C:\Python314\python.exe -m fkb verify`.
5. Generate the deliverables with `-m fka run <id>` and `-m fkb run <id>`.
6. Put a durable mathematical claim in a focused test and cite the thesis page
   or archived source log in its docstring.

Prefer a **generator** to a copied edge list when the instance has a stated
construction: `fka.instances` builds the Fano plane, the self-dualized Fano
benchmarks, and the matching and threshold families from their definitions, so
they are reproducible and carry no third-party redistribution question.

## Withdrawing an instance

An instance that turns out not to be a usable research input is moved to
`_archive/<YYYYMMDD>-<slug>/`, not deleted, and the directory gets a `README.md`
saying what and why. `fka.instances.load_archived` reads it back, so any finding
that depended on it stays executable; `list_ids` does not include it, so nothing
runs against it by accident. See
`_archive/20260804-verbatim-nontransversal/` for the worked example.

## Where things are

Two folders, and each holds all of its kind:

- **`data/`** — every research input. `instances/` is the `(G, H)` library,
  `reference/` the MATLAB FK-B authors' recorded test vectors, `baselines/` the
  node counts and split sequences printed in the thesis. No test embeds its own
  data; `fka.instances.load_data` reads it. Superseded inputs move to
  `_archive/` and return through `load_archived`.
- **`results/`** — every run of both algorithms and every variant of each.

## Result formats

- JSON is the complete machine-readable recursion tree and analysis record,
  written compact: a tree is regenerated wholesale and never hand-edited, and
  indenting one costs about 2.7x the bytes.
- HTML is a single-file interactive view supporting structure, automorphism, and
  property views; clicking a node reveals its exact subproblem and certificate.
- DOT is optional interchange for Graphviz. It is not required for HTML.

Files are `<instance>.<algorithm>.<variant>.<extension>` — for example
`fano.fk-a.modified.html` and `fano.fk-b.faithful.html` — so rerunning a
configuration replaces its own deterministic artifacts rather than creating
numbered copies, and neither the two algorithms nor their variants overwrite
each other. Regenerate everything with

```text
python -m fka run --all --all-variants --dot
python -m fkb run --all --all-variants --dot
```

`results/fk-a-vs-fk-b.{json,md}` is the cross-algorithm report, regenerated by
`experiments/compare_algorithms.py`.

## Annotation cost

The automorphism search and the primal-graph classification are exponential in
the worst case, and `AnalysisOptions` caps both. The defaults are set from
measurement: on `sdfp-sd-2` the group search costs about 19 seconds per distinct
node against about a millisecond on every thesis instance, so a cap that admits
it makes a thousand-node tree unannotatable. Raise the caps deliberately, per
run, when a specific group is the question — not as a default.

## Validation layers

- Unit tests cover representation, Sperner reduction, tree serialization,
  groups, graph classes, properties, and output safety.
- Randomized tests compare both algorithms, all variants, and all pivot and
  split rules against the independent transversal oracle, and verify that every
  certificate and conflicting assignment holds at the root.
- Cross-algorithm tests assert FK-A and FK-B never disagree.
- Thesis tests lock published node counts and, where available, full pivot
  sequences.
- Conformance tests replay the MATLAB FK-B reference's own recorded test
  vectors from `Unit_Tests_Equivalency.m` and `Unit_tests_Check_Conditions.m`.
- Sage-marked tests cross-check automorphism orders and graph classes using
  Sage's independent implementations.

The Sage tests are expected to skip in ordinary CPython. A release or quoted
cross-backend result should also be run in the Sage environment described in
`docs/sagemath.md`.
