# Architecture and research workflow

## Why the maintained model is different

The legacy code represented every hyperedge as a mutable NumPy row and mixed
preprocessing, recursion, analysis, plotting, and notebook state. The maintained
model stores an edge as an integer bitset and a hypergraph as an immutable,
canonical tuple. Intersection, subset, deletion, and contraction are exact
integer operations, while incidence matrices remain a lossless interchange
format for SageMath and MATLAB.

The components are deliberately separated:

1. `hypergraph.py` and `sperner.py` define exact data operations.
2. `algorithm.py` builds a serializable `RecursionTree`; it does not plot or
   calculate automorphism groups while recursing.
3. `transversal.py` supplies an independent Berge-method oracle.
4. `analysis.py` annotates completed nodes, memoizing repeated subproblems.
5. `report.py` and `dot.py` render the stored tree without changing it.
6. `data/instances/*.json` holds named inputs, provenance, notes, and baselines.

This separation is what makes a result reproducible and reviewable months later.

## Adding a research instance

1. Add one JSON file under `data/instances/`. Use a short stable `id`, 1-indexed
   edge lists, a precise source, and one of `verbatim`, `corrected`, or `derived`
   for provenance.
2. If a printed source needs repair, keep a separate `-verbatim` file and state
   exactly what changed and why in the corrected file's `notes`.
3. Recompute baselines with `C:\Python314\python.exe -m fka refresh` after the
   editable install. Review every changed expected value; never accept the diff
   mechanically.
4. Run `C:\Python314\python.exe -m pytest` and
   `C:\Python314\python.exe -m fka verify`.
5. Generate the deliverable with `C:\Python314\python.exe -m fka run <id>`.
6. Put a durable mathematical claim in a focused test and cite the thesis page
   or archived source log in its docstring.

## Result formats

- JSON is the complete machine-readable recursion tree and analysis record.
- HTML is a single-file interactive view. It supports structure,
  automorphism, and property views; clicking a node reveals its exact
  subproblem and certificate.
- DOT is optional interchange for Graphviz. It is not required for HTML.

Files use `<instance>.<variant>.<extension>`, so rerunning a configuration
replaces its own deterministic artifacts rather than creating numbered copies.

## Validation layers

- Unit tests cover representation, Sperner reduction, tree serialization,
  groups, graph classes, properties, and output safety.
- Randomized tests compare both variants and both pivot rules with the
  independent transversal oracle.
- Thesis tests lock published node counts and, where available, full pivot
  sequences.
- Sage-marked tests cross-check automorphism orders and graph classes using
  Sage's independent implementations.

The Sage tests are expected to skip in ordinary CPython. A release or quoted
cross-backend result should also be run in the Sage environment described in
`docs/sagemath.md`.
