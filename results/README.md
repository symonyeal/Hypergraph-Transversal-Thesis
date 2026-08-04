# Results

Regenerable research artifacts. Nothing here is authored by hand.

## Per-instance runs

`<instance>.<algorithm>.<variant>.<ext>`, for example `fano.fk-a.modified.html`
and `fano.fk-b.faithful.json`:

- `.json` — the complete machine-readable recursion tree, including each node's
  subproblem, verdict, certificate, and analysis;
- `.html` — a self-contained interactive explorer, no network or font
  dependency, with structure, automorphism, and property views;
- `.dot` — optional Graphviz interchange.

Regenerate one instance with `python -m fka run <id>` or `python -m fkb run <id>`,
or the whole library with `--all`. Backend attribution is embedded in annotated
JSON.

## Cross-algorithm report

- `fk-a-vs-fk-b.json` and `fk-a-vs-fk-b.md` — node counts, automorphism-group
  distributions, and symmetry retention for both algorithms on every instance,
  plus `sdfp-sd-3`, which is too large to keep in the library.

Regenerate with `python experiments/compare_algorithms.py`. Every row is checked
against the brute-force oracle before it is reported.

The reading of these numbers is in [`../docs/fk-a-vs-fk-b.md`](../docs/fk-a-vs-fk-b.md).
