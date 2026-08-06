# Data

Every research input, in one place. Nothing that either algorithm or the test
suite reads lives anywhere else.

| Folder | Holds | Read by |
| --- | --- | --- |
| `instances/` | the `(cnf, dnf)` library — 16 named problem instances | `fka.instances.load`, `load_all` |
| `reference/` | the MATLAB FK-B authors' own recorded test vectors | `tests/test_fkb_reference.py` |
| `baselines/` | node counts and split sequences printed in the thesis | `tests/test_thesis_reproduction.py` |

Read any of it with `fka.instances.load_data("baselines", "published-trees.json")`.
Superseded inputs are not deleted; they move to `../_archive/<date>-<slug>/`
with a note, and come back through `fka.instances.load_archived`.

`cnf` and `dnf` are the same two term sets whichever algorithm reads them; the
hypergraph reading is in [`../docs/dictionary.md`](../docs/dictionary.md).

## `instances/`

Each JSON file contains the input pair in 1-indexed variable notation, provenance,
a source citation, research notes, and independently refreshable expected values.

```json
{
  "id": "fano",
  "provenance": "verbatim" | "corrected" | "derived",
  "n_vertices": 7,
  "cnf": [[1, 2, 4], ...],
  "dnf": [[1, 2, 4], ...],
  "notes": "...",
  "expected": {
    "dual": true,                  from the brute-force oracle, never from FK-A or FK-B
    "epsilon": "3/7",              exact fraction
    "nC": 7, "nD": 7,
    "fka": {"faithful": {...}, "modified": {...}},
    "fkb": {"faithful": {...}, "multiple": {...}}
  }
}
```

The `expected` block is **not** hand-written. `experiments/build_instances.py`
recomputes it and `tests/test_thesis_reproduction.py` asserts the committed
values still hold, so it is a regression baseline rather than a claim: if a
change to either algorithm moves a node count, the diff says so.

## Provenance

| Value | Meaning |
| --- | --- |
| `verbatim` | Transcribed exactly as published. |
| `corrected` | A published instance with a stated, minimal repair. The correction and its evidence are in `notes`, and the as-printed form is preserved under `_archive/`. |
| `derived` | Constructed here — from a generator, or recovered from a run log. |

Two thesis instances are not transversal pairs as printed. Their repairs are
`6a` and `8ver`; the as-printed forms were withdrawn to
`../_archive/20260804-verbatim-nontransversal/` on 2026-08-04 and are still
re-derived by the test suite via `fka.instances.load_archived`.

Prefer a generator to a copied edge list when the instance has a stated
construction. `fka.instances` builds the Fano plane, the self-dualized Fano
benchmarks (`sdfp-sd-*`), and the matching and threshold families from their
definitions.

Do not store executable Python, notebook state, or rendered output here.
Generated trees belong in `results/`; programs that derive new data belong in
`experiments/` and should write reviewed JSON here only when explicitly run.
