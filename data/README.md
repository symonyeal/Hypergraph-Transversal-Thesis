# Data

`instances/` is the canonical instance library, shared by FK-A and FK-B. Each
JSON file contains the input pair in 1-indexed vertex notation, provenance, a
source citation, research notes, and independently refreshable expected values.

```json
{
  "id": "fano",
  "provenance": "verbatim" | "corrected" | "derived",
  "n_vertices": 7,
  "G": [[1, 2, 4], ...],
  "H": [[1, 2, 4], ...],
  "notes": "...",
  "expected": {
    "dual": true,                  from the brute-force oracle, never from FK-A or FK-B
    "epsilon": "3/7",              exact fraction
    "n_edges_G": 7, "n_edges_H": 7,
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
