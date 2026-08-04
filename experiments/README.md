# Experiments

Explicit, rerunnable research and migration programs that do not belong in the
importable library.

| Program | What it does |
| --- | --- |
| `build_instances.py` | Authors `data/instances/` and the withdrawn instances under `_archive/`. Every `expected` block is computed here, not typed. |
| `compare_algorithms.py` | Runs FK-A and FK-B on the same instances, annotates both trees, and writes `results/fk-a-vs-fk-b.{json,md}`. |
| `import_fkb_benchmarks.py` | Reads the local MATLAB FK-B snapshot's `*_CNF_DNF.mat` files and reports what is in them. Requires SciPy and the snapshot; neither is a project dependency. |

Core correctness claims belong in `tests/`, and normal runs belong in the `fka`
and `fkb` CLIs. Avoid adding notebooks or unstructured output logs.

Each program is deterministic and safe to rerun. Any that writes into
`data/instances/` prints what changed, so the diff can be reviewed rather than
accepted mechanically.
