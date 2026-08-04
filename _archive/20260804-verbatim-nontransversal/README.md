# Archived 2026-08-04 — thesis instances that are not transversal pairs

## What is here

The two instances of Islam (2024) Section 5.2.1 that, **as printed**, are not
transversal pairs:

| File | Thesis | Defect |
| --- | --- | --- |
| `6a-verbatim.json` | 6-A, p.49 | `Tr(G)` has five edges, not four. It also contains `{2,6}` and `{1,3,6}`. `H` contains `{2,4,6}`, a proper superset of `{2,6}`, so `H` is not an antichain of minimal transversals. |
| `8ver-verbatim.json` | 8-Ver, p.50 | `H` lists 19 edges; `Tr(G)` has 20. `H` omits `{1,2,3}` and `{1,6,7}`, and includes `{1,3,7}`, which is not a transversal of `G` at all — it misses the edge `{2,5,6}`. |

Both were verified non-dual by `fka.transversal.is_dual_oracle`, which builds
`Tr(H)` by Berge's method and compares, independently of FK-A and FK-B.

## Why they were archived rather than deleted

They are the published data. Deleting them would erase the evidence for the
repair; keeping them in `data/instances/` implied they were usable research
inputs, which they are not — every run against them answers "not dual", which is
a fact about the transcription, not about the algorithms.

The corrected instances **`6a`** and **`8ver`** remain live in
`data/instances/`. Each carries `"provenance": "corrected"` and states the
change and its evidence in its `notes`. The corrections are justified by the
thesis' own reported values: 6-A is quoted at `ε = 1/2` and 8-Ver at `ε = 2/5`
(p.55), and only the repaired forms reproduce those. Verbatim 8-Ver gives
`ε = 8/19`.

## The finding is still executable

`tests/test_thesis_reproduction.py::test_archived_verbatim_instances_are_not_transversal_pairs`
loads these two files from this directory and re-derives the defect on every
run. If a future edit "fixes" them in place, or the corrected instances drift
away from the thesis' reported ε, that test fails and says why.

Nothing else reads this directory. The instance library, the baselines, and both
algorithms' test suites are unaffected.
