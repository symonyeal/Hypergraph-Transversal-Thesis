# Results

Every run of either algorithm, in one place. Regenerable; nothing here is
authored by hand.

## Per-instance runs

`<instance>.<algorithm>.<variant>.<ext>` — the full cross product of 12
instances, 2 algorithms, and every variant of each:

| algorithm | variants |
| --- | --- |
| `fk-a` | `faithful` (FK96 as stated), `modified` (thesis §5.1.1 `D_{1,k}`) |
| `fk-b` | `faithful` (`FK_B.m`), `multiple` (`MFK_B.m`) |

Three formats each: `.json` is the complete machine-readable recursion tree —
every node's subproblem, verdict, certificate and analysis; `.html` is a
self-contained interactive explorer with structure, automorphism and property
views, needing no network or font; `.dot` is optional Graphviz interchange.

Regenerate the lot with:

```text
python -m fka run --all --all-variants --dot
python -m fkb run --all --all-variants --dot
```

or a single configuration with `python -m fka run <id> --variant faithful`.
Backend attribution is embedded in annotated JSON.

The JSON is written compact rather than pretty-printed. A tree is a machine
artifact regenerated wholesale and never hand-edited, and indenting one costs
about 2.7x the bytes — 7.4 MB against 2.7 MB for `sdfp-sd-2` under FK-A.

### On `fk-b.multiple`

`MFK_B.m` differs from `FK_B.m` only in how many conflicting assignments it
reports, never in control flow, so the two produce identical trees. Every live
instance is a transversal pair, so no node produces a conflicting assignment at
all, and `fk-b.multiple` therefore differs from `fk-b.faithful` **only in the
`variant` label** — all 367 node payloads are byte-identical. It is generated
for completeness; do not read a difference into it. The claim is asserted by
`tests/test_fkb.py::test_multiple_variant_leaves_the_tree_alone`.

## Cross-algorithm report

- `fk-a-vs-fk-b.json` and `fk-a-vs-fk-b.md` — node counts, automorphism-group
  distributions, and symmetry retention for both algorithms on every instance,
  plus `sdfp-sd-3`, which is too large to keep in the library.

Regenerate with `python experiments/compare_algorithms.py`. Every row is checked
against the brute-force oracle before it is reported.

The reading of these numbers is in [`../docs/fk-a-vs-fk-b.md`](../docs/fk-a-vs-fk-b.md).
