# Results

Every run of either algorithm, in one place. Regenerable; nothing here is
authored by hand.

## Per-instance runs

`<instance>.<algorithm>.<variant>.<ext>`, over 16 instances and both algorithms:

| algorithm | variant | committed | |
| --- | --- | :-: | --- |
| `fk-a` | `faithful` | yes | FK96 as stated |
| `fk-a` | `modified` | yes | thesis §5.1.1 `D_{1,k}`; a different tree on every instance |
| `fk-b` | `faithful` | yes | `FK_B.m` |
| `fk-b` | `multiple` | no | `MFK_B.m`; same tree, see below |

Three formats each: `.json` is the complete machine-readable recursion tree —
every node's subproblem, verdict, conflicting assignment and analysis; `.html` is a
self-contained interactive explorer with structure, automorphism and property
views, needing no network or font; `.dot` is optional Graphviz interchange.

Regenerate the committed set with:

```text
python -m fka run --all --all-variants --dot
python -m fkb run --all --dot
```

### `word-xy` and `word-yx` are structure only

Their artifacts carry no `analysis` block, and the `Annotated` chip on the HTML
says so. Both sides sit under `skip_group_above = 40` — 32 and 16 terms — so the
annotation pass attempts the automorphism group, but the group is past
`MAX_GROUP_ORDER` and every distinct node spends the full enumeration budget to
discover that and return nothing. Measured: about 40 s per distinct node against
under a second for every other instance, and the trees have 865 and 577 nodes.
The cap is on term count and the cost is driven by variable count and group
order, which is the gap; until that is closed these two are generated with
`--no-annotate`:

```text
python -m fka run word-xy word-yx --all-variants --dot --no-annotate
python -m fkb run word-xy word-yx --dot --no-annotate
```

`--all-variants` works on both and writes every variant; for FK-B that adds a
set this repository deliberately does not carry, for the reason below. A single
configuration is `python -m fka run <id> --variant faithful`. Backend
attribution is embedded in annotated JSON.

The JSON is written compact rather than pretty-printed. A tree is a machine
artifact regenerated wholesale and never hand-edited, and indenting one costs
about 2.7x the bytes — 7.4 MB against 2.7 MB for `sdfp-sd-2` under FK-A.

### Why `fk-b.multiple` is not committed

`MFK_B.m` differs from `FK_B.m` only in how many conflicting assignments it
reports, never in control flow, so the two produce identical trees. Every live
instance is a dual pair, so no node produces a conflicting assignment at all — leaving the `.dot` byte-identical, the `.html` differing by one header
chip, and the `.json` by one field. Committing it duplicated 1.4 MB and invited
readers to look for a difference that is not there.

Generate it any time with `python -m fkb run --all --variant multiple --dot`.
The variant keeps its baseline in every instance file and its own row in
`python -m fkb verify`. The withdrawal, and the measurement behind it, is
[`../_archive/20260804-redundant-fkb-multiple-artifacts/`](../_archive/20260804-redundant-fkb-multiple-artifacts/).

## Cross-algorithm report

- `fk-a-vs-fk-b.json` and `fk-a-vs-fk-b.md` — node counts, automorphism-group
  distributions, and symmetry retention for both algorithms on every instance,
  plus `sdfp-sd-3`, which is too large to keep in the library.
- `hard-case-search.json` — every vertex-transitive self-transversal hypergraph
  on at most 15 points, with its group, primitivity, frequency and both
  algorithms' node counts. Regenerate with
  `sage -python experiments/hard_case_search.py 3 15`; the reading is in
  [`../docs/hard-cases.md`](../docs/hard-cases.md).

Regenerate with `python experiments/compare_algorithms.py`. Every row is checked
against the brute-force oracle before it is reported.

The reading of these numbers is in [`../docs/fk-a-vs-fk-b.md`](../docs/fk-a-vs-fk-b.md).
