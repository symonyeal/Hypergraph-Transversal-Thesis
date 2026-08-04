# Withdrawn 2026-08-04 — the `fk-b.multiple` result artifacts

## What was removed

`results/<instance>.fk-b.multiple.{json,html,dot}` — 36 files, 1.4 MB, one set
per instance in the live library.

They were generated when `results/` was expanded to the full cross product of
instances, algorithms and variants. They carry no information that
`fk-b.faithful` does not.

## Why they are redundant

`MFK_B.m` differs from `FK_B.m` only in how many conflicting assignments it
reports when it finds one. Its `Multiple_Check_Conditions` and `MEasy_case`
return non-empty in exactly the same situations as their single-answer
counterparts, so control flow — and therefore the recursion tree — is identical.
That is asserted directly by
`tests/test_fkb.py::test_multiple_variant_leaves_the_tree_alone`, which compares
node count, every pivot and every path.

Every instance in the live library is a transversal pair, so **no node produces a
conflicting assignment at all**, and the one field that could have differed is
absent from both. What remains is a label:

- `.dot` — byte-identical; `fka.dot.to_dot` does not emit the variant;
- `.html` — one line differs, the `Variant` chip in the header;
- `.json` — one field differs, `"variant"`.

Measured across all 12 instances and all 3 formats on 2026-08-04:

| instance | .json | .html | .dot |
| --- | --- | --- | --- |
| `6a` | 1 line differs | 1 line differs | byte-identical |
| `6b` | 1 line differs | 1 line differs | byte-identical |
| `6c` | 1 line differs | 1 line differs | byte-identical |
| `6d` | 1 line differs | 1 line differs | byte-identical |
| `7ver` | 1 line differs | 1 line differs | byte-identical |
| `8ver` | 1 line differs | 1 line differs | byte-identical |
| `f2g2` | 1 line differs | 1 line differs | byte-identical |
| `fano` | 1 line differs | 1 line differs | byte-identical |
| `k3` | 1 line differs | 1 line differs | byte-identical |
| `sdfp-sd-1` | 1 line differs | 1 line differs | byte-identical |
| `sdfp-sd-2` | 1 line differs | 1 line differs | byte-identical |
| `trivial-aut-1` | 1 line differs | 1 line differs | byte-identical |

Independently, all 367 node payloads compare equal between the two variants.

## Nothing is lost

These are regenerable in one command, and the variant remains fully supported —
this is a decision about what is worth *committing*, not a removal of
capability:

```text
python -m fkb run --all --variant multiple --dot
```

The variant also keeps its baseline in every instance file
(`expected.fkb.multiple`) and its own row in `python -m fkb verify`, which runs
both variants against the oracle. Those are two lines of JSON and a test
respectively, they make `refresh_expected` uniform over the variant tuple, and
they are what would flag it if `multiple` ever stopped agreeing with `faithful`.
Only the bulk artifacts went.

## What was kept, and why it is not the same case

`fk-a.faithful` and `fk-a.modified` are both committed. They are *not*
redundant: the thesis §5.1.1 `D_{1,k}` rule changes the tree, and the two
variants disagree on node count for **every** instance in the library — 3382
nodes against 1686 in total, 2899 against 1427 on `sdfp-sd-2` alone.
