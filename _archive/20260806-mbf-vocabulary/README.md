# 2026-08-06 — one vocabulary, and the glossary it replaced

The repository used two vocabularies for one object: FK-A's hypergraph `G`/`H`
and FK-B's normal-form `cnf`/`dnf`, mixed inside the same layer, plus invented
third names that matched neither source. Everything now speaks the
monotone-Boolean register of the FK-B reference, and
[`../../docs/dictionary.md`](../../docs/dictionary.md) is the single place the
hypergraph reading is recorded.

## What is here

`glossary.md` — the previous four-vocabulary correspondence table. It is the one
file the change *deleted* rather than renamed, so it is kept here. It is
superseded in full by `docs/dictionary.md`, which adds the on-disk column, the
recursion-tree names, the rule that a concept has exactly one identifier, and a
change log that is only appended to after approval.

## What is not here, and why

Nothing else. Every other file in that change was renamed or edited in place, so
`git log` and `git show` carry the previous version exactly; copying them here
would duplicate what version control already holds. The renames are:

| was | now |
| --- | --- |
| `src/fka/sperner.py` | `src/fka/irredundant.py` (`Irredundant.m`) |
| `src/fka/threshold.py` | `src/fka/bound.py` (it collided with the `threshold` generator) |
| `RecursionTree`, `RecursionNode` | `Tree`, `Node` |
| `SpernerReduction` | `Reduction` |
| `HypergraphProperties` | `Properties` |

and the deletions, each a duplicate of something that survives:

| removed | kept |
| --- | --- |
| `fka.threshold.eps`, `fka.instances.exact_epsilon` | `fka.algorithm.eps_exact` |
| second `FANO_LINES` in `fka.instances` | `fka.selfdual.FANO_LINES` |
| `fka.properties.is_normal` | `is_conformal` — the same predicate |
| `remove_superset_pairwise` | `remove_superset`; the definition it was compared against now lives in the test that makes the comparison |
| `Hypergraph.contraction`, `.deletion` | `phi_x_0`, `phi_x_1` — the same operations |
| the incidence-graph automorphism search in `fka.selfdual` | `fka.automorphism.term_set_automorphisms`, which handles one term set or a pair |
| `report.ALGORITHM_STYLE`'s per-algorithm side names | one vocabulary left nothing to reconcile |

## The on-disk schema moved with the code

`data/instances/*.json` and `results/*.json` were regenerated: `G`/`H` became
`cnf`/`dnf`, `n_edges_G`/`n_edges_H` became `nC`/`nD`, `pivot` became
`split_var`, `certificate` became `CA`, `removed_G`/`removed_H` became
`cut_C`/`cut_D`. Node counts were checked against the committed baselines before
and after: all fifteen instances reproduce them exactly under both algorithms.

## One correctness fix went in at the same time

`fka.groups._quaternion8` was two undocumented permutation constants that do not
generate Q8 — they generate an abelian group, in fact C2 x C4. Q8 therefore came
back `unidentified`, and C2 x C4 was misnamed `Q8`. No committed result
contained either group, so nothing published was wrong; the bug was latent. It
is now built from `i^2 = j^2 = k^2 = ijk = -1`.

Separately, `_ALIASES` rewrote whole names only, so the same group printed as
`S3` standalone and `D3` inside a product. Aliasing is now per factor, which
**changes published group names**: `D3 wr C2` → `S3 wr C2`, `D3 x D4` →
`S3 x D4`, `D3 x S4` → `S3 x S4`.
