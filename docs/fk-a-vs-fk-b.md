# FK-A against FK-B: what the comparison shows

Findings from running both algorithms on the same instance library, with the
same hypergraph model and the same annotation pass. Every number here is
regenerable:

```text
python experiments/compare_algorithms.py   -> results/fk-a-vs-fk-b.{json,md}
python -m fkb benchmark                    -> the scaling families
python -m fkb compare --all                -> per-instance agreement
```

Duality is decided by the Berge-method oracle in every case before any statistic
is reported. **The two algorithms never disagree**, on the library, on the
benchmark families, or across 2,400 randomised pairs under every variant and
split rule. What differs is only how much work each does to get there.

---

## 1. The gap is not uniform, and it is not small

| | FK-A nodes | FK-B nodes | ratio |
| --- | ---: | ---: | ---: |
| `threshold(10)` | 9 | 27 | **0.33x** — FK-B worse |
| `8ver` | 21 | 21 | 1.00x |
| `matching(12)` | 125 | 65 | 1.92x |
| `w11-sd-1` | 955 | 719 | 1.33x |
| `w11` | 475 | 323 | **1.47x** |
| `fano` | 35 | 11 | 3.18x |
| `sdfp-sd-2` | 1,427 | 231 | 6.18x |
| `sdfp-sd-3` | 26,931 | 1,677 | **16.06x** |

FK-B is not uniformly better. It wins by a growing margin on symmetric,
self-dual instances and it *loses* by a factor of three on the threshold family.
Any claim of the form "FK-B is faster" needs the instance class attached.

`w11` is the qualification added on 2026-08-05, and it is the sharpest one: it
is symmetric *and* self-dual, the shape that should be FK-B's best case, yet the
margin collapses to 1.47x. Symmetry alone does not hand FK-B the win — see §3
and [hard-cases.md](hard-cases.md).

## 2. FK-A's blow-up is re-derived subproblems, and the rate predicts the gap

Sorting the instances by the fraction of FK-A nodes whose `(G, H)` pair already
appeared elsewhere in its own tree:

| instance | FK-A repeat rate | ratio |
| --- | ---: | ---: |
| `threshold(8..12)` | 0.0% | 0.30–0.39x |
| `6a`, `7ver` | 27–29% | 1.9–3.7x |
| `8ver` | 47.6% | 1.00x |
| `fano` | 60.0% | 3.18x |
| `f2g2` | 68.0% | 3.57x |
| `w11` | 70.5% | 1.47x |
| `sdfp-sd-1` | 78.7% | 2.03x |
| `sdfp-sd-2` | 91.2% | 6.18x |
| `sdfp-sd-3` | **98.3%** | **16.06x** |

On `sdfp-sd-3`, 26,460 of FK-A's 26,931 nodes are subproblems it has already
solved somewhere else in the same tree. The thesis noticed this in the small
case — "The resulting contraction subhypergraph at Node 8 is identical to that
at Node 17. This repetition of subproblems contributes to the enlargement of the
recursion tree" (p.53) — and at `k = 3` it is essentially the whole tree.

Where FK-A repeats nothing at all (the threshold family, 0%), FK-B has nothing
to save and its extra machinery is pure cost. **The repeat rate is the
predictor**, with one instructive exception below.

## 3. Why the repeats happen: symmetry

The repeats are not accidental. FK-A splits on a single vertex. On a
vertex-transitive instance every vertex looks the same, so the two subproblems
it produces are near-isomorphic to each other and to their parent — and FK-A
carries no memo, so each copy is explored from scratch.

The annotation pass makes this measurable. Over all annotated nodes of the
fourteen live instances, with the twelve this section originally reported on
kept alongside for comparison:

| | nodes analysed | trivial automorphism group | rate |
| --- | ---: | ---: | ---: |
| FK-A, the original twelve | 1,427 | **202** | 14.2% |
| FK-B, the original twelve | 363 | **20** | 5.5% |
| FK-A, all fourteen | 2,639 | 244 | 9.2% |
| FK-B, all fourteen | 1,395 | 144 | **10.3%** |

**Corrected 2026-08-05.** This table previously read `FK-B | 2,028 | 20` and
concluded that FK-A drives ten times as many nodes to the trivial group "on 30%
*fewer* analysed nodes". The 2,028 was not reproducible: FK-B's whole tree over
those twelve instances is 367 nodes, so it cannot have annotated 2,028 of them.
The correct figure is 363, and FK-A therefore analysed nearly **four times as
many** nodes, not fewer. What survives is the ratio: FK-A drives 14.2% of its
nodes to the trivial group against FK-B's 5.5%, a factor of 2.6 rather than 10.
On `sdfp-sd-2` alone the split is 186 against 2, which was and remains correct.

The fourteen-instance row is not the same claim, and the reversal is the
interesting part: on `w11` and `w11-sd-1` **FK-B** is the algorithm that
shreds the group (8.8% and 13.5% of its nodes reach the trivial group, against
3.5% for FK-A on both). Those two instances are the ones whose root group is
4-transitive, and there FK-B's per-term branch — the step that wins everywhere
else by removing a whole orbit representative at once — is what destroys the
symmetry fastest. See [hard-cases.md](hard-cases.md).

The same effect read the other way — mean `|Aut(G)|` over the tree, as a
fraction of the root's:

| instance | root `|Aut|` | FK-A retained | FK-B retained |
| --- | ---: | ---: | ---: |
| `fano` | 168 | 0.065 | **0.134** |
| `f2g2` | 128 | 0.102 | **0.179** |
| `sdfp-sd-1` | 336 | 0.066 | **0.072** |
| `6c` | 48 | 0.625 | **1.000** |
| `6d` | 8 | 0.446 | **0.583** |
| `8ver` | 16 | **0.632** | 0.283 |

FK-B's subproblems stay closer to the symmetry of the instance they came from,
on every instance except `8ver` — which is also the one instance in the library
where FK-B has no node-count advantage at all. Higher retained symmetry *and*
fewer nodes is the opposite of what one might expect, and the mechanism is the
point: FK-A shreds the group one vertex at a time and then explores every
resulting orbit member separately, while FK-B's per-term branches remove a whole
clause or monomial in one step and leave a residual instance that is closer to a
scaled copy of the original.

## 4. Why FK-B escapes: a shallower base case

FK-A's base case is `|G| · |H| ≤ 1`, plus the thesis Section 5.1.1 `D_{1,k}`
rule. FK-B's is `min(|C|, |D|) ≤ 2`, decided outright by enumerating the
minimal transversals of the two-term side.

That threshold difference is a whole binary subtree per leaf, and it shows in
where each algorithm's nodes terminate. On `sdfp-sd-3`:

| | leaf reasons |
| --- | --- |
| FK-A | `dual` 13,465 · `single_edge` 6,990 · `trivial` 6,476 |
| FK-B | `easy_case` 1,655 · `equivalent` 22 |

98.7% of FK-B's leaves are decided by the easy case. FK-A must grind each branch
down to a single edge on one side before it can stop, and half its nodes are
internal `dual` confirmations of subtrees that only exist because of that.

## 5. Why FK-B loses on the threshold family

The threshold family inverts the shape: `|G| = 6` against `|H| = 25` at
`v = 10`. FK-A's `single_edge` rule fires almost immediately — `|G| = 1` arises
after a couple of splits, and `D_{1,k}` then decides in `O(k)` — so it finishes
in 9 nodes with **no repeated subproblems at all**.

FK-B has no `D_{1,k}` equivalent. Its μ-test fires `mu_D` six times, and each
`mu_D` branch spawns one subproblem per clause. On an instance that was already
about to be solved, that fan-out is overhead: 27 nodes for the same answer.

The reading is that the two algorithms' cheap paths are aimed at different
things. FK-A's is a *shape* shortcut, for when one side collapses to a single
edge. FK-B's is a *frequency* shortcut, for when no single vertex is frequent
enough to make a two-way split productive. Neither subsumes the other.

## 6. The output-bound exception

`matching(v)` has `|Tr| = 2^(v/2)`: the answer itself is exponential in the
input. FK-A's repeat rate there is high — 90.4% at `v = 12` — but the ratio
stays pinned at about 1.9x across every size tested.

That is the check on the story in §2. A high repeat rate only converts into a
large gap when the repeats are *avoidable*. In the matching family both
algorithms must emit work proportional to an exponential output, so neither can
pull away, and the repeat rate stops being predictive. Symmetry-driven repeats
are avoidable; output-driven ones are not.

## 7. The thesis' `Aut(G) = Aut(H)` claim

Thesis p.51 states that "the input hypergraphs G and H (if they are transversals
of each other) have the same automorphism group". Across **4,034 annotated
nodes** — 2,639 FK-A and 1,395 FK-B, over all fourteen live instances — there
are **zero** nodes where the two sides' groups differ, in either algorithm.
(This section previously totalled 3,455 over twelve instances, carrying the same
mis-stated FK-B figure corrected in §3; the count of disagreements was zero
then and is zero now, over more nodes and two further instances.)

This is stronger than the thesis states it, in two ways: it holds at every node
of the recursion, not only at the root, and it holds for FK-B's decomposition,
which the thesis never examined.

The previously recorded counterexample was `6a` *as published*, where
`Aut(G) = C₂` and `Aut(H) = C₂ × C₂`. That instance is **not a transversal
pair** — which is exactly the hypothesis the claim carries — so it never
contradicted it. It was withdrawn to
`_archive/20260804-verbatim-nontransversal/` on 2026-08-04 and the distinction
is locked by `tests/test_symmetry.py`.

---

## Summary: four regimes

| regime | example | FK-A repeats | who wins | why |
| --- | --- | --- | --- | --- |
| **Symmetry-bound** | `sdfp-sd-*`, `fano`, `f2g2` | 60–98% | **FK-B**, growing with size | FK-A re-derives orbit-equivalent subproblems; FK-B's per-term branch takes a whole orbit representative at once |
| **Primitively symmetric** | `w11`, `w11-sd-1` | 70% | **neither**, ~1.5x | The root group is 4-transitive, so no vertex, pair, triple or quadruple is distinguishable. FK-B's cheap path stops firing: it needs 323 nodes at depth 9 where the Fano plane costs it 11 at depth 3. Added 2026-08-05; see [hard-cases.md](hard-cases.md) |
| **Output-bound** | `matching(v)` | high but unavoidable | neither, flat ~1.9x | `|Tr|` is exponential; both must pay it |
| **Shape-easy** | `threshold(v)`, `8ver` | 0–48% | **FK-A** | `D_{1,k}` settles it before FK-B's per-term fan-out earns anything |

The practical consequence for this project: **`ε(G, H)` and `|Aut|` at the root
predict which algorithm to reach for.** A high root automorphism order with a
self-dual structure is FK-A's worst case and FK-B's best; a small, lopsided pair
is the reverse. Neither dominates, and the instance library now contains
witnesses for both directions.
