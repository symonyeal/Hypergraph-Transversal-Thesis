# A new hard case: the Witt design S(4,5,11)

The thesis' two hard classes are the Fano plane (Section 5.2 — self-transversal,
`Aut = PSL(3,2)`, frequency `3/7`) and the alternating trees of Gurvich–
Khachiyan, Takata and Hertel (Section 4.2 — `Tr(f_k) = g_k`, frequency exactly
at the Fredman–Khachiyan floor). This note reports a third: the **Witt design
`S(4,5,11)`**, 66 blocks of size 5 on 11 points, `Aut = M11`. It is in the
library as `w11`, and its self-dualisation as `w11-sd-1`.

Everything below is regenerable:

```text
sage -python experiments/hard_case_search.py 3 15   -> results/hard-case-search.json
python experiments/build_instances.py               -> data/instances/w11*.json
python -m fkb compare w11                           -> both algorithms, one instance
pytest tests/test_hard_cases.py
```

---

## 1. What was searched, and why that class

A hypergraph is self-transversal exactly when it is a **maximal intersecting
family**: for every `S`, exactly one of `S`, `V\S` contains an edge. That is the
class with no excuses. The answer is always "dual", so every branch of the
recursion has to be closed; and `|Tr(H)| = |H|`, so no part of the tree can be
blamed on an exponential output the way `matching(v)`'s can.

Three facts turn the class into a finite search rather than a hunt.

**L1.** Fix an edge `e` of minimum size `s`. Every other edge meets `e`, so with
`D` the maximum degree, `m ≤ 1 + s(D − 1)`, hence

```text
ε(H) = D/m > 1/s.
```

**L2.** FK96 guarantees `ε ≥ 1/log₂(|G| + |H|)` for a dual pair, and thesis
Section 4.3 calls a family hard when `ε` sits on that floor. Writing

```text
T(H) = ε · log₂(2m)     ≥ 1, and = 1 exactly at the floor,
```

L1 gives `T > log₂(2m)/s`. **A tight instance needs edges of size about the log
of their own number** — which is exactly what the alternating trees do
(`|e| = 2^(k−1)`, `|E| = 2^(2^k − 1)`) and exactly why they attain the bound.

**L3.** If `H` is `s`-uniform and vertex-transitive then `D = ms/n`, and L1
collapses to `n ≤ s² − 1`. So each ground set admits only a few edge sizes, and
every candidate is a union of orbits of a transitive group on subsets.

`experiments/hard_case_search.py` enumerates that exactly: every transitive
group of degree ≤ 15 in SageMath's database, every orbit of `s`-subsets allowed
by L3, every intersecting pair of orbits, filtered by two cheap necessary
conditions and then decided exactly — over all `2ⁿ` subsets at once, never by
building `Tr(H)`.

## 2. What is out there

21 vertex-transitive self-transversal hypergraphs up to isomorphism on ≤ 15
points, sorted by FK-A nodes per unit of input:

| n | \|E\| | sizes | ε | T | FK-A | FK-B | A/2m | B/2m | primitive | Aut |
| ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | :---: | --- |
| 15 | 300 | 6 | 0.400 | 3.69 | 3687 | 2695 | 6.14 | 4.49 | no | (A5×A5×A5) : (C2×S4) |
| 15 | 270 | 6 | 0.400 | 3.63 | 3199 | 1659 | 5.92 | 3.07 | no | (C3⁵) : (C2×((C2⁴):S5)) |
| **11** | **66** | **5** | **5/11** | **3.20** | **483** | **321** | **3.66** | **2.43** | **yes** | **M11** |
| 7 | 7 | 3 | 3/7 | **1.63** | 35 | 11 | 2.50 | 0.79 | yes | PSL(3,2) — *the Fano plane* |
| 9 | 27 | 4 | 4/9 | 2.56 | 119 | 57 | 2.20 | 1.06 | no | S3 wr S3 |
| 11 | 121 | 5,6 | 0.496 | 3.93 | 479 | 317 | 1.98 | 1.31 | yes | PSL(2,11) |
| 10 | 51 | 4,5 | 0.471 | 3.14 | 179 | 109 | 1.75 | 1.07 | yes | A6 |
| 11 | 110 | 5 | 5/11 | 3.54 | 271 | 163 | 1.23 | 0.74 | yes | C11 : C10 |
| 6 | 10 | 3 | 1/2 | 2.16 | 19 | 11 | 0.95 | 0.55 | yes | A5 |
| 3, 5, 7, 9, 11, 13, 15 | C(n, (n+1)/2) | | ≈ 1/2 | 1.7–7.3 | | | ≤ 0.56 | ≤ 0.58 | yes | S_n — *majority* |

Two readings.

**The class splits by primitivity.** Every instance above the Fano plane other
than W11 has an *imprimitive* automorphism group, and each is a composition:
the 9-point one is `MAJ3[MAJ3]`, the two 15-point ones are `MAJ5[MAJ3]` and
`MAJ3[MAJ5]` (verified by isomorphism against the constructions). Composition
multiplies frequencies — `ε(f[g]) = ε(f)·ε(g)` — which is how they get to
`ε = 0.4`, and it is the same reason the alternating trees are hard: **FK does
not notice that the instance decomposes**, so it re-derives the same module
over and over. A decomposition-aware dualiser cracks all of them, and the
thesis says as much about the alternating trees, which are read-once
(Section 4.4).

**Among the primitive ones, W11 is the hardest, and the Fano plane is second.**
Those two are also the only primitive non-majority members with a 2-transitive
group. They are not a coincidence: they are the same construction on two
different perfect codes.

## 3. Fano and W11 are the same object one step apart

The blocks of both are the supports of the minimum-weight words of a perfect
code:

| code | design | n | \|E\| | \|e\| | ε | Aut |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Hamming [7,4,3] | S(2,3,7) — Fano | 7 | 7 | 3 | 3/7 | PSL(3,2), 168 |
| ternary Golay [11,6,5] | **S(4,5,11) — W11** | 11 | 66 | 5 | 5/11 | **M11, 7920** |
| binary Golay [23,12,7] | S(4,7,23) | 23 | 253 | 7 | 3/11 | M23, 10200960 |

`fka.instances.witt_11` builds W11 that way, from the cyclic generator
polynomial `g(x) = −1 + x² − x³ + x⁴ + x⁵` over GF(3), so nothing is a
transcribed block list.

**The family stops at 11.** S(4,7,23) is intersecting — its blocks meet in 1 or
3 points — but it is *not* self-transversal: 16744 of the 1352078 11-subsets of
its 23 points contain no block and have a complement containing no block, and
each is a thesis equation 2.1 certificate against it. One of them,
`{1,2,3,4,5,6,7,8,10,14,21}`, is asserted in `tests/test_hard_cases.py`. The
same test also shows why the search cannot be run at that size naively: the
certificate settles the question in `O(|E|)`, while building `Tr` of 253 blocks
on 23 points is not a bounded computation.

Note also that this is where the analogy to *projective planes* ends, and the
thesis' framing of the Fano plane as "the smallest projective plane" slightly
undersells it. PG(2,q) for `q ≥ 3` is **not** self-transversal — it has
non-trivial blocking sets, i.e. blocking sets containing no line, and such a set
and its complement are a certificate. The Fano plane is the only projective
plane in the class. What generalises it is not the plane, it is the perfect code.

## 4. What W11 costs the two algorithms

On the library instance `w11` (66 blocks, ternary-Golay labelling), against the
Fano plane:

| | Fano | W11 | |
| --- | ---: | ---: | --- |
| FK-A `faithful` | 59 nodes, 4.21/2m | 1015, **7.69/2m** | depth 6 → 10 |
| FK-A `modified` | 35 nodes, 2.50/2m | 475, **3.60/2m** | depth 5 → 9 |
| FK-B `faithful` | 11 nodes, 0.79/2m | 323, **2.45/2m** | depth 3 → 9 |
| FK-A repeat rate (`modified`) | 60.0% | 70.5% | |
| FK-B repeat rate | **0.0%** | **33.4%** | |
| FK-A / FK-B | 3.18× | 1.47× | |

The FK-B column is the interesting one. `docs/fk-a-vs-fk-b.md` puts the Fano
plane in the "symmetry-bound" regime where **FK-B wins comfortably** — it
decides the Fano plane in 11 nodes, at depth 3, with every leaf an `easy_case`
and not one repeated subproblem. W11 takes that away: FK-B needs 323 nodes at
depth 9, a third of them re-derivations of subproblems it has already solved.

That makes W11 the first instance in the library that is hard for *both*
algorithms at once. The regime table in `docs/fk-a-vs-fk-b.md` gains a row:

| regime | example | who wins | why |
| --- | --- | --- | --- |
| **Primitively symmetric** | `w11` | neither, ~1.5× | The group is 4-transitive, so no vertex, pair, triple or quadruple is distinguishable. FK-A's single-vertex split and FK-B's per-term branch are both choosing between indistinguishable options, and neither algorithm's cheap path fires. |

The mechanism is the one the thesis proposes in Chapter 3, at its extreme.
`M11` is 4-transitive on 11 points; the design is a 4-design. Any split FK-A
makes is, up to an automorphism, *the same split*, and any four vertices look
alike to it. The Fano plane's `PSL(3,2)` is only 2-transitive, so a Fano
subproblem starts distinguishing vertices one level sooner — which is what the
smaller repeat rate and the FK-B collapse to 11 nodes are measuring.

`w11-sd-1` is the same instance put through SDFP's self-dualisation (133 edges
on 13 vertices), for direct comparison with `sdfp-sd-1` (15 edges on 9). It
holds FK-A's rate (3.59/2m against 4.10) and roughly doubles FK-B's (2.70/2m
against 1.23). `k = 2` is not in the library: it would have `66² + 133 = 4489`
edges.

## 5. Frequency tightness and tree size are different axes

Worth stating plainly, because the search makes it measurable. The Fano plane is
the **tightest** instance in the library on the thesis' own criterion —
`T = 1.632`, closer to the FK floor than anything else including the alternating
family member `f2g2` (`T = 1.79`) and the SDFP benchmarks (`T = 2.6` and `5.5`).
W11 is looser (`T = 3.20`) and yet costs both algorithms more per edge. And
`sdfp-sd-2`, the loosest instance in the library at `T = 5.47`, has the largest
tree of all relative to its input (22.65 nodes/2m).

So `ε` bounds the *worst case per split* and predicts nothing about the tree on
its own; what predicts the tree is how much symmetry survives the split.
`ε(G,H)` and `|Aut|` at the root together are the pair of numbers to look at —
which is what `docs/fk-a-vs-fk-b.md` §3 already concluded from the other
direction, and W11 is the instance that separates the two axes cleanly.

## 6. Limits of the search

- Vertex-transitive only. A self-transversal hypergraph with a non-transitive
  group is not reachable this way; L3 does not apply to one, so it would need a
  different enumeration.
- Degrees ≤ 15, and pairs of orbits only up to degree 12. Unions of three or
  more orbits are not enumerated at any degree.
- Node counts depend on the vertex labelling, because the labelling moves the
  pivot choice. The search reports 483/321 for its own labelling of `S(4,5,11)`
  and the library instance reports 475/323 for the ternary-Golay labelling; they
  are the same design. `docs/fk-a-vs-fk-b.md` records the same caveat for the
  SDFP instances against the MATLAB reference.
