# Glossary: one object, four vocabularies

FK-A is stated in the language of hypergraphs and transversals. FK-B is stated
in the language of monotone Boolean functions and normal forms. **They are the
same objects.** A monotone Boolean function *is* a hypergraph, and its dual *is*
the transversal hypergraph. Everything in this repository is built on that
identification, which is why one `Hypergraph` type, one `RecursionTree`, one
instance library and one visualiser serve both algorithms.

## The correspondence

| hypergraph | monotone Boolean function | normal form | this repo |
| --- | --- | --- | --- |
| hypergraph `H` on `V` | monotone Boolean function `f` on `V` | — | `Hypergraph` |
| hyperedge `e` | minimal true point of `f` | monomial of the DNF | an `int` bitset |
| vertex set `S` | assignment setting `S` true | assignment | an `int` bitset |
| `S` contains an edge | `f(S) = 1` | DNF satisfied | `fkb.dnf_value` |
| `S` meets every edge | `S` is a transversal | CNF satisfied | `fkb.cnf_value`, `is_transversal` |
| **Sperner** (antichain) | irredundant | irredundant | `sperner_reduce`, FK-B's `Irredundant` |
| `Tr(H)`, minimal transversals | the **dual** `f^d` | the CNF of the same `f` | `transversal` |
| `G = Tr(H)` | `f^d = g`, i.e. `C ≡ D` | CNF and DNF agree | what both algorithms decide |

So the single question both algorithms answer can be written three ways, and
they mean the same thing:

```
G = Tr(H)          hypergraph dualization
f^d = g            monotone Boolean duality
C ≡ D              equivalence of a CNF and a DNF
```

## Which name each algorithm uses

| | FK-A | FK-B |
| --- | --- | --- |
| the two sides | `G`, `H` | `C` (clauses), `D` (monomials) |
| in this code | `G`, `H` | `cnf`, `dnf` — the same two `Hypergraph`s |
| a counterexample | **certificate** `S`, thesis eq. 2.1 | **conflicting assignment** `S`, `C(S) ≠ D(S)` |
| answer when they match | `dual = True` | `dual = True` (reported as `equivalent`) |

`G` is the CNF side and `H` is the DNF side. `fk_a(G, H)` and `fk_b(G, H)` take
their arguments in the same order and decide the same predicate — that is what
`python -m fkb compare --all` checks, and the test suite asserts they can never
disagree.

### The one asymmetry worth knowing

A certificate and a conflicting assignment are *nearly* the same thing, not
quite. Once FK-B's condition 1 holds — every clause meets every monomial — we
have `D ≤ C` pointwise, so the only possible conflict is `C(S) = 1, D(S) = 0`,
which is exactly equation 2.1: `S` meets every edge of `G` and contains no edge
of `H`.

Condition 1 *failing* is the one case that produces the other polarity,
`C(S) = 0, D(S) = 1`. So every FK-B witness is an FK-A certificate except those
from `cond_1_disjoint`. `fkb.is_conflict` accepts both;
`fka.verify_certificate` accepts only the first. Neither is wrong — they are
checking different, correctly-stated claims.

## Degenerate cases

Read as a CNF and as a DNF, the two empty-ish hypergraphs swap meaning. This
trips people up, so it is stated once here and relied on everywhere:

| hypergraph | as a CNF (`G` side) | as a DNF (`H` side) |
| --- | --- | --- |
| `{}` — no edges | constant **1** (nothing to satisfy) | constant **0** (nothing to contain) |
| `{{}}` — one empty edge | constant **0** (unsatisfiable clause) | constant **1** (empty monomial holds) |

Hence `Tr({}) = {{}}` and `Tr({{}}) = {}`, and the two are dual — which is the
convention FK96 uses and the one `fka.transversal` produces by construction
rather than by special case.

## How big is the space

`d(n)`, the *n*-th Dedekind number, counts the antichains on an `n`-set — so it
counts the Sperner hypergraphs, the monotone Boolean functions, and the
irredundant normal forms on `n` variables, all at once. That is the same
identification again, and it is why the vendored `Dedekind/` snapshot belongs in
this repository. At `n = 8`, where the thesis' largest instances sit, it is
5.6 × 10²².
