# Dictionary: hypergraph ↔ monotone Boolean function

One object, two vocabularies. FK-A is stated over hypergraphs and transversals,
FK-B over monotone Boolean functions and normal forms. A monotone Boolean
function **is** a hypergraph and its dual **is** the transversal hypergraph, so
one model serves both.

**The code speaks MBF.** Names come from the FK-B reference in
`FK- B/FK-master/src` and from *How to apply SAT-solving for the equivalence
test of monotone normal forms*, p.5. This file is the only place the hypergraph
reading is written down; nothing in `src/` carries a second name for a concept
that already has one here.

Rules that keep it that way:

1. One concept, one identifier. If a name appears in the reference, that name
   wins, transliterated to Python casing (`Check_Conditions` → `check_conditions`).
2. Thesis notation survives only where the reference has no counterpart, and
   only in the "hypergraph" column below.
3. Nothing is added to this table without approval. Appended entries go under
   [Change log](#change-log) with a date, never silently edited in place.

## The correspondence

| concept | hypergraph | MBF / normal form | code | on disk |
| --- | --- | --- | --- | --- |
| the object | hypergraph `H` on `V` | monotone Boolean function | `Hypergraph` | — |
| a value of it, role unfixed | hypergraph | `MBF` | `MBF` | — |
| the two sides | `G`, `H` | `C` (clauses), `D` (monomials) | `cnf`, `dnf` | `cnf`, `dnf` |
| a term | hyperedge `e` | clause / monomial | `int` bitset | variable list |
| all of them | `E(H)` | the term set | `.edges` | the `cnf`/`dnf` array |
| how many | `\|E(G)\|`, `\|E(H)\|` | `nC`, `nD` | `nC`, `nD` | `nC`, `nD` |
| a variable set | `S ⊆ V` | assignment setting `S` true | `S` | variable list |
| `f(S) = 1` | `S` contains an edge | DNF satisfied | `dnf_value` | — |
| `S` is a transversal | `S` meets every edge | CNF satisfied | `cnf_value` | — |
| antichain reduction | Sperner reduction | `Irredundant` | `irredundant` | — |
| the dual | `Tr(H)`, minimal transversals | `f^d`, the other normal form | `transversal` | — |
| the question | `G = Tr(H)` | `C ≡ D` | `fk_a`, `fk_b` | `dual` |
| a counterexample | certificate `S`, eq. 2.1 | conflicting assignment | `CA` | `CA` |
| the branch vertex | pivot | split variable | `split_var`, locally `x` | `split_var` |
| how it is chosen | pivot rule | `Choose_SplitVar` method | `split_rule` | `split_rule` |
| terms without `x` | deletion `H \ x` | `phi_x_1` | `phi_x_1` | — |
| terms with `x`, `x` dropped | contraction `H / x` | `phi_x_0` | `phi_x_0` | — |
| force a clause false | — | `A_c_x` | `A_c_x` | — |
| force a monomial true | — | `A_m_x` | `A_m_x` | — |
| frequency of `v` | `ε_v(G, H)` | — | `eps_v` | — |
| frequency of the pair | `ε(G, H)` | — | `eps` | `eps` |
| the µ bound | — | `mu_function` | `mu` | — |
| rare on one side | — | `mu_frequent_in_A` | `mu_frequent` | — |
| terms removed as redundant | non-Sperner edges | — | `cut_C`, `cut_D` | `cut_C`, `cut_D` |

`cnf` is the `G` side and `dnf` the `H` side. `fk_a(cnf, dnf)` and
`fk_b(cnf, dnf)` take their arguments in the same order and decide the same
predicate.

### The recursion tree

Neither source names these; they are the project's own, kept short and used the
same way by both algorithms.

| concept | code | on disk |
| --- | --- | --- |
| the whole trace | `Tree` | the file |
| one recursive call | `Node` | an entry of `nodes` |
| inputs before reduction | `cnf_in`, `dnf_in` | `cnf_in`, `dnf_in` |
| inputs after `Irredundant` | `cnf`, `dnf` | `cnf`, `dnf` |
| which subproblem this is | `step` — `L`/`R`, `x=0`/`x=1`, `c1..`, `m1..` | `path` |
| which split shape was taken | `branch` — `mu_D`, `mu_C`, `split` | `branch` |
| why the node concluded | `Verdict.reason` | `verdict.reason` |
| every assignment found here | `CAs` | `CAs` |

## Where the two readings genuinely differ

A certificate and a conflicting assignment are *nearly* the same thing. Once
FK-B's condition 1 holds — every clause meets every monomial — `D ≤ C`
pointwise, so the only possible conflict is `C(S) = 1, D(S) = 0`, which is
exactly thesis equation 2.1: `S` meets every edge of `G` and contains no edge of
`H`.

Condition 1 *failing* is the one case producing the other polarity,
`C(S) = 0, D(S) = 1`. Every FK-B witness is an FK-A certificate except those
from `cond_1_disjoint`. `is_conflict` accepts both; `verify_CA` accepts only the
first. Both are correct statements of different claims.

## Degenerate cases

Read as a CNF and as a DNF the two empty-ish families swap meaning:

| family | as `cnf` | as `dnf` |
| --- | --- | --- |
| `{}` — no terms | constant **1** (nothing to satisfy) | constant **0** (nothing to contain) |
| `{{}}` — one empty term | constant **0** (unsatisfiable clause) | constant **1** (empty monomial holds) |

Hence `Tr({}) = {{}}` and `Tr({{}}) = {}`, which is FK96's convention and what
`transversal` produces by construction rather than by special case.

## Names that are hypergraph-only

These have no normal-form counterpart and keep thesis notation. `MBF` is still
the parameter name; only the concept is hypergraph-side.

`primal_graph`, `automorphism_group` / `Aut`, orbit, block system, primitive,
conformal, α-acyclic, read-once, cograph, `support`, `degree`, isolated vertex,
`Tr`, Berge's method.

## Where a bare `G` or `H` is still correct

Three uses are *not* the hypergraph `G`/`H` and are left alone:

- **a group.** `|G|`, `[G, G]` in `fka.groups` — standard group notation.
- **a graph.** `G`, `G^c` in the cograph criterion — standard graph notation.
- **a verbatim quotation** of the thesis or a cited paper, which is reproduced
  as written and never silently renamed.

Anywhere else, a bare `G` or `H` in this repository is a leftover.

## Change log

Entries are appended, with a date, only after approval.

- **2026-08-06** — created. Adopted the MBF register repo-wide, replacing the
  mixed `G`/`H` + `cnf`/`dnf` vocabulary and the split `pivot` / `split_var`,
  `certificate` / `CA` pairs. Supersedes `docs/glossary.md`, archived under
  `_archive/20260806-mbf-vocabulary/`.
  Renamed to match the reference: `c_set_0`/`c_set_1`/`d_set_0`/`d_set_1` →
  `A_c_x`/`A_m_x`, `sperner_reduce` → `irredundant`, `pivot_rule` →
  `split_rule`, `fk_b(irredundant=)` → `fk_b(redundancy_check=)` after
  `FK_B.m`'s own `redundancychk`. Renamed for clarity, with no counterpart in
  either source: `RecursionTree`/`RecursionNode` → `Tree`/`Node`,
  `SpernerReduction` → `Reduction`, `HypergraphProperties` → `Properties`,
  `fka.threshold` → `fka.bound` (it collided with the `threshold` generator),
  `substitution.union`/`join` → `disjoint_vee`/`disjoint_wedge` (they are the
  variable-disjoint operators, unlike `Hypergraph.vee`). The on-disk JSON schema
  moved with the code.

  Deleted as duplicates of something already here: the second and third copies
  of `eps` (`fka.threshold.eps`, `fka.instances.exact_epsilon`, now the single
  `fka.algorithm.eps_exact`), the second `FANO_LINES` (`fka.selfdual` owns it),
  the `is_normal` alias of `is_conformal`, `remove_superset_pairwise` (the
  definition it existed to be compared against now lives in the test that makes
  the comparison), and `Hypergraph.contraction`/`.deletion` (they were
  `phi_x_0`/`phi_x_1` under another name).

- **2026-08-06** — group theory. `fka.selfdual`'s incidence-graph automorphism
  search was removed; `fka.automorphism.term_set_automorphisms` is now the one
  search, taking one term set (`Aut(MBF)`) or two (the pair group). Group naming
  gained `conjugate`, `commutator` and per-factor spelling, so the same group
  reads the same standalone and inside a product -- which **renames**
  `D3 wr C2` → `S3 wr C2`, `D3 x D4` → `S3 x D4`, `D3 x S4` → `S3 x S4` in
  published results. `_quaternion8` was rebuilt: the constants it used generate
  C2 x C4, not Q8, so Q8 was unnameable and C2 x C4 was misnamed `Q8`.
