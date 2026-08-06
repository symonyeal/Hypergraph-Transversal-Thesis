# FK-A against FK-B

Regenerate with `python experiments/compare_algorithms.py`. Every row is
checked against the brute-force transversal oracle before it is reported.

`retained` is the mean `|Aut(C)|` over the tree divided by the root's: how
much of the instance's symmetry the average subproblem still carries.

| instance | n | nC | nD | dual | root Aut | FK-A nodes | FK-B nodes | ratio | FK-A retained | FK-B retained |
| --- | ---: | ---: | ---: | :-: | --- | ---: | ---: | ---: | ---: | ---: |
| `6a` | 6 | 4 | 4 | yes | C2 x C2 (4) | 11 | 3 | 3.667x | 0.5714 | 0.6667 |
| `6b` | 6 | 8 | 3 | yes | C2 x S4 (48) | 13 | 6 | 2.167x | 0.175 | 0.2083 |
| `6c` | 6 | 8 | 2 | yes | C2 x S4 (48) | 5 | 1 | 5.0x | 0.625 | 1.0 |
| `6d` | 6 | 4 | 4 | yes | D4 (8) | 9 | 3 | 3.0x | 0.4464 | 0.5833 |
| `7ver` | 7 | 11 | 7 | yes | C2 (2) | 21 | 11 | 1.909x | 1.3333 | 1.3636 |
| `8ver` | 8 | 5 | 20 | yes | C2 x D4 (16) | 21 | 21 | 1.0x | 0.6324 | 0.2827 |
| `f2g2` | 8 | 8 | 4 | yes | D4 wr C2 (128) | 25 | 7 | 3.571x | 0.1016 | 0.1786 |
| `fano` | 7 | 7 | 7 | yes | PSL(3,2) (168) | 35 | 11 | 3.182x | 0.0648 | 0.1342 |
| `k3` | 3 | 3 | 3 | yes | S3 (6) | 3 | 3 | 1.0x | 0.5556 | 0.5556 |
| `sdfp-sd-1` | 9 | 15 | 15 | yes | order 336 (non-abelian, unidentified) (336) | 75 | 37 | 2.027x | 0.0658 | 0.0718 |
| `sdfp-sd-2` | 16 | 64 | 64 | yes | — (—) | 1427 | 231 | 6.177x | — | — |
| `trivial-aut-1` | 8 | 31 | 35 | yes | C1 (1) | 41 | 33 | 1.242x | 4.8462 | 5.3636 |
| `witt-s4511` | 11 | 66 | 66 | yes | — (—) | 483 | 321 | 1.505x | — | — |
| `word-xy` ‡ | 16 | 32 | 16 | yes | — (—) | 577 | 193 | 2.99x | — | — |
| `word-yx` ‡ | 16 | 16 | 32 | yes | — (—) | 577 | 223 | 2.587x | — | — |
| `sdfp-sd-3` ‡ | 23 | 365 | 365 | yes | — (—) | 26931 | 1677 | 16.059x | — | — |

‡ too large to annotate; node counts only. See `ANNOTATE_BELOW`.

## Automorphism groups seen inside each tree

Group of the `C` side, counted over every node whose group was
computed. A node above `skip_group_above` edges is not counted.

| instance | FK-A groups | FK-B groups |
| --- | --- | --- |
| `6a` | `C2`×6, `C2 x C2`×1 | `C2`×2, `C2 x C2`×1 |
| `6b` | `C2`×6, `D4`×3, `C2 x S4`×1 | `C1`×4, `C2 x S4`×1, `D4`×1 |
| `6c` | `S4`×3, `C2 x S4`×1 | `C2 x S4`×1 |
| `6d` | `C2`×3, `D4`×1, `C2 x C2`×1, `C1`×1 | `D4`×1, `C2`×1, `C2 x C2`×1 |
| `7ver` | `C2`×10, `C1`×4, `S3`×2, `C2 x C2`×1 | `C2`×6, `C1`×2, `C2 x C2`×2, `D4`×1 |
| `8ver` | `C2`×6, `S3`×4, `D4`×3, `S4`×2 | `C1`×7, `C2`×5, `S3`×3, `D4`×2 |
| `f2g2` | `C2`×9, `D4`×5, `C2 x D4`×3, `D4 wr C2`×1 | `C2`×4, `D4 wr C2`×1, `C2 x D4`×1, `D4`×1 |
| `fano` | `C2`×14, `S3`×4, `S4`×2, `D4`×2 | `C2`×4, `S4`×2, `D4`×2, `C2 x C2`×2 |
| `k3` | `C2`×2, `S3`×1 | `C2`×2, `S3`×1 |
| `sdfp-sd-1` | `C2`×28, `S3`×8, `PSL(3,2)`×4, `S4`×4 | `C2`×16, `C2 x C2`×8, `S4`×4, `PSL(3,2)`×2 |
| `sdfp-sd-2` | `C2`×266, `S4`×212, `C1`×186, `S3`×162 | `D4`×101, `S3`×78, `C2`×14, `C2 x D4`×7 |
| `trivial-aut-1` | `C2`×22, `S3`×6, `C1`×5, `S4`×3 | `C2`×15, `C1`×5, `S4`×3, `C2 x C2`×3 |
| `witt-s4511` | `C2`×202, `S3`×50, `D4`×48, `S4`×36 | `C2`×120, `S3`×74, `C2 x C2`×40, `C1`×30 |
| `word-xy` | not annotated | not annotated |
| `word-yx` | not annotated | not annotated |
| `sdfp-sd-3` | not annotated | not annotated |

## FK-B split branches

`mu_D` and `mu_C` are the per-clause and per-monomial branches, taken
when the split variable is at most mu-frequent on one side. `split` is
the plain two-way split, which is FK-A's step.

| instance | branches |
| --- | --- |
| `6a` | `split`×1 |
| `6b` | `mu_D`×1 |
| `6c` | — (no split) |
| `6d` | `split`×1 |
| `7ver` | `split`×5 |
| `8ver` | `mu_C`×1, `mu_D`×1, `split`×5 |
| `f2g2` | `mu_C`×2 |
| `fano` | `split`×5 |
| `k3` | `split`×1 |
| `sdfp-sd-1` | `mu_C`×1, `mu_D`×5, `split`×11 |
| `sdfp-sd-2` | `mu_C`×8, `mu_D`×18, `split`×7 |
| `trivial-aut-1` | `split`×16 |
| `witt-s4511` | `mu_C`×11, `mu_D`×11, `split`×127 |
| `word-xy` | `mu_C`×3, `mu_D`×28, `split`×1 |
| `word-yx` | `mu_C`×32, `mu_D`×5, `split`×1 |
| `sdfp-sd-3` | `mu_C`×10, `mu_D`×11, `split`×1 |
