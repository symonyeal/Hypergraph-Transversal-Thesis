# Dedekind numbers — third-party reference implementation

`DedekindNumbersShiftFourMatrixMultiplication-main/` is a verbatim snapshot of
[Christian Jäkel's](https://github.com/cjaekel) implementation accompanying
*A computation of the ninth Dedekind Number*
([arXiv:2304.00895](https://arxiv.org/abs/2304.00895), 2023).

Unlike the MATLAB FK-B snapshot in `FK- B/`, this **is** committed: it carries an
explicit MIT licence (`LICENSE.txt`, Copyright (c) 2023 Christian Jäkel), which
permits redistribution provided the copyright notice travels with it. That file
is included unmodified and is the notice.

The snapshot is unmodified apart from collapsing the doubled
`…-main/…-main/` directory that GitHub's zip export produces. Nothing in `src/`,
the instance library, or either test suite reads it.

## Why it is in this repository

`d(n)`, the *n*-th Dedekind number, counts the antichains on an `n`-element set
— which is exactly the number of Sperner hypergraphs on `n` vertices, and so the
size of the space FK-A and FK-B search:

| n | d(n) | |
| ---: | ---: | --- |
| 4 | 168 | |
| 5 | 7,581 | |
| 6 | 7,828,354 | |
| 7 | 2,414,682,040,998 | |
| 8 | 56,130,437,228,687,557,907,788 | the thesis' largest instances live at n = 8 |
| 9 | 286386577668298411128469151667598498812366 | Jäkel 2023, and independently Van Hirtum et al. |

That is the honest scale of the dualization problem. It is also why the
instance library is a curated dozen rather than a sample: at `n = 8` there is no
meaningful notion of sampling the space, and every instance has to earn its
place by being a published or constructed object.

`Reading/Yusun Dedekind.pdf` is the companion reading. See
`docs/fk-a-vs-fk-b.md` for what the two algorithms actually do inside that
space.

## Licensing

This snapshot is MIT and is **not** covered by the repository's
`GPL-3.0-or-later`. The two licences coexist here because MIT is
GPL-compatible and the code is vendored, not linked: nothing in `src/` imports
or compiles against it.

Building it needs Eigen, Intel TBB, and CUDA for the GPU path — see its own
`README.md`. None of that is a dependency of this project.
