"""Inspect the local MATLAB FK-B snapshot's benchmark ``.mat`` files.

    python experiments/import_fkb_benchmarks.py

Requires SciPy and the snapshot under ``FK- B/FK-master/``. Neither is a project
dependency: the snapshot is git-ignored and its redistribution terms are not
established, so this program *reports* on it and never copies it into the
repository.

What it found, run on 2026-08-04 against the twelve shipped files:

* Five hold a genuine transversal pair -- ``SDFP16``, ``SDFP23``, ``ms_33``,
  ``sms_33`` and ``BIOMD0000000034``.
* Seven store a ``cnf`` and a ``dnf`` with **different column counts**, so the
  two matrices are not over the same variable set and cannot be a pair at all.
  They appear to have been saved from different runs. They are reported and
  skipped rather than silently coerced.
* The two ``SDFP`` files are the self-dualisation of ``k = 2`` and ``k = 3``
  disjoint Fano planes. That construction is reproduced exactly by
  :func:`fka.instances.self_dual_fano`, which is committed as ``sdfp-sd-1`` and
  ``sdfp-sd-2`` in the instance library. Nothing needs to be imported for them.

So this program is documentation of provenance, not part of any pipeline. If a
pair here is wanted as a research instance, add a generator to
``fka.instances`` and an entry to ``build_instances.py``, as ``sdfp-sd-*`` did.
"""

from __future__ import annotations

import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fka.hypergraph import Hypergraph  # noqa: E402
from fka.sperner import sperner_reduce  # noqa: E402
from fka.transversal import is_dual_oracle  # noqa: E402

SNAPSHOT = (
    Path(__file__).resolve().parents[1]
    / "FK- B"
    / "FK-master"
    / "Modified FK-B Algorithm"
)


def main() -> int:
    try:
        import scipy.io as sio
    except ImportError:
        print("SciPy is not installed; this program only reads MATLAB v5 files.")
        return 2
    if not SNAPSHOT.is_dir():
        print(f"no MATLAB snapshot at {SNAPSHOT}; nothing to inspect.")
        return 2

    usable = 0
    for path in sorted(SNAPSHOT.glob("*_CNF_DNF.mat")):
        data = sio.loadmat(path)
        if "cnf" not in data or "dnf" not in data:
            print(f"{path.name:34} SKIP no cnf/dnf pair")
            continue
        cnf, dnf = data["cnf"], data["dnf"]
        if cnf.shape[1] != dnf.shape[1]:
            print(
                f"{path.name:34} SKIP cnf has {cnf.shape[1]} columns, "
                f"dnf has {dnf.shape[1]}: not the same variable set"
            )
            continue

        G = sperner_reduce(Hypergraph.from_incidence(cnf.tolist())).reduced
        H = sperner_reduce(Hypergraph.from_incidence(dnf.tolist())).reduced
        dual = is_dual_oracle(G, H)
        sizes = dict(sorted(collections.Counter(G.edge_sizes()).items()))
        usable += dual
        print(
            f"{path.name:34} n={G.n:3} |C|={len(G):4} |D|={len(H):4} "
            f"dual={str(dual):5} self-dual={str(G.edges == H.edges):5} "
            f"clause sizes={sizes}"
        )

    print(f"\n{usable} files hold a genuine transversal pair.")
    print("SDFP16/SDFP23 are reproduced by fka.instances.self_dual_fano(2) and (3);")
    print("they are committed as sdfp-sd-1 and sdfp-sd-2 and need no import.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
