# Legacy assessment and migration record

Assessment date: 2026-08-04.

The original notebooks and `fk_a_final.py` were treated as research evidence,
not as an executable specification. The thesis, the archived decomposition log,
the independent transversal oracle, and explicit tests are the maintained
contracts.

## Confirmed defects in the legacy code

1. Superset removal skipped an all-ones row because it searched only zero
   columns. A strict superset could therefore survive preprocessing.
2. The cross-edge intersection precondition attempted the truth value of a
   multi-element NumPy array and raised `ValueError`.
3. The recursion called `remove_duplicate` and `remove_superset`, while only
   `remove_duplicate_rows` and `remove_supersets` existed, causing `NameError`.
4. Base-case control flow returned before the recursive step could be reached
   on relevant paths.
5. Recursive return values were discarded and the parent could report `True`
   even after a child found a non-dual subproblem.
6. Module-level node dictionaries and counters made repeated runs and notebook
   execution order affect results.
7. The visualization regenerated Graphviz images for frames, producing large,
   static artifacts that could not expose the exact node data.

Each algorithmic defect above has a focused regression test under `tests/`.
Randomized oracle comparisons then test the repaired implementation through
both algorithm variants and both pivot rules.

## Preserved evidence

The original source, notebooks, GAP notes, and decomposition log are preserved
under `_archive/20260804-legacy-code/`. They were moved, not deleted or
rewritten. `LEGACY README.md` is the original project notice.

The separate MATLAB FK-B code remains under `FK- B/` as a reference project; it
was not mixed into the FK-A Python package.

## Data corrections

Two thesis instances as printed are not transversal pairs. The canonical data
library keeps both layers:

- `6a-verbatim` and `8ver-verbatim` reproduce the printed inputs;
- `6a` and `8ver` are separately marked `corrected`, with the repair and
  supporting epsilon/transversality evidence recorded in their JSON notes.

This distinction prevents a future cleanup from silently rewriting the source
evidence to fit the expected result.
