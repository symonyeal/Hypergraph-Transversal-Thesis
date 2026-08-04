# Data

`instances/` is the canonical FK-A instance library. Each JSON file contains
the input pair, 1-indexed vertex notation, provenance, a source citation,
research notes, and independently refreshable expected values.

Do not store executable Python, notebook state, or rendered output here.
Generated trees belong in `results/`; programs that derive new data belong in
`experiments/` and should write reviewed JSON here only when explicitly run.
