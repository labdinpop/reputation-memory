# Using the code and saved results

This repository accompanies *The Influence of Reputation Memory on Cooperation
in a Spatial Prisoner's Dilemma*. It contains simulation outputs, model
implementations, and scripts for calculating the quantities reported in the article.

## Requirements and layout

The scripts use Python with NumPy, Numba, Matplotlib, and Pillow, together with
the Python standard library. Keep `src/`, `config/`, `results/`, and `qa/` under
the repository root so that relative paths resolve correctly.

The data can be inspected directly as CSV and JSON files. NumPy arrays (`.npy`)
store saved lattice configurations.

## Reading existing results

- `src/model.py` implements the main model with a mixture of uniform and
  reputation-guided reference selection.
- `src/model_s3.py` implements proportional reference selection.
- `src/model_s3_past.py` includes the memory rule that uses only past actions.
- `src/analyze_e1.py` and `src/analyze_phase2.py` aggregate saved realizations
  and generate numerical summaries and figures.
- `src/analyze_s3_past.py` analyzes the comparison using only past actions.

For a parameter point with more than one transient duration, use the selected
output recorded in its metadata or the `source` column of its summary. Earlier
and extended runs should not be combined as additional independent realizations.

## Regenerating the supplementary comparison

From the repository root, run:

```bash
python reproduce_supplement.py
```

This reads completed simulations and generates the residual drift table,
numerical macros, and Figure S3. The command reports the generated counts and
uses the relative paths retained in the analysis script. It does not run new
simulations. `revision_numbers_sources.json` records the sources of these values.

The comparison using only past actions has `a=0`, `w=1`, `b=1.15`,
`epsilon=0`, `eta=1e-6`, `L=100`, and 30 realizations per memory length.
Its settings differ from the main comparisons with `w=0.5`.

## Running simulations

The `run_*.py` scripts launch simulations using the corresponding files in
`config/`. Inspect their parameter grids and output paths before starting a run.
For inspection and reanalysis of the published data, use the saved outputs.

## File integrity

`SHA256SUMS.txt` contains the SHA-256 digest of every other file distributed
with the dataset. It can be used to check that downloaded files are intact.
