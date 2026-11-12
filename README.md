# Reputation memory and spatial cooperation: research data

Data and code accompanying *The Influence of Reputation Memory on Cooperation
in a Spatial Prisoner's Dilemma* (Physics submission physics-4642818).

Repository: [labdinpop/reputation-memory](https://github.com/labdinpop/reputation-memory).

Authors: Pedro Henrique Triguis Schimit, Pedro Aranha, Marco A. R. Fontoura,
Wellington de Lima Barreto, and Abimael R. Sergio.

## Start here

- [Residual drift table](residual_drift.csv): the 14 parameter combinations
  discussed in the reply to Reviewer 2, comment 2. This is a supplementary
  data table derived from saved simulation results. The reply itself displays
  the table; this CSV provides an electronic copy with additional metadata.
- [Column definitions](DATA_DICTIONARY.md): meaning of each column and the
  criterion used to identify residual drift.
- [Configuration files](config/): simulation parameters.
- [Results](results/): individual realizations, summaries, and saved outputs.
- [Code](src/): model implementation, simulation drivers, and analysis.
- [Validation records](qa/): checks supplied with the data release.
- [Reproduction instructions](README_R2.txt): how the added numerical values
  and supplementary comparison were generated from completed runs.

The `source` column in the drift table is a path relative to this directory.
It links each row to the saved realizations used to calculate its count.
It is a provenance field, not a bibliographic reference. All 14 paths resolve
within this data collection.

## Interpreting the results

Each row in `residual_drift.csv` represents one distinct parameter combination
after the selected extension of the transient. The 14 rows are not 14
realizations. `R` gives the number of independent realizations, and
`drifting_realizations` gives the number that exceeded the drift criterion.

Some result directories contain successive transient durations. Use the
`source` fields in the tables to identify the selected output; do not combine
earlier and extended runs as independent replicates of the final estimate.
Empty cells in the drift table identify parameters that do not apply to the
memory rule in that row. They do not represent zero.

The lattice side length is `L`; the population contains `L * L` sites.
Absorption fractions describe the finite populations and simulation durations
reported in the manuscript.

## Reproducing the supplementary table

The script `submissão/03 Physics R2/revision_work/build_revision_assets.py`
reads saved results and regenerates the residual drift table, added numerical
macros, and revised Figure S3. See `README_R2.txt` for its directory layout and
dependencies. It does not launch new simulations. Simulation drivers named
`run_*.py` do launch simulations and are not needed to inspect these data.

`revision_numbers_sources.json` records the origin of the numerical additions.
The complete manuscript sources are distributed separately in
`manuscript_sources.zip`.

## File integrity

This release retains all files from the previous data archive and includes
the completed auxiliary comparison and supporting code used in revision 2.
`SHA256SUMS.txt` lists the SHA-256 digest of every other file in the release.
No new simulations were run to prepare the editorial revision.
