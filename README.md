# Reputation memory and spatial cooperation

Code and data accompanying *The Influence of Reputation Memory on Cooperation
in a Spatial Prisoner's Dilemma*.

Authors: Pedro Henrique Triguis Schimit, Pedro Aranha, Marco A. R. Fontoura,
Wellington de Lima Barreto, and Abimael R. Sergio.

## What the model studies

Individuals occupy a square lattice and play a weak prisoner's dilemma with
their nearest neighbors. Public reputation influences which neighbor they
consider imitating. The study compares finite-window, exponential, and
present–past memory rules to examine how the weight of the current action,
the length of the historical record, observation errors, and population size
affect cooperation.

## Repository contents

| Location | Contents |
| --- | --- |
| [src/](src/) | Model implementations, simulation drivers, and analysis scripts. |
| [config/](config/) | Parameter grids and simulation settings. |
| [results/](results/) | Individual realizations, summaries, and spatial configurations. |
| [qa/](qa/) | Saved numerical checks and diagnostic records. |
| [residual_drift.csv](residual_drift.csv) | Parameter combinations with residual temporal drift after extending the transient. |
| [DATA_DICTIONARY.md](DATA_DICTIONARY.md) | Definitions of parameters, measurements, summary columns, and drift diagnostics. |
| [REPRODUCIBILITY.md](REPRODUCIBILITY.md) | Dependencies and instructions for reading and reanalyzing the data. |

## Simulation groups

| Directory under results/ | Question |
| --- | --- |
| E1 | Cooperation as a function of temptation payoff and memory length. |
| E2 | Spatial configurations associated with selected parameter combinations. |
| E3 | Separate effects of the weight of the current action and the length of the past window. |
| E4 | Reputation gaps and temporal correlations extracted from the simulated dynamics. |
| E5 | Effects of observation errors. |
| E6 | Effects of lattice size. |
| S1 | Sensitivity to the reputation offset. |
| S2 | Comparison of synchronous and sequential strategy updates. |
| S3, S3_floors, S3_past | Proportional reference selection, reputation offsets, and the use of current versus past actions. |
| transient_refinement | Results from extended transient durations. |

## Reading the data

The lattice side length is `L`, so the population contains `L * L` individuals.
Each realization uses a recorded random seed. The outputs include mean
cooperation, absorption indicators, reputation gaps, and temporal correlations.
Summary tables report averages and standard errors across realizations.

Some parameter points have outputs for more than one transient duration.
The `source` column identifies the selected result file using a path relative
to the repository root. This preserves the connection between summaries and
individual realizations.

The residual drift table contains 14 distinct parameter combinations. `R` is
the number of realizations, and `drifting_realizations` counts those exceeding
the drift criterion. Empty parameter cells indicate values that do not apply
to that memory rule. Definitions are given in the data dictionary.

## Reanalysis

The CSV and JSON files can be inspected without running simulations.
Python scripts regenerate summaries and figures from the saved outputs.
See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the available commands.

`SHA256SUMS.txt` lists the SHA-256 digest of each distributed file.
