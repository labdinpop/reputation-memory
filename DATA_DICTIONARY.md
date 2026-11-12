# Supplementary data table: residual drift

File: `residual_drift.csv` (UTF-8, comma-separated, one header row).
The table has 14 rows, one per distinct parameter combination selected after
transient extension. It summarizes saved individual simulation results.

For each realization, the drift diagnostic is the absolute difference between
the mean cooperation fractions in the two halves of the measurement window.
A realization is counted when this difference exceeds 0.02. A parameter
combination appears in this table when the counted fraction exceeds 0.10.
These are diagnostic thresholds, not significance tests.

| Column | Meaning |
| --- | --- |
| `kernel` | Memory rule: `W` for finite window, `E` for exponential, `PP` for present and past weights. |
| `L` | Side length of the square lattice; population size is `L * L`. |
| `b` | Temptation payoff in the weak prisoner's dilemma. |
| `M` | Memory length for `W` and `E`; empty for `PP`. |
| `a` | Current-action weight for `PP`; empty for `W` and `E`. |
| `m` | Past window for `PP`; empty for `W` and `E`. |
| `epsilon` | Probability of an observation error. |
| `w` | Weight of reputation in reference selection. |
| `K` | Fermi sensitivity parameter for normalized payoffs. |
| `eta` | Reputation offset in reference selection. |
| `T0` | Number of discarded transient rounds in the selected output. |
| `T1` | Number of rounds in the measurement window. |
| `R` | Number of independent realizations at this parameter combination. |
| `drifting_realizations` | Number of realizations with absolute `half_difference` greater than 0.02. |
| `drifting_fraction` | `drifting_realizations / R`, a fraction between zero and one. |
| `experiments` | Internal experiment identifiers used to locate the combination. Reused combinations are counted once. |
| `source` | Relative path to the CSV containing the individual realizations used for this row. |

The `source` files include the realization index, random seed, parameters,
cooperation fraction, absorption indicators, reputation gap, and drift
diagnostic. Sampling occurs every 10 rounds for all combinations in this table.
The `source` field documents provenance; it does not cite an external study.
