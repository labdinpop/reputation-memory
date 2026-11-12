# Data dictionary

## Individual realizations

CSV files under `results/*/points/` contain results from individual realizations.
Each realization has its own random seed. Model parameters are identified below.

| Column | Meaning |
| --- | --- |
| `realization`, `seed` | Realization index and recorded random seed. |
| `T0`, `T1`, `sample_every` | Discarded transient, measurement duration, and sampling interval, in rounds. |
| `f_C` | Mean fraction of cooperators over the measurement window. |
| `all_D`, `all_C` | Indicators of a final population containing only defectors or only cooperators. |
| `G` | Reputation gap between current cooperators and defectors, pooling the sampled observations with coexistence. |
| `G_time` | Mean of the reputation gaps calculated separately at each sampled time; a diagnostic quantity. |
| `G_from_Cy` | Gap reconstructed from the weighted correlations of observed actions. |
| `G_approx` | Approximate gap calculated from true strategy histories and the observation error rate. |
| `relative_gap_discrepancy` | Relative difference `(G - G_approx) / G`, when defined. |
| `n_pairs`, `coexistence_samples` | Numbers of individual–time pairs and sampled times used for the gap. |
| `C_s_0`, `C_s_1`, ... | Normalized correlations between the current strategy and true past strategies at each lag. |
| `C_y_0`, `C_y_1`, ... | Corresponding correlations using observed past actions. |
| `half_difference` | Difference between mean cooperation in the two halves of the measurement window. |
| `elapsed_seconds` | Recorded runtime of the realization. |

A gap can be undefined when both strategies never coexist during measurement.
Such realizations remain part of the cooperation and absorption averages.

## Summary files and other formats

Files named `point_summary.csv` aggregate realizations at each parameter point.
For a measured quantity, a column ending in `_se` gives its standard error
across realizations, while `_n` gives the number of valid realizations used.
The `source` field points to the selected individual-result file.

JSON files record parameter settings, selected outputs, execution metadata,
and diagnostic information. NumPy files (`.npy`) contain saved spatial
configurations; strategies are encoded as zero
for defection and one for cooperation.

## Residual drift table

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
The `source` field connects the summary to the individual results used to calculate it.
