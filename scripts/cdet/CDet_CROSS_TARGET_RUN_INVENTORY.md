# CDet cross-target calibration run inventory

This inventory records the runs admitted to the cross-target timing-calibration workflow and their approved ECal ADC-time windows. The windows are analysis cuts used to define a stable calibration sample; they are not fitted run-by-run timing offsets.

## Included runs

| ECal ADC-time window (ns) | Runs |
|---|---|
| `-15 < t_ECal < 30` | 3573, 3575 |
| `0 < t_ECal < 35` | 3602, 3603, 3604, 3648, 3682, 3685, 3686, 3694, 3757, 3758, 3786, 3788, 3844, 3845, 3846, 4003, 4006 |
| `10 < t_ECal < 35` | 4344, 4345, 4400, 4722, 4735, 4926, 4929, 4985, 5066, 5068, 5290, 5292, 5294, 5423 |
| `10 < t_ECal < 35` | 5710, 5722, 5723, 5724, 5726, 5794, 5795 |
| `-15 < t_ECal < 15` | 5992 |

Runs 4344, 5710, and 5992 already have established projection configurations and run-timing files. They are reference cases and are not silently recalibrated by the bulk `shift_ns` script.

## Explicit exclusions

- 5711 and 6077 are LH2 runs, not cross-target runs.
- 5976 and 5979 are excluded because their timing behavior shows an HCal-carried trigger condition inconsistent with the ECal-singles cross-target calibration model. The run-sheet classification is therefore not trusted for this workflow.

The operational batch list is also hard-coded in `run_cross_target_shift_calibrations.sh`; changing this inventory alone does not add a run to a calibration batch.
