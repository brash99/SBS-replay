# Cross-target run-specific p1 qualification results

## Scope

This document records the frozen-sample run-specific ECal timing-slope (`p1`) qualification survey completed on 2026-09-20.

The survey covers all 41 accepted ECal-singles cross-target calibration runs. Runs 5711 and 6077 are excluded because they are LH2 runs. Runs 5976 and 5979 are excluded because their timing behavior is consistent with an HCal-carried trigger condition rather than the ECal-singles timing model.

The detector master value used for every test was:

```text
p1 = 0.810203 ns/ns
```

All run-specific `shift_ns` constants were already fitted and closed before this survey. The qualification procedure is read-only: no `CDet_run<run>.dat` file was changed. Post-run checksum comparison confirmed that all 41 files remained byte-for-byte identical.

## Method

Each run was tested with:

```cpp
Run_CDet_Qualify_RunSpecificP1_FrozenSample(
    run,
    -1,
    "CDet_run<run>_projection.conf",
    5,
    3.0,
    0.5,
    2.0,
    0.5,
    1.05);
```

The accepted timing sample is frozen once while using the detector master `p1`. Candidate slopes are evaluated on the same events, hits, pairs, and half-bars; observations are not reselected after changing the candidate slope. Five contiguous event folds test whether the apparent correction generalizes across the run.

A run qualifies only if all six gates pass:

1. Residual-slope significance is at least 3 sigma.
2. Layer 1 and Layer 2 residuals have the same sign and agree within the configured tolerance.
3. The residual produces at least a 0.5 ns effect across the central 80% ECal-time span.
4. The algebraic full-sample closure residual is negligible.
5. Every contiguous validation fold satisfies its pull, improvement, and candidate-stability requirements.
6. The candidate correction does not degrade the paired-time RMS by more than 5%.

## Survey result

Four of 41 runs passed all six gates:

| Run | Residual p1 (ns/ns) | Significance | Candidate p1 (ns/ns) | Central-80% timing effect (ns) | RMS ratio after/before |
|---:|---:|---:|---:|---:|---:|
| 3575 | +0.158236 +/- 0.002674 | 59.18 | 0.968439 | 2.921 | 0.886 |
| 4722 | -0.096859 +/- 0.014005 | 6.92 | 0.713344 | 0.695 | 0.993 |
| 5066 | -0.131426 +/- 0.012276 | 10.71 | 0.678777 | 1.028 | 0.991 |
| 5068 | -0.114273 +/- 0.012781 | 8.94 | 0.695930 | 0.890 | 0.994 |

The remaining 37 runs did not qualify. Counts of failed gates among those runs were:

| Gate | Non-qualified runs failing gate |
|---|---:|
| Contiguous-fold validation | 36 |
| Practical timing effect | 26 |
| Discovery significance | 13 |
| Layer agreement | 3 |
| Algebraic closure | 0 |
| Resolution preservation | 0 |

These counts overlap because a run may fail more than one gate.

Notable diagnostic cases include:

- Run 3573 has a residual of `+0.157312 +/- 0.002449 ns/ns`, nearly identical to qualified Run 3575, but it fails contiguous-fold validation. This supports the possibility of unstable conditions near the beginning of that running period rather than a clean detector-wide slope change.
- Run 4344 has a large residual (`-0.130010 +/- 0.015108 ns/ns`) and a 1.138 ns central-span effect, but fails fold validation.
- Run 3844 passes every gate except the practical-effect requirement: its effect is 0.391 ns, below the required 0.5 ns.
- Reference Runs 5710 and 5992 show no evidence for a run-specific slope: their residuals are `-0.001347 +/- 0.001693 ns/ns` and `+0.010425 +/- 0.013735 ns/ns`, respectively.

## Calibration decision

There is no strong detector-physics motivation for `p1` to change from run to run. The survey is broadly consistent with the 5710-derived master value, while many nominally significant run residuals do not generalize across contiguous event folds. Furthermore, the four qualifying candidates are not mutually consistent with a single alternative detector response.

Therefore:

1. Retain `p1 = 0.810203 ns/ns` as the detector master value for all cross-target runs.
2. Continue to use `shift_ns` as the routine run-dependent timing correction.
3. Treat a positive qualification result as a diagnostic flag for further run-stability investigation, not as automatic authorization to install a run-specific `p1`.
4. Do not install the four candidate values without independent evidence for a stable detector-response change.

## Authoritative artifacts

- `CDet_cross_target_p1_qualifications.tsv`: consolidated gate values and outcomes for all 41 runs.
- `CDet_run<run>_frozen_p1_qualification.txt`: complete numerical result for each run.
- `CDet_run<run>_frozen_p1_samples.tsv`: frozen observation set used for each run.
- `CDet_cross_target_p1_qualification_logs/`: complete ROOT output for each qualification.
- `CDet_CROSS_TARGET_RUN_INVENTORY.md`: accepted runs, ECal-time windows, and explicit exclusions.
