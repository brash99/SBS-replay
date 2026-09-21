# Run 5710 calibration-parameter evolution

## Purpose

This document records how the commissioned Run 5710 CDet timing constants were
created from a clean start.  It follows the actual executable sequence in
`Run_CDet_Calibrate_Run5710_FromScratch.C`,
`Run_CDet_Calibration_Run5710_Accepted.C`, and
`Run_CDet_Calibration_TwoPass_InSession_AllCross_Individually.C`.

The clean-start driver requires both `CDet_calibration_dt.dat` and
`CDet_run5710.dat` to be absent.  It also explicitly sets every in-memory pixel
offset, ECal timing constant, and time-walk coefficient to zero.  Consequently,
the first row below represents genuinely nonexistent calibration constants,
not legacy defaults left in memory.

## Evolution table

`written` means that a non-zero parameter family was fitted and saved at that
step.  The production workflow updates the same master file in place, so exact
numeric snapshots from intermediate fits were not retained.  Those cells are
therefore identified honestly as overwritten rather than reconstructed from
the final file.

| Step in the commissioned workflow | Analysis stage or operation | Pixel offsets | ECal `p0`, `p1`, `delta` | Time walk `p1_L1`, `p1_L2` | Run `shift_ns` | What changed |
|---:|---|---|---|---|---|---|
| 0 | Clean-start guard and global reset | **Absent; all 2,688 values forced to 0** | **Absent; 0, 0, 0** | **Absent; 0, 0** | **Absent** | Both output files must not exist before execution. |
| 1 | Stage 0: uncalibrated inspection | 0, not applied | 0, not applied | 0, not applied | absent | Establishes the raw timing state; writes no calibration. |
| 2 | Stage 1: first hierarchical pixel-offset fit | **Written, non-zero** | 0 | 0 | absent | Fits ECal−CDet peak differences pixel by pixel and writes the first master file. |
| 3 | Stage 2: apply pixel offsets | loaded and applied | 0 | 0 | absent | Checks the first pixel-alignment result; no new parameter family is introduced. |
| 4 | Stage 8: first absolute ECal fit | loaded and applied | **Written, non-zero** | 0 | absent | Fits the detector-wide ECal slope and intercept and establishes the fixed timing-origin `delta`. |
| 5 | Stage 6: first time-walk fit | loaded and applied | loaded and applied | **Written, non-zero** | absent | Fits separate Layer-1 and Layer-2 ToT coefficients over 5–25 ns. |
| 6 | Stage 1: second pixel-offset fit | **Updated** | retained | retained | absent | Refines residual pixel offsets with the existing calibration loaded. |
| 7 | Stage 3: second ECal fit | retained | **Updated** | retained | absent | Refines the ECal timing constants after the first offset/time-walk cycle. |
| 8 | Stage 6: second time-walk fit | retained | retained | **Updated** | absent | Refines the two layer-dependent time-walk coefficients. |
| 9 | Stage 7: full-closure pixel-offset pass | **Updated** | applied | applied | absent | Applies the complete correction chain and fits remaining hierarchical pixel residuals. |
| 10 | Reviewed 102-cut residual-offset pass | **Updated** | applied and closure-refitted | applied | absent | Reuses the preserved human-reviewed pixel-quality decisions. |
| 11 | Half-bar intercept alignment | **Updated in 16-pixel half-bar blocks** | applied and closure-refitted | applied | absent | Aligns accepted half-bar intercepts at a common ECal reference time. |
| 12 | Final Stage-7 master closure | **Final values written** | **Final values written** | **Final values written** | absent | Produces the accepted detector-wide master calibration. |
| 13 | Run-specific timing-origin fit | final master applied | final master applied | final master applied | **Written, non-zero** | Centers the accepted projected-half-bar pair time at 30 ns and writes `CDet_run5710.dat`. |

## Final commissioned values

The final master file is `CDet_calibration_dt.dat`:

| Parameter | Clean-start value | Final value | Units or interpretation |
|---|---:|---:|---|
| Pixel-offset table | nonexistent; 2,688 zeros in memory | 2,688 rows; 1,964 non-zero offsets | ns added to CDet LE/TE |
| Non-zero pixel-offset range | — | −17.990589 to +17.952877 | ns |
| ECal `p0` | 0 | 14.807420 | ns |
| ECal `p1` | 0 | 0.810203 | ns/ns |
| ECal `delta` | 0 | 31.680784 | ns |
| Time-walk `p1_L1` | 0 | 13.987991 | ns·√ns |
| Time-walk `p1_L2` | 0 | 15.764852 | ns·√ns |
| Time-walk `ToTref_L1` | nonexistent | 12.000000 | ns |
| Time-walk `ToTref_L2` | nonexistent | 12.000000 | ns |
| Time-walk fit interval | nonexistent | 5–25 | ns |

The final run-specific file is `CDet_run5710.dat`:

| Parameter | Clean-start value | Final value | Units |
|---|---:|---:|---|
| Run 5710 `shift_ns` | nonexistent | 1.195534 | ns |

The fact that only 1,681 pixel rows carry direct fit-hit counts while 1,964
offsets are non-zero is expected.  The later half-bar alignment adds a common
increment to all 16 pixels in an accepted half-bar, including pixels without a
direct individual fit count.

## Final correction chain

For an accepted hit, the commissioned timing calculation is summarized by

```text
t_pixel = t_raw + pixel_offset

t_ecal = t_pixel - (p0 + p1*t_ECal) + delta

t_walk = t_ecal
         - p1_layer*(1/sqrt(ToT) - 1/sqrt(ToTref_layer))

t_final = t_walk + shift_ns
```

Thus the final calibrated time evolves from having no correction parameters at
all to using a per-pixel offset, three detector-wide ECal timing constants, one
time-walk coefficient per layer, and one Run 5710 timing-origin shift.

## Closure evidence

The clean, configuration-authoritative reproduction processed all 2,951,891
Run 5710 events.  Its archived acceptance record reports:

| Closure observable | Final result |
|---|---:|
| Accepted pairs before/after the run-shift operation | 217,367 / 217,367 |
| Pair-time core centroid | 30.0000 ± 0.0042 ns |
| Residual fixed-effects ECal slope | 0.000005 ± 0.001681 ns/ns |

These closure quantities demonstrate that the non-zero constants do not merely
exist in the output files: together they reproduce the intended absolute time
origin and remove the measurable within-half-bar ECal-time dependence.

## Scope note

The later physics-motivated ECal-y propagation correction using `n = 1.59` is
an optional analysis correction configured in `CDet_run5710_projection.conf`.
It was developed after the commissioned timing calibration and is not one of
the constants created by the clean master-calibration sequence documented
above.
