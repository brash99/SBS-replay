# CDet Run 5710 Timing Calibration Record

## Summary

This document records the production of the new CDet timing calibration from
cross-target run 5710. The calibration was generated with the two-pass,
in-session cross-target driver and then checked by rerunning the analysis at
calibration stage 7.

The final calibration file is:

```text
CDet_calibration_dt.dat
```

It contains per-pixel timing offsets, an ECal-time correction, and separate
Layer-1 and Layer-2 time-walk corrections. The final analysis processed all
2,951,891 available events from run 5710.

## Files preserved for comparison

Before starting the calibration, the existing calibration was copied to:

```text
CDet_calibration_dt.before_run5710.dat
```

The result of the first successful sequence, which still used an 8 ns minimum
ToT for the Layer-1/Layer-2 pairing step, was preserved as:

```text
CDet_calibration_dt.after_8ns_pairing.dat
CDet_pixel_timing_fit_results.after_8ns_pairing.dat
tdcPlots/run5710_after_8ns_pairing/
```

These files distinguish the effect of the pairing-cut change from the much
larger difference between the new calibration and the pre-run-5710
calibration.

## Problems corrected before the production run

### Driver argument mismatch

The calibration drivers contained calls written for an older
`plotCDetLayersTimeComp` signature. The current signature has `pixelBase`
between `overwrite` and the histogram width. Without that argument, the old
call assigned a negative value to the width and stopped the sequence during
the first ECal fit.

The affected drivers were updated to pass `pixelBase = 416`. Their ECal ADC
time window was also updated to 10--35 ns and their ECal-minus-CDet plotting
range to -60--30 ns.

### ROOT loading order

The driver refers to functions and status globals declared in the master
macro. It therefore cannot be loaded first in a clean ROOT session. The
working order is:

```cpp
.L PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C+
.L Run_CDet_Calibration_TwoPass_InSession_AllCross_Individually.C
```

### Layer-pair ToT selection

The first completed calibration used a minimum ToT of 8 ns in the
Layer-1/Layer-2 pairing fits. The ToT-versus-ToT distribution showed that this
removed a visible population of otherwise useful hits. The driver calls were
changed from 8--40 ns to 4--40 ns for this pairing step.

This is separate from:

- the main analysis ToT selection of 4--50 ns; and
- the time-walk fit range of 5--30 ns stored in the calibration file.

The ToT plot axes were corrected to say ToT rather than LE. Four empty
per-bar diagnostic PDFs were also repaired by filling the `hBarGoodLe`
histograms from the retained per-bar data.

## Production command

From `scripts/cdet`, the successful production run was started in ROOT with:

```cpp
.L PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C+
.L Run_CDet_Calibration_TwoPass_InSession_AllCross_Individually.C

Run_CDet_Calibration_TwoPass_InSession_AllCross_Individually(
    5710,       // run number
    -1,         // all events
    0,          // group index; irrelevant in single-run mode
    -1, -1,     // all available segment files
    0.02, 60,   // CDet leading-edge limits, ns
    4, 50,      // main CDet ToT limits, ns
    1, 100,     // Layer-1 occupancy limits
    1, 100,     // Layer-2 occupancy limits
    0.10,       // ECal-CDet x-difference limit, m
    0.02,       // x offset, m
    0.1,        // y offset, m
    3,          // require both CDet layers
    false,      // do not suppress listed bad channels
    30,         // unused for a single run
    2,          // maximum stream
    1,          // first event; currently unused
    true        // remove the active calibration and start fresh
);
```

The final `true` is destructive with respect to the active
`CDet_calibration_dt.dat`; preserve that file before using this option.

## Calibration sequence

The driver ran the following sequence in one ROOT session:

```text
Stage 0: uncalibrated inspection
Stage 1: first per-pixel offset fit
Stage 2: apply and inspect the first pixel offsets
Stage 3: first ECal-time fit
Stage 6: first time-walk fit
Stage 1: second per-pixel offset fit
Stage 3: second ECal-time fit
Stage 6: second time-walk fit
Stage 7: full-correction pixel-offset closure fit
Stage 7: final calibrated analysis state
```

The calibration is iterative. A valid pixel fit adds a residual correction to
the existing pixel offset rather than replacing it. For the ECal correction,
the residual slope can likewise be added to the loaded `p1`. The fitted
intercept needs different treatment because the correction's mean-restoring
shift makes `p0` cancel from the final corrected time, as discussed below.

## Final calibration constants

After the automated sequence, the stage-7 event-display analysis still found
the following residual dependence:

```text
<t_CDet> = p0 + p1 * t_ECal
p0 = 27.0283 +/- 0.0380626 ns
p1 = 0.136368 +/- 0.00176527
```

Because this was a residual fit made after applying the existing calibration,
its slope was added to, not substituted for, the then-current slope:

```text
p1:  0.824191 +  0.136368 =  0.960559
```

An initial manual update also changed `p0` from 21.717814 to 48.746114. A
second stage-7 analysis then measured:

```text
<t_CDet> = p0 + p1 * t_ECal
p0 = 29.5188 +/- 0.0378898 ns
p1 = 0.0230989 +/- 0.00175628
```

The slope was not consistent with zero, so the second residual slope was
added:

```text
p1: 0.960559 + 0.0230989 = 0.9836579
```

The intercept of approximately 29.5 ns was not added. The applied ECal
correction is

```text
t_corrected = t - (p0 + p1 * t_ECal) + delta
```

where `delta` is recalculated to restore the configured target mean CDet
time. Since that recalculation includes the opposite dependence on `p0`, the
stored intercept cancels from `t_corrected`. The intercept returned by a fit
to already-corrected data is consequently near the restored mean; it is not a
residual intercept that should be accumulated. The driver's automatic
`p0 += fitted_p0` behavior should not be used for such a closure update.

The final global sections of `CDet_calibration_dt.dat` are:

```ini
[ECalTiming]
p0 48.746114
p1 0.9836579

[TimeWalk]
p1_L1 13.812842
p1_L2 16.733828
totref_L1 12.000000
totref_L2 12.000000
totmin 5.000000
totmax 30.000000
```

The final `p1` value should be checked with another stage-7 analysis. The
target result is a residual slope consistent with zero.

### Context for the large ECal intercept change

The pre-run-5710 calibration used `p0 = -43.152900 ns`, whereas the new file
uses `p0 = 48.746114 ns`, a change of approximately +91.9 ns. The ECal group
independently shifted the ECal timing calibration in the database by roughly
100 ns, and this analysis now uses that new ECal calibration. The scale and
direction of the CDet/ECal intercept change are therefore qualitatively
consistent with the known upstream ECal timing change rather than, by
themselves, indicating a comparable shift in the CDet hardware timing.

## Stage-7 validation and event display

The reusable analysis and event-display settings are in:

```text
CDet_run5710_event_display.conf
```

For validation, it must contain:

```text
analysis.calibration_stage: 7
```

Stage 7 applies all three detector corrections: pixel offsets, ECal timing,
and layer-dependent ToT time walk. A clean validation session is:

```cpp
.L PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C
PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget(
    "CDet_run5710_event_display.conf");
plotCDetLayersTimeComp("CDet_run5710_event_display.conf");
```

For the configured Layer-1 pixel 471, the display routine normalizes the
selection to bar 29, covering pixels 464--479. The validation analysis built
11,221 events containing accepted Layer-1/Layer-2 pairs for that bar.

Display-only parameters can be edited in the configuration file and the
display function rerun without repeating the main event analysis. Changes to
the main analysis selection or calibration stage require rerunning the main
analysis.

## Comparison with earlier calibrations

Changing the pairing minimum from 8 ns to 4 ns had only a small effect on the
result. Among the 790 pixels with valid fits in both runs, the median absolute
offset change was 0.043 ns and the RMS change was 0.153 ns.

The comparison to `CDet_calibration_dt.before_run5710.dat` was much larger.
For the 768 pixels with nonzero offsets in both files:

| Quantity | All common pixels | Layer 1 | Layer 2 |
|---|---:|---:|---:|
| Number compared | 768 | 354 | 414 |
| Median signed change | +2.389 ns | +2.340 ns | +2.480 ns |
| RMS change | 3.962 ns | 3.679 ns | 4.189 ns |

The old file contained 857 nonzero pixel offsets and the new file contained
813. The global corrections also changed substantially:

| Constant | Before run 5710 | Automated 4 ns result | Final result |
|---|---:|---:|---:|
| ECal `p0` | -43.152900 | 21.717814 | 48.746114 |
| ECal `p1` | 0.986600 | 0.824191 | 0.9836579 |
| Time-walk Layer 1 | 19.736000 | 13.812842 | 13.812842 |
| Time-walk Layer 2 | 18.037000 | 16.733828 | 16.733828 |

The pixel offsets should not be interpreted independently of these global
changes. The ECal and time-walk models changed at the same time, so this is a
comparison of two complete calibration systems rather than a direct
measurement of paddle drift.

## Related documentation

The detailed derivation of the per-pixel offset algorithm, sign conventions,
fit validation, and complete correction chain is in:

```text
CDet_CROSSTARGET_TIMING_CALIBRATION.md
```
