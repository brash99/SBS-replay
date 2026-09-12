# CDet LH2 Run Calibration from Scratch

This is the student-facing procedure for calibrating a new LH2 coincidence run
when no run-specific calibration file exists. It reproduces the sequence used
for Run 5711: determine the coincidence-time selections, establish the absolute
CDet timing origin, determine the run-specific ECal dependence, restore the
timing origin after the slope converges, and perform a final visual review.

## Meaning of “from scratch”

An LH2 run calibration does **not** redetermine the detector-wide pixel offsets,
time-walk parameters, half-bar alignment, or propagation corrections. Those
must already exist in the commissioned Run 5710 master file:

```text
CDet_calibration_dt.dat
```

“From scratch” here means that `CDet_runRUN.dat` is absent. If the master file
is also absent, first reproduce and approve the Run 5710 detector calibration
using `Run_CDet_Reproduce_5710_5992_FromScratch.C`; do not attempt to derive a
detector calibration from an LH2 coincidence run.

All commands below are entered from `scripts/cdet` in a fresh Analyzer/ROOT
session. Replace `RUN` with the run number and use the matching
`CDet_runRUN_projection.conf` throughout.

## 1. Prepare the authoritative run configuration

Create `CDet_runRUN_projection.conf` from the Run 5711 configuration and change
only values justified for the new run. Before starting any fit, verify at least:

- `analysis.run_number` is the intended run;
- `analysis.events: -1` selects the complete run;
- the ECal energy interval selects the intended electron sample;
- `analysis.ecal_time_min/max` select the physical ECal coincidence peak;
- `display.ecal_time_min/max` and `display.hcal_time_min/max` contain the
  intended coincidence population;
- the LH2 good-hit ToT selection is consistently `8 < ToT < 35 ns` in
  `analysis`, `diagnostics`, and `display`;
- the projected-x match remains `analysis.x_difference_max: 0.08` m unless a
  separate geometry study justifies changing it.

For Run 5711 the accepted working selections are:

```text
3.0 < E_ECal < 4.5 GeV
-10 < t_ECal < 10 ns
-10 < t_HCal < 10 ns
8 < ToT < 35 ns
|x_CDet - x_ECal projected| < 0.08 m
```

The ECal event cut used by the main analysis is
`earm.ecal.adctime`. The ECal timing-window proposal plots
`earm.ecal.a_time`; confirm visually that both identify the same physical
population before accepting a window.

To generate an ECal timing-window proposal without changing the configuration:

```cpp
.L Run_CDet_Select_ECalTimingWindow.C+
Run_CDet_Select_ECalTimingWindow(RUN)
```

Inspect `CDet_runRUN_ECalTimingWindow_PROPOSAL.png`. After a human chooses the
bounds, store them explicitly:

```cpp
Run_CDet_Select_ECalTimingWindow(RUN, true, ECAL_MIN, ECAL_MAX)
```

For Run 5711 the initially approved physics peak was around -1 ns. The final
working ECal interval was widened to -10 to 10 ns so the ECal-dependence fit had
enough lever arm.

## 2. Start with no run-specific file

Preserve any pre-existing file outside `scripts/cdet`, then ensure
`CDet_runRUN.dat` is absent before the first analysis. Do not delete or replace
`CDet_calibration_dt.dat`.

With no run file, Stage 7 automatically uses the accepted detector-master `p1`
and zero run shift. Generate the seed diagnostics:

```cpp
.L Run_CDet_Inspect_RunTimingShift.C+
Run_CDet_Inspect_RunTimingShift(
    RUN, -1, "CDet_runRUN_projection.conf", "seed")
```

This command does not write calibration constants. Inspect, in order:

1. `CDet_runRUN_shift_diagnostics_seed/CDet_projected_halfbar_timing.png`
2. `CDet_runRUN_shift_diagnostics_seed/CDet_layer_LE.png`
3. `CDet_runRUN_shift_diagnostics_seed/cCDetTvsECalTDiagnostics.png`
4. `CDet_runRUN_shift_diagnostics_seed/cCDetTvsHCalT.png`
5. `CDet_runRUN_shift_diagnostics_seed/CDet_bar30_layer1_LE_vs_ToT.png`
6. the two Bar 30 plots under
   `CDet_runRUN_shift_diagnostics_seed/bar_ecal_cdet_timing/run_RUN_stage_7/`

The lower-left projected-half-bar pair-mean histogram is the authoritative
visual estimator for the narrow LH2 CDet peak. The broad arithmetic mean is not
appropriate because the LH2 spectrum has a long upper-time tail.

## 3. Fit the initial run timing shift

Fit the modal projected-half-bar peak locally and move it to 30 ns:

```cpp
.L Run_CDet_Calibrate_RunTimingShift.C+
Run_CDet_Calibrate_RunTimingShift(
    RUN, -1, 30.0, 5.0, "CDet_runRUN_projection.conf", true)
```

The final `true` selects the projected-half-bar estimator. The driver writes
only `[GlobalTiming] shift_ns`, reruns the full analysis, requires the accepted
pair population to remain unchanged, and checks the Gaussian-core closure. Do
not accept the value merely because ROOT returns to the prompt: require the
printed `physical-fit gate` and closure gate to pass.

Regenerate inspection plots with a distinct tag:

```cpp
Run_CDet_Inspect_RunTimingShift(
    RUN, -1, "CDet_runRUN_projection.conf", "after_initial_shift")
```

Confirm visually that the narrow projected-half-bar peak is near 30 ns and that
the selected Bar 30 peaks remain physical.

## 4. Fit the run-specific ECal dependence

Once the timing peak is inside the useful analysis range, fit `p1`:

```cpp
.L Run_CDet_Calibrate_RunSpecificP1.C+
Run_CDet_Calibrate_RunSpecificP1(
    RUN, -1, false, "CDet_runRUN_projection.conf")
```

The first call adds `[ECalTiming] p1` while preserving `shift_ns`. The fitted
quantity is the combined half-bar fixed-effects slope, not the pooled profile
slope. Require the printed closure residual to be consistent with zero under
the driver's closure gate.

If the first update does not close, inspect the candidate before continuing:

```cpp
Run_CDet_Inspect_RunTimingShift(
    RUN, -1, "CDet_runRUN_projection.conf", "p1_iteration_1")
```

If the peak and selected populations remain physical, apply another additive
fixed-effects correction by rerunning with permission to replace the existing
run-specific key:

```cpp
Run_CDet_Calibrate_RunSpecificP1(
    RUN, -1, true, "CDet_runRUN_projection.conf")
```

Repeat the inspect/update pair only as needed. Each iteration must use the same
configuration and full event sample. Stop if the physical population changes,
the fit ceases to be stable, or the residual does not decrease coherently.

For an LH2 coincidence run, a final `p1` near zero can be physically sensible:
after tight ECal, HCal, energy, geometry, and layer-coincidence selections, the
accepted events can share nearly the same event start time. Treat this as an
interpretation supported by control plots, not as an assumption imposed on the
fit.

## 5. Refit `shift_ns` after `p1` closes

Changing `p1` can move the timing origin. Repeat the projected-half-bar shift
fit after the final slope iteration:

```cpp
Run_CDet_Calibrate_RunTimingShift(
    RUN, -1, 30.0, 5.0, "CDet_runRUN_projection.conf", true)
```

This preserves the fitted `p1` and updates only `shift_ns`. Require the
projected Gaussian-core centroid, population-invariance check, and final
closure gate to pass.

## 6. Produce the final review package

Start a fresh Analyzer/ROOT session and run:

```cpp
.L Run_CDet_Inspect_RunTimingShift.C+
Run_CDet_Inspect_RunTimingShift(
    RUN, -1, "CDet_runRUN_projection.conf", "final")
```

Review all saved plots, with particular attention to:

- the projected-half-bar pair-mean peak and its long-time tail;
- Layer-1 versus Layer-2 timing correlation;
- combined and per-layer CDet time versus ECal time;
- the CDet-versus-HCal control plot;
- the fully corrected Bar 30 LE-versus-ToT distribution;
- the 4-by-4 Bar 30 pixel timing spectra;
- the Bar 30 amalgamated and per-pixel ECal-CDet timing fits.

The final `CDet_runRUN.dat` must contain only run-dependent constants:

```ini
[ECalTiming]
p1 VALUE

[GlobalTiming]
shift_ns VALUE
```

Saved per-pixel LE-versus-ToT polygons are not applied in this run-specific
analysis. The configured global ToT interval is an actual good-hit filter, and
the accepted detector-master time-walk correction is applied to every hit that
passes that interval, even outside the narrower ToT range originally used to
fit the time-walk parameters.

## 7. Run 5711 reproduction target

Using `CDet_run5711_projection.conf` and the current commissioned master, the
procedure converged to:

```ini
[ECalTiming]
p1 -0.012423

[GlobalTiming]
shift_ns 3.422000
```

The final nominal-window inspection uses all 719,218 events, retains 80,645
good CDet events, and produces 47,235 projected-half-bar pairs. The combined
fixed-effects residual slopes are consistent with zero for both the ECal and
HCal controls. These numerical targets are a reproduction check for Run 5711,
not defaults to copy into another LH2 run.

## 8. Optional selection studies

Selection studies must use separate configuration files and output tags. Never
overwrite the authoritative run configuration merely to make a comparison.
For example, the Run 5711 `12 < ToT < 35 ns` study used:

```cpp
Run_CDet_Inspect_RunTimingShift(
    5711, -1, "CDet_run5711_projection_totmin12_test.conf", "totmin12")
```

That test retained about two-thirds of the nominal accepted pairs while only
modestly narrowing the Bar 30 timing peak. It remains diagnostic evidence, not
the accepted LH2 selection.
