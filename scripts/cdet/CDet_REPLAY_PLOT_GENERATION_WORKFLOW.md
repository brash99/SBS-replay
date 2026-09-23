# CDet replay plot-generation workflow

This document gives the complete procedure for generating the current CDet
replay diagnostics for one run after that run has been replayed with the
database-backed calibrated pulse and Layer-1/Layer-2 pair output enabled. It
also distinguishes the per-run plotting macros from calibration-survey and
simulation utilities that are not normally part of this workflow.

The examples use LH2 Run 6077. Replace `6077`, the input directory, and the
output-directory names for another run.

## Prerequisites

The replay ROOT files must contain the current branches under:

```text
earm.cdet.pulse.*
earm.cdet.pair.*
earm.cdet.timing.*
```

In particular, the plotting workflow requires calibrated pulse LE/TE/ToT,
calibration and quality flags, ECal-projection quantities, and the replay-level
pair source-pulse indices. Older ROOT files made before these branches were
implemented must be replayed again.

The files for one replay job may consist of a base file plus ROOT rollover
parts, for example:

```text
cdet_6077_stream0_2_seg0_5.root
cdet_6077_stream0_2_seg0_5_1.root
...
cdet_6077_stream0_2_seg0_5_5.root
```

The dataset helper selects one coherent set and includes all of its rollover
parts. Do not mix an older full-run replay with a newer limited-event replay in
the same selected dataset.

Run the commands from `SBS-replay/scripts/cdet` after sourcing the usual replay
environment. On the local Mac used for this study, the input directory is:

```text
/Users/brash/CDet_replay/sbs/Rootfiles
```

On the ifarm, pass the directory that directly contains the replay output
files. Do not assume that `OUT_DIR` and `OUT_DIR/Rootfiles` are interchangeable;
inspect the directory first.

## Important ROOT syntax

At an interactive ROOT prompt, compile or load a macro first and then call its
function:

```cpp
.L SomeMacro.C+
SomeMacro(arguments);
```

Do not enter `SomeMacro.C+(arguments)` directly at the ROOT prompt. That syntax
is for ROOT's command-line `-q` argument and is parsed differently when entered
interactively.

## Complete interactive workflow

Start ROOT from `scripts/cdet`:

```bash
root -l
```

Set `OUT_DIR` to the directory containing the replay ROOT files. This variable
is required by the timing-window proposal macro; the two main analysis macros
also accept the directory explicitly.

```cpp
gSystem->Setenv("OUT_DIR", "/Users/brash/CDet_replay/sbs/Rootfiles");
```

### 1. Review the ECal timing window

```cpp
.L Run_CDet_Select_ECalTimingWindow.C+
Run_CDet_Select_ECalTimingWindow(6077);
```

This creates:

```text
CDet_run6077_ECalTimingWindow_PROPOSAL.png
CDet_run6077_ECalTimingWindow_PROPOSAL.pdf
CDet_run6077_ECalTimingWindow_PROPOSAL.root
```

The top row contains the ECal and HCal ADC-time spectra. The lower panel shows
HCal versus ECal ADC time. Red vertical lines mark the proposed ECal bounds;
dashed horizontal lines show the same values on the HCal axis for comparison
and do not impose an HCal cut.

The default call is read-only. To preview explicit bounds without changing the
configuration:

```cpp
Run_CDet_Select_ECalTimingWindow(6077, false, -10.0, 10.0);
```

Only after human review, explicit approval can be written to the existing run
configuration with:

```cpp
Run_CDet_Select_ECalTimingWindow(6077, true, -10.0, 10.0);
```

The approval form changes `CDet_run6077_projection.conf`; it is not required
merely to generate plots.

### 2. Generate the calibrated-pulse and pair plots

```cpp
.L Plot_CDet_GoodPulseCandidates_AllTDC.C+
Plot_CDet_GoodPulseCandidates_AllTDC(
    6077,
    "/Users/brash/CDet_replay/sbs/Rootfiles",
    "CDet_run6077_good_pulse_tdc",
    1.0, 0.0, 60.0,
    0.0, 40.0,
    false,
    0.0, -26.0,
    0.020, 5.0,
    2.0);
```

The arguments are:

| Argument | Meaning | Adopted value |
| --- | --- | ---: |
| `runNumber` | Run number | `6077` |
| `inputDirectory` | Directory containing the coherent replay dataset | local `Rootfiles` path |
| `outputDirectory` | New directory for generated files | `CDet_run6077_good_pulse_tdc` |
| `binWidthNs` | LE/TE/ToT histogram bin width | `1.0 ns` |
| `leMinNs`, `leMaxNs` | Corrected-LE display range | `0`, `60 ns` |
| `totMinNs`, `totMaxNs` | ToT display range | `0`, `40 ns` |
| `recoveredOnly` | Plot only pulses absent from the legacy accepted-hit population | `false` |
| `pairResidualCenterM` | Center of pair trajectory residual | `0 m` |
| `pairTimingCenterNs` | Center of ECal-minus-pair-mean timing | `-26 ns` |
| `pairResidualScaleM` | Spatial normalization scale | `0.020 m` |
| `pairTimingScaleNs` | Timing normalization scale | `5 ns` |
| `pairCutRadius` | Radius in normalized trajectory-time space | `2` |

The normalized pair selection is

```text
R^2 = ((r_x - r_x0) / s_x)^2
    + ((delta_t_pair - delta_t0) / s_t)^2,
```

where

```text
r_x = (x_L2,corr - x_L1,corr)
      - (x_ECal / z_ECal) * (z_L2 - z_L1)

delta_t_pair = t_ECal - (t_L1,corr + t_L2,corr) / 2.
```

The adopted radius 2 corresponds to spatial and timing semiaxes of 4 cm and
10 ns. It is intended to generate candidates for a later global ECal-CDet-GEM-
HCal region-of-interest tracking algorithm, not to make the final track choice.

Principal outputs include:

| File | Content |
| --- | --- |
| `CDetGoodPulse_AllTDC.pdf` | Global calibrated LE, TE, and ToT spectra |
| `CDetGoodPulse_AllChannels.pdf` | Pulse pixel/bar occupancy and multiplicity |
| `CDetGoodPulse_Bars_*.pdf` | Per-bar corrected-LE spectra by layer and side |
| `CDetGoodPulse_Bar030_Amalgamated.{png,pdf}` | Historical Bar 30 ECal-CDet timing sequence |
| `CDetGoodPulse_PairTrajectoryDiagnostics.{png,pdf}` | All-pair trajectory residual, best-pair residual, and timing versus trajectory residual with the ellipse |
| `CDetGoodPulse_Detector_Amalgamated.{png,pdf}` | Detector-wide selection sequence and one-entry-per-pair mean timing |
| `CDetGoodPulse_SelectedPairXCorrelation.{png,pdf}` | Mean corrected CDet-pair x versus ECal x projected to the pair mean z, for final ellipse-selected pairs |
| `CDetGoodPulse_SelectedPairXResidual.{png,pdf}` | One-dimensional residual `<x_CDet,corr>_pair - x_ECal projected`, with one entry per final ellipse-selected pair |
| `CDetGoodPulse_SelectedPairXDiagnostics.{png,pdf}` | Side-by-side view of the selected-pair x correlation and its one-dimensional residual |
| `CDetGoodPulse_AllTDC.root` | Histogram objects for further interactive study; normally kept local |

The selected pair-mean timing histogram contains one entry per pair, not one
entry per member pulse. Individual member timing histograms remain available in
the ROOT output for diagnostics.

Because timing is one coordinate of the ellipse, the selected pair-mean timing
width is conditioned by the selection and is not an independent measurement of
the intrinsic detector timing resolution.

### 3. Scan the ECal timing and energy cuts

```cpp
.L Plot_CDet_PairYieldVsECalTimingWindow.C+
Plot_CDet_PairYieldVsECalTimingWindow(
    6077,
    "/Users/brash/CDet_replay/sbs/Rootfiles",
    "CDet_run6077_pair_timing_scan",
    0.5, 10.0, 0.5,
    3.0, 4.5,
    0.0, -26.0,
    0.020, 5.0,
    2.0);
```

The first three numerical scan arguments specify symmetric ECal ADC-time
half-widths from 0.5 to 10 ns in 0.5 ns steps. The next two specify the full
3.0--4.5 GeV ECal energy reference interval. The remaining five reproduce the
adopted radius-2 trajectory-time ellipse.

Although its name mentions timing, this macro produces both studies:

```text
CDetPairYieldVsECalTimingWindow.{png,pdf,csv,root}
CDetPairYieldVsECalEnergyWindow.{png,pdf,csv,root}
```

The energy scan uses symmetric windows about the midpoint of the supplied full
energy interval, 3.75 GeV for the adopted 3.0--4.5 GeV selection.

Each plot separates three quantities so their denominators remain explicit:

1. the fraction of events admitted by the trial cut that contain at least one
   selected CDet pair;
2. the number of selected pairs divided by the number of events admitted by
   the trial cut;
3. the absolute retained pair-event yield divided by the full reference
   population.

The CSV files preserve the numerator and denominator counts and should be used
for numerical comparisons rather than reading values from the plot.

## Non-interactive commands

The same workflow can be run directly from the shell. Set the input path to the
directory containing the ROOT files.

```bash
export OUT_DIR=/Users/brash/CDet_replay/sbs/Rootfiles

root -l -b -q 'Run_CDet_Select_ECalTimingWindow.C+(6077)'

root -l -b -q 'Plot_CDet_GoodPulseCandidates_AllTDC.C+(6077,"/Users/brash/CDet_replay/sbs/Rootfiles","CDet_run6077_good_pulse_tdc",1.0,0.0,60.0,0.0,40.0,false,0.0,-26.0,0.020,5.0,2.0)'

root -l -b -q 'Plot_CDet_PairYieldVsECalTimingWindow.C+(6077,"/Users/brash/CDet_replay/sbs/Rootfiles","CDet_run6077_pair_timing_scan",0.5,10.0,0.5,3.0,4.5,0.0,-26.0,0.020,5.0,2.0)'
```

On the ifarm, replace the local path with the actual replay-output directory.

## Meaning of the other new macros

The following files are useful, but they are not part of the standard per-run
plot-generation sequence.

### `Run_CDet_Extract_RunTimingMeans_FromList.C`

This is a multi-run cross-target survey driver. It reads run numbers from a
text file, invokes the cross-target calibration macro for every run, and writes
the run-by-run mean CDet LE and ECal ADC times to CSV files.

Its full signature is:

```cpp
Run_CDet_Extract_RunTimingMeans_FromList(
    runList, events, calibrationStage,
    minSegment, maxSegment,
    leMin, leMax, totMin, totMax,
    ecalTimeMin, ecalTimeMax,
    leCsv, ecalCsv);
```

It is not the standard way to analyze an individual LH2 run such as 5711 or
6077.

### `run_cdet_timing_means.C`

This is a convenience wrapper for the preceding multi-run survey. It loads the
cross-target calibration macro and processes `runs.txt` with default settings:

```cpp
.x run_cdet_timing_means.C
```

It should not be run merely to generate the calibrated-pulse and pair plots for
one run.

### `SimulateECalCDetTimes.C`

This is a toy Monte Carlo that independently samples ECal and CDet time
distributions from `ECalTime.csv` and `CDetTime.csv`. It is useful for studying
whether apparent two-dimensional timing structures can arise from independent
marginal distributions:

```cpp
.L SimulateECalCDetTimes.C+
SimulateECalCDetTimes("ECalTime.csv", "CDetTime.csv");
```

It does not read a replay run and is not part of routine run processing.

## Recommended production values for LH2 runs

The independent 100,000-event Run 5711 and Run 6077 studies support the same
starting policy:

```text
3.0 < E_ECal < 4.5 GeV
-10 < t_ECal < 10 ns
pair trajectory-time center = (0 m, -26 ns)
pair trajectory-time scales = (0.020 m, 5 ns)
normalized ellipse radius = 2
```

Narrowing the ECal timing window modestly increases the conditional pair
fraction but removes a large fraction of usable candidate events. Narrowing the
energy window provides essentially no pair enrichment. The current broad cuts
are therefore appropriate for candidate generation before global ROI tracking.

For the detailed Run 5711 derivation, numerical tables, and physics reasoning,
see [`CDet_RUN5711_PAIR_SELECTION_STUDY.md`](CDet_RUN5711_PAIR_SELECTION_STUDY.md).

## Troubleshooting

### No coherent dataset found

Confirm that the input argument names the directory containing the ROOT files,
that the filenames contain the requested run number, and that all rollover
parts belong to the same replay job. Do not point to the parent directory unless
the files are actually stored there.

### Missing `earm.cdet.pulse.*` or `earm.cdet.pair.*` branch

The ROOT files predate the new replay implementation or were produced without
the current `replay_CDet.C`. Rebuild against the current SBS-offline library and
replay the run again.

### ROOT reports that `Macro.C+` is not a function

At an interactive prompt, use `.L Macro.C+` followed by the C++ function call.
Reserve `Macro.C+(arguments)` for shell invocations using `root -q`.

### Unexpectedly doubled event count

The input directory likely contains multiple complete replay datasets for the
same run. Move obsolete copies elsewhere or pass a directory containing only
the intended coherent dataset. The dataset helper reports every selected file
at startup; review that list before trusting the plots.
