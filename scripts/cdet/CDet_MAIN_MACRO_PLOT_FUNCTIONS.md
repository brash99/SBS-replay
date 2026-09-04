# CDet main-analysis plotting-function reference

This document describes the plotting and interactive diagnostic functions in
`PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C`. It is
intended as a practical guide to what each function shows, when to use it, and
whether it can change the active timing calibration.

## Before using the plotting functions

First run `PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget(...)`
in the same ROOT session. The plot functions consume event and hit vectors held
in memory by the main macro; most do not reread the replay ROOT file. Loading the
macro without running the analysis does not populate those vectors.

For routine work, prefer the configuration-file overload of the layer-comparison
function:

```cpp
plotCDetLayersTimeComp("CDet_run5710_event_display.conf")
```

Pixel IDs in this document are logical CDet pixel IDs, from 0 through 2687. A
bar contains 16 logical pixels, two of which are normally uninstrumented. Layer
1 is IDs 0--1343 and Layer 2 is IDs 1344--2687.

### Calibration-safety warning

Most functions are read-only diagnostics. The following can write calibration
constants when their write/overwrite argument is enabled:

- `plotGoodLeVsTotByLayer(true, ...)` updates the layer time-walk coefficients.
- `plotCDetLayersTimeComp(true, ...)`, or a configuration with
  `display.overwrite = true`, updates the ECal timing parameters.
- `plotAllTDC(true, ...)` can run the older residual-offset calculation and
  update pixel offsets. It is not the current production offset method.
- `extractHierarchicalCDetPixelTimingOffsets(true, ...)` runs the current
  production per-pixel offset method and updates the active calibration file.

Use `false` for inspection unless an update is intentional.

## Quick reference

| Function | Main purpose | Writes calibration? |
|---|---|---|
| `plotECalHCalTimeComp` | Compare ECal, HCal, and selected CDet-pixel timing | No |
| `plotGoodLeVsTotByLayer` | Measure and inspect layer-dependent time walk | Optional |
| `plotNumAdjacent` | Compare raw and good adjacent-hit multiplicities | No |
| `plotAveTotPerPixel` | Show mean TOT versus pixel ID | No |
| `plotSingleTot` | Inspect all 16 TOT spectra in one bar | No |
| `printOccupancy` | Print one pixel's occupancy | No |
| `plotOccupancyVsID` | Plot occupancy versus ID in four detector sections | No |
| `plotRawSinglesRateVsID` | Plot raw singles rate versus ID | No |
| `plotBarRateHV` | Display the accumulated bar-rate/HV histogram | No |
| `plotPaddleTOT` | Tune per-pixel TOT selections for one bar | No |
| `plotPaddles` | Fit and display LE spectra for one bar | No |
| `plotAllPaddles` | Produce LE-fit pages for every bar | No |
| `plotAllTDC` | Broad TDC/channel/bar diagnostics; legacy offset route | Optional |
| `plot2DrefVsLE` | Compare reference time with good CDet LE | No |
| `plotCDetLayersTimeComp` | Pair layers, compare timing/position, fit ECal dependence, and build event display | Optional |
| `plotECalCDetTimeCutStudy` | Study an ECal-minus-CDet timing cut and selected-pixel fits | No |
| `plotCDetPixelOffsetMethodDifference` | Compare two calibration files' pixel offsets | No |
| `plotCDetPixelLeAndDtSpectra` | Inspect LE and ECal-CDet delta-t for one pixel | No |
| `plotECalCDetTimeComp` | General ECal/CDet timing and position comparison | No |
| `plotRawXCorrelation` | Plot raw CDet/ECal x correlation after a timing cut | No |
| `plotTimeECalVsCDet` | Plot ECal time versus good CDet LE | No |
| `plotXDiffSections` | Plot projected ECal-CDet x residuals by layer | No |
| `plotPMTRates` | Display channel rates for selected PMT/module sections | No |

## Timing-calibration and timing-quality plots

### `plotGoodLeVsTotByLayer(...)`

Builds LE-versus-TOT histograms separately for Layers 1 and 2, optionally draws
profiles, and fits the profiles with `[0] + [1]/sqrt(TOT)`. This is the primary
plot for checking the time-walk correction.

- With `fitTimeWalk = false`, it is only an LE/TOT diagnostic.
- With `fitTimeWalk = true`, it performs the fits.
- With both `fitTimeWalk = true` and `overwrite = true`, the fitted slope
  corrections are added to the active layer time-walk coefficients and written
  to the calibration file.

### `plotCDetLayersTimeComp(...)`

This is the main post-analysis timing and matching diagnostic. It forms
one-to-one Layer-1/Layer-2 hit pairs using timing, x-position, and y-position;
then displays layer timing, layer-to-layer time differences, hit-ID correlation,
position correlations, ECal/CDet timing, and selected-bar pixel spectra.

It fits

```text
<t_CDet> = p0 + p1 * t_ECal
```

and prints the fitted `p0` and `p1`. With `overwrite = true`, those residual
fit parameters are added to the existing `[ECalTiming]` constants and written
to the active calibration file. With `overwrite = false`, the fit is diagnostic.

The integer `pixelBase`/`display.pixel` may be any Layer-1 pixel in the desired
bar; the function normalizes it to the first pixel of that 16-pixel bar. In an
interactive session it also constructs the event browser described below. In
ROOT batch mode the event browser is intentionally omitted.

Two interfaces are available:

```cpp
plotCDetLayersTimeComp(false, 471, /* remaining numeric arguments */)
plotCDetLayersTimeComp("CDet_run5710_event_display.conf")
```

The configuration form is preferred because it avoids a long positional
argument list and rereads the file each time the function is called. You can
therefore change the selected bar, cuts, or display ranges and rerun this plot
without repeating the main analysis.

### `plotECalCDetTimeCutStudy(...)`

Studies the ECal-CDet quantity `t_ECal - t_CDet,LE` before and after a proposed
delta-t cut. It shows the uncut distribution versus pixel ID, the global and
selected-pixel distributions, post-cut LE/TE/TOT spectra, and ECal-versus-CDet
timing. It can optionally draw delta-t versus TOT for the selected pixel.

This function also exercises the local peak/background fit checks used during
development of the per-pixel timing method. It is diagnostic only; it does not
write constants.

### `plotECalHCalTimeComp(...)`

Compares ECal and HCal ADC timing with timing selections on two chosen CDet
pixels. It is useful for establishing detector-to-detector timing relationships,
testing an HCal timing gate, and studying ECal amplitude or position versus time.
The `cdet_shift` parameter provides an overlay shift for visual comparison; it
does not alter calibration constants.

This function automatically saves several PDFs in the current directory,
including `ECalTimeComparison.pdf`, `HCalTimeComparison.pdf`,
`CDetECalTimeOverlay.pdf`, `ECalVsCDetShifted.pdf`, and associated ECal
amplitude/position plots.

### `plotECalCDetTimeComp(...)`

Provides a general ECal/CDet comparison using all good CDet hits. It overlays
cut and uncut ECal-minus-CDet timing, fits the main uncut peak, and plots timing
and projected x-position correlations. This is a useful broad diagnostic, but
the layer-pairing function is preferable when a clean two-layer CDet track is
required.

### `plotTimeECalVsCDet(...)`

Makes the simple two-dimensional plot of ECal ADC time versus good CDet LE after
LE and TOT selections. Use it for a quick check of residual timing dependence;
use `plotCDetLayersTimeComp` for the paired-layer fit used to update ECal timing.

### `plot2DrefVsLE(...)`

Plots the stored raw reference time against good CDet leading-edge time. It is
mainly useful when diagnosing reference-time alignment and requires the
reference and good-hit vectors to have been populated consistently.

## Per-pixel, per-bar, and detector-wide plots

### `plotSingleTot(pixel_base, raw, width, totMin, totMax)`

Draws a 4-by-4 page of TOT spectra for the 16 logical pixels in one bar. Set
`raw = true` for raw-hit spectra and `false` for good-hit spectra. `pixel_base`
should be the first logical pixel in the bar (a multiple of 16). Uninstrumented
pixels are visibly marked.

### `plotPaddleTOT(...)`

Displays TOT and the corresponding LE spectrum after a per-pixel TOT selection
for all 16 pixels in a chosen bar. If external means and sigmas are not supplied,
the function fits sufficiently populated TOT spectra and uses approximately a
two-sigma window. An optional HCal timing window can further restrict the LE
sample. Use this to understand or tune TOT acceptance; it does not save a new
calibration.

### `plotPaddles(bar, ...)`

Displays one bar's 16 good-LE spectra with Gaussian fits and a companion set of
LE-versus-TOT plots. It automatically saves `PaddleOffsetBarN.pdf` and
`PaddleOffset2DBarN.pdf`. This is an older, convenient bar-level inspection
function; its Gaussian peak selection is not the current hierarchical offset
algorithm.

### `plotAllPaddles(...)`

Runs the bar-level LE display over all 168 bars and saves one 16-panel PDF per
bar under `CDetPaddleTimes`, optionally in a `saveTag` subdirectory. It is useful
for a fast detector-wide visual survey of missing channels, secondary peaks, and
poor fits. It does not write constants.

### `plotAllTDC(...)`

Builds a broad suite of raw and good LE, TE, TOT, channel-ID, and bar-level
plots, including per-pixel/per-bar fit summaries. With `savePlots = true`, it
writes the summary PDFs beneath `saveDir` (normally `tdcPlots`) and incorporates
`saveTag` in the names.

This function also contains the older residual pixel-offset calculation. Do not
set `overwrite = true` merely to save plots: `overwrite` authorizes calibration
changes. The full timing workflow now uses
`extractHierarchicalCDetPixelTimingOffsets(true, ...)` instead.

### `plotNumAdjacent(nbins)`

Draws raw and good adjacent-hit multiplicity histograms side by side. It is a
quick diagnostic for clustering, cross talk, and the effect of hit-quality cuts.

### `plotAveTotPerPixel()`

Plots average TOT versus logical pixel ID. It is useful for finding systematic
gain differences, dead/low-response pixels, and detector-section discontinuities.

## Occupancy and rate plots

### `printOccupancy(pixel, applyTotCut)`

Prints the accepted hits per analyzed event for one logical pixel. With
`applyTotCut = false` it uses raw occupancy; with `true` it counts hits passing
the macro's TOT selection.

### `plotOccupancyVsID(applyTotCut, savePdf)`

Plots occupancy versus logical pixel ID in four panels: Layer-1 left/right and
Layer-2 left/right. Error bars are included. With `savePdf = true`, it saves
`RawOccupancyVsID_runN.pdf` or `TotCutOccupancyVsID_runN.pdf`.

### `plotRawSinglesRateVsID(savePdf)`

Plots the raw singles rate in kHz versus pixel ID in the same four detector
sections. With `savePdf = true`, it saves `RawSinglesRateVsID_runN.pdf`.

### `plotBarRateHV()`

Draws the bar-rate-versus-HV histogram accumulated by the analysis, using a
logarithmic x axis. This helper is meaningful only when that histogram was
filled by the selected input data and analysis path.

### `plotPMTRates(mymodule, mylayer, choice)`

Displays channel-rate distributions for left and right PMT groups in a selected
module and layer. It subdivides the PMTs into 16-pixel panels, marks unused
pixels, and overlays expected-rate guide levels. This is primarily a hardware
rate and channel-health diagnostic.

## Position and correlation plots

### `plotRawXCorrelation(tDiffMin, tDiffMax)`

Plots projected ECal x versus raw CDet x for hits passing an ECal-CDet timing
window. Because the default window (`80`--`100` ns) is historical, explicitly
set a window appropriate to the currently loaded timing calibration.

### `plotXDiffSections(le_min, le_max, tDiffMin, tDiffMax)`

Plots the x residual between CDet and the ECal position projected to the CDet
plane, separately for Layers 1 and 2, plus the selected LE distribution. Its
default `80`--`100` ns timing window is also historical and should normally be
overridden.

## Event-display controls

`plotCDetLayersTimeComp` builds an in-memory list of events containing accepted
pairs for the selected Layer-1 bar. The following helpers operate on that list:

- `ShowCDetEvent(index)` displays a particular accepted event by display index.
- `NextCDetEvent()` and `PreviousCDetEvent()` navigate the event list.
- `PrintCDetEvent()` prints the currently displayed event and pair details.
- `SaveCDetEvent()` saves the current event-display canvas.
- `BuildCDetEventDisplay(...)` is the lower-level builder called by
  `plotCDetLayersTimeComp`; it is normally not called directly.

The event list contains only events surviving the pairing and final selection
cuts used in the most recent call to `plotCDetLayersTimeComp`.

## Pixel-offset comparison and hierarchical-fit diagnostics

### `plotCDetPixelOffsetMethodDifference(regularFile, dtFile)`

Reads the `[PixelOffsets]` sections of two calibration files and plots
`dt-method correction - regular correction` versus logical pixel ID. Pixels for
which either correction is exactly zero are skipped. This function reads files
directly and does not modify either one.

### `plotCDetPixelLeAndDtSpectra(pixel, ...)`

For one pixel, draws the good-hit LE spectrum and the ECal-CDet delta-t spectrum
after an ECal-energy selection. It marks the proposed delta-t fit interval and
prints that pixel's corrections from two calibration files. This is useful for
understanding individual outliers in the method-comparison plot.

### `extractHierarchicalCDetPixelTimingOffsetsDiagnostic(...)`

Although this is an extraction function rather than a `plot...` function, it is
the source of the most useful current offset-fit plots. It groups each PMT's
pixels into groups of eight, fits the combined unnormalized-count spectrum, and
then fits each sufficiently populated pixel in a local window around the group
centroid. If a group fit fails, it falls back to individual fits over the full
configured fit range. Failed individual fits fall back to the group result when
one is available.

It writes:

- a hierarchical ROOT file containing directories, histograms, fit functions,
  canvases, and result metadata;
- per-bar PDF and PNG summaries;
- a detailed text result table; and
- an optional candidate calibration file.

The ROOT file should be opened and its stored canvases inspected directly; avoid
drawing a stale canvas pointer retained from a previously closed file.

### `extractHierarchicalCDetPixelTimingOffsets(generateOffsets, outputTag)`

This is the tuned production wrapper. With `generateOffsets = true`, it uses the
current settings (fit range `-30` to `+10` ns, 35 minimum individual counts,
100 minimum group counts, significance threshold 1.5, and chi-square/NDF limit
15), creates tagged diagnostics under `tdcPlots`, and activates the resulting
pixel offsets in the current calibration file. With `false`, it only explains
that generation is disabled; it does not run a diagnostic pass.

The lower-level legacy functions `extractCDetBarPixelTimingOffsets(...)` and
`extractAllCDetPixelTimingOffsets(...)` also generate fit plots, but they use the
older independent-pixel approach. Retain them for comparison and debugging, not
as the preferred production calibration path.

## Practical recommendations

For a normal calibration review:

1. Run the main analysis or the complete calibration driver.
2. Use `plotAllPaddles` or the hierarchical per-bar pages for a detector-wide
   fit-quality survey.
3. Use `plotCDetPixelLeAndDtSpectra` and the hierarchical ROOT file for suspect
   individual pixels or cross-talk cases.
4. Use `plotGoodLeVsTotByLayer(false, ...)` to inspect time walk.
5. Use `plotCDetLayersTimeComp("CDet_run5710_event_display.conf")` to validate
   layer matching, the ECal timing dependence, and representative events.
6. Use the occupancy/rate functions to distinguish a timing-fit problem from a
   dead, low-rate, or unusually noisy channel.

When comparing results from different calls, remember that ROOT object names are
often reused. Close or rename old canvases if ROOT warns about replacing objects,
and treat the latest call's in-memory histograms as authoritative.
