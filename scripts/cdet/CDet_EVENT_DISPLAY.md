# CDet ROOT Event Display

This guide describes how to run the configuration-driven CDet analysis and use
the unified event display for cross-target or LH2 runs. The display operates on the event data retained by
`PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C`; it is not a
standalone replay program.

## 1. Start in the CDet scripts directory

```bash
cd /Users/brash/CDet_replay/git-repo/sbs_devel/SBS-replay/scripts/cdet
analyzer
```

Starting `analyzer` rather than plain `root` ensures that the Hall A Analyzer
and the installed SBS library are available.

## 2. Load the analysis macro

At the Analyzer prompt, load the macro:

```cpp
.L PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C
```

For a compiled ACLiC load, which is useful while developing or checking the
macro, use:

```cpp
.L PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C+
```

## 3. Run the main analysis

For the run-5710 example:

```cpp
PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget(
    "CDet_run5710_projection.conf"
)
```

Wait for this call to finish. It reads the replayed ROOT data, applies the
stage-7 analysis selections, and populates the per-event CDet and ECal vectors
used by the event display. The configuration applies the Run 5710 10--35 ns ECal
ADC timing cut and the Layer-1/Layer-2 hit-count limits without relying on
positional arguments. Unknown keys, invalid ranges, and unsupported
configuration versions are rejected before the analysis starts.

## 4. Select a CDet bar and build the display

The event browser is built by `plotCDetLayersTimeComp`. With the recommended
configuration-file interface, select the Layer-1 pixel with the `display.pixel`
key. The optional second function argument is `overwriteOverride`; it is not a
pixel ID. Each bar contains 16 pixels, and the macro normalizes the configured
pixel to the bar's first pixel using:

```text
pixelBase = 16 * barNumber
```

For example, the current configuration selects pixel 480, the base pixel of
Layer-1 bar 30. Run the configured display
with:

```cpp
plotCDetLayersTimeComp("CDet_run5710_projection.conf")
```

Only the legacy positional overload takes `pixelBase` as its second argument.
Do not mix that calling convention with the configuration-file overload.

Both calls read the same file, but each uses only its own `analysis.*` or
`display.*` settings. The complete file is
[`CDet_run5710_projection.conf`](CDet_run5710_projection.conf).

The main analysis and plotting calls must use the same configuration. The plot
function consumes the in-memory sample produced by the main call; it cannot
retroactively replace a sample loaded with the legacy positional interface.
Saved per-pixel LE-versus-ToT polygons are not applied by this event-display
path.

This call creates the aggregate timing and position plots, including the
selected-bar x-versus-z plot. In an interactive Analyzer session it also:

1. Finds events containing an accepted Layer-1/Layer-2 pair for bar 30.
2. Applies the same pair timing, x-difference, and ECal-to-CDet timing cuts as
   the aggregate analysis.
3. Opens the four-panel event-display canvas at the first accepted event.
4. Opens a small `CDet event display` control window.
5. Fits the selected bar's calibrated `t_ECal - t_CDet,corr` peak from the
   individual hits in its accepted Layer-1/Layer-2 pairs. This fit defines the
   optional best-hit display window.

The terminal reports how many events were placed in the browser:

```text
[CDet event display] Built N events containing accepted pairs for Layer-1 bar 30.
```

If `N` is zero, no event passed all the active pair-selection cuts. See
"Changing the cuts" below.

## 5. Read the four display panels

The upper-left panel is the principal x-z event view:

- Blue open circles are retained Layer-1 CDet hits.
- Green open triangles are retained Layer-2 CDet hits.
- The black star is the reconstructed ECal position.
- The dashed gray line runs from the nominal target position `(0,0)` to ECal.
- A thick red segment connects each accepted Layer-1/Layer-2 pair.
- Magenta vertical segments show the x residual between each selected CDet hit
  and the straight-line projection from the target to ECal.

All displayed and reported CDet x coordinates use a layer-specific scale and
offset correction:

```text
Layer 1: corrected x = GoodX * XCorr1 - CDetXOffset1
Layer 2: corrected x = GoodX * XCorr2 - CDetXOffset2
```

The macro currently defines `XCorr1 = 1.08`, `XCorr2 = 1.08`, and both offsets
as 0.03 m. The constants remain separate so either layer can be retuned
independently. ECal x is unchanged. Pairing and its x-difference use these
corrected CDet coordinates as well.

The upper-right panel shows the corresponding y-z projection. Because CDet has
almost no hit-position resolution along a paddle, a CDet hit is drawn as a
line spanning the paddle half-length on either side of its nominal y center,
rather than as a point-like y measurement. The macro currently uses a
half-length of 0.30 m. Accepted paddles are overlaid as thicker red lines.

The lower-left panel is an x-y face view. It uses the same paddle-length lines;
the red connector identifies an accepted Layer-1/Layer-2 pair. It also shows
the ECal `(x,y)` position projected from the target back to each CDet layer:
the blue star is the Layer-1 projection and the green star is the Layer-2
projection. The lower-right panel gives the run, tree entry, ECal quantities,
hit counts, pair timing and position differences, and projection residuals.

The display index counts only events accepted into the browser. The `tree
entry` in the title and information panel identifies the original analyzed
tree entry.

## 6. Navigate, print, and save

Use the buttons in the `CDet event display` control window:

- **Previous** displays the previous accepted event.
- **Next** displays the next accepted event.
- **All hits** displays every retained good CDet hit, matching the historical
  event-display behavior.
- **Best CDet hits** displays only hits whose calibrated timing satisfies
  `abs((t_ECal - t_CDet,corr) - peak_mean) <= 3*peak_sigma`. Accepted pair
  overlays are shown only when both hits pass this timing requirement.
- **Print** writes the current event and pair values to the Analyzer terminal.
- **Save PNG** saves the current four-panel canvas in the working directory.

The event list wraps around: advancing past the last event returns to the first
event, and moving backward from the first event goes to the last.

The same operations can be entered directly at the Analyzer prompt:

```cpp
ShowCDetEvent(0)       // First accepted event; indices are zero-based
ShowCDetEvent(25)      // Accepted display event 25
NextCDetEvent()
PreviousCDetEvent()
ShowAllCDetHits()
ShowBestCDetHits()
PrintCDetEvent()
SaveCDetEvent()
```

The canvas title and information panel identify the active hit mode. The
information panel also reports the fitted peak, sigma, and the number of
displayed hits and pairs relative to the retained totals. Switching modes does
not rerun the analysis, change the event list, alter calibration constants, or
modify any physics selection. It changes only which retained hits and accepted
pairs are drawn for the current event.

Saved images have names of the form:

```text
CDetEventDisplay_run5710_entry123456_bar30_all_hits.png
CDetEventDisplay_run5710_entry123456_bar30_best_hits.png
```

## Best-hit timing configuration

The default display mode remains `all`. The best-hit timing peak is fitted
automatically from the individual corrected times of both hits in every
accepted pair whose Layer-1 hit lies in the selected bar. The modal bin seeds a
local Gaussian fit. The default display window is three fitted standard
deviations on either side of the fitted mean.

The following optional `display.*` settings control this behavior:

```text
display.hit_mode: all
display.best_hit_nsigma: 3
```

Set `display.hit_mode: best` to open the browser initially in best-hit mode.
The control-window buttons can still switch modes afterward.

For a reviewed peak whose numerical values should be fixed rather than
refitted, specify both of the following values in nanoseconds:

```text
display.best_hit_peak_mean: -26.2
display.best_hit_peak_sigma: 3.7
```

The mean and sigma overrides must be supplied together, and sigma must be
positive. If they are absent, the automatic selected-bar fit is used. If that
fit is unavailable—for example, because the run has too few accepted hits—the
browser reports a warning, remains in all-hits mode, and disables the best-hit
view for that sample.

## Changing the cuts

The configuration overload uses the `display.*` values in the same file that
identified the analyzed sample. For the current Run 5710 configuration, the
principal pair/display ranges are:

```text
Layer timing difference:  -15 to +15 ns
Layer x difference:       -0.15 to +0.15 m
ECal-CDet timing:         -150 to +150 ns
```

For a display-only study, edit the corresponding `display.*` values in a copy
of the configuration and rerun `plotCDetLayersTimeComp` with that file. The
event browser automatically receives the same values, so its event selection
remains consistent with the aggregate plots. If any `analysis.*` value changes,
rerun the main analysis first with the same configuration.

The legacy positional interface remains available for development work. Its
relevant leading arguments are:

```cpp
plotCDetLayersTimeComp(
    overwrite,
    pixelBase,
    histogramWidth,
    layerDtMin,
    layerDtMax,
    layerDxMin,
    layerDxMax
    // remaining arguments retain their defaults unless explicitly supplied
)
```

For example, bar 30 with a narrower Layer-1/Layer-2 x-difference window:

```cpp
plotCDetLayersTimeComp(
    false, 480, 1,
    -15, 15,
    -0.08, 0.08
)
```

## Rebuilding the browser manually

Normally there is no reason to invoke `BuildCDetEventDisplay` directly because
`plotCDetLayersTimeComp` does it with matching cuts. If selection vectors or
pairing results are changed interactively, rerun `plotCDetLayersTimeComp` to
rebuild both the aggregate plots and the browser from a consistent state.

## Batch-mode behavior

The browser and control window are intentionally not constructed when ROOT is
running in batch mode (`-b`). Batch calibration and plot-production jobs can
therefore continue to call `plotCDetLayersTimeComp` without retaining the
event-display snapshots or opening additional canvases.

To use the browser, start `analyzer` interactively as shown at the beginning of
this guide.

## Starting over or changing runs

To analyze another run in the same session, call the main analysis function
with the new run and then call `plotCDetLayersTimeComp` again. The event-display
event list and its current index are cleared and rebuilt for the newly selected
data.

It is also safe to close the event-display canvas and control window, then call
`plotCDetLayersTimeComp` again for another bar. The canvas and controls will be
recreated for the new event list.

For the cleanest reset while changing analysis configurations, exit and start
a new Analyzer session:

```cpp
.q
```
