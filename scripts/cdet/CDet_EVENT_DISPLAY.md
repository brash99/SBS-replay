# CDet ROOT Event Display

This guide describes how to run the single-file cross-target analysis and use
the CDet event display. The display operates on the event data retained by
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
    "CDet_run5710_event_display.conf"
)
```

Wait for this call to finish. It reads the replayed ROOT data, applies the
stage-2 analysis selections, and populates the per-event CDet and ECal vectors
used by the event display. The configuration applies the new 10--35 ns ECal
ADC timing cut and the Layer-1/Layer-2 hit-count limits without relying on
positional arguments. Unknown keys, invalid ranges, and unsupported
configuration versions are rejected before the analysis starts.

## 4. Select a CDet bar and build the display

The event browser is built by `plotCDetLayersTimeComp`. Its second argument can
be any pixel ID in the selected Layer-1 bar. Each bar contains 16 pixels; the
macro normalizes the supplied ID to the bar's first pixel using:

```text
pixelBase = 16 * barNumber
```

For example, the configuration selects pixel 471, which is in Layer-1 bar 29
and is normalized to that bar's base pixel, 464. Run the configured display
with:

```cpp
plotCDetLayersTimeComp("CDet_run5710_event_display.conf")
```

Both calls read the same file, but each uses only its own `analysis.*` or
`display.*` settings. The complete file is
[`CDet_run5710_event_display.conf`](CDet_run5710_event_display.conf).

This call creates the aggregate timing and position plots, including the
selected-bar x-versus-z plot. In an interactive Analyzer session it also:

1. Finds events containing an accepted Layer-1/Layer-2 pair for bar 29.
2. Applies the same pair timing, x-difference, and ECal-to-CDet timing cuts as
   the aggregate analysis.
3. Opens the four-panel event-display canvas at the first accepted event.
4. Opens a small `CDet event display` control window.

The terminal reports how many events were placed in the browser:

```text
[CDet event display] Built N events containing accepted pairs for Layer-1 bar 29.
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
PrintCDetEvent()
SaveCDetEvent()
```

Saved images have names of the form:

```text
CDetEventDisplay_run5710_entry123456_bar29.png
```

## Changing the cuts

The simple call above uses the defaults of `plotCDetLayersTimeComp`, including:

```text
Layer timing difference:  -15 to +15 ns
Layer x difference:       -0.01 to +0.01 m
ECal-CDet timing:         -100 to +100 ns
```

If a study uses different cuts, pass them to `plotCDetLayersTimeComp`. The
event browser automatically receives the same values, so its event selection
remains consistent with the aggregate plots. The relevant leading arguments
are:

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

For example, bar 29 with a wider Layer-1/Layer-2 x-difference window:

```cpp
plotCDetLayersTimeComp(
    false, 464, 1,
    -15, 15,
    -0.02, 0.02
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
