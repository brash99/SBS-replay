# Hydrogen CDet Event Display

The hydrogen analysis now uses the same corrected CDet coordinate and event
display conventions as the cross-target analysis while retaining its own run
selection and calibration stages.

## Interactive workflow

From the CDet scripts directory:

```bash
cd /Users/brash/CDet_replay/git-repo/sbs_devel/SBS-replay/scripts/cdet
analyzer
```

At the Analyzer prompt:

```cpp
.L PlotElastic_Calibration_Master_stageflag_singlefile_hydrogen.C+
```

Run the hydrogen master function with the desired run, calibration stage, file
segments, and cuts. For example, the default argument set can be exercised
with:

```cpp
PlotElastic_Calibration_Master_stageflag_singlefile_hydrogen(6077,-1,3,0,0,5)
```

Then select a Layer-1 CDet bar through any of its pixel IDs. Each bar contains
16 pixels, so bar 29 begins at pixel 464:

```cpp
plotCDetLayersTimeComp(false,464,1,-15,15,-0.15,0.15)
```

The complete trailing plotting arguments match the cross-target macro:
`XBinWidth`, `XMin`, `XMax`, `ZBinWidth`, `ZMin`, and `ZMax`.

The second argument is the Layer-1 pixel ID; it is normalized to the beginning
of its bar. Thus pixels 464 through 479 all select bar 29, while pixel 480
selects bar 30.

The call creates the aggregate timing and selected-bar x-z plots and, in an
interactive session, opens the four-panel event display and control bar.

## Display controls

Use the control-window buttons or call the functions directly:

```cpp
ShowCDetEvent(0)
NextCDetEvent()
PreviousCDetEvent()
PrintCDetEvent()
SaveCDetEvent()
```

Closing the display and control windows is safe. A subsequent call to
`plotCDetLayersTimeComp` recreates them for the newly selected bar.

## Coordinate and timing conventions

The hydrogen macro uses:

```text
Layer 1 corrected x = GoodX * XCorr1 - CDetXOffset1
Layer 2 corrected x = GoodX * XCorr2 - CDetXOffset2
```

The current values are `XCorr1 = XCorr2 = 1.07` and
`CDetXOffset1 = CDetXOffset2 = 0.03 m`. The four constants remain independently
adjustable.

ECal is located at `z = 6.144 m`. ECal-CDet timing is consistently defined as:

```text
t_ECal - (t_CDet_L1 + t_CDet_L2)/2
```

CDet y hits are displayed as full paddle extents centered on the nominal
`GoodY` value because CDet has little position resolution along a paddle.
