# CDet Cross-Target Timing Calibration

## Purpose

This note documents how the current cross-target analysis calculates and
applies CDet timing offsets. The principal implementation is in
`PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C`, with the
automated sequence in `Run_CDet_Calibration_TwoPass_InSession_AllCross.C`.

The current workflow calculates one offset for each of the 2,688 logical CDet
pixel/channel IDs. It does not calculate only one offset for each 16-pixel
bar. Some variable and function names use "paddle" for a logical pixel, so the
logical ID should be checked when interpreting output.

## Timing value retained by the main analysis

For an accepted hit in logical pixel `i`, the stored CDet leading-edge time is

```text
t_CDet,i = 0.01 * t_TDC - t_reference + O_i(existing)
```

where:

- `0.01` converts the stored TDC value to ns;
- `t_reference` is the event reference-channel time when reference timing is
  enabled, and zero otherwise; and
- `O_i(existing)` is the pixel offset loaded from the active calibration file.

The offset is therefore an additive correction to both the CDet leading and
trailing edge. ToT is not shifted.

Before entering the calibration vectors, a hit must pass the main analysis
selection. This includes:

- CDet leading-edge and ToT limits;
- the ECal ADC timing window, currently 10--35 ns in the cross-target drivers;
- valid ECal reconstruction;
- ECal-to-CDet x and y matching;
- the configured Layer-1 and Layer-2 occupancy limits; and
- the requested `layer_choice`.

With the usual `layer_choice = 3`, the event must contain at least one accepted
hit in each CDet layer. The later per-pixel extraction does not itself require
that an individual hit be paired with an opposite-layer hit.

## Per-pixel ECal-CDet spectra

The all-cross drivers now use the production hierarchical extractor:

```cpp
extractHierarchicalCDetPixelTimingOffsets(true, "hierarchical_pass1_offsets")
```

For every retained hit in pixel `i`, this function fills

```text
delta_t_i = t_ECal - t_CDet,i
```

into a separate histogram for that logical pixel. It also applies an ECal
energy selection, which defaults to:

```text
1 GeV <= E_ECal <= 12 GeV
```

The current production settings are:

```text
Histogram range:  -60 to +30 ns
Group/broad range: -30 to +10 ns
Minimum group entries: 100
Minimum individual entries: 35
Allowed sigma:    0.5 to 8 ns
Maximum individual chi2/NDF: 15
Minimum signal significance: 1.5
```

## Peak fitting

Each 16-pixel bar is divided into contiguous groups 0--7 and 8--15. Actual
counts from the instrumented pixels in each group are summed without
normalization and fitted with a Gaussian signal plus linear background:

```text
f_i(delta_t) = Gaussian(A_i, mu_i, sigma_i) + a_i + b_i * delta_t
```

The ROOT expression is:

```cpp
gaus(0) + pol1(3)
```

The group centroid seeds an adaptive individual fit window with a half-width
between 4 and 8 ns. A valid constrained individual centroid is used directly.
If that fit fails, the group centroid is used. Group chi-square is diagnostic
only because the sum need not itself be Gaussian.

If the group fit is invalid, sufficiently populated pixels are instead fitted
independently across the full -30 to +10 ns interval. A valid result is labeled
`individual_broad_fallback`; otherwise the existing pixel offset is retained.

A fit is rejected when, among other checks, it has insufficient entries,
invalid parameters or uncertainties, a centroid too close to the fit boundary,
a sigma outside the configured limits, a nonpositive amplitude, an invalid
number of degrees of freedom, or excessive chi2/NDF. Known unused pixels are
also excluded.

## Detector-wide reference

The detector reference `mu_0` is the median of all valid constrained and broad
individual-fit centroids. If none are available, the median of valid group
centroids is used. This robust reference prevents a small number of unusually
precise or pathological fits from dominating the common timing origin.

The residual correction for pixel `i` is

```text
c_i = mu_i - mu_0
```

and the stored total offset becomes

```text
O_i(new) = O_i(existing) + c_i
```

This sign follows from the fact that the offset is added to the CDet time:

```text
delta_t_i(new)
  = t_ECal - (t_CDet,i + c_i)
  = mu_i - (mu_i - mu_0)
  = mu_0
```

Thus, after applying the residual correction, every valid pixel peak should be
aligned with the detector-wide reference. A rejected or low-statistics pixel
receives no new residual; any previously loaded offset is retained.

## Calibration stages and iteration

The intended combined cross-target sequence is:

1. Stage 0: inspect the uncalibrated timing.
2. Stage 1: fit ECal-CDet pixel offsets.
3. Stage 3: fit the dependence of the paired mean CDet time on ECal ADC time.
4. Stage 6: fit the layer-dependent ToT time walk.
5. Repeat stages 1, 3, and 6 with the existing corrections loaded.
6. Stage 7: apply all corrections and perform a final residual pixel-offset
   closure pass.
7. Run stage 7 again to inspect the final calibrated state.

The stages control which corrections are applied while constructing the
in-memory hit times:

```text
Stage 0: no calibration
Stage 1: fit pixel offsets; existing pixel offsets are applied when available
Stage 2: apply pixel offsets
Stage 3: fit ECal timing using pixel offsets
Stage 4: apply pixel offsets and ECal timing correction
Stage 5: inspect time walk after pixel and ECal corrections
Stage 6: fit time walk after pixel and ECal corrections
Stage 7: apply pixel, ECal, and time-walk corrections
Stage 8: refit ECal timing with pixel and time-walk corrections applied
```

Because each offset fit adds a residual to the previously loaded constant, the
sequence is iterative rather than a replacement calculation.

## Running the full workflow

For a single cross-target run such as 5710, start a fresh ROOT session and
compile the self-contained driver with ACLiC:

```cpp
.L Run_CDet_Calibration_TwoPass_InSession_AllCross_Individually.C+
Run_CDet_Calibration_TwoPass_InSession_AllCross_Individually(5710)
```

The default final argument is `removeExistingCalibrationFile=true`, which
starts from scratch by removing `CDet_calibration_dt.dat`. Set that final
argument to `false` only when deliberately continuing from the currently
active calibration. Each of the three pixel-offset passes now calls the
hierarchical production extractor and stops the sequence if extraction fails.

For a combined ECal cross-target group, load and call
`Run_CDet_Calibration_TwoPass_InSession_AllCross.C`; its calibration and
diagnostic output names include the group index.

## Other timing corrections

### ECal timing correlation

Accepted Layer-1 and Layer-2 hits are paired. For each pair,

```text
t_pair = (t_L1 + t_L2) / 2
```

A profile is fitted with

```text
<t_pair> = p0 + p1 * t_ECal
```

When the ECal correction is enabled, the macro applies

```text
t' = t - (p0 + p1 * t_ECal) + delta
```

where `delta` restores the configured target mean CDet leading-edge time.

### ToT time walk

Time walk is fitted separately for Layer 1 and Layer 2 using a dependence on
`1/sqrt(ToT)`. The applied correction is

```text
t' = t - p1_layer * (1/sqrt(ToT) - 1/sqrt(ToT_reference,layer))
```

The correction is zero at the layer's reference ToT.

### Run-dependent global shift

After the detector corrections, the macro can add a final run-dependent shift
from `CDet_run<runNumber>.dat`:

```ini
[GlobalTiming]
shift_ns <value>
```

The complete correction chain is therefore:

```text
raw TDC -> pixel offset -> ECal correction -> time walk -> run-global shift
```

## Output files

For a single cross-target run, the delta-t method normally writes:

```text
CDet_calibration_dt.dat
CDet_pixel_timing_fit_results.dat
```

For combined cross-target groups, it writes:

```text
CDet_calibration_dt_group<groupIndex>.dat
CDet_pixel_timing_fit_results_group<groupIndex>.dat
```

The master calibration file contains:

```ini
[PixelOffsets]
pixelID additive_offset_ns fit_entries

[ECalTiming]
p0 value
p1 value

[TimeWalk]
p1_L1 value
p1_L2 value
totref_L1 value
totref_L2 value
totmin value
totmax value
```

The detailed fit-results file records each pixel's fit status, centroid,
uncertainty, sigma, residual correction, total offset, detector reference,
chi2/NDF, and rejection reason.

## Older direct-LE method

The macro still contains an older method in `plotAllTDC(true)`. That method
fits the direct CDet leading-edge peak for each pixel and calculates

```text
c_i = mean_global_LE - mean_pixel_LE
```

It writes the legacy `CDet_calibration.dat`. The all-cross workflow instead
uses the ECal-CDet delta-t method and `CDet_calibration_dt*.dat`. Calibration
products from the two methods should not be mixed.

## Driver argument correction

The two-pass cross-target drivers previously contained stale positional calls
of the form:

```cpp
plotCDetLayersTimeComp(true, 1.0, -15, 15, ...)
```

The current function signature starts with:

```cpp
plotCDetLayersTimeComp(bool overwrite, int pixelBase, double Width, ...)
```

That stale call was interpreted as approximately:

```text
pixelBase = 1
Width = -15
```

The negative histogram width failed validation, stopping the automated
calibration sequence at its first stage-3 ECal fit. The affected combined-run,
individual-run, non-group two-pass, and basic drivers have been corrected by
inserting `pixelBase = 416` before the histogram width. Their stage-3 calls now
use the current 10--35 ns ECal ADC timing range and a -60--30 ns ECal-minus-CDet
timing range instead of the obsolete 95--125 ns and -104 to -60 ns ranges.

## Historical artifact caution

The checked-in `CDet_pixel_timing_fit_results.dat` was generated with an older
positive delta-t convention/range, including a 70--120 ns fit interval. It is
useful evidence of the algorithm's previous operation, but it is not evidence
of a calibration generated with the current negative default ranges.

## Hierarchical production method and standalone diagnostic

`extractHierarchicalCDetPixelTimingOffsets()` is the production method used by
the two full cross-target drivers. The underlying
`extractHierarchicalCDetPixelTimingOffsetsDiagnostic()` can still be run
without activating its candidate calibration for focused studies. The method divides each 16-pixel
bar into contiguous groups 0--7 and 8--15, sums the actual counts from all
eight pixels, and requires at least 100 total entries in the -30 to +10 ns
group-fit interval. Individual pixels are not normalized before summing. An
individual fit requires at least 35 entries in that same broad interval and a
fitted signal significance of at least 1.5.

The valid group centroid defines an adaptive individual fit window with a
half-width between 4 and 8 ns. It does not require the fitted individual
centroid to be within 2 ns of the group centroid. Individual fit failures use
the group centroid and are labeled `group_fallback`; other results are labeled
`individual_fit`, `retained_existing`, `unavailable`, or `unused`.

The group fit's chi-square per degree of freedom is retained as a diagnostic,
but it does not invalidate the entire group. A raw-count sum can be
non-Gaussian because its pixels have different offsets or cross-talk peaks;
the group fit is used only to provide a starting centroid for the constrained
individual fits. The chi-square quality requirement remains in force for each
individual fit, with a maximum chi-square per degree of freedom of 15.

If an eight-pixel group fit is invalid, pixels with at least 35 entries fall
back to independent fits over the full -30 to +10 ns interval. These fits use
the same individual quality requirements and are labeled
`individual_broad_fallback`. A pixel retains its existing offset if both its
group fit and broad individual fit are unusable.

The standalone diagnostic never activates its candidate calibration. A
complete stage-0 run for the currently selected problem bars can be launched with:

```cpp
.x Run_CDet_Hierarchical_Timing_Diagnostic.C
```

Its principal outputs are:

```text
CDet_calibration_dt_hierarchical_stage0_candidate.dat
CDet_pixel_timing_fit_results_hierarchical_stage0.dat
CDet_pixel_timing_hierarchical_stage0.root
tdcPlots/hierarchical_stage0/
```

In the comparison plots, the dashed blue curve is the old broad individual
fit, the red curve is the constrained individual fit, and the green line is
the eight-pixel group centroid. The ROOT diagnostic file stores group
histograms under `groups/`, pixel histograms under `pixels/`, and self-contained
copies of the selected-bar histograms with their fitted curves under
`fit_overlays/barNNN/`. Draw one of those fitted histograms directly; for
example:

```cpp
TFile *f = TFile::Open("CDet_pixel_timing_hierarchical_stage0.root");
TH1D *h = f->Get<TH1D>("fit_overlays/bar029/pixel_0472_with_fits");
h->Draw();
```

Keep the file open while using the histogram. The generated PDF and PNG retain
the complete 4-by-4 bar view, including the green group-centroid markers.
With the final diagnostic settings, the run-5710 stage-0 study found 879
accepted constrained individual fits, 6 accepted broad individual fallbacks,
715 group fallbacks, and 1,088 retained existing or unused pixels. The
detector reference was -8.893 ns. The interval was refined to -30 to +10 ns to
include bar 29 group 58's dominant timing structure without as much later-time
background.

During a full workflow, the three production invocations write the active
calibration file and update the in-memory offsets before the next stage. Their
auditable products are kept separately as:

```text
CDet_pixel_timing_fit_results_hierarchical_pass1_offsets.dat
CDet_pixel_timing_hierarchical_pass1_offsets.root
tdcPlots/hierarchical_pass1_offsets/

CDet_pixel_timing_fit_results_hierarchical_pass2_offsets.dat
CDet_pixel_timing_hierarchical_pass2_offsets.root
tdcPlots/hierarchical_pass2_offsets/

CDet_pixel_timing_fit_results_hierarchical_final_closure.dat
CDet_pixel_timing_hierarchical_final_closure.root
tdcPlots/hierarchical_final_closure/
```

The combined-run driver prefixes these output tags with `groupN_` so products
from different ECal cross-target groups do not overwrite one another.
