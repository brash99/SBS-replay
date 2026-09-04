# Run 5710 CDet--ECal timing: Simpson's-paradox snapshot

This directory preserves the CDet-versus-ECal timing state **before correcting
the newly identified within-half-bar dependence**. It is intentionally a
historical snapshot: do not replace these files when the calibration is rerun.

## Analysis state

- Run: 5710, full available event sample.
- Pixel timing calibration: run-5710 optimized v2
  (`CDet_calibration_dt_run5710_optimized_v2.dat`).
- ECal timing parameters applied: `p0 = 74.192800 ns` and
  `p1 = 1.120541 ns/ns`.
- Pixel-quality selection: 102 manual LE-versus-TOT polygons.
- Geometry: 84 half-bars per layer, arranged as six half-modules containing
  two banks of seven half-bars.

The ROOT file contains the 168 underlying two-dimensional histograms,
profiles, fits, and the four overview canvases. The text file contains one
fit-summary row per half-bar.

## Observed reversal

Of the 168 half-bars, 135 had successful profile fits:

- 130 of 135 fitted slopes were negative.
- 122 were negative by more than two standard errors.
- None were positive by more than two standard errors.
- The unweighted mean of the 135 fitted slopes was
  `-0.128439 ns/ns`.

Nevertheless, the profile formed after pooling hits from all half-bars was
nearly flat. A covariance decomposition of the binned ROOT histograms makes
the mechanism explicit:

| Layer | Pooled event-regression slope | Within-half-bar slope | Between-half-bar slope |
|---|---:|---:|---:|
| 1 | +0.014331 | -0.254412 | +0.922884 |
| 2 | +0.011965 | -0.255751 | +0.910988 |

For Layer 1, the within-half-bar covariance was `-582674`, while the
between-half-bar covariance was `+625204`, leaving total covariance
`+42529.8`. For Layer 2 the corresponding values were `-584841`, `+620350`,
and `+35509.3`. Thus a strong negative association within detector groups is
almost cancelled by a positive association among their group means. This is
a direct detector-calibration example of Simpson's paradox.

The pooled profile-fit slope shown in the diagnostic canvas is closer to zero
than the pooled event-regression slopes above because ROOT's profile fit uses
bin uncertainties as weights. That numerical distinction does not alter the
reversal.

## Preserved files

- `CDet_halfbar_ecal_timing_run5710_before_ecal_slope_correction.root`:
  self-contained histograms, profiles, fits, and four overview canvases.
- `CDet_halfbar_ecal_timing_run5710_before_ecal_slope_correction_summary.dat`:
  the 168 half-bar fit records.
- `layer1_bank0_halfbar_ecal_timing.jpg`
- `layer1_bank1_halfbar_ecal_timing.jpg`
- `layer2_bank0_halfbar_ecal_timing.jpg`
- `layer2_bank1_halfbar_ecal_timing.jpg`
- `pooled_and_selected_ecal_timing_diagnostics.jpg`: the pooled and selected
  half-bar comparison that initially exposed the contradiction.

## SHA-256 checksums

```text
72dad85b493960bed63da9b9e9dcfeae7fa41d0edcd9c0c34d1c3222433aedad  CDet_halfbar_ecal_timing_run5710_before_ecal_slope_correction.root
99a9a872a48a6b321e646cea4307ad2f5049c1322b669c910f212a6d5dcfbf47  CDet_halfbar_ecal_timing_run5710_before_ecal_slope_correction_summary.dat
c9d2bd045a7169338632d40e764f4458eefd0a9da0354949c54a22e0da078c5f  layer1_bank0_halfbar_ecal_timing.jpg
7785a52d740a56b51c3fcbfba84f4ba80e86e8e5ae386c2b151b4d5e0188ee05  layer1_bank1_halfbar_ecal_timing.jpg
16fa9792ccc770ced83710fa9708a4a4fb8a908feadaf8c7b624990d3f970721  layer2_bank0_halfbar_ecal_timing.jpg
5226f797e5f780c5c62174e435d849fb2d390f50f635c6e4d28440259ecc84b4  layer2_bank1_halfbar_ecal_timing.jpg
76c5c9f6e9f6238d6da831835596dac8cbda2dfd363234e0c09cfacdc1cbf81e  pooled_and_selected_ecal_timing_diagnostics.jpg
```

## Calibration consequence (not yet applied in this snapshot)

The appropriate ECal correction must be estimated from variation *within*
half-bars, rather than from the pooled profile. The two layers independently
suggest a residual slope near `-0.255 ns/ns`, implying a provisional common
`p1` near `1.120541 - 0.255 = 0.8655 ns/ns`. No such correction has been
applied to any artifact in this directory.
