# Definitive run-5710 CDet timing calibration

This directory is the matched archive of the final run-5710 timing
calibration after correcting the pooled-regression/Simpson's-paradox problem.
The files were copied together after the final hierarchical pixel-offset pass
and its subsequent closure analysis. Do not overwrite them with later studies.

## Final calibration state

```text
[ECalTiming]
p0 74.192800
p1 0.812267
```

The pixel-quality file contains 102 manually selected LE-versus-TOT polygon
cuts. The final hierarchical extraction reported:

```text
ECal-selected events                  261148 / 264337
manual LE-versus-TOT cuts applied    102
robust detector reference            -7.66232 ns
manual-population individual fits    56
individual constrained fits          984
broad individual fallbacks           2
eight-pixel group fallbacks           551
retained existing/unavailable         1095 / 0
```

## Final ECal-slope closure

The calibration slope is determined with an unbinned half-bar fixed-effects
regression, not the pooled timing profile. After the final pixel-offset
extraction, the closure results were:

| Sample | Residual slope (ns/ns) |
|---|---:|
| Layer 1 | +0.000330716 +/- 0.00245523 |
| Layer 2 | -0.00131549 +/- 0.00220603 |
| Combined | -0.000491710 +/- 0.00165043 |

The combined residual is `-0.30` standard deviations from zero. The pooled
slope remained `+0.215670 ns/ns`; that is the between-half-bar association
exposed by Simpson's paradox and is not a residual calibration error.

## Contents

- `CDet_calibration_dt_run5710_fixed_effects_final.dat`: definitive constants.
- `CDet_pixel_quality_cuts_run5710_fixed_effects_final.root`: the matched 102
  polygon selections.
- `CDet_pixel_timing_fit_results_run5710_fixed_effects_final.dat`: complete
  per-pixel extraction table.
- `CDet_pixel_timing_run5710_fixed_effects_final.root`: hierarchical fit
  diagnostics.
- `CDet_ECal_fixed_effects_final_closure.jpg`: final within-half-bar closure.
- `CDet_ECal_pooled_diagnostics_final.jpg`: pooled/component diagnostic retained
  to illustrate why the pooled slope must not drive calibration.
- `CDet_halfbar_ecal_timing_run5710_fixed_effects_final.root`: all 168 final
  half-bar histograms, profiles, fits, and the four overview canvases.
- `CDet_halfbar_ecal_timing_run5710_fixed_effects_final_summary.dat`: final
  168-row half-bar fit summary.
- `CDet_run5710_fixed_effects_final_L1_bank0.jpg`
- `CDet_run5710_fixed_effects_final_L1_bank1.jpg`
- `CDet_run5710_fixed_effects_final_L2_bank0.jpg`
- `CDet_run5710_fixed_effects_final_L2_bank1.jpg`: the four final 7-by-6
  half-bar overview canvases. Unlike the pre-correction set, their fitted
  slopes no longer have an overwhelmingly negative sign.
- `tdcPlots/run5710_fixed_effects_final/`: selected hierarchical fit canvases in
  PDF and PNG form.

The four pre-correction 7-by-6 half-bar overview canvases and their source ROOT
file remain preserved in the sibling
`simpsons_paradox_run5710_before_correction/` archive. They are intentionally
not relabeled as final-calibration products.

## SHA-256 checksums for primary products

```text
e1e2370edae459ee36cbad2a5ebbb6e8ea5b9be5e168ad9e5e4c0d83134220b2  CDet_calibration_dt_run5710_fixed_effects_final.dat
fcb05430884ce3e928e5712feececdfd47b47058ad1e8e4a9bb7ed9b33346cd6  CDet_pixel_quality_cuts_run5710_fixed_effects_final.root
fc752845399039e410c18c732b4dba06f3c0808b9ca8227fbef73a9c031c9b4b  CDet_pixel_timing_fit_results_run5710_fixed_effects_final.dat
2052da6999c88d905c3432e3dc2d7859b8cb94ef9443fef1ab38c66e7d4f72b8  CDet_pixel_timing_run5710_fixed_effects_final.root
896b5280214938109a0850f45e8268afc6d7c089c5a175236e29052387a68439  CDet_ECal_fixed_effects_final_closure.jpg
70674ad8a756cdea6672b32dfd7836ce7debfeb4b3c11221ec53805abc8709f1  CDet_ECal_pooled_diagnostics_final.jpg
d71828f0cdd015b9b230f2debb312ad9b362e57e0ab3b1c4ab55c92267bfd3e3  CDet_halfbar_ecal_timing_run5710_fixed_effects_final.root
56896856a4aa8245e7d1f0ca1ecb1be1c9bc1c53cfe3d6ab9451cdc450fe3d98  CDet_halfbar_ecal_timing_run5710_fixed_effects_final_summary.dat
e0d8f713865fed0995215b0d3036cbfae92a9dfac1467401c4f427bfe3e58a0f  CDet_run5710_fixed_effects_final_L1_bank0.jpg
e1a2527af72e77b47d07c2f2ef0cf5e4a7d5d88d04f8a08cdb93b8c6edec4d3b  CDet_run5710_fixed_effects_final_L1_bank1.jpg
eb570f9776fdc840e15c8bb94b3d812b939f48729cb4532c4727072de57b29e4  CDet_run5710_fixed_effects_final_L2_bank0.jpg
636bac3dc44cc87e153e75db5e4e8ee2a77b08507e38fd7d62760b01b0ff9c5a  CDet_run5710_fixed_effects_final_L2_bank1.jpg
```
