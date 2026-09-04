# Run 5710 half-bar-aligned final calibration

This is the definitive run-5710 CDet timing calibration after the common
unbinned ECal-slope correction, hierarchical pixel-offset extraction, and
half-bar intercept alignment.

## Final constants and closure

```text
[ECalTiming]
p0 74.192800
p1 0.817261
```

Final common within-half-bar slopes:

| Sample | Residual slope (ns/ns) |
|---|---:|
| Layer 1 | +0.000743058 +/- 0.00247422 |
| Layer 2 | -0.000738382 +/- 0.00221575 |
| Combined | +0.000000028 +/- 0.00166039 |

Final half-bar intercept closure at `t_ECal = 22 ns`:

```text
validated half-bars                    123 / 168
detector reference                     29.6304 ns
within-half-bar residual RMS            2.99597 ns
between-half-bar intercept RMS          0.03961 ns
total individual-hit RMS                2.99623 ns
predicted RMS after another alignment   2.99597 ns
```

The between-half-bar RMS was reduced from `1.61545 ns` to `0.03961 ns`.
The paired Layer-1/Layer-2 mean-time distribution has an observed width of
approximately `2.27 ns`. Its width is smaller than the individual-hit RMS
because the plotted time is the average of two layer measurements.

The remaining proposed alignment would improve the calculated individual-hit
RMS by only `0.00026 ns`; it must not be applied. The calibration is closed.

## Contents

- `CDet_calibration_dt_run5710_halfbar_aligned_final.dat`: definitive constants.
- `CDet_halfbar_intercept_corrections_applied.dat`: half-bar increments that
  produced the alignment.
- `CDet_halfbar_intercept_final_closure.dat`: final diagnostic-only proposal,
  retained to demonstrate that no further iteration is warranted.
- `CDet_pixel_quality_cuts_run5710_halfbar_aligned_final.root`: the matched 102
  manual LE-versus-TOT polygon cuts.
- `CDet_pixel_timing_fit_results_run5710_fixed_effects_final.dat` and
  `CDet_pixel_timing_run5710_fixed_effects_final.root`: the immediately
  preceding final hierarchical pixel-offset extraction.
- `CDet_ECal_fixed_effects_final_closure.jpg`: unbinned fixed-effects closure.
- `CDet_halfbar_intercept_final_closure.jpg`: post-alignment intercept closure.
- `CDet_pooled_timing_final.jpg` and `CDet_timing_components_final.jpg`: final
  pooled and component timing diagnostics.

The diagnostic ROOT writer was corrected after these results were produced;
the two original ~500-byte empty ROOT outputs are deliberately excluded. The
tables and rendered canvases above are valid. Future invocations write the
diagnostic objects into the ROOT file correctly.

## Primary SHA-256 checksums

```text
7256494534304a07bf0545c21c00f1da5aadb0898444556d1ece4778ffc4cade  CDet_calibration_dt_run5710_halfbar_aligned_final.dat
aed6abc0e2a4b427f25cd2bb2fde437f94e881fdbbd95a8cc212b0ee1447dfaa  CDet_halfbar_intercept_corrections_applied.dat
ee007129bedd67da85b112c4adec0086d988f8f7fcd4dec5f3bb7d12454ca4c7  CDet_halfbar_intercept_final_closure.dat
fcb05430884ce3e928e5712feececdfd47b47058ad1e8e4a9bb7ed9b33346cd6  CDet_pixel_quality_cuts_run5710_halfbar_aligned_final.root
fc752845399039e410c18c732b4dba06f3c0808b9ca8227fbef73a9c031c9b4b  CDet_pixel_timing_fit_results_run5710_fixed_effects_final.dat
2052da6999c88d905c3432e3dc2d7859b8cb94ef9443fef1ab38c66e7d4f72b8  CDet_pixel_timing_run5710_fixed_effects_final.root
ff0c4d7f1d8e0ad29beaebb0b17981f32ccc83c8aa143f13c1d9f2c796928f37  CDet_ECal_fixed_effects_final_closure.jpg
cc4a19ab74426cb0cf0e4e12b57937b0f57cae14852a0393d8dfa623f34d3cdd  CDet_halfbar_intercept_final_closure.jpg
c02392371617f2292e7c2334ab0c4ecbd7d47cc42e6d25f652da910454890662  CDet_pooled_timing_final.jpg
d476aaf05eca233f6be8710ef9c51e34e816f5b8a04e6783b7463d2d61431532  CDet_timing_components_final.jpg
```
