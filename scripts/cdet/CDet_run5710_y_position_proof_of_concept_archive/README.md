# Run 5710 CDet y-position proof of concept

This archive records the first reconstruction of approximate CDet y-position
from corrected CDet leading-edge timing. ECal y was used only to train and
validate the position response; it was not applied as a CDet timing correction,
and `CDet_calibration_dt.dat` was not modified by this study.

## Method

The accepted paired hits were divided by layer, readout side, and the
geometrically equivalent section groups 1+6, 2+5, and 3+4. A fixed-effects
regression removed the mean ECal y and mean CDet time independently for each
half-bar before extracting twelve within-half-bar propagation slopes. The
separate position-calibration file stores those slopes and a timing intercept
for each half-bar with at least 100 accepted hits.

For a valid half-bar, CDet y is reconstructed from CDet timing alone as

```text
y_CDet = (t_CDet - intercept_halfbar) / slope_layer,side,group
```

## Principal results

- Layer 1 `y_CDet - y_ECal`: mean `-0.00365 m`, RMS `0.2133 m`.
- Layer 2 `y_CDet - y_ECal`: mean `-0.00486 m`, RMS `0.2102 m`.
- Layer difference `y_L1 - y_L2`: mean `-0.00132 m`, RMS `0.3525 m`.
- Under equal, independent layer errors, `0.3525/sqrt(2) = 0.249 m` per layer.

The agreement of the two ECal-referenced widths and the independent positive
Layer-1/Layer-2 correlation establishes a real CDet timing-based y response.
The quoted RMS values include broad non-Gaussian tails and should not yet be
interpreted as intrinsic Gaussian detector resolutions.

## Contents

- `CDet_time_vs_ECal_y_by_geometry.jpg`: twelve physical-coordinate timing
  plots split by layer, side, and section group.
- `CDet_ECal_y_fixed_effects.jpg`: twelve half-bar-centered distributions and
  fixed-effects propagation slopes.
- `CDet_y_position_calibration_run5710.dat`: independent TEnv-style CDet
  y-position calibration parameters.
- `CDet_y_position_ECal_validation.jpg`: reconstructed CDet y versus ECal y and
  residuals for both layers.
- `CDet_y_position_layer_comparison.jpg`: CDet-only layer-to-layer correlation
  and y-difference distribution.

## SHA-256 checksums

```text
9692fc07c3adf6c13e6d1fce005937c5fbb5ee73189daa77ae280c125fb8ccba  CDet_ECal_y_fixed_effects.jpg
006e5580a9f9a766cd29145484342263ca6fffa9680df568ca355e0c40db53be  CDet_time_vs_ECal_y_by_geometry.jpg
8043192685db9f027a055f19599a2c12761d68769b65e2c43a870f15ca9bed20  CDet_y_position_ECal_validation.jpg
8072ae3c813ae60b777c67241a5c2d26b9a472fb5c7a56a525438d552f8bec77  CDet_y_position_calibration_run5710.dat
197111149c0fb0c91df21bd926cd709da758b843b3b49712b31358d62a9337e5  CDet_y_position_layer_comparison.jpg
```
