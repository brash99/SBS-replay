# CDet Run 5992 Timing Calibration Record

## Status

Run 5992 was conditionally accepted on 2026-09-11 for run-specific CDet timing
calibration. The accepted Run 5710 detector-wide master calibration remains
fixed. Only the Run 5992 ECal slope and absolute timing origin were fitted.

The authoritative selection and display settings are in
`CDet_run5992_projection.conf`. The resulting run constants are:

```ini
[ECalTiming]
p1 0.822678

[GlobalTiming]
shift_ns 4.840388
```

## Reproducible procedure

Preserve and remove an existing `CDet_run5992.dat`, retain the accepted
`CDet_calibration_dt.dat`, and run in a fresh ROOT session:

```cpp
.L Run_CDet_Calibrate_RunSpecificP1.C+
Run_CDet_Calibrate_RunSpecificP1(5992)

.L Run_CDet_Calibrate_RunTimingShift.C+
Run_CDet_Calibrate_RunTimingShift(5992)
```

For the final human review, start another fresh ROOT session and use the same
configuration for both calls:

```cpp
.L PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C+
PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget(
    "CDet_run5992_projection.conf");
plotCDetLayersTimeComp("CDet_run5992_projection.conf");
```

Do not load the event sample through the legacy positional run-number overload
and then plot with the configuration overload. The latter cannot rebuild an
already selected sample.

## Closure results

- Run-specific `p1`: `0.822678 ns/ns`
- Residual fixed-effects slope: `0.003155 +/- 0.012956 ns/ns`
- Run-specific `shift_ns`: `4.840388 ns`
- Final Gaussian-core centroid: `30.0325 +/- 0.0779 ns`
- Accepted pairs before/after the shift: `1856 / 1856`
- Final `plotCDetLayersTimeComp` visual review: conditionally accepted

Saved per-pixel LE-versus-ToT polygons are not applied by either run-specific
calibration driver or by `plotCDetLayersTimeComp`. They remain part of the
detector-wide pixel-offset calibration and review procedure.

## Run-quality limitation

Run 5992 has 1,132,123 events, but only 196,161 contain a reconstructed ECal
cluster and only 49,403 have a reconstructed cluster inside the approved
`-15 < earm.ecal.adctime < 15 ns` window. Raw ECal ADC data are present in
essentially every event, so the dominant loss is upstream ECal reconstruction.
This is outside the scope of the CDet calibration.

Among in-time reconstructed ECal events, 4.18% produce a fully accepted CDet
event, compared with 9.03% for Run 5710. The remaining relative factor of 2.16
was traced to projected ECal-to-CDet geometry (approximately a factor 1.52,
predominantly the projected-x requirement) and two-layer coincidence
(approximately a factor 1.50). LE, ToT, and multiplicity cuts do not create
the discrepancy. Run 5992 is dominated by trigger bit 4096, whereas Run 5710
is dominated by trigger bit 2048.

The calibration is therefore accepted for use with an explicit low-statistics
and upstream-reconstruction qualification. It is not an independent
detector-wide calibration and must not replace the Run 5710 master constants.
