# CDet Run 5992 Timing Calibration Record

## Status

Run 5992 was conditionally accepted on 2026-09-11 for run-specific CDet timing
calibration. The accepted Run 5710 detector-wide master calibration remains
fixed. Run 5992 retains the master ECal slope; only the absolute timing origin
is run-dependent.

The authoritative selection and display settings are in
`CDet_run5992_projection.conf`. The resulting run constants are:

```ini
[ECalTiming]
p1 0.810203

[GlobalTiming]
shift_ns 5.171528
```

## Reproducible procedure

For an exact clean reproduction of both Run 5710 and Run 5992, preserve and
move `CDet_calibration_dt.dat`, `CDet_run5710.dat`, and `CDet_run5992.dat` out
of `scripts/cdet`. In a fresh ROOT process started from that directory, run:

```cpp
.L Run_CDet_Reproduce_5710_5992_FromScratch.C+
Run_CDet_Reproduce_5710_5992_FromScratch()
```

This is the authoritative student procedure. It regenerates and verifies the
Run 5710 master before seeding Run 5992 with that generated master `p1`; it
then fits only Run 5992 `shift_ns` and verifies the final rounded constants.

For a shift-only rerun when an already verified Run 5710 master is present:

Preserve any existing `CDet_run5992.dat`, retain the accepted
`CDet_calibration_dt.dat`, set or inherit its master `p1 = 0.810203`, and run
the shift-only calibration in a fresh ROOT session:

```cpp
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

- Fixed Run 5710 master `p1`: `0.810203 ns/ns`
- Residual fixed-effects slope: `-0.0131452 +/- 0.0152886 ns/ns`
- Run-specific `shift_ns`: `5.171528 ns`
- Final Gaussian-core centroid: `30.0000 +/- 0.0888 ns`
- Accepted pairs before/after the shift: `1453 / 1453`
- Final `plotCDetLayersTimeComp` visual review: conditionally accepted

An exploratory direct refit at `shift_ns = 5.202920 ns` produced
`p1 = 0.820628`, and isolated fits at shifts of 4.203 and 6.203 ns produced
`p1 = 0.819331` and `0.831321`, respectively. These differences are smaller
than the individual statistical uncertainties and arise from shift-sensitive
event selection in the limited Run 5992 sample. They do not establish a
physical change in the detector-wide ECal slope. The commissioned choice is
therefore the high-statistics Run 5710 `p1` with only `shift_ns` fitted.

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
