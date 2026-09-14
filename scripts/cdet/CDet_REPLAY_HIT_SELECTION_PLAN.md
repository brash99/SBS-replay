# Plan: ECal-associated CDet hit selection during replay

Status: proposed implementation; no SBS-offline selection changes implemented.

## Objective

For LH2 events, approximately 800 CDet TDC hits contain only about 3–4 hits
associated with the primary electron defining the ECal shower cluster. These
occupancies are the working estimates supplied for this plan, not measured
performance results.

Move the calibrated hit-selection procedure from the analysis macro into
SBS-offline so replay ROOT files contain an ECal-associated CDet candidate
collection. Aim to retain the electron hits in a manageable collection, ideally
around 8–10 candidates per event. Treat that multiplicity as a performance
target, not a fixed output count or an unconditional top-ten truncation.

Preserve existing hit collections initially. Add candidate information and
selection diagnostics so efficiency and background rejection can be evaluated.

## Current implementation and implications

- `SBSGenericDetector::FindGoodHit()` chooses one TDC pulse per channel by
  proximity to `tdc.GoodTimeCut`. This happens without ECal association. With
  multiple pulses in one channel, an accidental pulse can displace the electron
  pulse from the generic good-hit collection.
- `SBSCDet::CoarseProcess()` constructs hits from that good-hit collection using
  broad LE and ToT cuts. Its cut values are currently initialized in the
  constructor. `SBSCDet::FineProcess()` is empty.
- The analysis macro applies pixel offsets, an ECal-dependent correction,
  layer-dependent ToT time-walk corrections and a run-dependent time shift.
  It also applies ECal-projection and other preselection requirements.
- The event display's “best hits” mode adds a window on ECal time minus corrected
  CDet time. Its peak mean and width may be supplied or fitted from the display
  sample. That sample is conditioned on accepted pairs and a selected bar.
- `SBSGEPEArm::CoarseReconstruct()` already accesses ECal following the base
  coarse reconstruction, providing a potential explicit integration point.

The replay association must examine all valid decoded LE/TE pulse pairs, not
only the generic one-pulse-per-channel good-hit collection. Preserve channel
and pulse indices throughout processing. This does not require changing the
generic good-hit behavior for other detectors.

## Event inputs and database responsibilities

The database stores calibration constants, geometry and selection policy.
Event-specific ECal information must come from the reconstructed detector:

- selected/main cluster index and validity;
- `earm.ecal.adctime`, the energy-weighted mean ADC time of the main cluster;
- cluster position and the geometry needed to project it to CDet;
- cluster energy or other quality information if used by the chosen eligibility
  criteria.

Use the same ECal cluster consistently for time, position and energy. Do not
substitute `earm.ecal.a_time` (block-level pulse times) or `earm.ecal.atimeblk`
(highest-energy block time) for `adctime`.

An absent or invalid ECal cluster must produce an explicit association status,
not an implicit zero-time or zero-position match. HCal is not required for the
initial electron-associated CDet selection; no shared ECal/HCal cut is assumed.

## Proposed earm.cdet database content

The names below are conceptual groups, not implemented database keys. Final
spelling, vector dimensions and defaults must be settled during implementation.
Reuse existing decoder and geometry keys where their meanings match.

| Information | Granularity | Purpose |
| --- | --- | --- |
| Channel mapping, physical layer and position | Per channel | Map pulses to physical detector elements |
| Active/bad-channel status | Per channel | Identify unusable channels |
| TDC conversion and reference-time convention | Per channel/reference as appropriate | Define base times in ns |
| Pixel timing offset and calibration-valid flag | Per pixel | Align times and distinguish uncalibrated from measured-zero offsets |
| Time-walk coefficient and reference ToT | Per physical layer initially | Apply the calibrated ToT correction |
| Time-walk validity range and out-of-range policy | Per layer or calibration set | Make correction behavior explicit |
| ECal correlation coefficients p0, p1 and fixed delta | Detector-wide with validity-dependent overrides | Reproduce the calibrated ECal dependence |
| Global timing-origin shift | Run/validity period | Align the run's timing origin |
| Coordinate alignment corrections | Per layer | Reproduce geometry used in ECal matching |
| Broad LE and ToT acceptance | Global or per layer | Reject unusable pulses |
| ECal-projection residual limits or resolutions | Per layer initially | Select the spatial region consistent with the electron |
| ECal–CDet residual mean, width and acceptance multiplier | Global initially; per layer if validated | Define deterministic best-hit timing acceptance |
| Optional cross-layer spatial/time consistency limits | Per layer pair | Identify coherent candidate groups |
| Selection enablement and calibration provenance/version | Per configuration/validity period | Support controlled rollout and reproducibility |

Existing `tdc.calib`, `tdc.offset`, reference-TDC settings and `xpos/ypos/zpos`
are starting points. Verify their units, signs and application stage before
reusing them for exported macro constants. In particular, the macro multiplies
input TDC times by 0.01 and adds its pixel offsets; blindly copying these into
decoder fields could change units, reverse a sign or apply a correction twice.

Validate vector lengths against the actual detector mapping, including any
reference channels. The macro's physical layer conventions must be mapped
explicitly; do not assume database row/layer indices have identical meanings.

Keep calibration validity separate from channel health. An uncalibrated channel
requires an explicit handling policy rather than silently treating its offset
as a valid zero or automatically declaring the channel bad.

## Timing convention to preserve

For a valid pulse i in physical layer l, reproduce the macro's correction chain:

```text
t_corr[i] = t_base[i] + offset[i]
            - (p0 + p1 * t_ecal) + delta
            - timewalk[layer](ToT[i]) + run_shift

timewalk[layer](ToT) = a[layer] *
    (1 / sqrt(ToT) - 1 / sqrt(ToT_ref[layer]))

residual[i] = t_ecal - t_corr[i]
accept_time[i] = abs(residual[i] - peak_mean) <= nsigma * peak_sigma
```

`t_base` includes the agreed TDC conversion and reference treatment. Document
all units, signs, enable flags and the order of any broad cuts relative to
corrections. Apply common timing shifts consistently to LE and TE so ToT is
not unintentionally changed. Handle invalid/nonpositive ToT explicitly.

Because `t_corr` already depends on ECal time, the residual calibration must
belong to this exact convention. Use fixed, validated constants during replay;
do not fit or recenter from the events processed in that replay invocation.

The current Run 6077 artifacts provide:

- master pixel offsets, ECal p0/delta and layer time-walk parameters in
  `CDet_calibration_dt.dat`;
- a run-specific p1 of 0.002520 and additive shift of 2.913449 ns in
  `CDet_run6077.dat`;
- an analysis ToT window of 8–35 ns and x-projection half-width of 0.08 m in
  `CDet_run6077_projection.conf`.

These describe the current local configuration, not universal defaults. The
remaining calibration input for deterministic best-hit selection is an approved
residual peak mean and width, validated beyond the selected-bar display sample.

## Reconstruction sequence

1. Decode and retain all valid CDet pulse pairs with stable channel/pulse indices.
2. Apply detector-local calibration and broad pulse-quality requirements.
3. Reconstruct and select the ECal cluster.
4. Obtain a consistent cluster time, position, index and validity state.
5. Apply ECal-dependent CDet corrections and spatial/timing association.
6. Optionally evaluate cross-layer consistency without prematurely removing
   isolated valid candidates.
7. Write candidates, residuals, selection flags and event association status.

Keep detector-specific calibration and association logic in `SBSCDet` or a
small CDet-owned helper. Decide whether `FineProcess()` can safely invoke it
after the required ECal cluster is finalized, or whether `SBSGEPEArm` should
call it explicitly after ECal reconstruction. Verify actual call order and
cluster-selection timing before choosing; detector registration order alone
should not be an implicit dependency.

Use the established projection model first, with documented coordinate frames,
layer alignment and origin/vertex assumption. Any improved vertex-aware
projection is a separate, validated change.

## Proposed ROOT output

Preserve existing branch meanings. Add a clearly named associated-candidate
collection containing:

- original channel ID and pulse index;
- physical layer and position;
- corrected LE/TE time and ToT;
- ECal spatial residuals and timing residual;
- associated ECal cluster index;
- selection/quality flags and calibration-valid status;
- optional group identifier and consistency metrics when grouping is implemented.

Add event-level counts at the broad-quality, spatial and timing selection
stages, plus association status. Record effective calibration/configuration
provenance in suitable run metadata. Retain enough information to distinguish
no candidates from missing ECal input or invalid calibration.

## Implementation stages

1. **Freeze conventions and database schema.** Trace units, references, physical
   channel mapping, correction order and reconstruction ordering. Specify
   calibration validity and missing-input policies. Establish the fixed residual
   mean/width and a reference sample.
2. **Port calibrated times.** Read validated database constants and calculate
   corrected CDet times while preserving current outputs. Compare against the
   macro for the same pulses before adding rejection.
3. **Add ECal association.** Use all valid pulses, apply the established spatial
   cuts and fixed residual timing window, and export candidates and diagnostics.
4. **Validate on LH2.** Measure retained multiplicities, signal retention and
   accidental contamination. Resolve disagreements with the macro explicitly,
   including expected differences from retaining multiple pulses per channel.
5. **Evaluate grouping.** Add cross-layer consistency or ranking only when it
   demonstrably improves selection without unacceptable signal loss.

Do not transfer display-only settings such as selected bar, plot ranges or
histogram binning into replay physics selection. Event eligibility and optional
energy cuts must be explicit rather than inherited accidentally from a display
sample.

## Validation and acceptance criteria

- **Calibration parity:** corrected times and residuals agree with the macro for
  identical pulses and constants within a defined numerical tolerance.
- **Multiple-pulse recovery:** an electron-compatible pulse remains eligible when
  another pulse in its channel is closer to the generic good-time setting.
- **Determinism:** splitting a run into files or changing event-processing order
  or count does not change an event's calibrated result.
- **Geometry correctness:** confirm layer mapping, coordinate transformations,
  projection signs and boundary behavior using known examples.
- **Boundary behavior:** test missing ECal clusters, invalid calibration,
  nonpositive ToT, bad channels, empty events and ambiguous candidates.
- **Output compatibility:** existing branches retain their meanings; new
  candidate indices map back to the original pulses; event state clears fully.
- **LH2 performance:** report the candidate-count distribution, including tails,
  and signal-retention/background estimates on a representative sample. Define
  the truth or independent reference method before claiming efficiency; the
  display's selected hits alone are not ground truth.
- **Practical replay cost:** measure additional CPU and output size on realistic
  high-occupancy events.

Numerical efficiency and contamination acceptance thresholds remain to be
agreed. Do not claim the 3–4 electron hits are retained merely because the
candidate collection averages 8–10 entries.

## Source references

- [SBSCDet implementation](../../../SBS-offline/SBSCDet.cxx) — broad hit
  construction, database reading and empty fine processing.
- [Generic detector implementation](../../../SBS-offline/SBSGenericDetector.cxx)
  — pulse selection and decoder calibration settings.
- [Electron-arm reconstruction](../../../SBS-offline/SBSGEPEArm.cxx) — ECal
  access following coarse reconstruction.
- [CDet database](../../DB/db_earm.cdet.dat) — current mapping and geometry.
- [Calibration and event-display macro](PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C)
  — correction chain, projection cuts and best-hit residual selection.
- [Master calibration](CDet_calibration_dt.dat),
  [Run 6077 timing constants](CDet_run6077.dat), and
  [Run 6077 analysis/display configuration](CDet_run6077_projection.conf).
