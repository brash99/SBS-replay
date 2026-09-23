# Plan: ECal-associated CDet hit selection during replay

Status: implementation and validation in progress on 2026-09-22. The all-pulse
input, identity preservation, database-backed calibrated timing, diagnostic
ECal projection, and one-to-one cross-layer pairing stages are implemented.
Run 5711 studies now also define and validate a macro-level trajectory-time
ellipse for producing candidates for a future global ROI tracker.

## Implemented first slice

`SBSCDet` now has an opt-in `pulse.*` output collection built from every decoded
channel's complete `SBSData::TDC::GetAllHits()` store. An analysis pulse is
required to contain a nonzero LE, a nonzero TE, and a positive derived ToT;
incomplete or invalid decoder slots are never copied into the candidate arrays.
Separate event-level counters record LE-only, TE-only, and invalid two-edge
slots for decoder diagnostics. The collection adds the stable
`(pmtnum, index-within-channel)` identity, row/column/layer, channel position,
calibrated LE/TE/ToT, and raw decoder values. The existing `earm.cdet.hit.*` and
generic `earm.cdet.hits.*` meanings are unchanged.

`replay_CDet.C` enables this collection explicitly with
`SetStorePulseCandidates(kTRUE)`; the `SBSCDet` class default is disabled for
backward compatibility with other replay configurations.

A 2,000-event Run 5710 investigation found 8,082 decoder slots: 5,696 complete
triplets, 891 LE-only slots, and 1,495 TE-only slots. Although 912 events had a
channel with more than one decoder slot, only 130 events had a channel with two
or more complete triplets. Those 130 events contained 136 such channel
occurrences. Of 137 adjacent complete-triplet pairs, 86 overlapped in time and
51 did not; the overlapping population requires decoder-level investigation.
This result is why incomplete slots are retained only as diagnostic counts and
why the stable within-channel index is preserved for accepted complete
triplets.

The detailed results and decoder association mechanism are documented in
`CDet_MULTIPULSE_INVESTIGATION.md`. CDet now uses an alternating-edge state
machine: the first LE opens a pulse, repeated LEs are ignored until the first TE
closes it, and repeated TEs are then ignored until the next LE. Thus
`LE,LE,TE,TE` produces one `(first LE, first TE)` pulse. Original LE and TE
decoder-slot indices remain in the output for auditability.

## Implemented second slice: calibrated times without selection

For Run 5710, the standard analyzer database file `DB/db_earm.cdet.dat`
contains the 2,688 pixel offsets, TDC conversion, ECal p0/p1/delta parameters,
layer-dependent time-walk parameters, and global run shift. `SBSCDet` loads
these keys through its normal `ReadDatabase()` method. Neither `SBSCDet` nor
`replay_CDet.C` opens a macro calibration file.

The first database validity interval begins at Run 5710's recorded start
timestamp, `2025-08-07 01:26:21`. At Run 5711's start timestamp,
`2025-08-07 02:06:54`, a minimal run-validity block overrides only `p1`,
`shift_ns`, and the ECal timing window. Pixel offsets and every other common
calibration, geometry, pulse-quality, and pairing parameter are inherited.
Subsequent runs should follow the same timestamp-valid override pattern.

After both CDet and ECal have completed coarse processing,
`SBSGEPEArm::CoarseReconstruct()` passes the selected ECal cluster's
energy-weighted mean ADC time and cluster index to `SBSCDet`. Each complete
physical-pixel pulse then receives the same ordered corrections as the macro:
the 0.01 ns/count conversion, pixel offset, ECal p0/p1/delta correction,
layer-dependent ToT time walk, and Run 5710 global shift. No pulse is rejected
by this operation. Corrected LE/TE, ToT in ns, ECal-minus-CDet residual, and a
per-pulse calibration-valid flag are written under `earm.cdet.pulse.*`.
Event-level `earm.cdet.timing.*` fields record status, ECal cluster index, and
the ECal time used. A missing ECal cluster leaves corrected quantities invalid
rather than substituting a zero time.

A fresh 2,000-event Run 5710 replay contained 5,589 physical-pixel pulses.
There were 1,989 events with an applied ECal calibration and 11 explicitly
marked as missing ECal; those 11 events contained 23 pulses. All pulse-array
sizes and validity states were consistent. An independent calculation from the
original source calibration files reproduced every database-backed corrected
LE with a maximum absolute difference of `1.42e-14 ns`; ToT and residual
differences were zero at the stored precision.

## Implemented third slice: diagnostic pulse and projection decisions

Run 5710's database block now also contains the macro's broad pulse and initial
ECal-projection parameters. `SBSCDet` evaluates these decisions for every
complete physical-pixel pulse but does not remove, rank, or pair any pulse.
The new parallel outputs are:

- `pulse.broad_quality_pass`, using raw LE in 0.02--60 ns and ToT in 4--30 ns;
- `pulse.ecal_eligible`, using the selected cluster's 10--35 ns time window and
  the macro's finite/nonzero x/y acceptance;
- `pulse.x_corr`, `pulse.ecal_x_proj`, and `pulse.ecal_y_proj`;
- `pulse.ecal_x_residual` and `pulse.ecal_y_residual`;
- `pulse.spatial_pass`, using `|x residual| <= 0.08 m` and
  `|y residual| <= 0.36 m` after the configured alignment offsets.

The ECal coordinates are projected from the target origin to each pulse's CDet
z coordinate using `earm.ecaldist`. The Run 5710 layer-x convention remains
`x_corr = 1.07*x - 0.03 m`; the y residual retains the macro's `0.10 m` offset.
These values are explicit analyzer database keys rather than constants embedded
in the reconstruction code.

In the same 2,000-event replay, the 5,589 complete pulses produced 4,480 broad
quality passes, 5,540 ECal-eligibility passes, 1,085 spatial passes, and 960
pulses passing all three decisions. An independent calculation found zero
decision mismatches and zero numerical differences in the aligned coordinate,
projections, and residuals at stored precision.

## Implemented fourth slice: cross-layer pairing

`SBSCDet` now forms Layer-1/Layer-2 candidates from calibrated pulses passing
the broad-quality, ECal-eligibility, and projected-spatial diagnostics. Candidate
formation and ranking use only CDet-to-CDet quantities:

```text
dt = t_L2 - t_L1
dx = x_L2 - x_L1
dy = y_L2 - y_L1
score = ((dt - dt_center)/time_scale)^2 + (dx/x_scale)^2
```

For Run 5710 the database supplies `dt_center = 0 ns`, `|dt| <= 15 ns`,
`|dx| <= 0.15 m`, `|dy| <= 0.08 m`, `time_scale = 1 ns`, and
`x_scale = 0.01 m`. Projected ECal x does not enter the score. Candidates are
sorted by score and accepted greedily with one-to-one pulse use; the database
currently permits multiple non-conflicting pairs per event, matching the macro
configuration.

The `pair.*` output preserves the two source pulse-array indices and pixel IDs,
corrected member and mean times, dt/dx/dy, score, and ECal-minus-pair-mean time.
The pulse collection remains unchanged.

In the 2,000-event Run 5710 replay, 178 events produced 184 accepted pairs. An
independent reconstruction reproduced every pair count, source-pulse identity,
dt, dx, dy, and score exactly. No pulse was reused in any accepted pair.

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
- The analysis macro forms its initial hit collection using broad LE and ToT
  limits, hit multiplicity, layer occupancy, bad-channel suppression, ECal
  cluster validity, ECal time/energy eligibility, and projected ECal x/y
  consistency. It then applies pixel offsets, the ECal-time correction,
  layer-dependent ToT time walk, and the run-dependent timing shift.
- The optional physics-motivated y-propagation correction (`n = 1.59` for Run
  5710) is currently evaluated only after pairs have been selected. It produces
  parallel timing diagnostics and does not change the accepted sample or any
  stored calibration constant.
- The event display's “best hits” mode adds a window on ECal time minus corrected
  CDet hit time. Its peak mean and width may be supplied or fitted from the
  display sample. That sample is conditioned on accepted pairs and a selected
  bar, so this is a display filter, not yet an authoritative replay-selection
  definition.
- `SBSGEPEArm::CoarseReconstruct()` already accesses ECal following the base
  coarse reconstruction, providing a potential explicit integration point.

The replay association must examine all valid decoded LE/TE pulse pairs, not
only the generic one-pulse-per-channel good-hit collection. Preserve channel
and pulse indices throughout processing. This does not require changing the
generic good-hit behavior for other detectors.

### Current macro pairing and selection order

The current macro has two distinct spatial decisions that must not be conflated:

1. During initial hit preselection, each CDet hit must be compatible with the
   ECal projection in x and y, using the configured projection tolerances and
   the layer-dependent CDet x alignment.
2. After corrected hit times are available, Layer-1/Layer-2 pair candidates are
   formed using only CDet-to-CDet quantities: `t_L2 - t_L1`, `x_L2 - x_L1`, and
   a same-side y consistency requirement. The ranking score contains only the
   layer-to-layer time and x differences. **ECal x is not used to rank candidate
   pairs.**

Candidates are sorted by that score and accepted greedily as one-to-one pairs,
so neither hit can be reused. The `display.allow_multiple_pairs` setting
controls whether more than one non-conflicting pair may survive in an event.
Final pair diagnostics then apply configured layer-time, layer-x, and
ECal-minus-pair-time windows. The ECal projection into the physical footprint of
both hit half-bars defines an additional diagnostic population; it is not part
of the pair-ranking score.

This ordering is the current parity reference, not automatically the optimal
replay design. In particular, the first replay implementation should expose
the broad-quality, ECal-spatial, cross-layer-pairing, and ECal-timing decisions
as separate flags so their efficiencies can be measured independently.

### Run 5711 trajectory-time candidate study

The subsequent Run 5711 study uses the ECal-to-target trajectory to compare the
expected and measured Layer-1/Layer-2 x displacement after replay-level pairs
have been formed. It then combines that trajectory residual with the
ECal-minus-pair-mean timing residual in an adjustable normalized ellipse. This
is a downstream candidate-ranking refinement; it does not change the current
replay-level CDet-only greedy pairing score.

The adopted macro working point has center `(0 m, -26 ns)`, scales
`(0.020 m, 5 ns)`, and normalized radius 2. In the 100,000-event Run 5711
sample it retains 10,361 pairs in 6,215 events. Energy- and timing-window scans
support retaining the production selections `3.0 < E_ECal < 4.5 GeV` and
`-10 < t_ECal < 10 ns` for efficiency-oriented ROI candidate generation.
Definitions, denominators, plots, and reproduction commands are recorded in
[`CDet_RUN5711_PAIR_SELECTION_STUDY.md`](CDet_RUN5711_PAIR_SELECTION_STUDY.md).

## Event inputs and database responsibilities

The database stores calibration constants, geometry and selection policy.

### LH2 run-dependent timing policy

For LH2 production runs, use exactly `p1 = 0`. The Run 5711 and Run 6077
studies both found slopes close to zero, and a zero ECal-time coefficient is
the physically motivated model for these coincidence data. The only
run-dependent calibrated timing constant is therefore `shift_ns`; the ECal
timing window may also vary by run as a selection parameter. Detector-wide
pixel offsets, ECal `p0`/`delta`, time-walk constants, propagation correction,
geometry, pulse-quality cuts, and pairing policy remain common.
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
| Optional y-propagation speed or refractive index and enable flag | Per calibration/validity period | Reproduce the physics-motivated event-level propagation correction without confusing it with empirical position-fit constants |
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

optional timing-only diagnostic:
y_projected[i] = y_ecal * z_cdet[i] / z_ecal
t_ycorr[i] = t_corr[i]
             - signed_slope[side] * (y_projected[i] - y_halfbar_center[i])
signed_slope[left]  = -1 / (c/n)
signed_slope[right] = +1 / (c/n)
```

`t_base` includes the agreed TDC conversion and reference treatment. Document
all units, signs, enable flags and the order of any broad cuts relative to
corrections. Apply common timing shifts consistently to LE and TE so ToT is
not unintentionally changed. Handle invalid/nonpositive ToT explicitly.

The active Run 5710 diagnostic uses `n = 1.59`, corresponding to
`c/n = 18.8549 cm/ns` and `|dt/dy| = 5.30367 ns/m`. Pair selection uses
`t_corr`, and `t_ycorr` is calculated only afterward on the identical accepted
sample. Do not promote `t_ycorr` into replay selection until its intended role,
validity range, and behavior for missing ECal y are explicitly approved.

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
5. Apply the fixed ECal-time, time-walk, and run-shift corrections to every
   eligible pulse without fitting from the replayed event sample.
6. Evaluate ECal spatial association and retain its result as an explicit flag.
7. Apply a fixed ECal-minus-CDet timing association and retain its result as a
   separate flag.
8. Optionally evaluate cross-layer consistency without prematurely removing
   isolated valid candidates. If macro-like pairs are produced, rank them only
   with CDet layer-to-layer quantities, never ECal x.
9. Optionally calculate the approved y-propagation-corrected time as a parallel
   quantity; do not silently substitute it for the primary corrected time.
10. Write candidates, residuals, selection flags and event association status.

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
- optional y-propagation-corrected time and the projected y used to calculate it;
- associated ECal cluster index;
- separate broad-quality, ECal-spatial, ECal-timing, and calibration-valid flags;
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
   macro for the same pulses before adding rejection. Keep the optional
   y-propagation time as a separately named output.
3. **Add ECal association.** Use all valid pulses, apply the established spatial
   cuts and fixed residual timing window, and export candidates and diagnostics.
   Preserve the component decisions as flags rather than immediately collapsing
   them into one opaque good-hit result.
4. **Validate on LH2.** Measure retained multiplicities, signal retention and
   accidental contamination. Resolve disagreements with the macro explicitly,
   including expected differences from retaining multiple pulses per channel.
5. **Evaluate grouping.** Add cross-layer consistency or ranking only when it
   demonstrably improves selection without unacceptable signal loss. For
   macro-parity studies, reproduce its greedy one-to-one matching and confirm
   explicitly that ECal x does not enter the pair-ranking score.

Do not transfer display-only settings such as selected bar, plot ranges or
histogram binning into replay physics selection. Event eligibility and optional
energy cuts must be explicit rather than inherited accidentally from a display
sample.

## Validation and acceptance criteria

- **Calibration parity:** corrected times and residuals agree with the macro for
  identical pulses and constants within a defined numerical tolerance.
- **Selection-stage parity:** broad-quality, ECal-spatial, cross-layer, and
  ECal-timing decisions can be compared independently with the macro rather than
  only through one final pass/fail bit.
- **Pairing parity:** when macro-like pairs are requested, their identities and
  ordering agree for the same preselected hits, including greedy one-to-one hit
  use and the absence of ECal x from the ranking score.
- **Y-correction parity:** when enabled, the projected y, side sign, and
  propagation-corrected time agree with the macro while the uncorrected
  candidate identity remains unchanged.
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
- [Event-display and pairing guide](CDet_EVENT_DISPLAY.md) — current pair
  construction, display-only best-hit behavior and configuration interface.
- [Y-position and propagation-correction study](CDet_Y_POSITION_RECONSTRUCTION.md)
  — empirical position response and the separate `n = 1.59` timing correction.
- [Master calibration](CDet_calibration_dt.dat),
  [Run 6077 timing constants](CDet_run6077.dat), and
  [Run 6077 analysis/display configuration](CDet_run6077_projection.conf).
