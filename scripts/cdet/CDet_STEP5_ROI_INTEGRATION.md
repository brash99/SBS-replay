# Step 5: integrating CDet candidates with the GEP ROI module

Status: Step 5A and the non-operative Step 5B diagnostic slice are
implemented. CDet-derived hypotheses are exported for study, but they do not
yet alter any GEM constraint.

This document tracks Step 5 of
`CDet_SBS_OFFLINE_ROI_IMPLEMENTATION.md`: connect the validated, public
`SBSCDet` candidate interface to `SBSGEPRegionOfInterestModule` without
duplicating CDet calibration or prematurely discarding alternate CDet
hypotheses.

## Scope of this first checkpoint

Before changing the ROI implementation, this checkpoint establishes:

- when the ROI module runs in the analyzer event sequence;
- which detector information it currently consumes;
- what it currently means by a region of interest;
- which GEM constraints it creates;
- when calibrated CDet candidates become available;
- the safe software boundary for adding CDet;
- the physics decisions that must be explicit before implementation.

The source files reviewed for this description are:

- `SBS-offline/SBSGEPRegionOfInterestModule.h` and `.cxx`;
- `SBS-offline/SBSGEPEArm.cxx`;
- `SBS-offline/SBSEArm.cxx`;
- `SBS-offline/SBSGEMSpectrometerTracker.cxx`;
- `SBS-offline/SBSCDet.h` and `.cxx`;
- `SBS-replay/replay/replay_gep.C`;
- `SBS-replay/DB/db_FTROI.dat`.

## What the existing ROI module is

`SBSGEPRegionOfInterestModule` is a Podd `InterStageModule`. It is not a
detector and it does not itself perform GEM tracking. Its present job is to
run between analyzer stages, use already reconstructed calorimeter
information to predict a family of elastic proton trajectories, and install
front/back constraint points in the two GEM trackers before constrained
tracking begins.

`replay_gep.C` constructs it as:

```cpp
new SBSGEPRegionOfInterestModule(
    "FTROI", "GEP region of interest calculation",
    THaAnalyzer::kCoarseRecon);
```

It is registered only when `dogems != 0`. The `kCoarseRecon` stage means its
`Process()` method runs after apparatus coarse reconstruction and before the
later tracking/reconstruction work that consumes the installed constraints.

In this code, "ROI" currently means **a set of allowed constraint rays for
GEM tracking**, not an event display region, a rectangular detector cut, or a
final physics-event decision.

## Event sequence relevant to CDet

The relevant event flow is:

```text
detector Decode
    |
apparatus CoarseReconstruct
    |
    +-- ECal CoarseProcess and cluster reconstruction
    |
    +-- CDet CoarseProcess
    |       build complete pulse triplets
    |
    +-- SBSGEPEArm::CoarseReconstruct
            call CDet::ApplyECalTimingCalibration(...)
            apply timing and position corrections
            build pair_candidate.*, pair.*, and single_candidate.*
            create the E-arm pseudo-track from the best ECal cluster
    |
FTROI inter-stage Process
    |       read already reconstructed detector information
    |       clear and install GEM constraint points
    |
fine/constrained GEM tracking and later reconstruction
```

This ordering is essential. CDet candidate construction cannot be moved into
`SBSGEPRegionOfInterestModule`, because the calibration and candidate
decisions belong to `SBSCDet` and are already complete before the ROI module
runs. Conversely, the ROI module is the first natural global location that
can combine those candidates with other detectors.

## Detector discovery and prerequisites

On every event, `SBSGEPRegionOfInterestModule::Process()` scans `gHaApps` and
finds objects using names loaded from `db_FTROI.dat` (or constructor defaults):

| Role | Default apparatus/detector | Required class |
|---|---|---|
| Electron arm | `earm` | `SBSGEPEArm` |
| Electron calorimeter | `earm.ecal` | `SBSECal` |
| Hadron arm | `sbs` | `SBSEArm` |
| Front GEM tracker | `sbs.gemFT` | `SBSGEMSpectrometerTracker` |
| Polarimeter/back GEM tracker | `sbs.gemFPP` | `SBSGEMPolarimeterTracker` |
| Hadron calorimeter | `sbs.hcal` | `SBSHCal` |

The current module returns with invalid data if any of these objects is
missing, or if the hadron arm is not in GEP tracking mode. It clears the
existing front- and back-tracker constraints before imposing event-specific
ones.

It then requires:

```text
at least one HCal cluster
AND at least one E-arm track.
```

The E-arm "track" used here is presently the pseudo-track constructed from
the best ECal cluster, not an independently fitted charged-particle track.
Therefore the current ROI seed is calorimeter driven.

## How the current ROI is calculated

### 1. Establish the measured calorimeter points

The module takes the best HCal cluster and expresses its position in the SBS
focal-plane coordinate system. It takes the first E-arm pseudo-track and uses
its ECal position and energy to build the global ECal cluster position.

The ECal position is also passed to the front GEM tracker through
`SetECALpos()`, where it may be used by the tracker's separate elastic
constraint logic.

### 2. Calculate a central elastic solution

For a nominal target vertex `z0targ`, the ECal direction defines the electron
scattering angles. With the run beam energy and elastic electron-proton
kinematics, the module calculates:

- expected scattered-electron energy;
- expected proton momentum;
- expected proton polar and azimuthal angles.

The predicted proton ray is transformed into SBS transport coordinates and
transported through the forward optics matrix to obtain a central focal-plane
position and direction:

```text
xfp0, yfp0, xpfp0, ypfp0.
```

These values are exported as diagnostic variables; they are not a reconstructed
GEM track.

### 3. Scan the extended target

`db_FTROI.dat` currently defines 24 target-z bins. For the later five-pass
period the nominal range is `-0.28 < z < 0.08 m`. For every assumed vertex,
the module repeats the elastic-kinematics and forward-optics calculation.

This produces a family of allowed proton rays rather than one point-target
prediction. The target scan is the present source of multiple front/back
constraint-point pairs.

### 4. Install GEM constraints

For the front tracker (`gemFT`):

- the front constraint is placed at its local `z = 0` plane;
- the back constraint is placed 5 cm beyond its last GEM layer;
- each point follows the predicted proton focal-plane ray, plus configured
  constraint offsets.

For the polarimeter/back tracker (`gemFPP`):

- the front constraint is placed at the analyzer midpoint;
- the back constraint is anchored to the measured HCal cluster position.

When a tracker is in multi-track mode, all vertex-bin hypotheses are added.
Otherwise, only the central hypothesis is installed. The tracker later uses
these points and its independently configured constraint widths to restrict
track finding.

## What the ROI module currently exports

The `FTROI.*` variables are diagnostics describing the calorimeter seed and
central elastic prediction:

- global ECal x, y, and z;
- ECal energy;
- predicted electron and proton angles;
- predicted scattered-electron energy and proton momentum;
- central predicted focal-plane position and slopes;
- the inherited `InterStageModule` data-valid state.

It presently exports no CDet identity, candidate count, candidate score,
candidate-to-constraint association, or ROI acceptance classification.

## CDet information available at the Step 5 boundary

By the time `FTROI.Process()` begins, `SBSCDet` provides four distinct levels
of information through public read-only accessors:

| Interface | Meaning |
|---|---|
| `GetPulseCandidate()` | Complete calibrated pulse with detector coordinates, ECal projection/residuals, and quality flags |
| `GetLayerPairCandidate()` | Every ranked two-layer hypothesis passing detector-local gates and the configured ECal trajectory-time ellipse; hypotheses may share pulses |
| `GetLayerPair()` | Final detector-local one-to-one greedy pairs retained for backward compatibility |
| `GetSingleLayerCandidate()` | Accepted candidates in exclusive one-layer events |

`GetROICandidateStatus()` summarizes whether the event has pair hypotheses,
Layer-1-only candidates, Layer-2-only candidates, or no accepted CDet
candidate.

For global ROI work, `GetLayerPairCandidate()` is the important two-layer
interface. The final `GetLayerPair()` collection has already discarded
alternatives through a CDet-local one-to-one choice. The purpose of the global
ROI is precisely to allow ECal, CDet, GEM, and HCal information to resolve
such ambiguities. The selected-pair collection should remain available for
comparison and backward compatibility, but should not silently replace the
candidate collection at this boundary.

## Safe software integration point

The smallest structurally correct Step 5 change is inside
`SBSGEPRegionOfInterestModule`, not inside `SBSCDet::FineProcess()`:

1. find `earm.cdet` alongside `earm.ecal` during detector discovery;
2. verify that CDet timing calibration completed for the event (there is not
   yet a public `GetTimingStatus()` accessor, so Step 5A must either add a
   read-only status accessor or define an equally explicit validity test);
3. read immutable pulse, pair-candidate, and single-layer snapshots through
   the public `SBSCDet` accessors;
4. translate accepted CDet hypotheses into explicit global-ROI hypotheses;
5. preserve the association between every ROI hypothesis and its source CDet
   candidate(s);
6. install or publish the resulting constraints without changing the
   detector-local collections.

This preserves the ownership boundary established in Steps 1-4: `SBSCDet`
calibrates and classifies CDet information; the ROI module combines detector
information.

## Physics decisions established for Step 5

The current ROI predicts the **proton** trajectory from ECal elastic
kinematics and uses HCal to constrain the back tracker. CDet measures the
electron arm. The following decisions now define how CDet enters that logic.

### The existing fallback is unconditional

If an event has no accepted CDet candidate, it must retain the current
ECal/HCal-only ROI calculation. Missing CDet information is not grounds for
rejecting an otherwise processable event. This is both a physics requirement
and the principal backward-compatibility rule for Step 5.

### ECal and CDet share one rigid detector frame

CDet and ECal are physically housed in the same rigid detector frame, and the
CDet database positions have been characterized in the ECal coordinate
system. CDet points must therefore use the same detector-to-global rotation
and translation used for ECal. Step 5 must not introduce an independent CDet
alignment frame or a second set of arm axes.

Concretely, the ROI code should first represent all ECal and CDet measurements
as points in their shared electron-arm frame, and then apply the existing
electron-arm-to-global transformation consistently to every point. A useful
validation invariant is that projecting the transformed CDet points back into
the electron-arm frame reproduces the stored CDet coordinates.

### All available three-dimensional detector points contribute

For a two-layer CDet hypothesis, the electron-trajectory calculation should
use all three measured points:

```text
(x, y, z)_ECal
(x, y, z)_CDet Layer 1
(x, y, z)_CDet Layer 2
```

For an exclusive single-layer hypothesis, it should use the ECal point and
the one available CDet point:

```text
(x, y, z)_ECal
(x, y, z)_CDet Layer 1 or Layer 2.
```

The pulse indices stored in a `LayerPair` identify the two
`PulseCandidate` snapshots containing the required coordinates. A
`SingleLayerCandidate` similarly identifies its source pulse. The ROI module
should obtain coordinates through those source identities rather than
reconstructing pixel geometry independently.

The two-layer case is an overconstrained straight-line measurement in x and
can provide both a best electron ray and a fit-consistency diagnostic. The
single-layer case supplies only two detector points, so it defines an x-z ray
but has less redundancy. Its uncertainty should follow from the ECal and CDet
position resolutions rather than from an arbitrary extra cut.

### CDet y is a coarse half-bar constraint, not a precision measurement

At present, CDet supplies essentially no useful independent continuous-y
measurement. The stored y coordinate identifies the relevant half-bar, but it
must not be assigned the same precision as the CDet x coordinate.

For tracking, each CDet point therefore carries:

```text
y_CDet = nominal half-bar center
uncertainty interval = y_CDet +/- (half-bar length)/2.
```

This interval expresses the actual detector information: the particle passed
somewhere along that half-bar. The half-bar length should come from the
detector geometry applicable to the source pulse rather than from a new
independent hard-coded tracking constant.

If the initial line fitter requires a Gaussian sigma rather than an interval,
Step 5 will conservatively use the half-width itself,
`sigma_y,CDet = (half-bar length)/2`. It will not convert a uniform interval
to the smaller RMS value `length/sqrt(12)`, because the present purpose is to
avoid assigning unjustified precision to CDet y. This should give the CDet y
points practically zero weight relative to ECal y while keeping their coarse
geometric information available for diagnostics and gross-consistency tests.

Consequently, the initial electron-ray model is deliberately asymmetric:

- x-z fitting uses the meaningful ECal and CDet layer positions and their
  measured resolutions;
- y-z fitting is driven primarily by ECal and the target/vertex hypothesis,
  with CDet contributing only a very broad half-bar compatibility interval.

The fit output should report the CDet y pull or interval compatibility so this
assumption can be verified from data. A future independent CDet y
reconstruction can replace this uncertainty model without changing the
candidate identities or the Step 5 interface.

### Initial CDet x uncertainty from Run 5710

The initial effective CDet x uncertainty is taken from the standard deviation
of the ECal-trajectory residual in the Run 5710 candidate comparison (row 2,
column 3 of `CDet_Run5710_CandidateBranchComparison.pdf`). The underlying
histogram contains:

```text
entries = 365,256
mean    = -0.00536588 m
sigma   =  0.0179730 m = 1.79730 cm.
```

The independently reconstructed and analyzer-native histograms have identical
statistics. Step 5 will therefore begin with:

```text
sigma_x,CDet = 0.017973 m
```

for each CDet layer point in the electron-ray fit. This is intentionally a
conservative **effective** per-point uncertainty, not a claim that the
intrinsic position resolution of one layer is 1.80 cm.

The measured residual is

```text
(x_L2 - x_L1) - x_ECal * (z_L2 - z_L1) / z_ECal,
```

so its variance contains both CDet layer uncertainties, the smaller
ECal-projection contribution, and any residual alignment, selection, or
non-Gaussian effects. If the two CDet layers had equal independent resolution
and every other contribution were negligible, the corresponding per-layer
value would instead be `1.797/sqrt(2) = 1.27 cm`. Using the full measured
1.80 cm for each point avoids overstating CDet precision in the initial ROI
fit. Pull distributions from Step 5B will show whether this conservative value
should later be refined.

### Provisional ECal x and y uncertainties

The ECal group reports position uncertainties of 6 mm in both transverse
coordinates. These values have not been independently established by the
present CDet study, and there is some concern that they may be optimistic.
Nevertheless, they are the appropriate collaboration inputs to use unless and
until data provide evidence for different values. The initial electron-ray fit
will therefore use:

```text
sigma_x,ECal = 0.006 m
sigma_y,ECal = 0.006 m.
```

These must be configurable inputs rather than immutable hard-coded constants.
With `sigma_x,CDet = 0.017973 m`, ECal will receive substantially greater
weight than either CDet layer in the x-z fit. In y, the 6 mm ECal uncertainty
is much smaller than the CDet half-bar uncertainty, so ECal should dominate
the y-z information as intended.

Step 5B must export normalized residuals or pulls separately for ECal x, ECal
y, and each CDet measurement. Their widths should be examined versus run,
detector position, ECal energy, and candidate class. Systematically non-unit
pull widths would be evidence that the assumed measurement uncertainties—or
the geometrical model—need revision. Until that validation is complete, the
6 mm values are working assumptions and must not be presented as results of
the CDet analysis.

### CDet supplements the target-z scan

A CDet hypothesis should initially supplement, not replace, the current
extended-target scan. The scan remains important because the interaction
vertex is not known a priori.

There is a genuine coupling that must be studied. The present CDet
trajectory-time ellipse is formed relative to an ECal-to-nominal-target ray,
effectively using `(x,y,z)_target = (0,0,0)` during candidate selection. The
true event vertex is distributed through the extended target. The large
target-to-detector distance should make this a modest effect, but that must be
demonstrated rather than assumed.

Two closely related formulations should be compared before the
physics-bearing constraint change:

1. fit the detector-only electron ray from the ECal and available CDet
   points, then compare its target intersection with each target-z hypothesis;
2. for every target-z bin, include the assumed target point together with the
   ECal/CDet measurements in a resolution-weighted constrained line fit and
   retain the fit quality.

This comparison will show whether CDet meaningfully narrows or reweights the
existing vertex scan and whether the nominal-origin candidate ellipse biases
the accepted vertex distribution.

### Resulting global interpretation

For every accepted CDet hypothesis, the intended global calculation is now:

1. build an electron-ray hypothesis from ECal and all available CDet layer
   points in their common frame;
2. combine that ray with the existing target-z scan;
3. recalculate the elastic proton kinematics for the resulting electron/vertex
   hypothesis;
4. create the corresponding proton-arm GEM constraint family;
5. retain multiple CDet-derived families when CDet is ambiguous, allowing the
   later global tracking information to resolve them.

If no CDet hypothesis exists, execute the present ECal/HCal-only calculation
unchanged.

## Physics and representation questions still open

The decisions above settle the basic geometry and fallback policy. The
following still require explicit choices or measurement:

- validation of the provisional 6 mm ECal x and y uncertainties using fitted
  residuals and pulls;
- the exact statistic used to rank electron-ray/target-z hypotheses;
- how much the current nominal-origin pair ellipse changes acceptance across
  the target length;
- whether single-layer hypotheses require broader downstream GEM constraints
  after normal uncertainty propagation;
- how duplicate or nearly identical constraint families are merged or capped;
- which source indices, fit residuals, fit quality, vertex-bin identity, and
  CDet scores must be exported for auditability.

These are physics-policy decisions, not merely software plumbing. They should
be resolved in this document before modifying the tracker or other code owned
outside `SBSCDet` and `SBSGEPRegionOfInterestModule`.

## Proposed incremental implementation sequence

### Step 5A: read-only CDet discovery and diagnostics

- Add the configured CDet detector name to the ROI database interface.
- Locate `SBSCDet` and read its status and candidate counts after coarse
  reconstruction.
- Export diagnostic counts and source classifications only.
- Do not change any GEM constraint.

This proves stage ordering and object access with essentially zero physics
risk.

Implemented in `SBSGEPRegionOfInterestModule`:

- `earm.cdet` is discovered using the configurable `FTROI.cdet_name`;
- `SBSCDet` timing status and detector-local ROI status are read through
  public, read-only accessors;
- pulse, pre-greedy pair-candidate, and exclusive single-layer counts are
  exported under `FTROI.cdet.*`;
- absence of CDet is non-fatal and leaves the established ROI path unchanged.

### Step 5B: construct auditable ROI hypotheses without applying them

- Build internal ROI-hypothesis records from CDet pair candidates and
  exclusive single-layer candidates.
- Transform the ECal and source-pulse CDet points through their common rigid
  electron-arm frame and form the two- or three-point electron-ray fit.
- Associate each electron-ray hypothesis with the existing target-z scan.
- Store source candidate/pulse indices, vertex-bin identity, all geometry
  inputs, fit residuals, and fit quality.
- Export enough information to compare the predicted constraints against the
  existing ECal/HCal-only constraints in replayed events.
- Still do not change GEM tracking.

The initial non-operative Step 5B slice is also implemented. For every
pre-greedy pair candidate and every exclusive single-layer candidate it:

- retrieves the source pulse coordinates by their event-local identities;
- performs independent resolution-weighted straight-line fits in x-z and
  y-z, using the configurable uncertainties documented above;
- transforms the fitted direction from the shared ECal/CDet frame into the
  global Hall frame using the existing E-arm axes;
- exports source type, source index, pulse indices, source score, number of
  fitted points, intercepts, slopes, chi-squares, NDF, global angles, and the
  fitted intercept at the nominal target z;
- projects each fitted ray through every bin of the existing target-z scan
  and exports a flattened, source-indexed association table.

The source-type encoding is `1 = pair candidate`, `2 = Layer-1-only`, and
`3 = Layer-2-only`. A two-point single-layer fit has zero fit NDF by
construction; its value is in the propagated ray and subsequent global
comparison, not an internal chi-square test.

No Step 5B value is passed to either tracker. The calls that install front and
back GEM constraints still consume only the original ECal/HCal calculation.

### Step 5C: optional application to GEM constraints

- Add a database switch, defaulting off, that allows CDet-derived hypotheses
  to affect tracker constraints.
- Retain the established ECal/HCal-only path as an explicit fallback and
  control sample.
- Apply multiple constraint families rather than forcing an early unique CDet
  choice.

### Step 5D: validation

- Verify no behavior change with the new switch disabled.
- Replay the established Run 5710, Run 5711, and Run 6077 samples.
- Audit candidate identity through the ROI hypothesis and resulting GEM
  constraint.
- Compare tracking efficiency, multiplicity, accidental rate, and elastic
  residuals against the unchanged baseline.

## Initial invariants

The Step 5 implementation must preserve the following unless a later,
explicitly reviewed physics decision changes them:

- CDet calibration remains entirely within `SBSCDet`.
- `pulse.*`, `pair_candidate.*`, `pair.*`, and `single_candidate.*` retain
  their current definitions.
- Every global hypothesis remains traceable to event-local CDet source
  indices.
- Multiple valid pair hypotheses remain available to the global algorithm.
- ECal and CDet use one common rigid detector coordinate frame and one common
  transformation into global Hall coordinates.
- Two-layer hypotheses use ECal plus both CDet layer points; exclusive
  single-layer hypotheses use ECal plus the available CDet layer point.
- CDet y is represented by its nominal half-bar center with an uncertainty
  interval of plus or minus one-half of the half-bar length. If a Gaussian
  fitter is used initially, that half-width itself is the conservative sigma.
- The initial ECal x and y uncertainties are both 0.006 m, treated as
  configurable, externally supplied working assumptions and tested with pull
  distributions before any physics-bearing ROI constraint is enabled.
- CDet supplements rather than removes the established target-z scan during
  the initial implementation and validation.
- Cross-target runs, for which hydrogen-specific CDet candidate production is
  disabled, do not acquire an accidental hydrogen ROI policy.
- Existing ROI/GEM behavior is bit-for-bit unchanged when CDet integration is
  disabled.
- Missing CDet information cannot invalidate an event that the current ROI
  would otherwise process; the current ECal/HCal-only path is mandatory.

## Current conclusion

The software timing and interface have now been exercised through the first
Step 5 implementation slice: calibrated CDet candidates are read after coarse
reconstruction, converted into auditable electron-ray and target-z diagnostic
hypotheses, and exported without exposing mutable detector storage.

The next task is empirical validation of these new branches in replayed Run
5710, Run 5711, and Run 6077 events. That study should determine the pull
widths, target-z behavior, hypothesis multiplicity, and useful ranking or
deduplication rules. Only after those diagnostics are understood should the
database-disabled Step 5C path be allowed to modify GEM constraints.
