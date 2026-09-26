# CDet candidate integration into the SBS-offline ROI algorithm

Status: implementation in progress. This document records the interface and
data-model work needed to make the calibrated CDet candidates available to the
event-level region-of-interest (ROI) algorithm in SBS-offline.

## Objective

The calibration and diagnostic macros read reconstructed quantities from ROOT
trees. The eventual ROI algorithm instead runs inside the analyzer and must
consume the same CDet information directly from `SBSCDet`, before an output
tree exists. Detector-local reconstruction should remain owned by `SBSCDet`;
the global ROI module should combine the resulting CDet hypotheses with ECal,
HCal, and tracking information.

The implementation must preserve the existing `hit.*`, `pulse.*`, and `pair.*`
ROOT output and must not silently change the current pulse or pair selection.

## Current state

| Capability | State at start | Current state |
|---|---|---|
| Build complete LE/TE/ToT pulse candidates | Implemented in `SBSCDet` | Implemented |
| Apply pixel, time-walk, run-shift, p1, and y corrections | Implemented | Implemented |
| Apply broad pulse-quality and ECal-projection decisions | Implemented | Implemented |
| Build and rank Layer-1/Layer-2 pairs | Implemented | Implemented |
| Store final greedy-selected, one-to-one pairs | Implemented as `pair.*` | Implemented |
| Public read-only C++ pulse/pair interface | Not implemented | **Implemented in Step 1** |
| Preserve all valid alternate pair hypotheses | Not implemented | **Implemented in Step 2** |
| Produce accepted exclusive single-layer candidates | Macro study only | **Implemented in Step 3** |
| Consume CDet candidates in `SBSGEPRegionOfInterestModule` | Not implemented | Not implemented |

The ROOT branches are suitable for post-replay analysis, but they are not an
appropriate in-process interface between analyzer modules. Before Step 1, the
underlying pulse and pair vectors were protected members of `SBSCDet`.

## Responsibility boundary

`SBSCDet` is responsible for:

- complete-pulse construction and stable pulse identity;
- timing calibration and detector-coordinate corrections;
- detector-local quality and projection decisions;
- Layer-1/Layer-2 candidate construction and detector-local scores;
- exposing event-local candidates through a read-only interface.

The global ROI module is responsible for:

- combining CDet candidates with ECal, HCal, and tracking information;
- resolving ambiguous detector-level hypotheses using global information;
- deciding which hypotheses belong in the event ROI.

The ROI module should not reproduce CDet calibration or reach into protected
`SBSCDet` storage.

## Implementation sequence

1. Add a stable, public, read-only C++ interface for the existing pulse
   candidates and greedy-selected pairs.
2. Preserve every pair hypothesis that passes the detector-local hard gates
   and ECal ellipse before greedy pulse-sharing resolution. Keep the existing
   `pair.*` collection unchanged for backward compatibility. **Implemented.**
3. Implement explicit exclusive Layer-1-only and Layer-2-only ROI candidates,
   using the selection validated by the hydrogen studies. **Implemented.**
4. Add the remaining database-backed policy needed by the new candidate
   collections and expose their classifications and scores. **Implemented.**
4a. Migrate the principal diagnostic macro to prefer the analyzer candidate
   collections while retaining its independent reconstruction for old files
   and regression tests. **Implemented and validated on a limited sample.**
5. Connect `SBSGEPRegionOfInterestModule` to `SBSCDet` and consume the public
   candidate interface after ECal-dependent CDet calibration has run.
6. Validate candidate identity, multiplicity, backward compatibility, and ROI
   behavior with the established Run 5710 and hydrogen samples.

## Step 1: public read-only pulse and pair interface

Status: implemented.

`SBSCDet.h` now defines two event-local snapshot types:

- `SBSCDet::PulseCandidate`, containing the complete pulse identity,
  coordinates, raw and calibrated timing values, ECal residuals, and the three
  existing selection flags;
- `SBSCDet::LayerPair`, containing the two source-pulse indices and pixel IDs,
  member and mean times, layer residuals, CDet and ECal scores, and y topology.

The corresponding accessors are:

```cpp
Int_t GetNumPulseCandidates() const;
Bool_t GetPulseCandidate(Int_t index, PulseCandidate& pulse) const;

Int_t GetNumLayerPairs() const;
Bool_t GetLayerPair(Int_t index, LayerPair& pair) const;
```

Each accessor copies one immutable snapshot and returns `false` for an invalid
index. Callers therefore cannot mutate `SBSCDet` storage or retain a reference
that becomes invalid when the next event is cleared. Pulse indices stored in a
`LayerPair` refer to the same event's pulse-candidate collection.

Step 1 intentionally exposes only the already existing collections. In
particular, `GetLayerPair()` returns the current greedy-selected pairs. It does
not yet provide alternate hypotheses or single-layer candidates; those remain
separate, explicitly testable implementation steps.

No ROOT branch names, calibration calculations, selection cuts, greedy-pairing
behavior, or default replay configuration are changed by Step 1.

### Step 1 verification

The SBS-offline ROOT dictionary and complete `libsbs` shared-library target
were rebuilt successfully after adding the nested snapshot types and accessors:

```text
cmake --build ../build --target sbs -- -j4
[100%] Built target sbs
```

Both SBS-offline and SBS-replay pass `git diff --check`. A replay-data parity
test is not required for this slice because no reconstruction path or output
definition changed; event-level parity will be checked when the candidate
collections themselves are extended in Steps 2 and 3.

## Step 2: preserve valid pair hypotheses before greedy assignment

Status: implemented.

`SBSCDet::BuildLayerPairs()` already constructs every Layer-1/Layer-2
combination from pulses passing the pulse-quality, ECal-eligibility, and
per-pulse spatial decisions. Each combination must then pass the configured
hard limits on `dt`, `dx`, and y topology. When ECal-informed ranking is
enabled, it must also lie inside the configured trajectory-time ellipse.

Previously, these valid combinations existed only in a temporary local vector.
After sorting, the greedy assignment retained a subset in `pair.*` and
discarded every alternative sharing either pulse with a better-ranked pair.

Step 2 preserves the complete sorted vector as a new parallel collection:

```text
earm.cdet.pair_candidate.*
```

It contains the same fields as `pair.*`: source pulse indices and pixel IDs,
corrected member and mean times, `dt`, `dx`, `dy`, the CDet-only score, ECal
timing residual, trajectory residual, ECal score, and y topology. Candidate
indices are ranks in the sorted hypothesis list. A source pulse may appear in
multiple candidates by design.

The in-process read-only interface is:

```cpp
Int_t GetNumLayerPairCandidates() const;
Bool_t GetLayerPairCandidate(Int_t index, LayerPair& pair) const;
```

The accessor reuses the Step 1 `LayerPair` snapshot because a pre-greedy
hypothesis and a selected pair carry the same physical quantities. Their
collection membership gives them different semantics:

- `GetLayerPairCandidate()` / `pair_candidate.*`: every ranked hypothesis
  passing detector-local gates, including hypotheses that share pulses;
- `GetLayerPair()` / `pair.*`: the existing one-to-one greedy winners.

The new collection is filled only after candidate scoring and sorting, and
immediately before the pre-existing greedy loop. That loop, its ordering, its
one-pulse-use rule, and `pairing.allow_multiple` behavior are unchanged.

### Step 2 verification

The complete ROOT dictionary and `libsbs` target rebuild successfully:

```text
cmake --build ../build --target sbs -- -j4
[100%] Built target sbs
```

The diff confirms that the existing candidate gates, scores, sorting
comparison, and greedy assignment were not modified. The only insertion in
that path copies the already-ranked temporary candidates into the new parallel
collection before greedy assignment. Both repositories pass
`git diff --check`.

A fresh 100-event Run 5710 replay was then run with the rebuilt library and a
freshly compiled replay macro using the matching source header. It completed
normally and wrote all `pair_candidate.*` branches. The sample contained:

| Check | Result |
|---|---:|
| Physics events | 100 |
| Events with at least one valid pair candidate | 9 |
| Events with more candidates than greedy-selected pairs | 2 |
| Events with `candidate count < selected-pair count` | 0 |
| Events where candidate rank 0 did not match selected pair 0 | 0 |

The last two checks confirm the expected subset and ordering relationships in
this validation sample. Because this was deliberately a short interface test,
larger-sample multiplicity and physics-distribution comparisons remain part of
the final validation step.

## Step 3: exclusive single-layer candidates

Status: implemented.

The hydrogen macro study defines this recovery sample at the event level. A
pulse is first required to pass all four existing analyzer decisions:

```text
calibration valid
AND broad pulse quality
AND ECal eligibility
AND projected spatial compatibility.
```

The event is eligible only when these fully selected pulses populate exactly
one CDet layer. Events with pulses in both layers are excluded even if no pair
survives pairing. This makes the collection exclusive from both the two-layer
candidate population and unresolved two-layer ambiguities.

For every eligible pulse, `SBSCDet` computes

```text
delta_x_single = x_CDet,corr - x_ECal,projected
delta_t_single = t_ECal - t_CDet,corr

R_single^2 = ((delta_x_single - x0) / sx)^2
             + ((delta_t_single - t0) / st)^2.
```

The validated hydrogen working point is `x0 = 0 m`, `sx = 0.040 m`,
`t0 = -26 ns`, `st = 5 ns`, and `R_single <= 2`. Every pulse inside the
ellipse is retained; no detector-local best-candidate choice is imposed.

The new output collection is:

```text
earm.cdet.single_candidate.*
```

It stores the collection index, source pulse-array index, pixel ID, layer,
corrected leading-edge time, ECal timing residual, projected-x residual, and
normalized ellipse score. The corresponding read-only C++ interface is:

```cpp
Int_t GetNumSingleLayerCandidates() const;
Bool_t GetSingleLayerCandidate(Int_t index,
                              SingleLayerCandidate& candidate) const;
```

The source pulse index connects each entry to the complete Step 1
`PulseCandidate` information without duplicating every pulse field.

### Configuration behavior

`SBSCDet` recognizes the following optional database keys:

```text
single_layer.enable
single_layer.residual_center
single_layer.timing_center
single_layer.residual_scale
single_layer.timing_scale
single_layer.radius
```

For backward compatibility with older database files, a missing
`single_layer.enable` inherits `pairing.ecal_rank_enable`, and missing numerical
keys use the validated hydrogen working point above. Step 4 makes all of these
values explicit in the authoritative database. Invalid enabled scales or
radius cause database initialization to fail rather than silently accepting an
ill-defined cut.

### Step 3 verification

The ROOT dictionary and complete `libsbs` target rebuild successfully, and
both repositories pass `git diff --check`.

A normal 100-event Run 5710 replay using the unchanged authoritative database
completed successfully and wrote the new branches. Since Run 5710 is a
cross-target block with ECal-informed pairing disabled, the inherited default
correctly left single-layer recovery disabled.

For a direct functional test, a temporary database copy explicitly enabled the
validated single-layer policy for a 2,000-event Run 5710 replay. This changed no
repository database file. The result was:

| Check | Result |
|---|---:|
| Physics events | 2,000 |
| Events with accepted single-layer candidates | 33 |
| Accepted single-layer candidates | 33 |
| Events containing both a single-layer candidate and selected pair | 0 |
| Single-candidate events containing candidates from mixed layers | 0 |
| Candidates outside the configured ellipse | 0 |
| Candidates with invalid source-pulse indices | 0 |

An independent tree-level implementation of the macro definition then
reconstructed the expected source pulse indices and scores from `pulse.*`.
Across all 2,000 events it found zero count mismatches, zero identity
mismatches, and zero score mismatches.

## Step 4: explicit policy and candidate classification

Status: implemented.

### Authoritative database policy

`DB/db_earm.cdet.dat` now defines the validated single-layer ellipse once in
the common calibration block:

```text
earm.cdet.single_layer.residual_center = 0.0
earm.cdet.single_layer.timing_center = -26.0
earm.cdet.single_layer.residual_scale = 0.040
earm.cdet.single_layer.timing_scale = 5.0
earm.cdet.single_layer.radius = 2.0
```

The enable switch is explicit in every target-policy validity block. It is
kept synchronized with the established ECal-informed hydrogen candidate mode:

- `single_layer.enable = 1` for non-cross-target/hydrogen intervals;
- `single_layer.enable = 0` for cross-target intervals.

All 47 `pairing.ecal_rank_enable` policy declarations have a matching
`single_layer.enable` declaration with the same state. This removes dependence
on the compatibility fallback for the current five-pass database while keeping
older database files readable.

### Pair-hypothesis classification

Every `pair_candidate.*` hypothesis now includes:

```text
greedy_selected
selected_pair_index
```

`greedy_selected` is one when the hypothesis also appears in the existing
one-to-one `pair.*` collection. `selected_pair_index` gives that collection's
index, or `-1` for an alternate hypothesis rejected only by greedy pulse-use
resolution. The `LayerPair` C++ snapshot exposes the same two fields. A
snapshot read through `GetLayerPair()` is always marked selected; a snapshot
read through `GetLayerPairCandidate()` reports its actual classification.

### Event-level classification

The scalar `earm.cdet.roi.status`, also available through
`GetROICandidateStatus()`, provides a compact mutually exclusive event result:

| Value | Meaning |
|---:|---|
| 0 | No accepted ROI candidate |
| 1 | At least one two-layer pair hypothesis |
| 2 | Accepted exclusive Layer-1-only candidate(s) |
| 3 | Accepted exclusive Layer-2-only candidate(s) |

Candidate-level scores remain available as `pair_candidate.score`,
`pair_candidate.ecal_score`, and `single_candidate.score`. The source pulse
indices retain the full connection to `pulse.*`.

### Step 4 verification

The ROOT dictionary and full `libsbs` target rebuilt successfully. A fresh
2,000-event Run 5710 replay with the authoritative database completed normally.
The cross-target policy produced no single-layer candidates, as configured.
For its pair hypotheses:

- every event had exactly as many `greedy_selected` flags as `pair.*` entries;
- all classification flags were either zero or one;
- every selected candidate had a nonnegative selected-pair index;
- every alternate candidate had selected-pair index `-1`;
- every event with pair hypotheses had ROI status 1.

The explicit-enabled 2,000-event functional sample from Step 3 was replayed
with the Step 4 classification code. Its 33 single-candidate events divided
into 25 Layer-1-only status-2 events and 8 Layer-2-only status-3 events, with
zero status/layer mismatches and zero pair-status mismatches.

An automated database audit found 47 paired target-policy declarations and no
state mismatch between `pairing.ecal_rank_enable` and
`single_layer.enable`. Both repositories pass `git diff --check`.

## Step 4.5: analyzer-native plotting and regression path

Status: implemented and validated on a limited Run 5710 sample.

`Plot_CDet_GoodPulseCandidates_AllTDC.C` now chooses its candidate source from
the branches actually present in the input tree:

- new trees use `pair_candidate.*` and `single_candidate.*` directly;
- older trees without those branches retain the existing `pair.*` and
  macro-level single-layer reconstruction;
- each run prints the selected source explicitly, so the fallback cannot be
  mistaken for analyzer-native operation.

The macro still reconstructs the single-layer selection independently from
`pulse.*` on new trees. It compares the analyzer and reconstructed source-pulse
indices event by event and reports count and identity mismatches. It also
checks each `pair_candidate.greedy_selected` classification against the
corresponding `pair.*` identity and checks `roi.status` against the candidate
collections. Thus the duplicated logic has become a regression oracle rather
than the nominal source of ROI candidates.

This first comparison exposed a bookkeeping defect before ROI integration:
the geometry element layer field was zero for pixels in both physical layers.
Pair construction was unaffected because it already uses the authoritative
pixel ranges, but `pulse.layer` and consequently `single_candidate.layer` were
incorrect for physical Layer 2. `SBSCDet` now defines the exported layer as
`pixel_id / 1344`, with 0 denoting Layer 1 and 1 denoting Layer 2.

After rebuilding `libsbs`, a fresh 2,000-event replay with single-layer
selection explicitly enabled gave:

| Regression check | Result |
|---|---:|
| ECal-energy-selected pair hypotheses | 49 |
| Pair classification count-mismatch events | 0 |
| Pair classification identity mismatches | 0 |
| Analyzer single-candidate events | 13 |
| Analyzer/reconstruction count mismatches | 0 |
| Analyzer/reconstruction identity mismatches | 0 |
| ROI-status consistency mismatches | 0 |
| Recovered Layer-1-only events/pulses | 12 / 12 |
| Recovered Layer-2-only events/pulses | 1 / 1 |

The macro also compiled successfully and automatically selected its legacy
fallback on the existing pre-candidate Run 5710 files. No full production
replay was performed for this checkpoint; broader Run 5710 and hydrogen
comparisons remain part of the final validation step.
