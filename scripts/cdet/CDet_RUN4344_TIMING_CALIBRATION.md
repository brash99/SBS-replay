# CDet Run 4344 Timing Calibration

## Purpose and present scope

This is the authoritative restart procedure for Run 4344. It uses the accepted
Run 5710 detector calibration in `CDet_calibration_dt.dat`, retains its
detector-wide ECal slope, and determines only the additive Run 4344 timing
origin.

The earlier iterative Run 4344 values are exploratory and must not be treated
as accepted constants. In particular, do not alternate
`Run_CDet_Calibrate_RunSpecificP1` and
`Run_CDet_Calibrate_RunTimingShift`. The current `p1` driver reconstructs its
sample after every slope change, so successive fits do not necessarily use the
same events, hits, pairs, or half-bars. A run-specific Run 4344 `p1` may be
introduced only by the frozen-sample qualification defined below. That
qualification is implemented by
`Run_CDet_Qualify_RunSpecificP1_FrozenSample.C`. The completed qualification
found that Run 4344 does **not** qualify for a run-specific `p1`: its apparent
residual slope failed the contiguous-fold validation gate. The accepted result
therefore retains the Run 5710 master `p1` and uses only the Run 4344
`shift_ns`.

Run every command from `scripts/cdet`. Use a fresh ROOT process for every
compiled macro command.

## Inputs that must remain fixed

- `CDet_calibration_dt.dat` must be the accepted Run 5710 master calibration.
- `CDet_run4344_projection.conf` is the authoritative Run 4344 selection.
- `analysis.run_number` must be `4344`.
- The ECal timing window must be explicitly reviewed and recorded in the
  configuration before calibration.
- The restart run file must not contain an `[ECalTiming] p1` key. Its omission
  makes the analysis inherit the master `p1` from `CDet_calibration_dt.dat`.

The similarly named `CDet_run4344_projections.conf` is not an input to this
procedure.

## 1. Preserve the exploratory result

Do not overwrite the existing exploratory artifacts. Create a new archive
directory and move the current run file into it. If the example directory name
already exists, choose a new explicit name rather than merging archives.

```bash
mkdir CDet_run4344_pre_restart_20260920
cp CDet_run4344_projection.conf CDet_run4344_pre_restart_20260920/
mv CDet_run4344.dat CDet_run4344_pre_restart_20260920/CDet_run4344_iterated.dat
```

Also copy any plots or terminal log that support the exploratory result into
that directory. Confirm that the active run file is absent:

```bash
test ! -e CDet_run4344.dat
```

Do not move or modify `CDet_calibration_dt.dat`.

## 2. Verify the master calibration

Inspect the master ECal constants:

```bash
sed -n '/\[ECalTiming\]/,/^\[/p' CDet_calibration_dt.dat
```

For the currently accepted Run 5710 master, this section is expected to contain:

```ini
[ECalTiming]
p0 14.807420
p1 0.810203
delta 31.680784
```

Stop if the file is absent or these are not the intended accepted constants.
Do not copy `p1` into `CDet_run4344.dat`; inheritance by omission is deliberate.

## 3. Review and approve the Run 4344 ECal window

Generate the advisory timing-window plots without modifying the configuration:

```bash
root -l -b -q 'Run_CDet_Select_ECalTimingWindow.C+(4344)'
```

Inspect `CDet_run4344_ECalTimingWindow_PROPOSAL.png`, its PDF, and its ROOT
file. Select the physical main-cluster bounds from `earm.ecal.adctime`; do not
derive the window from `earm.ecal.atimeblk` or block-level
`earm.ecal.a_time`.

The current Run 4344 configuration contains the candidate interval
`10--35 ns`. Preview that interval explicitly:

```bash
root -l -b -q 'Run_CDet_Select_ECalTimingWindow.C+(4344,false,10.0,35.0)'
```

If human review approves those bounds, record that approval atomically:

```bash
root -l -b -q 'Run_CDet_Select_ECalTimingWindow.C+(4344,true,10.0,35.0)'
```

If different bounds are approved, replace `10.0` and `35.0` in both commands
with the reviewed values.

Then verify that `CDet_run4344_projection.conf` contains the correct run and
explicit approved window:

```bash
rg '^analysis\.(run_number|ecal_time_min|ecal_time_max):' CDet_run4344_projection.conf
```

Do not continue until those three values are correct.

## 4. Fit only the Run 4344 timing origin

With `CDet_run4344.dat` absent, run the shift calibration once:

```bash
root -l -b -q 'Run_CDet_Calibrate_RunTimingShift.C+(4344,-1,30.0,5.0,"CDet_run4344_projection.conf")'
```

This command constructs the sample with the master `p1`, writes only
`[GlobalTiming] shift_ns`, reruns the analysis, and performs one deterministic
refinement internally if the first centroid closure misses tolerance. It must
report all of the following before the result is accepted:

- the baseline analysis and Gaussian-core fit succeed;
- the final centroid agrees with `30 ns` within
  `max(0.020 ns, 3*sigma_mu)`;
- the accepted-pair counts before and after the additive shift are identical;
- `gLastCalibrationSequenceSucceeded` is true; and
- the output file contains no run-specific `p1`.

Inspect the resulting file:

```bash
sed -n '1,120p' CDet_run4344.dat
```

Its calibration content must have this form:

```ini
# Run-specific timing calibration for cross-target run 4344.

[GlobalTiming]
shift_ns <fitted value>
```

If the command ends with `ERROR`, preserve the failed candidate for diagnosis,
remove it from the active path, and restart this step from an absent run file.
Do not repeatedly invoke the shift fitter on an unaccepted candidate.

## 5. Perform a fresh final analysis

Use a new ROOT process and the same authoritative configuration for loading and
plotting:

```cpp
root -l
.L PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C+
PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget(
    "CDet_run4344_projection.conf", 7, -1);
plotCDetLayersTimeComp("CDet_run4344_projection.conf", 0);
.q
```

Record at least:

- the loaded master `p1` and Run 4344 `shift_ns`;
- accepted event, hit, and pair counts;
- the Gaussian-core centroid and uncertainty;
- the combined and per-layer within-half-bar fixed-effects slopes;
- the fixed-effects entry and half-bar counts;
- the paired-time width;
- the ECal profile fit quality; and
- the generated timing diagnostic plots.

The pooled profile slope is diagnostic only. The within-half-bar fixed-effects
slope is the relevant estimator because it removes the separate half-bar
intercepts. A poor pooled-profile chi-square or sparse endpoint bins must not be
used as justification for changing `p1`.

## 6. Frozen-sample qualification for a possible Run 4344 `p1`

Do not use `Run_CDet_Calibrate_RunSpecificP1.C` for this decision. Start a
fresh ROOT process, compile the dedicated read-only qualifier, and invoke it
with the full argument list:

```cpp
root -l
.L Run_CDet_Qualify_RunSpecificP1_FrozenSample.C+
Run_CDet_Qualify_RunSpecificP1_FrozenSample(
    4344,
    -1,
    "CDet_run4344_projection.conf",
    5,
    3.0,
    0.5,
    2.0,
    0.5,
    1.05);
.q
```

At an interactive `root [n]` prompt, do **not** enter
`Run_CDet_Qualify_RunSpecificP1_FrozenSample.C+(...)`. The `.C+` suffix belongs
only to the preceding `.L` compilation command. The subsequent invocation is a
normal C++ function call whose name is
`Run_CDet_Qualify_RunSpecificP1_FrozenSample(...)`, as shown above.

The equivalent noninteractive shell command is shown below. Copy the entire
block, including both `ROOT` delimiter lines:

```bash
root -l -b <<'ROOT'
.L Run_CDet_Qualify_RunSpecificP1_FrozenSample.C+
Run_CDet_Qualify_RunSpecificP1_FrozenSample(
    4344,
    -1,
    "CDet_run4344_projection.conf",
    5,
    3.0,
    0.5,
    2.0,
    0.5,
    1.05);
.q
ROOT
```

The explicit arguments declare, in order:

1. five contiguous accepted-event folds;
2. minimum full-sample residual significance of `3 sigma`;
3. minimum correction across the central 80% ECal-time span of `0.5 ns`;
4. maximum held-out corrected-slope pull of `2 sigma`;
5. minimum held-out absolute-slope reduction of 50%; and
6. maximum candidate/master paired-time RMS ratio of `1.05`.

These are qualification thresholds, not fit starting values. Do not change
them in response to the Run 4344 result without a separately documented review.
HCal is not part of this ECal-singles qualification.

### Freeze the observations

Construct the accepted-pair sample once with:

- the accepted Run 5710 master `p1`;
- no run-specific ECal slope;
- the accepted Run 4344 timing shift already fitted in Step 4;
- exactly the ordinary Step-5 production selection coordinates and cuts;
- the approved Run 4344 configuration; and
- fixed event, hit, pair, layer, pixel, and half-bar identities.

The qualifier records the event, pair, fold, layer, pixel, half-bar, ECal time,
master-slope CDet time, candidate-slope CDet time, and ToT for every frozen hit
in:

```text
CDet_run4344_frozen_p1_samples.tsv
```

It records the full, per-layer, closure, fold-validation, practical-effect,
and resolution results in:

```text
CDet_run4344_frozen_p1_qualification.txt
```

Every master-versus-candidate comparison must use exactly these observations.
Changing `p1` must be an algebraic transformation of the stored CDet times, not
a new event-selection pass. With a fixed sample, the candidate is a one-step
calculation:

```text
p1_candidate = p1_master + residual_p1
```

A material residual after applying that candidate is a model or estimator
diagnostic. It is not a reason to iterate `p1`.

The qualifier deliberately uses the same ordinary population as Step 5. It
must not enable the special shift-invariant cut mode used internally by the
timing-shift fitter, because that would qualify a slope on a different sample
from the final production audit.

### Qualification gates and interpretation

A Run 4344 slope may be written only if the qualifier reports `QUALIFIED: YES`
and the summary shows all of the following on the frozen artifact:

1. the combined fixed-effects residual reaches the declared `3 sigma`
   threshold;
2. the Layer-1 and Layer-2 residuals have the same sign and are statistically
   compatible;
3. the correction across the central 80% ECal-time range is at least `0.5 ns`;
4. each training-fold estimate is within
   `max(0.03 ns/ns, 2*full_sample_error)` of the full estimate;
5. in every held-out contiguous-event fold, the corrected residual is no more
   than `2 sigma` from zero and its magnitude is reduced by at least 50%;
6. the candidate paired-time RMS is no more than 5% larger than the master
   result; and
7. the candidate closes on the frozen sample without a second slope update.

The command prints every fold and gate and ends by confirming that
`CDet_run4344.dat` was not modified. Verify that the run file still contains
only `[GlobalTiming] shift_ns`:

```bash
sed -n '1,120p' CDet_run4344.dat
```

If any gate fails, retain the master `p1`. Do not average candidates, remove a
failed fold, loosen a threshold, or iterate the slope.

If all gates pass, preserve the TSV, summary, terminal output, configuration,
master calibration, and current run file together. Do not add the candidate
manually and do not run:

```cpp
Run_CDet_Calibrate_RunSpecificP1(4344)
```

and do not add `[ECalTiming] p1` to `CDet_run4344.dat` manually.

## 7. Install a qualified candidate and refit the timing origin

Perform this step only when Step 6 reports `QUALIFIED: YES`. In a fresh ROOT
process, run:

```cpp
root -l
.L Run_CDet_Install_QualifiedP1.C+
Run_CDet_Install_QualifiedP1(
    4344,
    "CDet_run4344_frozen_p1_qualification.txt");
.q
```

The installer rereads the qualification artifact, requires its run to be 4344
and `qualified` to be 1, preserves the current run file as
`CDet_run4344_before_qualified_p1.dat`, and atomically installs exactly the
reported candidate. It refuses to overwrite an existing backup.

Inspect the active file and confirm that it contains both the qualified `p1`
and the still-provisional old `shift_ns`:

```bash
sed -n '1,120p' CDet_run4344.dat
```

Because changing `p1` changes the absolute timing origin, immediately refit
`shift_ns` once:

```bash
root -l -b -q 'Run_CDet_Calibrate_RunTimingShift.C+(4344,-1,30.0,5.0,"CDet_run4344_projection.conf")'
```

Require the timing-shift closure and exact before/after accepted-pair population
gates to pass. Do not refit `p1` during or after this operation.

Finally repeat Step 5 in a fresh ROOT process. Record the production-sample
fixed-effects residual and compare its accepted-pair population with the frozen
artifact. The frozen-sample slope is already closed by construction; this last
analysis is an audit for population change after installing the qualified
constant, not another slope-fitting iteration.

## Acceptance status

The restarted master-`p1` result is accepted after the ECal window review,
shift closure, population-invariance check, and fresh final diagnostics are
recorded. The completed read-only frozen-sample test reported `QUALIFIED: NO`
for Run 4344 because the candidate did not generalize across the contiguous
event folds. Retain the master `p1`; do not perform Step 7 for the accepted Run
4344 calibration. Step 7 remains documented only as the controlled procedure
that would apply if a future, separately approved qualification reported
`QUALIFIED: YES`.
