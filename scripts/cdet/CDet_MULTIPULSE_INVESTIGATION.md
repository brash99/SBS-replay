# CDet multiple-pulse and edge-pairing investigation

## Analysis requirement

Only a complete CDet TDC triplet is eligible for physics analysis. A complete
triplet has a nonzero leading edge (LE), a nonzero trailing edge (TE), and a
positive derived `ToT = TE - LE`. LE-only and TE-only decoder slots are not
pulses and must never enter calibration or event selection.

The new `earm.cdet.pulse.*` collection enforces this requirement. Incomplete
slots remain visible only through the event-level diagnostic counters
`pulse.n_le_only`, `pulse.n_te_only`, and `pulse.n_invalid_pair`.

## Run 5710 sample

The investigation used the two newly replayed segment-0 Run 5710 files:

- `cdet_5710_stream0_2_seg0_0.root`: 753,130 events
- `cdet_5710_stream0_2_seg0_0__1.root`: 292,340 events

Together they contain 1,045,470 events and 4,270,587 CDet decoder slots.

| Decoder-slot classification | Count | Fraction |
| --- | ---: | ---: |
| Complete LE+TE+positive-ToT triplet | 3,060,709 | 71.670% |
| LE only | 423,427 | 9.915% |
| TE only | 786,451 | 18.416% |

Incomplete edge slots are therefore common enough that treating every decoder
slot as a pulse would produce a substantial and avoidable error.

## Multiple complete triplets in one channel

There were 73,931 events (7.072%) in which at least one CDet channel contained
two or more complete triplets. These events contained 76,837 channel
occurrences with multiple complete triplets and 78,327 adjacent
complete-triplet pairs.

| Adjacent-pair class | Count | Fraction |
| --- | ---: | ---: |
| Overlapping: next LE occurs before previous TE | 44,860 | 57.273% |
| Non-overlapping | 33,467 | 42.727% |

Of the overlapping pairs, 31,769 (70.818%) repeat both their LE and TE within
1 ns. No additional pairs enter when that threshold is relaxed from 1 to 2 ns,
showing a pronounced separation between the close-repeat and broader
populations.

For the overlapping population, both input files independently give
approximately the same characteristic separations:

- median absolute LE separation: 0.205 ns;
- median absolute TE separation: 0.299 ns;
- median ToT difference: 0.094 ns.

The same channels dominate both files. The most prominent examples are channel
IDs 802, 1009, 1000, 1012, 1075, 1061, and 918. This stability argues against a
small-sample fluctuation.

## Decoder mechanism

The CDet VETROC decoder writes every hardware hit word with its edge bit and
time. `SBSGenericDetector::DecodeTDC()` then sorts those edges by time.
`SBSData::TDC::Process()` associates them by edge ordinal: the nth LE is paired
with the nth TE. It does not require the sorted sequence to alternate between
LE and TE.

As a result, a time-ordered hardware sequence such as

```text
LE0, LE1, TE0, TE1
```

is exported as two complete but temporally overlapping triplets. This explains
how a triplet can satisfy the basic completeness requirement without
necessarily representing an independent physical discriminator pulse.

The `suppress_hitsperchan` crate-map option only limits warning messages. It
does not suppress, merge, or otherwise alter these data.

The legacy `FindGoodHit()` has a separate problem: it chooses the slot whose LE
is closest to `GoodTimeCut` without first requiring a complete triplet. A
TE-only slot has an LE value of zero and can therefore displace a valid triplet;
the later `SBSCDet` cuts then reject that selected zero-LE slot. The new
all-complete-triplet collection avoids this loss mechanism.

## Adopted CDet pulse-building policy

Completeness is necessary but not sufficient for interpreting multiple entries
in one channel as independent pulses. The adopted CDet-specific policy is an
alternating-edge state machine operating on the time-sorted channel edges:

1. the first LE opens a pulse;
2. additional LEs are ignored while that pulse is open;
3. the first subsequent TE closes and records the pulse;
4. additional TEs are ignored until another LE opens a new pulse.

Therefore `LE0, LE1, TE0, TE1` records exactly one pulse using `(LE0, TE0)`.
This policy is implemented in `SBSCDet::BuildPulseCandidates()` rather than in
the generic decoder, because other detectors may depend on the existing generic
association behavior.

The output preserves separate original decoder-slot indices as
`pulse.le_index` and `pulse.te_index`. Ignored or unmatched edges are counted in
`pulse.n_le_only` and `pulse.n_te_only`; they never enter the pulse arrays.

On the 2,000-event Run 5710 smoke sample, the state machine produced 5,615
accepted pulses, 972 discarded/unmatched LEs, 1,576 discarded/unmatched TEs,
and zero invalid nonpositive-ToT pairs. This accounts exactly for all 13,778
nonzero input edges. The known event 25/channel 802 `LE,LE,TE,TE` example
produced one pulse using the first LE and first TE, as required.

The remaining investigation is to measure the effect of this adopted policy on
electron-hit efficiency and CDet multiplicity when the ECal association is
added.
