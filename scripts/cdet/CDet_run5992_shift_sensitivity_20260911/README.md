# Run 5992 shift sensitivity experiment

Each subdirectory contains an isolated copy of the current Run 5710 master
calibration and the authoritative Run 5992 configuration. The only seeded
run-specific value is the indicated `shift_ns`; no run-specific `p1` is present
before fitting.

Active-file SHA-256 values before the experiment:

```text
CDet_calibration_dt.dat  b39d3deda7bd70bfc4272942a368d333f3fcd527a894050b320709c0ef405600
CDet_run5992.dat         96e12db01bacb3391076ee0708fe4d46488493640d70463e2c6218d1e0f97172
```

The hashes were unchanged after both isolated fits.

## Results

| Seed `shift_ns` | Baseline residual `p1` | Fitted run `p1` | Closure residual `p1` |
|---:|---:|---:|---:|
| 4.203 ns | `0.00912792 +/- 0.0143342` | `0.819331` | `0.00380351 +/- 0.0145059` |
| 6.203 ns | `0.0211177 +/- 0.0130388` | `0.831321` | `0.00314037 +/- 0.0131414` |

Both closure gates passed. The fitted values differ by `0.011990 ns/ns`, less
than either individual statistical uncertainty. The baseline accepted sample
also changed slightly with the seeded shift. A constant additive shift cannot
change the ECal-time slope for an identical sample, so any observed dependence
comes from shift-sensitive event selection and finite-sample fluctuations, not
from minimizing the CDet timing width.

## Commissioning conclusion

This experiment is diagnostic, not the commissioned calibration. A final
shift-only fit retained the Run 5710 master `p1 = 0.810203` and obtained
`shift_ns = 5.171528 ns`. Its closure centroid was
`30.0000 +/- 0.0888 ns`, with `1453 / 1453` accepted pairs and a residual
fixed-effects slope of `-0.0131452 +/- 0.0152886 ns/ns`. Because that slope is
consistent with zero, Run 5992 provides no statistically compelling reason to
replace the high-statistics detector-wide `p1`.
