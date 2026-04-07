# CDet Timing Calibration Log

## Overview

This document summarizes the development and implementation of timing
corrections for CDet, including:

-   ECal-based timing correction
-   Time-walk correction using TOT
-   Iterative refinement of calibration parameters

------------------------------------------------------------------------

## 1. Initial Goal

Improve timing resolution and remove systematic dependencies in:

-   CDet LE time vs. TOT (time-walk)
-   CDet time vs. ECal time (residual slope)

------------------------------------------------------------------------

## 2. Time-Walk Investigation

### Method

Studied: LE (corrected for offsets + ECal) vs TOT

Separated by detector layer: - Layer 1 - Layer 2

### Observation

Clear dependence: - LE decreases with increasing TOT - Functional form
consistent with:

LE(TOT) = p0 + p1 / sqrt(TOT)

------------------------------------------------------------------------

## 3. Time-Walk Fit

Fit applied using:

f(TOT) = p0 + p1 / sqrt(TOT)

Independently for: - Layer 1 → p1_L1 - Layer 2 → p1_L2

------------------------------------------------------------------------

## 4. Time-Walk Correction Implementation

### Correction Formula

t_corr = t - p1 \* (1/sqrt(TOT) - 1/sqrt(TOT_ref))

### Code Implementation

``` cpp
double tw_corr = 0.0;

if (TOT > gTimeWalkTotMin && TOT < gTimeWalkTotMax) {
    if (layer == 1) {
        tw_corr = gTimeWalkP1_L1 * (1.0/sqrt(TOT) - 1.0/sqrt(gTimeWalkTotRef_L1));
    } else if (layer == 2) {
        tw_corr = gTimeWalkP1_L2 * (1.0/sqrt(TOT) - 1.0/sqrt(gTimeWalkTotRef_L2));
    }
}

LEcorr -= tw_corr;
TEcorr -= tw_corr;
```

------------------------------------------------------------------------

## 5. Secondary Effect: Residual ECal Dependence

After time-walk correction, examined:

CDet time vs ECal time

Fit result:

t_CDet = 7.425 + 0.19742 \* t_ECal

------------------------------------------------------------------------

## 6. ECal Correction Update

Original: P0 = -40.303 P1 = 0.63615

Updated: P1 = 0.83357

------------------------------------------------------------------------

## 7. Final Correction Chain

Raw → Offsets → ECal → Time-walk → Global shift

------------------------------------------------------------------------

## Summary

-   Time-walk removed using 1/sqrt(TOT)
-   Residual ECal slope corrected
-   Final P1 = 0.83357
