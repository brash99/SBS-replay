#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>
#include <limits>

// Complete commissioned Run 5710 timing-calibration recipe.
//
// Run this from scripts/cdet in a fresh ROOT process, starting with neither
// CDet_calibration_dt.dat nor CDet_run5710.dat present:
//
//   root -l
//   .L Run_CDet_Calibrate_Run5710_FromScratch.C+
//   Run_CDet_Calibrate_Run5710_FromScratch()
//
// The 10--35 ns ECal window was visually approved on 2026-09-11.  This macro
// records that reviewed choice; it does not replace the human review required
// when commissioning a different run or a changed input dataset.
void Run_CDet_Calibrate_Run5710_FromScratch(
    Int_t nevents = std::numeric_limits<Int_t>::min())
{
  const char *masterFile = "CDet_calibration_dt.dat";
  const char *runFile = "CDet_run5710.dat";
  if (!gSystem->AccessPathName(masterFile) ||
      !gSystem->AccessPathName(runFile)) {
    std::cerr
        << "[Run 5710 from scratch] ERROR: clean start required.\n"
        << "  Move both CDet_calibration_dt.dat and CDet_run5710.dat out of "
           "scripts/cdet before running this macro.\n";
    return;
  }

  // 1. Build the accepted detector-wide master calibration from raw Run 5710
  //    using the preserved 102-pixel human-reviewed cut artifact.
  gROOT->ProcessLine(".L Run_CDet_Calibration_Run5710_Accepted.C+");
  gROOT->ProcessLine(TString::Format(
      "Run_CDet_Calibration_Run5710_Accepted(%d);", nevents));
  if (!static_cast<bool>(gROOT->ProcessLine(
          "gLastCalibrationSequenceSucceeded"))) {
    std::cerr << "[Run 5710 from scratch] ERROR: accepted master workflow "
                 "failed.\n";
    return;
  }

  // 2. Regenerate the ECal-window proposal for the audit trail.  The approved
  //    10--35 ns selection is authoritative in CDet_run5710_projection.conf;
  //    calibration constants, not selection cuts, belong in the .dat file.
  gROOT->ProcessLine(".L Run_CDet_Select_ECalTimingWindow.C+");
  gROOT->ProcessLine("Run_CDet_Select_ECalTimingWindow(5710);");

  // 3. Determine s_5710 from the final Stage-7 pair-time distribution and
  //    require centroid, ECal-slope, and accepted-population closure.
  gROOT->ProcessLine(".L Run_CDet_Calibrate_RunTimingShift.C+");
  gROOT->ProcessLine(TString::Format(
      "Run_CDet_Calibrate_RunTimingShift(5710, %d, 30.0, 5.0);", nevents));
  if (!static_cast<bool>(gROOT->ProcessLine(
          "gLastCalibrationSequenceSucceeded"))) {
    std::cerr << "[Run 5710 from scratch] ERROR: timing-shift closure failed.\n";
    return;
  }

  std::cout
      << "\n[Run 5710 from scratch] COMPLETE\n"
      << "  master calibration: " << masterFile << "\n"
      << "  run-specific timing constants: " << runFile << "\n";
}
