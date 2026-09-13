#include "Run_CDet_Calibration_TwoPass_InSession_AllCross_Individually.C"

#include <TString.h>
#include <TSystem.h>

#include <iostream>

namespace {
bool Run5710AcceptedStage(const TString &configFile, Int_t nevents, Int_t stage)
{
  ResetCalibrationGlobals();
  PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget(
      configFile.Data(), stage, nevents);
  if (!gLastCalibrationStageSucceeded)
    std::cerr << "[Run 5710 accepted workflow] ERROR: stage " << stage
              << " failed.\n";
  return gLastCalibrationStageSucceeded;
}

bool FitRun5710ECalResidual(const TString &configFile)
{
  plotCDetLayersTimeComp(configFile.Data(), 1);
  return gLastCalibrationFitSucceeded;
}
} // namespace

// Reproducible Run 5710 workflow accepted by matched human review on
// 2026-09-11.  It starts with the clean automatic calibration, applies the
// preserved 102 reviewed pixel cuts, refits the within-half-bar ECal slope,
// and performs one half-bar intercept-alignment pass plus final closure.
void Run_CDet_Calibration_Run5710_Accepted(
    Int_t nevents = std::numeric_limits<Int_t>::min(),
    TString reviewedCutFile =
        "CDet_run5710_halfbar_aligned_final_archive/"
        "CDet_pixel_quality_cuts_run5710_halfbar_aligned_final.root",
    TString configFile = "CDet_run5710_projection.conf")
{
  gLastCalibrationSequenceSucceeded = false;
  if (gSystem->AccessPathName(reviewedCutFile)) {
    std::cerr << "[Run 5710 accepted workflow] ERROR: reviewed cut file not found: "
              << reviewedCutFile << "\n";
    return;
  }

  Run_CDet_Calibration_TwoPass_InSession_AllCross_Individually(
      configFile.Data(), nevents, true);
  if (!gLastCalibrationSequenceSucceeded) return;
  gLastCalibrationSequenceSucceeded = false;

  if (!Run5710AcceptedStage(configFile, nevents, 7) ||
      !FitRun5710ECalResidual(configFile)) return;

  if (!Run5710AcceptedStage(configFile, nevents, 7)) return;
  extractHierarchicalCDetPixelTimingOffsets(
      true, "run5710_reviewed_102_cut_pass", reviewedCutFile);
  if (!gLastCalibrationFitSucceeded) return;

  if (!Run5710AcceptedStage(configFile, nevents, 7) ||
      !FitRun5710ECalResidual(configFile)) return;

  if (!Run5710AcceptedStage(configFile, nevents, 7)) return;
  calibrateCDetHalfBarIntercepts(
      true, 22.0, 100, 1.0,
      "CDet_halfbar_intercept_diagnostics_applied.root",
      "CDet_halfbar_intercept_corrections_applied.dat");
  if (!gLastCalibrationFitSucceeded) return;

  if (!Run5710AcceptedStage(configFile, nevents, 7) ||
      !FitRun5710ECalResidual(configFile)) return;
  calibrateCDetHalfBarIntercepts(
      false, 22.0, 100, 1.0,
      "CDet_halfbar_intercept_final_closure.root",
      "CDet_halfbar_intercept_final_closure.dat");
  if (!gLastCalibrationFitSucceeded) return;

  reportCDetPairedTimeResolution(true, 25.0, 35.0);
  gLastCalibrationSequenceSucceeded = true;
  std::cout << "\n[Run 5710 accepted workflow] Sequence complete.\n"
            << "[Run 5710 accepted workflow] Final constants: "
            << gCalibrationFile << "\n";
}
