#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include "PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C"

// In-session two-pass CDet calibration driver.
//
// This version does NOT spawn fresh ROOT subprocesses.
// Instead, it resets the master-macro globals between stages and runs
// the full sequence in one ROOT session. This avoids environment/library
// mismatches between parent and child ROOT processes.
//
// Compile this driver with ACLiC (.L ...C+) so the master macro and driver are
// built into one self-contained library.
//
// Sequence:
//   pass 1: 0 -> 1 -> 3 -> 6
//   pass 2: 1 -> 3 -> 6
//   final : 7

void Run_CDet_Calibration_TwoPass_InSession_AllCross_Individually(
    Int_t RunNumber1 = 3575,
    Int_t nevents = -1,
    Int_t elastic = 0,
    Int_t minSeg = 0,
    Int_t maxSeg = 5,
    Double_t LeMin = 0.02,
    Double_t LeMax = 60.0,
    Double_t TotMin = 4.0,
    Double_t TotMax = 50.0,
    Int_t nhitcutlow1 = 1,
    Int_t nhitcuthigh1 = 100,
    Int_t nhitcutlow2 = 1,
    Int_t nhitcuthigh2 = 100,
    Double_t XDiffCut = 0.10,
    Double_t XOffset = 0.02,
    Double_t YOffset = 0.1,
    Int_t layer_choice = 3,
    bool suppress_bad = false,
    Int_t nruns = 30,
    Int_t maxstream = 2,
    Int_t firstevent = 1,
    bool removeExistingCalibrationFile = true,
    TString authoritativeConfig = ""
){
    gLastCalibrationSequenceSucceeded = false;
    const TString calibFile   = "CDet_calibration_dt.dat";

    if (removeExistingCalibrationFile && !gSystem->AccessPathName(calibFile)) {
        std::cout << "[Driver] Removing existing calibration file: " << calibFile << "\\n";
        gSystem->Unlink(calibFile);
    }
    if (removeExistingCalibrationFile) {
        // A missing file alone is insufficient: the master macro has legacy
        // in-memory defaults.  A clean run must invalidate every correction
        // explicitly so Stage 1 cannot write those defaults as if fitted.
        gPixelToffsetCorr.assign(NumCDetPaddles, 0.0);
        gPixelToffsetNhits.assign(NumCDetPaddles, 0);
        gPixelToffsetLoaded = false;
        gECalFitP0 = 0.0;
        gECalFitP1 = 0.0;
        gECalDeltaShift = 0.0;
        gECalParamsLoaded = false;
        gECalDeltaLoaded = false;
        gTimeWalkP1_L1 = 0.0;
        gTimeWalkP1_L2 = 0.0;
        gTimeWalkParamsLoaded = false;
    }

    auto runMain = [&](Int_t stage) -> bool {
        if (!authoritativeConfig.IsNull()) {
          PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget(
              authoritativeConfig.Data(), stage, nevents);
        } else {
          PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget(
            RunNumber1, nevents, stage,elastic, minSeg, maxSeg,
            LeMin, LeMax, TotMin, TotMax,
            10.0, 35.0,
            nhitcutlow1, nhitcuthigh1, nhitcutlow2, nhitcuthigh2,
            XDiffCut, XOffset, YOffset, layer_choice,
            suppress_bad, nruns, maxstream, firstevent, false
          );
        }
        if (!gLastCalibrationStageSucceeded)
            std::cerr << "[Driver] ERROR: stage " << stage << " failed; stopping sequence.\\n";
        return gLastCalibrationStageSucceeded;
    };

    auto stageBanner = [&](const char* name, Int_t stage) {
        std::cout << "\\n[Driver] ==================================================\\n";
        std::cout << "[Driver] " << name << "  (stage " << stage << ")\\n";
        std::cout << "[Driver] ==================================================\\n";
    };
    
    stageBanner("before_pixel_offsets", 0);
    ResetCalibrationGlobals();
    if (!runMain(0)) return;
    plotAllPaddles(1, 0, 60, 0, 60, 0, 60, TString::Format("%d", RunNumber1));
    runStats();
    // plotAllTDC(false, 1.0, 0.0, 60.0, true,
    //         "stage0_beforeOffsets",
    //        TString::Format("tdcPlots/run%d", RunNumber1));

    stageBanner("pass1_timeoffset_fit", 1);
    ResetCalibrationGlobals();
    if (!runMain(1)) return;
    extractHierarchicalCDetPixelTimingOffsets(true, "hierarchical_pass1_offsets");
    if (!gLastCalibrationFitSucceeded) return;

    // Apply offsets and then plot corrected spectra
    stageBanner("after_pixel_offsets_applied", 2);
    ResetCalibrationGlobals();
    if (!runMain(2)) return;
    plotAllTDC(false, 1.0, 0.0, 60.0, true,
           "stage2_afterOffsets",
           TString::Format("tdcPlots/run%d", RunNumber1));
    if (!gLastCalibrationFitSucceeded) return;

    stageBanner("pass1_absolute_ecal_fit", 8);
    ResetCalibrationGlobals();
    if (!runMain(8)) return;
    if (!authoritativeConfig.IsNull())
      plotCDetLayersTimeComp(authoritativeConfig.Data(), 1);
    else
      plotCDetLayersTimeComp(true, 416, 1.0, -15, 15, -0.1, 0.1, 20, 45, 4, 40, -15, 15, 0, 60, 0, 80, 10, 35, -60, 30, true, 0.005, -1.5, 1.5, 0.01, 0.0, 7.0);
    if (!gLastCalibrationFitSucceeded) return;

    stageBanner("pass1_timewalk_fit", 6);
    ResetCalibrationGlobals();
    if (!runMain(6)) return;
    plotGoodLeVsTotByLayer(true, 15, 45, 4, 30, 0.2, 0.5, true, true, 5.0, 25.0);
    if (!gLastCalibrationFitSucceeded) return;

    stageBanner("pass2_timeoffset_refit", 1);
    ResetCalibrationGlobals();
    if (!runMain(1)) return;
    extractHierarchicalCDetPixelTimingOffsets(true, "hierarchical_pass2_offsets");
    if (!gLastCalibrationFitSucceeded) return;

    stageBanner("pass2_ecal_refit", 3);
    ResetCalibrationGlobals();
    if (!runMain(3)) return;
    if (!authoritativeConfig.IsNull())
      plotCDetLayersTimeComp(authoritativeConfig.Data(), 1);
    else
      plotCDetLayersTimeComp(true, 416, 1.0, -15, 15, -0.1, 0.1, 20, 45, 4, 40, -15, 15, 0, 60, 0, 80, 10, 35, -60, 30, true, 0.005, -1.5, 1.5, 0.01, 0.0, 7.0);
    if (!gLastCalibrationFitSucceeded) return;

    stageBanner("pass2_timewalk_refit", 6);
    ResetCalibrationGlobals();
    if (!runMain(6)) return;
    plotGoodLeVsTotByLayer(true, 15, 45, 4, 30, 0.2, 0.5, true, true, 5.0, 25.0);
    if (!gLastCalibrationFitSucceeded) return;

    // -----------------------------
    // FINAL PIXEL-OFFSET CLOSURE PASS
    // -----------------------------
    stageBanner("pass3_fullclosure_offsets", 7);
    ResetCalibrationGlobals();
    if (!runMain(7)) return;
    extractHierarchicalCDetPixelTimingOffsets(true, "hierarchical_final_closure");
    if (!gLastCalibrationFitSucceeded) return;
    
    stageBanner("final_calibrated_state", 7);
    ResetCalibrationGlobals();
    if (!runMain(7)) return;
    plotAllTDC(false, 1.0, 0.0, 60.0, false, "", "tdcPlots");
    if (!authoritativeConfig.IsNull())
      plotCDetLayersTimeComp(authoritativeConfig.Data(), 0);
    else
      plotCDetLayersTimeComp(false, 416, 1.0, -15, 15, -0.1, 0.1, 20, 45, 4, 40, -15, 15, 0, 60, 0, 80, 10, 35, -60, 30, true, 0.005, -1.5, 1.5, 0.01, 0.0, 7.0);
    if (!gLastCalibrationFitSucceeded) return;
    plotGoodLeVsTotByLayer(false, 15, 45, 4, 30, 0.2, 0.5, true, false, 5.0, 25.0);
    if (!gLastCalibrationFitSucceeded) return;
    gLastCalibrationSequenceSucceeded = true;
    std::cout << "\\n[Driver] Two-pass in-session calibration sequence complete.\\n";
    std::cout << "[Driver] Final calibration file should be in: " << calibFile << "\\n";
}

// Configuration-authoritative production entry point.  Only the calibration
// stage and optional event limit are transient; every cut and display setting
// comes from configFile.
void Run_CDet_Calibration_TwoPass_InSession_AllCross_Individually(
    const char *configFile, Int_t nevents = std::numeric_limits<Int_t>::min(),
    bool removeExistingCalibrationFile = true)
{
    TEnv env;
    if (!LoadCDetConfiguration(env, configFile,
                               "CDet configuration calibration driver"))
      return;
    const Int_t runNumber = env.GetValue("analysis.run_number", -1);
    const Int_t configuredEvents = env.GetValue("analysis.events", -1);
    const Int_t effectiveEvents = nevents == std::numeric_limits<Int_t>::min()
                                      ? configuredEvents : nevents;
    if (runNumber <= 0) {
      std::cerr << "[Driver] ERROR: invalid analysis.run_number in "
                << configFile << ".\n";
      return;
    }
    Run_CDet_Calibration_TwoPass_InSession_AllCross_Individually(
        runNumber, effectiveEvents, 0, 0, 5, 0.02, 60.0, 4.0, 50.0,
        1, 100, 1, 100, 0.10, 0.02, 0.1, 3, false, 30, 2, 1,
        removeExistingCalibrationFile, configFile);
}
