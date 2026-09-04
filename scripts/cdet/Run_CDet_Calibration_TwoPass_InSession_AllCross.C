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
// RunNumber1 = 0 will add all runs from the .txt file runs.txt
void Run_CDet_Calibration_TwoPass_InSession_AllCross(
    Int_t groupIndex = 0,
    Int_t RunNumber1 = 0,
    Int_t nevents = -1,
    Int_t minSeg = 0,
    Int_t maxSeg = 20,
    Double_t LeMin = 10.0,
    Double_t LeMax = 60.0,
    Double_t TotMin = 4.0,
    Double_t TotMax = 50.0,
    Int_t nhitcutlow1 = 0,
    Int_t nhitcuthigh1 = 100,
    Int_t nhitcutlow2 = 0,
    Int_t nhitcuthigh2 = 100,
    Double_t XDiffCut = 0.10,
    Double_t XOffset = 0.02,
    Double_t YOffset = 0.1,
    Int_t layer_choice = 3,
    bool suppress_bad = false,
    Int_t nruns = 30,
    Int_t maxstream = 2,
    Int_t firstevent = 1,
    bool removeExistingCalibrationFile = true
){
    gLastCalibrationSequenceSucceeded = false;
    const TString calibFile   = TString::Format("CDet_calibration_dt_group%d.dat",groupIndex);

    if (removeExistingCalibrationFile && !gSystem->AccessPathName(calibFile)) {
        std::cout << "[Driver] Removing existing calibration file: " << calibFile << "\\n";
        gSystem->Unlink(calibFile);
    }

    auto runMain = [&](Int_t stage) -> bool {
        PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget(
            RunNumber1, nevents, stage,groupIndex, minSeg, maxSeg,
            LeMin, LeMax, TotMin, TotMax,
            10.0, 35.0,
            nhitcutlow1, nhitcuthigh1, nhitcutlow2, nhitcuthigh2,
            XDiffCut, XOffset, YOffset, layer_choice,
            suppress_bad, nruns, maxstream, firstevent, false
        );
        if (!gLastCalibrationStageSucceeded)
            std::cerr << "[Driver] ERROR: stage " << stage << " failed; stopping sequence.\\n";
        return gLastCalibrationStageSucceeded;
    };

    auto stageBanner = [&](const char* name, Int_t stage) {
        std::cout << "\\n[Driver] ==================================================\\n";
        std::cout << "[Driver] " << name << "  (stage " << stage << ")\\n";
        std::cout << "[Driver] ==================================================\\n";
    };

    stageBanner("pass1_raw_inspection", 0);
    ResetCalibrationGlobals();
    if (!runMain(0)) return;
    plotAllTDC(false, 1.0, 0.0, 60.0, false, "", "tdcPlots");
    if (!gLastCalibrationFitSucceeded) return;

    stageBanner("pass1_timeoffset_fit", 1);
    ResetCalibrationGlobals();
    if (!runMain(1)) return;
    extractHierarchicalCDetPixelTimingOffsets(true,
        TString::Format("group%d_hierarchical_pass1_offsets", groupIndex));
    if (!gLastCalibrationFitSucceeded) return;

    stageBanner("pass1_ecal_fit", 3);
    ResetCalibrationGlobals();
    if (!runMain(3)) return;
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
    extractHierarchicalCDetPixelTimingOffsets(true,
        TString::Format("group%d_hierarchical_pass2_offsets", groupIndex));
    if (!gLastCalibrationFitSucceeded) return;

    stageBanner("pass2_ecal_refit", 3);
    ResetCalibrationGlobals();
    if (!runMain(3)) return;
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
    extractHierarchicalCDetPixelTimingOffsets(true,
        TString::Format("group%d_hierarchical_final_closure", groupIndex));
    if (!gLastCalibrationFitSucceeded) return;
    
    stageBanner("final_calibrated_state", 7);
    ResetCalibrationGlobals();
    if (!runMain(7)) return;
    plotAllTDC(false, 1.0, 0.0, 60.0, false, "", "tdcPlots");
    plotCDetLayersTimeComp(false, 416, 1.0, -15, 15, -0.1, 0.1, 20, 45, 4, 40, -15, 15, 0, 60, 0, 80, 10, 35, -60, 30, true, 0.005, -1.5, 1.5, 0.01, 0.0, 7.0);
    plotGoodLeVsTotByLayer(false, 15, 45, 4, 30, 0.2, 0.5, true, false, 5.0, 25.0);


    gLastCalibrationSequenceSucceeded = true;
    std::cout << "\\n[Driver] Two-pass in-session calibration sequence complete.\\n";
    std::cout << "[Driver] Final calibration file should be in: " << calibFile << "\\n";
}
