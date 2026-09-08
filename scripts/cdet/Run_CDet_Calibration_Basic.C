#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>

// In-session two-pass CDet calibration driver.
//
// This version does NOT spawn fresh ROOT subprocesses.
// Instead, it resets the master-macro globals between stages and runs
// the full sequence in one ROOT session. This avoids environment/library
// mismatches between parent and child ROOT processes.
//
// Assumptions:
//   - The master macro file is named:
//       PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C
//   - The master macro provides a function:
//       void ResetCalibrationGlobals();
//     If that function does not yet exist, add it there.
//
// Sequence:
//   pass 1: 0 -> 1 -> 3 -> 6
//   pass 2: 1 -> 3 -> 6
//   final : 7

void Run_CDet_Calibration_Basic(
    Int_t RunNumber1 = 3575,
    Int_t nevents = -1,
    Int_t elastic = 0,
    Int_t minSeg = 0,
    Int_t maxSeg = 5,
    Double_t LeMin = 10.0,
    Double_t LeMax = 50.0,
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
    bool removeExistingCalibrationFile = true
){
    const TString masterMacro = "PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C";
    const TString calibFile   = "CDet_calibration_dt.dat";

    if (gSystem->AccessPathName(masterMacro)) {
        std::cerr << "[Driver] ERROR: Could not find master macro " << masterMacro << "\\n";
        return;
    }

    gROOT->ProcessLine(TString::Format(".L %s+", masterMacro.Data()));

    if (removeExistingCalibrationFile && !gSystem->AccessPathName(calibFile)) {
        std::cout << "[Driver] Removing existing calibration file: " << calibFile << "\\n";
        gSystem->Unlink(calibFile);
    }

    auto runMain = [&](Int_t stage) {
        PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget(
            RunNumber1, nevents, stage,elastic, minSeg, maxSeg,
            LeMin, LeMax, TotMin, TotMax,
            10.0, 35.0,
            nhitcutlow1, nhitcuthigh1, nhitcutlow2, nhitcuthigh2,
            XDiffCut, XOffset, YOffset, layer_choice,
            suppress_bad, nruns, maxstream, firstevent
        );
    };

    auto stageBanner = [&](const char* name, Int_t stage) {
        std::cout << "\\n[Driver] ==================================================\\n";
        std::cout << "[Driver] " << name << "  (stage " << stage << ")\\n";
        std::cout << "[Driver] ==================================================\\n";
    };

    stageBanner("final_calibrated_state", 8);
    ResetCalibrationGlobals();
    runMain(8);
    //plotAllTDC(false, 1.0, 0.0, 60.0);
    plotCDetLayersTimeComp(false, 416, 1.0, -15, 15, -0.1, 0.1, 20, 45, 4, 40, -15, 15, 0, 60, 0, 80, 10, 35, -60, 30);
    //plotGoodLeVsTotByLayer(false, 15, 45, 4, 30, 0.2, 0.5, true, false, 5.0, 25.0);


    std::cout << "\\n[Driver] Two-pass in-session calibration sequence complete.\\n";
    std::cout << "[Driver] Final calibration file should be in: " << calibFile << "\\n";
}
