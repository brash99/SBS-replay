#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>

// In-session hydrogen calibration/update driver.
//
// Intended workflow for hydrogen:
//   1) Fit UPDATED CDet/ECal timing in a stage where
//        - bar offsets are applied
//        - timewalk is applied
//        - existing ECal correction is NOT applied
//      i.e. stage 8
//   2) Run final fully-calibrated stage 7
//
// Assumptions:
//   - The master macro file is named:
//       PlotElastic_Calibration_Master_stageflag_singlefile_hydrogen.C
//   - The master macro provides:
//       void ResetCalibrationGlobals();
//   - Stage 8 exists in the master macro.

void Run_CDet_Calibration_Hydrogen(
    Int_t RunNumber1 = 6077,
    Int_t nevents = -1,
    Int_t elastic = 0,
    Int_t minSeg = 0,
    Int_t maxSeg = 5,
    Double_t LeMin = 1.0,
    Double_t LeMax = 60.0,
    Double_t TotMin = 8.0,
    Double_t TotMax = 30.0,
    Int_t nhitcutlow1 = 1,
    Int_t nhitcuthigh1 = 100,
    Int_t nhitcutlow2 = 1,
    Int_t nhitcuthigh2 = 100,
    Double_t XDiffCut = 0.04,
    Double_t XOffset = 0.02,
    Double_t YOffset = 0.1,
    Int_t layer_choice = 3,
    bool suppress_bad = false,
    Int_t nruns = 30,
    Int_t maxstream = 2,
    Int_t firstevent = 1,
    bool removeExistingCalibrationFile = false
){
    const TString masterMacro = "PlotElastic_Calibration_Master_stageflag_singlefile_hydrogen.C";
    const TString calibFile   = "CDet_calibration.dat";

    if (gSystem->AccessPathName(masterMacro)) {
        std::cerr << "[Driver] ERROR: Could not find master macro " << masterMacro << "\n";
        return;
    }

    gROOT->ProcessLine(TString::Format(".L %s+", masterMacro.Data()));

    if (removeExistingCalibrationFile && !gSystem->AccessPathName(calibFile)) {
        std::cout << "[Driver] Removing existing calibration file: " << calibFile << "\n";
        gSystem->Unlink(calibFile);
    }

    auto runMain = [&](Int_t stage) {
        PlotElastic_Calibration_Master_stageflag_singlefile_hydrogen(
            RunNumber1, nevents, stage, elastic, minSeg, maxSeg,
            LeMin, LeMax, TotMin, TotMax,
            nhitcutlow1, nhitcuthigh1, nhitcutlow2, nhitcuthigh2,
            XDiffCut, XOffset, YOffset, layer_choice,
            suppress_bad, nruns, maxstream, firstevent
        );
    };

    auto stageBanner = [&](const char* name, Int_t stage) {
        std::cout << "\n[Driver] ==================================================\n";
        std::cout << "[Driver] " << name << "  (stage " << stage << ")\n";
        std::cout << "[Driver] ==================================================\n";
    };


/*
    // First pass on hydrogen:
    // Fit CDet/ECal timing with bar offsets + timewalk applied, but NO ECal correction.
    stageBanner("hydrogen_ecal_fit", 8);
    ResetCalibrationGlobals();
    runMain(8);
    plotCDetLayersTimeComp(true, 1.0, -15, 15, -0.05, 0.05,
                           LeMin, LeMax, TotMin, TotMax, -15, 15, 0, 60, 0, 80,
                           70, 125, -160, -40);
*/

  
    // Final fully calibrated pass on hydrogen
    stageBanner("final_calibrated_state", 7);
    ResetCalibrationGlobals();
    runMain(7);
    plotAllTDC(false, 1.0, 1, 60);
    plotCDetLayersTimeComp(false, 1.0, -15, 15, -0.02, 0.02,
                           LeMin, LeMax, TotMin, TotMax, -15, 15, 0, 60, 0, 80,
                           85, 115, -160, -40);
    plotGoodLeVsTotByLayer(false, LeMin, LeMax, TotMin, TotMax, 0.2, 0.5, true, false, 5.0, 25.0);


    std::cout << "\n[Driver] Hydrogen calibration/update sequence complete.\n";
    std::cout << "[Driver] Final calibration file should be in: " << calibFile << "\n";

}
