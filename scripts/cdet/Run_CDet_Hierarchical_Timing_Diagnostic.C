#include <TROOT.h>
#include <TSystem.h>
#include <iostream>

// Rebuild the uncalibrated run-5710 event vectors and run the non-production
// hierarchical eight-pixel timing-offset diagnostic in the same ROOT process.
// The diagnostic writes only *_hierarchical_* products and never activates
// its candidate calibration.
void Run_CDet_Hierarchical_Timing_Diagnostic()
{
    if (gROOT->ProcessLine(
            ".L PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C+") != 0) {
        std::cerr << "[Hierarchical diagnostic driver] ERROR: master macro failed to load.\n";
        gSystem->Exit(2);
        return;
    }

    gROOT->ProcessLine(
        "PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget("
        "5710,-1,0,0,-1,-1,0.02,60.0,4.0,50.0,10.0,35.0,"
        "1,100,1,100,0.10,0.02,0.1,3,false,30,2,1,false);");

    gROOT->ProcessLine(
        "extractHierarchicalCDetPixelTimingOffsetsDiagnostic("
        "1.0,-60.0,30.0,-30.0,10.0,35,100,0.5,8.0,15.0,2.0,1.5,4.0,8.0,"
        "1.0,12.0,"
        "\"CDet_calibration_dt_hierarchical_stage0_candidate.dat\","
        "\"CDet_pixel_timing_fit_results_hierarchical_stage0.dat\","
        "\"CDet_pixel_timing_hierarchical_stage0.root\","
        "\"tdcPlots/hierarchical_stage0\","
        "\"29,79,104,118,139,148\");");
}
