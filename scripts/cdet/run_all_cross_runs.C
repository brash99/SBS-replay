//macro to run all cross runs in a single session
void run_all_cross_runs(Int_t nGroups) {

  // Compile heavy macro
  gROOT->ProcessLine(".L PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C+");

  // Load driver + wrapper (interpreted)
  gROOT->ProcessLine(".L Run_CDet_Calibration_TwoPass_InSession_AllCross.C");

  // Run everything
  for (Int_t i = 0; i < nGroups; i++) gROOT->ProcessLine(TString::Format("Run_CDet_Calibration_TwoPass_InSession_AllCross(%d);", i));
}