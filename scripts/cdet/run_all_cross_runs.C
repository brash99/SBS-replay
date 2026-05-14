
void run_all_cross_runs() {

  // Compile heavy macro
  gROOT->ProcessLine(".L PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C+");

  // Load driver + wrapper (interpreted)
  gROOT->ProcessLine(".L Run_CDet_Calibration_TwoPass_InSession_AllCross.C");

  // Run everything
  gROOT->ProcessLine("Run_CDet_Calibration_TwoPass_InSession_AllCross();");
}