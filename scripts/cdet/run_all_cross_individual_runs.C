void run_all_cross_individual_runs() {

  // Compile heavy macro
  gROOT->ProcessLine(".L PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C+");

  // Load driver + wrapper (interpreted)
  gROOT->ProcessLine(".L Run_CDet_Calibration_TwoPass_InSession_AllCross_Individually.C");
  gROOT->ProcessLine(".L Run_CDet_Calibration_FromList.C");

  // Run everything
  gROOT->ProcessLine("Run_CDet_Calibration_FromList(\"crossRuns.txt\", false);");
}