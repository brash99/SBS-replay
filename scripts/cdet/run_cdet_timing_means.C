#include <TROOT.h>

void run_cdet_timing_means() {
  gROOT->ProcessLine(
      ".L PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C+");
  gROOT->ProcessLine(".L Run_CDet_Extract_RunTimingMeans_FromList.C");
  gROOT->ProcessLine("Run_CDet_Extract_RunTimingMeans_FromList(\"runs.txt\");");
}
