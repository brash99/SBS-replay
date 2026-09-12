#include "PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C"

#include <TCanvas.h>
#include <TEnv.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>
#include <limits>

// Run a calibrated analysis without changing any constants and save the
// spectra needed for a human estimate of the run-specific timing shift.
void Run_CDet_Inspect_RunTimingShift(
    Int_t runNumber,
    Int_t nevents = std::numeric_limits<Int_t>::min(),
    TString configFile = "")
{
  if (runNumber <= 0) {
    std::cerr << "[Run timing-shift inspection] ERROR: runNumber must be positive.\n";
    return;
  }
  if (configFile.IsNull())
    configFile = TString::Format("CDet_run%d_projection.conf", runNumber);

  TEnv env;
  if (!LoadCDetConfiguration(env, configFile.Data(),
                             "Run timing-shift inspection")) return;
  if (env.GetValue("analysis.run_number", -1) != runNumber) {
    std::cerr << "[Run timing-shift inspection] ERROR: analysis.run_number in "
              << configFile << " does not match " << runNumber << ".\n";
    return;
  }

  ResetCalibrationGlobals();
  PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget(
      configFile.Data(), 7, nevents);
  if (!gLastCalibrationStageSucceeded) {
    std::cerr << "[Run timing-shift inspection] ERROR: analysis failed.\n";
    return;
  }
  plotCDetLayersTimeComp(configFile.Data(), 0);
  if (!gLastCalibrationFitSucceeded || gCDetAcceptedPairMeanTimes.empty()) {
    std::cerr << "[Run timing-shift inspection] ERROR: timing spectra are unavailable.\n";
    return;
  }

  const TString outputDirectory =
      TString::Format("CDet_run%d_shift_diagnostics", runNumber);
  gSystem->mkdir(outputDirectory, true);

  TCanvas *layerTimes = static_cast<TCanvas *>(gROOT->FindObject("cCDetLayerTimes"));
  if (layerTimes) {
    layerTimes->SaveAs(outputDirectory + "/CDet_layer_LE.png");
    layerTimes->SaveAs(outputDirectory + "/CDet_layer_LE.pdf");
  }
  if (gCDetProjectedHalfBarTimingCanvas) {
    gCDetProjectedHalfBarTimingCanvas->SaveAs(
        outputDirectory + "/CDet_projected_halfbar_timing.png");
    gCDetProjectedHalfBarTimingCanvas->SaveAs(
        outputDirectory + "/CDet_projected_halfbar_timing.pdf");
  }
  TCanvas *selectedBarLeVsTotSource =
      static_cast<TCanvas *>(gROOT->FindObject("cCDetLeVsTotBar"));
  TVirtualPad *selectedBarLayer1Pad =
      selectedBarLeVsTotSource ? selectedBarLeVsTotSource->GetPad(1) : nullptr;
  TH2D *selectedBarLayer1LeVsTot = selectedBarLayer1Pad
      ? dynamic_cast<TH2D *>(selectedBarLayer1Pad->GetPrimitive(
            "hCDet1BarLeVsTot"))
      : nullptr;
  if (selectedBarLayer1LeVsTot) {
    const int selectedBar = env.GetValue("display.pixel", 0) / 16;
    TCanvas selectedBarLeVsTotCanvas(
        "cCDetSelectedBarLayer1LeVsTotInspection",
        "Selected Layer-1 bar LE versus ToT", 1000, 800);
    selectedBarLayer1LeVsTot->Draw("COLZ");
    selectedBarLeVsTotCanvas.SaveAs(
        outputDirectory + TString::Format(
            "/CDet_bar%d_layer1_LE_vs_ToT.png", selectedBar));
    selectedBarLeVsTotCanvas.SaveAs(
        outputDirectory + TString::Format(
            "/CDet_bar%d_layer1_LE_vs_ToT.pdf", selectedBar));
  }

  TH1D pairedMean("hCDetPairedMeanTimeShiftInspection",
      TString::Format("Run %d with current constants;accepted-pair mean CDet time (ns);accepted pairs / 0.5 ns",
                      runNumber),
      240, -30.0, 90.0);
  for (double value : gCDetAcceptedPairMeanTimes) pairedMean.Fill(value);
  TCanvas pairedCanvas("cCDetPairedMeanTimeShiftInspection",
                       "CDet timing-shift inspection", 1500, 700);
  pairedCanvas.Divide(2, 1);
  pairedCanvas.cd(1);
  pairedMean.Draw("HIST");
  pairedCanvas.cd(2);
  gPad->SetLogy();
  pairedMean.SetMinimum(0.7);
  pairedMean.Draw("HIST");
  pairedCanvas.SaveAs(outputDirectory + "/CDet_paired_mean_time.png");
  pairedCanvas.SaveAs(outputDirectory + "/CDet_paired_mean_time.pdf");

  std::cout << "\n[Run timing-shift inspection]\n"
            << "  constants were not modified\n"
            << "  accepted pairs: " << gCDetAcceptedPairMeanTimes.size() << "\n"
            << "  review: " << outputDirectory
            << "/CDet_projected_halfbar_timing.png\n";
}
