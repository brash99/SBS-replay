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
    TString configFile = "",
    TString outputTag = "")
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

  TString outputDirectory =
      TString::Format("CDet_run%d_shift_diagnostics", runNumber);
  if (!outputTag.IsNull()) outputDirectory += "_" + outputTag;
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
  TCanvas *ecalTiming =
      static_cast<TCanvas *>(gROOT->FindObject("cCDetTvsECalT"));
  if (ecalTiming) {
    ecalTiming->SaveAs(outputDirectory + "/cCDetTvsECalT.png");
    ecalTiming->SaveAs(outputDirectory + "/cCDetTvsECalT.pdf");
  }
  TCanvas *ecalTimingDiagnostics = static_cast<TCanvas *>(
      gROOT->FindObject("cCDetTvsECalTDiagnostics"));
  if (ecalTimingDiagnostics) {
    ecalTimingDiagnostics->SaveAs(
        outputDirectory + "/cCDetTvsECalTDiagnostics.png");
    ecalTimingDiagnostics->SaveAs(
        outputDirectory + "/cCDetTvsECalTDiagnostics.pdf");
  }
  TCanvas *hcalTiming =
      static_cast<TCanvas *>(gROOT->FindObject("cCDetTvsHCalT"));
  if (hcalTiming) {
    hcalTiming->SaveAs(outputDirectory + "/cCDetTvsHCalT.png");
    hcalTiming->SaveAs(outputDirectory + "/cCDetTvsHCalT.pdf");
  }
  TCanvas *timingStructures = static_cast<TCanvas *>(
      gROOT->FindObject("cCDetSelectedBarTimingStructures"));
  if (timingStructures) {
    timingStructures->SaveAs(
        outputDirectory + "/CDet_bar30_vertical_vs_diagonal_LE_vs_ToT.png");
    timingStructures->SaveAs(
        outputDirectory + "/CDet_bar30_vertical_vs_diagonal_LE_vs_ToT.pdf");
  }
  TCanvas *allBarsTimingStructures = static_cast<TCanvas *>(
      gROOT->FindObject("cCDetAllBarsTimingStructures"));
  if (allBarsTimingStructures) {
    allBarsTimingStructures->SaveAs(
        outputDirectory + "/CDet_allbars_three_region_LE_vs_ToT.png");
    allBarsTimingStructures->SaveAs(
        outputDirectory + "/CDet_allbars_three_region_LE_vs_ToT.pdf");
  }
  const int selectedBar = env.GetValue("display.pixel", 0) / 16;
  TCanvas *selectedBarLeVsTotSource =
      static_cast<TCanvas *>(gROOT->FindObject("cCDetLeVsTotBar"));
  TVirtualPad *selectedBarLayer1Pad =
      selectedBarLeVsTotSource ? selectedBarLeVsTotSource->GetPad(1) : nullptr;
  TH2D *selectedBarLayer1LeVsTot = selectedBarLayer1Pad
      ? dynamic_cast<TH2D *>(selectedBarLayer1Pad->GetPrimitive(
            "hCDet1BarLeVsTot"))
      : nullptr;
  if (selectedBarLayer1LeVsTot) {
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

  const double correctedLeMin = env.GetValue("display.cdet_time_min", 0.0);
  const double correctedLeMax = env.GetValue("display.cdet_time_max", 60.0);
  const double correctedLeBinWidth =
      env.GetValue("display.histogram_width", 1.0);
  const int correctedLeBins = static_cast<int>(
      (correctedLeMax - correctedLeMin) / correctedLeBinWidth);
  if (selectedBar >= 0 && selectedBar < NumCDetPaddles / NumPaddles &&
      correctedLeBins > 0 &&
      vPaddleGoodLe.size() == static_cast<std::size_t>(NumCDetPaddles)) {
    TCanvas selectedBarCorrectedLeCanvas(
        "cCDetSelectedBarCorrectedLeInspection",
        TString::Format("Run %d, Layer 1 bar %d corrected LE", runNumber,
                        selectedBar),
        1600, 1000);
    selectedBarCorrectedLeCanvas.Divide(4, 4, 0.001, 0.001);
    std::vector<TH1D *> correctedLeHistograms;
    correctedLeHistograms.reserve(NumPaddles);
    for (int localPixel = 0; localPixel < NumPaddles; ++localPixel) {
      const int pixel = selectedBar * NumPaddles + localPixel;
      TH1D *histogram = new TH1D(
          TString::Format("hCDetRun%dBar%dPixel%dCorrectedLe", runNumber,
                          selectedBar, pixel),
          TString::Format("Pixel %d;Fully corrected LE time (ns);Good hits / %.2g ns",
                          pixel, correctedLeBinWidth),
          correctedLeBins, correctedLeMin, correctedLeMax);
      histogram->SetDirectory(nullptr);
      for (double value : vPaddleGoodLe[pixel]) histogram->Fill(value);
      selectedBarCorrectedLeCanvas.cd(localPixel + 1);
      histogram->Draw("HIST");
      correctedLeHistograms.push_back(histogram);
    }
    selectedBarCorrectedLeCanvas.SaveAs(
        outputDirectory + TString::Format(
            "/CDet_bar%d_layer1_corrected_LE_4x4.png", selectedBar));
    selectedBarCorrectedLeCanvas.SaveAs(
        outputDirectory + TString::Format(
            "/CDet_bar%d_layer1_corrected_LE_4x4.pdf", selectedBar));
  }

  // Reproduce the selective ECal-CDet timing views used during the earlier
  // hydrogen study.  This is diagnostic only: do not write a candidate table,
  // and deliberately disable saved pixel polygons so the configured global
  // LH2 ToT interval is the quality selection applied to every hit.
  const double ecalEnergyMin = env.GetValue("analysis.ecal_energy_min", 3.0);
  const double ecalEnergyMax = env.GetValue("analysis.ecal_energy_max", 4.5);
  const double acceptedTotMin = env.GetValue("analysis.tot_min", 8.0);
  const double acceptedTotMax = env.GetValue("analysis.tot_max", 35.0);
  extractCDetBarPixelTimingOffsets(
      selectedBar * NumPaddles,
      1.0,                         // histogram bin width (ns)
      -60.0, 30.0,                // displayed tECal-tCDet range (ns)
      -45.0, -10.0,               // fit range for the current candidate
      50, 0.5, 20.0, 10.0, 1.0,  // fit-quality requirements
      true, outputDirectory + "/bar_ecal_cdet_timing",
      false, "",                 // never write calibration candidates
      ecalEnergyMin, ecalEnergyMax,
      acceptedTotMin, acceptedTotMax,
      8.0, 2.5, -40.0, -15.0,
      true,
      acceptedTotMin, acceptedTotMax,
      "");                       // no polygon gate in run analysis

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
