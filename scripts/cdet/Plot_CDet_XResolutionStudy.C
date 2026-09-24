#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "CDetGoodPulseConfig.h"
#include "CDetRunDataset.h"

namespace {

constexpr double kECalZFromTargetM = 6.144;
constexpr int kCDetNPixels = 2688;
constexpr double kPairDeltaTimeMaxNs = 15.0;
constexpr double kPairDeltaXMaxM = 0.15;
constexpr double kPairDeltaYMaxM = 0.08;
TH1D *gCDetXResolutionSidebandTemplate = nullptr;

double SignalPlusSidebandTemplate(double *x, double *par) {
  const double sigma = par[2];
  if (!(sigma > 0.0) || !gCDetXResolutionSidebandTemplate)
    return 0.0;
  const double binWidth = gCDetXResolutionSidebandTemplate->GetBinWidth(1);
  const double gaussian =
      par[0] * binWidth * TMath::Gaus(x[0], par[1], sigma, true);
  const double background =
      par[3] * gCDetXResolutionSidebandTemplate->Interpolate(x[0]);
  return gaussian + background;
}

double SidebandTemplateOnly(double *x, double *par) {
  if (!gCDetXResolutionSidebandTemplate)
    return 0.0;
  return par[0] * gCDetXResolutionSidebandTemplate->Interpolate(x[0]);
}

double GaussianOnly(double *x, double *par) {
  if (!gCDetXResolutionSidebandTemplate || !(par[2] > 0.0))
    return 0.0;
  const double binWidth = gCDetXResolutionSidebandTemplate->GetBinWidth(1);
  return par[0] * binWidth * TMath::Gaus(x[0], par[1], par[2], true);
}

bool HasRequiredBranches(TChain &chain) {
  const char *required[] = {
      "earm.cdet.pulse.pmtnum", "earm.cdet.pulse.tdc_le_corr",
      "earm.cdet.pulse.ecal_residual", "earm.cdet.pulse.calib_valid",
      "earm.cdet.pulse.broad_quality_pass",
      "earm.cdet.pulse.ecal_eligible", "earm.cdet.pulse.spatial_pass",
      "earm.cdet.pulse.x_corr", "earm.cdet.pulse.y",
      "earm.cdet.pulse.z", "earm.ecal.adctime", "earm.ecal.e",
      "earm.ecal.x", "earm.ecal.y"};
  chain.LoadTree(0);
  for (const char *name : required) {
    if (!chain.GetBranch(name)) {
      std::cerr << "[CDet x resolution] Missing required branch: " << name
                << '\n';
      return false;
    }
  }
  return true;
}

} // namespace

// Estimate the inter-layer CDet x residual width without cutting on that
// residual. Analyzer-stored pairs have already passed the pulse-quality and
// hard inter-layer gates. The ECal-CDet pair timing selects a signal band and
// two sidebands; the sideband residual shape supplies the accidental template
// in a Gaussian-signal-plus-template fit.
void Plot_CDet_XResolutionStudy(
    const char *configFile = "CDet_run6077_projection.conf",
    const char *inputDirectory = nullptr,
    const char *outputDirectory = "CDet_run6077_x_resolution",
    double signalHalfWidthNs = 5.0, double sidebandInnerNs = 10.0,
    double sidebandOuterNs = 15.0, double residualMinM = -0.20,
    double residualMaxM = 0.20, int residualBins = 200,
    double fitMinM = -0.03, double fitMaxM = 0.03) {
  CDetGoodPulseConfig::Values config;
  if (!CDetGoodPulseConfig::Load(configFile, config,
                                 "CDet x-resolution study"))
    return;

  TString input(inputDirectory ? inputDirectory : "");
  if (input.IsNull()) {
    const char *outDir = gSystem->Getenv("OUT_DIR");
    if (outDir)
      input = outDir;
  }
  if (input.IsNull()) {
    std::cerr << "[CDet x resolution] An input directory or OUT_DIR is "
                 "required.\n";
    return;
  }
  if (!(signalHalfWidthNs > 0.0) ||
      !(sidebandOuterNs > sidebandInnerNs) ||
      !(sidebandInnerNs > signalHalfWidthNs) ||
      !(residualMaxM > residualMinM) || residualBins < 20 ||
      !(fitMaxM > fitMinM) || fitMinM < residualMinM ||
      fitMaxM > residualMaxM) {
    std::cerr << "[CDet x resolution] Invalid signal, sideband, histogram, or "
                 "fit interval.\n";
    return;
  }

  TChain chain("T");
  const int filesAdded =
      CDetRunDataset::AddToChain(&chain, config.runNumber, input.Data());
  if (filesAdded <= 0 || chain.GetEntries() <= 0) {
    std::cerr << "[CDet x resolution] No events found for run "
              << config.runNumber << " in " << input << '\n';
    return;
  }
  if (!HasRequiredBranches(chain))
    return;

  TH2D hTimingVsResidual(
      "hCDetXResolutionTimingVsResidual",
      "Hard-gate CDet pairs before trajectory-time ellipse;"
      "r_{x} = #Delta x_{CDet} - (x_{ECal}/z_{ECal})#Delta z (m);"
      "t_{ECal} - <t_{CDet,corr}> (ns)",
      residualBins, residualMinM, residualMaxM, 160,
      config.pairTimingCenterNs - 40.0, config.pairTimingCenterNs + 40.0);
  TH1D hSignal("hCDetXResolutionSignalRegion",
               "Timing signal region;r_{x} (m);Pairs", residualBins,
               residualMinM, residualMaxM);
  TH1D hSideband("hCDetXResolutionSidebands",
                 "Timing sidebands;r_{x} (m);Pairs", residualBins,
                 residualMinM, residualMaxM);

  TTreeReader reader(&chain);
  TTreeReaderArray<Double_t> pixel(reader, "earm.cdet.pulse.pmtnum");
  TTreeReaderArray<Double_t> le(reader, "earm.cdet.pulse.tdc_le_corr");
  TTreeReaderArray<Double_t> pulseZ(reader, "earm.cdet.pulse.z");
  TTreeReaderArray<Double_t> ecalResidual(
      reader, "earm.cdet.pulse.ecal_residual");
  TTreeReaderArray<Double_t> calibValid(reader,
                                        "earm.cdet.pulse.calib_valid");
  TTreeReaderArray<Double_t> broadQuality(
      reader, "earm.cdet.pulse.broad_quality_pass");
  TTreeReaderArray<Double_t> ecalEligible(
      reader, "earm.cdet.pulse.ecal_eligible");
  TTreeReaderArray<Double_t> spatialPass(
      reader, "earm.cdet.pulse.spatial_pass");
  TTreeReaderArray<Double_t> correctedX(reader, "earm.cdet.pulse.x_corr");
  TTreeReaderArray<Double_t> pulseY(reader, "earm.cdet.pulse.y");
  TTreeReaderValue<Double_t> ecalTime(reader, "earm.ecal.adctime");
  TTreeReaderValue<Double_t> ecalEnergy(reader, "earm.ecal.e");
  TTreeReaderValue<Double_t> ecalX(reader, "earm.ecal.x");
  TTreeReaderValue<Double_t> ecalY(reader, "earm.ecal.y");

  Long64_t eventCount = 0;
  Long64_t admittedEventCount = 0;
  Long64_t hardGatePairCount = 0;
  Long64_t signalPairCount = 0;
  Long64_t sidebandPairCount = 0;
  while (reader.Next()) {
    ++eventCount;
    if (!std::isfinite(*ecalEnergy) ||
        *ecalEnergy < config.ecalEnergyMinGeV ||
        *ecalEnergy > config.ecalEnergyMaxGeV ||
        !std::isfinite(*ecalTime) || *ecalTime < config.ecalTimeMinNs ||
        *ecalTime > config.ecalTimeMaxNs || !std::isfinite(*ecalX))
      continue;
    ++admittedEventCount;

    const size_t nPulses = std::min(
        {pixel.GetSize(), le.GetSize(), ecalResidual.GetSize(),
         calibValid.GetSize(), broadQuality.GetSize(), ecalEligible.GetSize(),
         spatialPass.GetSize(), correctedX.GetSize(), pulseY.GetSize(),
         pulseZ.GetSize()});
    std::vector<size_t> layer1;
    std::vector<size_t> layer2;
    for (size_t i = 0; i < nPulses; ++i) {
      if (!(calibValid[i] > 0.5 && broadQuality[i] > 0.5 &&
            ecalEligible[i] > 0.5 && spatialPass[i] > 0.5) ||
          !std::isfinite(pixel[i]) || !std::isfinite(le[i]) ||
          !std::isfinite(ecalResidual[i]) ||
          !std::isfinite(correctedX[i]) || !std::isfinite(pulseY[i]) ||
          !std::isfinite(pulseZ[i]))
        continue;
      const int pixelID = std::lround(pixel[i]);
      if (pixelID < 0 || pixelID >= kCDetNPixels)
        continue;
      (pixelID < kCDetNPixels / 2 ? layer1 : layer2).push_back(i);
    }

    for (size_t i1 : layer1) {
      for (size_t i2 : layer2) {
        const double dt = le[i2] - le[i1];
        const double dx = correctedX[i2] - correctedX[i1];
        const double dy = pulseY[i2] - pulseY[i1];
        const double meanZ = 0.5 * (pulseZ[i1] + pulseZ[i2]);
        const double alignedProjectedY =
            *ecalY * meanZ / kECalZFromTargetM + 0.10;
        const bool sameSide = std::fabs(dy) <= kPairDeltaYMaxM;
        const bool oppositeSide = config.oppositeSideEnabled &&
            std::fabs(std::fabs(dy) - config.oppositeDYCenterM) <=
                config.oppositeDYToleranceM &&
            std::fabs(alignedProjectedY -
                      config.oppositeProjectedYCenterM) <=
                config.oppositeProjectedYMaxM;
        if (std::fabs(dt) > kPairDeltaTimeMaxNs ||
            std::fabs(dx) > kPairDeltaXMaxM ||
            !(sameSide || oppositeSide))
          continue;

        const double residual =
            dx - (*ecalX / kECalZFromTargetM) * (pulseZ[i2] - pulseZ[i1]);
        const double timing =
            0.5 * (ecalResidual[i1] + ecalResidual[i2]);
        if (!std::isfinite(residual) || !std::isfinite(timing))
          continue;
        ++hardGatePairCount;
        hTimingVsResidual.Fill(residual, timing);
        const double timingOffset =
            std::fabs(timing - config.pairTimingCenterNs);
        if (timingOffset <= signalHalfWidthNs) {
          hSignal.Fill(residual);
          ++signalPairCount;
        } else if (timingOffset >= sidebandInnerNs &&
                   timingOffset <= sidebandOuterNs) {
          hSideband.Fill(residual);
          ++sidebandPairCount;
        }
      }
    }
  }

  if (hSignal.GetEntries() < 20 || hSideband.GetEntries() < 20) {
    std::cerr << "[CDet x resolution] Insufficient signal or sideband pairs "
              << "for fitting.\n";
    return;
  }

  const double signalTimingWidth = 2.0 * signalHalfWidthNs;
  const double sidebandTimingWidth =
      2.0 * (sidebandOuterNs - sidebandInnerNs);
  const double timingWidthScale = signalTimingWidth / sidebandTimingWidth;

  gCDetXResolutionSidebandTemplate = &hSideband;
  TF1 fit("fCDetXResolutionSignalPlusSideband", SignalPlusSidebandTemplate,
          fitMinM, fitMaxM, 4);
  fit.SetParNames("Gaussian yield", "Gaussian mean", "Gaussian sigma",
                  "Sideband scale");
  const double approximateBackground =
      timingWidthScale * hSideband.Integral(
          hSideband.FindBin(fitMinM), hSideband.FindBin(fitMaxM));
  const double approximateSignal = std::max(
      1.0, hSignal.Integral(hSignal.FindBin(fitMinM),
                            hSignal.FindBin(fitMaxM)) - approximateBackground);
  fit.SetParameters(approximateSignal, 0.0, 0.010, timingWidthScale);
  fit.SetParLimits(0, 0.0, 10.0 * hSignal.GetEntries());
  fit.SetParLimits(1, -0.03, 0.03);
  fit.SetParLimits(2, 0.0005, 0.05);
  fit.SetParLimits(3, 0.0, 10.0 * timingWidthScale);
  const int fitStatus = hSignal.Fit(&fit, "QRS");

  TF1 background("fCDetXResolutionSidebandComponent", SidebandTemplateOnly,
                 fitMinM, fitMaxM, 1);
  background.SetParameter(0, fit.GetParameter(3));
  background.SetLineColor(kBlue + 1);
  background.SetLineStyle(2);
  background.SetLineWidth(2);

  TF1 gaussian("fCDetXResolutionGaussianComponent", GaussianOnly, fitMinM,
               fitMaxM, 3);
  gaussian.SetParameters(fit.GetParameter(0), fit.GetParameter(1),
                         fit.GetParameter(2));
  gaussian.SetLineColor(kGreen + 2);
  gaussian.SetLineStyle(2);
  gaussian.SetLineWidth(2);

  TH1D hSubtracted(hSignal);
  hSubtracted.SetName("hCDetXResolutionSidebandSubtracted");
  hSubtracted.SetTitle(
      "Signal region minus width-scaled sidebands;r_{x} (m);Pairs");
  hSubtracted.Add(&hSideband, -timingWidthScale);
  TH1D hSidebandScaled(hSideband);
  hSidebandScaled.SetName("hCDetXResolutionSidebandsWidthScaled");
  hSidebandScaled.SetTitle(
      "Timing sidebands scaled to signal-window width;r_{x} (m);Pairs");
  hSidebandScaled.Scale(timingWidthScale);

  gSystem->mkdir(outputDirectory, true);
  gStyle->SetOptFit(1111);
  TCanvas canvas("cCDetXResolutionStudy", "CDet x-resolution study", 1500,
                 1100);
  canvas.Divide(2, 2);
  canvas.cd(1);
  gPad->SetRightMargin(0.14);
  hTimingVsResidual.SetStats(false);
  hTimingVsResidual.Draw("COLZ");
  const double timingLines[] = {
      config.pairTimingCenterNs - sidebandOuterNs,
      config.pairTimingCenterNs - sidebandInnerNs,
      config.pairTimingCenterNs - signalHalfWidthNs,
      config.pairTimingCenterNs + signalHalfWidthNs,
      config.pairTimingCenterNs + sidebandInnerNs,
      config.pairTimingCenterNs + sidebandOuterNs};
  for (int i = 0; i < 6; ++i) {
    TLine line(residualMinM, timingLines[i], residualMaxM, timingLines[i]);
    line.SetLineColor((i == 2 || i == 3) ? kRed + 1 : kBlue + 1);
    line.SetLineStyle(2);
    line.DrawClone("same");
  }

  canvas.cd(2);
  hSignal.SetLineColor(kBlack);
  hSignal.SetLineWidth(2);
  hSignal.Draw("E");
  hSidebandScaled.SetLineColor(kBlue + 1);
  hSidebandScaled.SetLineWidth(2);
  hSidebandScaled.Draw("HIST SAME");

  canvas.cd(3);
  hSignal.Draw("E");
  fit.SetLineColor(kRed + 1);
  fit.SetLineWidth(3);
  fit.Draw("same");
  background.Draw("same");
  gaussian.Draw("same");

  canvas.cd(4);
  hSubtracted.SetLineColor(kBlack);
  hSubtracted.SetLineWidth(2);
  hSubtracted.Draw("E");
  gaussian.Draw("same");
  canvas.SaveAs(Form("%s/CDetXResolutionStudy.pdf", outputDirectory));
  canvas.SaveAs(Form("%s/CDetXResolutionStudy.png", outputDirectory));

  TFile output(Form("%s/CDetXResolutionStudy.root", outputDirectory),
               "RECREATE");
  hTimingVsResidual.Write();
  hSignal.Write();
  hSideband.Write();
  hSidebandScaled.Write();
  hSubtracted.Write();
  fit.Write();
  background.Write();
  gaussian.Write();
  output.Close();

  std::ofstream summary(
      Form("%s/CDetXResolutionSummary.txt", outputDirectory));
  summary << std::setprecision(10)
          << "run = " << config.runNumber << '\n'
          << "files = " << filesAdded << '\n'
          << "events = " << eventCount << '\n'
          << "ecal_admitted_events = " << admittedEventCount << '\n'
          << "hard_gate_pairs = " << hardGatePairCount << '\n'
          << "signal_pairs = " << signalPairCount << '\n'
          << "sideband_pairs = " << sidebandPairCount << '\n'
          << "timing_center_ns = " << config.pairTimingCenterNs << '\n'
          << "signal_half_width_ns = " << signalHalfWidthNs << '\n'
          << "sideband_inner_ns = " << sidebandInnerNs << '\n'
          << "sideband_outer_ns = " << sidebandOuterNs << '\n'
          << "timing_width_scale = " << timingWidthScale << '\n'
          << "fit_status = " << fitStatus << '\n'
          << "gaussian_yield = " << fit.GetParameter(0) << '\n'
          << "gaussian_yield_error = " << fit.GetParError(0) << '\n'
          << "residual_mean_m = " << fit.GetParameter(1) << '\n'
          << "residual_mean_error_m = " << fit.GetParError(1) << '\n'
          << "residual_sigma_m = " << fit.GetParameter(2) << '\n'
          << "residual_sigma_error_m = " << fit.GetParError(2) << '\n'
          << "single_layer_equal_resolution_m = "
          << fit.GetParameter(2) / std::sqrt(2.0) << '\n'
          << "sideband_fit_scale = " << fit.GetParameter(3) << '\n'
          << "sideband_fit_scale_error = " << fit.GetParError(3) << '\n'
          << "chi2 = " << fit.GetChisquare() << '\n'
          << "ndf = " << fit.GetNDF() << '\n';

  std::cout << "[CDet x resolution] Files/events/ECal-admitted: "
            << filesAdded << '/' << eventCount << '/' << admittedEventCount
            << '\n'
            << "[CDet x resolution] Hard-gate/signal/sideband pairs: "
            << hardGatePairCount << '/' << signalPairCount << '/'
            << sidebandPairCount << '\n'
            << "[CDet x resolution] Fit status: " << fitStatus
            << "; residual mean = " << 1000.0 * fit.GetParameter(1)
            << " +/- " << 1000.0 * fit.GetParError(1)
            << " mm; sigma = " << 1000.0 * fit.GetParameter(2) << " +/- "
            << 1000.0 * fit.GetParError(2) << " mm\n"
            << "[CDet x resolution] Equal independent single-layer estimate: "
            << 1000.0 * fit.GetParameter(2) / std::sqrt(2.0) << " mm\n"
            << "[CDet x resolution] Sideband scale fitted/width expected: "
            << fit.GetParameter(3) << '/' << timingWidthScale << '\n'
            << "[CDet x resolution] Output directory: " << outputDirectory
            << std::endl;
}
