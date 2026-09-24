#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#include "CDetGoodPulseConfig.h"
#include "CDetRunDataset.h"

namespace {

constexpr int kCDetNPixelsDY = 2688;
constexpr double kECalZFromTargetMDY = 6.144;
constexpr double kPairDeltaTimeMaxNsDY = 15.0;
constexpr double kPairDeltaXMaxMDY = 0.15;
constexpr double kNominalDeltaYMaxMDY = 0.08;
constexpr double kSelectionYResidualOffsetMDY = 0.10;

bool HasDYStudyBranches(TChain &chain) {
  const char *required[] = {
      "earm.cdet.pulse.pmtnum", "earm.cdet.pulse.tdc_le_corr",
      "earm.cdet.pulse.ecal_residual", "earm.cdet.pulse.calib_valid",
      "earm.cdet.pulse.broad_quality_pass",
      "earm.cdet.pulse.ecal_eligible", "earm.cdet.pulse.spatial_pass",
      "earm.cdet.pulse.x_corr", "earm.cdet.pulse.y",
      "earm.cdet.pulse.z", "earm.ecal.e", "earm.ecal.adctime",
      "earm.ecal.x", "earm.ecal.y"};
  chain.LoadTree(0);
  for (const char *name : required) {
    if (!chain.GetBranch(name)) {
      std::cerr << "[CDet dy study] Missing branch: " << name << '\n';
      return false;
    }
  }
  return true;
}

struct DYCandidate {
  double absoluteDY = std::numeric_limits<double>::infinity();
  double dy = 0.0;
  double trajectoryDY = 0.0;
  double timingResidual = 0.0;
  double xResidual = 0.0;
  double yResidualL1 = 0.0;
  double yResidualL2 = 0.0;
  double projectedYAtPair = 0.0;
  double ellipseScore = std::numeric_limits<double>::infinity();
};

} // namespace

void Plot_CDet_DYGateStudy(
    const char *configFile = "CDet_run6077_projection.conf",
    const char *inputDirectory = nullptr,
    const char *outputDirectory = "CDet_run6077_dy_gate_study",
    double maxDYMaxM = 0.40, double stepM = 0.01) {
  CDetGoodPulseConfig::Values config;
  if (!CDetGoodPulseConfig::Load(configFile, config, "CDet dy study"))
    return;
  if (maxDYMaxM < kNominalDeltaYMaxMDY || stepM <= 0.0) {
    std::cerr << "[CDet dy study] Invalid scan range.\n";
    return;
  }

  TString input(inputDirectory ? inputDirectory : "");
  if (input.IsNull()) {
    const char *outDir = gSystem->Getenv("OUT_DIR");
    if (outDir)
      input = outDir;
  }
  if (input.IsNull()) {
    std::cerr << "[CDet dy study] An input directory or OUT_DIR is required.\n";
    return;
  }

  TChain chain("T");
  const int filesAdded =
      CDetRunDataset::AddToChain(&chain, config.runNumber, input.Data());
  if (filesAdded <= 0 || chain.GetEntries() <= 0 ||
      !HasDYStudyBranches(chain))
    return;

  std::vector<double> thresholds;
  for (double value = kNominalDeltaYMaxMDY;
       value <= maxDYMaxM + 0.5 * stepM;
       value += stepM)
    thresholds.push_back(value);
  std::vector<Long64_t> eventCounts(thresholds.size(), 0);

  TH1D hMinimumAbsDY(
      "hCDetMinimumAbsDY",
      "Best available pair in each ECal-admitted event;minimum |y_{L2}-y_{L1}| (m);Events",
      std::max(1, int(std::ceil(maxDYMaxM / 0.005))), 0.0, maxDYMaxM);
  TH1D hRecoveredDY(
      "hCDetDYRecoveredPairDY",
      "Candidates newly admitted beyond nominal 8 cm gate;y_{L2}-y_{L1} (m);Events",
      160, -maxDYMaxM, maxDYMaxM);
  TH1D hRecoveredTrajectoryDY(
      "hCDetDYRecoveredTrajectoryResidual",
      "Newly admitted candidates;#Delta y_{CDet}-(y_{ECal}/z_{ECal})#Delta z (m);Events",
      160, -maxDYMaxM, maxDYMaxM);
  TH1D hRecoveredTiming(
      "hCDetDYRecoveredTimingResidual",
      "Newly admitted candidates;t_{ECal}-<t_{CDet,corr}> (ns);Events",
      120, -60.0, 0.0);
  TH1D hRecoveredXResidual(
      "hCDetDYRecoveredXResidual",
      "Newly admitted candidates;<x_{CDet,corr}>-x_{ECal projected} (m);Events",
      160, -0.20, 0.20);
  TH2D hRecoveredLayerYResiduals(
      "hCDetDYRecoveredLayerYResiduals",
      "Newly admitted candidates;Layer 1 y residual (m);Layer 2 y residual (m)",
      160, -0.45, 0.45, 160, -0.45, 0.45);
  TH1D hRecoveredProjectedY(
      "hCDetDYRecoveredProjectedECalY",
      "Newly admitted candidates;ECal y projected to pair mean z (m);Events",
      160, -0.8, 0.8);

  TTreeReader reader(&chain);
  TTreeReaderArray<Double_t> pixel(reader, "earm.cdet.pulse.pmtnum");
  TTreeReaderArray<Double_t> le(reader, "earm.cdet.pulse.tdc_le_corr");
  TTreeReaderArray<Double_t> ecalResidual(
      reader, "earm.cdet.pulse.ecal_residual");
  TTreeReaderArray<Double_t> calibValid(
      reader, "earm.cdet.pulse.calib_valid");
  TTreeReaderArray<Double_t> broadQuality(
      reader, "earm.cdet.pulse.broad_quality_pass");
  TTreeReaderArray<Double_t> ecalEligible(
      reader, "earm.cdet.pulse.ecal_eligible");
  TTreeReaderArray<Double_t> spatialPass(
      reader, "earm.cdet.pulse.spatial_pass");
  TTreeReaderArray<Double_t> correctedX(reader, "earm.cdet.pulse.x_corr");
  TTreeReaderArray<Double_t> pulseY(reader, "earm.cdet.pulse.y");
  TTreeReaderArray<Double_t> pulseZ(reader, "earm.cdet.pulse.z");
  TTreeReaderValue<Double_t> ecalEnergy(reader, "earm.ecal.e");
  TTreeReaderValue<Double_t> ecalTime(reader, "earm.ecal.adctime");
  TTreeReaderValue<Double_t> ecalX(reader, "earm.ecal.x");
  TTreeReaderValue<Double_t> ecalY(reader, "earm.ecal.y");

  Long64_t admittedEvents = 0;
  Long64_t eventsWithDTDXEllipseCandidate = 0;
  while (reader.Next()) {
    if (!std::isfinite(*ecalEnergy) || !std::isfinite(*ecalTime) ||
        !std::isfinite(*ecalX) || !std::isfinite(*ecalY) ||
        *ecalEnergy < config.ecalEnergyMinGeV ||
        *ecalEnergy > config.ecalEnergyMaxGeV ||
        *ecalTime < config.ecalTimeMinNs ||
        *ecalTime > config.ecalTimeMaxNs)
      continue;
    ++admittedEvents;

    const size_t n = std::min(
        {pixel.GetSize(), le.GetSize(), ecalResidual.GetSize(),
         calibValid.GetSize(), broadQuality.GetSize(), ecalEligible.GetSize(),
         spatialPass.GetSize(), correctedX.GetSize(), pulseY.GetSize(),
         pulseZ.GetSize()});
    std::vector<size_t> layer1;
    std::vector<size_t> layer2;
    for (size_t i = 0; i < n; ++i) {
      if (!(calibValid[i] > 0.5 && broadQuality[i] > 0.5 &&
            ecalEligible[i] > 0.5 && spatialPass[i] > 0.5) ||
          !std::isfinite(pixel[i]) || !std::isfinite(le[i]) ||
          !std::isfinite(ecalResidual[i]) || !std::isfinite(correctedX[i]) ||
          !std::isfinite(pulseY[i]) || !std::isfinite(pulseZ[i]))
        continue;
      const int id = int(std::lround(pixel[i]));
      if (id < 0 || id >= kCDetNPixelsDY)
        continue;
      (id < kCDetNPixelsDY / 2 ? layer1 : layer2).push_back(i);
    }

    DYCandidate best;
    bool found = false;
    for (size_t i1 : layer1) {
      for (size_t i2 : layer2) {
        const double dt = le[i2] - le[i1];
        const double dx = correctedX[i2] - correctedX[i1];
        if (std::fabs(dt) > kPairDeltaTimeMaxNsDY ||
            std::fabs(dx) > kPairDeltaXMaxMDY)
          continue;
        const double deltaZ = pulseZ[i2] - pulseZ[i1];
        const double trajectoryX =
            dx - (*ecalX / kECalZFromTargetMDY) * deltaZ;
        const double timing =
            0.5 * (ecalResidual[i1] + ecalResidual[i2]);
        const double xPull =
            (trajectoryX - config.pairResidualCenterM) /
            config.pairResidualScaleM;
        const double timingPull =
            (timing - config.pairTimingCenterNs) /
            config.pairTimingScaleNs;
        const double score = xPull * xPull + timingPull * timingPull;
        if (score > config.pairCutRadius * config.pairCutRadius)
          continue;

        const double dy = pulseY[i2] - pulseY[i1];
        const double absoluteDY = std::fabs(dy);
        if (!found || absoluteDY < best.absoluteDY ||
            (absoluteDY == best.absoluteDY && score < best.ellipseScore)) {
          found = true;
          best.absoluteDY = absoluteDY;
          best.dy = dy;
          best.trajectoryDY =
              dy - (*ecalY / kECalZFromTargetMDY) * deltaZ;
          best.timingResidual = timing;
          const double meanZ = 0.5 * (pulseZ[i1] + pulseZ[i2]);
          const double meanX = 0.5 * (correctedX[i1] + correctedX[i2]);
          best.xResidual = meanX - *ecalX * meanZ / kECalZFromTargetMDY;
          best.yResidualL1 = pulseY[i1] -
              *ecalY * pulseZ[i1] / kECalZFromTargetMDY -
              kSelectionYResidualOffsetMDY;
          best.yResidualL2 = pulseY[i2] -
              *ecalY * pulseZ[i2] / kECalZFromTargetMDY -
              kSelectionYResidualOffsetMDY;
          best.projectedYAtPair =
              *ecalY * meanZ / kECalZFromTargetMDY;
          best.ellipseScore = score;
        }
      }
    }
    if (!found)
      continue;

    ++eventsWithDTDXEllipseCandidate;
    hMinimumAbsDY.Fill(best.absoluteDY);
    const size_t firstPassing = std::lower_bound(
        thresholds.begin(), thresholds.end(), best.absoluteDY) -
        thresholds.begin();
    for (size_t i = firstPassing; i < thresholds.size(); ++i)
      ++eventCounts[i];

    if (best.absoluteDY > kNominalDeltaYMaxMDY &&
        best.absoluteDY <= maxDYMaxM) {
      hRecoveredDY.Fill(best.dy);
      hRecoveredTrajectoryDY.Fill(best.trajectoryDY);
      hRecoveredTiming.Fill(best.timingResidual);
      hRecoveredXResidual.Fill(best.xResidual);
      hRecoveredLayerYResiduals.Fill(best.yResidualL1, best.yResidualL2);
      hRecoveredProjectedY.Fill(best.projectedYAtPair);
    }
  }

  gSystem->mkdir(outputDirectory, true);
  const auto nominalIterator = std::lower_bound(
      thresholds.begin(), thresholds.end(), kNominalDeltaYMaxMDY);
  const size_t nominalIndex = std::min<size_t>(
      nominalIterator - thresholds.begin(), thresholds.size() - 1);
  const Long64_t nominalEvents = eventCounts[nominalIndex];

  TGraph totalGraph(thresholds.size());
  TGraph addedGraph(thresholds.size());
  for (size_t i = 0; i < thresholds.size(); ++i) {
    totalGraph.SetPoint(i, thresholds[i], eventCounts[i]);
    addedGraph.SetPoint(i, thresholds[i], eventCounts[i] - nominalEvents);
  }
  totalGraph.SetName("gCDetPairEventsVsDYMax");
  totalGraph.SetTitle(
      "Pair-event yield versus inter-layer y gate;|#Delta y| maximum (m);Events with selected pair");
  addedGraph.SetName("gCDetAdditionalPairEventsVsDYMax");
  addedGraph.SetTitle(
      "Additional pair events beyond nominal 8 cm gate;|#Delta y| maximum (m);Additional events");
  for (TGraph *graph : {&totalGraph, &addedGraph}) {
    graph->SetMarkerStyle(20);
    graph->SetLineWidth(2);
  }

  TCanvas scanCanvas("cCDetDYGateScan", "CDet dy gate scan", 1800, 550);
  scanCanvas.Divide(3, 1);
  scanCanvas.cd(1);
  hMinimumAbsDY.Draw("HIST");
  scanCanvas.cd(2);
  totalGraph.Draw("APL");
  scanCanvas.cd(3);
  addedGraph.Draw("APL");
  scanCanvas.SaveAs(Form("%s/CDetDYGateScan.pdf", outputDirectory));
  scanCanvas.SaveAs(Form("%s/CDetDYGateScan.png", outputDirectory));

  TCanvas qualityCanvas(
      "cCDetDYRecoveredQuality", "Quality of dy-gate recoveries", 1800, 900);
  qualityCanvas.Divide(3, 2);
  qualityCanvas.cd(1); hRecoveredDY.Draw("HIST");
  qualityCanvas.cd(2); hRecoveredTrajectoryDY.Draw("HIST");
  qualityCanvas.cd(3); hRecoveredTiming.Draw("HIST");
  qualityCanvas.cd(4); hRecoveredXResidual.Draw("HIST");
  qualityCanvas.cd(5); hRecoveredLayerYResiduals.Draw("COLZ");
  qualityCanvas.cd(6); hRecoveredProjectedY.Draw("HIST");
  qualityCanvas.SaveAs(Form("%s/CDetDYRecoveredQuality.pdf", outputDirectory));
  qualityCanvas.SaveAs(Form("%s/CDetDYRecoveredQuality.png", outputDirectory));

  std::ofstream table(Form("%s/CDetDYGateScan.csv", outputDirectory));
  table << "dy_max_m,selected_pair_events,additional_events,conditional_fraction\n";
  for (size_t i = 0; i < thresholds.size(); ++i)
    table << thresholds[i] << ',' << eventCounts[i] << ','
          << eventCounts[i] - nominalEvents << ','
          << (admittedEvents > 0
                  ? static_cast<double>(eventCounts[i]) / admittedEvents
                  : 0.0)
          << '\n';

  TFile output(Form("%s/CDetDYGateStudy.root", outputDirectory), "RECREATE");
  hMinimumAbsDY.Write();
  hRecoveredDY.Write();
  hRecoveredTrajectoryDY.Write();
  hRecoveredTiming.Write();
  hRecoveredXResidual.Write();
  hRecoveredLayerYResiduals.Write();
  hRecoveredProjectedY.Write();
  totalGraph.Write();
  addedGraph.Write();
  output.Close();

  std::cout << "[CDet dy study] Files/events: " << filesAdded << '/'
            << chain.GetEntries() << "\n"
            << "[CDet dy study] ECal-admitted events: " << admittedEvents
            << "\n[CDet dy study] Events with an ellipse + dt/dx candidate: "
            << eventsWithDTDXEllipseCandidate
            << "\n[CDet dy study] Events passing nominal |dy| <= "
            << thresholds[nominalIndex] << " m: " << nominalEvents
            << "\n[CDet dy study] Events passing |dy| <= "
            << thresholds.back() << " m: " << eventCounts.back()
            << " (additional " << eventCounts.back() - nominalEvents << ")\n"
            << "[CDet dy study] Output directory: " << outputDirectory
            << std::endl;
}
