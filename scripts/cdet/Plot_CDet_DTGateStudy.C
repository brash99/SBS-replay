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

constexpr int kCDetNPixelsDT = 2688;
constexpr double kECalZFromTargetMDT = 6.144;
constexpr double kPairDeltaXMaxMDT = 0.15;
constexpr double kNominalDeltaTimeMaxNsDT = 15.0;
constexpr double kSelectionYResidualOffsetMDT = 0.10;

bool HasDTStudyBranches(TChain &chain) {
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
      std::cerr << "[CDet dt study] Missing branch: " << name << '\n';
      return false;
    }
  }
  return true;
}

struct DTCandidate {
  size_t indexL1 = 0;
  size_t indexL2 = 0;
  double absoluteDT = std::numeric_limits<double>::infinity();
  double dt = 0.0;
  double timingResidual = 0.0;
  double xResidual = 0.0;
  double trajectoryResidual = 0.0;
  double dy = 0.0;
  double alignedProjectedY = 0.0;
  double ellipseScore = std::numeric_limits<double>::infinity();
  int yTopology = -1;
};

bool CandidateOrder(const DTCandidate &left, const DTCandidate &right) {
  if (left.ellipseScore != right.ellipseScore)
    return left.ellipseScore < right.ellipseScore;
  if (left.indexL1 != right.indexL1)
    return left.indexL1 < right.indexL1;
  return left.indexL2 < right.indexL2;
}

std::vector<size_t> SelectOneToOne(const std::vector<DTCandidate> &candidates,
                                   double maximumAbsoluteDT,
                                   size_t pulseCount) {
  std::vector<bool> used(pulseCount, false);
  std::vector<size_t> selected;
  for (size_t i = 0; i < candidates.size(); ++i) {
    const DTCandidate &candidate = candidates[i];
    if (candidate.absoluteDT > maximumAbsoluteDT ||
        used[candidate.indexL1] || used[candidate.indexL2])
      continue;
    used[candidate.indexL1] = true;
    used[candidate.indexL2] = true;
    selected.push_back(i);
  }
  return selected;
}

} // namespace

void Plot_CDet_DTGateStudy(
    const char *configFile = "CDet_run6077_projection.conf",
    const char *inputDirectory = nullptr,
    const char *outputDirectory = "CDet_run6077_dt_gate_study",
    double minDTMaxNs = 5.0, double maxDTMaxNs = 30.0,
    double stepNs = 0.5) {
  CDetGoodPulseConfig::Values config;
  if (!CDetGoodPulseConfig::Load(configFile, config, "CDet dt study"))
    return;
  if (minDTMaxNs < 0.0 || maxDTMaxNs < kNominalDeltaTimeMaxNsDT ||
      maxDTMaxNs < minDTMaxNs || stepNs <= 0.0) {
    std::cerr << "[CDet dt study] Invalid scan range.\n";
    return;
  }

  TString input(inputDirectory ? inputDirectory : "");
  if (input.IsNull()) {
    const char *outDir = gSystem->Getenv("OUT_DIR");
    if (outDir)
      input = outDir;
  }
  if (input.IsNull()) {
    std::cerr << "[CDet dt study] An input directory or OUT_DIR is required.\n";
    return;
  }

  TChain chain("T");
  const int filesAdded =
      CDetRunDataset::AddToChain(&chain, config.runNumber, input.Data());
  if (filesAdded <= 0 || chain.GetEntries() <= 0 ||
      !HasDTStudyBranches(chain))
    return;

  std::vector<double> thresholds;
  for (double value = minDTMaxNs; value <= maxDTMaxNs + 0.5 * stepNs;
       value += stepNs)
    thresholds.push_back(value);
  if (std::none_of(thresholds.begin(), thresholds.end(), [](double value) {
        return std::fabs(value - kNominalDeltaTimeMaxNsDT) < 1e-9;
      })) {
    thresholds.push_back(kNominalDeltaTimeMaxNsDT);
    std::sort(thresholds.begin(), thresholds.end());
  }

  std::vector<Long64_t> eventCounts(thresholds.size(), 0);
  std::vector<Long64_t> pairCounts(thresholds.size(), 0);
  std::vector<double> recoveredTimingSum(thresholds.size(), 0.0);
  std::vector<double> recoveredTimingSum2(thresholds.size(), 0.0);
  std::vector<double> recoveredXSum(thresholds.size(), 0.0);
  std::vector<double> recoveredXSum2(thresholds.size(), 0.0);
  TH1D hMinimumAbsDT(
      "hCDetMinimumAbsDT",
      "Best available ellipse-compatible pair;minimum |t_{L2}-t_{L1}| (ns);Events",
      std::max(1, int(std::ceil(maxDTMaxNs / 0.25))), 0.0, maxDTMaxNs);
  TH1D hRecoveredDT(
      "hCDetDTRecoveredPairDT",
      "Events newly selected beyond nominal 15 ns gate;t_{L2}-t_{L1} (ns);Events",
      160, -maxDTMaxNs, maxDTMaxNs);
  TH1D hRecoveredTiming(
      "hCDetDTRecoveredTimingResidual",
      "Newly selected events;t_{ECal}-<t_{CDet,corr}> (ns);Events",
      120, -60.0, 0.0);
  TH1D hRecoveredXResidual(
      "hCDetDTRecoveredXResidual",
      "Newly selected events;<x_{CDet,corr}>-x_{ECal projected} (m);Events",
      160, -0.20, 0.20);
  TH1D hRecoveredTrajectoryResidual(
      "hCDetDTRecoveredTrajectoryResidual",
      "Newly selected events;inter-layer trajectory residual (m);Events",
      160, -0.08, 0.08);
  TH2D hRecoveredDTvsTiming(
      "hCDetDTRecoveredDTvsTiming",
      "Newly selected events; t_{L2}-t_{L1} (ns);t_{ECal}-<t_{CDet,corr}> (ns)",
      160, -maxDTMaxNs, maxDTMaxNs, 120, -60.0, 0.0);
  TH2D hRecoveredDTvsX(
      "hCDetDTRecoveredDTvsX",
      "Newly selected events;t_{L2}-t_{L1} (ns);projected-x residual (m)",
      160, -maxDTMaxNs, maxDTMaxNs, 160, -0.20, 0.20);
  TH1D hRecoveredYTopology(
      "hCDetDTRecoveredYTopology",
      "Newly selected events;y topology (0=same side, 1=opposite seam);Events",
      2, -0.5, 1.5);

  TTreeReader reader(&chain);
  TTreeReaderArray<Double_t> pixel(reader, "earm.cdet.pulse.pmtnum");
  TTreeReaderArray<Double_t> le(reader, "earm.cdet.pulse.tdc_le_corr");
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
  TTreeReaderArray<Double_t> pulseZ(reader, "earm.cdet.pulse.z");
  TTreeReaderValue<Double_t> ecalEnergy(reader, "earm.ecal.e");
  TTreeReaderValue<Double_t> ecalTime(reader, "earm.ecal.adctime");
  TTreeReaderValue<Double_t> ecalX(reader, "earm.ecal.x");
  TTreeReaderValue<Double_t> ecalY(reader, "earm.ecal.y");

  Long64_t admittedEvents = 0;
  Long64_t eventsWithEllipseDXYCandidate = 0;
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
      if (id < 0 || id >= kCDetNPixelsDT)
        continue;
      (id < kCDetNPixelsDT / 2 ? layer1 : layer2).push_back(i);
    }

    std::vector<DTCandidate> candidates;
    for (size_t i1 : layer1) {
      for (size_t i2 : layer2) {
        const double dt = le[i2] - le[i1];
        const double dx = correctedX[i2] - correctedX[i1];
        if (std::fabs(dx) > kPairDeltaXMaxMDT)
          continue;
        const double dy = pulseY[i2] - pulseY[i1];
        const double meanZ = 0.5 * (pulseZ[i1] + pulseZ[i2]);
        const double alignedProjectedY =
            *ecalY * meanZ / kECalZFromTargetMDT +
            kSelectionYResidualOffsetMDT;
        const bool sameSide = std::fabs(dy) <= 0.08;
        const bool oppositeSide = config.oppositeSideEnabled &&
            std::fabs(std::fabs(dy) - config.oppositeDYCenterM) <=
                config.oppositeDYToleranceM &&
            std::fabs(alignedProjectedY -
                      config.oppositeProjectedYCenterM) <=
                config.oppositeProjectedYMaxM;
        if (!sameSide && !oppositeSide)
          continue;

        const double deltaZ = pulseZ[i2] - pulseZ[i1];
        const double trajectoryResidual =
            dx - (*ecalX / kECalZFromTargetMDT) * deltaZ;
        const double timingResidual =
            0.5 * (ecalResidual[i1] + ecalResidual[i2]);
        const double residualPull =
            (trajectoryResidual - config.pairResidualCenterM) /
            config.pairResidualScaleM;
        const double timingPull =
            (timingResidual - config.pairTimingCenterNs) /
            config.pairTimingScaleNs;
        const double score = residualPull * residualPull +
                             timingPull * timingPull;
        if (score > config.pairCutRadius * config.pairCutRadius)
          continue;

        DTCandidate candidate;
        candidate.indexL1 = i1;
        candidate.indexL2 = i2;
        candidate.absoluteDT = std::fabs(dt);
        candidate.dt = dt;
        candidate.timingResidual = timingResidual;
        candidate.trajectoryResidual = trajectoryResidual;
        candidate.dy = dy;
        candidate.alignedProjectedY = alignedProjectedY;
        candidate.ellipseScore = score;
        candidate.yTopology = sameSide ? 0 : 1;
        const double meanX = 0.5 * (correctedX[i1] + correctedX[i2]);
        candidate.xResidual =
            meanX - *ecalX * meanZ / kECalZFromTargetMDT;
        candidates.push_back(candidate);
      }
    }
    if (candidates.empty())
      continue;
    ++eventsWithEllipseDXYCandidate;
    std::sort(candidates.begin(), candidates.end(), CandidateOrder);
    const auto bestDT = std::min_element(
        candidates.begin(), candidates.end(),
        [](const DTCandidate &left, const DTCandidate &right) {
          return left.absoluteDT < right.absoluteDT;
        });
    hMinimumAbsDT.Fill(bestDT->absoluteDT);

    std::vector<std::vector<size_t>> selectedByThreshold;
    selectedByThreshold.reserve(thresholds.size());
    for (size_t i = 0; i < thresholds.size(); ++i) {
      selectedByThreshold.push_back(
          SelectOneToOne(candidates, thresholds[i], n));
      if (!selectedByThreshold.back().empty())
        ++eventCounts[i];
      pairCounts[i] += selectedByThreshold.back().size();
    }

    const size_t nominalIndex = std::lower_bound(
        thresholds.begin(), thresholds.end(), kNominalDeltaTimeMaxNsDT) -
        thresholds.begin();
    const bool selectedNominal = !selectedByThreshold[nominalIndex].empty();
    if (!selectedNominal) {
      for (size_t i = nominalIndex + 1; i < thresholds.size(); ++i) {
        if (selectedByThreshold[i].empty())
          continue;
        const DTCandidate &best = candidates[selectedByThreshold[i].front()];
        recoveredTimingSum[i] += best.timingResidual;
        recoveredTimingSum2[i] += best.timingResidual * best.timingResidual;
        recoveredXSum[i] += best.xResidual;
        recoveredXSum2[i] += best.xResidual * best.xResidual;
      }
    }
    const std::vector<size_t> &selectedMaximum = selectedByThreshold.back();
    if (!selectedNominal && !selectedMaximum.empty()) {
      const DTCandidate &best = candidates[selectedMaximum.front()];
      hRecoveredDT.Fill(best.dt);
      hRecoveredTiming.Fill(best.timingResidual);
      hRecoveredXResidual.Fill(best.xResidual);
      hRecoveredTrajectoryResidual.Fill(best.trajectoryResidual);
      hRecoveredDTvsTiming.Fill(best.dt, best.timingResidual);
      hRecoveredDTvsX.Fill(best.dt, best.xResidual);
      hRecoveredYTopology.Fill(best.yTopology);
    }
  }

  gSystem->mkdir(outputDirectory, true);
  const size_t nominalIndex = std::lower_bound(
      thresholds.begin(), thresholds.end(), kNominalDeltaTimeMaxNsDT) -
      thresholds.begin();
  const Long64_t nominalEvents = eventCounts[nominalIndex];
  const Long64_t nominalPairs = pairCounts[nominalIndex];

  TGraph eventGraph(thresholds.size());
  TGraph addedEventGraph(thresholds.size());
  TGraph pairGraph(thresholds.size());
  for (size_t i = 0; i < thresholds.size(); ++i) {
    eventGraph.SetPoint(i, thresholds[i], eventCounts[i]);
    addedEventGraph.SetPoint(i, thresholds[i], eventCounts[i] - nominalEvents);
    pairGraph.SetPoint(i, thresholds[i], pairCounts[i]);
  }
  eventGraph.SetName("gCDetPairEventsVsDTMax");
  eventGraph.SetTitle(
      "Pair-event yield versus inter-layer timing gate;|#Delta t| maximum (ns);Events with selected pair");
  addedEventGraph.SetName("gCDetAdditionalPairEventsVsDTMax");
  addedEventGraph.SetTitle(
      "Additional pair events beyond nominal 15 ns gate;|#Delta t| maximum (ns);Additional events");
  pairGraph.SetName("gCDetPairsVsDTMax");
  pairGraph.SetTitle(
      "Selected-pair yield versus inter-layer timing gate;|#Delta t| maximum (ns);Selected pairs");
  for (TGraph *graph : {&eventGraph, &addedEventGraph, &pairGraph}) {
    graph->SetMarkerStyle(20);
    graph->SetLineWidth(2);
  }

  TCanvas scanCanvas("cCDetDTGateScan", "CDet dt gate scan", 1800, 900);
  scanCanvas.Divide(2, 2);
  scanCanvas.cd(1); hMinimumAbsDT.Draw("HIST");
  scanCanvas.cd(2); eventGraph.Draw("APL");
  scanCanvas.cd(3); addedEventGraph.Draw("APL");
  scanCanvas.cd(4); pairGraph.Draw("APL");
  scanCanvas.SaveAs(Form("%s/CDetDTGateScan.pdf", outputDirectory));
  scanCanvas.SaveAs(Form("%s/CDetDTGateScan.png", outputDirectory));

  TCanvas qualityCanvas(
      "cCDetDTRecoveredQuality", "Quality of dt-gate recoveries", 1800, 900);
  qualityCanvas.Divide(4, 2);
  qualityCanvas.cd(1); hRecoveredDT.Draw("HIST");
  qualityCanvas.cd(2); hRecoveredTiming.Draw("HIST");
  qualityCanvas.cd(3); hRecoveredXResidual.Draw("HIST");
  qualityCanvas.cd(4); hRecoveredTrajectoryResidual.Draw("HIST");
  qualityCanvas.cd(5); hRecoveredDTvsTiming.Draw("COLZ");
  qualityCanvas.cd(6); hRecoveredDTvsX.Draw("COLZ");
  qualityCanvas.cd(7); hRecoveredYTopology.Draw("HIST");
  qualityCanvas.SaveAs(Form("%s/CDetDTRecoveredQuality.pdf", outputDirectory));
  qualityCanvas.SaveAs(Form("%s/CDetDTRecoveredQuality.png", outputDirectory));

  std::ofstream table(Form("%s/CDetDTGateScan.csv", outputDirectory));
  table << "dt_max_ns,selected_pair_events,selected_pairs,additional_events,"
           "additional_pairs,conditional_fraction,pairs_per_selected_event,"
           "additional_pairs_per_additional_event,recovered_timing_mean_ns,"
           "recovered_timing_rms_ns,recovered_x_mean_m,recovered_x_rms_m\n";
  for (size_t i = 0; i < thresholds.size(); ++i) {
    const Long64_t additionalEvents = eventCounts[i] - nominalEvents;
    const Long64_t additionalPairs = pairCounts[i] - nominalPairs;
    const double timingMean = additionalEvents > 0
        ? recoveredTimingSum[i] / additionalEvents : 0.0;
    const double timingVariance = additionalEvents > 0
        ? recoveredTimingSum2[i] / additionalEvents - timingMean * timingMean
        : 0.0;
    const double xMean = additionalEvents > 0
        ? recoveredXSum[i] / additionalEvents : 0.0;
    const double xVariance = additionalEvents > 0
        ? recoveredXSum2[i] / additionalEvents - xMean * xMean : 0.0;
    table << thresholds[i] << ',' << eventCounts[i] << ',' << pairCounts[i]
          << ',' << additionalEvents << ',' << additionalPairs << ','
          << (admittedEvents > 0
                  ? static_cast<double>(eventCounts[i]) / admittedEvents
                  : 0.0)
          << ',' << (eventCounts[i] > 0
                  ? static_cast<double>(pairCounts[i]) / eventCounts[i] : 0.0)
          << ',' << (additionalEvents > 0
                  ? static_cast<double>(additionalPairs) / additionalEvents
                  : 0.0)
          << ',' << timingMean << ',' << std::sqrt(std::max(0.0, timingVariance))
          << ',' << xMean << ',' << std::sqrt(std::max(0.0, xVariance))
          << '\n';
  }

  TFile output(Form("%s/CDetDTGateStudy.root", outputDirectory), "RECREATE");
  hMinimumAbsDT.Write();
  hRecoveredDT.Write();
  hRecoveredTiming.Write();
  hRecoveredXResidual.Write();
  hRecoveredTrajectoryResidual.Write();
  hRecoveredDTvsTiming.Write();
  hRecoveredDTvsX.Write();
  hRecoveredYTopology.Write();
  eventGraph.Write();
  addedEventGraph.Write();
  pairGraph.Write();
  output.Close();

  std::cout << "[CDet dt study] Files/events: " << filesAdded << '/'
            << chain.GetEntries() << "\n"
            << "[CDet dt study] ECal-admitted events: " << admittedEvents
            << "\n[CDet dt study] Events with an ellipse + dx/y candidate: "
            << eventsWithEllipseDXYCandidate
            << "\n[CDet dt study] Nominal |dt| <= "
            << thresholds[nominalIndex] << " ns: " << nominalEvents
            << " events, " << nominalPairs << " pairs"
            << "\n[CDet dt study] Maximum |dt| <= " << thresholds.back()
            << " ns: " << eventCounts.back() << " events, "
            << pairCounts.back() << " pairs (additional "
            << eventCounts.back() - nominalEvents << " events)\n"
            << "[CDet dt study] Newly selected event quality at maximum gate: "
            << "timing mean/RMS=" << hRecoveredTiming.GetMean() << '/'
            << hRecoveredTiming.GetRMS() << " ns, x mean/RMS="
            << hRecoveredXResidual.GetMean() << '/'
            << hRecoveredXResidual.GetRMS() << " m\n"
            << "[CDet dt study] Output directory: " << outputDirectory
            << std::endl;
}
