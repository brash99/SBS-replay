#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TStyle.h>
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

#include "CDetRunDataset.h"

namespace {

struct PairCandidate {
  size_t pulse1;
  size_t pulse2;
  double score;
};

bool HasPairScanBranches(TChain &chain) {
  const char *required[] = {
      "earm.cdet.pulse.pmtnum", "earm.cdet.pulse.tdc_le_corr",
      "earm.cdet.pulse.tdc_tot_ns", "earm.cdet.pulse.calib_valid",
      "earm.cdet.pulse.broad_quality_pass", "earm.cdet.pulse.x_corr",
      "earm.cdet.pulse.y", "earm.cdet.pulse.z", "earm.ecal.e",
      "earm.ecal.adctime", "earm.ecal.x", "earm.ecal.y"};
  chain.LoadTree(0);
  for (const char *name : required) {
    if (!chain.GetBranch(name)) {
      std::cerr << "[CDet pair-window scan] Missing branch: " << name
                << std::endl;
      return false;
    }
  }
  return true;
}

} // namespace

// Scan symmetric ECal ADC-time windows about zero. Pairing is reconstructed
// from stored calibrated pulses so this study may extend beyond the ECal-time
// window that was active when the ROOT file was produced.
void Plot_CDet_PairYieldVsECalTimingWindow(
    int runNumber = 5711, const char *inputDirectory = nullptr,
    const char *outputDirectory = "CDet_pair_timing_window_scan",
    double minHalfWidthNs = 1.0, double maxHalfWidthNs = 30.0,
    double stepNs = 1.0, double ecalEnergyMinGeV = 3.0,
    double ecalEnergyMaxGeV = 4.5, double pairResidualCenterM = 0.0,
    double pairTimingCenterNs = -26.0, double pairResidualScaleM = 0.020,
    double pairTimingScaleNs = 5.0, double pairCutRadius = 2.0) {
  if (minHalfWidthNs <= 0.0 || maxHalfWidthNs < minHalfWidthNs ||
      stepNs <= 0.0 || ecalEnergyMaxGeV <= ecalEnergyMinGeV ||
      pairResidualScaleM <= 0.0 || pairTimingScaleNs <= 0.0 ||
      pairCutRadius <= 0.0) {
    std::cerr << "[CDet pair-window scan] Invalid scan parameters."
              << std::endl;
    return;
  }

  TString input(inputDirectory ? inputDirectory : "");
  if (input.IsNull()) {
    const char *outDir = gSystem->Getenv("OUT_DIR");
    if (outDir)
      input = outDir;
  }
  if (input.IsNull()) {
    std::cerr << "[CDet pair-window scan] An input directory or OUT_DIR is required."
              << std::endl;
    return;
  }

  TChain chain("T");
  const int filesAdded =
      CDetRunDataset::AddToChain(&chain, runNumber, input.Data());
  if (filesAdded <= 0 || chain.GetEntries() <= 0 ||
      !HasPairScanBranches(chain))
    return;

  std::vector<double> halfWidths;
  for (double width = minHalfWidthNs;
       width <= maxHalfWidthNs + 0.5 * stepNs; width += stepNs)
    halfWidths.push_back(width);
  std::vector<Long64_t> denominators(halfWidths.size(), 0);
  std::vector<Long64_t> eventsWithPairs(halfWidths.size(), 0);
  std::vector<Long64_t> selectedPairs(halfWidths.size(), 0);
  const double energyCenterGeV =
      0.5 * (ecalEnergyMinGeV + ecalEnergyMaxGeV);
  const double maxEnergyHalfWidthGeV =
      0.5 * (ecalEnergyMaxGeV - ecalEnergyMinGeV);
  const double energyStepGeV = 0.05;
  std::vector<double> energyHalfWidths;
  for (double width = energyStepGeV;
       width <= maxEnergyHalfWidthGeV + 0.5 * energyStepGeV;
       width += energyStepGeV)
    energyHalfWidths.push_back(width);
  std::vector<Long64_t> energyDenominators(energyHalfWidths.size(), 0);
  std::vector<Long64_t> energyEventsWithPairs(energyHalfWidths.size(), 0);
  std::vector<Long64_t> energySelectedPairs(energyHalfWidths.size(), 0);
  Long64_t energyOnlyReferenceEvents = 0;

  TTreeReader reader(&chain);
  TTreeReaderArray<Double_t> pixel(reader, "earm.cdet.pulse.pmtnum");
  TTreeReaderArray<Double_t> le(reader, "earm.cdet.pulse.tdc_le_corr");
  TTreeReaderArray<Double_t> tot(reader, "earm.cdet.pulse.tdc_tot_ns");
  TTreeReaderArray<Double_t> calibValid(reader,
                                        "earm.cdet.pulse.calib_valid");
  TTreeReaderArray<Double_t> broadQuality(
      reader, "earm.cdet.pulse.broad_quality_pass");
  TTreeReaderArray<Double_t> correctedX(reader, "earm.cdet.pulse.x_corr");
  TTreeReaderArray<Double_t> pulseY(reader, "earm.cdet.pulse.y");
  TTreeReaderArray<Double_t> pulseZ(reader, "earm.cdet.pulse.z");
  TTreeReaderValue<Double_t> ecalEnergy(reader, "earm.ecal.e");
  TTreeReaderValue<Double_t> ecalTime(reader, "earm.ecal.adctime");
  TTreeReaderValue<Double_t> ecalX(reader, "earm.ecal.x");
  TTreeReaderValue<Double_t> ecalY(reader, "earm.ecal.y");

  constexpr int kLayerBoundary = 1344;
  constexpr double kECalZFromTargetM = 6.144;
  Long64_t eventCount = 0;
  while (reader.Next()) {
    ++eventCount;
    if (!std::isfinite(*ecalEnergy) || !std::isfinite(*ecalTime) ||
        *ecalEnergy < ecalEnergyMinGeV || *ecalEnergy > ecalEnergyMaxGeV)
      continue;
    ++energyOnlyReferenceEvents;
    if (std::fabs(*ecalTime) > maxHalfWidthNs)
      continue;

    const size_t firstWindow = static_cast<size_t>(std::lower_bound(
        halfWidths.begin(), halfWidths.end(), std::fabs(*ecalTime)) -
        halfWidths.begin());
    for (size_t window = firstWindow; window < halfWidths.size(); ++window)
      ++denominators[window];
    const double energyDistance = std::fabs(*ecalEnergy - energyCenterGeV);
    const size_t firstEnergyWindow = static_cast<size_t>(std::lower_bound(
        energyHalfWidths.begin(), energyHalfWidths.end(), energyDistance) -
        energyHalfWidths.begin());
    for (size_t window = firstEnergyWindow;
         window < energyHalfWidths.size(); ++window)
      ++energyDenominators[window];

    const bool positionEligible =
        std::isfinite(*ecalX) && std::isfinite(*ecalY) &&
        *ecalX > -1.5 && *ecalX < 1.5 && *ecalY > -1.2 && *ecalY < 1.2 &&
        *ecalX != 0.0 && *ecalY != 0.0;
    if (!positionEligible)
      continue;

    const size_t nPulses = std::min(
        {pixel.GetSize(), le.GetSize(), tot.GetSize(), calibValid.GetSize(),
         broadQuality.GetSize(), correctedX.GetSize(), pulseY.GetSize(),
         pulseZ.GetSize()});
    std::vector<size_t> layer1;
    std::vector<size_t> layer2;
    for (size_t i = 0; i < nPulses; ++i) {
      if (!(calibValid[i] > 0.5 && broadQuality[i] > 0.5) ||
          !std::isfinite(pixel[i]) || !std::isfinite(le[i]) ||
          !std::isfinite(tot[i]) || !std::isfinite(correctedX[i]) ||
          !std::isfinite(pulseY[i]) || !std::isfinite(pulseZ[i]))
        continue;
      const double projectedX = *ecalX * pulseZ[i] / kECalZFromTargetM;
      const double projectedY = *ecalY * pulseZ[i] / kECalZFromTargetM;
      if (std::fabs(correctedX[i] - projectedX) > 0.08 ||
          std::fabs(pulseY[i] - projectedY - 0.10) > 0.36)
        continue;
      if (static_cast<int>(std::lround(pixel[i])) < kLayerBoundary)
        layer1.push_back(i);
      else
        layer2.push_back(i);
    }

    std::vector<PairCandidate> candidates;
    candidates.reserve(layer1.size() * layer2.size());
    for (size_t pulse1 : layer1) {
      for (size_t pulse2 : layer2) {
        const double dt = le[pulse2] - le[pulse1];
        const double dx = correctedX[pulse2] - correctedX[pulse1];
        const double dy = pulseY[pulse2] - pulseY[pulse1];
        if (std::fabs(dt) > 15.0 || std::fabs(dx) > 0.15 ||
            std::fabs(dy) > 0.08)
          continue;
        candidates.push_back(
            {pulse1, pulse2, dt * dt + (dx / 0.01) * (dx / 0.01)});
      }
    }
    std::stable_sort(candidates.begin(), candidates.end(),
                     [](const PairCandidate &left,
                        const PairCandidate &right) {
                       return left.score < right.score;
                     });

    std::vector<bool> used(nPulses, false);
    Long64_t eventSelectedPairs = 0;
    for (const PairCandidate &candidate : candidates) {
      if (used[candidate.pulse1] || used[candidate.pulse2])
        continue;
      used[candidate.pulse1] = true;
      used[candidate.pulse2] = true;
      const double deltaZ =
          pulseZ[candidate.pulse2] - pulseZ[candidate.pulse1];
      const double trajectoryResidual =
          correctedX[candidate.pulse2] - correctedX[candidate.pulse1] -
          (*ecalX / kECalZFromTargetM) * deltaZ;
      const double timingResidual =
          *ecalTime - 0.5 * (le[candidate.pulse1] + le[candidate.pulse2]);
      const double normalizedX =
          (trajectoryResidual - pairResidualCenterM) / pairResidualScaleM;
      const double normalizedTime =
          (timingResidual - pairTimingCenterNs) / pairTimingScaleNs;
      if (normalizedX * normalizedX + normalizedTime * normalizedTime <=
          pairCutRadius * pairCutRadius)
        ++eventSelectedPairs;
    }
    if (eventSelectedPairs <= 0)
      continue;
    for (size_t window = firstWindow; window < halfWidths.size(); ++window) {
      ++eventsWithPairs[window];
      selectedPairs[window] += eventSelectedPairs;
    }
    for (size_t window = firstEnergyWindow;
         window < energyHalfWidths.size(); ++window) {
      ++energyEventsWithPairs[window];
      energySelectedPairs[window] += eventSelectedPairs;
    }
  }

  gSystem->mkdir(outputDirectory, true);
  TGraphErrors eventFraction(halfWidths.size());
  TGraphErrors pairsPerEvent(halfWidths.size());
  TGraphErrors absoluteTimingRetention(halfWidths.size());
  for (size_t i = 0; i < halfWidths.size(); ++i) {
    const double denominator = denominators[i];
    const double fraction = denominator > 0.0
        ? static_cast<double>(eventsWithPairs[i]) / denominator
        : 0.0;
    const double fractionError = denominator > 0.0
        ? std::sqrt(fraction * (1.0 - fraction) / denominator)
        : 0.0;
    const double pairRatio = denominator > 0.0
        ? static_cast<double>(selectedPairs[i]) / denominator
        : 0.0;
    const double pairRatioError = denominator > 0.0
        ? std::sqrt(static_cast<double>(selectedPairs[i])) / denominator
        : 0.0;
    eventFraction.SetPoint(i, halfWidths[i], fraction);
    eventFraction.SetPointError(i, 0.0, fractionError);
    pairsPerEvent.SetPoint(i, halfWidths[i], pairRatio);
    pairsPerEvent.SetPointError(i, 0.0, pairRatioError);
    const double absoluteRetention = energyOnlyReferenceEvents > 0
        ? static_cast<double>(eventsWithPairs[i]) / energyOnlyReferenceEvents
        : 0.0;
    absoluteTimingRetention.SetPoint(i, halfWidths[i], absoluteRetention);
  }
  eventFraction.SetName("gCDetPairEventFractionVsECalTimeHalfWidth");
  eventFraction.SetTitle(
      "ECal-energy events with #geq1 selected CDet pair;"
      "Symmetric ECal timing half-width (ns);Event fraction");
  pairsPerEvent.SetName("gCDetPairsPerECalEventVsECalTimeHalfWidth");
  pairsPerEvent.SetTitle(
      "Selected CDet pair yield per ECal-energy event;"
      "Symmetric ECal timing half-width (ns);Pairs / ECal-energy event");
  absoluteTimingRetention.SetName(
      "gCDetAbsolutePairEventRetentionVsECalTimeHalfWidth");
  absoluteTimingRetention.SetTitle(Form(
      "Absolute retention relative to all %.1f--%.1f GeV events;"
      "Symmetric ECal timing half-width (ns);Pair events / %lld energy events",
      ecalEnergyMinGeV, ecalEnergyMaxGeV, energyOnlyReferenceEvents));
  for (TGraphErrors *graph : {&eventFraction, &pairsPerEvent,
                              &absoluteTimingRetention}) {
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlue + 1);
    graph->SetLineColor(kBlue + 1);
  }

  TCanvas canvas("cCDetPairYieldVsECalTimingWindow",
                 "CDet pair yield vs ECal timing window", 1800, 500);
  canvas.Divide(3, 1);
  canvas.cd(1);
  eventFraction.Draw("APL");
  canvas.cd(2);
  pairsPerEvent.Draw("APL");
  canvas.cd(3);
  absoluteTimingRetention.Draw("APL");
  canvas.SaveAs(Form("%s/CDetPairYieldVsECalTimingWindow.pdf",
                     outputDirectory));
  canvas.SaveAs(Form("%s/CDetPairYieldVsECalTimingWindow.png",
                     outputDirectory));

  TFile output(Form("%s/CDetPairYieldVsECalTimingWindow.root",
                    outputDirectory), "RECREATE");
  eventFraction.Write();
  pairsPerEvent.Write();
  absoluteTimingRetention.Write();
  output.Close();

  std::ofstream table(Form("%s/CDetPairYieldVsECalTimingWindow.csv",
                           outputDirectory));
  table << "half_width_ns,ecal_energy_events,events_with_selected_pair,"
           "selected_pairs,event_fraction,pairs_per_ecal_event\n";
  for (size_t i = 0; i < halfWidths.size(); ++i) {
    const double denominator = denominators[i];
    table << halfWidths[i] << ',' << denominators[i] << ','
          << eventsWithPairs[i] << ',' << selectedPairs[i] << ','
          << (denominator > 0.0 ? eventsWithPairs[i] / denominator : 0.0)
          << ','
          << (denominator > 0.0 ? selectedPairs[i] / denominator : 0.0)
          << '\n';
  }

  TGraphErrors energyEventFraction(energyHalfWidths.size());
  TGraphErrors energyPairsPerEvent(energyHalfWidths.size());
  TGraphErrors absoluteEnergyRetention(energyHalfWidths.size());
  const double fullEnergyTimingReference = energyDenominators.empty()
      ? 0.0
      : static_cast<double>(energyDenominators.back());
  for (size_t i = 0; i < energyHalfWidths.size(); ++i) {
    const double denominator = energyDenominators[i];
    const double fraction = denominator > 0.0
        ? static_cast<double>(energyEventsWithPairs[i]) / denominator
        : 0.0;
    const double fractionError = denominator > 0.0
        ? std::sqrt(fraction * (1.0 - fraction) / denominator)
        : 0.0;
    const double pairRatio = denominator > 0.0
        ? static_cast<double>(energySelectedPairs[i]) / denominator
        : 0.0;
    const double pairRatioError = denominator > 0.0
        ? std::sqrt(static_cast<double>(energySelectedPairs[i])) / denominator
        : 0.0;
    energyEventFraction.SetPoint(i, energyHalfWidths[i], fraction);
    energyEventFraction.SetPointError(i, 0.0, fractionError);
    energyPairsPerEvent.SetPoint(i, energyHalfWidths[i], pairRatio);
    energyPairsPerEvent.SetPointError(i, 0.0, pairRatioError);
    absoluteEnergyRetention.SetPoint(
        i, energyHalfWidths[i],
        fullEnergyTimingReference > 0.0
            ? energyEventsWithPairs[i] / fullEnergyTimingReference
            : 0.0);
  }
  energyEventFraction.SetName("gCDetPairEventFractionVsECalEnergyHalfWidth");
  energyEventFraction.SetTitle(Form(
      "Events with #geq1 selected CDet pair;"
      "ECal energy half-width about %.2f GeV;Pair events / events in window",
      energyCenterGeV));
  energyPairsPerEvent.SetName("gCDetPairsPerEventVsECalEnergyHalfWidth");
  energyPairsPerEvent.SetTitle(Form(
      "Selected CDet pair yield;"
      "ECal energy half-width about %.2f GeV;Pairs / events in window",
      energyCenterGeV));
  absoluteEnergyRetention.SetName(
      "gCDetAbsolutePairEventRetentionVsECalEnergyHalfWidth");
  absoluteEnergyRetention.SetTitle(Form(
      "Absolute retention relative to %.1f--%.1f GeV reference;"
      "ECal energy half-width about %.2f GeV;Pair events / reference events",
      ecalEnergyMinGeV, ecalEnergyMaxGeV, energyCenterGeV));
  for (TGraphErrors *graph : {&energyEventFraction, &energyPairsPerEvent,
                              &absoluteEnergyRetention}) {
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlue + 1);
    graph->SetLineColor(kBlue + 1);
  }
  TCanvas energyCanvas("cCDetPairYieldVsECalEnergyWindow",
                       "CDet pair yield vs ECal energy window", 1800, 500);
  energyCanvas.Divide(3, 1);
  energyCanvas.cd(1);
  energyEventFraction.Draw("APL");
  energyCanvas.cd(2);
  energyPairsPerEvent.Draw("APL");
  energyCanvas.cd(3);
  absoluteEnergyRetention.Draw("APL");
  energyCanvas.SaveAs(Form("%s/CDetPairYieldVsECalEnergyWindow.pdf",
                           outputDirectory));
  energyCanvas.SaveAs(Form("%s/CDetPairYieldVsECalEnergyWindow.png",
                           outputDirectory));

  TFile energyOutput(Form("%s/CDetPairYieldVsECalEnergyWindow.root",
                          outputDirectory), "RECREATE");
  energyEventFraction.Write();
  energyPairsPerEvent.Write();
  absoluteEnergyRetention.Write();
  energyOutput.Close();

  std::ofstream energyTable(Form("%s/CDetPairYieldVsECalEnergyWindow.csv",
                                 outputDirectory));
  energyTable << "energy_center_gev,half_width_gev,energy_min_gev,"
                 "energy_max_gev,ecal_events,events_with_selected_pair,"
                 "selected_pairs,event_fraction,pairs_per_event,"
                 "absolute_event_retention\n";
  for (size_t i = 0; i < energyHalfWidths.size(); ++i) {
    const double denominator = energyDenominators[i];
    energyTable << energyCenterGeV << ',' << energyHalfWidths[i] << ','
                << energyCenterGeV - energyHalfWidths[i] << ','
                << energyCenterGeV + energyHalfWidths[i] << ','
                << energyDenominators[i] << ',' << energyEventsWithPairs[i]
                << ',' << energySelectedPairs[i] << ','
                << (denominator > 0.0
                        ? energyEventsWithPairs[i] / denominator
                        : 0.0)
                << ','
                << (denominator > 0.0
                        ? energySelectedPairs[i] / denominator
                        : 0.0)
                << ','
                << (fullEnergyTimingReference > 0.0
                        ? energyEventsWithPairs[i] /
                              fullEnergyTimingReference
                        : 0.0)
                << '\n';
  }

  std::cout << "[CDet pair-window scan] Processed " << eventCount
            << " events from " << filesAdded << " files." << std::endl;
  std::cout << "[CDet pair-window scan] Output directory: "
            << outputDirectory << std::endl;
}
