#include <TCanvas.h>
#include <TChain.h>
#include <TEllipse.h>
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
#include <iostream>
#include <limits>
#include <unordered_set>
#include <vector>

#include "CDetRunDataset.h"

namespace {

const int kCDetNPixels = 2688;
const int kCDetNBars = 168;
const int kPixelsPerBar = 16;

bool HasRequiredBranches(TChain &chain) {
  const char *required[] = {
      "earm.cdet.pulse.pmtnum", "earm.cdet.pulse.tdc_le_corr",
      "earm.cdet.pulse.tdc_te_corr", "earm.cdet.pulse.tdc_tot_ns",
      "earm.cdet.pulse.ecal_residual",
      "earm.cdet.pulse.calib_valid",
      "earm.cdet.pulse.broad_quality_pass",
      "earm.cdet.pulse.ecal_eligible", "earm.cdet.pulse.spatial_pass",
      "earm.cdet.pulse.x_corr", "earm.cdet.pulse.z",
      "earm.cdet.pair.pulse_index_l1", "earm.cdet.pair.pulse_index_l2"};

  chain.LoadTree(0);
  for (const char *name : required) {
    if (!chain.GetBranch(name)) {
      std::cerr << "[good-pulse TDC] Missing required branch: " << name
                << std::endl;
      return false;
    }
  }
  return true;
}

void DrawBarPage(TCanvas *canvas, const std::vector<TH1D *> &histograms,
                 int firstBar, int lastBar) {
  canvas->Clear();
  canvas->Divide(7, 6, 0.001, 0.001);
  for (int bar = firstBar; bar <= lastBar; ++bar) {
    canvas->cd(bar - firstBar + 1);
    TH1D *hist = histograms[bar];
    hist->SetLineColor(kBlue + 1);
    hist->SetStats(false);
    hist->Draw();

    if (hist->GetEntries() >= 20 && hist->GetRMS() > 0.0) {
      const double fitLow = std::max(hist->GetXaxis()->GetXmin(),
                                     hist->GetMean() - 2.0 * hist->GetRMS());
      const double fitHigh = std::min(hist->GetXaxis()->GetXmax(),
                                      hist->GetMean() + 2.0 * hist->GetRMS());
      if (fitHigh > fitLow) {
        TF1 fit(Form("fGoodPulseBar%d", bar), "gaus", fitLow, fitHigh);
        fit.SetLineColor(kRed + 1);
        hist->Fit(&fit, "QNR");
        fit.DrawCopy("same");
      }
    }
  }
  canvas->Modified();
  canvas->Update();
}

} // namespace

// Reproduce the plotAllTDC timing views using only calibrated good pulse
// candidates.  A pulse is accepted when all four analyzer flags are true:
// calib_valid, broad_quality_pass, ecal_eligible, and spatial_pass.  An event
// contributes only if its accepted candidates include both CDet layers.
//
// Example:
// root -l -b -q 'Plot_CDet_GoodPulseCandidates_AllTDC.C+(5710,"/path/to/Rootfiles","CDet_run5710_good_pulse_tdc")'
void Plot_CDet_GoodPulseCandidates_AllTDC(
    int runNumber = 5710, const char *inputDirectory = nullptr,
    const char *outputDirectory = "CDet_good_pulse_tdc",
    double binWidthNs = 1.0, double leMinNs = 0.0, double leMaxNs = 60.0,
    double totMinNs = 0.0, double totMaxNs = 40.0,
    bool recoveredOnly = false, double pairResidualCenterM = 0.0,
    double pairTimingCenterNs = -26.0, double pairResidualScaleM = 0.020,
    double pairTimingScaleNs = 5.0, double pairCutRadius = 2.0) {
  TString resolvedInputDirectory(inputDirectory ? inputDirectory : "");
  if (resolvedInputDirectory.IsNull()) {
    const char *outDir = gSystem->Getenv("OUT_DIR");
    if (outDir)
      resolvedInputDirectory = outDir;
  }
  if (resolvedInputDirectory.IsNull()) {
    std::cerr << "[good-pulse TDC] An input directory or OUT_DIR is required."
              << std::endl;
    return;
  }
  if (binWidthNs <= 0.0 || leMaxNs <= leMinNs || totMaxNs <= totMinNs ||
      pairResidualScaleM <= 0.0 || pairTimingScaleNs <= 0.0 ||
      pairCutRadius <= 0.0) {
    std::cerr << "[good-pulse TDC] Invalid histogram limits." << std::endl;
    return;
  }

  TChain chain("T");
  const int filesAdded = CDetRunDataset::AddToChain(
      &chain, runNumber, resolvedInputDirectory.Data());
  if (filesAdded <= 0 || chain.GetEntries() <= 0) {
    std::cerr << "[good-pulse TDC] No events found for run " << runNumber
              << " in " << resolvedInputDirectory << std::endl;
    return;
  }
  if (!HasRequiredBranches(chain))
    return;

  const int nLEBins = std::max(1, int(std::ceil((leMaxNs - leMinNs) /
                                                binWidthNs)));
  const double teMaxNs = leMaxNs + totMaxNs;
  const int nTEBins = std::max(1, int(std::ceil((teMaxNs - leMinNs) /
                                                binWidthNs)));
  const int nToTBins = std::max(1, int(std::ceil((totMaxNs - totMinNs) /
                                                 binWidthNs)));
  const int bar30FirstPixel = 30 * kPixelsPerBar;
  const double bar30DeltaTMinNs = -60.0;
  const double bar30DeltaTMaxNs = 30.0;
  const int nBar30DeltaTBins = std::max(
      1, int(std::ceil((bar30DeltaTMaxNs - bar30DeltaTMinNs) / binWidthNs)));

  const char *populationTitle = recoveredOnly
      ? "Recovered CDet pulses absent from legacy hAllGoodLe"
      : "Good CDet pulse candidates, two-layer events";
  TH1D hGoodLE("hCDetGoodPulseLE",
               Form("%s;Corrected LE time (ns);Pulses", populationTitle),
               nLEBins, leMinNs, leMaxNs);
  TH1D hGoodTE("hCDetGoodPulseTE",
               Form("%s;Corrected TE time (ns);Pulses", populationTitle),
               nTEBins, leMinNs, teMaxNs);
  TH1D hGoodToT("hCDetGoodPulseToT",
                Form("%s;Corrected ToT (ns);Pulses", populationTitle),
                nToTBins, totMinNs, totMaxNs);
  TH1D hGoodPixel("hCDetGoodPulsePixel",
                  "Good CDet pulse candidates;Pixel ID;Pulses",
                  kCDetNPixels, -0.5, kCDetNPixels - 0.5);
  TH1D hGoodBar("hCDetGoodPulseBar",
                "Good CDet pulse candidates;Bar ID;Pulses", kCDetNBars,
                -0.5, kCDetNBars - 0.5);
  TH1D hGoodMultiplicity(
      "hCDetGoodPulseMultiplicity",
      "Good CDet pulse candidates per event;Candidate pulses;Events", 101,
      -0.5, 100.5);
  TH1D hBar30ECalCut(
      "hCDetBar30ECalMinusCDet_ECalCut",
      "After ECal energy cut, CDet bar 30 (all instrumented pixels);"
      "t_{ECal} - t_{CDet,corr} (ns);Selected hits",
      nBar30DeltaTBins, bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH1D hBar30Projected(
      "hCDetBar30ECalMinusCDet_Projected",
      "ECal projection in CDet bar 30;"
      "t_{ECal} - t_{CDet,corr} (ns);Trajectory-matched hits",
      nBar30DeltaTBins, bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH1D hBar30ProjectedQuality(
      "hCDetBar30ECalMinusCDet_ProjectedQuality",
      "ECal projection + pulse-quality selection, CDet bar 30;"
      "t_{ECal} - t_{CDet,corr} (ns);Selected hits",
      nBar30DeltaTBins, bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH2D hBar30ProjectedQualityVsToT(
      "hCDetBar30ECalMinusCDetVsToT_ProjectedQuality",
      "ECal projection + pulse-quality selection, CDet bar 30;"
      "CDet ToT (ns);t_{ECal} - t_{CDet,corr} (ns)",
      nToTBins, totMinNs, totMaxNs, nBar30DeltaTBins,
      bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH1D hDetectorECalCut(
      "hCDetDetectorECalMinusCDet_ECalCut",
      "ECal energy cut, events with an accepted CDet pair;"
      "t_{ECal} - t_{CDet,corr} (ns);Selected hits",
      nBar30DeltaTBins, bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH1D hDetectorProjected(
      "hCDetDetectorECalMinusCDet_Projected",
      "Projection-matched pulses, events with an accepted CDet pair;"
      "t_{ECal} - t_{CDet,corr} (ns);Trajectory-matched hits",
      nBar30DeltaTBins, bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH1D hDetectorProjectedQuality(
      "hCDetDetectorECalMinusCDet_ProjectedQuality",
      "Trajectory-time selected Layer-1/Layer-2 pair members, all CDet;"
      "t_{ECal} - t_{CDet,corr} (ns);Paired hits",
      nBar30DeltaTBins, bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH2D hDetectorProjectedQualityVsToT(
      "hCDetDetectorECalMinusCDetVsToT_ProjectedQuality",
      "Trajectory-time selected Layer-1/Layer-2 pair members, all CDet;"
      "CDet ToT (ns);t_{ECal} - t_{CDet,corr} (ns)",
      nToTBins, totMinNs, totMaxNs, nBar30DeltaTBins,
      bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH1D hSelectedPairMeanResidual(
      "hCDetSelectedPairMeanECalResidual",
      "Trajectory-time selected pair mean, all CDet;"
      "t_{ECal} - <t_{CDet,corr}>_{pair} (ns);Selected pairs",
      nBar30DeltaTBins, bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH2D hSelectedPairMeanResidualVsMeanToT(
      "hCDetSelectedPairMeanECalResidualVsMeanToT",
      "Trajectory-time selected pair mean, all CDet;"
      "<CDet ToT>_{pair} (ns);t_{ECal} - <t_{CDet,corr}>_{pair} (ns)",
      nToTBins, totMinNs, totMaxNs, nBar30DeltaTBins,
      bar30DeltaTMinNs, bar30DeltaTMaxNs);
  TH2D hSelectedPairXCorrelation(
      "hCDetSelectedPairXVsProjectedECalX",
      "Trajectory-time selected pairs;"
      "x_{ECal} projected to pair mean z (m);"
      "<x_{CDet,corr}>_{pair} (m)",
      160, -1.6, 1.6, 160, -1.6, 1.6);
  TH1D hPairTrajectoryResidual(
      "hCDetPairTrajectoryResidual",
      "All accepted pairs;#Delta x_{pair} - (x_{ECal}/z_{ECal})#Delta z (m);Pairs",
      160, -0.08, 0.08);
  TH1D hBestPairTrajectoryResidual(
      "hCDetBestPairTrajectoryResidual",
      "Best trajectory-matched pair per event;#Delta x_{pair} - (x_{ECal}/z_{ECal})#Delta z (m);Events",
      160, -0.08, 0.08);
  TH2D hPairTimingVsTrajectoryResidual(
      "hCDetPairTimingVsTrajectoryResidual",
      "Accepted pairs;#Delta x_{pair} - (x_{ECal}/z_{ECal})#Delta z (m);"
      "t_{ECal} - <t_{CDet}>_{pair} (ns)",
      160, -0.08, 0.08, nBar30DeltaTBins, bar30DeltaTMinNs,
      bar30DeltaTMaxNs);

  std::vector<TH1D *> pixelLE(kCDetNPixels, nullptr);
  for (int pixel = 0; pixel < kCDetNPixels; ++pixel) {
    pixelLE[pixel] =
        new TH1D(Form("hCDetGoodPulseLE_pixel%04d", pixel),
                 Form("Good pulse LE, pixel %d;Corrected LE time (ns);Pulses",
                      pixel),
                 nLEBins, leMinNs, leMaxNs);
  }

  std::vector<TH1D *> barLE(kCDetNBars, nullptr);
  for (int bar = 0; bar < kCDetNBars; ++bar) {
    barLE[bar] =
        new TH1D(Form("hCDetGoodPulseLE_bar%03d", bar),
                 Form("Good pulse LE, bar %d;Corrected LE time (ns);Pulses",
                      bar),
                 nLEBins, leMinNs, leMaxNs);
  }

  TTreeReader reader(&chain);
  TTreeReaderArray<Double_t> pixel(reader, "earm.cdet.pulse.pmtnum");
  TTreeReaderArray<Double_t> le(reader, "earm.cdet.pulse.tdc_le_corr");
  TTreeReaderArray<Double_t> te(reader, "earm.cdet.pulse.tdc_te_corr");
  TTreeReaderArray<Double_t> tot(reader, "earm.cdet.pulse.tdc_tot_ns");
  TTreeReaderArray<Double_t> ecalResidual(
      reader, "earm.cdet.pulse.ecal_residual");
  TTreeReaderArray<Double_t> calibValid(reader,
                                        "earm.cdet.pulse.calib_valid");
  TTreeReaderArray<Double_t> broadQuality(
      reader, "earm.cdet.pulse.broad_quality_pass");
  TTreeReaderArray<Double_t> ecalEligible(reader,
                                          "earm.cdet.pulse.ecal_eligible");
  TTreeReaderArray<Double_t> spatialPass(reader,
                                         "earm.cdet.pulse.spatial_pass");
  TTreeReaderArray<Double_t> correctedX(reader, "earm.cdet.pulse.x_corr");
  TTreeReaderArray<Double_t> pulseZ(reader, "earm.cdet.pulse.z");
  TTreeReaderArray<Double_t> pairPulseIndexL1(
      reader, "earm.cdet.pair.pulse_index_l1");
  TTreeReaderArray<Double_t> pairPulseIndexL2(
      reader, "earm.cdet.pair.pulse_index_l2");
  TTreeReaderArray<Double_t> legacyPMT(reader, "earm.cdet.hit.pmtnum");
  TTreeReaderArray<Double_t> legacyLE(reader, "earm.cdet.hit.tdc_le");
  TTreeReaderArray<Double_t> legacyToT(reader, "earm.cdet.hit.tdc_tot");
  TTreeReaderArray<Double_t> legacyX(reader, "earm.cdet.hit.xhit");
  TTreeReaderArray<Double_t> legacyY(reader, "earm.cdet.hit.yhit");
  TTreeReaderArray<Double_t> legacyZ(reader, "earm.cdet.hit.zhit");
  TTreeReaderArray<Double_t> legacyMultiplicity(reader, "earm.cdet.tdc_mult");
  TTreeReaderValue<Double_t> ecalX(reader, "earm.ecal.x");
  TTreeReaderValue<Double_t> ecalY(reader, "earm.ecal.y");
  TTreeReaderValue<Double_t> ecalTime(reader, "earm.ecal.adctime");
  TTreeReaderValue<Double_t> ecalEnergy(reader, "earm.ecal.e");

  Long64_t eventCount = 0;
  Long64_t pulseCount = 0;
  Long64_t goodPulseCount = 0;
  Long64_t twoLayerEventCount = 0;
  Long64_t plottedEventCount = 0;
  Long64_t malformedEvents = 0;
  Long64_t ecalEnergyEventCount = 0;
  Long64_t detectorBaselineEventCount = 0;
  Long64_t detectorProjectedEventCount = 0;
  Long64_t detectorProjectedQualityEventCount = 0;
  Long64_t acceptedPairCount = 0;
  Long64_t trajectoryTimeSelectedPairCount = 0;
  constexpr double kECalZFromTargetM = 6.144;

  while (reader.Next()) {
    ++eventCount;
    const std::vector<size_t> sizes = {
        pixel.GetSize(),        le.GetSize(),          te.GetSize(),
        tot.GetSize(),          ecalResidual.GetSize(), calibValid.GetSize(),
        broadQuality.GetSize(), ecalEligible.GetSize(), spatialPass.GetSize(),
        correctedX.GetSize(), pulseZ.GetSize()};
    const size_t nPulses = *std::min_element(sizes.begin(), sizes.end());
    const size_t maxPulses = *std::max_element(sizes.begin(), sizes.end());
    if (nPulses != maxPulses)
      ++malformedEvents;
    pulseCount += nPulses;

    // Analyzer-native reproduction of the historical Bar 30 timing canvas.
    // The saved residual has the desired historical sign, ECal minus CDet.
    if (std::isfinite(*ecalEnergy) && *ecalEnergy >= 3.0 &&
        *ecalEnergy <= 4.5) {
      ++ecalEnergyEventCount;
      std::unordered_set<size_t> pairedPulseIndices;
      std::unordered_set<size_t> trajectoryTimeSelectedPulseIndices;
      const size_t nPairs =
          std::min(pairPulseIndexL1.GetSize(), pairPulseIndexL2.GetSize());
      double bestTrajectoryAbsResidual = std::numeric_limits<double>::infinity();
      double bestTrajectoryResidual = std::numeric_limits<double>::quiet_NaN();
      for (size_t pair = 0; pair < nPairs; ++pair) {
        if (!std::isfinite(pairPulseIndexL1[pair]) ||
            !std::isfinite(pairPulseIndexL2[pair]))
          continue;
        const Long64_t indexL1 = std::llround(pairPulseIndexL1[pair]);
        const Long64_t indexL2 = std::llround(pairPulseIndexL2[pair]);
        if (indexL1 < 0 || indexL2 < 0 ||
            indexL1 >= static_cast<Long64_t>(nPulses) ||
            indexL2 >= static_cast<Long64_t>(nPulses))
          continue;
        pairedPulseIndices.insert(static_cast<size_t>(indexL1));
        pairedPulseIndices.insert(static_cast<size_t>(indexL2));
        const double deltaZ = pulseZ[indexL2] - pulseZ[indexL1];
        const double trajectoryResidual =
            (correctedX[indexL2] - correctedX[indexL1]) -
            (*ecalX / kECalZFromTargetM) * deltaZ;
        const double pairTimingResidual =
            0.5 * (ecalResidual[indexL1] + ecalResidual[indexL2]);
        if (std::isfinite(trajectoryResidual) &&
            std::isfinite(pairTimingResidual)) {
          hPairTrajectoryResidual.Fill(trajectoryResidual);
          hPairTimingVsTrajectoryResidual.Fill(trajectoryResidual,
                                               pairTimingResidual);
          if (std::fabs(trajectoryResidual) < bestTrajectoryAbsResidual) {
            bestTrajectoryAbsResidual = std::fabs(trajectoryResidual);
            bestTrajectoryResidual = trajectoryResidual;
          }
          const double normalizedResidual =
              (trajectoryResidual - pairResidualCenterM) /
              pairResidualScaleM;
          const double normalizedTiming =
              (pairTimingResidual - pairTimingCenterNs) /
              pairTimingScaleNs;
          if (normalizedResidual * normalizedResidual +
                  normalizedTiming * normalizedTiming <=
              pairCutRadius * pairCutRadius) {
            trajectoryTimeSelectedPulseIndices.insert(
                static_cast<size_t>(indexL1));
            trajectoryTimeSelectedPulseIndices.insert(
                static_cast<size_t>(indexL2));
            hSelectedPairMeanResidual.Fill(pairTimingResidual);
            const double pairMeanToT =
                0.5 * (tot[indexL1] + tot[indexL2]);
            if (std::isfinite(pairMeanToT))
              hSelectedPairMeanResidualVsMeanToT.Fill(pairMeanToT,
                                                      pairTimingResidual);
            const double pairMeanZ =
                0.5 * (pulseZ[indexL1] + pulseZ[indexL2]);
            const double pairMeanX =
                0.5 * (correctedX[indexL1] + correctedX[indexL2]);
            const double projectedECalX =
                *ecalX * pairMeanZ / kECalZFromTargetM;
            if (std::isfinite(pairMeanX) && std::isfinite(projectedECalX))
              hSelectedPairXCorrelation.Fill(projectedECalX, pairMeanX);
            ++trajectoryTimeSelectedPairCount;
          }
        }
        ++acceptedPairCount;
      }
      if (std::isfinite(bestTrajectoryResidual))
        hBestPairTrajectoryResidual.Fill(bestTrajectoryResidual);
      const bool hasAcceptedPair = !pairedPulseIndices.empty();
      bool hasDetectorBaseline = false;
      bool hasDetectorProjected = false;
      for (size_t i = 0; i < nPulses; ++i) {
        if (!(calibValid[i] > 0.5 && ecalEligible[i] > 0.5) ||
            !std::isfinite(pixel[i]) || !std::isfinite(ecalResidual[i]) ||
            !std::isfinite(tot[i]))
          continue;
        const int pixelID = int(std::lround(pixel[i]));
        if (pixelID < 0 || pixelID >= kCDetNPixels)
          continue;
        if (hasAcceptedPair) {
          hDetectorECalCut.Fill(ecalResidual[i]);
          hasDetectorBaseline = true;
          if (spatialPass[i] > 0.5) {
            hDetectorProjected.Fill(ecalResidual[i]);
            hasDetectorProjected = true;
          }
          if (trajectoryTimeSelectedPulseIndices.count(i) != 0U) {
            hDetectorProjectedQuality.Fill(ecalResidual[i]);
            hDetectorProjectedQualityVsToT.Fill(tot[i], ecalResidual[i]);
          }
        }
        if (pixelID < bar30FirstPixel ||
            pixelID >= bar30FirstPixel + kPixelsPerBar)
          continue;
        hBar30ECalCut.Fill(ecalResidual[i]);
        if (!(spatialPass[i] > 0.5))
          continue;
        hBar30Projected.Fill(ecalResidual[i]);
        if (!(broadQuality[i] > 0.5))
          continue;
        hBar30ProjectedQuality.Fill(ecalResidual[i]);
        hBar30ProjectedQualityVsToT.Fill(tot[i], ecalResidual[i]);
      }
      if (hasDetectorBaseline)
        ++detectorBaselineEventCount;
      if (hasDetectorProjected)
        ++detectorProjectedEventCount;
      if (!trajectoryTimeSelectedPulseIndices.empty())
        ++detectorProjectedQualityEventCount;
    }

    std::vector<size_t> acceptedIndices;
    bool hasLayer1 = false;
    bool hasLayer2 = false;
    for (size_t i = 0; i < nPulses; ++i) {
      if (!(calibValid[i] > 0.5 && broadQuality[i] > 0.5 &&
            ecalEligible[i] > 0.5 && spatialPass[i] > 0.5))
        continue;
      if (!std::isfinite(pixel[i]) || !std::isfinite(le[i]) ||
          !std::isfinite(te[i]) || !std::isfinite(tot[i]) ||
          !std::isfinite(ecalResidual[i]))
        continue;

      const int pixelID = int(std::lround(pixel[i]));
      if (pixelID < 0 || pixelID >= kCDetNPixels)
        continue;
      acceptedIndices.push_back(i);
      hasLayer1 |= pixelID < kCDetNPixels / 2;
      hasLayer2 |= pixelID >= kCDetNPixels / 2;
    }

    if (!(hasLayer1 && hasLayer2)) {
      hGoodMultiplicity.Fill(0);
      continue;
    }

    std::unordered_set<int> legacyAcceptedIDs;
    if (recoveredOnly) {
      const size_t nLegacy = std::min(
          {legacyPMT.GetSize(), legacyLE.GetSize(), legacyToT.GetSize(),
           legacyX.GetSize(), legacyY.GetSize(), legacyZ.GetSize(),
           legacyMultiplicity.GetSize()});
      std::vector<bool> legacyBasic(nLegacy, false);
      int legacyLayer1 = 0;
      int legacyLayer2 = 0;
      for (size_t i = 0; i < nLegacy; ++i) {
        const int id = static_cast<int>(legacyPMT[i]);
        const int layer = id < kCDetNPixels / 2 ? 0 : 1;
        const double correctedX = legacyX[i] * 1.08 - 0.03;
        const bool pass = *ecalY > -1.2 && *ecalY < 1.2 &&
            *ecalX > -1.5 && *ecalX < 1.5 && *ecalX != 0.0 &&
            *ecalY != 0.0 && legacyLE[i] * 0.01 >= 0.02 &&
            legacyLE[i] * 0.01 <= 60.0 && legacyToT[i] * 0.01 >= 4.0 &&
            legacyToT[i] * 0.01 <= 30.0 && legacyMultiplicity[i] < 100.0 &&
            std::fabs(correctedX - (*ecalX) * legacyZ[i] / 6.144) <= 0.08 &&
            std::fabs(legacyY[i] - (*ecalY) * legacyZ[i] / 6.144 - 0.10) <=
                0.36 &&
            *ecalTime > 10.0 && *ecalTime < 35.0;
        legacyBasic[i] = pass;
        if (pass) {
          if (layer == 0)
            ++legacyLayer1;
          else
            ++legacyLayer2;
        }
      }
      if (legacyLayer1 >= 1 && legacyLayer1 <= 100 && legacyLayer2 >= 1 &&
          legacyLayer2 <= 100) {
        for (size_t i = 0; i < nLegacy; ++i) {
          if (legacyBasic[i])
            legacyAcceptedIDs.insert(static_cast<int>(legacyPMT[i]));
        }
      }
    }

    ++twoLayerEventCount;
    int plottedThisEvent = 0;
    for (size_t i : acceptedIndices) {
      const int pixelID = int(std::lround(pixel[i]));
      if (recoveredOnly && legacyAcceptedIDs.count(pixelID))
        continue;
      const int barID = pixelID / kPixelsPerBar;
      ++plottedThisEvent;
      ++goodPulseCount;
      hGoodLE.Fill(le[i]);
      hGoodTE.Fill(te[i]);
      hGoodToT.Fill(tot[i]);
      hGoodPixel.Fill(pixelID);
      hGoodBar.Fill(barID);
      pixelLE[pixelID]->Fill(le[i]);
      barLE[barID]->Fill(le[i]);
    }
    if (plottedThisEvent > 0)
      ++plottedEventCount;
    hGoodMultiplicity.Fill(plottedThisEvent);
  }

  if (gSystem->mkdir(outputDirectory, true) != 0 &&
      gSystem->AccessPathName(outputDirectory)) {
    std::cerr << "[good-pulse TDC] Cannot create output directory "
              << outputDirectory << std::endl;
    return;
  }

  gStyle->SetOptStat(1110);
  TCanvas cAll("cCDetGoodPulseAllTDC", "Good pulse candidate timing", 1500,
               900);
  cAll.Divide(2, 2);
  cAll.cd(1);
  hGoodLE.Draw();
  cAll.cd(2);
  hGoodTE.Draw();
  cAll.cd(3);
  hGoodToT.Draw();
  cAll.cd(4);
  hGoodMultiplicity.Draw();
  cAll.SaveAs(Form("%s/CDetGoodPulse_AllTDC.pdf", outputDirectory));

  TCanvas cChannels("cCDetGoodPulseChannels", "Good pulse channels", 1500,
                    700);
  cChannels.Divide(2, 1);
  cChannels.cd(1);
  hGoodPixel.Draw();
  cChannels.cd(2);
  hGoodBar.Draw();
  cChannels.SaveAs(Form("%s/CDetGoodPulse_AllChannels.pdf", outputDirectory));

  auto drawBar30Peak = [](TH1D &hist, const char *fitName) {
    hist.SetLineColor(kBlack);
    hist.SetLineWidth(2);
    hist.Draw();
    if (hist.GetEntries() < 20)
      return;
    const double fitMin = -55.0;
    const double fitMax = -10.0;
    int peakBin = hist.FindBin(-40.0);
    for (int bin = peakBin + 1; bin <= hist.FindBin(-15.0); ++bin) {
      if (hist.GetBinContent(bin) > hist.GetBinContent(peakBin))
        peakBin = bin;
    }
    const double peak = hist.GetBinCenter(peakBin);
    const double background =
        0.5 * (hist.GetBinContent(hist.FindBin(fitMin)) +
               hist.GetBinContent(hist.FindBin(fitMax)));
    TF1 fit(fitName, "gaus(0)+pol1(3)", fitMin, fitMax);
    fit.SetParameters(std::max(1.0, hist.GetBinContent(peakBin) - background),
                      peak, 3.0, background, 0.0);
    fit.SetLineColor(kRed + 1);
    hist.Fit(&fit, "RQN");
    fit.DrawCopy("same");
  };

  TCanvas cBar30("cCDetGoodPulseBar30Amalgamated",
                 "CDet bar 30 amalgamated ECal-CDet timing", 1200, 800);
  cBar30.Divide(2, 2);
  cBar30.cd(1);
  drawBar30Peak(hBar30ECalCut, "fCDetBar30ECalCut");
  cBar30.cd(2);
  drawBar30Peak(hBar30Projected, "fCDetBar30Projected");
  cBar30.cd(3);
  drawBar30Peak(hBar30ProjectedQuality, "fCDetBar30ProjectedQuality");
  cBar30.cd(4);
  hBar30ProjectedQualityVsToT.SetStats(false);
  hBar30ProjectedQualityVsToT.Draw("COLZ");
  cBar30.SaveAs(Form("%s/CDetGoodPulse_Bar030_Amalgamated.pdf",
                     outputDirectory));
  cBar30.SaveAs(Form("%s/CDetGoodPulse_Bar030_Amalgamated.png",
                     outputDirectory));

  TCanvas cDetector("cCDetGoodPulseDetectorAmalgamated",
                    "Whole-detector amalgamated ECal-CDet timing", 1200, 800);
  cDetector.Divide(2, 2);
  cDetector.cd(1);
  drawBar30Peak(hDetectorECalCut, "fCDetDetectorECalCut");
  cDetector.cd(2);
  drawBar30Peak(hDetectorProjected, "fCDetDetectorProjected");
  cDetector.cd(3);
  drawBar30Peak(hSelectedPairMeanResidual,
                "fCDetSelectedPairMeanResidual");
  cDetector.cd(4);
  hSelectedPairMeanResidualVsMeanToT.SetStats(false);
  hSelectedPairMeanResidualVsMeanToT.Draw("COLZ");
  cDetector.SaveAs(Form("%s/CDetGoodPulse_Detector_Amalgamated.pdf",
                        outputDirectory));
  cDetector.SaveAs(Form("%s/CDetGoodPulse_Detector_Amalgamated.png",
                        outputDirectory));

  TCanvas cPairSlope("cCDetGoodPulsePairSlopeDiagnostics",
                     "CDet pair trajectory diagnostics", 1200, 500);
  cPairSlope.Divide(3, 1);
  cPairSlope.cd(1);
  hPairTrajectoryResidual.Draw();
  cPairSlope.cd(2);
  hBestPairTrajectoryResidual.Draw();
  cPairSlope.cd(3);
  hPairTimingVsTrajectoryResidual.SetStats(false);
  hPairTimingVsTrajectoryResidual.Draw("COLZ");
  TEllipse pairSelectionEllipse(pairResidualCenterM, pairTimingCenterNs,
                                pairResidualScaleM * pairCutRadius,
                                pairTimingScaleNs * pairCutRadius);
  pairSelectionEllipse.SetFillStyle(0);
  pairSelectionEllipse.SetLineColor(kRed + 1);
  pairSelectionEllipse.SetLineWidth(3);
  pairSelectionEllipse.Draw("same");
  cPairSlope.SaveAs(Form("%s/CDetGoodPulse_PairTrajectoryDiagnostics.pdf",
                         outputDirectory));
  cPairSlope.SaveAs(Form("%s/CDetGoodPulse_PairTrajectoryDiagnostics.png",
                         outputDirectory));

  TCanvas cPairX("cCDetGoodPulseSelectedPairXCorrelation",
                 "Selected CDet pair x versus projected ECal x", 800, 700);
  cPairX.SetRightMargin(0.14);
  hSelectedPairXCorrelation.SetStats(false);
  hSelectedPairXCorrelation.Draw("COLZ");
  TLine pairXDiagonal(-1.6, -1.6, 1.6, 1.6);
  pairXDiagonal.SetLineColor(kRed + 1);
  pairXDiagonal.SetLineWidth(3);
  pairXDiagonal.Draw("same");
  cPairX.SaveAs(Form("%s/CDetGoodPulse_SelectedPairXCorrelation.pdf",
                     outputDirectory));
  cPairX.SaveAs(Form("%s/CDetGoodPulse_SelectedPairXCorrelation.png",
                     outputDirectory));

  const char *pageNames[4] = {"Layer1_Left", "Layer1_Right",
                              "Layer2_Left", "Layer2_Right"};
  TCanvas cBars("cCDetGoodPulseBars", "Good pulse bar timing", 1800, 1200);
  for (int page = 0; page < 4; ++page) {
    const int firstBar = page * 42;
    DrawBarPage(&cBars, barLE, firstBar, firstBar + 41);
    cBars.SaveAs(Form("%s/CDetGoodPulse_Bars_%s.pdf", outputDirectory,
                      pageNames[page]));
  }

  TFile output(Form("%s/CDetGoodPulse_AllTDC.root", outputDirectory),
               "RECREATE");
  hGoodLE.Write();
  hGoodTE.Write();
  hGoodToT.Write();
  hGoodPixel.Write();
  hGoodBar.Write();
  hGoodMultiplicity.Write();
  hBar30ECalCut.Write();
  hBar30Projected.Write();
  hBar30ProjectedQuality.Write();
  hBar30ProjectedQualityVsToT.Write();
  hDetectorECalCut.Write();
  hDetectorProjected.Write();
  hDetectorProjectedQuality.Write();
  hDetectorProjectedQualityVsToT.Write();
  hSelectedPairMeanResidual.Write();
  hSelectedPairMeanResidualVsMeanToT.Write();
  hSelectedPairXCorrelation.Write();
  hPairTrajectoryResidual.Write();
  hBestPairTrajectoryResidual.Write();
  hPairTimingVsTrajectoryResidual.Write();
  for (TH1D *hist : pixelLE)
    hist->Write();
  for (TH1D *hist : barLE)
    hist->Write();
  output.Close();

  std::cout << "\n[good-pulse TDC] Files/events/pulses: " << filesAdded << "/"
            << eventCount << "/" << pulseCount << std::endl;
  std::cout << "[good-pulse TDC] Accepted calibrated good pulse candidates: "
            << goodPulseCount << std::endl;
  std::cout << "[good-pulse TDC] Population mode: "
            << (recoveredOnly ? "recovered-only" : "all new candidates")
            << std::endl;
  std::cout << "[good-pulse TDC] Events with accepted candidates in both layers: "
            << twoLayerEventCount << std::endl;
  std::cout << "[good-pulse TDC] Events contributing to plotted population: "
            << plottedEventCount << std::endl;
  std::cout << "[good-pulse TDC] Events with inconsistent pulse-array sizes: "
            << malformedEvents << std::endl;
  std::cout << "[good-pulse TDC] ECal-energy-selected events: "
            << ecalEnergyEventCount << std::endl;
  std::cout << "[good-pulse TDC] ECal-energy events with >=1 accepted Layer-1/Layer-2 pair: "
            << detectorBaselineEventCount << std::endl;
  std::cout << "[good-pulse TDC] Pair-containing events with >=1 projected CDet pulse: "
            << detectorProjectedEventCount << std::endl;
  std::cout << "[good-pulse TDC] ECal-energy events contributing trajectory-time selected pair members: "
            << detectorProjectedQualityEventCount << std::endl;
  std::cout << "[good-pulse TDC] Accepted one-to-one Layer-1/Layer-2 pairs in ECal-energy events: "
            << acceptedPairCount << std::endl;
  std::cout << "[good-pulse TDC] Trajectory-time ellipse: center=("
            << pairResidualCenterM << " m, " << pairTimingCenterNs
            << " ns), scales=(" << pairResidualScaleM << " m, "
            << pairTimingScaleNs << " ns), radius=" << pairCutRadius
            << std::endl;
  std::cout << "[good-pulse TDC] Pairs inside trajectory-time ellipse: "
            << trajectoryTimeSelectedPairCount << std::endl;
  std::cout << "[good-pulse TDC] Selected pair-mean timing entries: "
            << static_cast<Long64_t>(hSelectedPairMeanResidual.GetEntries())
            << std::endl;
  std::cout << "[good-pulse TDC] Output directory: " << outputDirectory
            << std::endl;

  for (TH1D *hist : pixelLE)
    delete hist;
  for (TH1D *hist : barLE)
    delete hist;
}
