#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
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
      "earm.cdet.pulse.calib_valid",
      "earm.cdet.pulse.broad_quality_pass",
      "earm.cdet.pulse.ecal_eligible", "earm.cdet.pulse.spatial_pass"};

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
    bool recoveredOnly = false) {
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
  if (binWidthNs <= 0.0 || leMaxNs <= leMinNs || totMaxNs <= totMinNs) {
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
  TTreeReaderArray<Double_t> calibValid(reader,
                                        "earm.cdet.pulse.calib_valid");
  TTreeReaderArray<Double_t> broadQuality(
      reader, "earm.cdet.pulse.broad_quality_pass");
  TTreeReaderArray<Double_t> ecalEligible(reader,
                                          "earm.cdet.pulse.ecal_eligible");
  TTreeReaderArray<Double_t> spatialPass(reader,
                                         "earm.cdet.pulse.spatial_pass");
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

  Long64_t eventCount = 0;
  Long64_t pulseCount = 0;
  Long64_t goodPulseCount = 0;
  Long64_t twoLayerEventCount = 0;
  Long64_t plottedEventCount = 0;
  Long64_t malformedEvents = 0;

  while (reader.Next()) {
    ++eventCount;
    const std::vector<size_t> sizes = {
        pixel.GetSize(),        le.GetSize(),          te.GetSize(),
        tot.GetSize(),          calibValid.GetSize(),  broadQuality.GetSize(),
        ecalEligible.GetSize(), spatialPass.GetSize()};
    const size_t nPulses = *std::min_element(sizes.begin(), sizes.end());
    const size_t maxPulses = *std::max_element(sizes.begin(), sizes.end());
    if (nPulses != maxPulses)
      ++malformedEvents;
    pulseCount += nPulses;

    std::vector<size_t> acceptedIndices;
    bool hasLayer1 = false;
    bool hasLayer2 = false;
    for (size_t i = 0; i < nPulses; ++i) {
      if (!(calibValid[i] > 0.5 && broadQuality[i] > 0.5 &&
            ecalEligible[i] > 0.5 && spatialPass[i] > 0.5))
        continue;
      if (!std::isfinite(pixel[i]) || !std::isfinite(le[i]) ||
          !std::isfinite(te[i]) || !std::isfinite(tot[i]))
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
  std::cout << "[good-pulse TDC] Output directory: " << outputDirectory
            << std::endl;

  for (TH1D *hist : pixelLE)
    delete hist;
  for (TH1D *hist : barLE)
    delete hist;
}
