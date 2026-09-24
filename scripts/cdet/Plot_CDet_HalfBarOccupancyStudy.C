#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TString.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "CDetRunDataset.h"

namespace {

constexpr int kCDetPixelsOccupancy = 2688;
constexpr int kCDetPixelsPerHalfBar = 16;
constexpr int kCDetInstrumentedPerHalfBar = 14;
constexpr int kCDetHalfBars =
    kCDetPixelsOccupancy / kCDetPixelsPerHalfBar;
constexpr int kCDetHalfBarsPerLayer = kCDetHalfBars / 2;

bool LoadInstrumentedMask(const char *databaseFile,
                          std::array<bool, kCDetPixelsOccupancy> &instrumented,
                          std::array<double, kCDetPixelsOccupancy> &xPosition) {
  std::ifstream input(databaseFile);
  if (!input) {
    std::cerr << "[CDet occupancy] Cannot open database file " << databaseFile
              << '\n';
    return false;
  }

  bool readingX = false;
  std::vector<double> values;
  std::string line;
  while (std::getline(input, line)) {
    const size_t comment = line.find('#');
    const std::string content = line.substr(0, comment);
    if (!readingX) {
      if (content.find("earm.cdet.xpos") != std::string::npos)
        readingX = true;
      continue;
    }
    if (content.find("earm.cdet.ypos") != std::string::npos)
      break;
    std::istringstream tokens(content);
    double value = 0.0;
    while (tokens >> value)
      values.push_back(value);
  }

  if (values.size() < kCDetPixelsOccupancy) {
    std::cerr << "[CDet occupancy] Database supplies only " << values.size()
              << " xpos values; expected at least " << kCDetPixelsOccupancy
              << ".\n";
    return false;
  }

  instrumented.fill(false);
  for (int pixel = 0; pixel < kCDetPixelsOccupancy; ++pixel) {
    xPosition[pixel] = values[pixel];
    instrumented[pixel] = std::isfinite(values[pixel]) &&
                          std::fabs(values[pixel]) < 900.0;
  }
  for (int halfBar = 0; halfBar < kCDetHalfBars; ++halfBar) {
    int count = 0;
    for (int local = 0; local < kCDetPixelsPerHalfBar; ++local)
      count += instrumented[halfBar * kCDetPixelsPerHalfBar + local] ? 1 : 0;
    if (count != kCDetInstrumentedPerHalfBar) {
      std::cerr << "[CDet occupancy] Half-bar " << halfBar << " has " << count
                << " instrumented xpos entries; expected "
                << kCDetInstrumentedPerHalfBar << ".\n";
      return false;
    }
  }
  return true;
}

double Median(std::vector<double> values) {
  if (values.empty())
    return 0.0;
  std::sort(values.begin(), values.end());
  const size_t middle = values.size() / 2;
  return values.size() % 2 ? values[middle]
                           : 0.5 * (values[middle - 1] + values[middle]);
}

struct HalfBarOccupancy {
  int global = -1;
  int layer = -1;
  int inLayer = -1;
  double completeMin = 0.0;
  double completeMedian = 0.0;
  double completeMax = 0.0;
  double goodMin = 0.0;
  double goodMedian = 0.0;
  double goodMax = 0.0;
  double completeMaxFraction = 0.0;
  double goodMaxFraction = 0.0;
};

} // namespace

void Plot_CDet_HalfBarOccupancyStudy(
    int runNumber = 6077, const char *inputDirectory = nullptr,
    const char *outputDirectory = "CDet_run6077_halfbar_occupancy",
    const char *databaseFile = "../../DB/db_earm.cdet.dat",
    double coherentLowFraction = 0.10) {
  if (runNumber <= 0 || coherentLowFraction <= 0.0) {
    std::cerr << "[CDet occupancy] Invalid run or threshold.\n";
    return;
  }

  std::array<bool, kCDetPixelsOccupancy> instrumented;
  std::array<double, kCDetPixelsOccupancy> xPosition;
  if (!LoadInstrumentedMask(databaseFile, instrumented, xPosition))
    return;

  TString input(inputDirectory ? inputDirectory : "");
  if (input.IsNull()) {
    const char *outDir = gSystem->Getenv("OUT_DIR");
    if (outDir)
      input = outDir;
  }
  if (input.IsNull()) {
    std::cerr << "[CDet occupancy] An input directory or OUT_DIR is required.\n";
    return;
  }

  TChain chain("T");
  const int filesAdded =
      CDetRunDataset::AddToChain(&chain, runNumber, input.Data());
  if (filesAdded <= 0 || chain.GetEntries() <= 0)
    return;
  chain.LoadTree(0);
  const char *required[] = {"earm.cdet.pulse.pmtnum",
                            "earm.cdet.pulse.calib_valid",
                            "earm.cdet.pulse.broad_quality_pass"};
  for (const char *branch : required) {
    if (!chain.GetBranch(branch)) {
      std::cerr << "[CDet occupancy] Missing branch " << branch << '\n';
      return;
    }
  }

  std::array<Long64_t, kCDetPixelsOccupancy> completeCounts{};
  std::array<Long64_t, kCDetPixelsOccupancy> goodCounts{};
  TTreeReader reader(&chain);
  TTreeReaderArray<Double_t> pixel(reader, "earm.cdet.pulse.pmtnum");
  TTreeReaderArray<Double_t> calibValid(reader,
                                        "earm.cdet.pulse.calib_valid");
  TTreeReaderArray<Double_t> broadQuality(
      reader, "earm.cdet.pulse.broad_quality_pass");
  Long64_t events = 0;
  Long64_t ignoredReferenceOrInvalid = 0;
  Long64_t pulsesInSentinelSlots = 0;
  while (reader.Next()) {
    ++events;
    const size_t n = std::min(
        {pixel.GetSize(), calibValid.GetSize(), broadQuality.GetSize()});
    for (size_t i = 0; i < n; ++i) {
      if (!std::isfinite(pixel[i]))
        continue;
      const int id = int(std::lround(pixel[i]));
      if (id < 0 || id >= kCDetPixelsOccupancy) {
        ++ignoredReferenceOrInvalid;
        continue;
      }
      if (!instrumented[id]) {
        ++pulsesInSentinelSlots;
        continue;
      }
      if (calibValid[i] > 0.5)
        ++completeCounts[id];
      if (calibValid[i] > 0.5 && broadQuality[i] > 0.5)
        ++goodCounts[id];
    }
  }

  std::array<double, 2> layerCompleteMedian{};
  std::array<double, 2> layerGoodMedian{};
  for (int layer = 0; layer < 2; ++layer) {
    std::vector<double> complete;
    std::vector<double> good;
    const int begin = layer * kCDetPixelsOccupancy / 2;
    const int end = (layer + 1) * kCDetPixelsOccupancy / 2;
    for (int id = begin; id < end; ++id) {
      if (!instrumented[id])
        continue;
      complete.push_back(completeCounts[id]);
      good.push_back(goodCounts[id]);
    }
    layerCompleteMedian[layer] = Median(complete);
    layerGoodMedian[layer] = Median(good);
  }

  std::vector<HalfBarOccupancy> halfBars;
  halfBars.reserve(kCDetHalfBars);
  for (int halfBar = 0; halfBar < kCDetHalfBars; ++halfBar) {
    std::vector<double> complete;
    std::vector<double> good;
    for (int local = 0; local < kCDetPixelsPerHalfBar; ++local) {
      const int id = halfBar * kCDetPixelsPerHalfBar + local;
      if (!instrumented[id])
        continue;
      complete.push_back(completeCounts[id]);
      good.push_back(goodCounts[id]);
    }
    HalfBarOccupancy result;
    result.global = halfBar;
    result.layer = halfBar / kCDetHalfBarsPerLayer;
    result.inLayer = halfBar % kCDetHalfBarsPerLayer;
    const auto completeRange = std::minmax_element(complete.begin(), complete.end());
    const auto goodRange = std::minmax_element(good.begin(), good.end());
    result.completeMin = *completeRange.first;
    result.completeMedian = Median(complete);
    result.completeMax = *completeRange.second;
    result.goodMin = *goodRange.first;
    result.goodMedian = Median(good);
    result.goodMax = *goodRange.second;
    result.completeMaxFraction = layerCompleteMedian[result.layer] > 0.0
        ? result.completeMax / layerCompleteMedian[result.layer] : 0.0;
    result.goodMaxFraction = layerGoodMedian[result.layer] > 0.0
        ? result.goodMax / layerGoodMedian[result.layer] : 0.0;
    halfBars.push_back(result);
  }

  gSystem->mkdir(outputDirectory, true);
  TH2D hCompleteMap(
      "hCDetCompleteTripletOccupancyMap",
      "Complete calibrated triplets;Half-bar within layer;Local electronic slot",
      kCDetHalfBarsPerLayer, -0.5, kCDetHalfBarsPerLayer - 0.5,
      kCDetPixelsPerHalfBar, -0.5, kCDetPixelsPerHalfBar - 0.5);
  TH2D hGoodMap(
      "hCDetGoodPulseOccupancyMap",
      "Broad-quality good pulses;Half-bar within layer;Local electronic slot",
      kCDetHalfBarsPerLayer, -0.5, kCDetHalfBarsPerLayer - 0.5,
      kCDetPixelsPerHalfBar, -0.5, kCDetPixelsPerHalfBar - 0.5);
  // The two layer maps are placed consecutively along x in the display below.
  TH2D hCompleteBoth(
      "hCDetCompleteTripletOccupancyBothLayers",
      "Complete calibrated triplets;Global half-bar (0-83 L1, 84-167 L2);Local electronic slot",
      kCDetHalfBars, -0.5, kCDetHalfBars - 0.5,
      kCDetPixelsPerHalfBar, -0.5, kCDetPixelsPerHalfBar - 0.5);
  TH2D hGoodBoth(
      "hCDetGoodPulseOccupancyBothLayers",
      "Broad-quality good pulses;Global half-bar (0-83 L1, 84-167 L2);Local electronic slot",
      kCDetHalfBars, -0.5, kCDetHalfBars - 0.5,
      kCDetPixelsPerHalfBar, -0.5, kCDetPixelsPerHalfBar - 0.5);
  for (int id = 0; id < kCDetPixelsOccupancy; ++id) {
    if (!instrumented[id])
      continue;
    const int halfBar = id / kCDetPixelsPerHalfBar;
    const int local = id % kCDetPixelsPerHalfBar;
    hCompleteBoth.SetBinContent(halfBar + 1, local + 1, completeCounts[id]);
    hGoodBoth.SetBinContent(halfBar + 1, local + 1, goodCounts[id]);
  }

  TGraph completeMaximum(kCDetHalfBars);
  TGraph goodMaximum(kCDetHalfBars);
  int lowCompleteBars = 0;
  int lowGoodBars = 0;
  for (int i = 0; i < kCDetHalfBars; ++i) {
    completeMaximum.SetPoint(i, halfBars[i].global,
                             halfBars[i].completeMaxFraction);
    goodMaximum.SetPoint(i, halfBars[i].global, halfBars[i].goodMaxFraction);
    lowCompleteBars +=
        halfBars[i].completeMaxFraction < coherentLowFraction ? 1 : 0;
    lowGoodBars += halfBars[i].goodMaxFraction < coherentLowFraction ? 1 : 0;
  }
  completeMaximum.SetName("gCDetHalfBarCompleteMaximumFraction");
  completeMaximum.SetTitle(
      "Coherent half-bar occupancy;Global half-bar;maximum physical-pixel count / layer median");
  goodMaximum.SetName("gCDetHalfBarGoodMaximumFraction");
  goodMaximum.SetTitle(
      "Coherent good-pulse occupancy;Global half-bar;maximum physical-pixel count / layer median");
  for (TGraph *graph : {&completeMaximum, &goodMaximum}) {
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.6);
  }

  TCanvas mapCanvas("cCDetHalfBarOccupancyMaps",
                    "CDet half-bar occupancy maps", 1800, 900);
  mapCanvas.Divide(1, 2);
  mapCanvas.cd(1); gPad->SetLogz(); hCompleteBoth.Draw("COLZ");
  mapCanvas.cd(2); gPad->SetLogz(); hGoodBoth.Draw("COLZ");
  mapCanvas.SaveAs(Form("%s/CDetHalfBarOccupancyMaps.pdf", outputDirectory));
  mapCanvas.SaveAs(Form("%s/CDetHalfBarOccupancyMaps.png", outputDirectory));

  TCanvas summaryCanvas("cCDetHalfBarOccupancySummary",
                        "CDet coherent low-occupancy summary", 1800, 850);
  summaryCanvas.Divide(1, 2);
  summaryCanvas.cd(1); gPad->SetLogy(); completeMaximum.Draw("AP");
  TLine completeThreshold(-0.5, coherentLowFraction,
                          kCDetHalfBars - 0.5, coherentLowFraction);
  completeThreshold.SetLineColor(kRed + 1); completeThreshold.SetLineWidth(2);
  completeThreshold.Draw("SAME");
  summaryCanvas.cd(2); gPad->SetLogy(); goodMaximum.Draw("AP");
  TLine goodThreshold(-0.5, coherentLowFraction,
                      kCDetHalfBars - 0.5, coherentLowFraction);
  goodThreshold.SetLineColor(kRed + 1); goodThreshold.SetLineWidth(2);
  goodThreshold.Draw("SAME");
  summaryCanvas.SaveAs(Form("%s/CDetHalfBarOccupancySummary.pdf", outputDirectory));
  summaryCanvas.SaveAs(Form("%s/CDetHalfBarOccupancySummary.png", outputDirectory));

  std::ofstream csv(Form("%s/CDetHalfBarOccupancy.csv", outputDirectory));
  csv << "global_halfbar,layer,halfbar_in_layer,instrumented_pixels,"
         "complete_min,complete_median,complete_max,complete_max_fraction,"
         "good_min,good_median,good_max,good_max_fraction,"
         "coherent_low_complete,coherent_low_good\n";
  for (const HalfBarOccupancy &bar : halfBars)
    csv << bar.global << ',' << bar.layer + 1 << ',' << bar.inLayer << ','
        << kCDetInstrumentedPerHalfBar << ',' << bar.completeMin << ','
        << bar.completeMedian << ',' << bar.completeMax << ','
        << bar.completeMaxFraction << ',' << bar.goodMin << ','
        << bar.goodMedian << ',' << bar.goodMax << ',' << bar.goodMaxFraction
        << ',' << (bar.completeMaxFraction < coherentLowFraction ? 1 : 0)
        << ',' << (bar.goodMaxFraction < coherentLowFraction ? 1 : 0) << '\n';

  TFile output(Form("%s/CDetHalfBarOccupancyStudy.root", outputDirectory),
               "RECREATE");
  hCompleteBoth.Write();
  hGoodBoth.Write();
  completeMaximum.Write();
  goodMaximum.Write();
  output.Close();

  std::cout << "[CDet occupancy] Files/events: " << filesAdded << '/'
            << events << "\n[CDet occupancy] Instrumented detector pixels: "
            << kCDetHalfBars * kCDetInstrumentedPerHalfBar
            << " (14 in each of 168 half-bars)"
            << "\n[CDet occupancy] Layer complete-triplet pixel medians: "
            << layerCompleteMedian[0] << ", " << layerCompleteMedian[1]
            << "\n[CDet occupancy] Layer good-pulse pixel medians: "
            << layerGoodMedian[0] << ", " << layerGoodMedian[1]
            << "\n[CDet occupancy] Coherently low half-bars at maximum-pixel fraction < "
            << coherentLowFraction << ": complete=" << lowCompleteBars
            << ", good=" << lowGoodBars
            << "\n[CDet occupancy] Ignored reference/out-of-range pulse records: "
            << ignoredReferenceOrInvalid
            << "\n[CDet occupancy] Pulse records in geometry-sentinel slots: "
            << pulsesInSentinelSlots
            << "\n[CDet occupancy] Output directory: " << outputDirectory
            << std::endl;
}
