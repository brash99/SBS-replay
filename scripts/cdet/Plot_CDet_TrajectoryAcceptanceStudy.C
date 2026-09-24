#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "CDetGoodPulseConfig.h"
#include "CDetRunDataset.h"

namespace {

constexpr int kAcceptancePixels = 2688;
constexpr int kAcceptancePixelsPerHalfBar = 16;
constexpr int kAcceptanceHalfBars = 168;
constexpr int kAcceptanceHalfBarsPerLayer = 84;
constexpr double kAcceptanceECalZ = 6.144;
constexpr double kAcceptanceXScale = 1.07;
constexpr double kAcceptanceXAlignment = 0.03;
constexpr double kAcceptanceYAlignment = 0.10;

enum class AcceptanceIntersection { kActive, kSuppressed, kGap, kOutside };

const char *IntersectionName(AcceptanceIntersection value) {
  switch (value) {
  case AcceptanceIntersection::kActive: return "active";
  case AcceptanceIntersection::kSuppressed: return "suppressed";
  case AcceptanceIntersection::kGap: return "gap";
  case AcceptanceIntersection::kOutside: return "outside";
  }
  return "unknown";
}

bool ReadGeometryVector(const char *databaseFile, const std::string &key,
                        std::array<double, kAcceptancePixels> &values,
                        bool useLastOccurrence = false) {
  std::ifstream input(databaseFile);
  if (!input)
    return false;
  bool reading = false;
  std::vector<double> current;
  std::vector<double> selected;
  std::string line;
  while (std::getline(input, line)) {
    const size_t comment = line.find('#');
    const std::string content = line.substr(0, comment);
    if (content.find(key) != std::string::npos) {
      if (current.size() >= values.size()) {
        selected.assign(current.begin(), current.begin() + values.size());
        if (!useLastOccurrence)
          break;
      }
      current.clear();
      reading = true;
      continue;
    }
    if (!reading)
      continue;
    std::istringstream tokens(content);
    double value = 0.0;
    while (tokens >> value) {
      current.push_back(value);
      if (current.size() == values.size()) {
        selected = current;
        reading = false;
        if (!useLastOccurrence)
          break;
      }
    }
    if (!useLastOccurrence && selected.size() == values.size())
      break;
  }
  if (current.size() >= values.size())
    selected.assign(current.begin(), current.begin() + values.size());
  if (selected.size() != values.size())
    return false;
  std::copy(selected.begin(), selected.end(), values.begin());
  return true;
}

struct PixelGeometry {
  int id = -1;
  int halfBar = -1;
  int layer = -1;
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  bool suppressed = false;
};

struct LayerGeometry {
  std::vector<PixelGeometry> pixels;
  double z = 0.0;
  double xMin = std::numeric_limits<double>::infinity();
  double xMax = -std::numeric_limits<double>::infinity();
  double yMin = std::numeric_limits<double>::infinity();
  double yMax = -std::numeric_limits<double>::infinity();
};

struct LayerResult {
  AcceptanceIntersection category = AcceptanceIntersection::kOutside;
  int halfBar = -1;
  int pixel = -1;
  double x = 0.0;
  double y = 0.0;
};

LayerResult Classify(const LayerGeometry &geometry, double projectedX,
                     double projectedY, double halfWidthX,
                     double halfLengthY) {
  LayerResult result;
  result.x = projectedX;
  result.y = projectedY;
  double bestDistance = std::numeric_limits<double>::infinity();
  bool foundSuppressed = false;
  for (const PixelGeometry &pixel : geometry.pixels) {
    const double dx = std::fabs(projectedX - pixel.x);
    const double dy = std::fabs(projectedY - pixel.y);
    if (dx > halfWidthX || dy > halfLengthY)
      continue;
    const double distance = (dx / halfWidthX) * (dx / halfWidthX) +
                            (dy / halfLengthY) * (dy / halfLengthY);
    if (!pixel.suppressed) {
      if (result.category != AcceptanceIntersection::kActive || distance < bestDistance) {
        result.category = AcceptanceIntersection::kActive;
        result.halfBar = pixel.halfBar;
        result.pixel = pixel.id;
        bestDistance = distance;
      }
    } else if (result.category != AcceptanceIntersection::kActive &&
               (!foundSuppressed || distance < bestDistance)) {
      foundSuppressed = true;
      result.category = AcceptanceIntersection::kSuppressed;
      result.halfBar = pixel.halfBar;
      result.pixel = pixel.id;
      bestDistance = distance;
    }
  }
  if (result.category == AcceptanceIntersection::kActive ||
      result.category == AcceptanceIntersection::kSuppressed)
    return result;
  const bool insideEnvelope = projectedX >= geometry.xMin &&
      projectedX <= geometry.xMax && projectedY >= geometry.yMin &&
      projectedY <= geometry.yMax;
  if (!insideEnvelope) {
    result.category = AcceptanceIntersection::kOutside;
    return result;
  }

  // Ordinary inter-paddle and inter-module seams are much smaller than the
  // ECal-to-CDet projection resolution. Within the detector envelope, assign
  // a point in such a seam to the nearest physical paddle whose half-bar y
  // extent contains the projection. Preserve only a genuinely uncovered y
  // topology as kGap.
  bestDistance = std::numeric_limits<double>::infinity();
  double bestYDistance = std::numeric_limits<double>::infinity();
  for (const PixelGeometry &pixel : geometry.pixels) {
    const double dy = std::fabs(projectedY - pixel.y);
    if (dy > halfLengthY)
      continue;
    const double dx = std::fabs(projectedX - pixel.x);
    if (dx > bestDistance ||
        (dx == bestDistance && dy >= bestYDistance))
      continue;
    bestDistance = dx;
    bestYDistance = dy;
    result.category = pixel.suppressed ? AcceptanceIntersection::kSuppressed
                                       : AcceptanceIntersection::kActive;
    result.halfBar = pixel.halfBar;
    result.pixel = pixel.id;
  }
  if (!std::isfinite(bestDistance))
    result.category = AcceptanceIntersection::kGap;
  return result;
}

int CategoryIndex(AcceptanceIntersection layer1,
                  AcceptanceIntersection layer2) {
  if (layer1 == AcceptanceIntersection::kActive &&
      layer2 == AcceptanceIntersection::kActive)
    return 0;
  if (layer1 == AcceptanceIntersection::kOutside ||
      layer2 == AcceptanceIntersection::kOutside)
    return 3;
  if (layer1 == AcceptanceIntersection::kGap ||
      layer2 == AcceptanceIntersection::kGap)
    return 2;
  return 1;
}

} // namespace

void Plot_CDet_TrajectoryAcceptanceStudy(
    const char *configFile = "CDet_run6077_projection.conf",
    const char *inputDirectory = nullptr,
    const char *outputDirectory = "CDet_run6077_trajectory_acceptance",
    const char *databaseFile = "../../DB/db_earm.cdet.dat",
    double paddleWidthM = 0.005, double halfBarLengthM = 0.60) {
  CDetGoodPulseConfig::Values config;
  if (!CDetGoodPulseConfig::Load(configFile, config,
                                  "CDet trajectory acceptance"))
    return;
  if (paddleWidthM <= 0.0 || halfBarLengthM <= 0.0) {
    std::cerr << "[CDet acceptance] Invalid physical dimensions.\n";
    return;
  }

  std::array<double, kAcceptancePixels> xRaw{};
  std::array<double, kAcceptancePixels> yRaw{};
  std::array<double, kAcceptancePixels> zRaw{};
  if (!ReadGeometryVector(databaseFile, "earm.cdet.xpos", xRaw) ||
      !ReadGeometryVector(databaseFile, "earm.cdet.ypos", yRaw) ||
      !ReadGeometryVector(databaseFile, "earm.cdet.zpos", zRaw, true)) {
    std::cerr << "[CDet acceptance] Could not read 2688-channel geometry from "
              << databaseFile << ".\n";
    return;
  }

  const std::array<int, 6> suppressedHalfBars{{7, 24, 70, 136, 137, 155}};
  std::array<LayerGeometry, 2> layers;
  std::array<int, 2> instrumented{{0, 0}};
  for (int id = 0; id < kAcceptancePixels; ++id) {
    if (!std::isfinite(xRaw[id]) || !std::isfinite(yRaw[id]) ||
        !std::isfinite(zRaw[id]) || std::fabs(xRaw[id]) >= 900.0 ||
        std::fabs(yRaw[id]) >= 900.0 || std::fabs(zRaw[id]) >= 900.0)
      continue;
    PixelGeometry pixel;
    pixel.id = id;
    pixel.halfBar = id / kAcceptancePixelsPerHalfBar;
    pixel.layer = id / (kAcceptancePixels / 2);
    pixel.x = xRaw[id] * kAcceptanceXScale - kAcceptanceXAlignment;
    pixel.y = yRaw[id];
    pixel.z = zRaw[id];
    pixel.suppressed = std::find(suppressedHalfBars.begin(),
                                 suppressedHalfBars.end(), pixel.halfBar) !=
                       suppressedHalfBars.end();
    LayerGeometry &layer = layers[pixel.layer];
    layer.pixels.push_back(pixel);
    layer.z += pixel.z;
    layer.xMin = std::min(layer.xMin, pixel.x);
    layer.xMax = std::max(layer.xMax, pixel.x);
    layer.yMin = std::min(layer.yMin, pixel.y);
    layer.yMax = std::max(layer.yMax, pixel.y);
    ++instrumented[pixel.layer];
  }
  const double halfWidthX = 0.5 * paddleWidthM * kAcceptanceXScale;
  const double halfLengthY = 0.5 * halfBarLengthM;
  for (int layer = 0; layer < 2; ++layer) {
    if (instrumented[layer] != 1176) {
      std::cerr << "[CDet acceptance] Layer " << layer + 1 << " has "
                << instrumented[layer] << " physical pixels; expected 1176.\n";
      return;
    }
    layers[layer].z /= instrumented[layer];
    layers[layer].xMin -= halfWidthX;
    layers[layer].xMax += halfWidthX;
    layers[layer].yMin -= halfLengthY;
    layers[layer].yMax += halfLengthY;
  }

  TString input(inputDirectory ? inputDirectory : "");
  if (input.IsNull()) {
    const char *outDir = gSystem->Getenv("OUT_DIR");
    if (outDir)
      input = outDir;
  }
  if (input.IsNull()) {
    std::cerr << "[CDet acceptance] An input directory or OUT_DIR is required.\n";
    return;
  }
  TChain chain("T");
  const int filesAdded =
      CDetRunDataset::AddToChain(&chain, config.runNumber, input.Data());
  if (filesAdded <= 0 || chain.GetEntries() <= 0)
    return;
  chain.LoadTree(0);
  const char *required[] = {"earm.ecal.e", "earm.ecal.adctime",
                            "earm.ecal.x", "earm.ecal.y",
                            "earm.cdet.pair.ecal_score"};
  for (const char *branch : required) {
    if (!chain.GetBranch(branch)) {
      std::cerr << "[CDet acceptance] Missing branch " << branch << '\n';
      return;
    }
  }

  TH2D hLayer1Active("hCDetAcceptanceL1Active",
      "Layer 1: active physical pixel;x projected at CDet (m);aligned y projected at CDet (m)",
      180, -1.8, 1.8, 100, -0.7, 0.7);
  TH2D hLayer1Rejected("hCDetAcceptanceL1Rejected",
      "Layer 1: suppressed/gap/outside;x projected at CDet (m);aligned y projected at CDet (m)",
      180, -1.8, 1.8, 100, -0.7, 0.7);
  TH2D hLayer2Active("hCDetAcceptanceL2Active",
      "Layer 2: active physical pixel;x projected at CDet (m);aligned y projected at CDet (m)",
      180, -1.8, 1.8, 100, -0.7, 0.7);
  TH2D hLayer2Rejected("hCDetAcceptanceL2Rejected",
      "Layer 2: suppressed/gap/outside;x projected at CDet (m);aligned y projected at CDet (m)",
      180, -1.8, 1.8, 100, -0.7, 0.7);
  TH1D hEventCategory("hCDetTrajectoryAcceptanceCategory",
      "ECal-admitted trajectory acceptance;Category;Events", 4, -0.5, 3.5);
  hEventCategory.GetXaxis()->SetBinLabel(1, "active both");
  hEventCategory.GetXaxis()->SetBinLabel(2, "suppressed HB");
  hEventCategory.GetXaxis()->SetBinLabel(3, "gap");
  hEventCategory.GetXaxis()->SetBinLabel(4, "outside");
  TH1D hProjectedHalfBar("hCDetProjectedHalfBar",
      "Best physical half-bar intersection;Global half-bar;Layer intersections",
      kAcceptanceHalfBars, -0.5, kAcceptanceHalfBars - 0.5);

  TTreeReader reader(&chain);
  TTreeReaderValue<Double_t> ecalEnergy(reader, "earm.ecal.e");
  TTreeReaderValue<Double_t> ecalTime(reader, "earm.ecal.adctime");
  TTreeReaderValue<Double_t> ecalX(reader, "earm.ecal.x");
  TTreeReaderValue<Double_t> ecalY(reader, "earm.ecal.y");
  TTreeReaderArray<Double_t> pairScore(reader, "earm.cdet.pair.ecal_score");
  std::array<std::array<Long64_t, 4>, 2> layerCounts{};
  std::array<Long64_t, 4> eventCounts{};
  std::array<Long64_t, 4> selectedEventCounts{};
  std::array<Long64_t, kAcceptanceHalfBars> projectedHalfBarCounts{};
  Long64_t admittedEvents = 0;
  while (reader.Next()) {
    if (!std::isfinite(*ecalEnergy) || !std::isfinite(*ecalTime) ||
        !std::isfinite(*ecalX) || !std::isfinite(*ecalY) ||
        *ecalEnergy < config.ecalEnergyMinGeV ||
        *ecalEnergy > config.ecalEnergyMaxGeV ||
        *ecalTime < config.ecalTimeMinNs ||
        *ecalTime > config.ecalTimeMaxNs)
      continue;
    ++admittedEvents;
    std::array<LayerResult, 2> result;
    for (int layer = 0; layer < 2; ++layer) {
      const double projectedX =
          *ecalX * layers[layer].z / kAcceptanceECalZ;
      const double projectedY =
          *ecalY * layers[layer].z / kAcceptanceECalZ +
          kAcceptanceYAlignment;
      result[layer] = Classify(layers[layer], projectedX, projectedY,
                               halfWidthX, halfLengthY);
      const int category = static_cast<int>(result[layer].category);
      ++layerCounts[layer][category];
      if (result[layer].halfBar >= 0) {
        ++projectedHalfBarCounts[result[layer].halfBar];
        hProjectedHalfBar.Fill(result[layer].halfBar);
      }
      TH2D &histogram = layer == 0
          ? (result[layer].category == AcceptanceIntersection::kActive
                 ? hLayer1Active : hLayer1Rejected)
          : (result[layer].category == AcceptanceIntersection::kActive
                 ? hLayer2Active : hLayer2Rejected);
      histogram.Fill(projectedX, projectedY);
    }
    const int eventCategory = CategoryIndex(result[0].category,
                                            result[1].category);
    ++eventCounts[eventCategory];
    if (pairScore.GetSize() > 0)
      ++selectedEventCounts[eventCategory];
    hEventCategory.Fill(eventCategory);
  }

  gSystem->mkdir(outputDirectory, true);
  TCanvas mapCanvas("cCDetTrajectoryAcceptanceMaps",
                    "CDet trajectory acceptance maps", 1800, 1000);
  mapCanvas.Divide(2, 2);
  mapCanvas.cd(1); hLayer1Active.Draw("COLZ");
  mapCanvas.cd(2); hLayer1Rejected.Draw("COLZ");
  mapCanvas.cd(3); hLayer2Active.Draw("COLZ");
  mapCanvas.cd(4); hLayer2Rejected.Draw("COLZ");
  mapCanvas.SaveAs(Form("%s/CDetTrajectoryAcceptanceMaps.pdf", outputDirectory));
  mapCanvas.SaveAs(Form("%s/CDetTrajectoryAcceptanceMaps.png", outputDirectory));

  TCanvas summaryCanvas("cCDetTrajectoryAcceptanceSummary",
                        "CDet trajectory acceptance summary", 1600, 700);
  summaryCanvas.Divide(2, 1);
  summaryCanvas.cd(1); hEventCategory.Draw("HIST TEXT0");
  summaryCanvas.cd(2); hProjectedHalfBar.Draw("HIST");
  summaryCanvas.SaveAs(
      Form("%s/CDetTrajectoryAcceptanceSummary.pdf", outputDirectory));
  summaryCanvas.SaveAs(
      Form("%s/CDetTrajectoryAcceptanceSummary.png", outputDirectory));

  std::ofstream summary(Form("%s/CDetTrajectoryAcceptanceSummary.csv",
                             outputDirectory));
  summary << "scope,category,count,fraction,events_with_stored_pair,"
             "stored_pair_event_fraction\n";
  for (int layer = 0; layer < 2; ++layer) {
    for (int category = 0; category < 4; ++category)
      summary << "layer" << layer + 1 << ','
              << IntersectionName(
                     static_cast<AcceptanceIntersection>(category)) << ','
              << layerCounts[layer][category] << ','
              << (admittedEvents > 0
                      ? double(layerCounts[layer][category]) / admittedEvents
                      : 0.0) << ",,\n";
  }
  const char *eventNames[] = {"active_both", "suppressed_halfbar",
                              "gap", "outside"};
  for (int category = 0; category < 4; ++category)
    summary << "event," << eventNames[category] << ','
            << eventCounts[category] << ','
            << (admittedEvents > 0
                    ? double(eventCounts[category]) / admittedEvents : 0.0)
            << ',' << selectedEventCounts[category] << ','
            << (eventCounts[category] > 0
                    ? double(selectedEventCounts[category]) /
                          eventCounts[category]
                    : 0.0) << '\n';

  std::ofstream bars(Form("%s/CDetProjectedHalfBarCounts.csv",
                          outputDirectory));
  bars << "global_halfbar,layer,halfbar_in_layer,suppressed,projected_intersections\n";
  for (int halfBar = 0; halfBar < kAcceptanceHalfBars; ++halfBar) {
    const bool suppressed = std::find(suppressedHalfBars.begin(),
                                      suppressedHalfBars.end(), halfBar) !=
                            suppressedHalfBars.end();
    bars << halfBar << ',' << halfBar / kAcceptanceHalfBarsPerLayer + 1
         << ',' << halfBar % kAcceptanceHalfBarsPerLayer << ','
         << (suppressed ? 1 : 0) << ','
         << projectedHalfBarCounts[halfBar] << '\n';
  }

  TFile output(Form("%s/CDetTrajectoryAcceptanceStudy.root", outputDirectory),
               "RECREATE");
  hLayer1Active.Write();
  hLayer1Rejected.Write();
  hLayer2Active.Write();
  hLayer2Rejected.Write();
  hEventCategory.Write();
  hProjectedHalfBar.Write();
  output.Close();

  std::cout << "[CDet acceptance] Files/events: " << filesAdded << '/'
            << chain.GetEntries() << "\n[CDet acceptance] ECal-admitted events: "
            << admittedEvents << "\n[CDet acceptance] Layer z positions: "
            << layers[0].z << ", " << layers[1].z << " m"
            << "\n[CDet acceptance] Physical rectangle: "
            << 2.0 * halfWidthX << " m aligned x width x "
            << 2.0 * halfLengthY << " m y length"
            << "\n[CDet acceptance] Event categories: active both="
            << eventCounts[0] << ", suppressed half-bar=" << eventCounts[1]
            << ", gap=" << eventCounts[2] << ", outside=" << eventCounts[3]
            << "\n[CDet acceptance] Active-both fraction: "
            << (admittedEvents > 0 ? double(eventCounts[0]) / admittedEvents
                                   : 0.0)
            << "\n[CDet acceptance] Stored-pair events by category: active both="
            << selectedEventCounts[0] << ", suppressed="
            << selectedEventCounts[1] << ", gap=" << selectedEventCounts[2]
            << ", outside=" << selectedEventCounts[3]
            << "\n[CDet acceptance] Output directory: " << outputDirectory
            << std::endl;
}
