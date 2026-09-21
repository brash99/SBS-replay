#include "PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C"

#include <TEnv.h>
#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {
struct FrozenFitResult {
  bool valid = false;
  long long entries = 0;
  int groups = 0;
  long long ndf = 0;
  double slope = std::numeric_limits<double>::quiet_NaN();
  double error = std::numeric_limits<double>::quiet_NaN();
};

std::string TrimFrozenLine(const std::string &line)
{
  const std::string whitespace = " \t\r\n";
  const size_t first = line.find_first_not_of(whitespace);
  if (first == std::string::npos) return "";
  return line.substr(first, line.find_last_not_of(whitespace) - first + 1);
}

bool ReadFrozenP1(const TString &filename, double &value)
{
  std::ifstream input(filename.Data());
  if (!input) return false;
  std::string line, section;
  while (std::getline(input, line)) {
    const std::string trimmed = TrimFrozenLine(line);
    if (!trimmed.empty() && trimmed.front() == '[') {
      section = trimmed;
      continue;
    }
    std::istringstream parser(trimmed);
    std::string key;
    parser >> key;
    if (section == "[ECalTiming]" && key == "p1") {
      parser >> value;
      return std::isfinite(value);
    }
  }
  return false;
}

FrozenFitResult FitFrozenSamples(
    const std::vector<CDetFrozenECalTimingSample> &samples,
    int firstLayer, int lastLayer, int excludedFold,
    const std::map<size_t, int> &eventFold, double appliedResidual)
{
  FrozenFitResult result;
  std::map<int, std::vector<std::pair<double, double>>> groups;
  for (const auto &sample : samples) {
    if (sample.layer < firstLayer || sample.layer >= lastLayer) continue;
    const auto fold = eventFold.find(sample.eventIndex);
    if (excludedFold >= 0 && fold != eventFold.end() &&
        fold->second == excludedFold) continue;
    groups[sample.halfBar].push_back(
        {sample.ecalTime,
         sample.cdetTime - appliedResidual*sample.ecalTime});
  }

  double withinXY = 0.0, withinXX = 0.0;
  std::vector<std::pair<int, std::pair<double, double>>> means;
  for (const auto &group : groups) {
    if (group.second.size() < 3) continue;
    double sumX = 0.0, sumY = 0.0;
    for (const auto &sample : group.second) {
      sumX += sample.first;
      sumY += sample.second;
    }
    const double meanX = sumX/group.second.size();
    const double meanY = sumY/group.second.size();
    means.push_back({group.first, {meanX, meanY}});
    result.entries += group.second.size();
    for (const auto &sample : group.second) {
      withinXX += std::pow(sample.first - meanX, 2);
      withinXY += (sample.first - meanX)*(sample.second - meanY);
    }
  }
  result.groups = means.size();
  result.ndf = result.entries - result.groups - 1;
  if (!(withinXX > 0.0) || result.ndf <= 0) return result;
  result.slope = withinXY/withinXX;

  double residualSquares = 0.0;
  for (const auto &mean : means) {
    const auto group = groups.find(mean.first);
    if (group == groups.end()) continue;
    for (const auto &sample : group->second) {
      const double dx = sample.first - mean.second.first;
      const double dy = sample.second - mean.second.second;
      residualSquares += std::pow(dy - result.slope*dx, 2);
    }
  }
  result.error = std::sqrt((residualSquares/result.ndf)/withinXX);
  result.valid = std::isfinite(result.slope) &&
                 std::isfinite(result.error) && result.error > 0.0;
  return result;
}

FrozenFitResult FitHeldOutFold(
    const std::vector<CDetFrozenECalTimingSample> &samples, int selectedFold,
    const std::map<size_t, int> &eventFold, double appliedResidual)
{
  std::vector<CDetFrozenECalTimingSample> heldOut;
  for (const auto &sample : samples) {
    const auto fold = eventFold.find(sample.eventIndex);
    if (fold != eventFold.end() && fold->second == selectedFold)
      heldOut.push_back(sample);
  }
  const std::map<size_t, int> noExclusions;
  return FitFrozenSamples(heldOut, 0, 2, -1, noExclusions, appliedResidual);
}

double Quantile(std::vector<double> values, double fraction)
{
  if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
  std::sort(values.begin(), values.end());
  const double position = fraction*(values.size() - 1);
  const size_t lower = static_cast<size_t>(std::floor(position));
  const size_t upper = static_cast<size_t>(std::ceil(position));
  const double weight = position - lower;
  return values[lower]*(1.0 - weight) + values[upper]*weight;
}
} // namespace

// Read-only Run-specific p1 qualification. This freezes the accepted sample
// once, derives the candidate algebraically, and validates it on contiguous
// event folds. It never writes CDet_run<run>.dat.
void Run_CDet_Qualify_RunSpecificP1_FrozenSample(
    Int_t runNumber,
    Int_t nevents = std::numeric_limits<Int_t>::min(),
    TString configFile = "",
    Int_t folds = 5,
    double minimumSignificance = 3.0,
    double minimumTimingSpan = 0.5,
    double maximumValidationPull = 2.0,
    double minimumValidationImprovement = 0.5,
    double maximumResolutionRatio = 1.05)
{
  if (runNumber <= 0 || folds < 3 || minimumSignificance <= 0.0 ||
      minimumTimingSpan < 0.0 || maximumValidationPull <= 0.0 ||
      minimumValidationImprovement < 0.0 ||
      minimumValidationImprovement >= 1.0 || maximumResolutionRatio < 1.0) {
    std::cerr << "[Frozen p1] ERROR: invalid arguments.\n";
    return;
  }
  if (configFile.IsNull())
    configFile = TString::Format("CDet_run%d_projection.conf", runNumber);
  TEnv env;
  if (!LoadCDetConfiguration(env, configFile.Data(), "Frozen p1")) return;
  if (env.GetValue("analysis.run_number", -1) != runNumber) {
    std::cerr << "[Frozen p1] ERROR: configuration run number does not match.\n";
    return;
  }

  if (gSystem->AccessPathName("CDet_calibration_dt.dat")) {
    std::cerr << "[Frozen p1] ERROR: CDet_calibration_dt.dat is missing.\n";
    return;
  }
  double masterFileP1 = std::numeric_limits<double>::quiet_NaN();
  if (!ReadFrozenP1("CDet_calibration_dt.dat", masterFileP1)) {
    std::cerr << "[Frozen p1] ERROR: master p1 is missing or invalid.\n";
    return;
  }
  const TString runFile = TString::Format("CDet_run%d.dat", runNumber);
  double runFileP1 = std::numeric_limits<double>::quiet_NaN();
  if (ReadFrozenP1(runFile, runFileP1) &&
      std::fabs(runFileP1 - masterFileP1) > 5.0e-7) {
    std::cerr << "[Frozen p1] ERROR: " << runFile << " contains p1="
              << runFileP1 << ", which differs from master p1="
              << masterFileP1 << ". Restore the master-p1 state before "
                 "qualification.\n";
    return;
  }

  ResetCalibrationGlobals();
  PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget(
      configFile.Data(), 7, nevents);
  if (!gLastCalibrationStageSucceeded) {
    std::cerr << "[Frozen p1] ERROR: Stage-7 analysis failed.\n";
    return;
  }
  const double masterP1 = masterFileP1;
  plotCDetLayersTimeComp(configFile.Data(), 0);
  if (!gLastCalibrationFitSucceeded || !gLastECalFixedEffectsValid ||
      gCDetFrozenECalTimingSamples.empty()) {
    std::cerr << "[Frozen p1] ERROR: frozen timing samples are unavailable.\n";
    return;
  }

  const std::vector<CDetFrozenECalTimingSample> frozen =
      gCDetFrozenECalTimingSamples;
  std::set<size_t> eventSet;
  std::vector<double> ecalTimes;
  for (const auto &sample : frozen) {
    eventSet.insert(sample.eventIndex);
    ecalTimes.push_back(sample.ecalTime);
  }
  if (eventSet.size() < static_cast<size_t>(folds)) {
    std::cerr << "[Frozen p1] ERROR: fewer accepted events than folds.\n";
    return;
  }
  std::map<size_t, int> eventFold;
  size_t eventRank = 0;
  for (size_t event : eventSet) {
    int fold = static_cast<int>((eventRank*folds)/eventSet.size());
    if (fold >= folds) fold = folds - 1;
    eventFold[event] = fold;
    ++eventRank;
  }

  const FrozenFitResult full = FitFrozenSamples(
      frozen, 0, 2, -1, eventFold, 0.0);
  const FrozenFitResult layer1 = FitFrozenSamples(
      frozen, 0, 1, -1, eventFold, 0.0);
  const FrozenFitResult layer2 = FitFrozenSamples(
      frozen, 1, 2, -1, eventFold, 0.0);
  const FrozenFitResult closure = FitFrozenSamples(
      frozen, 0, 2, -1, eventFold, full.slope);
  if (!full.valid || !layer1.valid || !layer2.valid || !closure.valid) {
    std::cerr << "[Frozen p1] ERROR: a required fixed-effects fit failed.\n";
    return;
  }

  const double candidateP1 = masterP1 + full.slope;
  const double significance = std::fabs(full.slope)/full.error;
  const double layerDifferenceError = std::hypot(layer1.error, layer2.error);
  const double layerDifferencePull =
      std::fabs(layer1.slope - layer2.slope)/layerDifferenceError;
  const bool layerSignAgreement = layer1.slope*layer2.slope > 0.0;
  const double centralECalSpan = Quantile(ecalTimes, 0.90) - Quantile(ecalTimes, 0.10);
  const double timingEffect = std::fabs(full.slope)*centralECalSpan;

  std::map<std::pair<size_t, size_t>, std::vector<const CDetFrozenECalTimingSample *>> pairs;
  for (const auto &sample : frozen)
    pairs[{sample.eventIndex, sample.pairIndex}].push_back(&sample);
  std::vector<double> pairBefore, pairAfter;
  for (const auto &pair : pairs) {
    if (pair.second.size() != 2) continue;
    const double before = 0.5*(pair.second[0]->cdetTime + pair.second[1]->cdetTime);
    const double ecal = pair.second[0]->ecalTime;
    pairBefore.push_back(before);
    pairAfter.push_back(before - full.slope*ecal);
  }
  auto rms = [](const std::vector<double> &values) {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    double mean = 0.0;
    for (double value : values) mean += value;
    mean /= values.size();
    double variance = 0.0;
    for (double value : values) variance += std::pow(value - mean, 2);
    return std::sqrt(variance/values.size());
  };
  const double rmsBefore = rms(pairBefore);
  const double rmsAfter = rms(pairAfter);
  const double resolutionRatio = rmsAfter/rmsBefore;

  const TString sampleFile = TString::Format(
      "CDet_run%d_frozen_p1_samples.tsv", runNumber);
  const TString summaryFile = TString::Format(
      "CDet_run%d_frozen_p1_qualification.txt", runNumber);
  std::ofstream sampleOutput(sampleFile.Data());
  std::ofstream summary(summaryFile.Data());
  if (!sampleOutput || !summary) {
    std::cerr << "[Frozen p1] ERROR: cannot create output artifacts.\n";
    return;
  }
  sampleOutput << "event_index\tpair_index\tfold\tlayer\tpixel\thalf_bar"
                  "\tecal_time_ns\tcdet_time_master_ns\tcdet_time_candidate_ns\ttot_ns\n";
  sampleOutput << std::setprecision(12);
  for (const auto &sample : frozen) {
    sampleOutput << sample.eventIndex << "\t" << sample.pairIndex << "\t"
                 << eventFold[sample.eventIndex] << "\t" << sample.layer + 1
                 << "\t" << sample.pixel << "\t" << sample.halfBar << "\t"
                 << sample.ecalTime << "\t" << sample.cdetTime << "\t"
                 << sample.cdetTime - full.slope*sample.ecalTime << "\t"
                 << sample.tot << "\n";
  }

  bool foldsPass = true;
  summary << std::setprecision(12)
          << "run " << runNumber << "\nconfig " << configFile << "\n"
          << "master_p1 " << masterP1 << "\nresidual_p1 " << full.slope
          << "\nresidual_p1_error " << full.error << "\ncandidate_p1 "
          << candidateP1 << "\nfull_significance " << significance
          << "\nlayer1_residual " << layer1.slope << "\nlayer1_error "
          << layer1.error << "\nlayer2_residual " << layer2.slope
          << "\nlayer2_error " << layer2.error << "\nlayer_difference_pull "
          << layerDifferencePull << "\ncentral_ecal_span_ns " << centralECalSpan
          << "\ntiming_effect_ns " << timingEffect << "\nfull_closure_residual "
          << closure.slope << "\nfull_closure_error " << closure.error
          << "\npair_rms_before_ns " << rmsBefore << "\npair_rms_after_ns "
          << rmsAfter << "\nresolution_ratio " << resolutionRatio << "\n";

  std::cout << "\n[Frozen p1 qualification]\n"
            << "  frozen hits/events/pairs: " << frozen.size() << "/"
            << eventSet.size() << "/" << pairBefore.size() << "\n"
            << "  master p1: " << masterP1 << " ns/ns\n"
            << "  residual: " << full.slope << " +/- " << full.error
            << " ns/ns (" << significance << " sigma)\n"
            << "  one-step candidate: " << candidateP1 << " ns/ns\n"
            << "  algebraic full-sample closure: " << closure.slope
            << " +/- " << closure.error << " ns/ns\n"
            << "  central-80% ECal span/effect: " << centralECalSpan << "/"
            << timingEffect << " ns\n"
            << "  pair RMS before/after: " << rmsBefore << "/" << rmsAfter
            << " ns\n";

  for (int fold = 0; fold < folds; ++fold) {
    const FrozenFitResult training = FitFrozenSamples(
        frozen, 0, 2, fold, eventFold, 0.0);
    const FrozenFitResult validationBefore = FitHeldOutFold(
        frozen, fold, eventFold, 0.0);
    const FrozenFitResult validationAfter = FitHeldOutFold(
        frozen, fold, eventFold, training.slope);
    const double pull = validationAfter.valid
        ? std::fabs(validationAfter.slope)/validationAfter.error
        : std::numeric_limits<double>::infinity();
    const double improvement = validationBefore.valid &&
        std::fabs(validationBefore.slope) > 0.0
        ? 1.0 - std::fabs(validationAfter.slope)/std::fabs(validationBefore.slope)
        : -std::numeric_limits<double>::infinity();
    const double candidateDeviation = training.valid
        ? std::fabs(training.slope - full.slope)
        : std::numeric_limits<double>::infinity();
    const double stabilityTolerance = std::max(0.03, 2.0*full.error);
    const bool pass = training.valid && validationBefore.valid &&
        validationAfter.valid && pull <= maximumValidationPull &&
        improvement >= minimumValidationImprovement &&
        candidateDeviation <= stabilityTolerance;
    foldsPass = foldsPass && pass;
    summary << "fold " << fold << " training_residual " << training.slope
            << " validation_before " << validationBefore.slope
            << " validation_before_error " << validationBefore.error
            << " validation_after " << validationAfter.slope
            << " validation_after_error " << validationAfter.error
            << " validation_pull " << pull << " improvement " << improvement
            << " candidate_deviation " << candidateDeviation
            << " pass " << pass << "\n";
    std::cout << "  fold " << fold << ": train=" << training.slope
              << ", held-out " << validationBefore.slope << " -> "
              << validationAfter.slope << " +/- " << validationAfter.error
              << ", pull=" << pull << ", improvement=" << improvement
              << ", " << (pass ? "PASS" : "FAIL") << "\n";
  }

  const bool discoveryPass = significance >= minimumSignificance;
  const bool layersPass = layerSignAgreement && layerDifferencePull <= 2.0;
  const bool practicalPass = timingEffect >= minimumTimingSpan;
  const bool closurePass = std::fabs(closure.slope) <= 1.0e-10;
  const bool resolutionPass = std::isfinite(resolutionRatio) &&
                              resolutionRatio <= maximumResolutionRatio;
  const bool qualified = discoveryPass && layersPass && practicalPass &&
                         closurePass && foldsPass && resolutionPass;
  summary << "gate_discovery " << discoveryPass << "\ngate_layers "
          << layersPass << "\ngate_practical " << practicalPass
          << "\ngate_closure " << closurePass << "\ngate_folds " << foldsPass
          << "\ngate_resolution " << resolutionPass << "\nqualified "
          << qualified << "\n";
  summary.close();
  sampleOutput.close();

  std::cout << "  gates: discovery=" << discoveryPass
            << " layers=" << layersPass << " practical=" << practicalPass
            << " closure=" << closurePass << " folds=" << foldsPass
            << " resolution=" << resolutionPass << "\n"
            << "  QUALIFIED: " << (qualified ? "YES" : "NO") << "\n"
            << "  sample artifact: " << sampleFile << "\n"
            << "  summary: " << summaryFile << "\n"
            << "  CDet_run" << runNumber << ".dat was not modified.\n";
}
