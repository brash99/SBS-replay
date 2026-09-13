#include "PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C"

#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
std::string TrimRunShiftLine(const std::string &line)
{
  const std::string whitespace = " \t\r\n";
  const std::string::size_type first = line.find_first_not_of(whitespace);
  if (first == std::string::npos) return "";
  const std::string::size_type last = line.find_last_not_of(whitespace);
  return line.substr(first, last - first + 1);
}

bool WriteRunTimingShift(const TString &filename, Int_t runNumber, double shift)
{
  std::vector<std::string> lines;
  if (!gSystem->AccessPathName(filename)) {
    std::ifstream input(filename.Data());
    if (!input) {
      std::cerr << "[Run timing-shift calibration] ERROR: cannot read "
                << filename << ".\n";
      return false;
    }
    std::string line;
    while (std::getline(input, line)) lines.push_back(line);
  } else {
    lines.push_back(TString::Format(
        "# Run-specific timing calibration for cross-target run %d.",
        runNumber).Data());
  }

  std::ostringstream value;
  value << "shift_ns " << std::fixed << std::setprecision(6) << shift;
  bool inGlobal = false;
  bool globalSeen = false;
  bool shiftWritten = false;
  std::vector<std::string> updated;
  for (const std::string &line : lines) {
    const std::string trimmed = TrimRunShiftLine(line);
    const bool isSection = trimmed.size() >= 2 && trimmed.front() == '[' &&
                           trimmed.back() == ']';
    if (isSection) {
      if (inGlobal && !shiftWritten) {
        updated.push_back(value.str());
        shiftWritten = true;
      }
      inGlobal = trimmed == "[GlobalTiming]";
      if (inGlobal) globalSeen = true;
      updated.push_back(line);
      continue;
    }
    std::istringstream parser(trimmed);
    std::string key;
    parser >> key;
    if (inGlobal && key == "shift_ns") {
      if (!shiftWritten) updated.push_back(value.str());
      shiftWritten = true;
    } else {
      updated.push_back(line);
    }
  }
  if (inGlobal && !shiftWritten) {
    updated.push_back(value.str());
    shiftWritten = true;
  }
  if (!globalSeen) {
    if (!updated.empty() && !updated.back().empty()) updated.push_back("");
    updated.push_back("[GlobalTiming]");
    updated.push_back(value.str());
  }

  const TString temporary = filename + ".tmp";
  std::ofstream output(temporary.Data());
  if (!output) {
    std::cerr << "[Run timing-shift calibration] ERROR: cannot write "
              << temporary << ".\n";
    return false;
  }
  for (const std::string &line : updated) output << line << "\n";
  output.flush();
  if (!output) {
    output.close();
    std::remove(temporary.Data());
    std::cerr << "[Run timing-shift calibration] ERROR: failed while writing "
              << temporary << ".\n";
    return false;
  }
  output.close();
  if (std::rename(temporary.Data(), filename.Data()) != 0) {
    std::remove(temporary.Data());
    std::cerr << "[Run timing-shift calibration] ERROR: cannot replace "
              << filename << ".\n";
    return false;
  }
  return true;
}

bool AnalyzeRunForTimingShift(const TString &configFile, Int_t nevents,
                              double gaussianFitHalfWidth,
                              bool useProjectedHalfBarPairs)
{
  ResetCalibrationGlobals();
  PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget(
      configFile.Data(), 7, nevents);
  if (!gLastCalibrationStageSucceeded) return false;
  gShiftInvariantPairTimingCuts = true;
  plotCDetLayersTimeComp(configFile.Data(), 0);
  const std::vector<double> &timingSample = useProjectedHalfBarPairs
      ? gCDetProjectedHalfBarPairMeanTimes : gCDetAcceptedPairMeanTimes;
  if (!gCDetPairedMeanTimeVsECal || timingSample.empty()) {
    gShiftInvariantPairTimingCuts = false;
    return false;
  }
  double sampleMean = 0.0;
  for (double value : timingSample) sampleMean += value;
  sampleMean /= static_cast<double>(timingSample.size());
  if (!std::isfinite(sampleMean)) {
    gShiftInvariantPairTimingCuts = false;
    return false;
  }
  if (useProjectedHalfBarPairs) {
    // The LH2 projected sample has a narrow physical peak followed by a long
    // upper-time tail. Locate the peak on a translation-invariant grid, then
    // fit locally; the ordinary sample mean is not a suitable fit center.
    const double binWidth = 0.5;
    double maxDeviation = 0.0;
    for (double value : timingSample)
      maxDeviation = std::max(maxDeviation, std::fabs(value - sampleMean));
    const int halfBins = std::max(20, static_cast<int>(
        std::ceil(maxDeviation/binWidth)) + 2);
    TH1D peakHistogram("hCDetProjectedHalfBarShiftFit",
        "Projected-half-bar pair mean time;Corrected pair mean time (ns);Pairs",
        2*halfBins, sampleMean - halfBins*binWidth,
        sampleMean + halfBins*binWidth);
    peakHistogram.SetDirectory(nullptr);
    for (double value : timingSample) peakHistogram.Fill(value);
    const double peakCenter = peakHistogram.GetBinCenter(
        peakHistogram.GetMaximumBin());
    TF1 coreFit("fCDetProjectedHalfBarShiftCore", "gaus",
                peakCenter - gaussianFitHalfWidth,
                peakCenter + gaussianFitHalfWidth);
    const TFitResultPtr fitResult = peakHistogram.Fit(&coreFit, "QRSN");
    gLastPairedTimeEntries = static_cast<double>(timingSample.size());
    gLastPairedCoreFitValid = int(fitResult) == 0;
    if (gLastPairedCoreFitValid) {
      gLastPairedCoreMean = coreFit.GetParameter(1);
      gLastPairedCoreMeanError = coreFit.GetParError(1);
      gLastPairedCoreSigma = coreFit.GetParameter(2);
      gLastPairedCoreSigmaError = coreFit.GetParError(2);
      gLastPairedCoreFitValid =
          std::isfinite(gLastPairedCoreMean) &&
          std::isfinite(gLastPairedCoreMeanError) &&
          std::isfinite(gLastPairedCoreSigma) &&
          gLastPairedCoreSigma > 0.0 &&
          gLastPairedCoreSigma <= gaussianFitHalfWidth &&
          std::fabs(gLastPairedCoreMean - peakCenter) <= gaussianFitHalfWidth;
    }
    std::cout << "\n[Run timing-shift projected-half-bar fit]\n"
              << "  entries: " << timingSample.size() << "\n"
              << "  modal-bin center: " << peakCenter << " ns\n"
              << "  Gaussian core mean: " << gLastPairedCoreMean << " +/- "
              << gLastPairedCoreMeanError << " ns\n"
              << "  Gaussian core sigma: " << gLastPairedCoreSigma << " +/- "
              << gLastPairedCoreSigmaError << " ns\n"
              << "  physical-fit gate: "
              << (gLastPairedCoreFitValid ? "PASS" : "FAIL") << "\n";
  } else {
    reportCDetPairedTimeResolution(
        false, sampleMean - gaussianFitHalfWidth,
        sampleMean + gaussianFitHalfWidth);
  }
  const bool analysisSucceeded =
      gLastCalibrationFitSucceeded && gLastECalFixedEffectsValid &&
      gLastPairedCoreFitValid;
  gShiftInvariantPairTimingCuts = false;
  return analysisSucceeded;
}
} // namespace

// Set the run-specific final timing origin from the Gaussian core centroid of
// either all accepted pairs or, when explicitly requested, the ECal-projected
// half-bar pairs. Existing run-file keys (notably a fixed ECal p1) are
// preserved. The projected estimator locates the narrow modal peak before its
// local fit, which avoids bias from the long upper-time tail in LH2 data.
void Run_CDet_Calibrate_RunTimingShift(
    Int_t runNumber,
    Int_t nevents = std::numeric_limits<Int_t>::min(),
    double targetMean = 30.0,
    double gaussianFitHalfWidth = 5.0,
    TString configFile = "",
    bool useProjectedHalfBarPairs = false)
{
  gLastCalibrationSequenceSucceeded = false;
  if (runNumber <= 0 || !std::isfinite(targetMean) ||
      !(gaussianFitHalfWidth > 0.0)) {
    std::cerr << "[Run timing-shift calibration] ERROR: invalid arguments.\n";
    return;
  }
  if (configFile.IsNull())
    configFile = TString::Format("CDet_run%d_projection.conf", runNumber);
  TEnv env;
  if (!LoadCDetConfiguration(env, configFile.Data(),
                             "Run timing-shift calibration")) return;
  if (env.GetValue("analysis.run_number", -1) != runNumber) {
    std::cerr << "[Run timing-shift calibration] ERROR: analysis.run_number in "
              << configFile << " does not match " << runNumber << ".\n";
    return;
  }
  if (!env.Defined("analysis.ecal_time_min") ||
      !env.Defined("analysis.ecal_time_max") ||
      env.GetValue("analysis.ecal_time_min", 0.0) >=
          env.GetValue("analysis.ecal_time_max", 0.0)) {
    std::cerr << "[Run timing-shift calibration] ERROR: " << configFile
              << " must explicitly define a valid analysis ECal-time window.\n";
    return;
  }
  if (gSystem->AccessPathName("CDet_calibration_dt.dat")) {
    std::cerr << "[Run timing-shift calibration] ERROR: "
                 "CDet_calibration_dt.dat is missing.\n";
    return;
  }

  const TString runFile = TString::Format("CDet_run%d.dat", runNumber);
  LoadRunTimingConstants(runFile.Data());
  if (!AnalyzeRunForTimingShift(configFile, nevents,
                                gaussianFitHalfWidth,
                                useProjectedHalfBarPairs)) {
    std::cerr << "[Run timing-shift calibration] ERROR: baseline analysis or "
                 "Gaussian-core fit failed.\n";
    return;
  }
  const double baselineCentroid = gLastPairedCoreMean;
  const double baselineCentroidError = gLastPairedCoreMeanError;
  const double baselinePairEntries = gLastPairedTimeEntries;
  const double previousShift = gGlobalTimingLoaded ? gGlobalTimingShift : 0.0;
  const double calibratedShift = previousShift + targetMean - baselineCentroid;
  if (!std::isfinite(calibratedShift) ||
      !WriteRunTimingShift(runFile, runNumber, calibratedShift)) return;

  std::cout << "\n[Run timing-shift calibration]\n"
            << "  estimator: "
            << (useProjectedHalfBarPairs ? "projected-half-bar pairs"
                                         : "all accepted pairs") << "\n"
            << "  baseline Gaussian-core centroid: " << baselineCentroid
            << " +/- " << baselineCentroidError << " ns\n"
            << "  previous shift_ns: " << previousShift << " ns\n"
            << "  target centroid: " << targetMean << " ns\n"
            << "  stored shift_ns: " << calibratedShift << " ns\n"
            << "  output: " << runFile << "\n";

  if (!AnalyzeRunForTimingShift(configFile, nevents,
                                gaussianFitHalfWidth,
                                useProjectedHalfBarPairs)) {
    std::cerr << "[Run timing-shift calibration] ERROR: closure analysis failed.\n";
    return;
  }
  double closureTolerance =
      std::max(0.02, 3.0 * gLastPairedCoreMeanError);
  double closureResidual = gLastPairedCoreMean - targetMean;
  std::cout << "  closure Gaussian-core centroid: " << gLastPairedCoreMean
            << " +/- " << gLastPairedCoreMeanError << " ns\n"
            << "  closure residual: " << closureResidual << " ns\n"
            << "  accepted-pair entries before/after: "
            << baselinePairEntries << "/" << gLastPairedTimeEntries << "\n"
            << "  ECal-slope closure: " << gLastECalFixedEffectsSlope
            << " +/- " << gLastECalFixedEffectsSlopeError << " ns/ns\n";
  if (gLastPairedTimeEntries != baselinePairEntries) {
    std::cerr << "[Run timing-shift calibration] ERROR: the additive timing "
                 "shift changed the accepted-pair population.\n";
    return;
  }
  if (std::fabs(closureResidual) > closureTolerance) {
    const double refinedShift = calibratedShift - closureResidual;
    std::cout << "  first closure exceeds " << closureTolerance
              << " ns; refining shift_ns to " << refinedShift << " ns.\n";
    if (!WriteRunTimingShift(runFile, runNumber, refinedShift) ||
        !AnalyzeRunForTimingShift(configFile, nevents,
                                  gaussianFitHalfWidth,
                                  useProjectedHalfBarPairs)) {
      std::cerr << "[Run timing-shift calibration] ERROR: refined closure "
                   "analysis failed.\n";
      return;
    }
    closureTolerance = std::max(0.02, 3.0 * gLastPairedCoreMeanError);
    closureResidual = gLastPairedCoreMean - targetMean;
    std::cout << "  refined closure Gaussian-core centroid: "
              << gLastPairedCoreMean << " +/- "
              << gLastPairedCoreMeanError << " ns\n"
              << "  refined closure residual: " << closureResidual << " ns\n"
              << "  refined accepted-pair entries before/after: "
              << baselinePairEntries << "/" << gLastPairedTimeEntries << "\n"
              << "  refined ECal-slope closure: "
              << gLastECalFixedEffectsSlope << " +/- "
              << gLastECalFixedEffectsSlopeError << " ns/ns\n";
    if (gLastPairedTimeEntries != baselinePairEntries ||
        std::fabs(closureResidual) > closureTolerance) {
      std::cerr << "[Run timing-shift calibration] ERROR: refined closure "
                   "did not satisfy the population and centroid gates.\n";
      return;
    }
  }
  gLastCalibrationSequenceSucceeded = true;
}
