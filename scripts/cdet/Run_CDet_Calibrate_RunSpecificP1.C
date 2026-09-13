#include "PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget.C"

#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
std::string TrimRunP1Line(const std::string &line)
{
  const std::string whitespace = " \t\r\n";
  const size_t first = line.find_first_not_of(whitespace);
  if (first == std::string::npos) return "";
  const size_t last = line.find_last_not_of(whitespace);
  return line.substr(first, last - first + 1);
}

bool RunFileHasP1(const TString &filename)
{
  std::ifstream input(filename.Data());
  std::string line, section;
  while (std::getline(input, line)) {
    const std::string trimmed = TrimRunP1Line(line);
    if (!trimmed.empty() && trimmed.front() == '[') {
      section = trimmed;
      continue;
    }
    std::istringstream parser(trimmed);
    std::string key;
    parser >> key;
    if (section == "[ECalTiming]" && key == "p1") return true;
  }
  return false;
}

bool WriteRunSpecificP1(const TString &filename, Int_t runNumber, double p1)
{
  std::vector<std::string> lines;
  if (!gSystem->AccessPathName(filename)) {
    std::ifstream input(filename.Data());
    if (!input) return false;
    std::string line;
    while (std::getline(input, line)) lines.push_back(line);
  } else {
    lines.push_back(TString::Format(
        "# Run-specific timing calibration for cross-target run %d.", runNumber).Data());
  }
  std::ostringstream value;
  value << "p1 " << std::fixed << std::setprecision(6) << p1;
  bool inTiming = false, sectionSeen = false, valueWritten = false;
  std::vector<std::string> updated;
  for (const std::string &line : lines) {
    const std::string trimmed = TrimRunP1Line(line);
    const bool isSection = trimmed.size() >= 2 && trimmed.front() == '[' &&
                           trimmed.back() == ']';
    if (isSection) {
      if (inTiming && !valueWritten) {
        updated.push_back(value.str());
        valueWritten = true;
      }
      inTiming = trimmed == "[ECalTiming]";
      if (inTiming) sectionSeen = true;
      updated.push_back(line);
      continue;
    }
    std::istringstream parser(trimmed);
    std::string key;
    parser >> key;
    if (inTiming && key == "p1") {
      if (!valueWritten) updated.push_back(value.str());
      valueWritten = true;
    } else {
      updated.push_back(line);
    }
  }
  if (inTiming && !valueWritten) {
    updated.push_back(value.str());
    valueWritten = true;
  }
  if (!sectionSeen) {
    if (!updated.empty() && !updated.back().empty()) updated.push_back("");
    updated.push_back("[ECalTiming]");
    updated.push_back(value.str());
  }

  const TString temporary = filename + ".tmp";
  std::ofstream output(temporary.Data());
  if (!output) {
    std::cerr << "[Run p1 calibration] ERROR: cannot write " << temporary << ".\n";
    return false;
  }
  for (const std::string &line : updated) output << line << "\n";
  output.flush();
  if (!output) {
    output.close();
    std::remove(temporary.Data());
    std::cerr << "[Run p1 calibration] ERROR: failed while writing "
              << temporary << ".\n";
    return false;
  }
  output.close();
  if (std::rename(temporary.Data(), filename.Data()) != 0) {
    std::remove(temporary.Data());
    std::cerr << "[Run p1 calibration] ERROR: cannot replace " << filename << ".\n";
    return false;
  }
  return true;
}

bool AnalyzeRunForP1(const TString &configFile, Int_t nevents)
{
  ResetCalibrationGlobals();
  PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget(
      configFile.Data(), 7, nevents);
  if (!gLastCalibrationStageSucceeded) return false;
  plotCDetLayersTimeComp(configFile.Data(), 0);
  return gLastCalibrationFitSucceeded && gLastECalFixedEffectsValid;
}
} // namespace

// Fit only the run-specific ECal timing slope.  The input master calibration
// is read but never overwritten. The approved [ECalSelection] section is
// mandatory and is preserved when p1 is added or replaced.
void Run_CDet_Calibrate_RunSpecificP1(
    Int_t runNumber,
    Int_t nevents = std::numeric_limits<Int_t>::min(),
    bool overwriteExistingRunFile = false,
    TString configFile = "")
{
  gLastCalibrationSequenceSucceeded = false;
  if (runNumber <= 0) {
    std::cerr << "[Run p1 calibration] ERROR: runNumber must be positive.\n";
    return;
  }
  if (configFile.IsNull())
    configFile = TString::Format("CDet_run%d_projection.conf", runNumber);
  TEnv env;
  if (!LoadCDetConfiguration(env, configFile.Data(),
                             "Run p1 calibration")) return;
  if (env.GetValue("analysis.run_number", -1) != runNumber) {
    std::cerr << "[Run p1 calibration] ERROR: analysis.run_number in "
              << configFile << " does not match " << runNumber << ".\n";
    return;
  }
  if (!env.Defined("analysis.ecal_time_min") ||
      !env.Defined("analysis.ecal_time_max") ||
      env.GetValue("analysis.ecal_time_min", 0.0) >=
          env.GetValue("analysis.ecal_time_max", 0.0)) {
    std::cerr << "[Run p1 calibration] ERROR: " << configFile
              << " must explicitly define a valid analysis ECal-time window.\n";
    return;
  }
  if (gSystem->AccessPathName("CDet_calibration_dt.dat")) {
    std::cerr << "[Run p1 calibration] ERROR: CDet_calibration_dt.dat is missing.\n";
    return;
  }

  const TString runFile = TString::Format("CDet_run%d.dat", runNumber);
  LoadRunTimingConstants(runFile.Data());
  if (RunFileHasP1(runFile) && !overwriteExistingRunFile) {
    std::cerr << "[Run p1 calibration] ERROR: " << runFile
              << " already contains p1. Pass overwriteExistingRunFile=true "
                 "to replace that key while preserving the rest of the file.\n";
    return;
  }

  // With no p1 key present, Stage 7 loads the detector-wide master p1.
  if (!AnalyzeRunForP1(configFile, nevents)) {
    std::cerr << "[Run p1 calibration] ERROR: baseline fixed-effects fit failed.\n";
    return;
  }
  const double masterP1 = gECalFitP1;
  const double residualP1 = gLastECalFixedEffectsSlope;
  const double residualError = gLastECalFixedEffectsSlopeError;
  const double runP1 = masterP1 + residualP1;
  if (!std::isfinite(runP1)) {
    std::cerr << "[Run p1 calibration] ERROR: candidate p1 is not finite.\n";
    return;
  }
  if (!WriteRunSpecificP1(runFile, runNumber, runP1)) return;

  std::cout << "\n[Run p1 calibration]\n"
            << "  detector master p1: " << masterP1 << " ns/ns\n"
            << "  run residual p1: " << residualP1 << " +/- "
            << residualError << " ns/ns\n"
            << "  stored run p1: " << runP1 << " ns/ns\n"
            << "  output: " << runFile << "\n";

  // Reload with the newly written override and require fixed-effects closure.
  if (!AnalyzeRunForP1(configFile, nevents)) {
    std::cerr << "[Run p1 calibration] ERROR: closure analysis failed.\n";
    return;
  }
  std::cout << "  closure residual p1: " << gLastECalFixedEffectsSlope
            << " +/- " << gLastECalFixedEffectsSlopeError << " ns/ns\n";
  const double closureTolerance =
      std::max(0.01, 3.0*gLastECalFixedEffectsSlopeError);
  if (!std::isfinite(gLastECalFixedEffectsSlope) ||
      !std::isfinite(gLastECalFixedEffectsSlopeError) ||
      std::fabs(gLastECalFixedEffectsSlope) > closureTolerance) {
    std::cerr << "[Run p1 calibration] ERROR: closure residual exceeds "
              << closureTolerance << " ns/ns.\n";
    return;
  }
  std::cout << "  closure gate: PASS (|residual| <= "
            << closureTolerance << " ns/ns)\n";
  gLastCalibrationSequenceSucceeded = true;
}
