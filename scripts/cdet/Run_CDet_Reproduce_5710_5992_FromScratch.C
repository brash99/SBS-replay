#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

namespace {
bool ReadCDetIniValue(const char *filename, const std::string &section,
                      const std::string &key, double &value)
{
  std::ifstream input(filename);
  if (!input) return false;
  std::string line;
  bool inSection = false;
  while (std::getline(input, line)) {
    const std::string::size_type first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos || line[first] == '#') continue;
    const std::string trimmed = line.substr(first);
    if (trimmed.front() == '[') {
      inSection = trimmed == "[" + section + "]";
      continue;
    }
    if (!inSection) continue;
    std::istringstream parser(trimmed);
    std::string candidate;
    if (parser >> candidate >> value && candidate == key) return true;
  }
  return false;
}

bool RequireCDetValue(const char *filename, const char *section,
                      const char *key, double expected,
                      double tolerance = 0.0000005)
{
  double observed = std::numeric_limits<double>::quiet_NaN();
  if (!ReadCDetIniValue(filename, section, key, observed)) {
    std::cerr << "[5710/5992 reproduction] ERROR: missing [" << section
              << "] " << key << " in " << filename << ".\n";
    return false;
  }
  if (!std::isfinite(observed) || std::fabs(observed - expected) > tolerance) {
    std::cerr << std::setprecision(9)
              << "[5710/5992 reproduction] ERROR: " << filename << " ["
              << section << "] " << key << " = " << observed
              << "; expected " << expected << " +/- " << tolerance << ".\n";
    return false;
  }
  return true;
}

bool WriteRun5992MasterP1Seed(double p1)
{
  std::ofstream output("CDet_run5992.dat");
  if (!output) return false;
  output << "# Run-specific timing calibration for cross-target run 5992.\n\n"
         << "[ECalTiming]\n"
         << "p1 " << std::fixed << std::setprecision(6) << p1 << "\n";
  output.flush();
  return static_cast<bool>(output);
}
} // namespace

// Reproduce the commissioned cross-target calibration in one fresh ROOT
// process. Run from scripts/cdet with all three output files absent:
//
//   root -l
//   .L Run_CDet_Reproduce_5710_5992_FromScratch.C+
//   Run_CDet_Reproduce_5710_5992_FromScratch()
//
// The authoritative configurations and reviewed Run 5710 polygon artifact are
// inputs. CDet_calibration_dt.dat, CDet_run5710.dat, and CDet_run5992.dat are
// outputs. No run-specific p1 fit is performed for Run 5992.
void Run_CDet_Reproduce_5710_5992_FromScratch(
    Int_t nevents = std::numeric_limits<Int_t>::min())
{
  const char *outputs[] = {
      "CDet_calibration_dt.dat", "CDet_run5710.dat", "CDet_run5992.dat"};
  for (const char *output : outputs) {
    if (!gSystem->AccessPathName(output)) {
      std::cerr << "[5710/5992 reproduction] ERROR: clean start required; move "
                << output << " out of scripts/cdet and restart ROOT.\n";
      return;
    }
  }

  gROOT->ProcessLine(".L Run_CDet_Calibrate_Run5710_FromScratch.C+");
  gROOT->ProcessLine(TString::Format(
      "Run_CDet_Calibrate_Run5710_FromScratch(%d);", nevents));
  if (!static_cast<bool>(
          gROOT->ProcessLine("gLastCalibrationSequenceSucceeded"))) {
    std::cerr << "[5710/5992 reproduction] ERROR: Run 5710 failed.\n";
    return;
  }

  bool run5710Matches = true;
  run5710Matches &= RequireCDetValue(
      "CDet_calibration_dt.dat", "ECalTiming", "p0", 14.807420);
  run5710Matches &= RequireCDetValue(
      "CDet_calibration_dt.dat", "ECalTiming", "p1", 0.810203);
  run5710Matches &= RequireCDetValue(
      "CDet_calibration_dt.dat", "ECalTiming", "delta", 31.680784);
  run5710Matches &= RequireCDetValue(
      "CDet_calibration_dt.dat", "TimeWalk", "p1_L1", 13.987991);
  run5710Matches &= RequireCDetValue(
      "CDet_calibration_dt.dat", "TimeWalk", "p1_L2", 15.764852);
  run5710Matches &= RequireCDetValue(
      "CDet_run5710.dat", "GlobalTiming", "shift_ns", 1.195534);
  if (!run5710Matches) {
    std::cerr << "[5710/5992 reproduction] ERROR: Run 5710 output does not "
                 "match the commissioned calibration; Run 5992 was not run.\n";
    return;
  }

  // Store the generated master p1 explicitly in the Run 5992 file so the
  // artifact documents the fixed-slope decision and matches production.
  double masterP1 = 0.0;
  if (!ReadCDetIniValue(
          "CDet_calibration_dt.dat", "ECalTiming", "p1", masterP1) ||
      !WriteRun5992MasterP1Seed(masterP1)) {
    std::cerr << "[5710/5992 reproduction] ERROR: could not seed Run 5992 "
                 "with the generated Run 5710 p1.\n";
    return;
  }

  gROOT->ProcessLine(".L Run_CDet_Calibrate_RunTimingShift.C+");
  gROOT->ProcessLine(TString::Format(
      "Run_CDet_Calibrate_RunTimingShift(5992, %d, 30.0, 5.0, "
      "\"CDet_run5992_projection.conf\");", nevents));
  if (!static_cast<bool>(
          gROOT->ProcessLine("gLastCalibrationSequenceSucceeded"))) {
    std::cerr << "[5710/5992 reproduction] ERROR: Run 5992 shift fit failed.\n";
    return;
  }

  bool run5992Matches = true;
  run5992Matches &= RequireCDetValue(
      "CDet_run5992.dat", "ECalTiming", "p1", 0.810203);
  run5992Matches &= RequireCDetValue(
      "CDet_run5992.dat", "GlobalTiming", "shift_ns", 5.171528);
  if (!run5992Matches) {
    std::cerr << "[5710/5992 reproduction] ERROR: Run 5992 output does not "
                 "match the commissioned calibration.\n";
    return;
  }

  std::cout << "\n[5710/5992 reproduction] COMPLETE AND EXACT\n"
            << "  CDet_calibration_dt.dat: commissioned Run 5710 master\n"
            << "  CDet_run5710.dat: shift_ns = 1.195534 ns\n"
            << "  CDet_run5992.dat: p1 = 0.810203, shift_ns = 5.171528 ns\n"
            << "  Run 5992 p1 was inherited from Run 5710, not refitted.\n";
}
