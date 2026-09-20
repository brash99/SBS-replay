#include <TString.h>
#include <TSystem.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {
std::string TrimInstallLine(const std::string &line)
{
  const std::string whitespace = " \t\r\n";
  const size_t first = line.find_first_not_of(whitespace);
  if (first == std::string::npos) return "";
  return line.substr(first, line.find_last_not_of(whitespace) - first + 1);
}

bool CopyInstallFile(const TString &source, const TString &destination)
{
  std::ifstream input(source.Data(), std::ios::binary);
  std::ofstream output(destination.Data(), std::ios::binary);
  if (!input || !output) return false;
  output << input.rdbuf();
  return input.good() || input.eof() ? output.good() : false;
}

bool WriteQualifiedP1(const TString &filename, double p1)
{
  std::ifstream input(filename.Data());
  if (!input) return false;
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(input, line)) lines.push_back(line);

  std::ostringstream value;
  value << "p1 " << std::fixed << std::setprecision(6) << p1;
  bool inSection = false, sectionSeen = false, written = false;
  std::vector<std::string> updated;
  for (const std::string &original : lines) {
    const std::string trimmed = TrimInstallLine(original);
    const bool section = trimmed.size() >= 2 && trimmed.front() == '[' &&
                         trimmed.back() == ']';
    if (section) {
      if (inSection && !written) {
        updated.push_back(value.str());
        written = true;
      }
      inSection = trimmed == "[ECalTiming]";
      sectionSeen = sectionSeen || inSection;
      updated.push_back(original);
      continue;
    }
    std::istringstream parser(trimmed);
    std::string key;
    parser >> key;
    if (inSection && key == "p1") {
      if (!written) updated.push_back(value.str());
      written = true;
    } else {
      updated.push_back(original);
    }
  }
  if (inSection && !written) updated.push_back(value.str());
  if (!sectionSeen) {
    if (!updated.empty() && !updated.back().empty()) updated.push_back("");
    updated.push_back("[ECalTiming]");
    updated.push_back(value.str());
  }

  const TString temporary = filename + ".tmp";
  std::ofstream output(temporary.Data());
  if (!output) return false;
  for (const std::string &entry : updated) output << entry << "\n";
  output.close();
  if (!output) {
    std::remove(temporary.Data());
    return false;
  }
  if (std::rename(temporary.Data(), filename.Data()) != 0) {
    std::remove(temporary.Data());
    return false;
  }
  return true;
}
} // namespace

// Install only a candidate from a successful frozen-sample qualification.
// The current run file is preserved once before the atomic update.
void Run_CDet_Install_QualifiedP1(Int_t runNumber,
    TString qualificationFile = "")
{
  if (runNumber <= 0) {
    std::cerr << "[Install qualified p1] ERROR: invalid run number.\n";
    return;
  }
  if (qualificationFile.IsNull())
    qualificationFile = TString::Format(
        "CDet_run%d_frozen_p1_qualification.txt", runNumber);
  std::ifstream input(qualificationFile.Data());
  if (!input) {
    std::cerr << "[Install qualified p1] ERROR: cannot read "
              << qualificationFile << ".\n";
    return;
  }
  int artifactRun = -1, qualified = 0;
  double candidate = std::numeric_limits<double>::quiet_NaN();
  std::string key;
  while (input >> key) {
    std::string rest;
    if (key == "run") input >> artifactRun;
    else if (key == "candidate_p1") input >> candidate;
    else if (key == "qualified") input >> qualified;
    else std::getline(input, rest);
  }
  if (artifactRun != runNumber || qualified != 1 || !std::isfinite(candidate)) {
    std::cerr << "[Install qualified p1] ERROR: artifact is not a qualified "
                 "result for run " << runNumber << ".\n";
    return;
  }

  const TString runFile = TString::Format("CDet_run%d.dat", runNumber);
  const TString backup = TString::Format(
      "CDet_run%d_before_qualified_p1.dat", runNumber);
  if (gSystem->AccessPathName(runFile)) {
    std::cerr << "[Install qualified p1] ERROR: missing " << runFile << ".\n";
    return;
  }
  if (!gSystem->AccessPathName(backup)) {
    std::cerr << "[Install qualified p1] ERROR: backup already exists: "
              << backup << ". Preserve or rename it before retrying.\n";
    return;
  }
  if (!CopyInstallFile(runFile, backup) ||
      !WriteQualifiedP1(runFile, candidate)) {
    std::cerr << "[Install qualified p1] ERROR: installation failed. The "
                 "preserved input is " << backup << ".\n";
    return;
  }
  std::cout << "[Install qualified p1] Installed p1=" << candidate
            << " in " << runFile << ".\n"
            << "  preserved prior file: " << backup << "\n"
            << "  next required operation: refit shift_ns once.\n";
}
