#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>

#include "CDetGoodPulseConfig.h"

namespace {

TString CDetRootStringLiteral(const char *value) {
  TString escaped(value ? value : "");
  escaped.ReplaceAll("\\", "\\\\");
  escaped.ReplaceAll("\"", "\\\"");
  return TString::Format("\"%s\"", escaped.Data());
}

} // namespace

// Load both plotting macros and run them from one authoritative per-run
// configuration. Invoke this macro from scripts/cdet so the relative macro
// names resolve in both interactive and batch ROOT sessions.
void Run_CDet_GoodPulseDiagnostics(
    const char *configFile, const char *inputDirectory = nullptr,
    const char *goodPulseOutputDirectory = nullptr,
    const char *pairScanOutputDirectory = nullptr) {
  CDetGoodPulseConfig::Values config;
  if (!CDetGoodPulseConfig::Load(
          configFile, config, "CDet good-pulse diagnostic suite"))
    return;

  TString input(inputDirectory ? inputDirectory : "");
  if (input.IsNull()) {
    const char *outDir = gSystem->Getenv("OUT_DIR");
    if (outDir)
      input = outDir;
  }
  if (input.IsNull()) {
    std::cerr << "[CDet good-pulse diagnostic suite] An input directory or "
                 "OUT_DIR is required.\n";
    return;
  }

  const TString goodOutput =
      goodPulseOutputDirectory && goodPulseOutputDirectory[0]
          ? TString(goodPulseOutputDirectory)
          : TString::Format("CDet_run%d_good_pulse_tdc", config.runNumber);
  const TString scanOutput =
      pairScanOutputDirectory && pairScanOutputDirectory[0]
          ? TString(pairScanOutputDirectory)
          : TString::Format("CDet_run%d_pair_timing_scan", config.runNumber);

  gROOT->ProcessLine(".L Plot_CDet_GoodPulseCandidates_AllTDC.C+");
  gROOT->ProcessLine(TString::Format(
      "Plot_CDet_GoodPulseCandidates_AllTDC(%s,%s,%s);",
      CDetRootStringLiteral(configFile).Data(),
      CDetRootStringLiteral(input.Data()).Data(),
      CDetRootStringLiteral(goodOutput.Data()).Data()));

  gROOT->ProcessLine(".L Plot_CDet_PairYieldVsECalTimingWindow.C+");
  gROOT->ProcessLine(TString::Format(
      "Plot_CDet_PairYieldVsECalTimingWindow(%s,%s,%s);",
      CDetRootStringLiteral(configFile).Data(),
      CDetRootStringLiteral(input.Data()).Data(),
      CDetRootStringLiteral(scanOutput.Data()).Data()));
}
