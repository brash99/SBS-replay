#include <TDirectory.h>
#include <TH1.h>
#include <TList.h>
#include <TObject.h>
#include <TSystem.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

namespace {
struct CDetRunTimingMean {
  int run;
  int calibrationStage;
  double meanLe;
  double meanECalAdcTime;
  std::vector<double> groupMeanLe;
  std::size_t leHits;
  std::size_t ecalEvents;
};

std::unordered_set<TObject*> SnapshotDirectoryObjects() {
  std::unordered_set<TObject*> objects;
  if (!gDirectory || !gDirectory->GetList()) return objects;

  TIter next(gDirectory->GetList());
  while (TObject *object = next()) objects.insert(object);
  return objects;
}

void DeleteNewRunHistograms(const std::unordered_set<TObject*>& objectsBefore) {
  if (!gDirectory || !gDirectory->GetList()) return;

  std::vector<TH1*> histograms;
  TIter next(gDirectory->GetList());
  while (TObject *object = next()) {
    if (objectsBefore.count(object) == 0) {
      if (TH1 *histogram = dynamic_cast<TH1*>(object)) histograms.push_back(histogram);
    }
  }
  for (TH1 *histogram : histograms) delete histogram;
}

bool ReadCDetMeanRunList(const char *runList, std::vector<int>& runs) {
  std::ifstream input(runList);
  if (!input) {
    std::cerr << "[Run means] ERROR: could not open " << runList << ".\n";
    return false;
  }

  std::string line;
  while (std::getline(input, line)) {
    const std::size_t first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos || line[first] == '#' || line[first] == '[')
      continue;

    std::istringstream parser(line.substr(first));
    int run = -1;
    if (parser >> run && run > 0) runs.push_back(run);
  }

  if (runs.empty()) {
    std::cerr << "[Run means] ERROR: no run numbers found in " << runList << ".\n";
    return false;
  }
  return true;
}

bool WriteCDetRunMeanCsvFiles(const std::vector<CDetRunTimingMean>& results,
                              const char *leFile,
                              const char *ecalFile) {
  const std::string leTemporary = std::string(leFile) + ".tmp";
  const std::string ecalTemporary = std::string(ecalFile) + ".tmp";

  std::ofstream leOutput(leTemporary);
  std::ofstream ecalOutput(ecalTemporary);
  if (!leOutput || !ecalOutput) {
    std::cerr << "[Run means] ERROR: could not create temporary CSV files.\n";
    std::remove(leTemporary.c_str());
    std::remove(ecalTemporary.c_str());
    return false;
  }

  leOutput << "run_number,calibration_stage,mean_le_ns,mean_le_pixels_1472_1487_ns,"
              "mean_le_pixels_1200_1215_ns,mean_le_pixels_480_495_ns,n_le_hits\n";
  ecalOutput << "run_number,calibration_stage,mean_ecal_adctime_ns,n_ecal_events\n";
  for (const CDetRunTimingMean& result : results) {
    leOutput << result.run << "," << result.calibrationStage
             << "," << result.meanLe;
    for (int group = 0; group < 3; ++group)
      leOutput << "," << result.groupMeanLe[group];
    leOutput << "," << result.leHits << "\n";

    ecalOutput << result.run << "," << result.calibrationStage
               << "," << result.meanECalAdcTime
               << "," << result.ecalEvents << "\n";
  }
  leOutput.close();
  ecalOutput.close();

  if (!leOutput || !ecalOutput) {
    std::cerr << "[Run means] ERROR: failed while writing temporary CSV files.\n";
    std::remove(leTemporary.c_str());
    std::remove(ecalTemporary.c_str());
    return false;
  }

  if (std::rename(leTemporary.c_str(), leFile) != 0 ||
      std::rename(ecalTemporary.c_str(), ecalFile) != 0) {
    std::cerr << "[Run means] ERROR: could not replace the final CSV files.\n";
    std::remove(leTemporary.c_str());
    std::remove(ecalTemporary.c_str());
    return false;
  }
  return true;
}
} // namespace

void Run_CDet_Extract_RunTimingMeans_FromList(
    const char *runList = "runs.txt", Int_t events = -1,
    Int_t calibrationStage = 2,
    Int_t minSegment = -1, Int_t maxSegment = -1,
    Double_t leMin = 0.02, Double_t leMax = 100.0,
    Double_t totMin = 0.02, Double_t totMax = 150.0,
    Double_t ecalTimeMin = -10.0, Double_t ecalTimeMax = 50.0,
    const char *leCsv = "cdet_le_means.csv",
    const char *ecalCsv = "ecal_adctime_means.csv") {
  std::vector<int> runs;
  if (!ReadCDetMeanRunList(runList, runs)) return;
  if (calibrationStage < 0 || calibrationStage > 8) {
    std::cerr << "[Run means] ERROR: calibration stage must be between 0 and 8.\n";
    return;
  }

  std::vector<CDetRunTimingMean> results;
  results.reserve(runs.size());
  const bool previousDisableRunTimingConstants = gDisableRunTimingConstants;
  const bool previousAddDirectory = TH1::AddDirectoryStatus();
  gDisableRunTimingConstants = true;
  TH1::AddDirectory(kTRUE);

  for (std::size_t index = 0; index < runs.size(); ++index) {
    const int run = runs[index];
    std::cout << "\n[Run means] Processing run " << run << " ("
              << index + 1 << "/" << runs.size() << ")\n";
    const std::unordered_set<TObject*> objectsBefore = SnapshotDirectoryObjects();

    PlotElastic_Calibration_Master_stageflag_singlefile_crosstarget(
        run, events, calibrationStage, 0, minSegment, maxSegment,
        leMin, leMax, totMin, totMax, ecalTimeMin, ecalTimeMax,
        1, 100, 0, 100, 0.05, 0.0, 0.1, 3,
        false, 30, 2, 1, false);

    if (!gLastCalibrationStageSucceeded ||
        !std::isfinite(gLastRunMeanGoodLe) ||
        !std::isfinite(gLastRunMeanECalAdcTime) ||
        gLastRunGoodLeCount == 0 || gLastRunGoodECalEventCount == 0) {
      DeleteNewRunHistograms(objectsBefore);
      std::cerr << "[Run means] ERROR: run " << run
                << " did not produce valid timing means. Existing CSV files "
                   "were left unchanged.\n"
                << "  calibration stage: " << calibrationStage << "\n"
                << "  events processed: " << gNumEventsInRun << "\n"
                << "  good LE hits: " << gLastRunGoodLeCount << "\n"
                << "  good ECal events: " << gLastRunGoodECalEventCount << "\n"
                << "  mean good LE: " << gLastRunMeanGoodLe << " ns\n"
                << "  mean ECal ADC time: " << gLastRunMeanECalAdcTime << " ns\n";
      gDisableRunTimingConstants = previousDisableRunTimingConstants;
      TH1::AddDirectory(previousAddDirectory);
      return;
    }

    results.push_back({run, calibrationStage,
                       gLastRunMeanGoodLe, gLastRunMeanECalAdcTime,
                       gLastRunGroupMeanGoodLe, gLastRunGoodLeCount,
                       gLastRunGoodECalEventCount});
    DeleteNewRunHistograms(objectsBefore);
  }

  gDisableRunTimingConstants = previousDisableRunTimingConstants;
  TH1::AddDirectory(previousAddDirectory);
  if (!WriteCDetRunMeanCsvFiles(results, leCsv, ecalCsv)) return;
  std::cout << "\n[Run means] Completed all " << results.size()
            << " runs. Replaced " << leCsv << " and " << ecalCsv << ".\n";
}
