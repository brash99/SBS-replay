#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

/* Wrapper script to call PlotElastic routine many times for different run numbers in .txt file given in constructor */

void Run_CDet_Calibration_FromList(const char* runlist = "crossRuns.txt", bool removeCalib = true){
    gROOT->SetBatch(kTRUE);
    TH1::AddDirectory(kFALSE);

    const TString calibFile = "CDet_calibration.dat";
    const TString outDir    = "calibrationFiles";

    // Make output directory if needed
    if (gSystem->AccessPathName(outDir)) {
        std::cout << "[Wrapper] Creating directory: " << outDir << "\n";
        gSystem->mkdir(outDir, true);
    }

    std::ifstream infile(runlist);
    if (!infile) {
        std::cerr << "[Wrapper] ERROR: Could not open " << runlist << "\n";
        return;
    }

    std::vector<int> runs;
    std::string line;

    while (std::getline(infile, line)) {
        // Trim leading whitespace
        size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;

        line = line.substr(first);

        // Skip comments
        if (line[0] == '#') continue;

        // Skip group headers like [Group 1], [Group 2], etc.
        if (line[0] == '[') continue;

        // Read run number from line
        std::stringstream ss(line);
        int run = -1;

        if (ss >> run) {
            runs.push_back(run);
        }
    }

    if (runs.empty()) {
        std::cerr << "[Wrapper] ERROR: No runs found in " << runlist << "\n";
        return;
    }

    std::cout << "[Wrapper] Found " << runs.size() << " runs\n";

    for (size_t i = 0; i < runs.size(); i++) {
        int thisRun = runs[i];

        std::cout << "\n========================================\n";
        std::cout << "[Wrapper] Processing run " << thisRun
                  << " (" << i+1 << "/" << runs.size() << ")\n";
        std::cout << "========================================\n";

        Run_CDet_Calibration_TwoPass_InSession_AllCross_Individually(
            thisRun,     // RunNumber1
            -1,          // nevents
            0,           // elastic
            0, 5,        // segments
            10.0, 50.0,  // LE min, max
            4.0, 50.0,   // TOT min, max
            1, 100,
            1, 100,
            0.10,
            0.02,
            0.1,
            3,
            false,
            30,
            2,
            1,
            removeCalib
        );

        if (!gLastCalibrationSequenceSucceeded) {
            std::cerr << "[Wrapper] ERROR: calibration failed for run "
                      << thisRun << "; no output will be archived.\n";
            continue;
        }

        // After the run finishes, save a copy of the calibration file
        if (gSystem->AccessPathName(calibFile)) {
            std::cerr << "[Wrapper] WARNING: " << calibFile
                      << " was not found after run " << thisRun << "\n";
            continue;
        }

        TString outFile = TString::Format("%s/%s.%d",
                                          outDir.Data(),
                                          calibFile.Data(),
                                          thisRun);

        int copyStatus = gSystem->CopyFile(calibFile, outFile, true);
        if (copyStatus != 0) {
            std::cerr << "[Wrapper] WARNING: Failed to copy "
                      << calibFile << " to " << outFile << "\n";
        } else {
            std::cout << "[Wrapper] Saved calibration file to "
                      << outFile << "\n";
        }
    }

    std::cout << "\n[Wrapper] All runs completed.\n";
}
