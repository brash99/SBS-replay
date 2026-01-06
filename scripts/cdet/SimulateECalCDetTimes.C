#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "TH2D.h"
#include "TCanvas.h"
#include "TRandom3.h"

// Helper to load a 1D distribution from CSV: columns x_center,content
bool LoadTimeDistribution(const char* filename,
                          std::vector<double>& centers,
                          std::vector<double>& weights,
                          double& totalContent)
{
    centers.clear();
    weights.clear();
    totalContent = 0.0;

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return false;
    }

    std::string line;
    bool firstLine = true;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        std::stringstream ss(line);
        std::string col1, col2;

        if (!std::getline(ss, col1, ',')) continue;
        if (!std::getline(ss, col2, ',')) continue;

        // Skip header if present
        if (firstLine) {
            firstLine = false;
            if (!std::isdigit(col1[0]) && col1.find("x_center") != std::string::npos)
                continue;
        }

        double x = std::atof(col1.c_str());
        double w = std::atof(col2.c_str());

        centers.push_back(x);
        weights.push_back(w);
        totalContent += w;
    }

    if (centers.empty()) {
        std::cerr << "Error: no data read from " << filename << std::endl;
        return false;
    }

    return true;
}

// Build CDF from weights
void BuildCDF(const std::vector<double>& weights,
              std::vector<double>& cdf,
              double totalContent)
{
    cdf.resize(weights.size());
    double running = 0.0;
    for (std::size_t i = 0; i < weights.size(); ++i) {
        running += weights[i] / totalContent;
        cdf[i] = running;
    }
    // Ensure last element is exactly 1.0 (avoid rounding issues)
    if (!cdf.empty()) cdf.back() = 1.0;
}

// Sample from discrete distribution defined by centers + CDF
double SampleFromDistribution(TRandom3& rng,
                              const std::vector<double>& centers,
                              const std::vector<double>& cdf)
{
    double u = rng.Rndm(); // uniform in [0,1)
    auto it = std::upper_bound(cdf.begin(), cdf.end(), u);
    std::size_t idx = std::distance(cdf.begin(), it);
    if (idx >= centers.size()) idx = centers.size() - 1; // safety
    return centers[idx];
}

/**
 * Simulate ECal vs CDet timing correlations.
 *
 * - Reads ECalTime.csv and CDetTime.csv (columns: x_center,content)
 * - Number of events = sum(content) from the files
 * - Fills a TH2D with x = cdet_time, y = ecal_time
 * - Draws a 900x700 canvas with COLZ
 *
 * You can call this from any ROOT macro:
 *    root [0] .L YourMacro.C+
 *    root [1] SimulateECalCDetTimes();
 */
void SimulateECalCDetTimes(const char* ecalFile = "ECalTime.csv",
                           const char* cdetFile = "CDetTime.csv")
{
    // --- Load distributions ---
    std::vector<double> ecalCenters, ecalWeights;
    std::vector<double> cdetCenters, cdetWeights;
    double ecalTotal = 0.0, cdetTotal = 0.0;

    if (!LoadTimeDistribution(ecalFile, ecalCenters, ecalWeights, ecalTotal)) return;
    if (!LoadTimeDistribution(cdetFile, cdetCenters, cdetWeights, cdetTotal)) return;

    if (std::fabs(ecalTotal - cdetTotal) > 1e-6) {
        std::cerr << "Warning: total contents differ between files: "
                  << "ECal=" << ecalTotal << ", CDet=" << cdetTotal << std::endl;
    }

    Long64_t nEvents = static_cast<Long64_t>(std::llround(ecalTotal));
    std::cout << "Simulating " << nEvents << " events." << std::endl;

    // --- Build CDFs for sampling ---
    std::vector<double> ecalCDF, cdetCDF;
    BuildCDF(ecalWeights, ecalCDF, ecalTotal);
    BuildCDF(cdetWeights, cdetCDF, cdetTotal);

    // --- Define 2D histogram binning from centers ---
    int nXBins = static_cast<int>(cdetCenters.size());
    int nYBins = static_cast<int>(ecalCenters.size());

    double dx = (nXBins > 1) ? (cdetCenters[1] - cdetCenters[0]) : 1.0;
    double dy = (nYBins > 1) ? (ecalCenters[1] - ecalCenters[0]) : 1.0;

    double xMin = cdetCenters.front() - 0.5 * dx;
    double xMax = cdetCenters.back()  + 0.5 * dx;
    double yMin = ecalCenters.front() - 0.5 * dy;
    double yMax = ecalCenters.back()  + 0.5 * dy;

    TH2D* hECalVsCDetSim = new TH2D("hECalVsCDetSim",
                                 "ECal time vs CDet time (Sim);CDet time [ns];ECal time [ns]",
                                 nXBins, xMin, xMax,
                                 nYBins, yMin, yMax);

    // --- Random generator ---
    TRandom3 rng(0); // seed with 0 => random seed

    // --- Event loop ---
    for (Long64_t i = 0; i < nEvents; ++i) {
        double ecal_time = SampleFromDistribution(rng, ecalCenters, ecalCDF);
        double cdet_time = SampleFromDistribution(rng, cdetCenters, cdetCDF);
	if (ecal_time-cdet_time>=70.0 && ecal_time-cdet_time <= 115.0) {
        	hECalVsCDetSim->Fill(cdet_time, ecal_time);
	}
    }

    // --- Draw the 2D histogram ---
    TCanvas* c = new TCanvas("cECalVsCDet", "ECal vs CDet time", 900, 700);
    c->cd();
    hECalVsCDetSim->Draw("COLZ");
    c->Update();
}
