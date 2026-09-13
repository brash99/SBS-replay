#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLine.h>
#include <TString.h>
#include <TStyle.h>
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
std::string TrimECalSelectionLine(const std::string &line)
{
  const std::string whitespace = " \t\r\n";
  const size_t first = line.find_first_not_of(whitespace);
  if (first == std::string::npos) return "";
  const size_t last = line.find_last_not_of(whitespace);
  return line.substr(first, last - first + 1);
}

bool WriteApprovedECalWindow(const TString &filename, Int_t,
                             double timeMin, double timeMax)
{
  std::vector<std::string> lines;
  std::ifstream input(filename.Data());
  if (!input) return false;
  std::string line;
  while (std::getline(input, line)) lines.push_back(line);

  std::ostringstream minimum, maximum;
  minimum << "analysis.ecal_time_min: " << std::fixed
          << std::setprecision(3) << timeMin;
  maximum << "analysis.ecal_time_max: " << std::fixed
          << std::setprecision(3) << timeMax;
  bool minimumWritten = false, maximumWritten = false;
  std::vector<std::string> updated;
  for (const std::string &line : lines) {
    const std::string trimmed = TrimECalSelectionLine(line);
    if (trimmed.rfind("analysis.ecal_time_min:", 0) == 0) {
      if (!minimumWritten) updated.push_back(minimum.str());
      minimumWritten = true;
    } else if (trimmed.rfind("analysis.ecal_time_max:", 0) == 0) {
      if (!maximumWritten) updated.push_back(maximum.str());
      maximumWritten = true;
    } else {
      updated.push_back(line);
    }
  }
  if (!minimumWritten || !maximumWritten) return false;

  const TString temporary = filename + ".tmp";
  std::ofstream output(temporary.Data());
  if (!output) return false;
  for (const std::string &line : updated) output << line << '\n';
  output.close();
  if (!output || std::rename(temporary.Data(), filename.Data()) != 0) {
    std::remove(temporary.Data());
    return false;
  }
  return true;
}

std::pair<double, double> ProposeWindow(TH1D &input)
{
  TH1D smoothed(input);
  smoothed.SetDirectory(nullptr);
  smoothed.Smooth(5);
  int peakBin = smoothed.GetMaximumBin();
  const int sideDistance = std::max(10, static_cast<int>(25.0/input.GetBinWidth(1)));
  const int leftReference = std::max(1, peakBin - sideDistance);
  const int rightReference = std::min(smoothed.GetNbinsX(), peakBin + sideDistance);
  const double leftBackground = smoothed.GetBinContent(leftReference);
  const double rightBackground = smoothed.GetBinContent(rightReference);
  const double peak = smoothed.GetBinContent(peakBin);
  auto background = [&](int bin) {
    const double fraction = static_cast<double>(bin - leftReference)/
        std::max(1, rightReference - leftReference);
    return leftBackground + fraction*(rightBackground - leftBackground);
  };
  auto threshold = [&](int bin) {
    const double localBackground = background(bin);
    return localBackground + 0.02*(peak - background(peakBin));
  };
  int lowBin = peakBin, highBin = peakBin;
  while (lowBin > leftReference && smoothed.GetBinContent(lowBin) > threshold(lowBin)) --lowBin;
  while (highBin < rightReference && smoothed.GetBinContent(highBin) > threshold(highBin)) ++highBin;
  const double low = std::floor(2.0*input.GetXaxis()->GetBinLowEdge(lowBin))/2.0;
  const double high = std::ceil(2.0*input.GetXaxis()->GetBinUpEdge(highBin))/2.0;
  return {low, high};
}
}

// First call with approve=false. Inspect the PNG/PDF, choose the bounds, then
// call again with approve=true and the human-approved min/max. A proposal is
// never silently promoted into CDet_run<run>_projection.conf.
void Run_CDet_Select_ECalTimingWindow(Int_t runNumber,
                                      bool approve = false,
                                      double approvedMin = 0.0,
                                      double approvedMax = 0.0)
{
  if (runNumber <= 0) {
    std::cerr << "[ECal window] ERROR: runNumber must be positive.\n";
    return;
  }
  if (approve) {
    if (!std::isfinite(approvedMin) || !std::isfinite(approvedMax) ||
        approvedMin >= approvedMax) {
      std::cerr << "[ECal window] ERROR: approval requires valid explicit bounds.\n";
      return;
    }
    const TString configFile = Form("CDet_run%d_projection.conf", runNumber);
    if (!WriteApprovedECalWindow(configFile, runNumber, approvedMin, approvedMax)) {
      std::cerr << "[ECal window] ERROR: could not update existing "
                << configFile << ".\n";
      return;
    }
    std::cout << "[ECal window] Stored HUMAN-APPROVED window " << approvedMin
              << " < time < " << approvedMax << " ns in " << configFile << ".\n";
    return;
  }

  const char *outDirectory = gSystem->Getenv("OUT_DIR");
  if (!outDirectory || !*outDirectory) {
    std::cerr << "[ECal window] ERROR: OUT_DIR is not defined. Source the replay setup first.\n";
    return;
  }
  TChain chain("T");
  const TString pattern = Form("%s/cdet_%d_stream*_seg*_*.root", outDirectory, runNumber);
  const int files = chain.Add(pattern);
  if (files <= 0 || chain.GetEntries() <= 0) {
    std::cerr << "[ECal window] ERROR: no input trees match " << pattern << ".\n";
    return;
  }

  TH1D spectrum(Form("hECalATime%d", runNumber),
      Form("Run %d: earm.ecal.a_time;earm.ecal.a_time (ns);ECal block hits / 0.5 ns", runNumber),
      800, -200.0, 200.0);
  chain.Draw(Form("earm.ecal.a_time>>hECalATime%d", runNumber), "", "goff");
  const auto proposal = ProposeWindow(spectrum);

  TCanvas canvas(Form("cECalATime%d", runNumber), "ECal timing-window review", 1500, 700);
  canvas.Divide(2, 1);
  for (int pad = 1; pad <= 2; ++pad) {
    canvas.cd(pad);
    TH1D *copy = static_cast<TH1D *>(spectrum.Clone(Form("hECalATime%dPad%d", runNumber, pad)));
    copy->GetXaxis()->SetRangeUser(-120.0, 160.0);
    if (pad == 2) {
      gPad->SetLogy();
      copy->SetMinimum(0.7);
    }
    gPad->SetGridx();
    copy->Draw("HIST");
    gPad->Update();
    const double lineMinimum = pad == 2 ? 0.7 : 0.0;
    const double lineMaximum = 1.05*copy->GetMaximum();
    TLine *low = new TLine(proposal.first, lineMinimum, proposal.first, lineMaximum);
    TLine *high = new TLine(proposal.second, lineMinimum, proposal.second, lineMaximum);
    low->SetLineColor(kRed + 1); high->SetLineColor(kRed + 1);
    low->SetLineWidth(3); high->SetLineWidth(3);
    low->Draw(); high->Draw();
  }
  const TString base = Form("CDet_run%d_ECalTimingWindow_PROPOSAL", runNumber);
  canvas.SaveAs(base + ".png");
  canvas.SaveAs(base + ".pdf");
  TFile output(base + ".root", "RECREATE");
  spectrum.Write();
  output.Close();

  std::cout << "\n[ECal window] AUTOMATIC PROPOSAL ONLY: " << proposal.first
            << " < earm.ecal.a_time < " << proposal.second << " ns\n"
            << "  Inspect " << base << ".png before approval.\n"
            << "  To record explicit human approval:\n"
            << "  Run_CDet_Select_ECalTimingWindow(" << runNumber
            << ", true, APPROVED_MIN, APPROVED_MAX)\n";
}
