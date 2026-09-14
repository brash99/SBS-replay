#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TPad.h>
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
// pass explicit min/max with approve=false to preview that window and zoom the
// plots to those bounds. Call with approve=true and the human-approved min/max.
// A proposal is
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
  const bool explicitWindow = approvedMin != 0.0 || approvedMax != 0.0;
  if ((approve || explicitWindow) &&
      (!std::isfinite(approvedMin) || !std::isfinite(approvedMax) ||
       approvedMin >= approvedMax)) {
    std::cerr << "[ECal window] ERROR: explicit bounds must be finite with min < max.\n";
    return;
  }
  if (approve) {
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

  TH1D spectrum(Form("hECalAdcTime%d", runNumber),
      Form("Run %d: earm.ecal.adctime;earm.ecal.adctime (ns);Events / 0.5 ns", runNumber),
      800, -200.0, 200.0);
  if (chain.Draw(Form("earm.ecal.adctime>>hECalAdcTime%d", runNumber), "", "goff") < 0) {
    std::cerr << "[ECal window] ERROR: could not read earm.ecal.adctime for ECal timing.\n";
    return;
  }
  const auto proposal = explicitWindow
      ? std::make_pair(approvedMin, approvedMax) : ProposeWindow(spectrum);

  TH1D hcalSpectrum(Form("hHCalAdcTime%d", runNumber),
      Form("Run %d: sbs.hcal.adctime (ECal bounds for comparison);sbs.hcal.adctime (ns);Events / 0.5 ns", runNumber),
      800, -200.0, 200.0);
  if (chain.Draw(Form("sbs.hcal.adctime>>hHCalAdcTime%d", runNumber), "", "goff") < 0) {
    std::cerr << "[ECal window] ERROR: could not read sbs.hcal.adctime for HCal comparison.\n";
    return;
  }

  TH2D correlation(Form("hHCalVsECalAdcTime%d", runNumber),
      Form("Run %d: HCal vs. ECal ADC time;earm.ecal.adctime (ns);sbs.hcal.adctime (ns);Events / (0.5 ns #times 0.5 ns)", runNumber),
      800, -200.0, 200.0, 800, -200.0, 200.0);
  if (chain.Draw(Form("sbs.hcal.adctime:earm.ecal.adctime>>hHCalVsECalAdcTime%d", runNumber), "", "goff") < 0) {
    std::cerr << "[ECal window] ERROR: could not read paired ECal/HCal ADC times.\n";
    return;
  }

  const double plotMin = explicitWindow ? approvedMin : -120.0;
  const double plotMax = explicitWindow ? approvedMax : 160.0;
  TCanvas canvas(Form("cECalATime%d", runNumber), "ECal and HCal timing-window review", 1500, 1200);
  canvas.Divide(1, 2);
  TPad *top = static_cast<TPad *>(canvas.cd(1));
  top->Divide(2, 1);
  for (int pad = 1; pad <= 2; ++pad) {
    top->cd(pad);
    TH1D &source = pad == 1 ? spectrum : hcalSpectrum;
    TH1D *copy = static_cast<TH1D *>(source.Clone(Form("%sPad%d", source.GetName(), pad)));
    copy->GetXaxis()->SetRangeUser(plotMin, plotMax);
    gPad->SetGridx();
    copy->Draw("HIST");
    gPad->Update();
    const double lineMinimum = 0.0;
    const double lineMaximum = std::max(1.0, 1.05*copy->GetMaximum());
    TLine *low = new TLine(proposal.first, lineMinimum, proposal.first, lineMaximum);
    TLine *high = new TLine(proposal.second, lineMinimum, proposal.second, lineMaximum);
    low->SetLineColor(kRed + 1); high->SetLineColor(kRed + 1);
    low->SetLineWidth(3); high->SetLineWidth(3);
    low->Draw(); high->Draw();
  }
  canvas.cd(2);
  gPad->SetRightMargin(0.16);
  correlation.SetStats(false);
  correlation.GetXaxis()->SetRangeUser(plotMin, plotMax);
  correlation.GetYaxis()->SetRangeUser(plotMin, plotMax);
  correlation.Draw("COLZ");
  for (double bound : {proposal.first, proposal.second}) {
    TLine *ecalBound = new TLine(bound, plotMin, bound, plotMax);
    TLine *hcalReference = new TLine(plotMin, bound, plotMax, bound);
    ecalBound->SetLineColor(kRed + 1);
    hcalReference->SetLineColor(kRed + 1);
    ecalBound->SetLineWidth(2);
    hcalReference->SetLineWidth(2);
    hcalReference->SetLineStyle(2);
    ecalBound->Draw();
    hcalReference->Draw();
  }
  const TString base = Form("CDet_run%d_ECalTimingWindow_PROPOSAL", runNumber);
  canvas.SaveAs(base + ".png");
  canvas.SaveAs(base + ".pdf");
  TFile output(base + ".root", "RECREATE");
  spectrum.Write();
  hcalSpectrum.Write();
  correlation.Write();
  output.Close();

  std::cout << "\n[ECal window] "
            << (explicitWindow ? "EXPLICIT WINDOW PREVIEW ONLY: " : "AUTOMATIC PROPOSAL ONLY: ")
            << proposal.first
            << " < earm.ecal.adctime < " << proposal.second << " ns\n"
            << "  Inspect " << base << ".png before approval.\n"
            << "  HCal ADC-time panels show the same bounds for comparison, not an applied HCal cut.\n"
            << "  To record explicit human approval:\n"
            << "  Run_CDet_Select_ECalTimingWindow(" << runNumber
            << ", true, APPROVED_MIN, APPROVED_MAX)\n";
}
