#include "CDetPlotStyle.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TProfile.h>
#include <TString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace {

void DrawPairAndSingle(TH1D &pair, TH1D &single, const char *xTitle)
{
  if (pair.Integral() > 0.0)
    pair.Scale(1.0 / pair.Integral());
  if (single.Integral() > 0.0)
    single.Scale(1.0 / single.Integral());
  pair.SetLineColor(kBlue + 1);
  pair.SetLineWidth(2);
  single.SetLineColor(kOrange + 7);
  single.SetLineWidth(2);
  pair.GetXaxis()->SetTitle(xTitle);
  pair.GetYaxis()->SetTitle("Fraction of hypotheses");
  const double ymax = 1.15 * std::max(pair.GetMaximum(), single.GetMaximum());
  pair.SetMaximum(ymax > 0.0 ? ymax : 1.0);
  pair.Draw("hist");
  single.Draw("hist same");
  TLegend *legend = new TLegend(0.62, 0.72, 0.87, 0.87);
  legend->SetBorderSize(0);
  legend->AddEntry(&pair, "Two-layer pair", "l");
  legend->AddEntry(&single, "Single layer", "l");
  legend->Draw();
}

void SetProfileDisplayRange(TProfile &profile,
                            double reference =
                                std::numeric_limits<double>::quiet_NaN(),
                            double fractionalPadding = 0.15)
{
  double ymin = std::numeric_limits<double>::infinity();
  double ymax = -std::numeric_limits<double>::infinity();
  for (int bin = 1; bin <= profile.GetNbinsX(); ++bin) {
    if (profile.GetBinEntries(bin) <= 0.0)
      continue;
    ymin = std::min(ymin, profile.GetBinContent(bin));
    ymax = std::max(ymax, profile.GetBinContent(bin));
  }
  if (!std::isfinite(ymin) || !std::isfinite(ymax))
    return;
  if (std::isfinite(reference)) {
    ymin = std::min(ymin, reference);
    ymax = std::max(ymax, reference);
  }
  const double span = std::max(ymax - ymin, 0.02 * std::max(std::abs(ymin), 1.0));
  profile.SetMinimum(ymin - fractionalPadding * span);
  profile.SetMaximum(ymax + fractionalPadding * span);
}

} // namespace

// Plot the Step-5 diagnostic relationship among ECal, CDet, and the existing
// target-z scan. This macro reads only exported ROOT-tree quantities and does
// not rerun candidate selection or modify the GEM search.
void Plot_CDet_FTROI_TargetZDiagnostics(
    const char *inputFile,
    const char *outputStem = "CDet_FTROI_TargetZDiagnostics",
    double centralEArmAngleDegrees = 27.0)
{
  ApplyCDetPlotStyle();
  gStyle->SetOptStat(0);

  TFile input(inputFile, "READ");
  if (input.IsZombie()) {
    std::cerr << "[CDet FTROI target-z] Cannot open " << inputFile << std::endl;
    return;
  }
  TTree *tree = dynamic_cast<TTree *>(input.Get("T"));
  if (!tree) {
    std::cerr << "[CDet FTROI target-z] Tree T is missing from " << inputFile
              << std::endl;
    return;
  }

  const char *required[] = {
      "earm.ecal.x", "earm.ecal.y", "earm.cdet.pulse.x_corr",
      "FTROI.cdet.hyp.n", "FTROI.cdet.hyp.source_type",
      "FTROI.cdet.hyp.pulse_index_l1", "FTROI.cdet.hyp.pulse_index_l2",
      "FTROI.cdet.vertex.n", "FTROI.cdet.vertex.hyp_index",
      "FTROI.cdet.vertex.z", "FTROI.cdet.vertex.xchi2",
      "FTROI.cdet.vertex.xndf", "FTROI.cdet.vertex.xslope",
      "FTROI.cdet.vertex.yslope", "FTROI.cdet.vertex.yresidual_l1",
      "FTROI.cdet.vertex.yresidual_l2", "FTROI.cdet.vertex.ycompatible",
      "FTROI.cdet.vertex.theta_global"};
  for (const char *name : required) {
    if (!tree->GetBranch(name)) {
      std::cerr << "[CDet FTROI target-z] Required branch is missing: "
                << name << std::endl;
      return;
    }
  }

  // The exported target positions are bin centers. Reconstruct the exact
  // database scan edges so that each scan point occupies one display bin.
  // Adding visual margins by widening the histogram while retaining 24 bins
  // would merge selected adjacent scan points and create false hot spots.
  const int targetBins = 24;
  const double targetCenterMin = tree->GetMinimum("FTROI.cdet.vertex.z");
  const double targetCenterMax = tree->GetMaximum("FTROI.cdet.vertex.z");
  const double targetBinWidth =
      (targetCenterMax - targetCenterMin) / double(targetBins - 1);
  const double targetMin = targetCenterMin - 0.5 * targetBinWidth;
  const double targetMax = targetCenterMax + 0.5 * targetBinWidth;
  if (!(targetBinWidth > 0.0)) {
    std::cerr << "[CDet FTROI target-z] Cannot determine target-z scan binning"
              << std::endl;
    return;
  }

  TTreeReader reader(tree);
  TTreeReaderValue<Double_t> ecalX(reader, "earm.ecal.x");
  TTreeReaderValue<Double_t> ecalY(reader, "earm.ecal.y");
  TTreeReaderArray<Double_t> pulseX(reader, "earm.cdet.pulse.x_corr");
  TTreeReaderArray<Double_t> hypSourceType(reader, "FTROI.cdet.hyp.source_type");
  TTreeReaderArray<Double_t> hypPulseL1(reader, "FTROI.cdet.hyp.pulse_index_l1");
  TTreeReaderArray<Double_t> hypPulseL2(reader, "FTROI.cdet.hyp.pulse_index_l2");
  TTreeReaderArray<Double_t> vertexHyp(reader, "FTROI.cdet.vertex.hyp_index");
  TTreeReaderArray<Double_t> vertexZ(reader, "FTROI.cdet.vertex.z");
  TTreeReaderArray<Double_t> vertexXChi2(reader, "FTROI.cdet.vertex.xchi2");
  TTreeReaderArray<Double_t> vertexXNDF(reader, "FTROI.cdet.vertex.xndf");
  TTreeReaderArray<Double_t> vertexXSlope(reader, "FTROI.cdet.vertex.xslope");
  TTreeReaderArray<Double_t> vertexYSlope(reader, "FTROI.cdet.vertex.yslope");
  TTreeReaderArray<Double_t> vertexYResidualL1(reader, "FTROI.cdet.vertex.yresidual_l1");
  TTreeReaderArray<Double_t> vertexYResidualL2(reader, "FTROI.cdet.vertex.yresidual_l2");
  TTreeReaderArray<Double_t> vertexYCompatible(reader, "FTROI.cdet.vertex.ycompatible");
  TTreeReaderArray<Double_t> vertexTheta(reader, "FTROI.cdet.vertex.theta_global");

  TH2D hECalVsL1("hECalVsL1",
                 "Layer 1 CDet x versus ECal x;ECal x (m);CDet Layer 1 x (m)",
                 120, -1.2, 1.2, 120, -1.2, 1.2);
  TH2D hECalVsL2("hECalVsL2",
                 "Layer 2 CDet x versus ECal x;ECal x (m);CDet Layer 2 x (m)",
                 120, -1.2, 1.2, 120, -1.2, 1.2);
  TH2D hXChi2VsZ("hXChi2VsZ",
                 "Target-constrained x fit;Target z (m);#chi^{2}/NDF",
                 targetBins, targetMin, targetMax, 160, 0.0, 40.0);
  TProfile pXChi2VsZ("pXChi2VsZ",
                     "Mean target-constrained x fit quality;Target z (m);Mean #chi^{2}/NDF",
                     targetBins, targetMin, targetMax);
  TH1D hBestZPair("hBestZPair", "", targetBins, targetMin, targetMax);
  TH1D hBestZSingle("hBestZSingle", "", targetBins, targetMin, targetMax);
  TH1D hMinChi2Pair("hMinChi2Pair", "", 160, 0.0, 40.0);
  TH1D hMinChi2Single("hMinChi2Single", "", 160, 0.0, 40.0);
  TProfile pYCompatibleVsZ(
      "pYCompatibleVsZ",
      "CDet half-bar y compatibility;Target z (m);Compatible fraction",
      targetBins, targetMin, targetMax, 0.0, 1.0);
  TH2D hYResidualL1VsZ(
      "hYResidualL1VsZ",
      "Layer 1 y compatibility residual;Target z (m);y_{CDet}-y_{ECal ray} (m)",
      targetBins, targetMin, targetMax, 160, -0.8, 0.8);
  TH2D hYResidualL2VsZ(
      "hYResidualL2VsZ",
      "Layer 2 y compatibility residual;Target z (m);y_{CDet}-y_{ECal ray} (m)",
      targetBins, targetMin, targetMax, 160, -0.8, 0.8);
  TProfile pThetaVsZ(
      "pThetaVsZ",
      "Global electron scattering angle;Target z (m);Mean #theta_{lab} (deg)",
      targetBins, targetMin, targetMax);
  TProfile pOutOfPlaneVsZ(
      "pOutOfPlaneVsZ",
      "Transport out-of-plane angle;Target z (m);Mean tan^{-1}(x') (deg)",
      targetBins, targetMin, targetMax);
  TProfile pInPlaneVsZ(
      "pInPlaneVsZ",
      "Transport in-plane deviation;Target z (m);Mean tan^{-1}(y') (deg)",
      targetBins, targetMin, targetMax);

  const double radiansToDegrees = 180.0 / TMath::Pi();

  Long64_t eventsWithHypotheses = 0;
  Long64_t hypotheses = 0;
  while (reader.Next()) {
    const int nHyp = static_cast<int>(hypSourceType.GetSize());
    if (nHyp <= 0)
      continue;
    ++eventsWithHypotheses;
    hypotheses += nHyp;

    std::vector<double> bestChi2NDF(nHyp,
        std::numeric_limits<double>::infinity());
    std::vector<double> bestZ(nHyp,
        std::numeric_limits<double>::quiet_NaN());

    for (int ih = 0; ih < nHyp; ++ih) {
      const int iL1 = static_cast<int>(std::lround(hypPulseL1[ih]));
      const int iL2 = static_cast<int>(std::lround(hypPulseL2[ih]));
      if (iL1 >= 0 && iL1 < static_cast<int>(pulseX.GetSize()))
        hECalVsL1.Fill(*ecalX, pulseX[iL1]);
      if (iL2 >= 0 && iL2 < static_cast<int>(pulseX.GetSize()))
        hECalVsL2.Fill(*ecalX, pulseX[iL2]);
    }

    const int nVertex = static_cast<int>(vertexZ.GetSize());
    for (int iv = 0; iv < nVertex; ++iv) {
      const int ih = static_cast<int>(std::lround(vertexHyp[iv]));
      if (ih < 0 || ih >= nHyp)
        continue;
      const double z = vertexZ[iv];
      const double ndf = vertexXNDF[iv];
      const double chi2NDF = ndf > 0.0 ? vertexXChi2[iv] / ndf :
          std::numeric_limits<double>::quiet_NaN();
      if (std::isfinite(chi2NDF)) {
        hXChi2VsZ.Fill(z, chi2NDF);
        pXChi2VsZ.Fill(z, chi2NDF);
        if (chi2NDF < bestChi2NDF[ih]) {
          bestChi2NDF[ih] = chi2NDF;
          bestZ[ih] = z;
        }
      }
      if (std::isfinite(vertexYResidualL1[iv]))
        hYResidualL1VsZ.Fill(z, vertexYResidualL1[iv]);
      if (std::isfinite(vertexYResidualL2[iv]))
        hYResidualL2VsZ.Fill(z, vertexYResidualL2[iv]);
      if (vertexYCompatible[iv] == 0.0 || vertexYCompatible[iv] == 1.0)
        pYCompatibleVsZ.Fill(z, vertexYCompatible[iv]);
      if (std::isfinite(vertexTheta[iv]))
        pThetaVsZ.Fill(z, vertexTheta[iv] * radiansToDegrees);
      if (std::isfinite(vertexXSlope[iv]))
        pOutOfPlaneVsZ.Fill(z, std::atan(vertexXSlope[iv]) *
                                radiansToDegrees);
      if (std::isfinite(vertexYSlope[iv]))
        pInPlaneVsZ.Fill(z, std::atan(vertexYSlope[iv]) *
                             radiansToDegrees);
    }

    for (int ih = 0; ih < nHyp; ++ih) {
      if (!std::isfinite(bestZ[ih]) || !std::isfinite(bestChi2NDF[ih]))
        continue;
      if (static_cast<int>(std::lround(hypSourceType[ih])) == 1) {
        hBestZPair.Fill(bestZ[ih]);
        hMinChi2Pair.Fill(bestChi2NDF[ih]);
      } else {
        hBestZSingle.Fill(bestZ[ih]);
        hMinChi2Single.Fill(bestChi2NDF[ih]);
      }
    }
  }

  TCanvas cGeometry("cCDetFTROIGeometry", "CDet/ECal geometry", 1500, 1000);
  cGeometry.Divide(2, 2);
  cGeometry.cd(1); hECalVsL1.Draw("colz");
  cGeometry.cd(2); hECalVsL2.Draw("colz");
  cGeometry.cd(3); hXChi2VsZ.Draw("colz");
  cGeometry.cd(4); pXChi2VsZ.SetLineColor(kRed + 1);
  pXChi2VsZ.SetLineWidth(2); SetProfileDisplayRange(pXChi2VsZ);
  pXChi2VsZ.SetMinimum(0.0);
  pXChi2VsZ.SetMaximum(5.0);
  pXChi2VsZ.Draw();

  TCanvas cTarget("cCDetFTROITarget", "CDet target-z diagnostics", 1500, 1000);
  cTarget.Divide(2, 2);
  cTarget.cd(1); DrawPairAndSingle(hBestZPair, hBestZSingle, "Best target z (m)");
  cTarget.cd(2); DrawPairAndSingle(hMinChi2Pair, hMinChi2Single, "Minimum #chi^{2}/NDF");
  cTarget.cd(3); pYCompatibleVsZ.SetMinimum(0.0); pYCompatibleVsZ.SetMaximum(1.0);
  pYCompatibleVsZ.SetLineColor(kBlue + 1); pYCompatibleVsZ.SetLineWidth(2);
  pYCompatibleVsZ.Draw();
  cTarget.cd(4); pThetaVsZ.SetLineColor(kBlue + 1); pThetaVsZ.SetLineWidth(2);
  SetProfileDisplayRange(pThetaVsZ, centralEArmAngleDegrees);
  pThetaVsZ.Draw();
  TLine centralTheta(targetMin, centralEArmAngleDegrees,
                     targetMax, centralEArmAngleDegrees);
  centralTheta.SetLineColor(kGray + 2); centralTheta.SetLineStyle(2);
  centralTheta.SetLineWidth(2); centralTheta.Draw("same");

  TCanvas cY("cCDetFTROIY", "CDet y and angle diagnostics", 1500, 1000);
  cY.Divide(2, 2);
  cY.cd(1); hYResidualL1VsZ.Draw("colz");
  cY.cd(2); hYResidualL2VsZ.Draw("colz");
  cY.cd(3); pOutOfPlaneVsZ.SetLineColor(kBlue + 1);
  pOutOfPlaneVsZ.SetLineWidth(2); SetProfileDisplayRange(pOutOfPlaneVsZ);
  pOutOfPlaneVsZ.Draw();
  cY.cd(4); pInPlaneVsZ.SetLineColor(kMagenta + 2);
  pInPlaneVsZ.SetLineWidth(2); SetProfileDisplayRange(pInPlaneVsZ);
  pInPlaneVsZ.Draw();

  const TString pdf = TString::Format("%s.pdf", outputStem);
  cGeometry.Print(pdf + "[");
  cGeometry.Print(pdf);
  cTarget.Print(pdf);
  cY.Print(pdf);
  cY.Print(pdf + "]");
  cGeometry.Print(TString::Format("%s_geometry.png", outputStem));
  cTarget.Print(TString::Format("%s_target_z.png", outputStem));
  cY.Print(TString::Format("%s_y_angles.png", outputStem));

  std::cout << "[CDet FTROI target-z] events with hypotheses: "
            << eventsWithHypotheses << std::endl;
  std::cout << "[CDet FTROI target-z] hypotheses: " << hypotheses << std::endl;
  std::cout << "[CDet FTROI target-z] wrote " << pdf << std::endl;
}
