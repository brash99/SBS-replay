#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TMarker.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAxis.h>

#include <vector>
#include <iostream>
#include <cmath>

// ----------------------------------------------------------------------
// Minimal geometry / constants (copied from your analysis macro)
// ----------------------------------------------------------------------
static const int NumPaddles   = 16;
static const int NumBars      = 14;
static const int NumLayers    = 2;
static const int NumSides     = 2;
static const int NumModules   = 3;
static const int NumHalfModules = NumModules * NumSides * NumLayers;

static const int NumCDetPaddles   = NumHalfModules * NumBars * NumPaddles; // 2688
static const int nRef             = 4;
static const int NumRefPaddles    = 4;
static const int nTdc             = NumCDetPaddles + NumRefPaddles;        // 2704

// ----------------------------------------------------------------------
// Minimal TCDet namespace: only what we need for the event display
// ----------------------------------------------------------------------
namespace TCDet {
  // Vector sizes
  Int_t   NdataGoodX;
  Int_t   NdataGoodY;
  Int_t   NdataGoodLayer;

  // Hit positions and layer
  Double_t GoodX[nTdc*2];
  Double_t GoodY[nTdc*2];
  Double_t GoodLayer[nTdc*2];

  // ECal reconstructed position
  Double_t GoodECalX;
  Double_t GoodECalY;
  Double_t GoodECalE;

  // Hit counts
  Double_t nhits;
  Double_t ngoodhits;
  Double_t ngoodTDChits;
}

// ----------------------------------------------------------------------
// Global chain and replay directory
// ----------------------------------------------------------------------
TChain *T = nullptr;

// You can swap this to use OUT_DIR if you prefer:
// const TString REPLAYED_DIR = gSystem->Getenv("OUT_DIR");
const TString REPLAYED_DIR = "/work/brash/CDet_replay/sbs/Rootfiles";

// ----------------------------------------------------------------------
// Setup the TChain and branch addresses for the event display
// ----------------------------------------------------------------------
void CDet_SetupChainForDisplay(Int_t RunNumber1,
                               Int_t neventsr,
                               Int_t nruns)
{
  if (T) return; // already set up

  T = new TChain("T");

  TString subfile, sInFile;

  // First "main" file (same pattern as your analysis macro)
  subfile = TString::Format("cdet_%d_%d", RunNumber1, neventsr);
  sInFile = REPLAYED_DIR + "/" + subfile + ".root";
  std::cout << "[CDet display] Adding " << sInFile << std::endl;
  T->Add(sInFile);

  // Additional files
  std::cout << "[CDet display] Adding " << nruns << " files..." << std::endl;
  for (Int_t i = 1; i <= nruns; i++) {
    subfile = TString::Format(
      "cdet_%d_stream_0_2_seg0_9_firstevent1_nevent%d_%d",
      RunNumber1, neventsr, i
    );
    sInFile = REPLAYED_DIR + "/" + subfile + ".root";
    std::cout << "  " << sInFile << std::endl;
    T->Add(sInFile);
  }

  // Branch status and addresses (only what we actually need)
  T->SetBranchStatus("*", 0);

  // Enable all CDet and ECal branches (easiest)
  T->SetBranchStatus("earm.cdet.*", 1);
  T->SetBranchStatus("earm.ecal.*", 1);

  // CDet good hit positions and layer
  T->SetBranchAddress("earm.cdet.hit.xhit",  TCDet::GoodX);
  T->SetBranchAddress("earm.cdet.hit.yhit",  TCDet::GoodY);
  T->SetBranchAddress("earm.cdet.hit.layer", TCDet::GoodLayer);

  // Hit counts
  T->SetBranchAddress("earm.cdet.nhits",        &TCDet::nhits);
  T->SetBranchAddress("earm.cdet.ngoodhits",    &TCDet::ngoodhits);
  T->SetBranchAddress("earm.cdet.ngoodTDChits", &TCDet::ngoodTDChits);

  // ECal reconstructed position
  T->SetBranchAddress("earm.ecal.x", &TCDet::GoodECalX);
  T->SetBranchAddress("earm.ecal.y", &TCDet::GoodECalY);
  T->SetBranchAddress("earm.ecal.e", &TCDet::GoodECalE);

  // Ndata branches (vector sizes)
  T->SetBranchAddress("Ndata.earm.cdet.hit.xhit",  &TCDet::NdataGoodX);
  T->SetBranchAddress("Ndata.earm.cdet.hit.yhit",  &TCDet::NdataGoodY);
  T->SetBranchAddress("Ndata.earm.cdet.hit.layer", &TCDet::NdataGoodLayer);

  std::cout << "[CDet display] Chain setup complete." << std::endl;
}

// ----------------------------------------------------------------------
// Draw a single event on a canvas with 2 pads (Layer 1 / Layer 2)
// ----------------------------------------------------------------------
void CDet_DrawEvent(Long64_t iev,
                    TH2F *frameL1,
                    TH2F *frameL2,
                    TCanvas *c)
{
  T->GetEntry(iev);

  std::vector<double> xL1, yL1;
  std::vector<double> xL2, yL2;

  // Split hits into layer 1 and layer 2
  for (Int_t ih = 0; ih < TCDet::NdataGoodX; ++ih) {
    double x = TCDet::GoodX[ih];
    double y = TCDet::GoodY[ih];
    int    layer = static_cast<int>(TCDet::GoodLayer[ih]);

    if (!std::isfinite(x) || !std::isfinite(y)) continue;

    if (layer == 1) {
      xL1.push_back(x);
      yL1.push_back(y);
    } else if (layer == 2) {
      xL2.push_back(x);
      yL2.push_back(y);
    }
  }

  // ---------------- Pad 1: Layer 1 ----------------
  c->cd(1);
  frameL1->Draw("axis");

  TGraph *gL1 = nullptr;
  if (!xL1.empty()) {
    gL1 = new TGraph((Int_t)xL1.size(), xL1.data(), yL1.data());
    gL1->SetMarkerStyle(20);
    gL1->SetMarkerColor(kRed+1);
    gL1->SetMarkerSize(1.1);
    gL1->Draw("P SAME");
  }

  // ECal marker overlaid on both pads (if finite)
  TMarker *mE1 = nullptr;
  if (std::isfinite(TCDet::GoodECalX) && std::isfinite(TCDet::GoodECalY)) {
    mE1 = new TMarker(TCDet::GoodECalX, TCDet::GoodECalY, 29);
    mE1->SetMarkerColor(kGreen+2);
    mE1->SetMarkerSize(2.0);
    mE1->Draw("SAME");
  }

  TLegend *leg1 = new TLegend(0.12, 0.80, 0.55, 0.93);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.035);
  leg1->SetHeader("Layer 1", "L");
  if (gL1) leg1->AddEntry(gL1, "CDet L1 hits", "p");
  if (mE1) leg1->AddEntry(mE1, "ECal (X,Y)", "p");
  leg1->Draw();

  // Event label
  {
    TLatex lab;
    lab.SetNDC();
    lab.SetTextSize(0.035);
    lab.DrawLatex(0.15, 0.94,
      Form("Event %lld   nhits=%.0f, ngood=%.0f",
           iev, TCDet::nhits, TCDet::ngoodhits));
  }

  // ---------------- Pad 2: Layer 2 ----------------
  c->cd(2);
  frameL2->Draw("axis");

  TGraph *gL2 = nullptr;
  if (!xL2.empty()) {
    gL2 = new TGraph((Int_t)xL2.size(), xL2.data(), yL2.data());
    gL2->SetMarkerStyle(21);
    gL2->SetMarkerColor(kBlue+1);
    gL2->SetMarkerSize(1.1);
    gL2->Draw("P SAME");
  }

  TMarker *mE2 = nullptr;
  if (std::isfinite(TCDet::GoodECalX) && std::isfinite(TCDet::GoodECalY)) {
    mE2 = new TMarker(TCDet::GoodECalX, TCDet::GoodECalY, 29);
    mE2->SetMarkerColor(kGreen+2);
    mE2->SetMarkerSize(2.0);
    mE2->Draw("SAME");
  }

  TLegend *leg2 = new TLegend(0.12, 0.80, 0.55, 0.93);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.035);
  leg2->SetHeader("Layer 2", "L");
  if (gL2) leg2->AddEntry(gL2, "CDet L2 hits", "p");
  if (mE2) leg2->AddEntry(mE2, "ECal (X,Y)", "p");
  leg2->Draw();

  c->Update();
}

// ----------------------------------------------------------------------
// Main driver: side-by-side event display with simple keyboard control
// ----------------------------------------------------------------------
void CDet_EventDisplay(Int_t RunNumber1 = 5811,
                       Int_t neventsr   = 103000,
                       Int_t nruns      = 30)
{
  CDet_SetupChainForDisplay(RunNumber1, neventsr, nruns);

  Long64_t Nev = T->GetEntries();
  std::cout << "[CDet display] Tree has " << Nev << " events." << std::endl;

  if (Nev <= 0) {
    std::cout << "[CDet display] No events found!" << std::endl;
    return;
  }

  // ------------------------------------------------------------------
  // Scan to determine global XY range (both layers + ECal)
  // ------------------------------------------------------------------
  double xmin =  1e9, xmax = -1e9;
  double ymin =  1e9, ymax = -1e9;

  Long64_t nScan = std::min<Long64_t>(Nev, 5000);
  for (Long64_t iev = 0; iev < nScan; ++iev) {
    T->GetEntry(iev);

    for (Int_t ih = 0; ih < TCDet::NdataGoodX; ++ih) {
      double x = TCDet::GoodX[ih];
      double y = TCDet::GoodY[ih];
      if (!std::isfinite(x) || !std::isfinite(y)) continue;
      xmin = std::min(xmin, x);
      xmax = std::max(xmax, x);
      ymin = std::min(ymin, y);
      ymax = std::max(ymax, y);
    }

    if (std::isfinite(TCDet::GoodECalX) && std::isfinite(TCDet::GoodECalY)) {
      xmin = std::min(xmin, TCDet::GoodECalX);
      xmax = std::max(xmax, TCDet::GoodECalX);
      ymin = std::min(ymin, TCDet::GoodECalY);
      ymax = std::max(ymax, TCDet::GoodECalY);
    }
  }

  // Fallback if nothing sensible found
  if (!(xmin < xmax) || !(ymin < ymax)) {
    xmin = -1.0; xmax = 1.0;
    ymin = -1.0; ymax = 1.0;
  } else {
    double dx = xmax - xmin;
    double dy = ymax - ymin;
    xmin -= 0.1 * dx;  xmax += 0.1 * dx;
    ymin -= 0.1 * dy;  ymax += 0.1 * dy;
  }

  // ------------------------------------------------------------------
  // Canvas + frames
  // ------------------------------------------------------------------
  TCanvas *c = new TCanvas("cCDetEvent",
                           "CDet Event Display (Layer 1 / Layer 2)",
                           1200, 600);
  c->Divide(2, 1);

  TH2F *frameL1 = new TH2F("hCDetL1Frame",
                           "CDet Layer 1;X;Y",
                           100, xmin, xmax,
                           100, ymin, ymax);
  TH2F *frameL2 = new TH2F("hCDetL2Frame",
                           "CDet Layer 2;X;Y",
                           100, xmin, xmax,
                           100, ymin, ymax);

  // ------------------------------------------------------------------
  // Interactive loop
  // ------------------------------------------------------------------
  Long64_t iev = 0;
  char cmd;

  std::cout << "Commands: n = next, p = previous, g = goto, q = quit" << std::endl;

  // Draw first event
  CDet_DrawEvent(iev, frameL1, frameL2, c);

  while (true) {
    std::cout << "Event " << iev << " / " << (Nev-1) << "  -> command? ";
    std::cin >> cmd;

    if (!std::cin) {
      std::cout << "\n[stdin EOF/bad] Exiting display loop." << std::endl;
      break;
    }

    if      (cmd == 'q') break;
    else if (cmd == 'n') { if (iev < Nev-1) ++iev; }
    else if (cmd == 'p') { if (iev > 0)    --iev; }
    else if (cmd == 'g') {
      Long64_t target;
      std::cout << "  Go to event: ";
      std::cin >> target;
      if (!std::cin) {
        std::cin.clear();
        std::cin.ignore(9999, '\n');
        std::cout << "  [input error]" << std::endl;
      } else if (target >= 0 && target < Nev) {
        iev = target;
      } else {
        std::cout << "  [invalid event number]" << std::endl;
      }
    } else {
      std::cout << "  [unknown command: " << cmd << "]" << std::endl;
      continue;
    }

    CDet_DrawEvent(iev, frameL1, frameL2, c);
  }

  std::cout << "[CDet display] Finished." << std::endl;
}
