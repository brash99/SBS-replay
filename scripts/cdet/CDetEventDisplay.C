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
  Int_t   NdataGoodZ;
  Int_t   NdataGoodLayer;
  Int_t   NdataGoodElID;

  // Hit positions and layer
  Double_t GoodX[nTdc*2];
  Double_t GoodY[nTdc*2];
  Double_t GoodZ[nTdc*2];
  Double_t GoodLayer[nTdc*2];
  Double_t GoodElID[nTdc*2];

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
                               Int_t nruns,
			       Int_t seg_start,
			       Int_t seg_end)
{
  if (T) return; // already set up

  T = new TChain("T");

  TString subfile, sInFile;

  // First "main" file (same pattern as your analysis macro)
  subfile = TString::Format("cdet_%d_stream_0_2_seg%d_%d_firstevent1_nevent%d", RunNumber1, seg_start, seg_end, neventsr);
  sInFile = REPLAYED_DIR + "/" + subfile + ".root";
  std::cout << "[CDet display] Adding " << sInFile << std::endl;
  T->Add(sInFile);

  // Additional files
  std::cout << "[CDet display] Adding " << nruns << " files..." << std::endl;
  for (Int_t i = 1; i <= nruns; i++) {
    subfile = TString::Format(
      "cdet_%d_stream_0_2_seg%d_%d_firstevent1_nevent%d_%d",
      RunNumber1, seg_start, seg_end, neventsr, i
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
  T->SetBranchAddress("earm.cdet.hit.zhit",  TCDet::GoodZ);
  T->SetBranchAddress("earm.cdet.hit.layer", TCDet::GoodLayer);
  T->SetBranchAddress("earm.cdet.hit.pmtnum", TCDet::GoodElID);

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
  T->SetBranchAddress("Ndata.earm.cdet.hit.zhit",  &TCDet::NdataGoodZ);
  T->SetBranchAddress("Ndata.earm.cdet.hit.layer", &TCDet::NdataGoodLayer);
  T->SetBranchAddress("Ndata.earm.cdet.hit.pmtnum", &TCDet::NdataGoodElID);

  std::cout << "[CDet display] Chain setup complete." << std::endl;
}

// ----------------------------------------------------------------------
// Draw a single event on a canvas with 2 pads (Layer 1 / Layer 2)
// ----------------------------------------------------------------------
void CDet_DrawEvent(Long64_t iev,
                    TH2F *frameL1,
                    TH2F *frameL2,
                    TCanvas *c,
		    double XOffset,
		    double XDiffCut)
{
  T->GetEntry(iev);

  std::vector<double> xL1, yL1;
  std::vector<double> xL2, yL2;
	
  double ecal_dist = 6.6;
  double cdet_dist_offset = 2.0;
  double cdet_y_half_length = 9.30;
  double zlayer1 = 0;
  double zlayer2 = 0;

  // Split hits into layer 1 and layer 2
  for (Int_t ih = 0; ih < TCDet::NdataGoodX; ++ih) {
        bool good_ecal_reconstruction = TCDet::GoodECalY > -1.2 && TCDet::GoodECalY < 1.2 &&
		TCDet::GoodECalX > -1.5 && TCDet::GoodECalX < 1.5 &&
		TCDet::GoodECalX != 0.00 && TCDet::GoodECalY != 0.00;
	bool good_cdet_X = TCDet::GoodX[ih] < 998;
	bool good_ecal_diff_x = (TCDet::GoodX[ih]-(TCDet::GoodECalX*(TCDet::GoodZ[ih]-cdet_dist_offset)/ecal_dist)-XOffset) <= XDiffCut && 
					(TCDet::GoodX[ih]-(TCDet::GoodECalX*(TCDet::GoodZ[ih]-cdet_dist_offset)/ecal_dist)-XOffset) >= -1.0*XDiffCut; 
	bool good_ecal_diff_y = (TCDet::GoodY[ih]-(TCDet::GoodECalY*(TCDet::GoodZ[ih]-cdet_dist_offset)/ecal_dist)) <= cdet_y_half_length && 
					(TCDet::GoodY[ih]-(TCDet::GoodECalY*(TCDet::GoodZ[ih]-cdet_dist_offset)/ecal_dist)) >= -1.0*cdet_y_half_length; 
	bool good_CDet_event = good_ecal_reconstruction && good_ecal_diff_x && good_ecal_diff_y && good_cdet_X;
   
	//cout << "good_ecal_reconstruction = " << good_ecal_reconstruction <<
        //		" good_cdet_X = " << good_cdet_X <<
	//	" good_ecal_diff_x = " << good_ecal_diff_x <<
	//	" good_ecal_diff_y = " << good_ecal_diff_y <<
	//	" good_CDet_event = " << good_CDet_event << endl;

        double x = TCDet::GoodX[ih];
        double y = TCDet::GoodY[ih];
        int    layer = static_cast<int>(TCDet::GoodLayer[ih]);
    

        if (!std::isfinite(x) || !std::isfinite(y)) continue;

	if (TCDet::GoodElID[ih] < 1344) {
		layer = 1;
		zlayer1 = TCDet::GoodZ[ih];
	} else {
		layer = 2;
		zlayer2 = TCDet::GoodZ[ih];
	}

	if (good_CDet_event) {
		//cout << "cdet hit " << ih << " layer = " << layer << " x = " << x << " y = " << y << endl;
        	if (layer == 1) {
          		xL1.push_back(x);
          		yL1.push_back(y);
        	} else if (layer == 2) {
          		xL2.push_back(x);
          		yL2.push_back(y);
        	}
	}
  }

  // ---------------- Pad 1: Layer 1 ----------------
  c->cd(1);
  //frameL1->Draw("axis");


  // draw rectangles for each hit
  double half_dx = 0.005;   // total width = 0.01
  double half_dy = 0.25;    // total height = 0.50
  //cout << "xL1 size = " << xL1.size() << endl;

    for (size_t i = 0; i < xL1.size(); ++i) {
    double x = xL1[i];
    double y = yL1[i];
    //cout << "cdet layer 1 hit " << i+1 << " x = " << x << " y = " << y << endl;
    TBox *box = new TBox(x - half_dx, y - half_dy,
                       x + half_dx, y + half_dy);
    box->SetFillColorAlpha(kRed+1, 0.35);   // semi-transparent red fill
    box->SetLineColor(kRed+2);
    box->SetLineWidth(2);
    box->Draw("same");
  }

  // ECal marker overlaid on both pads (if finite)
  TMarker *mE1 = nullptr;
  if (std::isfinite(TCDet::GoodECalX) && std::isfinite(TCDet::GoodECalY)) {
    mE1 = new TMarker(TCDet::GoodECalX*(zlayer1-cdet_dist_offset)/ecal_dist, TCDet::GoodECalY*(zlayer1-cdet_dist_offset)/ecal_dist, 29);
    mE1->SetMarkerColor(kGreen+2);
    mE1->SetMarkerSize(2.0);
    mE1->Draw("SAME");
  }
/*
  TLegend *leg1 = new TLegend(0.12, 0.80, 0.55, 0.93);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.035);
  leg1->SetHeader("Layer 1", "L");
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
*/

  // ---------------- Pad 2: Layer 2 ----------------
  c->cd(2);
  //frameL2->Draw("axis");


  for (size_t i = 0; i < xL2.size(); ++i) {
    double x = xL2[i];
    double y = yL2[i];
    //cout << "cdet layer 2 hit " << i+1 << " x = " << x << " y = " << y << endl;
    TBox *box = new TBox(x - half_dx, y - half_dy,
                       x + half_dx, y + half_dy);
    box->SetFillColorAlpha(kRed+1, 0.35);   // semi-transparent red fill
    box->SetLineColor(kRed+2);
    box->SetLineWidth(2);
    box->Draw("same");
    
  }


  TMarker *mE2 = nullptr;
  if (std::isfinite(TCDet::GoodECalX) && std::isfinite(TCDet::GoodECalY)) {
    mE2 = new TMarker(TCDet::GoodECalX*(zlayer1-cdet_dist_offset)/ecal_dist, TCDet::GoodECalY*(zlayer1-cdet_dist_offset)/ecal_dist, 29);
    mE2->SetMarkerColor(kGreen+2);
    mE2->SetMarkerSize(2.0);
    mE2->Draw("SAME");
  }
/*
  TLegend *leg2 = new TLegend(0.12, 0.80, 0.55, 0.93);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.035);
  leg2->SetHeader("Layer 2", "L");
  if (mE2) leg2->AddEntry(mE2, "ECal (X,Y)", "p");
  leg2->Draw();
*/
  c->Update();
}

// ----------------------------------------------------------------------
// Main driver: side-by-side event display with simple keyboard control
// ----------------------------------------------------------------------
void CDet_EventDisplay(Int_t RunNumber1  = 5811,
                       Int_t neventsr    = 300000,
		       Double_t XOffset  = 0.02,
		       Double_t XDiffCut = 0.04,
                       Int_t nruns       = 30,
		       Int_t seg_start   = 0,
		       Int_t seg_end     = 14)
{
  CDet_SetupChainForDisplay(RunNumber1, neventsr, nruns, seg_start,seg_end);

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
      xmax = -1.0*xmin;
      ymax = std::max(ymax, y);
      ymin = -1.0*ymax;
      //std::cout << "xmin = " << xmin << " xmax = " << xmax << " ymin = " << ymin << " ymax = " << ymax << std::endl;
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
  CDet_DrawEvent(iev, frameL1, frameL2, c, XOffset, XDiffCut);

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
    
    c->cd(1);
    frameL1->Draw("axis");
    c->cd(2);
    frameL2->Draw("axis");

    for (int ii=iev; ii<iev+1000; ii++) {
      cout << "Event = " << ii << endl;
      CDet_DrawEvent(ii, frameL1, frameL2, c, XOffset, XDiffCut);
    }

  }

  std::cout << "[CDet display] Finished." << std::endl;
}
