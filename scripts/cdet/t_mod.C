#include <TROOT.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>      // for sscanf
#include <algorithm>   // for std::sort
#include <TMath.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TLatex.h>
#include <vector>

using namespace std;

const TString REPLAYED_DIR = TString(gSystem->Getenv("OUT_DIR")) + "/rootfiles";

// const TString ANALYSED_DIR = gSystem->Getenv("ANALYSED_DIR");
//const TString REPLAYED_DIR = "/volatile/halla/sbs/btspaude/cdet/rootfiles";
const TString ANALYSED_DIR = "/work/halla/sbs/btspaude/sbs/Rootfiles/cdetFiles";

// Parse the "segX_Y" part: returns true and fills firstSeg/lastSeg if found.
bool GetSegRange(const TString& fname, int& firstSeg, int& lastSeg) {
  // Find "_seg"
  Ssiz_t pos = fname.Index("_seg");
  if (pos == kNPOS) return false;

  // Tail looks like "9_9.root" or "9_9_1.root"
  TString tail = fname(pos + 4, fname.Length() - (pos + 4));

  // Extract first two ints; ignore any further suffix
  int a = -1, b = -1;
  if (sscanf(tail.Data(), "%d_%d", &a, &b) == 2) {
    firstSeg = a; lastSeg = b;
    return true;
  }
  return false;
}

void AddRunFilesToChain(TChain *chain, const char *dir, int runnum, int segMin = -1, int segMax = -1) {
  TString prefix = dir;
  std::vector<TString> runfiles;

  TSystemDirectory directory(prefix, prefix);
  TList *files = directory.GetListOfFiles();

  if (files) {
    TIter next(files);
    TSystemFile *f;
    while ((f = (TSystemFile*) next())) {
      if (f->IsDirectory()) continue; // skip dirs like "." and ".."

      TString fname = f->GetName();
      if (!fname.BeginsWith(Form("cdet_%d_", runnum))) continue;
      if (!fname.EndsWith(".root")) continue;

      // Range filtering enabled only if segMin/segMax are set
      if (segMin >= 0 || segMax >= 0) {
        if (segMin < 0) segMin = segMax;
        if (segMax < 0) segMax = segMin;
        if (segMin > segMax) std::swap(segMin, segMax);

        int firstSeg = -1, lastSeg = -1;
        if (!GetSegRange(fname, firstSeg, lastSeg)) continue;

        // accept if [firstSeg,lastSeg] overlaps [segMin,segMax]
        if (lastSeg < segMin || firstSeg > segMax) continue;
      }

      runfiles.push_back(prefix + "/" + fname);
    }
  }

  std::sort(runfiles.begin(), runfiles.end());

  std::cout << "Adding " << runfiles.size() << " files for run " << runnum << "...\n";
  for (auto &file : runfiles) {
    std::cout << "  " << file << "\n";
    chain->Add(file);
  }
}

std::vector<int> getLocation(int pixelID) {
  // Check valid range
  if (pixelID < 0 || pixelID > 2687) {
      std::cerr << "Error: pixelID must be in the range 0 to 2687.\n";
      return {};  // return empty vector to signal error
  }

  int layerNum      = pixelID / 1344; //1344 pixels per layer
  pixelID               %= 1344;

  int sideNum       = pixelID / 672;  //672 pixels per side
  pixelID               %= 672;

  int submoduleNum  = pixelID / 224; //224 pixels per side of module
  pixelID               %= 224;

  int pmtNum        = pixelID / 16; //16 pixels per bar
  pixelID               %= 16;
  int pixelNum      = pixelID % 16;

  return {layerNum, sideNum, submoduleNum, pmtNum, pixelNum};
}


TChain *T = 0;

static const double TDC_calib_to_ns = 0.01;

double TDCBinLow;
double TDCBinHigh;
int NTDCBins;
Double_t RefNTDCBins;
Double_t RefLeMin;
Double_t RefLeMax;

std::vector<double> vRefRawLe;
std::vector<double> vRefRawTe;
std::vector<double> vRefRawTot;
std::vector<int>    vRefRawPMT;
std::vector<double> vAllRawLe;
std::vector<double> vAllRawLeNoRef;
std::vector<int>    vAllRawPMT;
std::vector<double> vT_mod;
std::vector<double> vAllRefForLe;

std::vector<std::vector<double>> vCDetPaddleRawTot;
std::vector<std::vector<double>> vCDetPaddleCutTot;

std::vector<int> rawRate(2688, 0); 
int rateEvTrack = 0;
std::vector<double> chanRates(2688,0);
std::vector<int> cutRate(2688, 0); 
int cutRateEvTrack = 0;
std::vector<double> cutChanRates(2688,0);

int NTotBins = 200;
double TotBinLow = 1.;
double TotBinHigh = 51.;

void t_mod(int runnum = 5811, Int_t neventsr=500000, Int_t minSeg = -1, Int_t maxSeg = -1, Double_t LeMin = 0.02, Double_t LeMax = 60){
    RefLeMin = 0.0;
    RefLeMax = 252.0;
    RefNTDCBins = (RefLeMax-RefLeMin)/4;
    Double_t RefTotMin = 1.0;
    Double_t RefTotMax = 252.0;

    NTDCBins = 2*(LeMax-LeMin)/.0160167; // 4 ns is the trigger time, 0.018 ns is the expected time resolution, if we use a reference TDC ? 
                    // 4 ns resolution is the best we can hope for, I think, using only the module trigger time.
    TDCBinLow = LeMin;
    TDCBinHigh = LeMax;

    if (!T) {
        T = new TChain("T");
        //int onlySegment = -1; // set to >=0 to pick just one
        AddRunFilesToChain(T, REPLAYED_DIR.Data(), runnum, minSeg, maxSeg);
    }

    TTreeReader reader(T);

    TTreeReaderArray<double> RawElID   (reader, "earm.cdet.hits.TDCelemID");
    TTreeReaderArray<double> RawElLE   (reader, "earm.cdet.hits.t");
    TTreeReaderArray<double> RawElTE   (reader, "earm.cdet.hits.t_te");
    TTreeReaderArray<double> RawElTot  (reader, "earm.cdet.hits.t_tot");
    
    TTreeReaderArray<double> GoodElID  (reader, "earm.cdet.hit.pmtnum");
    TTreeReaderArray<double> GoodElLE  (reader, "earm.cdet.hit.tdc_le");
    TTreeReaderArray<double> GoodElTE  (reader, "earm.cdet.hit.tdc_te");
    TTreeReaderArray<double> GoodElTot (reader, "earm.cdet.hit.tdc_tot");

    vCDetPaddleRawTot.assign(2688, std::vector<double>{});
    vCDetPaddleCutTot.assign(2688, std::vector<double>{});

    Int_t Nev = T->GetEntries();
    cout << "N entries in tree is " << Nev << endl;
    Int_t NEventsAnalysis;// = Nev;
    if(neventsr==-1) NEventsAnalysis = Nev;
    else NEventsAnalysis = neventsr;
    cout << "Running analysis for " << NEventsAnalysis << " events" << endl;

    //================================================================= Event Loop
    // variables outside event loop
    Int_t EventCounter = 0;
    cout << "Starting Event Loop" << endl;

    // event loop start
    Int_t event = 0;
    while(reader.Next()){
        event++;
        event = event - 1;
        EventCounter++;
        // Only stop early if neventsr > 0
        if (neventsr > 0 && EventCounter > neventsr) {
            break;
        }
        
        if (EventCounter % 1000 == 0) {
        cout << EventCounter << "/" << NEventsAnalysis << endl;
        }

        // First pass through hits:  purpose is to get reference LE TDC Value for this event
      
      double event_ref_tdc = 0.0;
      double ref_int = 0;
      double ref_corr = 0;
      // std::cout << "RawElID Size = " << " " << RawElID.GetSize() << std::endl;
      // for (auto val : RawElLE) {
      //   std::cout << " rawLE array= "<< " " << val <<std::endl; 
      // }
        for(Int_t el=0; el<RawElID.GetSize(); el++) {
          if ((Int_t)RawElID[el] == 2696) {  // only look at ref PMT 
            bool good_ref_le_time = RawElLE[el] > 0.0/TDC_calib_to_ns && RawElLE[el] <= 252.0/TDC_calib_to_ns;
            bool good_ref_event = good_ref_le_time;
            if ( good_ref_event ) {
              //std::cout << "event = " << " " << EventCounter << " " << "ref time = " << RawElLE[el] << std::endl;
              //std::cout << "ref Tot = " << RawElTot[el] << std::endl;
              if ( (Int_t)RawElID[el] == 2696 && (Int_t)RawElLE[el] > 0 ) {
                  vRefRawLe.push_back(RawElLE[el] * TDC_calib_to_ns);
                  vRefRawTe.push_back(RawElTE[el] * TDC_calib_to_ns);
                  vRefRawTot.push_back(RawElTot[el] * TDC_calib_to_ns);
                  vRefRawPMT.push_back((int)RawElID[el]);

                  event_ref_tdc = RawElLE[el]*TDC_calib_to_ns - 48;
                  ref_int = std::floor(event_ref_tdc);
                  ref_corr = event_ref_tdc - ref_int;
              }
            }
          }
        }// end ref TDC loop
        int nID  = RawElID.GetSize();
        int nTot = RawElTot.GetSize();
        int n    = std::min(nID, nTot);  
        rateEvTrack++;
        bool eventHasCutHit = false;
        for(Int_t el = 0; el < n; el++) {
          int idx = RawElID[el];
          if (0 <= idx && idx < 2688) {
            double tot_ns = RawElTot[el]*TDC_calib_to_ns;
            rawRate[idx]++;
            vCDetPaddleRawTot[idx].push_back(tot_ns);
            if (tot_ns >= 10.0){
              cutRate[idx]++;
              vCDetPaddleCutTot[idx].push_back(tot_ns);
              eventHasCutHit = true;
            }
          }
        }
        if (eventHasCutHit) cutRateEvTrack++;

        for(Int_t el=0; el<RawElID.GetSize(); el++){

        bool good_raw_le_time = RawElLE[el] >= LeMin/TDC_calib_to_ns && RawElLE[el] <= LeMax/TDC_calib_to_ns;
        bool good_raw_event = good_raw_le_time;

        //if ((Int_t)RawElID[el] > 1000) cout << "el = " << el << " Hit ID = " << (Int_t)RawElID[el] << "    TDC = " << RawElLE[el]*TDC_calib_to_ns << endl;
        //cout << "Raw ID = " << RawElID[el] << " raw le = " << RawElLE[el] << " raw te = " << RawElTE[el] << " raw tot = " << RawElTot[el] << endl;
            if ( good_raw_event ) {
                if ( (Int_t)RawElID[el] < 2688 ) {
                    //fill all hits vectors
                    double t_ns = RawElLE[el]*TDC_calib_to_ns - event_ref_tdc;//ref_corr;
                    double wrapped = fmod(t_ns,2.0);

                    vAllRawLe.push_back(RawElLE[el]*TDC_calib_to_ns - event_ref_tdc);//ref_corr);
                    vAllRawLeNoRef.push_back(RawElLE[el]*TDC_calib_to_ns);//ref_corr);
                    vAllRawPMT.push_back(RawElID[el]);
                    vT_mod.push_back(wrapped);
                    vAllRefForLe.push_back(event_ref_tdc);
                }
            }
        }// all raw tdc hit loop
    }//end event loop 
    std::cout << "nevents = " << rateEvTrack << std::endl;
    for (int i = 0; i < 2688; i++){
      chanRates[i] = (double)rawRate[i] / (rateEvTrack);
      cutChanRates[i] = (double)cutRate[i] / cutRateEvTrack;
      //std::cout << "triggered Rate in Pixel " << 417 + i << " = " << chanRates << " & with time window Rate = " << chanRates / winWidth <<std::endl;
    }

    std::cout << "someone cooked here - Walter White" << std::endl;
}//end main

void plotSingleTot(int pixel_base = 0, double width = 1, double totMin=1, double totMax=80){
  TH1::AddDirectory(kFALSE);
  if (pixel_base % 16 != 0) {
    Error("plotSingleTot", "pixel_base = %d is not a multiple of 16", pixel_base);
    return;
  }

  const int nPlots = 16;
  int TDCBinNum = (int)((totMax-totMin)/width);

  // Decode pixel → {layer, side, submodule, pmt, pixel}
  auto info = getLocation(pixel_base);

  int layer     = info[0] + 1;  // display as 1-based
  int side      = info[1];      // 0=L, 1=R
  int submodule = info[2] + 1;
  int bar       = info[3] + 1;

  TString sideStr = (side == 0) ? "L" : "R";

  TString canvasTitle = Form("Layer %d | %s | Module %d | Bar %d", layer, sideStr.Data(), submodule, bar);
  // Canvas with 4x4 pads
  TCanvas* cTot = new TCanvas("cTot", canvasTitle, 1200, 1000);
  cTot->Divide(4, 4, 0.001, 0.001);

  // Histogram array
  TH1D* hTot[nPlots];

  for (int i = 0; i < nPlots; i++) {

    int pixel = pixel_base + i;

    TString hname  = Form("hTot_pix%d", pixel);
    TString htitle = Form("Pixel %d;TOT (ns);Counts", pixel);

    hTot[i] = new TH1D(hname, htitle, TDCBinNum, totMin, totMax);

    // Fill histogram
    for (const auto& x : vCDetPaddleRawTot[pixel]) {
      hTot[i]->Fill(x);
    }

    // Draw
    cTot->cd(i + 1);
    hTot[i]->Draw();
  }

  cTot->Update();
}

void getRate(int pixel, bool cut = false){
  if (!cut) std:: cout << "Rate in Pixel " << pixel << " = " << chanRates[pixel] << std::endl;
  if (cut) std:: cout << "Rate in Pixel " << pixel << " (with Tot Cut) = " << cutChanRates[pixel] << std::endl;
}

void plotRateVsID(bool raw = true){
  TH1::AddDirectory(kFALSE);
    // --- constants ---
  const int NCHAN_TOTAL = 2688;
  const int NCHAN_LAYER = 1344;
  const int NCHAN_SIDE  = 672;
  const int NMOD        = 3;
  const int NCHAN_MOD   = NCHAN_SIDE / NMOD; // 224

  // Helper: create one segment histogram and fill from chanRates
  auto MakeRateHist = [&](const char* hname,
                          const char* htitle,
                          int idStart, int idEnd,
                          double yMax = 1.5) -> TH1D* {
    const int nbins = idEnd - idStart + 1; // inclusive
    TH1D* h = new TH1D(hname, htitle, nbins, idStart, idEnd + 1); // [start, end+1)
    h->SetStats(0);
    h->SetMinimum(0.0);
    h->SetMaximum(yMax);

    for (int id = idStart; id <= idEnd; id++) {
      const int bin = h->FindBin(id);
      if (raw){
        h->SetBinContent(bin, chanRates[id]);
      }
      if (!raw){
        h->SetBinContent(bin, cutChanRates[id]);
      }
    }
    return h;
  };

  struct Seg { int layer; int mod; const char* side; int start; int end; };

  std::vector<Seg> segs;
  segs.reserve(12);

  auto AddLayerSegs_LeftThenRight = [&](int layer, int base) {
    const int L0 = base + 0;
    const int R0 = base + NCHAN_SIDE;

    // Left side: M1, M2, M3
    for (int m = 0; m < NMOD; m++) {
      int s = L0 + m*NCHAN_MOD;
      int e = s + NCHAN_MOD - 1;
      segs.push_back({layer, m+1, "L", s, e});
    }
    // Right side: M1, M2, M3
    for (int m = 0; m < NMOD; m++) {
      int s = R0 + m*NCHAN_MOD;
      int e = s + NCHAN_MOD - 1;
      segs.push_back({layer, m+1, "R", s, e});
    }
  };

  // Layer 1: IDs 0..1343
  AddLayerSegs_LeftThenRight(1, 0);
  // Layer 2: IDs 1344..2687
  AddLayerSegs_LeftThenRight(2, 1344);

  // --- build histograms ---
  TH1D* hRateSeg[12] = {nullptr};

  for (int i = 0; i < 12; i++) {
    const auto& s = segs[i];
    if (raw) {
      TString name  = Form("hRateVsIDL%dM%d%s", s.layer, s.mod, s.side);
      TString title = Form("CDet L%d %s M%d Rate vs Paddle ID;Paddle ID;Rate",
                          s.layer, s.side, s.mod);
      hRateSeg[i] = MakeRateHist(name.Data(), title.Data(), s.start, s.end, 1.5);
    }
    if (!raw){
      TString name  = Form("hRateVsIDL%dM%d%s", s.layer, s.mod, s.side);
      TString title = Form("CDet L%d %s M%d Rate w/Cut vs Paddle ID;Paddle ID;Rate",
                          s.layer, s.side, s.mod);
      hRateSeg[i] = MakeRateHist(name.Data(), title.Data(), s.start, s.end, 1.5);
    }
  }

  // --- Draw: Layer 1 canvas (Left M1-3 then Right M1-3) ---
  TCanvas* cRateL1 = new TCanvas("cRateL1", "CDet Rate vs ID (Layer 1)", 1400, 800);
  cRateL1->Divide(3,2); // top row: left M1-3, bottom row: right M1-3

  int pad = 1;
  for (int i = 0; i < 12; i++) {
    if (segs[i].layer != 1) continue;
    cRateL1->cd(pad++);
    hRateSeg[i]->Draw("HIST");
  }

  // --- Draw: Layer 2 canvas (Left M1-3 then Right M1-3) ---
  TCanvas* cRateL2 = new TCanvas("cRateL2", "CDet Rate vs ID (Layer 2)", 1400, 800);
  cRateL2->Divide(3,2);

  pad = 1;
  for (int i = 0; i < 12; i++) {
    if (segs[i].layer != 2) continue;
    cRateL2->cd(pad++);
    hRateSeg[i]->Draw("HIST");
  }
}

TCanvas *plotAllTDC(double TDCBinLow, double TDCBinHigh){
    //define histograms
    TH1F *hAllRawLe = new TH1F(TString::Format("hRawLe"),
            TString::Format("hRawLe"),
            NTDCBins, TDCBinLow, TDCBinHigh);
    TH1F *hT_mod = new TH1F("hT_mod","hT_mod",
            200, 0, 2.5);
    //fill necessary histograms from vectors
    for (double x : vAllRawLe) hAllRawLe->Fill(x);
    for (double x : vT_mod) hT_mod->Fill(x);

    TCanvas *c = new TCanvas("c", "Tmod (0,2)",800,800);
    // c->Divide(2,1);

    c->cd(1);
    hAllRawLe->Draw();

    // c->cd(2);
    // hT_mod->Draw("Fill");
    return c;
}

TCanvas *plotNoRef(double TDCBinLow, double TDCBinHigh){
    //define histograms
    TH1F *hAllRawLeNoRef = new TH1F("hAllRawLeNoRef",
            "hAllRawLeNoRef",
            NTDCBins, TDCBinLow, TDCBinHigh);
    //fill necessary histograms from vectors
    for (double x : vAllRawLeNoRef) hAllRawLeNoRef->Fill(x);

    TCanvas *c = new TCanvas("c", "LE no Ref",800,800);
    hAllRawLeNoRef->Draw();
    return c;
}

TCanvas *plotRef(Double_t cutTotMin = 90, Double_t cutTotMax = 150, Int_t RefNTDCBins = 252, Double_t RefLeMin = 0.0, Double_t RefLeMax = 52.0, Int_t RefNTotBins= 252, Double_t RefTotMin = 1, Double_t RefTotMax = 250 ){
    //define histograms
    TH1F *hRefLe = new TH1F("hRefLe",
            "hRefLe",
            RefNTDCBins, RefLeMin, RefLeMax);
    TH1F *hRefTe = new TH1F("hRefTe",
            "hRefTe",
            RefNTDCBins, RefLeMin+RefTotMin, RefLeMax+RefTotMax);
    TH1F *hRefTot = new TH1F("hRefTot",
            "hRefTot",
            RefNTotBins, RefTotMin, RefTotMax);
  
    size_t N = vRefRawLe.size();
    for (size_t i = 0; i < N; i++){
      if (vRefRawTot[i] >= cutTotMin && vRefRawTot[i] <= cutTotMax) {
        hRefLe->Fill(vRefRawLe[i]);
        hRefTe->Fill(vRefRawTe[i]);
        hRefTot->Fill(vRefRawTot[i]);
      }
    }

    TCanvas *c = new TCanvas("c", "Ref Signal",800,800);
    c->Divide(1,3);

    c->cd(1);
    //gPad->SetLogy();
    hRefLe->Draw("HIST");

    c->cd(2);
    //gPad->SetLogy();
    hRefTe->Draw("HIST");

    c->cd(3);
    //gPad->SetLogy();
    hRefTot->Draw("HIST");

    return c;
}

TCanvas *plotRefvsLE(){
  TH2F *hRefVsLE = new TH2F("hRefVsLE", "hRefvsLE;LE (ns);Ref time (ns)",NTDCBins,TDCBinLow,TDCBinHigh,RefNTDCBins, RefLeMin, RefLeMax);

  for (size_t i = 0; i < vAllRawLe.size(); ++i){
    hRefVsLE->Fill(vAllRawLe[i],vAllRefForLe[i]);
  }
  TCanvas* c1 = new TCanvas("c2", "Ref vs LE", 800,800);
  hRefVsLE->Draw("COLZ");
  return c1;
}
