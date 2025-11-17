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

const TString REPLAYED_DIR = TString(gSystem->Getenv("OUT_DIR"));// + "/rootfiles";

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

void AddRunFilesToChain(TChain *chain, const char *dir, int runnum, int onlySegment = -1) {
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

      if (onlySegment >= 0) {
        int firstSeg = -1, lastSeg = -1;
        if (!GetSegRange(fname, firstSeg, lastSeg)) {
          // If there is no _seg part, skip when filtering by segment
          continue;
        }
        // Accept if onlySegment is within [firstSeg, lastSeg]
        if (!(onlySegment >= firstSeg && onlySegment <= lastSeg)) continue;
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

int NTotBins = 200;
double TotBinLow = 1.;
double TotBinHigh = 51.;

void t_mod(int runnum = 5811, Int_t neventsr=500000, Int_t onlySegment = -1, Double_t LeMin = 0.02, Double_t LeMax = 60){
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
        AddRunFilesToChain(T, REPLAYED_DIR.Data(), runnum, onlySegment);
    }

    TTreeReader reader(T);

    TTreeReaderArray<double> RawElID   (reader, "earm.cdet.hits.TDCelemID");
    TTreeReaderArray<double> RawElLE   (reader, "earm.cdet.hits.t");
    TTreeReaderArray<double> RawElTE   (reader, "earm.cdet.hits.t_te");
    TTreeReaderArray<double> RawElTot  (reader, "earm.cdet.hits.t_tot");

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
    std::cout << "someone cooked here" << std::endl;
}//end main

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
