#include <TROOT.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <iostream>
#include <fstream>
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

static const int TDCmult_cut = 100;
static const double xcut = 998.0;
//static const int nhitcutlow = 2;
//static const int nhitcuthigh = 20;
static const double TDC_calib_to_ns = 0.01;
static const double HotChannelRatio = .01;

static const int NumPaddles = 16;
static const int NumBars = 14;
static const int NumLayers = 2;
static const int NumSides = 2;
static const int NumModules = 3;
static const int NumHalfModules = NumModules*NumSides*NumLayers;

static const int NumCDetPaddles = NumHalfModules*NumBars*NumPaddles;
static const int nRef = 2;
static const int NumRefPaddles = NumPaddles*nRef;
static const int nTdc = NumCDetPaddles+NumRefPaddles;

static const int nBarsADC = 0;
static const double ADCCUT = 150.;//100.0;

const TString REPLAYED_DIR = "/adaqfs/home/a-onl/sbs/Rootfiles";
const TString ANALYSED_DIR = "/adaqfs/home/a-onl/sbs/Rootfiles/cdetFiles/cdet_histfiles";

// // for local analysis at uog (please leave in comments)
// TString REPLAYED_DIR = "/w/work0/home/rachel/HallA/BB_Hodo/FallRun2021/Replayed";
//TString REPLAYED_DIR = "/w/work2/jlab/halla/sbs_hodo/Rootfiles";
//TString ANALYSED_DIR = "/w/work2/jlab/halla/sbs_hodo/Rootfiles/bbhodo_hist";
//TString ANALYSED_DIR = "/w/work0/home/rachel/HallA/BB_Hodo/FallRun2021/Analysed";

namespace TCDet {
  Int_t NdataMult;
  Double_t TDCmult[nTdc*2];

  Int_t NdataRawElID;
  Double_t RawElID[nTdc*2];
  Int_t NdataRawElLE;
  Double_t RawElLE[nTdc*2];
  Int_t NdataRawElTE;
  Double_t RawElTE[nTdc*2];
  Int_t NdataRawElTot;
  Double_t RawElTot[nTdc*2];
  
  Int_t NdataGoodRow;
  Double_t GoodRow[nTdc*2];
  Int_t NdataGoodCol;
  Double_t GoodCol[nTdc*2];
  Int_t NdataGoodLayer;
  Double_t GoodLayer[nTdc*2];

  Int_t NdataGoodElID;
  Double_t GoodElID[nTdc*2];
  Int_t NdataGoodElLE;
  Double_t GoodElLE[nTdc*2];
  Int_t NdataGoodElTE;
  Double_t GoodElTE[nTdc*2];
  Int_t NdataGoodElTot;
  Double_t GoodElTot[nTdc*2];

  Int_t NdataGoodX;
  Double_t GoodX[nTdc*2];
  Int_t NdataGoodY;
  Double_t GoodY[nTdc*2];
  Int_t NdataGoodZ;
  Double_t GoodZ[nTdc*2];

  Double_t GoodECalX;
  Double_t GoodECalY;
  Double_t nhits;
  Double_t ngoodhits;
  Double_t ngoodTDChits;

  Double_t ngoodTDChits_paddles[nTdc*2];
  Double_t ngoodTDCpaddles;

};

TChain *T = 0;
  
//===================================================== Histogram Declarations
// number of histo bins
const int NTotBins = 200;
const double TotBinLow = 1.;
const double TotBinHigh = 201.;

//const int num_bad = 0;
const int num_bad = 12;

//const int bad_channels[] = {
//	45, 207, 604, 728, 1174, 1488, 1684, 1685, 1712, 1726, 1912, 2131, 2132, 
//	2410, 2411, 2417, 2425, 2426, 2427, 2430, 2466, 2470, 2471, 2472, 
//	2673, 2675
//};

const int bad_channels[] = {
	604, 673, 674, 675, 676, 677, 678, 679, 680, 681, 682, 728
};
	

/*
const int bad_channels[] = {
		45,   224, 226, 234, 241, 284, 240, 604, 690, 698, 702,
		864, 1174,1437,1488,1726,1912,2131,2132,2213,2214,
		2411,2417,2425,2426,2427,2434,2436,2437,2438,2448, 2457,
		2466,2467,2468,2469,2470,2471,2472,2473,2474,2475,
		2476,2477,2478,2479,2502,2517,2518,2519,2520,2522,
		2592,2607,2613,2614,2616,2620,2624,2628,2632,2633,
		2634,2656,2661,2664,2665,2666,2667,2673,2675,2681,
		2683,2687};
*/

// Raw hits ie all hits
TH1F *hRawLe[nTdc];
TH1F *hRawTe[nTdc];
TH1F *hRawTot[nTdc];
TH1F *hGoodLe[nTdc];
TH1F *hGoodTe[nTdc];
TH1F *hGoodTot[nTdc];

TH1F *hAllRawLe;
TH1F *hAllRawTe;
TH1F *hAllRawTot;
TH1F *hAllRawPMT;
TH1F *hAllGoodLe;
TH1F *hAllGoodTe;
TH1F *hAllGoodTot;
TH1F *hAllGoodPMT;

TH2F *h2AllGoodLe;
TH2F *h2AllGoodTe;
TH2F *h2AllGoodTot;

TH2F *h2TDCTOTvsLE;
TH2F *h2CDetX1vsX2;

TH1F *hRefRawLe;
TH1F *hRefRawTe;
TH1F *hRefRawTot;
TH1F *hRefRawPMT;
TH1F *hRefGoodLe;
TH1F *hRefGoodTe;
TH1F *hRefGoodTot;
TH1F *hRefGoodPMT;

TH1F *hMultiplicityL[nTdc];
TH1F *hMultiplicity;


// hit channel id
TH1F *hHitPMT;
TH1F *hRow;
TH1F *hRowLayer1Side1;
TH1F *hRowLayer1Side2;
TH1F *hRowLayer2Side1;
TH1F *hRowLayer2Side2;
TH1F *hLayer;
TH1F *hCol;

TH1F *hnhits1;
TH1F *hngoodhits1;
TH1F *hngoodTDChits1;
TH1F *hnhits2;
TH1F *hngoodhits2;
TH1F *hngoodTDChits2;
TH1F *hnhits_ev;
TH1F *hngoodhits_ev;
TH1F *hngoodTDChits_ev;

TH1F *hngoodTDCpaddles;

TH1F *hHitX;
TH1F *hHitY;
TH1F *hHitZ;

TH2F *hHitXY1;
TH2F *hHitXY2;

TH1F *hXECal;
TH1F *hYECal;

TH2F *hXECalCDet1;
TH2F *hXECalCDet2;
TH2F *hYECalCDet1;
TH2F *hYECalCDet2;

TH2F *hXYECal;

TH1F *hXDiffECalCDet1;
TH1F *hXPlusECalCDet1;
TH1F *hXDiffECalCDet2;
TH1F *hXPlusECalCDet2;

TH2F *hXCDet1CDet2;
  
// 2D histograms
TH2F* h2d_RawLE;
TH2F* h2d_RawTE;
TH2F* h2d_RawTot;

TH2F* h2d_GoodLE;
TH2F* h2d_GoodTE;
TH2F* h2d_GoodTot;

TH2F* h2d_Mult;

using namespace std;

bool check_bad(int pmt, bool suppress_bad) {
	bool flag = false;
	if (!suppress_bad) return flag;
	for (int i=0;i<num_bad;i++) {
		if (pmt == bad_channels[i])  flag = true;
	}
	return flag;
}

void PlotRawTDC2DCut(Int_t RunNumber1=2845, Int_t nevents=60000, Int_t neventsr=300000, 
	const TString InFilePrefix="cdet",
	Double_t LeMin = 0.0, Double_t LeMax = 4000.0,
	Double_t TotMin = 2.0, Double_t TotMax = 60.0, 
	Int_t nhitcutlow1 = 1, Int_t nhitcuthigh1 = 50,
	Int_t nhitcutlow2 = 0, Int_t nhitcuthigh2 = 50,
	Double_t XDiffCut = 0.2, Double_t XOffset = 0.0,
        Int_t layer_choice=1,	
	bool suppress_bad = false, Int_t nruns=30){

	int NTDCBins = (LeMax-LeMin)/4; // 4 ns is the trigger time 
					// resolution, so this is the best we can hope for, I think.

  // InFile is the input file without absolute path and without .root suffix
  // nevents is how many events to analyse, -1 for all
  
  // To execute
  // root -l
  // .L PlotRawTDC2D.C+
  // PlotRawTDC2D("filename", -1)
  double TDCBinLow = LeMin;
  double TDCBinHigh = LeMax;
  
  // hit channel id
  
  hHitPMT = new TH1F("hHitPMT","hHitPMT",nTdc,0,nTdc);
  hnhits1 = new TH1F("hnhits1","hnhits1",100,1,101);
  hngoodhits1 = new TH1F("hngoothits1","hngoodhits1",100,1,101);
  hngoodTDChits1 = new TH1F("hngoodTDChits1","hngoodTDChits1",100,1,101);
  hnhits2 = new TH1F("hnhits2","hnhits2",100,1,101);
  hngoodhits2 = new TH1F("hngoothits2","hngoodhits2",100,1,101);
  hngoodTDChits2 = new TH1F("hngoodTDChits2","hngoodTDChits2",100,1,101);
  hnhits_ev = new TH1F("hnhits_ev","hnhits_ev",500,0,200000);
  hngoodhits_ev = new TH1F("hngoothits_ev","hngoodhits_ev",500,0,200000);
  hngoodTDChits_ev = new TH1F("hngoodTDChits_ev","hngoodTDChits_ev",500,0,200000);

  hngoodTDCpaddles = new TH1F("hngoodTDCpaddles","hngoodTDCpaddles",50,0,50);

  hRow = new TH1F("RowNumber","RowNumber",680, 0, 680);
  hRowLayer1Side1 = new TH1F("RowNumberL1S1","RowNumberL1S1",680, 0, 680);
  hRowLayer1Side2 = new TH1F("RowNumberL1S2","RowNumberL1S2",680, 0, 680);
  hRowLayer2Side1 = new TH1F("RowNumberL2S1","RowNumberL2S1",680, 0, 680);
  hRowLayer2Side2 = new TH1F("RowNumberL2S2","RowNumberL2S2",680, 0, 680);
  hLayer = new TH1F("LayerNumber","LayerNumber",3, 0, 3);
  hCol = new TH1F("ColNumber","ColNumber",3, 0, 3);

  hHitX = new TH1F("HitXposition","HitXPosition",1000,-2.0,2.0);
  hHitY = new TH1F("HitYposition","HitYPosition",200,-0.5,0.5);
  hHitZ = new TH1F("HitZposition","HitZPosition",200,7.5,8.0);
  
  hHitXY1 = new TH2F("HitXY1position","HitXY1Position",9,-1.0,1.0,800,-2.0,2.0);
  hHitXY2 = new TH2F("HitXY2position","HitXY2Position",9,-1.0,1.0,800,-2.0,2.0);
  
  hXECal = new TH1F("XEcal","XEcal",200,-1.5,1.5);
  hYECal = new TH1F("YEcal","YEcal",200,-1.0,1.0);
  
  hXECalCDet1 = new TH2F("XECalCDet1","XECalCDet1",100,-2.0,2.0,100,-2.0,2.0);
  hXECalCDet2 = new TH2F("XECalCDet2","XECalCDet2",100,-2.0,2.0,100,-2.0,2.0);
  hYECalCDet1 = new TH2F("YECalCDet1","YECalCDet1",100,-1.0,1.0,9,-1.0,1.0);
  hYECalCDet2 = new TH2F("YECalCDet2","YECalCDet2",100,-1.0,1.0,9,-1.0,1.0);
  
  hXYECal = new TH2F("XYECal","XYECal",200,-2.0,2.0,200,-2.0,2.0);
  
  hXDiffECalCDet1 = new TH1F("XDiffECalCDet1","XDiffECalCDet1",200,-3.0,3.0);
  hXPlusECalCDet1 = new TH1F("XPlusECalCDet1","XPlusECalCDet1",200,-3.0,3.0);
  hXDiffECalCDet2 = new TH1F("XDiffECalCDet2","XDiffECalCDet2",200,-3.0,3.0);
  hXPlusECalCDet2 = new TH1F("XPlusECalCDet2","XPlusECalCDet2",200,-3.0,3.0);
  
  hXCDet1CDet2 = new TH2F("XCDet1CDet2","XCDet1CDet2",200,-0.5,0.5,200,-0.5,0.5);
  
  // 2D histograms
  h2d_RawLE  = new TH2F("h2d_RawLE","", NTDCBins,TDCBinLow,TDCBinHigh,nTdc+1,0,nTdc+1);
  h2d_RawTE  = new TH2F("h2d_RawTE","", NTDCBins,TDCBinLow,TDCBinHigh,nTdc+1,0,nTdc+1);
  h2d_RawTot = new TH2F("h2d_RawTot","", NTotBins,TotBinLow,TotBinHigh,nTdc+1,0,nTdc+1);
  h2d_Mult   = new TH2F("h2d_Mult","", 100,0,100,nTdc+1,0,nTdc+1);
  
  hMultiplicity = new TH1F("hMultiplicity","hMultiplicity",20,0,20);
  
  for(Int_t tdc=0; tdc<nTdc; tdc++){
    hMultiplicityL[tdc] =  new TH1F(TString::Format("hMultiplicity_Bar%d",tdc),
		      TString::Format("hMultiplicity_Bar%d",tdc),
		      10, 0, 10);
  }// element loop

  hAllRawLe = new TH1F(TString::Format("hRawLe"),
            TString::Format("hRawLe"),
            NTDCBins, TDCBinLow, TDCBinHigh);
  hAllRawTe = new TH1F(TString::Format("hRawTe"),
            TString::Format("hRawTe"),
            NTDCBins, TDCBinLow, TDCBinHigh);
  hAllRawTot = new TH1F(TString::Format("hRawTot"),
            TString::Format("hRawTot"),
            NTotBins, TotBinLow, TotBinHigh);
  hAllRawPMT = new TH1F(TString::Format("hRawPMT"),
            TString::Format("hRawPMT"),
            nTdc, 0, nTdc);
  hAllGoodLe = new TH1F(TString::Format("hAllGoodLe"),
            TString::Format("hAllGoodLe"),
            NTDCBins, TDCBinLow, TDCBinHigh);
  hAllGoodTe = new TH1F(TString::Format("hAllGoodTe"),
            TString::Format("hAllGoodTe"),
            NTDCBins, TDCBinLow, TDCBinHigh);
  hAllGoodTot = new TH1F(TString::Format("hAllGoodTot"),
            TString::Format("hAllGoodTot"),
            NTotBins, TotBinLow, TotBinHigh);
  hAllGoodPMT = new TH1F(TString::Format("hAllGoodPMT"),
            TString::Format("hAllGoodPMT"),
            nTdc, 0, nTdc);
  h2AllGoodLe = new TH2F(TString::Format("h2AllGoodLe"),
            TString::Format("h2AllGoodLe"),nTdc,0,nTdc,
            NTDCBins, TDCBinLow, TDCBinHigh);
  h2AllGoodTe = new TH2F(TString::Format("h2AllGoodTe"),
            TString::Format("h2AllGoodTe"),nTdc,0,nTdc,
            NTDCBins, TDCBinLow, TDCBinHigh);
  h2AllGoodTot = new TH2F(TString::Format("h2AllGoodTot"),
            TString::Format("h2AllGoodTot"),nTdc,0,nTdc,
            NTotBins, TotBinLow, TotBinHigh);

  h2TDCTOTvsLE = new TH2F(TString::Format("h2TDCTOTvsLE"),
            TString::Format("h2TDCTOTvsLE"),NTotBins,TotBinLow,TotBinHigh,
            NTDCBins, TDCBinLow, TDCBinHigh);
  h2CDetX1vsX2 = new TH2F(TString::Format("h2CDetX1vsX2"),
            TString::Format("h2CDetX1vsX2"),1000, -2.0, 2.0, 
            1000, -2.0, 2.0);

  hRefRawLe = new TH1F(TString::Format("hRefRawLe"),
            TString::Format("hRefRawLe"),
            NTDCBins, TDCBinLow, TDCBinHigh);
  hRefRawTe = new TH1F(TString::Format("hRefRawTe"),
            TString::Format("hRefRawTe"),
            NTDCBins, TDCBinLow, TDCBinHigh);
  hRefRawTot = new TH1F(TString::Format("hRefRawTot"),
            TString::Format("hRefRawTot"),
            NTotBins, TotBinLow, TotBinHigh);
  hRefRawPMT = new TH1F(TString::Format("hRefRawPMT"),
            TString::Format("hRefRawPMT"),
            2720, 0, 2720);
  hRefGoodLe = new TH1F(TString::Format("hRefGoodLe"),
            TString::Format("hRefGoodLe"),
            NTDCBins, TDCBinLow, TDCBinHigh);
  hRefGoodTe = new TH1F(TString::Format("hRefGoodTe"),
            TString::Format("hRefGoodTe"),
            NTDCBins, TDCBinLow, TDCBinHigh);
  hRefGoodTot = new TH1F(TString::Format("hRefGoodTot"),
            TString::Format("hRefGoodTot"),
            NTotBins, TotBinLow, TotBinHigh);
  hRefGoodPMT = new TH1F(TString::Format("hRefGoodPMT"),
            TString::Format("hRefGoodPMT"),
            2720, 0, 2720);
  
  
  for(Int_t bar=0; bar<(nTdc); bar++){
    // raw hits
    // leading edge
    hRawLe[bar] = new TH1F(TString::Format("hRawLe_Bar%d",bar),
 	    TString::Format("hRawLe_Bar%d",bar),
	    NTDCBins, TDCBinLow, TDCBinHigh);
    // trailing edge 
    hRawTe[bar] = new TH1F(TString::Format("hRawTe_Bar%d",bar),
	    TString::Format("hRawTe_Bar%d",bar),
	    NTDCBins, TDCBinLow, TDCBinHigh);
    // tot 
    hRawTot[bar] = new TH1F(TString::Format("hRawTot_Bar%d",bar),
			     TString::Format("hRawTot_Bar%d",bar),
			     NTotBins, TotBinLow, TotBinHigh);
    // good hits
    // leading edge
    hGoodLe[bar] = new TH1F(TString::Format("hGoodLe_Bar%d",bar),
 	    TString::Format("hGoodLe_Bar%d",bar),
	    NTDCBins, TDCBinLow, TDCBinHigh);
    // trailing edge 
    hGoodTe[bar] = new TH1F(TString::Format("hGoodTe_Bar%d",bar),
	    TString::Format("hGoodTe_Bar%d",bar),
	    NTDCBins, TDCBinLow, TDCBinHigh);
    // tot 
    hGoodTot[bar] = new TH1F(TString::Format("hGoodTot_Bar%d",bar),
			     TString::Format("hGoodTot_Bar%d",bar),
			     NTotBins, TotBinLow, TotBinHigh);
  }// bar loop


  //========================================================= Get data from tree
  if(!T) { 
    // TString sInFile = REPLAYED_DIR + "/" + InFile + ".root";
    T = new TChain("T");

    TString subfile, sInFile;

    subfile = TString::Format("_%d_%d",RunNumber1,neventsr);
    sInFile = REPLAYED_DIR + "/" + InFilePrefix + subfile + ".root";
    cout << "Input ROOT file = " << sInFile << endl;
    cout << "Adding " << sInFile << endl;
    T->Add(sInFile);
    cout << "Adding " << nruns << " files ... " << endl;
    for (Int_t i=1; i<=nruns; i++) {
    	subfile = TString::Format("_%d_%d_%d",RunNumber1,neventsr,i);
    	//subfile = TString::Format("_%d_1000000_%d",RunNumber1,i);
    	sInFile = REPLAYED_DIR + "/" + InFilePrefix + subfile + ".root";
    	cout << "Input ROOT file = " << sInFile << endl;
    	cout << "Adding " << sInFile << endl;
    	T->Add(sInFile);
    }
    
    
    // disable all branches
    T->SetBranchStatus("*",0);
    // enable branches
    T->SetBranchStatus("earm.cdet.*",1);
    T->SetBranchStatus("earm.ecal.*",1);

    T->SetBranchAddress("earm.cdet.tdc_mult",TCDet::TDCmult);
    
    T->SetBranchAddress("earm.cdet.hits.TDCelemID",TCDet::RawElID);
    T->SetBranchAddress("earm.cdet.hits.t",TCDet::RawElLE);
    T->SetBranchAddress("earm.cdet.hits.t_te",TCDet::RawElTE);
    T->SetBranchAddress("earm.cdet.hits.t_tot",TCDet::RawElTot);

    T->SetBranchAddress("earm.cdet.hit.pmtnum",TCDet::GoodElID);
    T->SetBranchAddress("earm.cdet.hit.tdc_le",TCDet::GoodElLE);
    T->SetBranchAddress("earm.cdet.hit.tdc_te",TCDet::GoodElTE);
    T->SetBranchAddress("earm.cdet.hit.tdc_tot",TCDet::GoodElTot);
    
    T->SetBranchAddress("earm.cdet.hit.xhit",TCDet::GoodX);
    T->SetBranchAddress("earm.cdet.hit.yhit",TCDet::GoodY);
    T->SetBranchAddress("earm.cdet.hit.zhit",TCDet::GoodZ);
    
    T->SetBranchAddress("earm.cdet.hit.row",TCDet::GoodCol);
    T->SetBranchAddress("earm.cdet.hit.col",TCDet::GoodRow);
    T->SetBranchAddress("earm.cdet.hit.layer",TCDet::GoodLayer);

        
    T->SetBranchAddress("earm.cdet.nhits",&TCDet::nhits);
    T->SetBranchAddress("earm.cdet.ngoodhits",&TCDet::ngoodhits);
    T->SetBranchAddress("earm.cdet.ngoodTDChits",&TCDet::ngoodTDChits);
    T->SetBranchAddress("earm.ecal.x",&TCDet::GoodECalX);
    T->SetBranchAddress("earm.ecal.y",&TCDet::GoodECalY);

    // enable vector size branches
    T->SetBranchAddress("Ndata.earm.cdet.tdc_mult",&TCDet::NdataMult); 
    T->SetBranchAddress("Ndata.earm.cdet.hits.TDCelemID",&TCDet::NdataRawElID); 
    T->SetBranchAddress("Ndata.earm.cdet.hits.t",&TCDet::NdataRawElLE); 
    T->SetBranchAddress("Ndata.earm.cdet.hits.t_te",&TCDet::NdataRawElTE); 
    T->SetBranchAddress("Ndata.earm.cdet.hits.t_tot",&TCDet::NdataRawElTot); 
    
    T->SetBranchAddress("Ndata.earm.cdet.hit.pmtnum",&TCDet::NdataGoodElID);
    T->SetBranchAddress("Ndata.earm.cdet.hit.tdc_le",&TCDet::NdataGoodElLE);
    T->SetBranchAddress("Ndata.earm.cdet.hit.tdc_te",&TCDet::NdataGoodElTE);
    T->SetBranchAddress("Ndata.earm.cdet.hit.tdc_tot",&TCDet::NdataGoodElTot);
    
    T->SetBranchAddress("Ndata.earm.cdet.hit.xhit",&TCDet::NdataGoodX);
    T->SetBranchAddress("Ndata.earm.cdet.hit.yhit",&TCDet::NdataGoodY);
    T->SetBranchAddress("Ndata.earm.cdet.hit.zhit",&TCDet::NdataGoodZ);
    
    T->SetBranchAddress("Ndata.earm.cdet.hit.row",&TCDet::NdataGoodCol);
    T->SetBranchAddress("Ndata.earm.cdet.hit.col",&TCDet::NdataGoodRow);
    T->SetBranchAddress("Ndata.earm.cdet.hit.layer",&TCDet::NdataGoodLayer);
    

  }//setting tree
  
  //========================================================= Check no of events
  Int_t Nev = T->GetEntries();
  cout << "N entries in tree is " << Nev << endl;
  Int_t NEventsAnalysis;
  if(nevents==-1) NEventsAnalysis = Nev;
  else NEventsAnalysis = nevents;
  cout << "Running analysis for " << NEventsAnalysis << " events" << endl;
  

  
  //==================================================== Create output root file
  // root file for viewing fits
  
  TString outrootfile = ANALYSED_DIR + "/RawTDC_" + InFilePrefix + ".root";
  TFile *f = new TFile(outrootfile, "RECREATE");



  //================================================================= Event Loop
  // variables outside event loop
  Int_t EventCounter = 0;
  cout << "Starting Event Loop" << endl;

  // event loop start
  for(Int_t event=0; event<NEventsAnalysis; event++){
    
    T->GetEntry(event);
    EventCounter++;
    if (EventCounter % 100 == 0) {
	cout << EventCounter << "/" << NEventsAnalysis << "/ Nhits = " << (Int_t)TCDet::nhits << endl;
    	for (Int_t nfill=0; nfill<(Int_t)TCDet::nhits; nfill++) {hnhits_ev->Fill(EventCounter);}
    	for (Int_t nfill=0; nfill<(Int_t)TCDet::ngoodhits; nfill++) {hngoodhits_ev->Fill(EventCounter);}
    	for (Int_t nfill=0; nfill<(Int_t)TCDet::ngoodTDChits; nfill++) {hngoodTDChits_ev->Fill(EventCounter);}
    }


    //cout << "Raw TDC hit loop: " << TCDet::NdataRawElID << endl;
    
    
    for(Int_t el=0; el<TCDet::NdataRawElID; el++){
	//if ((Int_t)TCDet::RawElID[el] > 0) cout << "el = " << el << " Hit ID = " << (Int_t)TCDet::RawElID[el] << "    TDC = " << TCDet::RawElLE[el]*TDC_calib_to_ns << endl;
	//cout << "Raw ID = " << TCDet::RawElID[el] << " raw le = " << TCDet::RawElLE[el] << " raw te = " << TCDet::RawElTE[el] << " raw tot = " << TCDet::RawElTot[el] << endl;
	if (TCDet::RawElLE[el] >= LeMin/TDC_calib_to_ns && TCDet::RawElLE[el] <= LeMax/TDC_calib_to_ns &&
		TCDet::RawElTot[el] >= TotMin/TDC_calib_to_ns && TCDet::RawElTot[el] <= TotMax/TDC_calib_to_ns &&
		TCDet::TDCmult[el] < TDCmult_cut ) {
	  //cout << "el = " << el << " Raw ID = " << TCDet::RawElID[el] << " raw le = " << 
	//	TCDet::RawElLE[el] << " raw te = " << TCDet::RawElTE[el] << " raw tot = " << 
	//	TCDet::RawElTot[el] << " CDet X = " << TCDet::GoodX[el] << " ECal X = " << TCDet::GoodECalX << endl;

	  if ( !check_bad(TCDet::RawElID[el],suppress_bad) ) {
	   //cout << " el = " << el << endl;
	   //cout << " tdc = " << TCDet::RawElLE[el]*TDC_calib_to_ns << endl;
	   if ( TCDet::RawElID[el] < NumCDetPaddles ) {
	    //if ((Int_t)TCDet::RawElID[el] > nTdc) cout << " CDet ID = " << (Int_t)TCDet::RawElID[el] << "    TDC = " << TCDet::RawElLE[el]*TDC_calib_to_ns << endl;
	    hRawLe[(Int_t)TCDet::RawElID[el]]->Fill(TCDet::RawElLE[el]*TDC_calib_to_ns);
	    hRawTe[(Int_t)TCDet::RawElID[el]]->Fill(TCDet::RawElTE[el]*TDC_calib_to_ns);
	    hRawTot[(Int_t)TCDet::RawElID[el]]->Fill(TCDet::RawElTot[el]*TDC_calib_to_ns);
 	    hAllRawLe->Fill(TCDet::RawElLE[el]*TDC_calib_to_ns);
	    hAllRawTe->Fill(TCDet::RawElTE[el]*TDC_calib_to_ns);
	    hAllRawTot->Fill(TCDet::RawElTot[el]*TDC_calib_to_ns);
	    hAllRawPMT->Fill(TCDet::RawElID[el]);

	    h2d_RawLE->Fill(TCDet::RawElLE[el]*TDC_calib_to_ns, (Int_t)TCDet::RawElID[el]);
	    h2d_RawTE->Fill(TCDet::RawElTE[el]*TDC_calib_to_ns, (Int_t)TCDet::RawElID[el]);
	    h2d_RawTot->Fill(TCDet::RawElTot[el]*TDC_calib_to_ns, (Int_t)TCDet::RawElID[el]);

	   } else {
	    if ((Int_t)TCDet::RawElID[el] >= NumCDetPaddles) cout << " Ref  ID = " << (Int_t)TCDet::RawElID[el] << "    TDC = " << TCDet::RawElLE[el]*TDC_calib_to_ns << endl;
 	    hRefRawLe->Fill(TCDet::RawElLE[el]*TDC_calib_to_ns);
	    hRefRawTe->Fill(TCDet::RawElTE[el]*TDC_calib_to_ns);
	    hRefRawTot->Fill(TCDet::RawElTot[el]*TDC_calib_to_ns);
	    hRefRawPMT->Fill(TCDet::RawElID[el]);
	   }
	  }
	}

    }// all raw tdc hit loop
    
    int nhitsc1 = 0;
    int nhitsc2 = 0;
    int ngoodhitsc1 = 0;
    int ngoodhitsc2 = 0;
    int ngoodTDChitsc1 = 0;
    int ngoodTDChitsc2 = 0;
    for (int j=0; j<nTdc; j++) {
	TCDet::ngoodTDChits_paddles[j]=0;
    }
    TCDet::ngoodTDCpaddles=0;
    
    for(Int_t el=0; el<TCDet::NdataGoodElID; el++){
	if (TCDet::GoodElLE[el] >= LeMin/TDC_calib_to_ns && TCDet::GoodElLE[el] <= LeMax/TDC_calib_to_ns 
		&& TCDet::GoodElTot[el] >= TotMin/TDC_calib_to_ns && TCDet::GoodElTot[el] <= TotMax/TDC_calib_to_ns
		&& TCDet::GoodX[el] < xcut && TCDet::TDCmult[el] < TDCmult_cut
		&& (TCDet::GoodX[el]-(TCDet::GoodECalX)-XOffset) <= XDiffCut && (TCDet::GoodX[el]-(TCDet::GoodECalX)-XOffset) >= -1.0*XDiffCut
		&& TCDet::GoodECalX > -1.5 && TCDet::GoodECalX < 1.5  
		&& TCDet::GoodECalY > -0.8 && TCDet::GoodECalY < 0.8  
	) {

	  if ( !check_bad(TCDet::GoodElID[el], suppress_bad) ) {
	   if ( TCDet::GoodElID[el] < NumCDetPaddles )  {

	    int sbsrown = (Int_t)TCDet::GoodRow[el];
	    int sbscoln = (Int_t)TCDet::GoodCol[el];
	    int mylayern = sbscoln/2;
	    int mypaddlen = sbscoln*672 + sbsrown;
	    //cout << "Hit number " << el << ":    Paddle = " << mypaddlen << " Row = " << sbsrown  << " Col = " << sbscoln  << " hits = " << TCDet::ngoodTDChits_paddles[mypaddlen] << endl;
	    //cout << "el = " << el << " Good ID = " << TCDet::GoodElID[el] << " Good le = " << 
	//	TCDet::GoodElLE[el] << " Good te = " << TCDet::GoodElTE[el] << " Good tot = " << 
	//	TCDet::GoodElTot[el] << " CDet X = " << TCDet::GoodX[el] << " ECal X = " << TCDet::GoodECalX << endl;
	    if (mylayern == 0) {
		nhitsc1++;
	        ngoodhitsc1++;
		ngoodTDChitsc1++;
	    } else {
		nhitsc2++;
	        ngoodhitsc2++;
		ngoodTDChitsc2++;
	    }
	    TCDet::ngoodTDChits_paddles[mypaddlen]++;
	   }
	  }
	}
    }
    for (int j=0; j<nTdc; j++) {
	if (TCDet::ngoodTDChits_paddles[j] > 0) {
		 TCDet::ngoodTDCpaddles++;
		 //cout << "Paddle = " << j <<  "  nhits = " << TCDet::ngoodTDChits_paddles[j] << endl;
	}
    }
    //cout << "event " << event << endl;
    //cout << "Number of good layer 1 hits: " << ngoodTDChitsc1 << endl;
    //cout << "Number of good layer 2 hits: " << ngoodTDChitsc2 << endl;
    //cout << "Layer 1 Hit Cut " << nhitcutlow1 << " " << nhitcuthigh1 << endl;
    //cout << "Layer 2 Hit Cut " << nhitcutlow2 << " " << nhitcuthigh2 << endl;

    hngoodTDCpaddles->Fill(TCDet::ngoodTDCpaddles);

    
    for(Int_t el=0; el<TCDet::NdataGoodElID; el++){
	if (TCDet::GoodElLE[el] >= LeMin/TDC_calib_to_ns && TCDet::GoodElLE[el] <= LeMax/TDC_calib_to_ns 
		&& TCDet::GoodElTot[el] >= TotMin/TDC_calib_to_ns && TCDet::GoodElTot[el] <= TotMax/TDC_calib_to_ns 
		&& TCDet::GoodX[el] < xcut && TCDet::TDCmult[el] < TDCmult_cut 
		&& ngoodTDChitsc1 >= nhitcutlow1  && ngoodTDChitsc2 >= nhitcutlow2 
		&& ngoodTDChitsc1 <= nhitcuthigh1 && ngoodTDChitsc2 <= nhitcuthigh2 
		&& (TCDet::GoodX[el]-(TCDet::GoodECalX)-XOffset) <= XDiffCut && (TCDet::GoodX[el]-(TCDet::GoodECalX)-XOffset) >= -1.0*XDiffCut
		&& TCDet::GoodECalX > -1.5 && TCDet::GoodECalX < 1.5  
		&& TCDet::GoodECalY > -0.8 && TCDet::GoodECalY < 0.8  
	) {

	  if ( !check_bad(TCDet::GoodElID[el], suppress_bad) ) {
	   if ( TCDet::GoodElID[el] < NumCDetPaddles )  {
    	    //cout << "event " << event << endl;
	    //cout << "el = " << el << " Good ID = " << TCDet::GoodElID[el] << " Good le = " << 
		//TCDet::GoodElLE[el] << " Good te = " << TCDet::GoodElTE[el] << " Good tot = " << 
		//TCDet::GoodElTot[el] << " CDet X = " << TCDet::GoodX[el] << " ECal X = " << TCDet::GoodECalX << endl;

	    //cout << "Filling good timing histos ... " << ngoodTDChitsc1 << " " << endl;
	    
	    //std::cout << "Layer = " << (Int_t)TCDet::GoodLayer[el] << " Side = " << (Int_t)TCDet::GoodCol[el] << std::endl;
	    
	    int sbscol = (Int_t)TCDet::GoodCol[el];
	    int myside = sbscol%2;
	    int mylayer = sbscol/2;
	    
	    if ((layer_choice == 1 && mylayer == 0) || (layer_choice == 2 && mylayer == 1) || layer_choice == 3) 
		{
	    	hGoodLe[(Int_t)TCDet::GoodElID[el]]->Fill(TCDet::GoodElLE[el]*TDC_calib_to_ns);
	    	hGoodTe[(Int_t)TCDet::GoodElID[el]]->Fill(TCDet::GoodElTE[el]*TDC_calib_to_ns);
	    	hGoodTot[(Int_t)TCDet::GoodElID[el]]->Fill(TCDet::GoodElTot[el]*TDC_calib_to_ns);
	    	hAllGoodLe->Fill(TCDet::GoodElLE[el]*TDC_calib_to_ns);
	    	hAllGoodTe->Fill(TCDet::GoodElTE[el]*TDC_calib_to_ns);
	    	hAllGoodTot->Fill(TCDet::GoodElTot[el]*TDC_calib_to_ns);
	    	hAllGoodPMT->Fill(TCDet::GoodElID[el]);
	    	h2AllGoodLe->Fill(TCDet::GoodElID[el],TCDet::GoodElLE[el]*TDC_calib_to_ns);
	    	h2AllGoodTe->Fill(TCDet::GoodElID[el],TCDet::GoodElTE[el]*TDC_calib_to_ns);
	    	h2AllGoodTot->Fill(TCDet::GoodElID[el],TCDet::GoodElTot[el]*TDC_calib_to_ns);

	    	h2TDCTOTvsLE->Fill(TCDet::GoodElTot[el]*TDC_calib_to_ns,TCDet::GoodElLE[el]*TDC_calib_to_ns);
	  
	    	hHitPMT->Fill((Int_t)TCDet::GoodElID[el]);
	    	hRow->Fill((Int_t)TCDet::GoodRow[el]);
	    }
	
	    if (myside == 0) {
		if (mylayer == 0) {
			hRowLayer1Side1->Fill((Int_t)TCDet::GoodRow[el]);
		} else {
			hRowLayer2Side1->Fill((Int_t)TCDet::GoodRow[el]);
		}
	    } else {
		if(mylayer == 0) {
			hRowLayer1Side2->Fill((Int_t)TCDet::GoodRow[el]);
		} else {
			hRowLayer2Side2->Fill((Int_t)TCDet::GoodRow[el]);
		}
	    }	
			
	    hCol->Fill(myside);
	    hLayer->Fill(mylayer);

	    hHitX->Fill(TCDet::GoodX[el]);
	    hHitY->Fill(TCDet::GoodY[el]);
	    hHitZ->Fill(TCDet::GoodZ[el]);
	    if (mylayer==0) {
		hHitXY1->Fill(TCDet::GoodY[el],TCDet::GoodX[el]);
		hXECalCDet1->Fill(TCDet::GoodX[el],TCDet::GoodECalX);
		hYECalCDet1->Fill(TCDet::GoodY[el],TCDet::GoodECalY);
	    	hXDiffECalCDet1->Fill(TCDet::GoodX[el]-TCDet::GoodECalX);
	    	hXPlusECalCDet1->Fill(TCDet::GoodX[el]+TCDet::GoodECalX);
	    } else {
		hHitXY2->Fill(TCDet::GoodY[el],TCDet::GoodX[el]);
		hXECalCDet2->Fill(TCDet::GoodX[el],TCDet::GoodECalX);
		hYECalCDet2->Fill(TCDet::GoodY[el],TCDet::GoodECalY);
	    	hXDiffECalCDet2->Fill(TCDet::GoodX[el]-TCDet::GoodECalX);
	    	hXPlusECalCDet2->Fill(TCDet::GoodX[el]+TCDet::GoodECalX);
	    }

	    hXYECal->Fill(TCDet::GoodECalY,TCDet::GoodECalX);
	    hXECal->Fill(TCDet::GoodECalX);
	    hYECal->Fill(TCDet::GoodECalY);

	   } else {
	    hRefGoodLe->Fill(TCDet::GoodElLE[el]*TDC_calib_to_ns);
	    hRefGoodTe->Fill(TCDet::GoodElTE[el]*TDC_calib_to_ns);
	    hRefGoodTot->Fill(TCDet::GoodElTot[el]*TDC_calib_to_ns);
	    hRefGoodPMT->Fill(TCDet::GoodElID[el]*TDC_calib_to_ns);
	   }
	  }
	}


    }// all good tdc hit loop
    
    hnhits1->Fill(nhitsc1);
    hngoodhits1->Fill(ngoodhitsc1);
    hngoodTDChits1->Fill(ngoodTDChitsc1);
    hnhits2->Fill(nhitsc2);
    hngoodhits2->Fill(ngoodhitsc2);
    hngoodTDChits2->Fill(ngoodTDChitsc2);

    //cout << "Element loop: " << TCDet::NdataMult << endl;
    for(Int_t tdc=0; tdc<TCDet::NdataMult; tdc++){
      if (!check_bad(TCDet::RawElID[tdc],suppress_bad)) {
      hMultiplicity->Fill(TCDet::TDCmult[tdc]);
	hMultiplicityL[(Int_t)TCDet::RawElID[tdc]]->Fill(TCDet::TDCmult[tdc]);
	if( TCDet::TDCmult[tdc] != 0 )
	  h2d_Mult->Fill(TCDet::TDCmult[tdc], (Int_t)TCDet::RawElID[tdc] );
      }
    }// element loop


  }// event loop

  for (Int_t b=0; b<NumCDetPaddles; b++) {
	if (hRawLe[b]->GetEntries() > EventCounter/HotChannelRatio) {
		int myhotlayer = b/1344 + 1;
		int myhotside = (b%1344)/672 + 1;
		int myhotmodule = (b%672)/244 + 1;
		int myhotbar = (b%672)%244/16 + 1;
		int myhotpaddle = ((b%672)%244)%16 + 1;
		std::cout << "Hot PMT!! ID = " << b << "  layer = " << myhotlayer <<
		"   side = " << myhotside << "   module = " << myhotmodule <<
		"   bar = " << myhotbar << "   paddle_PMT = " << myhotpaddle << std::endl;
	}
  }
  
    
    //========================================================== Write histos
  for(Int_t b=0; b<nTdc; b++){
    // hRawLe[b]->GetXaxis()->SetLabelSize(0.06);
    hRawLe[b]->GetXaxis()->SetTitle("time (ns)");
    // hRawLe[b]->GetXaxis()->SetTitleSize(0.05);
    hRawLe[b]->Write();
    // hRawLe[b]->GetXaxis()->SetLabelSize(0.06);
    hRawTe[b]->GetXaxis()->SetTitle("time (ns)");
    // hRawTe[b]->GetXaxis()->SetTitleSize(0.05);
    hRawTe[b]->Write();
    // hRawTe[b]->GetXaxis()->SetLabelSize(0.06);
    hRawTot[b]->GetXaxis()->SetTitle("tot (ns)");
    // hRawTot[b]->GetXaxis()->SetTitleSize(0.05);
    hRawTot[b]->Write();
    // hMultiplicityL[b]->GetXaxis()->SetLabelSize(0.06);
    hMultiplicityL[b]->GetXaxis()->SetTitle("tdc ref hit mult");
    hMultiplicityL[b]->Write();
  }
  // hMultiplicity->GetXaxis()->SetLabelSize(0.06);
  hMultiplicity->GetXaxis()->SetTitle("tdc hit multiplicity");
  hMultiplicity->Write();
  // hHitPMT->GetXaxis()->SetLabelSize(0.06);
  hHitPMT->GetXaxis()->SetTitle("Bar ID of Left PMT Hit");
  hHitPMT->SetTitle("");
  hHitPMT->Write();

  // 2D histograms
  h2d_RawLE->GetXaxis()->SetTitle("TDC Leading Edge Time [ns]");
  h2d_RawLE->GetYaxis()->SetTitle("PMT number (Left)");
  h2d_RawLE->SetTitle("");
  h2d_RawLE->Write();
  h2d_RawTE->GetXaxis()->SetTitle("TDC Trailing Edge Time [ns]");
  h2d_RawTE->GetYaxis()->SetTitle("PMT number (Left)");
  h2d_RawTE->SetTitle("");
  h2d_RawTE->Write();

  h2d_RawTot->GetXaxis()->SetTitle("TDC Time-over-threshold [ns]");
  h2d_RawTot->GetYaxis()->SetTitle("PMT number (Left)");
  h2d_RawTot->SetTitle("");
  h2d_RawTot->Write();

  h2d_Mult->GetXaxis()->SetTitle("TDC Multiplicity [ns]");
  h2d_Mult->GetYaxis()->SetTitle("PMT number (Left)");
  h2d_Mult->SetTitle("");
  h2d_Mult->Write();

  //========================================================== Close output file
  f->Close();



  //================================================================== End Macro
}// end main

TCanvas *plotAllTDC(){

  TCanvas *caa = new TCanvas("all", "all", 50,50,1000,1000);
  caa->Divide(2,4,0.01,0.01,0);

  caa->cd(1);
  hAllRawLe->SetFillColor(kRed);
  hAllRawLe->SetMinimum(0.0);
  hAllRawLe->Draw();
  
  caa->cd(2);
  hAllGoodLe->SetFillColor(kBlue);
  hAllGoodLe->SetMinimum(0.0);
  hAllGoodLe->Draw();
  
  caa->cd(3);
  hAllRawTe->SetFillColor(kRed);
  hAllRawTe->SetMinimum(0.0);
  hAllRawTe->Draw();
  
  caa->cd(4);
  hAllGoodTe->SetFillColor(kBlue);
  hAllGoodTe->SetMinimum(0.0);
  hAllGoodTe->Draw();
  

  caa->cd(5);
  //gPad->SetLogy();
  hAllRawTot->SetFillColor(kRed);
  hAllRawTot->Draw();

  caa->cd(6);
  //gPad->SetLogy();
  hAllGoodTot->SetFillColor(kBlue);
  hAllGoodTot->Draw();

  caa->cd(7);
  //gPad->SetLogy();
  //hs4->SetMinimum(0.);
  hAllRawPMT->SetFillColor(kRed);
  hAllRawPMT->Draw();
  
  caa->cd(8);
  //gPad->SetLogy();
  //hs4->SetMinimum(0.);
  hAllGoodPMT->SetFillColor(kBlue);
  hAllGoodPMT->Draw();

  return caa;
}

TCanvas *plotGoodTDC2D(){

  TCanvas *cac = new TCanvas("all2d", "all2d", 50,50,800,800);
  cac->Divide(2,2,0.01,0.01,0);
  
  cac->cd(1);
  gPad->SetLogz();
  h2AllGoodLe->Draw("colz");
  cac->cd(2);
  gPad->SetLogz();
  h2AllGoodTe->Draw("colz");
  cac->cd(3);
  gPad->SetLogz();
  h2AllGoodTot->Draw("colz");

  return cac;
}

TCanvas *plotRefTDC() {

  TCanvas *cbb = new TCanvas("ref", "ref", 850,50, 1200,800);
  cbb->Divide(2,2,0.01,0.01,0);

  cbb->cd(1);
  hRefRawLe->Draw();
  cbb->cd(2);
  hRefRawTe->Draw();
  cbb->cd(3);
  gPad->SetLogy();
  hRefRawTot->Draw();
  cbb->cd(4);
  gPad->SetLogy();
  hRefRawPMT->Draw();

  return cbb;
}


TCanvas *plotCDetTDC(){

  TCanvas *canvas[NumHalfModules];
  for (int cmodule=1;cmodule<=NumModules;cmodule++) {
   for (int layer=1;layer<=NumLayers;layer++) {
    for (int side=1;side<=NumSides;side++){
     int xposition1 = 100 + (layer-1)*1000 + (side-1)*500;
     int yposition1 = 100 + (NumModules-cmodule)*400;
     int xposition2 = 150 + (layer-1)*1000 + (side-1)*500;
     int yposition2 = 150 + (NumModules-cmodule)*400;

     int side_group_plot = 6*(layer-1) + 2*(cmodule-1) + (side-1);
     int elemID_start = (layer-1)*NumModules*NumBars*NumPaddles*NumLayers + (side-1)*NumModules*NumBars*NumPaddles + (cmodule-1)*NumBars*NumPaddles;

     TString cname;
     cname.Form("c1_%d_%d_%d",cmodule,layer,side);
     TCanvas *c1 = new TCanvas(cname, cname, xposition1,yposition1,380,280);
     c1->Divide(NumPaddles,NumBars, 0.01, 0.01, 0);
     canvas[side_group_plot] = c1;

     for (int ii = 0; ii < NumPaddles*NumBars; ii++) {

        c1->cd(ii+1);
        if (ii == 1) {
		cout << "side_group_plot = " << side_group_plot << endl;
	}
	
	hRawTe[elemID_start + ii + 1]->Draw();
	
     }

     cname.Form("c2_%d_%d_%d",cmodule,layer,side);
     TCanvas *c2 = new TCanvas(cname, cname, xposition2,yposition2,380,280);
     c2->Divide(NumPaddles,NumBars, 0.01, 0.01, 0);
     canvas[side_group_plot] = c2;

     for (int ii = 0; ii < NumPaddles*NumBars; ii++) {

        c2->cd(ii+1);
        if (ii == 1) {
                cout << "side_group_plot = " << side_group_plot << endl;
        }

        hRawLe[elemID_start + ii + 1]->Draw();

     }


    }
   }
  }

  return canvas[0];

}

TCanvas *plotHalfModule(int cmodule = 1, int side = 1, int layer = 1){


  TCanvas *c4 = new TCanvas("c4", "c4", 50,50,1000,1200);
  c4->Divide(NumPaddles,NumBars, 0.01, 0.01, 0);
  TCanvas *c4b = new TCanvas("c4b", "c4b", 1050,50,1000,1200);
  c4b->Divide(NumPaddles,NumBars, 0.01, 0.01, 0);

  int side_group_plot = 6*(layer-1) + 2*(cmodule-1) + (side-1);
  int elemID_start = (layer-1)*NumModules*NumBars*NumPaddles*NumLayers + (side-1)*NumModules*NumBars*NumPaddles + (cmodule-1)*NumBars*NumPaddles;


  for (int ii = 0; ii < NumPaddles*NumBars; ii++) {

        c4->cd(ii+1);
  	hRawTe[elemID_start + ii + 1]->Draw();
	c4b->cd(ii+1);
  	hRawLe[elemID_start + ii + 1]->Draw();
  }

  return c4;

}

TCanvas *plotBarTDC(int bar = 39, int side = 1, int layer = 1){

  int mymodule = (bar-1)/NumBars+1; //mymodule in this case represnts top, middle, or bottom as opposed todetector labels whe
  int paddle_start = (bar - (mymodule-1)*NumBars - 1)*NumPaddles;
  std::cout << "mymodule = " << mymodule << "  paddle_start = " << paddle_start << std::endl;

  int side_group = 6*(layer-1) + 2*(mymodule-1) + (side-1);
  int elemID_start = (layer-1)*NumModules*NumBars*NumPaddles*NumLayers + (side-1)*NumModules*NumBars*NumPaddles + (mymodule-1)*NumBars*NumPaddles + paddle_start;


  std::cout << "side_group = " << side_group << std::endl;

  TCanvas *c3 = new TCanvas("c3", "c3", 150,150,600,450);
  c3->Divide(4,4, 0.01, 0.01, 0);
  TCanvas *c333 = new TCanvas("c333", "c333", 750,150,600,450);
  c333->Divide(4,4, 0.01, 0.01, 0);
  TCanvas *c3a = new TCanvas("c3a", "c3a", 150,650,600,450);
  c3a->Divide(4,4, 0.01, 0.01, 0);
  TCanvas *c333a = new TCanvas("c333a", "c333a", 750,650,600,450);
  c333a->Divide(4,4, 0.01, 0.01, 0);

  for (int ii = 1; ii <= NumPaddles; ii++) {

        c3->cd(ii);
  	hRawLe[elemID_start + ii ]->Draw();
  }
  for (int ii = 1; ii <= NumPaddles; ii++) {

        c333->cd(ii);
  	hGoodLe[elemID_start + ii ]->Draw();
  }
  for (int ii = 1; ii <= NumPaddles; ii++) {

        c3a->cd(ii);
  	hRawTe[elemID_start + ii ]->Draw();
  }
  for (int ii = 1; ii <= NumPaddles; ii++) {

        c333a->cd(ii);
  	hGoodTe[elemID_start + ii ]->Draw();
  }

  return c3;

}

TCanvas *plotTOTvsLE(){

  TCanvas *c123 = new TCanvas("c123", "c123", 50,50,1000,1000);

  c123->cd();
  h2TDCTOTvsLE->Draw("colz");

  return c123;

}

TCanvas *plotRowColLayer(){


  TCanvas *c5 = new TCanvas("c5", "c5", 50,50,800,800);
  c5->Divide(4,2, 0.01, 0.01, 0);

  c5->cd(1);
  gPad->SetLogy();
  hRowLayer1Side1->Draw();
  c5->cd(2);
  gPad->SetLogy();
  hRowLayer1Side2->Draw();
  c5->cd(3);
  gPad->SetLogy();
  hRowLayer2Side1->Draw();
  c5->cd(4);
  gPad->SetLogy();
  hRowLayer2Side2->Draw();
  c5->cd(5);
  gPad->SetLogy();
  hRow->Draw();
  c5->cd(6);
  gPad->SetLogy();
  hCol->Draw();
  c5->cd(7);
  gPad->SetLogy();
  hLayer->Draw();
  c5->cd(8);
  gPad->SetLogy();
  hHitPMT->Draw();

  return c5;

}

TCanvas *plotNhits(){
  TCanvas *c55 = new TCanvas("c55", "c5", 50,50,800,800);
  c55->Divide(3,3, 0.01, 0.01, 0);

  c55->cd(1);
  hnhits1->Draw();
  c55->cd(2);
  hngoodhits1->Draw();
  c55->cd(3);
  hngoodTDChits1->Draw();
  //c55->cd(4);
  //hnhits2->Draw();
  //c55->cd(5);
  //hngoodhits2->Draw();
  //c55->cd(6);
  //hngoodTDChits2->Draw();
  c55->cd(4);
  hngoodTDCpaddles->Draw();
  c55->cd(7);
  hnhits_ev->Draw();
  c55->cd(8);
  hngoodhits_ev->Draw();
  c55->cd(9);
  hngoodTDChits_ev->Draw();

  return c55;

}

TCanvas *plotTDC2d(){


  TCanvas *c6 = new TCanvas("c6", "c6", 50,50,800,800);
  c6->Divide(2,2, 0.01, 0.01, 0);

  c6->cd(1);
  h2d_RawLE->Draw();
  c6->cd(2);
  h2d_RawTE->Draw();
  c6->cd(3);
  h2d_RawTot->Draw();
  c6->cd(4);
  h2d_Mult->Draw();

  return c6;

}

auto plotXYZ(){

   TCanvas *c7 = new TCanvas("c7", "c7", 800,800);
   c7->Draw();
   TPad *p1 = new TPad("p1","p1",0.05,0.0,0.45,1.0);
   p1->Draw();
   p1->Divide(1,3);

   p1->cd(1);
   gPad->SetLogy();
   hHitX->Draw();

   p1->cd(2);
   gPad->SetLogy();
   hHitY->Draw();
   
   p1->cd(3);
   gPad->SetLogy();
   hHitZ->Draw();

   c7->cd(0);
   TPad *p2 = new TPad("p1","p1",0.55,0.0,0.95,1.0);
   p2->Draw();
   p2->Divide(2,1);

   p2->cd(1);
   gPad->SetLogz();
   hHitXY1->Draw("colz");
   p2->cd(2);
   gPad->SetLogz();
   hHitXY2->Draw("colz");


  return c7;

}

auto plotXYECalCDet(){

   TCanvas *c8 = new TCanvas("c8", "c7", 1200,1200);
   c8->Divide(3,4);

   c8->cd(1);
   hXECalCDet1->Draw("colz");

   c8->cd(2);
   hXECalCDet2->Draw("colz");
   
   c8->cd(3);
   hXECal->Draw();
   
   c8->cd(4);
   hXECal->Draw();
   
   c8->cd(5);
   hXDiffECalCDet1->Draw();

   c8->cd(6);
   hXDiffECalCDet2->Draw();

   c8->cd(7);
   hXPlusECalCDet1->Draw();

   c8->cd(8);
   hXPlusECalCDet2->Draw();
   
   c8->cd(9);
   hYECalCDet1->Draw("colz");
   
   c8->cd(10);
   hYECalCDet2->Draw("colz");

   
   c8->cd(11);
   hXYECal->Draw("colz");
   

  return c8;
}
