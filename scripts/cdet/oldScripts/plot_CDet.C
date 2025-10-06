#include <TSystem.h>
#include <TChain.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>


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

const int num_bad = 12;

const int bad_channels[] = {
        604, 673, 674, 675, 676, 677, 678, 679, 680, 681, 682, 728
};


using namespace std;

bool check_bad(int pmt, bool suppress_bad) {
        bool flag = false;
        if (!suppress_bad) return flag;
        for (int i=0;i<num_bad;i++) {
                if (pmt == bad_channels[i])  flag = true;
        }
        return flag;
}

void plot_CDet(Int_t nrun=3311,Int_t nev=-1,
       Double_t LeMin = 0.0, Double_t LeMax = 60.0,
        Double_t TotMin = 2.0, Double_t TotMax = 60.0,
        Int_t nhitcutlow1 = 1, Int_t nhitcuthigh1 = 100,
        Int_t nhitcutlow2 = 1, Int_t nhitcuthigh2 = 100,
        Double_t XDiffCut = 0.073,
        bool suppress_bad = false){

  gStyle->SetPalette(1,0);
  gStyle->SetOptStat(1000011);
  gStyle->SetOptFit(11);
  gStyle->SetTitleOffset(1.,"Y");
  gStyle->SetTitleOffset(.7,"X");
  gStyle->SetLabelSize(0.04,"XY");
  gStyle->SetTitleSize(0.06,"XY");
  gStyle->SetPadLeftMargin(0.12);
   //
  TString basename="gep5_replayed_nogems";
  TString basename2="50k_events.root";
  TString Chainroot;
  Chainroot = Form("~/sbs/Rootfiles/%s_%d_%s",basename.Data(),nrun,basename2.Data());
  TChain *fchain = new TChain("T");
  fchain->Add(Chainroot);
  cout << "check file " << Chainroot << endl;
  TString fullname = Form("%s_%d_%s",basename.Data(),nrun,basename2.Data());
  TString outputhist;
  outputhist= "hist/"+fullname+"_hist.root";
  TObjArray HList(0);
  static const Int_t MaxHit=10000;
  fchain->SetBranchStatus("*",0);
   Int_t ecalHits;
   fchain->SetBranchAddress("Ndata.earm.ecal.a_amp_p",&ecalHits) ;
   fchain->SetBranchStatus("Ndata.earm.ecal.a_amp_p",1) ;
    Double_t ecal_amp[MaxHit]; // raw adc
     fchain->SetBranchAddress("earm.ecal.a_amp_p",&ecal_amp) ;
     fchain->SetBranchStatus("earm.ecal.a_amp_p",1) ;
    Double_t ecal_int[MaxHit]; // raw adc
     fchain->SetBranchAddress("earm.ecal.a_p",&ecal_int) ;
     fchain->SetBranchStatus("earm.ecal.a_p",1) ;
   Double_t ecal_id[MaxHit]; // raw adc
     fchain->SetBranchAddress("earm.ecal.adcelemID",&ecal_id) ;
     fchain->SetBranchStatus("earm.ecal.adcelemID",1) ;
   Double_t ecal_row[MaxHit]; // raw adc
     fchain->SetBranchAddress("earm.ecal.adcrow",&ecal_row) ;
     fchain->SetBranchStatus("earm.ecal.adcrow",1) ;
    Double_t ecal_col[MaxHit]; // raw adc
     fchain->SetBranchAddress("earm.ecal.adccol",&ecal_col) ;
     fchain->SetBranchStatus("earm.ecal.adccol",1) ;
   Double_t ecal_x[MaxHit]; // raw adc
     fchain->SetBranchAddress("earm.ecal.adcxpos",&ecal_x) ;
     fchain->SetBranchStatus("earm.ecal.adcxpos",1) ;
   Double_t ecal_y[MaxHit]; // raw adc
     fchain->SetBranchAddress("earm.ecal.adcypos",&ecal_y) ;
     fchain->SetBranchStatus("earm.ecal.adcypos",1) ;
   Double_t ecal_xpos; // raw adc
     fchain->SetBranchAddress("earm.ecal.x",&ecal_xpos) ;
     fchain->SetBranchStatus("earm.ecal.x",1) ;
   Double_t ecal_ypos; // raw adc
     fchain->SetBranchAddress("earm.ecal.y",&ecal_ypos) ;
     fchain->SetBranchStatus("earm.ecal.y",1) ;
//
    // enable branches
    fchain->SetBranchStatus("earm.cdet.*",1);
    fchain->SetBranchStatus("earm.ecal.*",1);

    fchain->SetBranchAddress("earm.cdet.tdc_mult",TCDet::TDCmult);

    fchain->SetBranchAddress("earm.cdet.hits.TDCelemID",TCDet::RawElID);
    fchain->SetBranchAddress("earm.cdet.hits.t",TCDet::RawElLE);
    fchain->SetBranchAddress("earm.cdet.hits.t_te",TCDet::RawElTE);
    fchain->SetBranchAddress("earm.cdet.hits.t_tot",TCDet::RawElTot);

    fchain->SetBranchAddress("earm.cdet.hit.pmtnum",TCDet::GoodElID);
    fchain->SetBranchAddress("earm.cdet.hit.tdc_le",TCDet::GoodElLE);
    fchain->SetBranchAddress("earm.cdet.hit.tdc_te",TCDet::GoodElTE);
    fchain->SetBranchAddress("earm.cdet.hit.tdc_tot",TCDet::GoodElTot);

    fchain->SetBranchAddress("earm.cdet.hit.xhit",TCDet::GoodX);
    fchain->SetBranchAddress("earm.cdet.hit.yhit",TCDet::GoodY);
    fchain->SetBranchAddress("earm.cdet.hit.zhit",TCDet::GoodZ);

    fchain->SetBranchAddress("earm.cdet.hit.row",TCDet::GoodCol);
    fchain->SetBranchAddress("earm.cdet.hit.col",TCDet::GoodRow);
    fchain->SetBranchAddress("earm.cdet.hit.layer",TCDet::GoodLayer);


    fchain->SetBranchAddress("earm.cdet.nhits",&TCDet::nhits);
    fchain->SetBranchAddress("earm.cdet.ngoodhits",&TCDet::ngoodhits);
    fchain->SetBranchAddress("earm.cdet.ngoodTDChits",&TCDet::ngoodTDChits);
    //fchain->SetBranchAddress("earm.ecal.x",&TCDet::GoodECalX);
    //fchain->SetBranchAddress("earm.ecal.y",&TCDet::GoodECalY);

    // enable vector size branches
    fchain->SetBranchAddress("Ndata.earm.cdet.tdc_mult",&TCDet::NdataMult);
    fchain->SetBranchAddress("Ndata.earm.cdet.hits.TDCelemID",&TCDet::NdataRawElID);
    fchain->SetBranchAddress("Ndata.earm.cdet.hits.t",&TCDet::NdataRawElLE);
    fchain->SetBranchAddress("Ndata.earm.cdet.hits.t_te",&TCDet::NdataRawElTE);
    fchain->SetBranchAddress("Ndata.earm.cdet.hits.t_tot",&TCDet::NdataRawElTot);

    fchain->SetBranchAddress("Ndata.earm.cdet.hit.pmtnum",&TCDet::NdataGoodElID);
    fchain->SetBranchAddress("Ndata.earm.cdet.hit.tdc_le",&TCDet::NdataGoodElLE);
    fchain->SetBranchAddress("Ndata.earm.cdet.hit.tdc_te",&TCDet::NdataGoodElTE);
    fchain->SetBranchAddress("Ndata.earm.cdet.hit.tdc_tot",&TCDet::NdataGoodElTot);

    fchain->SetBranchAddress("Ndata.earm.cdet.hit.xhit",&TCDet::NdataGoodX);
    fchain->SetBranchAddress("Ndata.earm.cdet.hit.yhit",&TCDet::NdataGoodY);
    fchain->SetBranchAddress("Ndata.earm.cdet.hit.zhit",&TCDet::NdataGoodZ);

    fchain->SetBranchAddress("Ndata.earm.cdet.hit.row",&TCDet::NdataGoodCol);
    fchain->SetBranchAddress("Ndata.earm.cdet.hit.col",&TCDet::NdataGoodRow);
    fchain->SetBranchAddress("Ndata.earm.cdet.hit.layer",&TCDet::NdataGoodLayer);




   //
     TH2D* h_ecalxy;
     h_ecalxy = new TH2D("","",60,-0.8,0.8,300,-1.5,1.5);
      TCanvas *c_ecalxy= new TCanvas("c_ecalxy","c_ecalxy",1000,700);
     TH2D* h_ecalrc;
     h_ecalrc = new TH2D("","",28,-1,27,69,0,69);
      c_ecalxy->Divide(1,1);
      
     TH2D* h_ecal_amp_col[69];
     for (Int_t i=0;i<69;i++) {
       h_ecal_amp_col[i] = new TH2D(Form("h_ecal_amp_col_r%d",i),Form("Row %d; Col; Amp",i),28,-1,27,200,0,200);
       HList.Add(h_ecal_amp_col[i]);
     }
      
//
      Int_t debug=1;
    vector<Int_t> goodhits;
  Long64_t nentries = fchain->GetEntries();
   cout << " Total Entry = " << nentries << endl;
  for (int ii = 0; ii < nentries; ii++) {
    fchain->GetEntry(ii);
    if (ii%50000==0) cout << " Entry = " << ii << endl;
    goodhits.clear();

        for (Int_t i=0;i<ecalHits;i++) {
	  if(ecal_row[i] == 0 || ecal_row[i] == 1 || ecal_row[i] == 2 || 
		ecal_row[i] == 6 || ecal_row[i] == 7 || ecal_row[i] == 8 || ecal_row[i] == 12 || 
		ecal_row[i] == 13 || ecal_row[i] == 14 || ecal_row[i] == 65 || ecal_row[i] == 64 || 
		ecal_row[i] == 63 || ecal_row[i] == 59 || ecal_row[i] == 58 || ecal_row[i] == 57) ecal_col[i]--;
	}

	//
        for (Int_t i=0;i<ecalHits;i++) {
	  goodhits.push_back(0); 
            for (Int_t j=0;j<ecalHits;j++) {
	      if (i!=j&&TMath::Abs( ecal_row[i]-ecal_row[j]) <=1 && TMath::Abs( ecal_col[i]-ecal_col[j])<=1) {
		 goodhits[i] = 1;
		 //		 cout << i << " " << j << " " << TMath::Abs( ecal_row[i]-ecal_row[j]) << " " << TMath::Abs( ecal_col[i]-ecal_col[j])<< endl;
	      }
	    }
	}
	//
    if (debug ==1 ) {
      Int_t multhits=0;
      for (UInt_t j=0;j<goodhits.size()-1;j++) {
	if (goodhits[j] == 1 && goodhits[j+1] == 1) multhits++;
      }
    if (goodhits.size() > 10 && multhits > 10) {
    for (Int_t i=0;i<ecalHits;i++) {
      //cout << i << " id =  " << ecal_id[i] << " row =  " << ecal_row[i] << " col =  " << ecal_col[i] << " x =  " << ecal_x[i] << " y =  " << ecal_y[i] << " amp = " << ecal_amp[i]<< " int = " << ecal_int[i] << " good hits = " << goodhits[i] << endl;
      //cout << "ECal: " << ecal_xpos << " " << ecal_ypos << endl;
      //if ( goodhits[i] ==1) h_ecalxy->Fill(ecal_y[i],-ecal_x[i],ecal_amp[i]);
      if ( goodhits[i] ==1) h_ecalxy->Fill(ecal_ypos,ecal_xpos,ecal_amp[i]);
      if ( goodhits[i] ==1) h_ecalrc->Fill(ecal_col[i],ecal_row[i],ecal_amp[i]);
    }
      //c_ecalxy->cd(1);
      //h_ecalxy->Draw("colz");
      //	c_ecalxy->cd(2);
      //h_ecalrc->Draw("colz");
      //c_ecalxy->Update();
      //gPad->WaitPrimitive();
      //h_ecalxy->Reset();
      //h_ecalrc->Reset();
      //    Int_t t;
      //cout << " Hit space bar for next event" << endl;
      //cin >> debug;
    }}
    
    if (goodhits.size() > 9) {
     for (UInt_t i=0;i<goodhits.size();i++) {
       Int_t row = ecal_row[i];
       //      if ( goodhits[i] ==1) h_ecal_amp_col[row]->Fill(ecal_col[i],ecal_amp[i]);
     }
    }

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
                && (TCDet::GoodX[el]-ecal_xpos) <= XDiffCut && (TCDet::GoodX[el]-ecal_xpos) >= -1.0*XDiffCut
                && ecal_xpos > -1.5 && ecal_xpos < 1.5
                && ecal_ypos > -0.8 && ecal_ypos < 0.8
	) {

          if ( !check_bad(TCDet::GoodElID[el], suppress_bad) ) {
           if ( TCDet::GoodElID[el] < NumCDetPaddles )  {

            int sbsrown = (Int_t)TCDet::GoodRow[el];
            int sbscoln = (Int_t)TCDet::GoodCol[el];
            int mylayern = sbscoln/2;
            int mypaddlen = sbscoln*672 + sbsrown;
            cout << "Hit number " << el << ":    Paddle = " << mypaddlen << " Row = " << sbsrown  << " Col = " << sbscoln  << " hits = " << TCDet::ngoodTDChits_paddles[mypaddlen] << endl;
            cout << "el = " << el << " Good ID = " << TCDet::GoodElID[el] << " Good le = " << 
              TCDet::GoodElLE[el] << " Good te = " << TCDet::GoodElTE[el] << " Good tot = " << 
              TCDet::GoodElTot[el] << " ECal X = " << ecal_xpos << " CDet X = " << TCDet::GoodX[el] << endl;
	    //cout << "Ecal: " << ecal_xpos << " " << ecal_ypos << endl;

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


   for(Int_t el=0; el<TCDet::NdataGoodElID; el++){
        if (TCDet::GoodElLE[el] >= LeMin/TDC_calib_to_ns && TCDet::GoodElLE[el] <= LeMax/TDC_calib_to_ns
                && TCDet::GoodElTot[el] >= TotMin/TDC_calib_to_ns && TCDet::GoodElTot[el] <= TotMax/TDC_calib_to_ns
                && TCDet::GoodX[el] < xcut && TCDet::TDCmult[el] < TDCmult_cut
                && ngoodTDChitsc1 >= nhitcutlow1  && ngoodTDChitsc2 >= nhitcutlow2
                && ngoodTDChitsc1 <= nhitcuthigh1 && ngoodTDChitsc2 <= nhitcuthigh2
                && (TCDet::GoodX[el]-ecal_xpos) <= XDiffCut && (TCDet::GoodX[el]-ecal_xpos) >= -1.0*XDiffCut
                && ecal_xpos > -1.5 && ecal_xpos < 1.5
                && ecal_ypos > -0.8 && ecal_ypos < 0.8
        ) {

          if ( !check_bad(TCDet::GoodElID[el], suppress_bad) ) {
           if ( TCDet::GoodElID[el] < NumCDetPaddles )  {
            int sbscol = (Int_t)TCDet::GoodCol[el];
            int myside = sbscol%2;
            int mylayer = sbscol/2;
	    h_ecalxy->Fill(TCDet::GoodY[el],TCDet::GoodX[el]);

	   } // good TDC channel
	  } // good channel 
	} // good cdet hit
    } // cdet hit loop

      c_ecalxy->cd(1);
      h_ecalxy->Draw("colz");
      //        c_ecalxy->cd(2);
      //h_ecalrc->Draw("colz");
      c_ecalxy->Update();
      gPad->WaitPrimitive();
      h_ecalxy->Reset();
      h_ecalrc->Reset();
      //    Int_t t;
      cout << " Hit space bar for next event" << endl;


  } // nentries loop
  //
  TFile hsimc(outputhist,"recreate");
  HList.Write();
  //
}
