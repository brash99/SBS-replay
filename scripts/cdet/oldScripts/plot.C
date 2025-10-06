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
using namespace std;

void plot(Int_t nrun=2885,Int_t nev=-1){
  gStyle->SetPalette(1,0);
  gStyle->SetOptStat(1000011);
  gStyle->SetOptFit(11);
  gStyle->SetTitleOffset(1.,"Y");
  gStyle->SetTitleOffset(.7,"X");
  gStyle->SetLabelSize(0.04,"XY");
  gStyle->SetTitleSize(0.06,"XY");
  gStyle->SetPadLeftMargin(0.12);
   //
  TString basename="cdet";
  TString Chainroot;
  Chainroot = Form("~/sbs/Rootfiles/%s_%d_%d.root",basename.Data(),nrun,nev);
  TChain *fchain = new TChain("T");
  fchain->Add(Chainroot);
  cout << "check file " << Chainroot << endl;
  Int_t npart = 1;
  Chainroot = Form("~/sbs/Rootfiles/%s_%d_%d_%d.root",basename.Data(),nrun,nev,npart);
  while (gSystem->FindFile(".",Chainroot)) { 
   cout << " add file " << Chainroot << endl;
  fchain->Add(Chainroot);
  npart++;
  Chainroot = Form("~/sbs/Rootfiles/%s_%d_%d_%d.root",basename.Data(),nrun,nev,npart);
  }
  TString fullname = Form("%s_%d",basename.Data(),nrun);
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
   //
     TH2D* h_ecalxy;
     h_ecalxy = new TH2D("","",60,-1.2,1.2,80,-1.6,1.6);
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
	  if(ecal_row[i] == 0 || ecal_row[i] == 1 || ecal_row[i] == 2 || ecal_row[i] == 6 || ecal_row[i] == 7 || ecal_row[i] == 8 || ecal_row[i] == 12 || ecal_row[i] == 13 || ecal_row[i] == 14 || ecal_row[i] == 65 || ecal_row[i] == 64 || ecal_row[i] == 63 || ecal_row[i] == 59 || ecal_row[i] == 58 || ecal_row[i] == 57) ecal_col[i]--;
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
      //      cout << i << " id =  " << ecal_id[i] << " row =  " << ecal_row[i] << " col =  " << ecal_col[i] << " x =  " << ecal_x[i] << " y =  " << ecal_y[i] << " amp = " << ecal_amp[i]<< " int = " << ecal_int[i] << " good hits = " << goodhits[i] << endl;
      if ( goodhits[i] ==1) h_ecalxy->Fill(ecal_y[i],-ecal_x[i],ecal_amp[i]);
      if ( goodhits[i] ==1) h_ecalrc->Fill(ecal_col[i],ecal_row[i],ecal_amp[i]);
    }
	c_ecalxy->cd(1);
      h_ecalxy->Draw("colz");
      //	c_ecalxy->cd(2);
      //h_ecalrc->Draw("colz");
      c_ecalxy->Update();
      gPad->WaitPrimitive();
      h_ecalxy->Reset();
      h_ecalrc->Reset();
      //    Int_t t;
      cout << " Hit space bar for next event" << endl;
      //cin >> debug;
    }}
    
    if (goodhits.size() > 9) {
     for (UInt_t i=0;i<goodhits.size();i++) {
       Int_t row = ecal_row[i];
       //      if ( goodhits[i] ==1) h_ecal_amp_col[row]->Fill(ecal_col[i],ecal_amp[i]);
     }
    }
    
  }
  //
   TFile hsimc(outputhist,"recreate");
  HList.Write();
  //
}
