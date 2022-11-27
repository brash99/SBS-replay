#include <TSystem.h>
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
#include <TText.h>
#include <TF1.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
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

void fit_bpm() {
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(0);
 gStyle->SetTitleOffset(1.2,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.16);
 TString DataFile="bpmx.dat";
 vector<Double_t> bpm;
 vector<Double_t> rb_ax;
 vector<Double_t> rb_bx;
 //
     ifstream file_Data(DataFile);
      if (file_Data.is_open()) {
	Int_t nr;
	Int_t ni;
	file_Data >> nr;
	file_Data >> ni;
        cout << " number of rows= " << nr << " number of items per row = " << ni << endl;
	for (Int_t ny=0;ny<nr;ny++)  {
		  Double_t temp;
	for (Int_t nn=0;nn<ni;nn++)  {
	  file_Data >>  temp;
	  if (nn==0) bpm.push_back(temp);
	  if (nn==1) rb_ax.push_back(1000*temp);
	  if (nn==2) rb_bx.push_back(1000*temp);
	    }
	cout << ny << " " << bpm[ny] << " " << rb_ax[ny] << " " << rb_bx[ny] << endl; 
	}
      }
        TGraph *gAx; 
        TGraph *gBx; 
      	gAx= new TGraph(rb_ax.size(), &rb_ax[0], &bpm[0]);
	gAx->SetMarkerStyle(23);
	gAx->SetTitle(" Bpmx vs Ax");
	gAx->GetXaxis()->SetTitle("Ax");
	gAx->GetYaxis()->SetTitle("BPMx");
      	gBx= new TGraph(rb_bx.size(), &rb_bx[0], &bpm[0]);
	gBx->SetMarkerStyle(23);
	gBx->SetTitle(" Bpmx vs Bx");
	gBx->GetXaxis()->SetTitle("Bx");
	gBx->GetYaxis()->SetTitle("BPMx");
	//
	TCanvas *cBPMx = new TCanvas("bpmx","bpmx",700,700);
	cBPMx->Divide(1,2);
	cBPMx->cd(1);
	gAx->Draw("AP");
	gAx->Fit("pol1","Q","");
	 TF1 *fit = gAx->GetFunction("pol1");
	 Double_t slope = fit->GetParameter(1);
	 Double_t inter = fit->GetParameter(0);
	 TString ttemp=Form(" Slope  = %3.2f  Int  = %3.2f ",slope,inter);
	 gAx->SetTitle(ttemp);
	cBPMx->cd(2);
	gBx->Draw("AP");
	gBx->Fit("pol1","Q","");
	 fit = gBx->GetFunction("pol1");
	 slope = fit->GetParameter(1);
	 inter = fit->GetParameter(0);
	 ttemp=Form(" Slope  = %3.2f  Int  = %3.2f ",slope,inter);
	 gBx->SetTitle(ttemp);

	//


}
