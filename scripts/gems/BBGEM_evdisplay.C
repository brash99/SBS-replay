#include "gmn_tree_allstripoutputs.C"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMarker.h"
#include "TCanvas.h"
#include "TPad.h"
#include <iostream>
#include <fstream>
#include "TStyle.h"
#include "TLegend.h"
#include "TAttLine.h"
#include "TAttMarker.h"
#include "TLine.h"
#include <vector>
#include <set>
#include <map>
#include "TClonesArray.h"
#include "TMath.h"
#include "TSystem.h"

void BBGEM_evdisplay( const char *rootfilename, const char *outputfilename="temp.root", int BBGEMconfig=3 ){

  //Now what do we want?
  // 1. We want 1D histograms of U/V strips fired, with height given by ADC sums on the strips. We want to do this with various cuts
  // We want to be able to indicate the "region of interest" on these plots along the direction measured by the strips (somewhat complicated) 
  // 2. We want 2D histograms of strips fired, color-coded by ADC sum
  // 3. We want reconstructed hit positions of 2D hits on tracks
  // 4. With our baseline strip and cluster outputs, we probably can't show the candidate 2D hit positions (future project for NIM paper)

  //BigBite GEM configuration 3:

  gStyle->SetTitleFont(42);
  gStyle->SetTitleSize(0.07);
  gStyle->SetTitleSize(0.07,"XYZ");
  gStyle->SetLabelSize(0.07,"XYZ");
  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetLabelFont(42,"XYZ");
  gROOT->ForceStyle();
  
  double lmargin = 0.18;
  double bmargin = 0.2;
  double rmargin = 0.05;
  double tmargin = 0.1;
  
  //This list of hardcoded parameters SHOULD give us what we need for an event display for GEM configuration 3:
  int nlayers=5, nmodules=8;
  
  vector<int> nstripu_module(nmodules), nstripv_module(nmodules);
  vector<double> uangle_mod(nmodules), vangle_mod(nmodules);
  vector<double> uoff_mod(nmodules), voff_mod(nmodules);
  vector<double> Lx_mod(nmodules), Ly_mod(nmodules);
  for( int i=0; i<4; i++ ){
    nstripu_module[i] = 3840;
    nstripv_module[i] = 3840;
    uangle_mod[i] = 150.0 * TMath::DegToRad();
    vangle_mod[i] = -150.0 * TMath::DegToRad();
    uoff_mod[i] = 0.0108;
    voff_mod[i] = 0.0108;

    Lx_mod[i] = 150.0;
    Ly_mod[i] = 40.0;

    nstripu_module[i+4] = 1280;
    nstripv_module[i+4] = 1536;
    uangle_mod[i+4] = 0.0;
    vangle_mod[i+4] = -90.0*TMath::DegToRad();
    uoff_mod[i+4] = 0.0;
    voff_mod[i+4] = 0.0;
    Lx_mod[i+4] = 51.2;
    Ly_mod[i+4] = 61.44;
  }

  TClonesArray *hStripU1D_by_module_all = new TClonesArray("TH1D",nmodules);
  TClonesArray *hStripU1D_by_module_keep = new TClonesArray("TH1D",nmodules);
  TClonesArray *hStripU1D_by_module_keep_lmax = new TClonesArray("TH1D",nmodules);
  TClonesArray *hStripU1D_by_module_ontrack = new TClonesArray("TH1D",nmodules);

  TClonesArray *hStripV1D_by_module_all = new TClonesArray("TH1D",nmodules);
  TClonesArray *hStripV1D_by_module_keep = new TClonesArray("TH1D",nmodules);
  TClonesArray *hStripV1D_by_module_keep_lmax = new TClonesArray("TH1D",nmodules);
  TClonesArray *hStripV1D_by_module_ontrack = new TClonesArray("TH1D",nmodules);
  
  for( int imod=0; imod<nmodules; imod++ ){

    cout << "module, nstripsu, nstripsv = " << imod << ", " << nstripu_module[imod] << ", " << nstripv_module[imod] << endl;
    TString histname;
    new( (*hStripU1D_by_module_all)[imod] ) TH1D( histname.Format("hStripU1D_all_mod%d", imod), histname.Format("Module %d all; U strip ; ADC sum",imod), nstripu_module[imod], -0.5, nstripu_module[imod]-0.5 );
    new( (*hStripU1D_by_module_keep)[imod] ) TH1D( histname.Format("hStripU1D_keep_mod%d", imod), "Passed zero sup ; strip ; ADC sum", nstripu_module[imod], -0.5, nstripu_module[imod]-0.5 );
    new( (*hStripU1D_by_module_keep_lmax)[imod] ) TH1D( histname.Format("hStripU1D_keep_lmax_mod%d", imod), "Passed local max threshold ; strip ; ADC sum", nstripu_module[imod], -0.5, nstripu_module[imod]-0.5 );
    new( (*hStripU1D_by_module_ontrack)[imod] ) TH1D( histname.Format("hStripU1D_ontrack_mod%d", imod), histname.Format("Module %d on track ; U strip ; ADC sum",imod), nstripu_module[imod], -0.5, nstripu_module[imod]-0.5 );

    new( (*hStripV1D_by_module_all)[imod] ) TH1D( histname.Format("hStripV1D_all_mod%d", imod), histname.Format("Module %d all ; V strip ; ADC sum",imod), nstripv_module[imod], -0.5, nstripv_module[imod]-0.5 );
    new( (*hStripV1D_by_module_keep)[imod] ) TH1D( histname.Format("hStripV1D_keep_mod%d", imod), "Passed zero sup ; strip ; ADC sum", nstripv_module[imod], -0.5, nstripv_module[imod]-0.5 );
    new( (*hStripV1D_by_module_keep_lmax)[imod] ) TH1D( histname.Format("hStripV1D_keep_lmax_mod%d", imod), "Passed local max threshold ; strip ; ADC sum", nstripv_module[imod], -0.5, nstripv_module[imod]-0.5 );
    new( (*hStripV1D_by_module_ontrack)[imod] ) TH1D( histname.Format("hStripV1D_ontrack_mod%d", imod), histname.Format("Module %d on track ; V strip ; ADC sum",imod), nstripv_module[imod], -0.5, nstripv_module[imod]-0.5 );

    ( (TH1D*) (*hStripU1D_by_module_all)[imod] )->SetLineColor(1);
    ( (TH1D*) (*hStripV1D_by_module_all)[imod] )->SetLineColor(1);

    ( (TH1D*) (*hStripU1D_by_module_ontrack)[imod] )->SetLineColor(2);
    ( (TH1D*) (*hStripV1D_by_module_ontrack)[imod] )->SetLineColor(2);

    ( (TH1D*) (*hStripU1D_by_module_ontrack)[imod] )->SetFillColor(2);
    ( (TH1D*) (*hStripV1D_by_module_ontrack)[imod] )->SetFillColor(2);
  }
  

  TChain *C = new TChain("T");

  C->Add(rootfilename);

  gmn_tree_allstripoutputs *T = new gmn_tree_allstripoutputs(C);

  long nevent = 0;

  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  //TCanvas *c1v = new TCanvas("c1v","c1v",1200,900);

  c1->Divide(2,nlayers,.001,.001);
  //c1->Divide(2,nlayers,.001,.001);
  
  c1->cd(9)->Divide(4,1,.001,.001);
  c1->cd(10)->Divide(4,1,.001,.001);

  bool keepgoing = true;
  while( C->GetEntry( nevent++ ) && keepgoing ){
    if( T->bb_tr_n > 0 && T->bb_gem_track_nhits[0] > 4 && fabs(T->bb_tr_vz[0])<0.07 && T->bb_gem_track_chi2ndf[0]<=15.0 && fabs(T->bb_tr_r_x[0]-0.9*T->bb_tr_r_th[0])<=0.3 &&
	T->bb_ps_e>0.15 && fabs(T->bb_etot_over_p[0]-1.0)<=0.25 ){
      //reset all histograms:
      for( int imod=0; imod<nmodules; imod++ ){
	( (TH1D*) (*hStripU1D_by_module_all)[imod] )->Reset();
	( (TH1D*) (*hStripV1D_by_module_all)[imod] )->Reset();
	( (TH1D*) (*hStripU1D_by_module_all)[imod] )->GetYaxis()->SetRangeUser(15,4000);
	( (TH1D*) (*hStripV1D_by_module_all)[imod] )->GetYaxis()->SetRangeUser(15,4000);
	( (TH1D*) (*hStripU1D_by_module_all)[imod] )->GetXaxis()->SetRangeUser(-10000,10000);
	( (TH1D*) (*hStripV1D_by_module_all)[imod] )->GetXaxis()->SetRangeUser(-10000,10000);
	( (TH1D*) (*hStripU1D_by_module_ontrack)[imod] )->Reset();
	( (TH1D*) (*hStripV1D_by_module_ontrack)[imod] )->Reset();
      }

      //      int ngoodU_mod(nmodules,0);
      // int ngoodV_mod(nmodules,0);
      
      for( int istrip=0; istrip<T->bb_gem_m0_strip_nstripsfired; istrip++ ){
	if( T->bb_gem_m0_strip_IsU[istrip] ){
	  ( (TH1D*) (*hStripU1D_by_module_all)[0] )->Fill( T->bb_gem_m0_strip_istrip[istrip], T->bb_gem_m0_strip_ADCmax[istrip] );
	  if( T->bb_gem_m0_strip_ontrack[istrip] && T->bb_gem_m0_strip_itrack[istrip]==0 ){
	    ( (TH1D*) (*hStripU1D_by_module_ontrack)[0] )->Fill( T->bb_gem_m0_strip_istrip[istrip], T->bb_gem_m0_strip_ADCmax[istrip] );
	  }
	} else {
	  ( (TH1D*) (*hStripV1D_by_module_all)[0] )->Fill( T->bb_gem_m0_strip_istrip[istrip], T->bb_gem_m0_strip_ADCmax[istrip] );
	  if( T->bb_gem_m0_strip_ontrack[istrip] && T->bb_gem_m0_strip_itrack[istrip]==0 ){
	    ( (TH1D*) (*hStripV1D_by_module_ontrack)[0] )->Fill( T->bb_gem_m0_strip_istrip[istrip], T->bb_gem_m0_strip_ADCmax[istrip] );
	  }
	}
      }

      for( int istrip=0; istrip<T->bb_gem_m1_strip_nstripsfired; istrip++ ){
	if( T->bb_gem_m1_strip_IsU[istrip] ){
	  ( (TH1D*) (*hStripU1D_by_module_all)[1] )->Fill( T->bb_gem_m1_strip_istrip[istrip], T->bb_gem_m1_strip_ADCmax[istrip] );
	  if( T->bb_gem_m1_strip_ontrack[istrip] &&  T->bb_gem_m1_strip_itrack[istrip]==0){
	    ( (TH1D*) (*hStripU1D_by_module_ontrack)[1] )->Fill( T->bb_gem_m1_strip_istrip[istrip], T->bb_gem_m1_strip_ADCmax[istrip] );
	  }
	} else {
	  ( (TH1D*) (*hStripV1D_by_module_all)[1] )->Fill( T->bb_gem_m1_strip_istrip[istrip], T->bb_gem_m1_strip_ADCmax[istrip] );
	  if( T->bb_gem_m1_strip_ontrack[istrip] &&  T->bb_gem_m1_strip_itrack[istrip]==0){
	    ( (TH1D*) (*hStripV1D_by_module_ontrack)[1] )->Fill( T->bb_gem_m1_strip_istrip[istrip], T->bb_gem_m1_strip_ADCmax[istrip] );
	  }
	}
      }

      for( int istrip=0; istrip<T->bb_gem_m2_strip_nstripsfired; istrip++ ){
	if( T->bb_gem_m2_strip_IsU[istrip] ){
	  ( (TH1D*) (*hStripU1D_by_module_all)[2] )->Fill( T->bb_gem_m2_strip_istrip[istrip], T->bb_gem_m2_strip_ADCmax[istrip] );
	  if( T->bb_gem_m2_strip_ontrack[istrip] &&  T->bb_gem_m2_strip_itrack[istrip]==0){
	    ( (TH1D*) (*hStripU1D_by_module_ontrack)[2] )->Fill( T->bb_gem_m2_strip_istrip[istrip], T->bb_gem_m2_strip_ADCmax[istrip] );
	  }
	} else {
	  ( (TH1D*) (*hStripV1D_by_module_all)[2] )->Fill( T->bb_gem_m2_strip_istrip[istrip], T->bb_gem_m2_strip_ADCmax[istrip] );
	  if( T->bb_gem_m2_strip_ontrack[istrip] &&  T->bb_gem_m2_strip_itrack[istrip]==0){
	    ( (TH1D*) (*hStripV1D_by_module_ontrack)[2] )->Fill( T->bb_gem_m2_strip_istrip[istrip], T->bb_gem_m2_strip_ADCmax[istrip] );
	  }
	}
      }

      for( int istrip=0; istrip<T->bb_gem_m3_strip_nstripsfired; istrip++ ){
	if( T->bb_gem_m3_strip_IsU[istrip] ){
	  ( (TH1D*) (*hStripU1D_by_module_all)[3] )->Fill( T->bb_gem_m3_strip_istrip[istrip], T->bb_gem_m3_strip_ADCmax[istrip] );
	  if( T->bb_gem_m3_strip_ontrack[istrip] &&  T->bb_gem_m3_strip_itrack[istrip]==0){
	    ( (TH1D*) (*hStripU1D_by_module_ontrack)[3] )->Fill( T->bb_gem_m3_strip_istrip[istrip], T->bb_gem_m3_strip_ADCmax[istrip] );
	  }
	} else {
	  ( (TH1D*) (*hStripV1D_by_module_all)[3] )->Fill( T->bb_gem_m3_strip_istrip[istrip], T->bb_gem_m3_strip_ADCmax[istrip] );
	  if( T->bb_gem_m3_strip_ontrack[istrip] &&  T->bb_gem_m3_strip_itrack[istrip]==0){
	    ( (TH1D*) (*hStripV1D_by_module_ontrack)[3] )->Fill( T->bb_gem_m3_strip_istrip[istrip], T->bb_gem_m3_strip_ADCmax[istrip] );
	  }
	}
      }

      for( int istrip=0; istrip<T->bb_gem_m4_strip_nstripsfired; istrip++ ){
	if( T->bb_gem_m4_strip_IsU[istrip] ){
	  ( (TH1D*) (*hStripU1D_by_module_all)[4] )->Fill( T->bb_gem_m4_strip_istrip[istrip], T->bb_gem_m4_strip_ADCmax[istrip] );
	  if( T->bb_gem_m4_strip_ontrack[istrip] &&  T->bb_gem_m4_strip_itrack[istrip]==0){
	    ( (TH1D*) (*hStripU1D_by_module_ontrack)[4] )->Fill( T->bb_gem_m4_strip_istrip[istrip], T->bb_gem_m4_strip_ADCmax[istrip] );
	  }
	} else {
	  ( (TH1D*) (*hStripV1D_by_module_all)[4] )->Fill( T->bb_gem_m4_strip_istrip[istrip], T->bb_gem_m4_strip_ADCmax[istrip] );
	  if( T->bb_gem_m4_strip_ontrack[istrip] &&  T->bb_gem_m4_strip_itrack[istrip]==0){
	    ( (TH1D*) (*hStripV1D_by_module_ontrack)[4] )->Fill( T->bb_gem_m4_strip_istrip[istrip], T->bb_gem_m4_strip_ADCmax[istrip] );
	  }
	}
      }

      for( int istrip=0; istrip<T->bb_gem_m5_strip_nstripsfired; istrip++ ){
	if( T->bb_gem_m5_strip_IsU[istrip] ){
	  ( (TH1D*) (*hStripU1D_by_module_all)[5] )->Fill( T->bb_gem_m5_strip_istrip[istrip], T->bb_gem_m5_strip_ADCmax[istrip] );
	  if( T->bb_gem_m5_strip_ontrack[istrip] &&  T->bb_gem_m5_strip_itrack[istrip]==0){
	    ( (TH1D*) (*hStripU1D_by_module_ontrack)[5] )->Fill( T->bb_gem_m5_strip_istrip[istrip], T->bb_gem_m5_strip_ADCmax[istrip] );
	  }
	} else {
	  ( (TH1D*) (*hStripV1D_by_module_all)[5] )->Fill( T->bb_gem_m5_strip_istrip[istrip], T->bb_gem_m5_strip_ADCmax[istrip] );
	  if( T->bb_gem_m5_strip_ontrack[istrip] &&  T->bb_gem_m5_strip_itrack[istrip]==0){
	    ( (TH1D*) (*hStripV1D_by_module_ontrack)[5] )->Fill( T->bb_gem_m5_strip_istrip[istrip], T->bb_gem_m5_strip_ADCmax[istrip] );
	  }
	}
      }

      for( int istrip=0; istrip<T->bb_gem_m6_strip_nstripsfired; istrip++ ){
	if( T->bb_gem_m6_strip_IsU[istrip] ){
	  ( (TH1D*) (*hStripU1D_by_module_all)[6] )->Fill( T->bb_gem_m6_strip_istrip[istrip], T->bb_gem_m6_strip_ADCmax[istrip] );
	  if( T->bb_gem_m6_strip_ontrack[istrip] &&  T->bb_gem_m6_strip_itrack[istrip]==0){
	    ( (TH1D*) (*hStripU1D_by_module_ontrack)[6] )->Fill( T->bb_gem_m6_strip_istrip[istrip], T->bb_gem_m6_strip_ADCmax[istrip] );
	  }
	} else {
	  ( (TH1D*) (*hStripV1D_by_module_all)[6] )->Fill( T->bb_gem_m6_strip_istrip[istrip], T->bb_gem_m6_strip_ADCmax[istrip] );
	  if( T->bb_gem_m6_strip_ontrack[istrip] &&  T->bb_gem_m6_strip_itrack[istrip]==0){
	    ( (TH1D*) (*hStripV1D_by_module_ontrack)[6] )->Fill( T->bb_gem_m6_strip_istrip[istrip], T->bb_gem_m6_strip_ADCmax[istrip] );
	  }
	}
      }

      for( int istrip=0; istrip<T->bb_gem_m7_strip_nstripsfired; istrip++ ){
	if( T->bb_gem_m7_strip_IsU[istrip] ){
	  ( (TH1D*) (*hStripU1D_by_module_all)[7] )->Fill( T->bb_gem_m7_strip_istrip[istrip], T->bb_gem_m7_strip_ADCmax[istrip] );
	  if( T->bb_gem_m7_strip_ontrack[istrip] &&  T->bb_gem_m7_strip_itrack[istrip]==0){
	    ( (TH1D*) (*hStripU1D_by_module_ontrack)[7] )->Fill( T->bb_gem_m7_strip_istrip[istrip], T->bb_gem_m7_strip_ADCmax[istrip] );
	  }
	} else {
	  ( (TH1D*) (*hStripV1D_by_module_all)[7] )->Fill( T->bb_gem_m7_strip_istrip[istrip], T->bb_gem_m7_strip_ADCmax[istrip] );
	  if( T->bb_gem_m7_strip_ontrack[istrip] &&  T->bb_gem_m7_strip_itrack[istrip]==0){
	    ( (TH1D*) (*hStripV1D_by_module_ontrack)[7] )->Fill( T->bb_gem_m7_strip_istrip[istrip], T->bb_gem_m7_strip_ADCmax[istrip] );
	  }
	}
      }


      
      for( int imod=0; imod<4; imod++ ){

	lmargin = 0.11;
	c1->cd(1+2*imod)->SetLogy(1);

	gPad->SetLeftMargin(lmargin);
	gPad->SetRightMargin(rmargin);
	gPad->SetBottomMargin(bmargin);
	gPad->SetTopMargin(tmargin);

	( (TH1D*) (*hStripU1D_by_module_all)[imod] )->Draw("hist");

	double Umean = ( (TH1D*) (*hStripU1D_by_module_ontrack)[imod] )->GetMean();

	if( Umean > 0.1 && ( (TH1D*) (*hStripU1D_by_module_ontrack)[imod] )->GetEntries() > 0 ){
	  ( (TH1D*) (*hStripU1D_by_module_all)[imod] )->GetXaxis()->SetRangeUser(Umean-100,Umean+100);
	}
	
	( (TH1D*) (*hStripU1D_by_module_ontrack)[imod] )->Draw("histsame");
      
	c1->cd(2+2*imod)->SetLogy(1);

	gPad->SetLeftMargin(lmargin);
	gPad->SetRightMargin(rmargin);
	gPad->SetBottomMargin(bmargin);
	gPad->SetTopMargin(tmargin);
	
	( (TH1D*) (*hStripV1D_by_module_all)[imod] )->Draw("hist");
	
	
	double Vmean = ( (TH1D*) (*hStripV1D_by_module_ontrack)[imod] )->GetMean();

	if( Vmean > 0.1 && ( (TH1D*) (*hStripV1D_by_module_ontrack)[imod] )->GetEntries() > 0 ){
	  ( (TH1D*) (*hStripV1D_by_module_all)[imod] )->GetXaxis()->SetRangeUser(Vmean-100,Vmean+100);
	}

	( (TH1D*) (*hStripV1D_by_module_ontrack)[imod] )->Draw("histsame");
	
	lmargin = 0.24;
	
	c1->cd(9)->cd(imod+1)->SetLogy(1);

	gPad->SetLeftMargin(lmargin);
	gPad->SetRightMargin(rmargin);
	gPad->SetBottomMargin(bmargin);
	gPad->SetTopMargin(tmargin);
	
	( (TH1D*) (*hStripU1D_by_module_all)[imod+4] )->Draw("hist");

	Umean = ( (TH1D*) (*hStripU1D_by_module_ontrack)[imod+4] )->GetMean();

	if( Umean > 0.1 && ( (TH1D*) (*hStripU1D_by_module_ontrack)[imod+4] )->GetEntries() > 0 ){
	  ( (TH1D*) (*hStripU1D_by_module_all)[imod+4] )->GetXaxis()->SetRangeUser(Umean-25,Umean+25);
	}
	
	( (TH1D*) (*hStripU1D_by_module_ontrack)[imod+4] )->Draw("histsame");
	c1->cd(10)->cd(imod+1)->SetLogy(1);

	gPad->SetLeftMargin(lmargin);
	gPad->SetRightMargin(rmargin);
	gPad->SetBottomMargin(bmargin);
	gPad->SetTopMargin(tmargin);
	
	( (TH1D*) (*hStripV1D_by_module_all)[imod+4] )->Draw("hist");

	Vmean = ( (TH1D*) (*hStripV1D_by_module_ontrack)[imod+4] )->GetMean();

	if( Vmean > 0.1 && ( (TH1D*) (*hStripV1D_by_module_ontrack)[imod+4] )->GetEntries() > 0 ){
	  ( (TH1D*) (*hStripV1D_by_module_all)[imod+4] )->GetXaxis()->SetRangeUser(Vmean-25,Vmean+25);
	}
	
	( (TH1D*) (*hStripV1D_by_module_ontrack)[imod+4] )->Draw("histsame");
	
	
      }
      gPad->Modified();
      //c1->Update();
      c1->Update();
      gSystem->ProcessEvents();

      std::cout << "next event? (hit enter to continue or enter q to quit):";
      TString reply;
      reply.ReadLine(cin,kFALSE);

      if( reply.BeginsWith("q") ) keepgoing = false;
    }    

  }

}
