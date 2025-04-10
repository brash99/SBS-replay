#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTreeFormula.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "genrp_tree_SBS_straight_through.C"

void DrawHist(TH1D *hist, vector<double> value, int flag = 0){

  //flag = 0 for calorimeter plots
  //flag = 1 for ADC plots
  //flag = 2 for time plots

  hist->Draw();

  TPaveText *pt;
  if(flag == 0){
    pt = new TPaveText(0.12,0.75,0.35,0.88,"ndc");
    pt->AddText(Form("Center = %g",value[0]));
    pt->AddText(Form("Width = %g",value[1]));
  }
  else if(flag == 1){
    pt = new TPaveText(0.55,0.75,0.75,0.88,"ndc");
    pt->AddText(Form("Threshold = %g",value[0]));
  }
  else if(flag == 2){
    pt = new TPaveText(0.12,0.75,0.35,0.88,"ndc");
    if(value.size() == 2){
      pt->AddText(Form("Mean = %g",value[0]));
      pt->AddText(Form("#sigma = %g",value[1]));
    }
    if(value.size() == 1){
      pt->AddText(Form("Cut = %g",value[0]));
    }
 }
  pt->SetFillColor(0);
  pt->Draw("same");


}

void FitGaus_FWHM( TH1D *htest, double thresh=0.5 ){
  int binmax = htest->GetMaximumBin();
  int binlow = binmax, binhigh = binmax;

  double max = htest->GetBinContent(binmax);

  while( htest->GetBinContent(binlow) >= thresh*max && binlow > 1 ){binlow--;}
  while( htest->GetBinContent(binhigh) >= thresh*max && binhigh < htest->GetNbinsX() ){ binhigh++; }

  double xlow = htest->GetBinLowEdge(binlow);
  double xhigh = htest->GetBinLowEdge(binhigh+1);

  htest->Fit("gaus","q0S","",xlow, xhigh);
}

void GetTrackingCutsFastSBS_StraightThrough( const char *configfilename, const char *outfilename="SBSGEMtrackingcuts_temp.root", int nmodules=26, double thresh=0.003, double nsig_tstrip=4.5, double nsig_dt=5.0 ){

  ifstream infile(configfilename);
  
  TChain *C = new TChain("T");

  TString currentline;
  while( infile >> currentline && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline.Data());
    }
  }

  TCut globalcut = "";
  
  while( infile >> currentline && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline.Data();
    }
  }


  C->SetBranchStatus("*",0);
  
  double EPS, ESH;

  //Unclear whether we should care about these variables for present purposes:
  C->SetBranchStatus("bb.etot_over_p",1);
  C->SetBranchStatus("bb.ps.e",1);
  C->SetBranchStatus("bb.sh.e",1);
  
  C->SetBranchStatus("bb.grinch_tdc.nclus",1);
  C->SetBranchStatus("bb.grinch_tdc.clus.*",1);


  //Set up trigger branches:
  C->SetBranchStatus("Ndata.sbs.hcal.Ref.*",1);
  C->SetBranchStatus("sbs.hcal.Ref.*",1);

  C->SetBranchStatus("sbs.hcal.*",1);
  
  C->SetBranchStatus("Ndata.bb.tdctrig.*",1);
  C->SetBranchStatus("bb.tdctrig.*",1);
 
  C->SetBranchStatus("sbs.gem.track.*",1);
  C->SetBranchStatus("sbs.tr.*",1);
  C->SetBranchStatus("sbs.gem.hit.*",1);
  C->SetBranchStatus("bb.gem.track.*",1);
  C->SetBranchStatus("bb.gem.hit.*",1);
  C->SetBranchStatus("bb.tr.*",1);

  C->SetBranchStatus("sbs.x_fcp",1);
  C->SetBranchStatus("sbs.y_fcp",1);
  C->SetBranchStatus("sbs.z_fcp",1);
  C->SetBranchStatus("sbs.x_bcp",1);
  C->SetBranchStatus("sbs.y_bcp",1);
  C->SetBranchStatus("sbs.z_bcp",1);
  
  
  C->SetBranchStatus("bb.x_fcp",1);
  C->SetBranchStatus("bb.y_fcp",1);
  C->SetBranchStatus("bb.z_fcp",1);
  C->SetBranchStatus("bb.x_bcp",1);
  C->SetBranchStatus("bb.y_bcp",1);
  C->SetBranchStatus("bb.z_bcp",1);
  
 
  TFile *fout = new TFile(outfilename,"RECREATE");
  
  //Next: initialize histograms and output file:
  TH1D *hdxfcp = new TH1D("hdxfcp", "Front Constraint X;x_{track}-x_{fcp} (m);", 500, -0.15, 0.15 );
  TH1D *hdyfcp = new TH1D("hdyfcp", "Front Constraint Y;y_{track}-y_{fcp} (m);", 500, -0.15, 0.15 );
  TH1D *hdxbcp = new TH1D("hdxbcp", "Back Constraint X;x_{track}+x'_{track}z_{bcp}-x_{bcp} (m);", 500, -0.15, 0.15 );
  TH1D *hdybcp = new TH1D("hdybcp", "Back Constraint Y;y_{track}+y'_{track}z_{bcp}-y_{bcp} (m);", 500, -0.15, 0.15 );

  TH2D *hdxdyfcp = new TH2D("hdxdyfcp", ";y_{track}-y_{fcp} (m);x_{track}-x_{fcp} (m)", 150, -0.15, 0.15, 150, -0.15, 0.15 );
  TH2D *hdxdybcp = new TH2D("hdxdybcp", ";y_{track}+y'_{track}z_{bcp}-y_{bcp} (m);x_{track}+x'_{track}z_{bcp}-x_{bcp} (m)", 150, -0.15, 0.15, 150, -0.15, 0.15 );

  TH1D *hdthcp = new TH1D("hdthcp", "Track #theta Constraint; x'_{track}-x'_{constraint};", 300, -0.1, 0.1 );
  TH1D *hdphcp = new TH1D("hdphcp", "Track #phi Constraint; y'_{track}-y'_{constraint};", 300, -0.1, 0.1 );

  //Forget TClonesArrays, use TH2D for module-dependent stuff:
  TH2D *hADCmaxsamp_vs_module = new TH2D("hADCmaxsamp_vs_module",";module;max strip max ADC sample", nmodules, -0.5, nmodules-0.5, 1500,0,3000);
  TH2D *hADCmaxstrip_vs_module = new TH2D("hADCmaxstrip_vs_module",";module;max strip ADC sum", nmodules, -0.5, nmodules-0.5, 1500,0,15000);

  TH2D *hADCclust_vs_module = new TH2D("hADCclust_vs_module",";module;cluster sum", nmodules, -0.5, nmodules-0.5, 1500,0,30000);


  TH2D *hdeltat_mod = new TH2D("hdeltat_mod", "Standard hit times ;module;#Delta t (ns)", nmodules, -0.5, nmodules-0.5, 200,-30,30);
  TH2D *hdeltat_deconv_mod = new TH2D("hdeltat_deconv_mod", "Deconvoluted hit times ;module;#Delta t (ns)", nmodules, -0.5, nmodules-0.5, 200,-80,80);
  TH2D *hdeltat_fit_mod = new TH2D("hdeltat_fit_mod", "Fitted hit times ;module;#Delta t (ns)", nmodules, -0.5, nmodules-0.5, 200,-30,30);

  TH2D *hADCasym_mod = new TH2D("hADCasym_mod", ";module; ADC asymmetry", nmodules,-0.5,nmodules-0.5, 200,-0.8,0.8);
  TH2D *hADCasym_deconv_mod = new TH2D("hADCasym_deconv_mod", ";module; ADC asymmetry", nmodules,-0.5,nmodules-0.5, 200,-1.0,1.0);

  TH2D *hADCratio_mod = new TH2D("hADCratio_mod", ";module; ADC V/U ratio", nmodules,-0.5,nmodules-0.5, 200,0.0,2.0);
  TH2D *hADCratio_deconv_mod = new TH2D("hADCratio_deconv_mod", ";module; ADC V/U ratio", nmodules,-0.5,nmodules-0.5, 200,0.0,2.0);

  TH2D *hmaxstrip_tU_mod = new TH2D("hmaxstrip_tU_mod", "Max U strip ;module; time (ns)", nmodules, -0.5, nmodules-0.5, 150,0,150); 
  TH2D *hmaxstrip_tV_mod = new TH2D("hmaxstrip_tV_mod", "Max V strip ;module; time (ns)", nmodules, -0.5, nmodules-0.5, 150,0,150); 

  TH2D *hmaxstrip_tU_deconv_mod = new TH2D("hmaxstrip_tU_deconv_mod", "Max U strip ;module; deconv time (ns)", nmodules, -0.5, nmodules-0.5, 150,-100,200); 
  TH2D *hmaxstrip_tV_deconv_mod = new TH2D("hmaxstrip_tV_deconv_mod", "Max V strip ;module; deconv time (ns)", nmodules, -0.5, nmodules-0.5, 150,-100,200); 

  TH2D *hmaxstrip_tU_fit_mod = new TH2D("hmaxstrip_tU_fit_mod", "Max U strip ;module; fit time (ns)", nmodules, -0.5, nmodules-0.5, 150,-100,200); 
  TH2D *hmaxstrip_tV_fit_mod = new TH2D("hmaxstrip_tV_fit_mod", "Max V strip ;module; fit time (ns)", nmodules, -0.5, nmodules-0.5, 150,-100,200); 

  TH2D *htU_mod = new TH2D("htU_mod", "U clusters;module; time (ns)", nmodules,-0.5,nmodules-0.5, 150,0,150);
  TH2D *htV_mod = new TH2D("htV_mod", "V clusters;module; time (ns)", nmodules,-0.5,nmodules-0.5, 150,0,150);

  TH2D *htU_deconv_mod = new TH2D("htU_deconv_mod", "U clusters;module; time deconv (ns)", nmodules,-0.5,nmodules-0.5, 150,-100,200);
  TH2D *htV_deconv_mod = new TH2D("htV_deconv_mod", "V clusters;module; time deconv (ns)", nmodules,-0.5,nmodules-0.5, 150,-100,200);
  
  TH2D *htU_fit_mod = new TH2D("htU_fit_mod", "U clusters;module; fit time (ns)", nmodules,-0.5,nmodules-0.5, 150,-100,200);
  TH2D *htV_fit_mod = new TH2D("htV_fit_mod", "V clusters;module; fit time (ns)", nmodules,-0.5,nmodules-0.5, 150,-100,200);

  TH2D *hdt_corr_mod = new TH2D("hdt_corr_mod", "Hit average corrected time minus track t0;T_{corr}-t_{0}^{track} (ns);", nmodules,-0.5,nmodules-0.5, 100,-50,50);

  TH2D *htavg_corr_vs_ttrig_allhits = new TH2D("hTavg_corr_vs_ttrig_allhits", "All hits; t_{trig} (ns); t_{GEM} (ns)", 300, 200, 500, 300,-75,75);

  TH1D *htrackt0 = new TH1D("htrackt0", ";track t_{0} (ns);", 250,-25,25);
  
  // For the GEM time versus trigger time correlation, we want a TClonesArray(TH2D):
  TClonesArray *htavg_corr_vs_ttrig_by_module = new TClonesArray("TH2D",nmodules);
 

  for( int imod=0; imod<nmodules; imod++ ){
    TString hname,htitle;
    new( (*htavg_corr_vs_ttrig_by_module)[imod] ) TH2D( hname.Format("hTavg_corr_vs_ttrig_mod%d",imod), htitle.Format("Module %d;t_{trig} (ns);t_{GEM} (ns)",imod), 300,200,500, 300,-75,75); 
  }


  
  TString fname_db = outfilename;
  fname_db.ReplaceAll(".root",".dat");
  
  ofstream dbfile(fname_db.Data());

  fname_db.ReplaceAll(".dat", "_timecuts.dat");
  ofstream dbfile_tcuts(fname_db.Data() );

  long nevent=0;
  
  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut, C );

  genrp_tree_SBS_straight_through *T = new genrp_tree_SBS_straight_through(C);
  
  int treenum=0, currenttreenum=0;

  while( C->GetEntry( nevent++ ) && nevent){
    currenttreenum = C->GetTreeNumber();
    if( nevent == 1 || currenttreenum != treenum ){
      std::cout << "new file, name = " << C->GetFile()->GetName() << ", nevent = " << nevent << std::endl;
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();
    }

    if( nevent % 1000 == 0 ) cout << nevent << endl;
    
    bool passedcut = GlobalCut->EvalInstance(0) != 0;

    // cout << "passed cut, ntracks, ngoodhits = " << passedcut << ", "
    // 	 << ntracks << ", " << ngoodhits << endl;


    
    if( passedcut && T->sbs_gem_track_ntrack >= 1){

      //Grab trigger TDC time: 
      double ttrig = 0.0;
      for( int itrig=0; itrig<T->Ndata_bb_tdctrig_tdc; itrig++ ){
	if( int(T->bb_tdctrig_tdcelemID[itrig]) == 5 ){
	  ttrig = T->bb_tdctrig_tdcelemID[itrig];
	}
      }

      double xfcp = T->sbs_x_fcp[0];
      double yfcp = T->sbs_y_fcp[0];
      double zfcp = T->sbs_z_fcp[0];

      double xbcp = T->sbs_x_bcp[0];
      double ybcp = T->sbs_y_bcp[0];
      double zbcp = T->sbs_z_bcp[0];
      
      
      double thcp = (T->sbs_x_bcp[0]-T->sbs_x_fcp[0])/(T->sbs_z_bcp[0]-T->sbs_z_fcp[0]);
      double phcp = (T->sbs_y_bcp[0]-T->sbs_y_fcp[0])/(T->sbs_z_bcp[0]-T->sbs_z_fcp[0]);
      //if( EPS >= 0.2 ){ //tracking constraint histos always require preshower regardless what the user defined?
      hdxfcp->Fill( T->sbs_tr_x[0]-xfcp );
      hdyfcp->Fill( T->sbs_tr_y[0]-yfcp );
      hdxbcp->Fill( T->sbs_tr_x[0]+T->sbs_tr_th[0]*(zbcp-zfcp)-xbcp );
      hdybcp->Fill( T->sbs_tr_y[0]+T->sbs_tr_ph[0]*(zbcp-zfcp)-ybcp );
	
      hdxdyfcp->Fill( T->sbs_tr_y[0]-yfcp, T->sbs_tr_x[0] - xfcp );
      hdxdybcp->Fill( T->sbs_tr_y[0]+T->sbs_tr_ph[0]*(zbcp-zfcp) - ybcp,
		      T->sbs_tr_x[0]+T->sbs_tr_th[0]*(zbcp-zfcp) - xbcp );
       	
      hdthcp->Fill( T->sbs_tr_th[0] - thcp );
      hdphcp->Fill( T->sbs_tr_ph[0] - phcp );

      htrackt0->Fill( T->sbs_gem_track_t0[0] );
      
      
      for( int ihit=0; ihit<int(T->sbs_gem_track_nhits[0]); ihit++ ){
	if( T->sbs_gem_hit_trackindex[ihit] == 0 ){
	  htavg_corr_vs_ttrig_allhits->Fill( ttrig, T->sbs_gem_hit_Tavg_corr[ihit] );

	  ( (TH2D *) (*htavg_corr_vs_ttrig_by_module)[T->sbs_gem_hit_module[ihit]] )->Fill( ttrig, T->sbs_gem_hit_Tavg_corr[ihit] );
	  
	  hADCmaxsamp_vs_module->Fill( T->sbs_gem_hit_module[ihit], 0.5*(T->sbs_gem_hit_ADCmaxsampU[ihit]+T->sbs_gem_hit_ADCmaxsampV[ihit]) );
	  hADCmaxstrip_vs_module->Fill( T->sbs_gem_hit_module[ihit], 0.5*(T->sbs_gem_hit_ADCmaxstripU[ihit]+T->sbs_gem_hit_ADCmaxstripV[ihit]) );

	  hADCclust_vs_module->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_ADCavg[ihit] );

	  //For timing cuts we require hits to pass higher ADC thresholds:
	  if( T->sbs_gem_hit_ADCavg[ihit] >= 1500.0 ){
	    hADCasym_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_ADCasym[ihit] );
	    hADCasym_deconv_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_ADCasym_deconv[ihit] );
	    hADCratio_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_ADCV[ihit]/T->sbs_gem_hit_ADCU[ihit] );
	    hADCratio_deconv_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_ADCV_deconv[ihit]/T->sbs_gem_hit_ADCU_deconv[ihit] );
	    
	    hdeltat_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_deltat[ihit] );
	    hdeltat_deconv_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_deltat_deconv[ihit] );
	    hdeltat_fit_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_deltat_fit[ihit] );
	  }

	  if( T->sbs_gem_hit_ADCmaxsampU[ihit] >= 100. && T->sbs_gem_hit_ADCmaxsampV[ihit] >= 100. &&
	      T->sbs_gem_hit_ADCmaxstripU[ihit] >= 250. && T->sbs_gem_hit_ADCmaxstripV[ihit] >= 250. &&
	      T->sbs_gem_hit_ADCU[ihit] >= 600. && T->sbs_gem_hit_ADCV[ihit] >= 600. &&
	      T->sbs_gem_hit_ADCavg[ihit] >= 750. ){
	    hmaxstrip_tU_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_UtimeMaxStrip[ihit] );
	    hmaxstrip_tV_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_VtimeMaxStrip[ihit] );
	    
	    hmaxstrip_tU_deconv_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_UtimeMaxStripDeconv[ihit] );
	    hmaxstrip_tV_deconv_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_VtimeMaxStripDeconv[ihit] );

	    hmaxstrip_tU_fit_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_UtimeMaxStripFit[ihit] );
	    hmaxstrip_tV_fit_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_VtimeMaxStripFit[ihit] );
	    
	    htU_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_Utime[ihit] );
	    htV_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_Vtime[ihit] );
	    
	    htU_deconv_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_UtimeDeconv[ihit] );
	    htV_deconv_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_VtimeDeconv[ihit] );
	    
	    htU_fit_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_UtimeFit[ihit] );
	    htV_fit_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_VtimeFit[ihit] );

	    hdt_corr_mod->Fill( T->sbs_gem_hit_module[ihit], T->sbs_gem_hit_Tavg_corr[ihit]-T->sbs_gem_track_t0[0] );
	  }
	}
      }
    }
  }

  //// Here we create a pdf which will be filled in the following loop ////

  TString outfilepdf = outfilename;
  outfilepdf.ReplaceAll(".root",".pdf");

  TString outfilepdf_time = outfilename;
  outfilepdf_time.ReplaceAll(".root","_timecuts.pdf");

  TCanvas *c1 = new TCanvas("c1","",1600,1200);
  c1->Divide(2,3);
  c1->cd(1);
  DrawHist(hdxfcp,{hdxfcp->GetMean(),hdxfcp->GetRMS()*4.5});
  c1->cd(2);
  DrawHist(hdyfcp,{hdyfcp->GetMean(),hdyfcp->GetRMS()*4.5});
  c1->cd(3);
  DrawHist(hdxbcp,{hdxbcp->GetMean(),hdxbcp->GetRMS()*4.5});
  c1->cd(4);
  DrawHist(hdybcp,{hdybcp->GetMean(),hdybcp->GetRMS()*4.5});
  c1->cd(5);
  DrawHist(hdthcp,{hdthcp->GetMean(),hdthcp->GetRMS()*4.5});
  c1->cd(6);
  DrawHist(hdphcp,{hdphcp->GetMean(),hdphcp->GetRMS()*4.5});

  c1->Print(outfilepdf + "(");  //Open the pdf and make the first page

  ///////////////////////////////////////////////////////////////////
  double nsigma = nsig_tstrip;
  
  double maxstrip_t0[nmodules][2];
  double maxstrip_tsigma[nmodules][2];
  //double maxstrip_tcut[nmodules][2];

  double maxstrip_t0_deconv[nmodules][2];
  double maxstrip_tsigma_deconv[nmodules][2];
  //double maxstrip_tcut_deconv[nmodules][2];

  double maxstrip_t0_fit[nmodules][2];
  double maxstrip_tsigma_fit[nmodules][2];

  double tmean[nmodules][6];
  double tsigma[nmodules][6];
 
  //double maxstrip_tcut_fit[nmodules][2];
  
  //Let's start with timing cuts:
  for( int imod=0; imod<nmodules; imod++ ){

    TH1D *hADCtemp = hADCmaxsamp_vs_module->ProjectionY("hADCtemp",imod+1,imod+1);
    ///Create canvas which will hold the ADC plots //////
    TCanvas *c= new TCanvas("c","",1200,1200);
    c->Divide(1,3);

    if( hADCtemp->GetEntries() >= 300 ){

      int binlow=1;

      while( hADCtemp->Integral(1,binlow) < thresh*hADCtemp->GetEntries() ){binlow++;}

      vector<double> threshsamp;
      threshsamp.push_back(hADCtemp->GetBinLowEdge(binlow));
      
      TString dbline;
      dbfile << dbline.Format("sbs.gem.m%d.threshold_sample = %12.5g", imod, threshsamp[0] ) << endl;

      ///// Add this histogram to the pdf output
      TH1D *htemp = hADCtemp;
      htemp->SetName(Form("hsample_m%i",imod));
      htemp->SetTitle(Form("Module %i Sample ADC",imod));
      c->cd(1);
      DrawHist(htemp,threshsamp,1);
    }

    hADCtemp = hADCmaxstrip_vs_module->ProjectionY("hADCtemp",imod+1,imod+1);

    if( hADCtemp->GetEntries() >= 300 ){
      int binlow=1;
      
      while( hADCtemp->Integral(1,binlow) < thresh*hADCtemp->GetEntries() ){binlow++;}
      
      vector<double> threshstrip;
      threshstrip.push_back(hADCtemp->GetBinLowEdge(binlow));
      
      TString dbline;
      dbfile << dbline.Format("sbs.gem.m%d.threshold_stripsum = %12.5g", imod, threshstrip[0] ) << endl;

      ///// Add this histogram to the pdf output
      TH1D *htemp = hADCtemp;
      htemp->SetName(Form("hstrip_m%i",imod));
      htemp->SetTitle(Form("Module %i Strip Sum ADC",imod));
      c->cd(2);
      DrawHist(htemp,threshstrip,1);
    }

    hADCtemp = hADCclust_vs_module->ProjectionY("hADCtemp",imod+1,imod+1);

    if( hADCtemp->GetEntries() >= 300 ){
      int binlow=1;
      
      while( hADCtemp->Integral(1,binlow) < thresh*hADCtemp->GetEntries() ){binlow++;}
      
      vector<double> thresh;
      thresh.push_back(hADCtemp->GetBinLowEdge(binlow));
      
      TString dbline;
      dbfile << dbline.Format("sbs.gem.m%d.threshold_clustersum = %12.5g", imod, thresh[0] ) << endl << endl;

      ///// Add this histogram to the pdf output
      TH1D *htemp = hADCtemp;
      htemp->SetName(Form("hclust_m%i",imod));
      htemp->SetTitle(Form("Module %i Cluster Sum ADC",imod));
      c->cd(3);
      DrawHist(htemp,thresh,1);
    }
    

    TCanvas *c2= new TCanvas("c2","",1600,1200);
    c2->Divide(2,3);

    TString hname;
    TH1D *htempU, *htempV;

    htempU = hmaxstrip_tU_mod->ProjectionY( "htempU", imod+1, imod+1 );
    htempV = hmaxstrip_tV_mod->ProjectionY( "htempV", imod+1, imod+1 );
    
    if( htempU->GetEntries() >= 300 ){
      
      FitGaus_FWHM( htempU, 0.3 );
      FitGaus_FWHM( htempV, 0.3 );
      
      TF1 *fitfuncU = (TF1*) (htempU->GetListOfFunctions()->FindObject("gaus"));
      TF1 *fitfuncV = (TF1*) (htempV->GetListOfFunctions()->FindObject("gaus"));
      
      maxstrip_t0[imod][0] = fitfuncU->GetParameter("Mean");
      maxstrip_t0[imod][1] = fitfuncV->GetParameter("Mean");
      maxstrip_tsigma[imod][0] = fitfuncU->GetParameter("Sigma");
      maxstrip_tsigma[imod][1] = fitfuncV->GetParameter("Sigma");

      TString dbline;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.maxstrip_t0 = %10.4g %10.4g",imod,maxstrip_t0[imod][0], maxstrip_t0[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.maxstrip_tsigma = %10.4g %10.4g",imod,maxstrip_tsigma[imod][0], maxstrip_tsigma[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.maxstrip_tcut = %10.4g %10.4g", imod, nsigma, nsigma ) << endl << endl;
      
      ///// Now we create the histograms for the pdf results//////
      vector<double> resultsU, resultsV;
      resultsU.push_back(maxstrip_t0[imod][0]);
      resultsU.push_back(maxstrip_tsigma[imod][0]);
      resultsV.push_back(maxstrip_t0[imod][1]);
      resultsV.push_back(maxstrip_tsigma[imod][1]);

      TH1D *htU = htempU;
      htU->SetName(Form("htmaxstripU_m%i",imod));
      htU->SetTitle(Form("Module %i U Time;time (ns)",imod));
      c2->cd(1);
      DrawHist(htU,resultsU,2);
      fitfuncU->Draw("same");

      TH1D *htV = htempV;
      htV->SetName(Form("htmaxstripV_m%i",imod));
      htV->SetTitle(Form("Module %i V Time;time (ns)",imod));
      c2->cd(2);
      DrawHist(htV,resultsV,2);
      fitfuncV->Draw("same");
    }

    htempU = hmaxstrip_tU_deconv_mod->ProjectionY("htempU", imod+1,imod+1);
    htempV = hmaxstrip_tV_deconv_mod->ProjectionY("htempV", imod+1,imod+1);
    
    if( htempU->GetEntries() >= 300 ){
      
      FitGaus_FWHM( htempU, 0.3 );
      FitGaus_FWHM( htempV, 0.3 );
      
      TF1 *fitfuncU = (TF1*) (htempU->GetListOfFunctions()->FindObject("gaus"));
      TF1 *fitfuncV = (TF1*) (htempV->GetListOfFunctions()->FindObject("gaus"));
      
      maxstrip_t0_deconv[imod][0] = fitfuncU->GetParameter("Mean");
      maxstrip_t0_deconv[imod][1] = fitfuncV->GetParameter("Mean");
      maxstrip_tsigma_deconv[imod][0] = fitfuncU->GetParameter("Sigma");
      maxstrip_tsigma_deconv[imod][1] = fitfuncV->GetParameter("Sigma");

      TString dbline;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.maxstrip_t0_deconv = %10.4g %10.4g",imod,maxstrip_t0_deconv[imod][0], maxstrip_t0_deconv[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.maxstrip_tsigma_deconv = %10.4g %10.4g",imod,maxstrip_tsigma_deconv[imod][0], maxstrip_tsigma_deconv[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.maxstrip_tcut_deconv = %10.4g %10.4g", imod, nsigma, nsigma ) << endl << endl;

      ///// Now we create the histograms for the pdf results//////
      vector<double> resultsU, resultsV;
      resultsU.push_back(maxstrip_t0_deconv[imod][0]);
      resultsU.push_back(maxstrip_tsigma_deconv[imod][0]);
      resultsV.push_back(maxstrip_t0_deconv[imod][1]);
      resultsV.push_back(maxstrip_tsigma_deconv[imod][1]);
      
      TH1D *htU = htempU;
      htU->SetName(Form("htstripdeconU_m%i",imod));
      htU->SetTitle(Form("Module %i U Deconvolution;deconv time (ns)",imod));
      c2->cd(3);
      DrawHist(htU,resultsU,2);
      fitfuncU->Draw("same");

      TH1D *htV = htempV;
      htV->SetName(Form("htstripdeconV_m%i",imod));
      htV->SetTitle(Form("Module %i V Deconvolution;deconv time (ns)",imod));
      c2->cd(4);
      DrawHist(htV,resultsV,2);
      fitfuncV->Draw("same");
    }
 
    htempU = hmaxstrip_tU_fit_mod->ProjectionY("htempU", imod+1,imod+1);
    htempV = hmaxstrip_tV_fit_mod->ProjectionY("htempV", imod+1,imod+1);
    
    if( htempU->GetEntries() >= 300 ){
      
      FitGaus_FWHM( htempU, 0.3 );
      FitGaus_FWHM( htempV, 0.3 );
      
      TF1 *fitfuncU = (TF1*) (htempU->GetListOfFunctions()->FindObject("gaus"));
      TF1 *fitfuncV = (TF1*) (htempV->GetListOfFunctions()->FindObject("gaus"));
      
      maxstrip_t0_fit[imod][0] = fitfuncU->GetParameter("Mean");
      maxstrip_t0_fit[imod][1] = fitfuncV->GetParameter("Mean");
      maxstrip_tsigma_fit[imod][0] = fitfuncU->GetParameter("Sigma");
      maxstrip_tsigma_fit[imod][1] = fitfuncV->GetParameter("Sigma");

      TString dbline;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.maxstrip_t0_fit = %10.4g %10.4g",imod,maxstrip_t0_fit[imod][0], maxstrip_t0_fit[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.maxstrip_tsigma_fit = %10.4g %10.4g",imod,maxstrip_tsigma_fit[imod][0], maxstrip_tsigma_fit[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.maxstrip_tcut_fit = %10.4g %10.4g", imod, nsigma, nsigma ) << endl << endl;


      ///// Now we create the histograms for the pdf results//////
      vector<double> resultsU, resultsV;
      resultsU.push_back(maxstrip_t0_fit[imod][0]);
      resultsU.push_back(maxstrip_tsigma_fit[imod][0]);
      resultsV.push_back(maxstrip_t0_fit[imod][1]);
      resultsV.push_back(maxstrip_tsigma_fit[imod][1]);
      
      TH1D *htU = htempU;
      htU->SetName(Form("htstripfitU_m%i",imod));
      htU->SetTitle(Form("Module %i U Fit Time;deconv time (ns)",imod));
      c2->cd(5);
      DrawHist(htU,resultsU,2);
      fitfuncU->Draw("same");

      TH1D *htV = htempV;
      htV->SetName(Form("htstripfitV_m%i",imod));
      htV->SetTitle(Form("Module %i V Fit Time;deconv time (ns)",imod));
      c2->cd(6);
      DrawHist(htV,resultsV,2);
      fitfuncV->Draw("same");
    }

    TCanvas *c3= new TCanvas("c3","",1200,800);
    c3->Divide(2,2);


    TH1D *htemp;

    htemp = hdeltat_mod->ProjectionY("htemp", imod+1,imod+1 );

    if( htemp->GetEntries() >= 300 ){
      FitGaus_FWHM( htemp, 0.3 );

      TF1 *fitfunc = (TF1*) htemp->GetListOfFunctions()->FindObject("gaus");

      double mean = fitfunc->GetParameter("Mean");
      double sigma = fitfunc->GetParameter("Sigma");

      TString dbline;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.deltat_sigma = %12.5g", imod, sigma ) << endl;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.deltat_cut = %12.5g", imod, nsig_dt*sigma ) << endl;

      vector<double> results = {nsig_dt*sigma};

      TH1D *hplot = htemp;
      hplot->SetName(Form("hdeltat_m%i",imod));
      hplot->SetTitle(Form("Module %i #Deltat; #Deltat (ns)",imod));
      c3->cd(1);
      DrawHist(hplot,results,2);
      fitfunc->Draw("same");
    }

    htemp = hdeltat_deconv_mod->ProjectionY("htemp", imod+1,imod+1 );

    if( htemp->GetEntries() >= 300 ){
      FitGaus_FWHM( htemp, 0.3 );

      TF1 *fitfunc = (TF1*) htemp->GetListOfFunctions()->FindObject("gaus");

      double mean = fitfunc->GetParameter("Mean");
      double sigma = fitfunc->GetParameter("Sigma");

      TString dbline;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.deltat_sigma_deconv = %12.5g", imod, sigma ) << endl;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.deltat_cut_deconv = %12.5g", imod, nsig_dt*sigma ) << endl;

      vector<double> results = {nsig_dt*sigma};

      TH1D *hplot = htemp;
      hplot->SetName(Form("hdeltatdeconv_m%i",imod));
      hplot->SetTitle(Form("Module %i Deconv #Deltat;deconv #Deltat(ns)",imod));
      c3->cd(2);
      DrawHist(hplot,results,2);
      fitfunc->Draw("same");
    }

    htemp = hdeltat_fit_mod->ProjectionY("htemp", imod+1,imod+1 );

    if( htemp->GetEntries() >= 300 ){
      FitGaus_FWHM( htemp, 0.3 );

      TF1 *fitfunc = (TF1*) htemp->GetListOfFunctions()->FindObject("gaus");

      double mean = fitfunc->GetParameter("Mean");
      double sigma = fitfunc->GetParameter("Sigma");

      TString dbline;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.deltat_sigma_fit = %12.5g", imod, sigma ) << endl;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.deltat_cut_fit = %12.5g", imod, nsig_dt*sigma ) << endl << endl;

      vector<double> results = {nsig_dt*sigma};

      TH1D *hplot = htemp;
      hplot->SetName(Form("hdeltatfit_m%i",imod));
      hplot->SetTitle(Form("Module %i Fit #Deltat;Fit #Deltat (ns)",imod));
      c3->cd(3);
      DrawHist(hplot,results,2);
      fitfunc->Draw("same");
    }


    TCanvas *c4= new TCanvas("c4","",1600,1200);
    c4->Divide(2,3);

    htempU = htU_mod->ProjectionY("htempU",imod+1,imod+1);
    htempV = htV_mod->ProjectionY("htempV",imod+1,imod+1);

    if( htempU->GetEntries() >= 300 ){
      FitGaus_FWHM( htempU, 0.3 );
      FitGaus_FWHM( htempV, 0.3 );
      TF1 *fitfuncU = (TF1*) (htempU->GetListOfFunctions()->FindObject("gaus"));
      TF1 *fitfuncV = (TF1*) (htempV->GetListOfFunctions()->FindObject("gaus"));

      tmean[imod][0] = fitfuncU->GetParameter("Mean");
      tmean[imod][1] = fitfuncV->GetParameter("Mean");
      tsigma[imod][0] = fitfuncU->GetParameter("Sigma");
      tsigma[imod][1] = fitfuncV->GetParameter("Sigma");

      TString dbline;

      dbfile_tcuts << dbline.Format("sbs.gem.m%d.HitTimeMean = %12.5g %12.5g", imod, tmean[imod][0], tmean[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.HitTimeSigma = %12.5g %12.5g", imod, tsigma[imod][0], tsigma[imod][1] ) << endl;

      ///// Now we create the histograms for the pdf results//////
      vector<double> resultsU, resultsV;
      resultsU.push_back(tmean[imod][0]);
      resultsU.push_back(tsigma[imod][0]);
      resultsV.push_back(tmean[imod][1]);
      resultsV.push_back(tsigma[imod][1]);

      TH1D *htU = htempU;
      htU->SetName(Form("htclustU_m%i",imod));
      htU->SetTitle(Form("Module %i U Cluster Time;deconv time (ns)",imod));
      c4->cd(1);
      DrawHist(htU,resultsU,2);
      fitfuncU->Draw("same");

      TH1D *htV = htempV;
      htV->SetName(Form("htclustV_m%i",imod));
      htV->SetTitle(Form("Module %i V Cluster Time;deconv time (ns)",imod));
      c4->cd(2);
      DrawHist(htV,resultsV,2);
      fitfuncV->Draw("same");
      
    }

    htempU = htU_deconv_mod->ProjectionY("htempU",imod+1,imod+1);
    htempV = htV_deconv_mod->ProjectionY("htempV",imod+1,imod+1);

    if( htempU->GetEntries() >= 300 ){
      FitGaus_FWHM( htempU, 0.3 );
      FitGaus_FWHM( htempV, 0.3 );
      TF1 *fitfuncU = (TF1*) (htempU->GetListOfFunctions()->FindObject("gaus"));
      TF1 *fitfuncV = (TF1*) (htempV->GetListOfFunctions()->FindObject("gaus"));

      tmean[imod][2] = fitfuncU->GetParameter("Mean");
      tmean[imod][3] = fitfuncV->GetParameter("Mean");
      tsigma[imod][2] = fitfuncU->GetParameter("Sigma");
      tsigma[imod][3] = fitfuncV->GetParameter("Sigma");

      TString dbline;

      dbfile_tcuts << dbline.Format("sbs.gem.m%d.HitTimeMeanDeconv = %12.5g %12.5g", imod, tmean[imod][2], tmean[imod][3] ) << endl;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.HitTimeSigmaDeconv = %12.5g %12.5g", imod, tsigma[imod][2], tsigma[imod][3] ) << endl;

      ///// Now we create the histograms for the pdf results//////
      vector<double> resultsU, resultsV;
      resultsU.push_back(tmean[imod][0]);
      resultsU.push_back(tsigma[imod][0]);
      resultsV.push_back(tmean[imod][1]);
      resultsV.push_back(tsigma[imod][1]);

      TH1D *htU = htempU;
      htU->SetName(Form("htclustdeconU_m%i",imod));
      htU->SetTitle(Form("Module %i U Cluster Deconv Time;deconv time (ns)",imod));
      c4->cd(3);
      DrawHist(htU,resultsU,2);
      fitfuncU->Draw("same");

      TH1D *htV = htempV;
      htV->SetName(Form("htclustdeconV_m%i",imod));
      htV->SetTitle(Form("Module %i V Cluster Deconv Time;deconv time (ns)",imod));
      c4->cd(4);
      DrawHist(htV,resultsV,2);
      fitfuncV->Draw("same");
      
    }

    htempU = htU_fit_mod->ProjectionY("htempU",imod+1,imod+1);
    htempV = htV_fit_mod->ProjectionY("htempV",imod+1,imod+1);

    if( htempU->GetEntries() >= 300 ){
      FitGaus_FWHM( htempU, 0.3 );
      FitGaus_FWHM( htempV, 0.3 );
      TF1 *fitfuncU = (TF1*) (htempU->GetListOfFunctions()->FindObject("gaus"));
      TF1 *fitfuncV = (TF1*) (htempV->GetListOfFunctions()->FindObject("gaus"));

      tmean[imod][4] = fitfuncU->GetParameter("Mean");
      tmean[imod][5] = fitfuncV->GetParameter("Mean");
      tsigma[imod][4] = fitfuncU->GetParameter("Sigma");
      tsigma[imod][5] = fitfuncV->GetParameter("Sigma");

      TString dbline;

      dbfile_tcuts << dbline.Format("sbs.gem.m%d.HitTimeMeanFit = %12.5g %12.5g", imod, tmean[imod][4], tmean[imod][5] ) << endl;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.HitTimeSigmaFit = %12.5g %12.5g", imod, tsigma[imod][4], tsigma[imod][5] ) << endl << endl;

      ///// Now we create the histograms for the pdf results//////
      vector<double> resultsU, resultsV;
      resultsU.push_back(tmean[imod][0]);
      resultsU.push_back(tsigma[imod][0]);
      resultsV.push_back(tmean[imod][1]);
      resultsV.push_back(tsigma[imod][1]);

      TH1D *htU = htempU;
      htU->SetName(Form("htclustfitU_m%i",imod));
      htU->SetTitle(Form("Module %i U Cluster Fit Time;deconv time (ns)",imod));
      c4->cd(5);
      DrawHist(htU,resultsU,2);
      fitfuncU->Draw("same");

      TH1D *htV = htempV;
      htV->SetName(Form("htclustfitV_m%i",imod));
      htV->SetTitle(Form("Module %i V Cluster Fit Time;deconv time (ns)",imod));
      c4->cd(6);
      DrawHist(htV,resultsV,2);
      fitfuncV->Draw("same");
      
    }

    TCanvas *c5= new TCanvas("c5","",1200,800);
    c5->Divide(2,2);

    htemp = hADCasym_mod->ProjectionY("htemp",imod+1,imod+1);

    if( htemp->GetEntries() >= 300 ){
      FitGaus_FWHM( htemp, 0.4 );
      TF1 *fitfunc = (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus"));
      double mean = fitfunc->GetParameter("Mean");
      double sigma = fitfunc->GetParameter("Sigma");

      TString dbline;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.ADCasym_sigma = %12.5g", imod, sigma ) << endl;
      
      ///// Now we create the histograms for the pdf results//////
      vector<double> results = {mean,sigma};
      
      TH1D *hplot = htemp;
      hplot->SetName(Form("hADCasym_m%i",imod));
      hplot->SetTitle(Form("Module %i ADC Asymmetry;ADC Asym",imod));
      c5->cd(1);
      DrawHist(hplot,results,2);
      fitfunc->Draw("same");
      
    }

    htemp = hADCratio_mod->ProjectionY("htemp",imod+1,imod+1);

    if( htemp->GetEntries() >= 300 ){
      FitGaus_FWHM( htemp, 0.4 );

      TF1 *fitfunc = (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus"));
      double mean = fitfunc->GetParameter("Mean");
      double sigma = fitfunc->GetParameter("Sigma");

      TString dbline;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.ADCratio_sigma = %12.5g", imod, sigma ) << endl;
      
      ///// Now we create the histograms for the pdf results//////
      vector<double> results = {mean,sigma};
      
      TH1D *hplot = htemp;
      hplot->SetName(Form("hADCratio_m%i",imod));
      hplot->SetTitle(Form("Module %i ADC Ratio;Clust ADC V/U",imod));
      c5->cd(2);
      DrawHist(hplot,results,2);
      fitfunc->Draw("same");
    }

    htemp = hADCratio_deconv_mod->ProjectionY("htemp",imod+1,imod+1);

    if( htemp->GetEntries() >= 300 ){
      FitGaus_FWHM( htemp, 0.4 );

      TF1 *fitfunc = (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus"));
      double mean = fitfunc->GetParameter("Mean");
      double sigma = fitfunc->GetParameter("Sigma");
      
      TString dbline;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.ADCratio_sigma_deconv = %12.5g", imod, sigma ) << endl << endl;

      ///// Now we create the histograms for the pdf results//////
      vector<double> results = {mean,sigma};
      
      TH1D *hplot = htemp;
      hplot->SetName(Form("hADCratiodecon_m%i",imod));
      hplot->SetTitle(Form("Module %i Deconvolution ADC Ratio;Clust V ADC / Clust U ADC",imod));
      c5->cd(3);
      DrawHist(hplot,results,2);
      fitfunc->Draw("same");
      
    }

    htemp = hdt_corr_mod->ProjectionY("htemp",imod+1,imod+1);

    if( htemp->GetEntries() >= 300 ){
      FitGaus_FWHM( htemp, 0.3 );

      TF1 *fitfunc = (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus"));
      double mean = fitfunc->GetParameter("Mean");
      double sigma = fitfunc->GetParameter("Sigma");

      TString dbline;
      dbfile_tcuts << dbline.Format("sbs.gem.m%d.sigma_tcorr = %12.5g", imod, sigma ) << endl << endl;
      
      ///// Now we create the histograms for the pdf results//////
      vector<double> results = {mean,sigma};
      
      TH1D *hplot = htemp;
      hplot->SetName(Form("hADCratiodecon_m%i",imod));
      hplot->SetTitle(Form("Module %i Hit Corrected Time Difference;t_{corr}-t_{0}^{track} (ns);",imod));
      c5->cd(4);
      DrawHist(hplot,results,2);
      fitfunc->Draw("same");
    }
    
    ///// Now we will plot the ADC threshold histograms /////
    if(imod == nmodules - 1) c->Print(outfilepdf + ")");
    else c->Print(outfilepdf);

    if(imod == 0) {
      c2->Print(outfilepdf_time + "(");
      c3->Print(outfilepdf_time);
      c4->Print(outfilepdf_time);
      c5->Print(outfilepdf_time);
    }
    else if(imod == nmodules - 1){
      c2->Print(outfilepdf_time);
      c3->Print(outfilepdf_time);
      c4->Print(outfilepdf_time);
      c5->Print(outfilepdf_time + ")");
    }
    else{
      c2->Print(outfilepdf_time);
      c3->Print(outfilepdf_time);
      c4->Print(outfilepdf_time);
      c5->Print(outfilepdf_time);
    }
    /////////////////////////////////////////////////////////////
    
  }

  FitGaus_FWHM( htrackt0, 0.25 );

  TF1 *fitfunc = (TF1*) (htrackt0->GetListOfFunctions()->FindObject("gaus"));

  double sigma = fitfunc->GetParameter("Sigma");

  TString dbline;
  dbfile_tcuts << endl << dbline.Format("sbs.gem.sigmatrackt0 = %12.5g", sigma ) << endl << endl;
  
  
  //Tracking constraints:
  dbfile << "sbs.frontconstraint_x0 = " << hdxfcp->GetMean() << endl;
  dbfile << "sbs.frontconstraint_y0 = " << hdyfcp->GetMean() << endl;
  dbfile << "sbs.backconstraint_x0 = " << hdxbcp->GetMean() << endl;
  dbfile << "sbs.backconstraint_y0 = " << hdybcp->GetMean() << endl << endl;

  dbfile << "sbs.frontconstraintwidth_x = " << hdxfcp->GetRMS() * 4.5 << endl;
  dbfile << "sbs.frontconstraintwidth_y = " << hdyfcp->GetRMS() * 4.5 << endl;
  dbfile << "sbs.backconstraintwidth_x = " << hdxbcp->GetRMS() * 4.5 << endl;
  dbfile << "sbs.backconstraintwidth_y = " << hdybcp->GetRMS() * 4.5 << endl << endl;

  dbfile << "sbs.gem.constraintwidth_theta = " << hdthcp->GetRMS() * 4.5 << endl;
  dbfile << "sbs.gem.constraintwidth_phi = " << hdphcp->GetRMS() * 4.5 << endl;
  
  //hdeltat_mod->Draw("colz");

  //Is this line necessary? NO: it leads to duplicate copies of these histograms
  //htavg_corr_vs_ttrig_by_module->Write();

  fout->Write();

  
}
