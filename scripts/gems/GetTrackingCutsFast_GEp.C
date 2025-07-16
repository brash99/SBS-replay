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

void GetTrackingCutsFast_GEp( const char *configfilename, const char *outfilename="TrackingCuts_GEp.root", int nmodules=14, double thresh=0.003, double nsig_tstrip=4.5, double nsig_dt=5.0 ){

  ifstream infile(configfilename);
  
  TChain *C = new TChain("T");
  cout << "Ok, let's get started!" << endl;
  TString currentline;


  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline);
    }
  }

  TCut globalcut = "";
  
  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline.Data();
    }
  }

  //Set parameters for detector type
//  int nmodules=14;
  TString detname = "sbs.gemFT";

  //Parses onfigfile and update detector type
  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endconfig") ){
    if( currentline.BeginsWith("detname") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      if( tokens->GetEntries() >= 2 ){
        detname = ( (TObjString*) (*tokens)[1] )->GetString();
      }
    }
    if( currentline.BeginsWith("nmodules") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      if( tokens->GetEntries() >= 2 ){
        nmodules = ( (TObjString*) (*tokens)[1] )->GetString().Atoi();
       // nstripx_mod.resize(nmodules);
       // nstripy_mod.resize(nmodules);
      }
    }
  }


  cout << "number of modules = " << nmodules << endl;
  //Initialize branches:
  
  //what branches do we need? Tracking stuff:
  //int ntracks=0;
  const int MAXNTRACKS=100;
  //int nhits=0;
  const int MAXNHITS = 1000;
  const int MAXNCLUST = 1;
  const int MAXNTRIG = 10;
  
// List of variables to cut on
  double ntracks;
  double tracknhits[MAXNTRACKS];
  double trackngoodhits[MAXNTRACKS];
  double trackchi2ndf[MAXNTRACKS];
  double trackchi2ndf_hitquality[MAXNTRACKS];
  double trackt0[MAXNTRACKS];
  double xfp[MAXNTRACKS], yfp[MAXNTRACKS], thfp[MAXNTRACKS], phfp[MAXNTRACKS];
  double xtar[MAXNTRACKS], ytar[MAXNTRACKS], thtar[MAXNTRACKS], phtar[MAXNTRACKS];
  double rxfp[MAXNTRACKS], ryfp[MAXNTRACKS], rthfp[MAXNTRACKS], rphfp[MAXNTRACKS];
  double dxfp[MAXNTRACKS], dyfp[MAXNTRACKS], dthfp[MAXNTRACKS], dphfp[MAXNTRACKS];
  double p[MAXNTRACKS],px[MAXNTRACKS],py[MAXNTRACKS],pz[MAXNTRACKS],vz[MAXNTRACKS]; 
  double xfcp[MAXNTRACKS],yfcp[MAXNTRACKS],zfcp[MAXNTRACKS],xbcp[MAXNTRACKS],ybcp[MAXNTRACKS],zbcp[MAXNTRACKS];

  int Ntrig;
  double trig_elemID[MAXNTRIG];
  double trig_tdc[MAXNTRIG];

  // Get all these variables in a vector to easily enable their branches
  TString branchname;
  map<TString,TString> branchnames;
  vector<TString> varnames;
  varnames.push_back("hit.ngoodhits");
  varnames.push_back("hit.trackindex");
  varnames.push_back("hit.module");
  varnames.push_back("hit.layer");
  varnames.push_back("hit.nstripu");
  varnames.push_back("hit.nstripv");
  varnames.push_back("hit.ustripmax");
  varnames.push_back("hit.ustriplo");
  varnames.push_back("hit.ustriphi");
  varnames.push_back("hit.vstripmax");
  varnames.push_back("hit.vstriplo");
  varnames.push_back("hit.vstriphi");
  varnames.push_back("hit.ADCU");
  varnames.push_back("hit.ADCV");
  varnames.push_back("hit.Ugain");
  varnames.push_back("hit.Vgain");
  varnames.push_back("hit.ADCavg");
  varnames.push_back("hit.ADCasym");
 
  varnames.push_back("hit.ADCmaxsampU");
  varnames.push_back("hit.ADCmaxsampV");
  varnames.push_back("hit.ADCmaxstripU");
  varnames.push_back("hit.ADCmaxstripV");
 
  varnames.push_back("hit.ccor_clust");
  varnames.push_back("hit.ccor_strip");
  varnames.push_back("hit.deltat");
 
  varnames.push_back("hit.Tavg_corr");
  varnames.push_back("hit.UtimeMaxStrip");
 
  varnames.push_back("hit.UtimeMaxStrip");
  varnames.push_back("hit.VtimeMaxStrip");
  varnames.push_back("hit.UtimeMaxStripDeconv");
  varnames.push_back("hit.VtimeMaxStripDeconv");
  varnames.push_back("hit.UtimeMaxStripFit");
  varnames.push_back("hit.VtimeMaxStripFit");
 
  varnames.push_back("hit.Utime");
  varnames.push_back("hit.Vtime");
  varnames.push_back("hit.UtimeDeconv");
  varnames.push_back("hit.VtimeDeconv");
  varnames.push_back("hit.UtimeFit");
  varnames.push_back("hit.VtimeFit");
 
  varnames.push_back("hit.deltat_deconv");
  varnames.push_back("hit.deltat_fit");
  varnames.push_back("hit.ADCU_deconv");
  varnames.push_back("hit.ADCV_deconv");
  varnames.push_back("hit.ADCasym_deconv");
  
  varnames.push_back("track.chi2ndf_hitquality");
  varnames.push_back("track.t0");
  varnames.push_back("track.ngoodhits");
  varnames.push_back("track.ntrack");
  varnames.push_back("track.nhits");
  varnames.push_back("track.besttrack");
  varnames.push_back("track.chi2ndf");
  

  //Why are the branches disabled here? To make it run FASTER by only activating the ones you need!
  cout << "disabling all branches...";
  //the * applies it to all branches, the 0 disables those branches. to enable would need to make 1
  C->SetBranchStatus("*",0);

  cout << "done." << endl;

  for( int i=0; i<varnames.size(); i++ ){
    //What is actually going on in this for loop. It looks like branch and variable names are getting printed. But also setting branches up in the TChain
    branchnames[varnames[i]] = branchname.Format("%s.%s",detname.Data(),varnames[i].Data());
//    cout << "Branch " << i << " name = " << branchnames[varnames[i]] << endl;
    C->SetBranchStatus( branchnames[varnames[i]].Data(), 1 );
    if (varnames[i].IsNull() || varnames[i].Length() ==0) {
      cout << "error in varnames i = " << i << endl;
    }
  }
  //Populating data in the TChain branches
  cout << "Setting branch addresses: " << endl;

  //Set up trigger branches:
 // C->SetBranchAddress("Ndata.bb.tdctrig.tdc",&Ntrig);
 // C->SetBranchAddress("sbs.tdctrig.tdcelemID",trig_elemID);
 // C->SetBranchAddress("sbs.tdctrig.tdc",trig_tdc);

  C->SetBranchStatus("heep.dxECAL",1);
  C->SetBranchStatus("heep.dyECAL",1);
  C->SetBranchStatus("heep.dpp",1);
  C->SetBranchStatus("heep.dt_ADC",1);

  C->SetBranchStatus("sbs.tr.x",1);
  C->SetBranchStatus("sbs.tr.y",1);
  C->SetBranchStatus("sbs.tr.th",1);
  C->SetBranchStatus("sbs.tr.ph",1);
///  C->SetBranchStatus("sbs.tr.tg_x",1);
  C->SetBranchStatus("sbs.tr.tg_y",1);
  C->SetBranchStatus("sbs.tr.tg_th",1);
  C->SetBranchStatus("sbs.tr.tg_ph",1);
  C->SetBranchStatus("sbs.tr.r_x",1);
  C->SetBranchStatus("sbs.tr.r_y",1);
  C->SetBranchStatus("sbs.tr.r_th",1);
  C->SetBranchStatus("sbs.tr.r_ph",1);
  C->SetBranchStatus("sbs.tr.d_x",1);
  C->SetBranchStatus("sbs.tr.d_y",1);
  C->SetBranchStatus("sbs.tr.d_th",1);
  C->SetBranchStatus("sbs.tr.d_ph",1);
  C->SetBranchStatus("sbs.tr.n",1);
  C->SetBranchStatus("sbs.tr.p",1);
  C->SetBranchStatus("sbs.tr.px",1);
  C->SetBranchStatus("sbs.tr.py",1);
  C->SetBranchStatus("sbs.tr.pz",1);
  C->SetBranchStatus("sbs.tr.vz",1);
  C->SetBranchStatus("sbs.x_fcp",1);
  C->SetBranchStatus("sbs.y_fcp",1);
  C->SetBranchStatus("sbs.z_fcp",1);
  C->SetBranchStatus("sbs.x_bcp",1);
  C->SetBranchStatus("sbs.y_bcp",1);
  C->SetBranchStatus("sbs.z_bcp",1);
  
  C->SetBranchAddress( branchnames["track.nhits"].Data(),tracknhits);
  C->SetBranchAddress( branchnames["track.ngoodhits"].Data(),trackngoodhits);
  C->SetBranchAddress( branchnames["track.chi2ndf"].Data(),trackchi2ndf);
  C->SetBranchAddress( branchnames["track.chi2ndf_hitquality"].Data(),trackchi2ndf_hitquality);
  C->SetBranchAddress( branchnames["track.t0"].Data(), trackt0);
  C->SetBranchAddress("sbs.tr.x",xfp);
  C->SetBranchAddress("sbs.tr.y",yfp);
  C->SetBranchAddress("sbs.tr.th",thfp);
  C->SetBranchAddress("sbs.tr.ph",phfp);
///  C->SetBranchAddress("sbs.tr.tg_x",xtar);
  C->SetBranchAddress("sbs.tr.tg_y",ytar);
  C->SetBranchAddress("sbs.tr.tg_th",thtar);
  C->SetBranchAddress("sbs.tr.tg_ph",phtar);
  C->SetBranchAddress("sbs.tr.r_x",rxfp);
  C->SetBranchAddress("sbs.tr.r_y",ryfp);
  C->SetBranchAddress("sbs.tr.r_th",rthfp);
  C->SetBranchAddress("sbs.tr.r_ph",rphfp);
  C->SetBranchAddress("sbs.tr.d_x",dxfp);
  C->SetBranchAddress("sbs.tr.d_y",dyfp);
  C->SetBranchAddress("sbs.tr.d_th",dthfp);
  C->SetBranchAddress("sbs.tr.d_ph",dphfp);
  C->SetBranchAddress("sbs.tr.n",&ntracks); // **
  C->SetBranchAddress("sbs.tr.p",p);
  C->SetBranchAddress("sbs.tr.px",px);
  C->SetBranchAddress("sbs.tr.py",py);
  C->SetBranchAddress("sbs.tr.pz",pz);
  C->SetBranchAddress("sbs.tr.vz",vz);
  C->SetBranchAddress("sbs.x_fcp",xfcp);
  C->SetBranchAddress("sbs.y_fcp",yfcp);
  C->SetBranchAddress("sbs.z_fcp",zfcp);
  C->SetBranchAddress("sbs.x_bcp",xbcp);
  C->SetBranchAddress("sbs.y_bcp",ybcp);
  C->SetBranchAddress("sbs.z_bcp",zbcp);

  double ngoodhits;
  //Track hit variables:
  double module[MAXNHITS];
  double nstripu[MAXNHITS], nstripv[MAXNHITS];
  double trackindex[MAXNHITS];

  C->SetBranchAddress( branchnames["hit.ngoodhits"].Data(), &ngoodhits); //**
  C->SetBranchAddress( branchnames["hit.trackindex"].Data(), trackindex);
  C->SetBranchAddress( branchnames["hit.module"].Data(), module);
  C->SetBranchAddress( branchnames["hit.nstripu"].Data(), nstripu);
  C->SetBranchAddress( branchnames["hit.nstripv"].Data(), nstripv);
  
  double ADCmaxsampU[MAXNHITS], ADCmaxsampV[MAXNHITS], ADCmaxstripU[MAXNHITS], ADCmaxstripV[MAXNHITS], ADCU[MAXNHITS], ADCV[MAXNHITS], ADCavg[MAXNHITS], ADCasym[MAXNHITS], ADCasymDeconv[MAXNHITS],DeconvADCU[MAXNHITS], DeconvADCV[MAXNHITS];
 
  C->SetBranchAddress( branchnames["hit.ADCmaxsampU"].Data(), ADCmaxsampU);
  C->SetBranchAddress( branchnames["hit.ADCmaxsampV"].Data(), ADCmaxsampV);
  C->SetBranchAddress( branchnames["hit.ADCmaxstripU"].Data(), ADCmaxstripU);
  C->SetBranchAddress( branchnames["hit.ADCmaxstripV"].Data(), ADCmaxstripV);
  C->SetBranchAddress( branchnames["hit.ADCU"].Data(), ADCU);
  C->SetBranchAddress( branchnames["hit.ADCV"].Data(), ADCV);
  C->SetBranchAddress( branchnames["hit.ADCavg"].Data(), ADCavg);
  C->SetBranchAddress( branchnames["hit.ADCasym"].Data(), ADCasym);
  C->SetBranchAddress( branchnames["hit.ADCU_deconv"].Data(), DeconvADCU);
  C->SetBranchAddress( branchnames["hit.ADCV_deconv"].Data(), DeconvADCV);
  C->SetBranchAddress( branchnames["hit.ADCasym_deconv"].Data(), ADCasymDeconv);
 
  double UtimeMaxStrip[MAXNHITS],VtimeMaxStrip[MAXNHITS];
  double UtimeMaxStripDeconv[MAXNHITS],VtimeMaxStripDeconv[MAXNHITS];
  double UtimeMaxStripFit[MAXNHITS],VtimeMaxStripFit[MAXNHITS];
 
  double Utime[MAXNHITS],Vtime[MAXNHITS];
  double UtimeDeconv[MAXNHITS],VtimeDeconv[MAXNHITS];
  double UtimeFit[MAXNHITS],VtimeFit[MAXNHITS];
 
  double HitTavg_corr[MAXNHITS];
 
  C->SetBranchAddress( branchnames["hit.Tavg_corr"].Data(), &HitTavg_corr[0]);
  C->SetBranchAddress( branchnames["hit.UtimeMaxStrip"].Data(), UtimeMaxStrip);
  C->SetBranchAddress( branchnames["hit.VtimeMaxStrip"].Data(), VtimeMaxStrip);
  C->SetBranchAddress( branchnames["hit.UtimeMaxStripDeconv"].Data(), UtimeMaxStripDeconv);
  C->SetBranchAddress( branchnames["hit.VtimeMaxStripDeconv"].Data(), VtimeMaxStripDeconv);
  C->SetBranchAddress( branchnames["hit.UtimeMaxStripFit"].Data(), UtimeMaxStripFit);
  C->SetBranchAddress( branchnames["hit.VtimeMaxStripFit"].Data(), VtimeMaxStripFit);
 
  C->SetBranchAddress( branchnames["hit.Utime"].Data(), Utime);
  C->SetBranchAddress( branchnames["hit.Vtime"].Data(), Vtime);
  C->SetBranchAddress( branchnames["hit.UtimeDeconv"].Data(), UtimeDeconv);
  C->SetBranchAddress( branchnames["hit.VtimeDeconv"].Data(), VtimeDeconv);
  C->SetBranchAddress( branchnames["hit.UtimeFit"].Data(), UtimeFit);
  C->SetBranchAddress( branchnames["hit.VtimeFit"].Data(), VtimeFit);
  
  double deltat[MAXNHITS], deltatFit[MAXNHITS], deltatDeconv[MAXNHITS];
 
  C->SetBranchAddress( branchnames["hit.deltat"].Data(), deltat);
  C->SetBranchAddress( branchnames["hit.deltat_deconv"].Data(), deltatDeconv);
  C->SetBranchAddress( branchnames["hit.deltat_fit"].Data(), deltatFit);

  if (detname == "sbs.gemFPP") {
    outfilename="TrackingCuts_GEp_GEMFPP.root";
    C->SetBranchStatus("sbs.gemFT.track.nhits",1);
  } else {
    outfilename="TrackingCuts_GEp_GEMFT.root";
  }

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

//  globalcut = "1";
  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut.GetTitle(), C );

  int treenum=0, currenttreenum=0;

  cout << "Cut = " << globalcut << endl;
  cout << "Cut Dim = " << GlobalCut->GetNdim() << endl;
  cout << "Cut ??? = " << GlobalCut << endl;

  // Here we start doing our main loops
  while( C->GetEntry( nevent++ ) && nevent){
    currenttreenum = C->GetTreeNumber();
    if( nevent == 1 || currenttreenum != treenum ){
      treenum = currenttreenum;

      if (GlobalCut) {
        delete GlobalCut;
        GlobalCut = nullptr;
      }
      GlobalCut = new TTreeFormula( "GlobalCut", globalcut.GetTitle(), C );
      GlobalCut->UpdateFormulaLeaves();
    }
   
    if( nevent % 100000 == 0 ) cout << nevent << endl;
    bool passedcut = GlobalCut->EvalInstance(0) != 0;

//     cout << "passed cut, ntracks, ngoodhits = " << passedcut << ", "
//     	 << ntracks << ", " << ngoodhits << endl;


    
    if( passedcut && int(ntracks) >= 1){

      //Grab trigger TDC time: 
      double ttrig = 0.0;
      for( int itrig=0; itrig<Ntrig; itrig++ ){
	if( int(trig_elemID[itrig]) == 5 ){
	  ttrig = trig_tdc[itrig];
	}
      }
//      cout << "Let's try filling our first batch of histograms" << endl;

      double thcp = (xbcp[0]-xfcp[0])/(zbcp[0]-zfcp[0]);
      double phcp = (ybcp[0]-yfcp[0])/(zbcp[0]-zfcp[0]);
//      if( EPS >= 0.2 ){ //tracking constraint histos always require preshower regardless what the user defined?
      hdxfcp->Fill( xfp[0] - xfcp[0] );
      hdyfcp->Fill( yfp[0] - yfcp[0] );
      hdxbcp->Fill( xfp[0]+thfp[0]*(zbcp[0]-zfcp[0]) - xbcp[0] );
      hdybcp->Fill( yfp[0]+phfp[0]*(zbcp[0]-zfcp[0]) - ybcp[0] );
  
      hdxdyfcp->Fill( yfp[0]-yfcp[0], xfp[0]-xfcp[0] );
      hdxdybcp->Fill( yfp[0]+phfp[0]*(zbcp[0]-zfcp[0]) - ybcp[0],
      		xfp[0]+thfp[0]*(zbcp[0]-zfcp[0]) - xbcp[0] );
 
      hdthcp->Fill( thfp[0]-thcp );
      hdphcp->Fill( phfp[0]-phcp );

      htrackt0->Fill( trackt0[0] );
//      }
      
      for( int ihit=0; ihit<int(ngoodhits); ihit++ ){
	if( int(trackindex[ihit]) == 0 && nstripu[ihit]>1&&nstripv[ihit]>1 ){
	  htavg_corr_vs_ttrig_allhits->Fill( ttrig, HitTavg_corr[ihit] );

	  ( (TH2D *) (*htavg_corr_vs_ttrig_by_module)[module[ihit]] )->Fill( ttrig, HitTavg_corr[ihit] );
	  
	  hADCmaxsamp_vs_module->Fill( module[ihit], 0.5*(ADCmaxsampU[ihit]+ADCmaxsampV[ihit]) );
	  hADCmaxstrip_vs_module->Fill( module[ihit], 0.5*(ADCmaxstripU[ihit]+ADCmaxstripV[ihit]) );

	  hADCclust_vs_module->Fill( module[ihit], ADCavg[ihit] );

	  //For timing cuts we require hits to pass higher ADC thresholds:
	  if( ADCavg[ihit] >= 1500.0 ){
	    hADCasym_mod->Fill( module[ihit], ADCasym[ihit] );
	    hADCasym_deconv_mod->Fill( module[ihit], ADCasymDeconv[ihit] );
	    hADCratio_mod->Fill( module[ihit], ADCV[ihit]/ADCU[ihit] );
	    hADCratio_deconv_mod->Fill( module[ihit], DeconvADCV[ihit]/DeconvADCU[ihit] );
	    
	    hdeltat_mod->Fill( module[ihit], deltat[ihit] );
	    hdeltat_deconv_mod->Fill( module[ihit], deltatDeconv[ihit] );
	    hdeltat_fit_mod->Fill( module[ihit], deltatFit[ihit] );
	  }

	  if( ADCmaxsampU[ihit] >= 100. && ADCmaxsampV[ihit] >= 100. &&
	      ADCmaxstripU[ihit] >= 250. && ADCmaxstripV[ihit] >= 250. &&
	      ADCU[ihit] >= 600. && ADCV[ihit] >= 600. &&
	      ADCavg[ihit] >= 750. ){
	    hmaxstrip_tU_mod->Fill( module[ihit], UtimeMaxStrip[ihit] );
	    hmaxstrip_tV_mod->Fill( module[ihit], VtimeMaxStrip[ihit] );
	    
	    hmaxstrip_tU_deconv_mod->Fill( module[ihit], UtimeMaxStripDeconv[ihit] );
	    hmaxstrip_tV_deconv_mod->Fill( module[ihit], VtimeMaxStripDeconv[ihit] );

	    hmaxstrip_tU_fit_mod->Fill( module[ihit], UtimeMaxStripFit[ihit] );
	    hmaxstrip_tV_fit_mod->Fill( module[ihit], VtimeMaxStripFit[ihit] );
	    
	    htU_mod->Fill( module[ihit], Utime[ihit] );
	    htV_mod->Fill( module[ihit], Vtime[ihit] );
	    
	    htU_deconv_mod->Fill( module[ihit], UtimeDeconv[ihit] );
	    htV_deconv_mod->Fill( module[ihit], VtimeDeconv[ihit] );
	    
	    htU_fit_mod->Fill( module[ihit], UtimeFit[ihit] );
	    htV_fit_mod->Fill( module[ihit], VtimeFit[ihit] );

	    hdt_corr_mod->Fill( module[ihit], HitTavg_corr[ihit]-trackt0[0] );
	  }
	}
      }
    }
  }
  cout << "done debugging the first major while loop" << endl;
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
  cout << "Let's look at timing cuts now" << endl; 
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
      dbfile << dbline.Format("%s.m%d.threshold_sample = %12.5g", detname.Data(), imod, threshsamp[0] ) << endl;

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
      dbfile << dbline.Format("%s.m%d.threshold_stripsum = %12.5g", detname.Data(), imod, threshstrip[0] ) << endl;

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
      dbfile << dbline.Format("%s.m%d.threshold_clustersum = %12.5g", detname.Data(), imod, thresh[0] ) << endl << endl;

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
      dbfile_tcuts << dbline.Format("%s.m%d.maxstrip_t0 = %10.4g %10.4g", detname.Data(), imod,maxstrip_t0[imod][0], maxstrip_t0[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("%s.m%d.maxstrip_tsigma = %10.4g %10.4g", detname.Data(), imod,maxstrip_tsigma[imod][0], maxstrip_tsigma[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("%s.m%d.maxstrip_tcut = %10.4g %10.4g", detname.Data(), imod, nsigma, nsigma ) << endl << endl;
      
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
      dbfile_tcuts << dbline.Format("%s.m%d.maxstrip_t0_deconv = %10.4g %10.4g", detname.Data(), imod,maxstrip_t0_deconv[imod][0], maxstrip_t0_deconv[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("%s.m%d.maxstrip_tsigma_deconv = %10.4g %10.4g", detname.Data(), imod,maxstrip_tsigma_deconv[imod][0], maxstrip_tsigma_deconv[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("%s.m%d.maxstrip_tcut_deconv = %10.4g %10.4g", detname.Data(), imod, nsigma, nsigma ) << endl << endl;

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
      dbfile_tcuts << dbline.Format("%s.m%d.maxstrip_t0_fit = %10.4g %10.4g", detname.Data(), imod,maxstrip_t0_fit[imod][0], maxstrip_t0_fit[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("%s.m%d.maxstrip_tsigma_fit = %10.4g %10.4g", detname.Data(), imod,maxstrip_tsigma_fit[imod][0], maxstrip_tsigma_fit[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("%s.m%d.maxstrip_tcut_fit = %10.4g %10.4g", detname.Data(), imod, nsigma, nsigma ) << endl << endl;


      ///// Now we create the histograms for the pdf results//////
      vector<double> resultsU, resultsV;
      resultsU.push_back(maxstrip_t0_fit[imod][0]);
      resultsU.push_back(maxstrip_tsigma_fit[imod][0]);
      resultsV.push_back(maxstrip_t0_fit[imod][1]);
      resultsV.push_back(maxstrip_tsigma_fit[imod][1]);
      
      TH1D *htU = htempU;
      htU->SetName(Form("htstripfitU_m%i",imod));
      htU->SetTitle(Form("Module %i U Fit Time;fit time (ns)",imod));
      c2->cd(5);
      DrawHist(htU,resultsU,2);
      fitfuncU->Draw("same");

      TH1D *htV = htempV;
      htV->SetName(Form("htstripfitV_m%i",imod));
      htV->SetTitle(Form("Module %i V Fit Time;fit time (ns)",imod));
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
      dbfile_tcuts << dbline.Format("%s.m%d.deltat_sigma = %12.5g", detname.Data(), imod, sigma ) << endl;
      dbfile_tcuts << dbline.Format("%s.m%d.deltat_cut = %12.5g", detname.Data(), imod, nsig_dt*sigma ) << endl;

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
      dbfile_tcuts << dbline.Format("%s.m%d.deltat_sigma_deconv = %12.5g", detname.Data(), imod, sigma ) << endl;
      dbfile_tcuts << dbline.Format("%s.m%d.deltat_cut_deconv = %12.5g", detname.Data(), imod, nsig_dt*sigma ) << endl;

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
      dbfile_tcuts << dbline.Format("%s.m%d.deltat_sigma_fit = %12.5g", detname.Data(), imod, sigma ) << endl;
      dbfile_tcuts << dbline.Format("%s.m%d.deltat_cut_fit = %12.5g", detname.Data(), imod, nsig_dt*sigma ) << endl << endl;

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

      dbfile_tcuts << dbline.Format("%s.m%d.HitTimeMean = %12.5g %12.5g", detname.Data(), imod, tmean[imod][0], tmean[imod][1] ) << endl;
      dbfile_tcuts << dbline.Format("%s.m%d.HitTimeSigma = %12.5g %12.5g", detname.Data(), imod, tsigma[imod][0], tsigma[imod][1] ) << endl;

      ///// Now we create the histograms for the pdf results//////
      vector<double> resultsU, resultsV;
      resultsU.push_back(tmean[imod][0]);
      resultsU.push_back(tsigma[imod][0]);
      resultsV.push_back(tmean[imod][1]);
      resultsV.push_back(tsigma[imod][1]);

      TH1D *htU = htempU;
      htU->SetName(Form("htclustU_m%i",imod));
      htU->SetTitle(Form("Module %i U Cluster Time;time (ns)",imod));
      c4->cd(1);
      DrawHist(htU,resultsU,2);
      fitfuncU->Draw("same");

      TH1D *htV = htempV;
      htV->SetName(Form("htclustV_m%i",imod));
      htV->SetTitle(Form("Module %i V Cluster Time;time (ns)",imod));
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

      dbfile_tcuts << dbline.Format("%s.m%d.HitTimeMeanDeconv = %12.5g %12.5g", detname.Data(), imod, tmean[imod][2], tmean[imod][3] ) << endl;
      dbfile_tcuts << dbline.Format("%s.m%d.HitTimeSigmaDeconv = %12.5g %12.5g", detname.Data(), imod, tsigma[imod][2], tsigma[imod][3] ) << endl;

      ///// Now we create the histograms for the pdf results//////
      vector<double> resultsU, resultsV;
      resultsU.push_back(tmean[imod][2]);
      resultsU.push_back(tsigma[imod][2]);
      resultsV.push_back(tmean[imod][3]);
      resultsV.push_back(tsigma[imod][3]);

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

      dbfile_tcuts << dbline.Format("%s.m%d.HitTimeMeanFit = %12.5g %12.5g", detname.Data(), imod, tmean[imod][4], tmean[imod][5] ) << endl;
      dbfile_tcuts << dbline.Format("%s.m%d.HitTimeSigmaFit = %12.5g %12.5g", detname.Data(), imod, tsigma[imod][4], tsigma[imod][5] ) << endl << endl;

      ///// Now we create the histograms for the pdf results//////
      vector<double> resultsU, resultsV;
      resultsU.push_back(tmean[imod][4]);
      resultsU.push_back(tsigma[imod][4]);
      resultsV.push_back(tmean[imod][5]);
      resultsV.push_back(tsigma[imod][5]);

      TH1D *htU = htempU;
      htU->SetName(Form("htclustfitU_m%i",imod));
      htU->SetTitle(Form("Module %i U Cluster Fit Time;fit time (ns)",imod));
      c4->cd(5);
      DrawHist(htU,resultsU,2);
      fitfuncU->Draw("same");

      TH1D *htV = htempV;
      htV->SetName(Form("htclustfitV_m%i",imod));
      htV->SetTitle(Form("Module %i V Cluster Fit Time;fit time (ns)",imod));
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
      dbfile_tcuts << dbline.Format("%s.m%d.ADCasym_sigma = %12.5g", detname.Data(), imod, sigma ) << endl;
      
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
      dbfile_tcuts << dbline.Format("%s.m%d.ADCratio_sigma = %12.5g", detname.Data(), imod, sigma ) << endl;
      
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
      cout << "did I make it to the htemp loop?" << endl;
      TF1 *fitfunc = (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus"));
      double mean = fitfunc->GetParameter("Mean");
      double sigma = fitfunc->GetParameter("Sigma");
      
      TString dbline;
      dbfile_tcuts << dbline.Format("%s.m%d.ADCratio_sigma_deconv = %12.5g", detname.Data(), imod, sigma ) << endl << endl;

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
      dbfile_tcuts << dbline.Format("%s.m%d.sigma_tcorr = %12.5g", detname.Data(), imod, sigma ) << endl << endl;
      
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
  dbfile_tcuts << endl << dbline.Format("%s.sigmatrackt0 = %12.5g", detname.Data(), sigma ) << endl << endl;
  
  
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
