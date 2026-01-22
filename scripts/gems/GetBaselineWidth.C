#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "TCanvas.h"
#include "TPad.h"
#include "TSystem.h"
#include "TGraph.h"

void FitGaus_FWHM( TH1D *htest, double thresh_low=0.2, double thresh_high=0.5 ){
  int binmax = htest->GetMaximumBin();
  int binlow = binmax, binhigh = binmax;

  double max = htest->GetBinContent(binmax);

  while( htest->GetBinContent(binlow) >= thresh_low*max && binlow > 1 ){binlow--;}
  while( htest->GetBinContent(binhigh) >= thresh_high*max && binhigh < htest->GetNbinsX() ){ binhigh++; }

  double xlow = htest->GetBinLowEdge(binlow);
  double xhigh = htest->GetBinLowEdge(binhigh+1);

  htest->Fit("gaus","qS","",xlow, xhigh);
}

void GetBaselineWidth(const char *rootfilename, const char *outfilename="BaselineWidthTemp.root", const char *detname="sbs.gemFT", int nmodules=14){
  
  TFile *fout = new TFile(outfilename,"RECREATE");
  
  TFile *fin = new TFile(rootfilename,"READ");

  TString sdetname = detname;
  sdetname.ReplaceAll(".","_");

  TCanvas *c1 = new TCanvas("c1","c1",1600,800);
  c1->Divide(2,1,0.001,0.001);

  double default_width = 20.0;
  
  vector<double> widthU(nmodules, default_width),widthV(nmodules, default_width);
  vector<double> width_avg(nmodules,default_width), module(nmodules);

  TString datfilename = outfilename;
  datfilename.ReplaceAll(".root",".dat");
  
  ofstream resultsfile(datfilename.Data());

  TString dbfilename = outfilename;
  dbfilename.ReplaceAll(".root",".dat");
  dbfilename.Prepend("db_");

  ofstream dbfile(dbfilename.Data());
  
  for( int imod=0; imod<nmodules; imod++ ){

    module[imod] = imod;
    
    TString hnameU = Form( "hADCpedsubU_allstrips_%s_m%d_goodCM", sdetname.Data(), imod );
    TString hnameV = Form( "hADCpedsubV_allstrips_%s_m%d_goodCM", sdetname.Data(), imod );

    TH1D *htempU,*htempV;

    fin->GetObject( hnameU.Data(), htempU );
    fin->GetObject( hnameV.Data(), htempV );

    if( htempU != nullptr && htempV != nullptr ){
    
      if( htempU->GetEntries() > 1000 ){
	
	c1->cd(1)->SetLogy();
	htempU->Draw();
	FitGaus_FWHM( htempU, 0.3, 0.5 );
	htempU->Draw();

	widthU[imod] = ( (TF1*) (htempU->GetListOfFunctions()->FindObject("gaus")) )->GetParameter("Sigma");
	
	gPad->Modified();
	c1->Update();
	gSystem->ProcessEvents();
      }
      
      if( htempV->GetEntries() > 1000 ){
	c1->cd(2)->SetLogy();
	htempV->Draw();
	FitGaus_FWHM( htempV, 0.3, 0.5 );
	htempV->Draw();

	widthV[imod] = ( (TF1*) (htempV->GetListOfFunctions()->FindObject("gaus")) )->GetParameter("Sigma");
	
	gPad->Modified();
	c1->Update();
	gSystem->ProcessEvents();
      }

      width_avg[imod] = 0.5*(widthU[imod]+widthV[imod]);
      
      TString pdfname = outfilename;
      pdfname.ReplaceAll(".root",".pdf");
      
      if( imod == 0 ) pdfname += "(";
      //if( imod+1 == nmodules ) pdfname += ")";

      TString datline = Form(" %15.6g %15.6g %15.6g ", widthU[imod], widthV[imod], width_avg[imod] );

      resultsfile << datline << endl;

      TString dbline = Form("%s.m%d.threshold_sample = %15.5g", detname, imod, 3.0*width_avg[imod]);
      dbfile << dbline << endl;
      dbline = Form( "%s.m%d.threshold_stripsum = %15.5g", detname, imod, 2.5*3.0*width_avg[imod]);
      dbfile << dbline << endl;
      dbline = Form( "%s.m%d.threshold_clustersum = %15.5g", detname, imod, 4.0*3.0*width_avg[imod]);
      dbfile << dbline << endl;
      
      c1->Print(pdfname.Data(),"pdf");
    }
  }

  TGraph *gwU = new TGraph(nmodules,&(module[0]),&(widthU[0]));
  TGraph *gwV = new TGraph(nmodules,&(module[0]),&(widthV[0]));
  TGraph *gwavg = new TGraph(nmodules,&(module[0]),&(width_avg[0]));

  gwU->SetMarkerStyle(24);
  gwU->SetMarkerColor(2);
  gwU->SetLineColor(2);

  gwV->SetMarkerStyle(25);
  gwV->SetMarkerColor(4);
  gwV->SetLineColor(4);

  gwavg->SetMarkerStyle(20);
  gwavg->SetMarkerColor(1);
  gwavg->SetLineColor(1);

  c1->Clear();

  c1->cd();

  gwavg->Draw("AP");
  gwavg->GetYaxis()->SetTitle("Baseline width (ADC, individual samples)");
  gwavg->GetXaxis()->SetTitle(Form("%s module",detname));
  gwavg->GetYaxis()->SetRangeUser(0.0,150.0);

  TString title = outfilename;
  title.ReplaceAll(".root","");
  
  
  gwavg->SetTitle(title.Data());
  
  gwU->Draw("PSAME");
  gwV->Draw("PSAME");
  
  TString pdfname = outfilename;
  pdfname.ReplaceAll(".root",".pdf");

  pdfname += ")";

  c1->Print(pdfname.Data(),"pdf");
    
  fout->cd();
  gwU->Write("gwU");
  gwV->Write("gwV");
  gwavg->Write("gwavg");
  
  

}
