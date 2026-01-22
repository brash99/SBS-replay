#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TPad.h"
#include "TSystem.h"
#include <vector>
#include "TSpectrum.h"

using namespace std;

void GetSaturationMinByAPV(const char *rootfilename, const char *detname="sbs.gemFT", int nmodules=14, const char *outfilename="SaturationMinTemp.root", int minevents = 10){
  TFile *fin = new TFile(rootfilename,"READ");

  if( !fin ) return;

  TFile *fout = new TFile(outfilename,"RECREATE");
  
  //  TString histname;

  TString sdetname = detname;
  sdetname.ReplaceAll(".","_");

  TString dbfilename;
  dbfilename.Form("db_%s_satmin.dat", sdetname.Data() ); 
  
  ofstream dbfile(dbfilename.Data());

  TCanvas *c1 = new TCanvas("c1","c1",1600,1200);

  c1->Divide(2,2,0.001,0.001);
  
  //  c1->SetLogy();

  TLine L;
  L.SetLineWidth(2);
  L.SetLineColor(6);


  
  //TH1D *projtemp;
  
  for( int mod=0; mod<nmodules; mod++ ){

    cout << "About to process module " << mod << endl;
    TString histnameU;
    histnameU.Form( "hrawADCs_by_APV_U_nopedsub_%s_m%d_goodCM", sdetname.Data(), mod );

    TH2D *htempU, *htempV;
    fin->GetObject( histnameU.Data(), htempU );

    c1->cd(1)->SetLogz();
    htempU->SetTitle( Form("Module %d U strips; APV num (by pos) ; raw ADC (no ped. sub)",mod ) );
    htempU->DrawCopy("colz");


    
    int nAPVsU = htempU->GetNbinsX();
    vector<double> satmin_by_APV_U(nAPVsU,0.0);

    double minval = 0.0;

    //   c1->cd();

    c1->cd(3)->SetLogy();
    
    for( int iAPV=1; iAPV<=nAPVsU; iAPV++ ){
      TString projname;
      projname.Form( "hrawADC_nopedsub_%s_m%d_APV_U_%d_goodCM", sdetname.Data(), mod, iAPV-1 );

      c1->cd(3)->SetLogy();
      
      TH1D *projtemp = htempU->ProjectionY( projname.Data(), iAPV, iAPV );
      
      if( !projtemp ) continue;

      projtemp->SetTitle( Form( "Module %d U strip APV%d;raw ADC (no ped sub)", mod, iAPV-1 ) );
      
      if( projtemp->GetEntries() < 1000 ) continue;
      
      int binmax = projtemp->GetMaximumBin();

      int firstbinwithhits = 1;

      while( firstbinwithhits < binmax && projtemp->GetBinContent(firstbinwithhits) < 5 ){
	firstbinwithhits++;
      }

      minval = projtemp->GetBinCenter( firstbinwithhits );
      
      TSpectrum refspec;
      int npeaksref = refspec.Search( projtemp, 3.0, "", 0.002 );

      double *peakpos = refspec.GetPositionX();
      double *peakheight = refspec.GetPositionY();

      double firstpeakpos = 0.0;
      //Find the peak with the lowest ADC value that is within ~125 ADC of firstbinwithhits:
      int ipeak_low = -1;
      double xpeak_low = 5000.0;
	
      if( npeaksref > 1 ){
	
	for( int ipeak=0; ipeak<npeaksref; ipeak++ ){
	  if( ipeak_low < 0 ||
	      (peakpos[ipeak] < xpeak_low )){
	    xpeak_low = peakpos[ipeak];
	    ipeak_low = ipeak;
	  }
	}
	
      }

      bool satminpeakfound = false;
      
      if( ipeak_low >= 0 && xpeak_low > minval && xpeak_low < minval + 200.0 && peakheight[ipeak_low] >= minevents && peakheight[ipeak_low] < projtemp->GetBinContent(binmax) ){
	firstpeakpos = xpeak_low;
	satminpeakfound = true;
      }

      double ADCmin = minval - 50.0;
      
      if( satminpeakfound ){
	
	int binfirst = projtemp->FindBin( firstpeakpos );

	//double firstpeakheight = projtemp->GetBinContent( binfirst );
	//Let's find the lowest point between the maximum and the negative saturation peak position (if applicable);
	
	int binmin = binfirst;
	int testbin = binfirst;
	
	double mincontent = projtemp->GetBinContent(binfirst);
	
	while( testbin < binmax ){
	  double testbincontent = projtemp->GetBinContent(testbin);
	  if( testbincontent < mincontent ){
	    binmin = testbin;
	    mincontent = testbincontent;
	  }
	  testbin++;
	}

	//	minval = projtemp->GetBinCenter(firstbinwithhits);
      
	//      double ADCmin = projtemp->GetBinLowEdge(binmin) + 150.0;
	ADCmin = projtemp->GetBinCenter( binmin );
      }
      
      c1->cd(3);

      bool okay = false;
      while( !okay ){
      
	projtemp->Draw();
	
	L.DrawLine( ADCmin, projtemp->GetYaxis()->GetXmin(), ADCmin, projtemp->GetMaximum() ); 
	
	c1->Update();
	gPad->Modified();
	gSystem->ProcessEvents();
	
	cout << "Saturation min okay? (y/n/s)";
	
	TString reply;
	reply.ReadLine(cin,kFALSE);

	okay = true;
	
	if( reply.BeginsWith("n") ){ //enter manually:
	  okay = false;
	  cout << "Enter new sat min value (enter to use minval)";
	  reply.ReadLine(cin,kFALSE);
	  ADCmin = reply.IsFloat() ? reply.Atof() : minval - 50.0;
	}

	//default to zero:
	if( reply.BeginsWith("s") ) ADCmin = 0.0;
      }
      
      satmin_by_APV_U[iAPV-1] = ADCmin;
      c1->cd(1);

      L.DrawLine( iAPV-1.5,ADCmin, iAPV-0.5,ADCmin );
      cout << "done" << endl;

      fout->cd();
      projtemp->Write();
      //      projtemp->Delete();
    }

    TString histnameV;
    histnameV.Form( "hrawADCs_by_APV_V_nopedsub_%s_m%d_goodCM", sdetname.Data(), mod );

    fin->GetObject( histnameV.Data(), htempV );

    c1->cd(2)->SetLogz();
    htempV->SetTitle( Form("Module %d V strips; APV num (by pos) ; raw ADC (no ped. sub)",mod ) );
    htempV->DrawCopy("colz");
    
    int nAPVsV = htempV->GetNbinsX();
    vector<double> satmin_by_APV_V(nAPVsV,0.0);

    minval = 0.0;

    c1->cd(4);
    
    for( int iAPV=1; iAPV<=nAPVsV; iAPV++ ){
      TString projname;
      cout << "About to project APV " << iAPV << endl;
      
      projname.Form( "hrawADC_nopedsub_%s_m%d_APV_V_%d_goodCM", sdetname.Data(), mod, iAPV-1 );

      cout << "Projecting histogram \"" << projname << "\" onto projtemp" << endl;

      c1->cd(4)->SetLogy();
      
      TH1D *projtemp = htempV->ProjectionY( projname.Data(), iAPV, iAPV );
      
      if( !projtemp ) continue;

      projtemp->SetTitle( Form( "Module %d V strip APV%d;raw ADC (no ped sub)", mod, iAPV-1 ) );
      
      int binmin = 1;

      if( projtemp->GetEntries() < 1000 ) continue;

      int binmax = projtemp->GetMaximumBin();
      
      int firstbinwithhits = 1;

      while ( firstbinwithhits < binmax && projtemp->GetBinContent(firstbinwithhits) < 5 ){
	firstbinwithhits++; 
      }

      minval = projtemp->GetBinCenter( firstbinwithhits );
      
      TSpectrum refspec;
      int npeaksref = refspec.Search( projtemp, 3.0, "", 0.002 );

      double *peakpos = refspec.GetPositionX();
      double *peakheight = refspec.GetPositionY();

      double firstpeakpos = 0.0;
      //Find the peak with the lowest ADC value
      int ipeak_low = -1;
      double xpeak_low = 5000.0;
	
      if( npeaksref > 1 ){
	
	for( int ipeak=0; ipeak<npeaksref; ipeak++ ){
	  if( ipeak_low < 0 ||
	      (peakpos[ipeak] < xpeak_low )){
	    xpeak_low = peakpos[ipeak];
	    ipeak_low = ipeak;
	  }
	}
	
      }

      bool satminpeakfound = false;
      if( ipeak_low >= 0 && xpeak_low > minval && xpeak_low < minval + 200.0 && peakheight[ipeak_low] >= minevents && peakheight[ipeak_low] < projtemp->GetBinContent(binmax) ){
	firstpeakpos = xpeak_low;
	satminpeakfound = true;
      }

      double ADCmin = minval - 50.0;

      cout << "Found sat min peak = " << satminpeakfound << endl;
      
      if( satminpeakfound ){
	
	int binfirst = projtemp->FindBin( firstpeakpos );
	
	int binmin = binfirst;
	int testbin = binfirst;
	double mincontent = projtemp->GetBinContent(binfirst);

	while( testbin < binmax ){
	  double testbincontent = projtemp->GetBinContent(testbin);
	  if( testbincontent < mincontent ){
	    binmin = testbin;
	    mincontent = testbincontent;
	  }
	  testbin++;
	}

	ADCmin = projtemp->GetBinCenter( binmin );
      }

      
      
      c1->cd(4);
      
      bool okay = false;
      while( !okay ){
	
	projtemp->Draw();
	
	L.DrawLine( ADCmin, projtemp->GetYaxis()->GetXmin(), ADCmin, projtemp->GetMaximum() ); 
	
	c1->Update();
	gPad->Modified();
	gSystem->ProcessEvents();
	
	cout << "Saturation min okay? (y/n/s)";
	
	TString reply;
	reply.ReadLine(cin,kFALSE);
	
	okay = true;
	
	if( reply.BeginsWith("n") ){ //enter manually:
	  okay = false;
	  cout << "Enter new sat min value (enter to use minval or s to skip)";
	  reply.ReadLine(cin,kFALSE);
	  ADCmin = reply.IsFloat() ? reply.Atof() : minval - 50.0;
	}
	
	//default to zero:
	if( reply.BeginsWith("s") ) ADCmin = 0.0;
      }
      
      satmin_by_APV_V[iAPV-1] = ADCmin;

      c1->cd(2);
      L.DrawLine( iAPV-1.5, ADCmin, iAPV-0.5, ADCmin );

      c1->Update();
      gPad->Modified();
      gSystem->ProcessEvents();
      
      cout << "done" << endl;
      fout->cd();
      projtemp->Write();
      //      projtemp->Delete();
      cout << "Wrote projtemp to fout" << endl;
    }

    //Before we delete the temporary histograms, print the page to pdf

    TString pdffilename = outfilename;
    pdffilename.ReplaceAll(".root",".pdf");

    if( mod== 0 ) pdffilename += "(";

    if( mod+1 == nmodules ) pdffilename += ")";

    c1->Print(pdffilename.Data());
    
    htempU->Delete();
    htempV->Delete();
    
    cout << "module " << mod << " done" << endl;

    TString sdbkey;
    //The database entries will be ordered by "pos":
    sdbkey.Form( "%s.m%d.rawADCminU = ", detname, mod );

    dbfile << sdbkey << endl;
    
    for( int iAPV=0; iAPV<satmin_by_APV_U.size(); iAPV++ ){
      TString smin;

      smin.Form( " %10.6g ", satmin_by_APV_U[iAPV] );
      dbfile << smin;
      if( iAPV > 0 && (iAPV+1) % 8 == 0 ) dbfile << endl;
    }
    dbfile << endl << endl;
    
    sdbkey.Form( "%s.m%d.rawADCminV = ", detname, mod );

    dbfile << sdbkey << endl;
    
    for( int iAPV=0; iAPV<satmin_by_APV_V.size(); iAPV++ ){
      TString smin;

      smin.Form( " %10.6g ", satmin_by_APV_V[iAPV] );
      dbfile << smin;
      if( iAPV > 0 && (iAPV+1) % 8 == 0 ) dbfile << endl;
    }
    
    dbfile << endl << endl;
    
  }
  
}
