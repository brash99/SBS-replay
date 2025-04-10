#include "gmn_tree.C"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCut.h"
#include <iostream>
#include "TString.h"
#include "TCanvas.h"
#include <vector>
#include "TPad.h"
#include "TSystem.h"
#include "TLine.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TMinuit.h"
#include "TSpectrum.h"
#include "TClonesArray.h"


void GetRFoffsets( const char *configfilename, const char *outputfilename="HodoRFoffsets_temp.root", double RF_freq = 249.5e6, int refID = 40, double reft0 = 1.9895, double etof0=15.1667 ){

  TFile *fout = new TFile(outputfilename,"RECREATE");
  
  ifstream configfile(configfilename);
  if( !configfile ) return;

  TString currentline;

  TChain *C = new TChain("T");
  
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add( currentline.Data() );
    }
  }
  
  //TFile *fout = new TFile(outputfilename,"RECREATE");
  
  
  //  C->Add(rootfilename);

  gmn_tree *T = new gmn_tree(C);

  double bunchspace = 1.0e9/RF_freq; // ns 

  //Total time range is going to be 3*bunchspace, we want about 8 bins/ns so 
  
  //int nbins = ;
  int nbins=100;
  TH2D *hdT_hodo_ID = new TH2D("hdT_hodo_ID", "; hodoscope bar ID ; t_{mean}^{raw}-t_{RF} (ns)", 90,-0.5,89.5, nbins, -1.5*bunchspace, 1.5*bunchspace);

  long nevent = 0;

  C->SetBranchStatus("*",0);
  C->SetBranchStatus("bb.tr.*",1);
  C->SetBranchStatus("bb.hodotdc.*",1);
  C->SetBranchStatus("bb.tdctrig.*",1);
  C->SetBranchStatus("bb.ps.*",1);
  C->SetBranchStatus("bb.sh.*",1);

  //double wL[90];
  //double wR[90];
  //double t0[90];
  //double vscint[90];

  double walk_default = (0.2113+0.2794)/2.0;
  double vscint_default = 0.1718;
  
  vector<double> wL(90,walk_default);
  vector<double> wR(90,walk_default);
  vector<double> t0(90,0.0);
  vector<double> vscint(90,vscint_default);
  
  //Read parameters from initial fit:
  // ifstream wLfile("hodowL_temp.dat");
  // ifstream wRfile("hodowR_temp.dat");
  // Let's modify this so that it reads all the walk coefficients from a single file:
  ifstream walkfile("hodowalkcorr_default.dat");
  ifstream vscintfile("hodovscint_default.dat");
  ifstream t0file("hodomeantime_default.dat");

  if( walkfile && vscintfile && t0file ){
    for( int i=0; i<90; i++ ){
      vscintfile >> vscint[i];
      t0file >> t0[i];
      walkfile >> wL[i];
      if( wL[i] < 0 ) wL[i] *= -1.; //The convention in this code requires the walk correction slope to be positive
    }

    for( int i=0; i<90; i++ ){
      walkfile >> wR[i];
      if( wR[i] < 0 ) wR[i] *= -1.; 
    }
    
  }
 
  
  
  // Determine mean time offsets in an (RF-corrected) way;
  // Meantimes are already corrected for time walk, which itself contains a redundant offset.
  // The idea is to look at (raw) mean-time differences between neighboring bars in the same cluster.
  // Minimize sum of squared corrected time differences between 
  // chi^2 = sum_{i=1}^Nevent [sum_{j=1}^{Nhits_i-1} [sum_{k=j+1}^N_hits_i (t_j^i-t_k^i)^2]]
  // Where t_j^i,corr = t_mean_j^i - t0ji
  // Then dchi^2 / dt0_m = 0 = 2 * sum_i,j,k ( tmean_j^i - tmean_k^i - t0ji + t0ki ) * ( delta mk - delta mj )
  // = 2 * [ sum_i,j (tji-tmi + -t0ji + t0mi) ] - sum_i,k (t_mi - t_ki - t0mi + t0ki ) ] = 0
  // coefficient of 
  // TMatrixD M(90,90);
  // TVectorD b(90);

  // for( int i=0; i<90; i++ ){
  //   b(i) = 0.0;
  //   for( int j=0; j<90; j++ ){
  //     M(i,j) = 0.0;
  //     M(j,i) = 0.0;
  //   }
  // }

  double zhodo = 1.854454; //meters
  double Lbar_hodo = 0.6; //meters

  //We're going to do this a little bit differently now; 
  
  while( C->GetEntry(nevent++)){

    if( nevent % 10000 == 0 ) cout << "event " << nevent << endl;
    
    if( T->bb_tr_n >= 1. && T->bb_hodotdc_clus_trackindex[0]==0. && T->bb_ps_e>0.2 && T->bb_tdctrig_rftime != 0. && T->bb_tdctrig_rftime < 1.e10 && T->bb_tdctrig_trigtime < 1.e10 ){
      double RFtime = T->bb_tdctrig_rftime + T->bb_tr_vz[0]/0.299792458;

      //RFtime = fmod( RFtime, 160.321 );
      double trigtime = T->bb_tdctrig_trigtime;

      //grab needed track parameters:
      double vz = T->bb_tr_vz[0];
      double pathl = T->bb_tr_pathl[0];
      double etof = pathl/0.299792458-etof0; //electron TOF from vertex to hodo.
      double yhodo = T->bb_tr_y[0]+zhodo*T->bb_tr_ph[0];
      
      for( int ibar=0; ibar<T->bb_hodotdc_clus_size[0]; ibar++ ){
	//double tHODO = T->bb_hodotdc_clus_bar_tdc_meantime[ibar];
	double tleft = T->bb_hodotdc_clus_bar_tdc_tleft[ibar];
	double tright = T->bb_hodotdc_clus_bar_tdc_tright[ibar];
	double totleft = T->bb_hodotdc_clus_bar_tdc_totleft[ibar];
	double totright = T->bb_hodotdc_clus_bar_tdc_totright[ibar];

	double IDHODO = T->bb_hodotdc_clus_bar_tdc_id[ibar];

	double dLEFT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 - yhodo));
	double dRIGHT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 + yhodo));
	
	double tleft_corr = tleft - etof + wL[IDHODO]*totleft - t0[IDHODO] - dLEFT/vscint[IDHODO];
	double tright_corr = tright - etof + wR[IDHODO]*totright - t0[IDHODO] - dRIGHT/vscint[IDHODO];
	
	double tHODO = 0.5*(tleft_corr + tright_corr);

	double deltat = tHODO - (RFtime);
	deltat -= bunchspace * ( floor(deltat/bunchspace) + 0.5 );
	//double deltat = fmod(tHODO - (RFtime - trigtime - bunchspace * floor( ( RFtime-trigtime ) / bunchspace ) ), bunchspace);

	hdT_hodo_ID->Fill( IDHODO, deltat );
				 
	hdT_hodo_ID->Fill( IDHODO, deltat - bunchspace );
	hdT_hodo_ID->Fill( IDHODO, deltat + bunchspace );
	// for( int jbar=ibar+1; jbar<T->bb_hodotdc_clus_size[0]; jbar++ ){
	//   double tj = T->bb_hodotdc_clus_bar_tdc_meantime[jbar];
	//   double IDj = T->bb_hodotdc_clus_bar_tdc_id[jbar];

	//   double deltat_ij = tHODO - tj;

	//   //coefficient of t0i in equation for t0i
	//   // two cases: in both we have (t_ji - t_ki = t_0ji - t_0ki)
	//   // if m = j, coefficient of t_0ji in m equation is +1
	//   // OTOH if m = k, coefficient of t_0ki in m equation is -1
	//   M(IDHODO,IDj) += 1.0;
	//   M(IDj,IDHODO) += 1.0;
	//   M(IDHODO,IDHODO) += -1.0;
	//   M(IDj,IDj) += -1.0;
	//   b(IDHODO) += deltat_ij;
	//   b(IDj) += deltat_ij;	  
	// }
      }
    }
  }

  // TDecompSVD A(M);
  // A.Solve(b);

  // b.Print();
  
  // TString histname;

  // vector<double> T0RF;
  // vector<double> IDRF;
  // vector<double> dT0RF;
  // vector<double> dIDRF;

  TCanvas *c1 = new TCanvas("c1","c1",1200,900);

  TCanvas *c2 = new TCanvas("c2","c2",1200,900);


  c1->cd();
  hdT_hodo_ID->Draw("colz");

  TString histname;

  double toffset_by_bar[90],dtoffset_by_bar[90];
  TClonesArray *dthodo_by_bar = new TClonesArray("TH1D",90);

  for( int i=0; i<90; i++ ){
    TH1D *htemp = hdT_hodo_ID->ProjectionY( histname.Format( "hdThodo_bar%d", i ), i+1, i+1 );
    new( (*dthodo_by_bar)[i] ) TH1D(*htemp);
  }
  
  c2->cd();
  
  double t0temp = reft0;
  TH1D *href = ( (TH1D*) dthodo_by_bar->At(refID) );
    //hdT_hodo_ID->ProjectionY( histname.Format( "hdThodo_bar%d", refID ), refID+1, refID+1 );

  href->Draw();
  
  TSpectrum refspec;
  int npeaksref = refspec.Search( href );

  double *peakpos = refspec.GetPositionX();

  int ipeak_best = -1;
  double mindiff = 1000.0;
  double xpeak_best = t0temp;
  for( int ipeak=0; ipeak<npeaksref; ipeak++ ){
    //    double a = peakpos[ipeak];
    // int bin = 1 + int(a+0.5);
    //double pos = href->GetBinCenter(bin);
    double pos = peakpos[ipeak];
    //cout << "ref bar peak (i,bin,raw pos, pos)=(" << ipeak << ", " << bin << ", " << a << ", " << pos << ")" << endl;
    if( ipeak_best < 0 || fabs(pos - t0temp)<mindiff ){
      ipeak_best = ipeak;
      mindiff = fabs( pos - t0temp );
      xpeak_best = pos;
    }
  }

  href->Fit("gaus","SQ","",xpeak_best-1.0,xpeak_best+1.0);

  gPad->Modified();
  c2->Update();
  gSystem->ProcessEvents();

  TF1 *fitfunc = ( (TF1*) (href->GetListOfFunctions()->FindObject("gaus")) );
  
  reft0 = fitfunc->GetParameter("Mean");
  
  toffset_by_bar[refID] = reft0;
  dtoffset_by_bar[refID] = fitfunc->GetParError(fitfunc->GetParNumber("Mean"));

  double xtest = reft0;
  
  for( int bar=refID+1; bar<90; bar++ ){
    TH1D *htemp = ( (TH1D*) dthodo_by_bar->At(bar) );
    htemp->Draw();

    if( htemp->GetEntries() >= 500.0 ){
      
      gPad->Modified();
      c2->Update();
      gSystem->ProcessEvents();
      
      //double xtest = toffset_by_bar[bar-1];
      
      TSpectrum spectemp;
      int npeakstemp = spectemp.Search( htemp );
      double *xpeakstemp = spectemp.GetPositionX();
      ipeak_best = -1;
      mindiff = 1000.0;
      
      for( int ipeak=0; ipeak<npeakstemp; ipeak++ ){
	//int bin = 1 + int(peakpos[ipeak]+0.5);
	//double pos = htemp->GetBinCenter(bin);
	double pos = xpeakstemp[ipeak];
	if( ipeak_best < 0 || fabs( pos - xtest ) < mindiff ){
	  ipeak_best = ipeak;
	  mindiff = fabs( pos - xtest );
	  xpeak_best = pos;
	}
      }
      
      htemp->Fit("gaus","SQ","",xpeak_best-1.0,xpeak_best+1.0);
      gPad->Modified();
      c2->Update();
      gSystem->ProcessEvents();
      
      TF1* fitfunc = ( (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus")) );
      
      toffset_by_bar[bar] = fitfunc->GetParameter("Mean");
      dtoffset_by_bar[bar] = fitfunc->GetParError(fitfunc->GetParNumber("Mean"));

      xtest = toffset_by_bar[bar];
      
    } else {
      toffset_by_bar[bar] = 0.0;
      dtoffset_by_bar[bar] = 0.5;
    }
  }

  xtest = toffset_by_bar[refID];
  
  for( int bar=refID-1; bar>=0; bar-- ){
    TH1D *htemp = ( (TH1D*) dthodo_by_bar->At(bar) );
    htemp->Draw();

    if( htemp->GetEntries() >= 500.0 ){
      
      gPad->Modified();
      c2->Update();
      gSystem->ProcessEvents();
      
      //      double xtest = toffset_by_bar[bar+1];
      
      TSpectrum spectemp;
      int npeakstemp = spectemp.Search( htemp );
      double *xpeakstemp = spectemp.GetPositionX();
      ipeak_best = -1;
      mindiff = 1000.0;
      
      for( int ipeak=0; ipeak<npeakstemp; ipeak++ ){
	//int bin = 1 + int(peakpos[ipeak]+0.5);
	//double pos = htemp->GetBinCenter(bin);
	double pos = xpeakstemp[ipeak];
	if( ipeak_best < 0 || fabs( pos - xtest ) < mindiff ){
	  ipeak_best = ipeak;
	  mindiff = fabs( pos - xtest );
	  xpeak_best = pos;
	}
      }
      
      htemp->Fit("gaus","SQ","",xpeak_best-1.0,xpeak_best+1.0);
      gPad->Modified();
      c2->Update();
      gSystem->ProcessEvents();
      
      TF1* fitfunc = ( (TF1*) (htemp->GetListOfFunctions()->FindObject("gaus")) );
      
      toffset_by_bar[bar] = fitfunc->GetParameter("Mean");
      dtoffset_by_bar[bar] = fitfunc->GetParError(fitfunc->GetParNumber("Mean"));

      xtest = toffset_by_bar[bar];
      
    } else {
      toffset_by_bar[bar] = 0.0;
      dtoffset_by_bar[bar] = 0.5;
    }
  }


  double gbarID[90];
  double dgbarID[90];

  

  for( int bar=0; bar<90; bar++ ){
    gbarID[bar]=bar;
    dgbarID[bar] = 0.0;
  }

  TGraphErrors *gt0RF = new TGraphErrors( 90, gbarID, toffset_by_bar, dgbarID, dtoffset_by_bar );

  gt0RF->SetMarkerStyle(20);
  c1->cd();
  gt0RF->Draw("PSAME");

  gt0RF->Write("gT0RF");

  TString offsetsfilename = outputfilename;
  offsetsfilename.ReplaceAll(".root",".dat");

  ofstream offsetsfile(offsetsfilename.Data());

  //Let's write out 10 rows of 9 numbers each:
  int ibar=0;
  for( int irow=0; irow<10; irow++ ){
    for( int icol=0; icol<9; icol++ ){
      TString st0;
      offsetsfile << st0.Format("%12.6g",toffset_by_bar[ibar++]) << " ";
    }
    offsetsfile << endl;
  }
      
   
  
  // TF1 *fitfunc = new TF1("fitfunc","[0]*exp(-0.5*pow((x-[1])/[2],2))",-120,120);

  //c2->cd();
  
  // for( int ibar=0; ibar<90; ibar++ ){
  //   TH1D *htemp = hdT_hodo_ID->ProjectionY( histname.Format( "hdThodo_bar%d", ibar ), ibar+1, ibar+1 );

  //   c2->cd();
  //   htemp->Draw();

  //   gPad->Modified();
  //   c2->Update();
  //   gSystem->ProcessEvents();
    
  //   double Tmean = htemp->GetMean();

  //   cout << "bar " << ibar << ", nevents = " << htemp->GetEntries() << ", Tmean = " << Tmean << " ns" << endl;
    
  //   double startpar[3] = {20.0, Tmean, 0.4};

  //   fitfunc->SetParameters(startpar);

  //   //   fitfunc->SetParLimits(2,0.01,1.0);
  //   //fitfunc->SetParLimits(1,Tmean-2.004,Tmean+2.004);
    
    
  // }

  //  TGraphErrors *gT0RF = new TGraphErrors( IDRF.size(), &(IDRF[0]), &(T0RF[0]), &(dIDRF[0]), &(dT0RF[0]));
  //gT0RF->SetMarkerStyle(20);
  //gT0RF->SetMarkerColor(6);
  //gT0RF->SetLineColor(6);

  //c1->cd();
  //gT0RF->Draw("PSAME");

  hdT_hodo_ID->Write();
  //gT0RF->Write("gT0RF");
  
}
