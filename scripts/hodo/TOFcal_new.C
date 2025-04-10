#include "gmn_tree.C"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TCut.h"
#include "TTreeFormula.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include <cmath>
#include <vector>
#include "TCanvas.h"
#include "TPad.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TMinuit.h"

const double Mp=0.938272;

// long NEVENTSFIT;
// vector<int> NHITS_EVENT;
// vector<vector<int> > HITBAR; //BAR ID OF HITS
// vector<vector<double> > HIT_TLEFT; // "RAW" TLEFT
// vector<vector<double> > HIT_TRIGHT; // "RAW" TRIGHT
// vector<vector<double> > HIT_TOTLEFT; // TOT LEFT
// vector<vector<double> > 
// void CHI2_FCN( int &npar, double *gin, double &f, double *par, int flag ){

// }

void TOFcal_new(const char *inputfilename, const char *outputfilename="TOFcal_temp.root"){

  ifstream configfile(inputfilename);

  if( !configfile ) return;

  TFile *fout = new TFile(outputfilename,"RECREATE");
  
  TChain *C = new TChain("T");
  
  TString currentline;
  while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline.Data());
    }
  }

  TCut globalcut = "";
  
  while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline.Data();
    }
  }

  
  
  //Default reference PMTs for aligning all other PMTs:
  //NOTE: in the current procedure these definitions are meaningless; since we align all the corrected mean times to zero!
  int hodorefID = 43; // reference BAR (not PMT); we'll set the LEFT PMT to the reference by convention
  int hcalrefID = 199; //row 16 column 6

  double zhodo = 1.854454; //meters
  double Lbar_hodo = 0.6; //meters
  double etof0 = (1.96+3.0)/0.299792458; // "central" TOF value:

  double thetaHCAL = 29.4;
  double dHCAL = 11.0;
  double dSBS = 2.25;

  double W2min=0.6, W2max=1.2;
  
  double sbsmaxfield = 1.26;
  double sbsfield_default = 0.0; //default to zero-field assumption for proton TOF calculations
  // When combining different field settings for calculating proton
  // TOF, we need to associate field values with run ranges

  //Need some hit quality checks:
  double totmin = 6.0, totmax = 35.0;
  double totdiff_max = 10.0;
  double tdiff_max = 10.0; //max time difference with respect to cluster maximum
  double ydiff_max = 0.15; // 0.15 m/ (0.18 m/ns) 
  double xdiff_max = 0.06; // vertical position difference of 0.06 m max wrt track projection 
  
  
  int mineventsperbar = 100;
  
  int nfieldsettings;
  vector<int> rmin,rmax;
  vector<double> sbsfield;
  vector<double> ptof_central;
  
  const double Dgap = 48.0*2.54/100.0;
  
  double Ebeam = 5.9826;

  double RF_freq = 249.5e6;
  double bunch_spacing_ns = 1.0e9/RF_freq;
  
  cout << "bunch spacing = " << bunch_spacing_ns << " ns" << endl;

  TString fname_RFoffsets="";

  TString fname_meantime_offsets="";

  TString fname_vscint_old=""; 

  TString fname_walkcorr_old="";
  
  bool useRFcorr = false;
  bool fixoffsets = false;

  double vscint_default = 0.16; //m/ns
  double wL_default = 0.3;
  double wR_default = 0.21;
  
  while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endconfig") ){
    if( !currentline.BeginsWith("#") ){
      TObjArray *currentline_tokens = currentline.Tokenize(" ");
      if( currentline_tokens->GetEntries() >= 2 ){
	TString skey = ( (TObjString*) (*currentline_tokens)[0] )->GetString();

	TString sval = ( (TObjString*) (*currentline_tokens)[1] )->GetString();

	if( skey.BeginsWith("totmin") ){
	  totmin = sval.Atof();
	}

	if( skey.BeginsWith("totmax") ){
	  totmax = sval.Atof();
	}

	if( skey.BeginsWith("tdiff_max") ){
	  tdiff_max = sval.Atof();
	}

	if( skey.BeginsWith("ydiff_max") ){
	  ydiff_max = sval.Atof();
	}

	if( skey.BeginsWith("xdiff_max") ){
	  xdiff_max = sval.Atof();
	}

	if( skey.BeginsWith("totdiff_max") ){
	  totdiff_max = sval.Atof();
	}
	
	if( skey.BeginsWith("mineventsperbar") ){
	  mineventsperbar = sval.Atoi();
	}
	
	if( skey.BeginsWith("hodorefID") ){
	  hodorefID = sval.Atoi();
	}

	if( skey.BeginsWith("hcalrefID") ){
	  hcalrefID = sval.Atoi();
	}

	if( skey.BeginsWith("zhodo") ){
	  zhodo = sval.Atof();
	}

	if( skey.BeginsWith("Lbar_hodo") ){
	  Lbar_hodo = sval.Atof();
	}

	if( skey.BeginsWith("etof0") ){
	  etof0 = sval.Atof();
	}

	if( skey.BeginsWith("bunch_spacing_ns") ){
	  bunch_spacing_ns = sval.Atof();
	}

	if( skey.BeginsWith( "RFoffsets" ) ){
	  fname_RFoffsets = sval;
	}

	if( skey.BeginsWith( "userfcorr" ) ){
	  int val = sval.Atoi();
	  useRFcorr = val > 0;
	}

	if( skey.BeginsWith("ebeam") ){
	  Ebeam = sval.Atof();
	}

	if( skey.BeginsWith("sbsmaxfield") ){
	  sbsmaxfield = sval.Atof();
	}
	if( skey.BeginsWith("sbsdist") ){
	  dSBS = sval.Atof();
	}
	if( skey.BeginsWith("hcaldist") ){
	  dHCAL = sval.Atof();
	}
	if( skey.BeginsWith("thetaHCAL") ){
	  thetaHCAL = sval.Atof();
	}
	if( skey.BeginsWith("W2min") ){
	  W2min = sval.Atof();
	}
	if( skey.BeginsWith("W2max") ){
	  W2max = sval.Atof();
	}

	if( skey.BeginsWith("meantime_offsets_old") ){
	  fname_meantime_offsets = sval;
	}

	if( skey.BeginsWith("vscint_old") ){
	  fname_vscint_old = sval;
	}

	if( skey.BeginsWith("walkcorr_old") ){ //walk correction slopes
	  fname_walkcorr_old = sval;
	}
      }
    }
  }

  nfieldsettings = 0;
  
  //Parsing field configuration info to calculate proton TOF corrections:
  while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endfieldconfig") ){
    if( !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      //Here we expect three values per line, rmin, rmax, field (as fraction of maximum):
      if( tokens->GetEntries() >= 3 ){
	TString srmin = ( (TObjString*) (*tokens)[0] )->GetString();
	TString srmax = ( (TObjString*) (*tokens)[1] )->GetString();
	TString sfieldval = ( (TObjString*) (*tokens)[2] )->GetString();

	nfieldsettings++;
	rmin.push_back( srmin.Atoi() );
	rmax.push_back( srmax.Atoi() );
	sbsfield.push_back( sfieldval.Atof());
      }
    }
  }

  std::cout << "Number of field settings = " << nfieldsettings << endl;
  
  
  // Mean time offsets for start time:
  vector<double> RFoffsets(90,0.0);

  ifstream offsetsfile(fname_RFoffsets);
  if( offsetsfile ){
    int noffsets=0;
    while( offsetsfile >> RFoffsets[noffsets] ){
      cout << "bar, RF offset = " << noffsets << ", " << RFoffsets[noffsets] << endl;
      noffsets++;
    }
    cout << "Number of RF offsets loaded = " << noffsets << endl;
    
    if( noffsets != 90 ){
      cout << "Error: RF offsets file does not contain 90 entries, abort" << endl;
      useRFcorr = false;
      exit(-1);
    }
  }

  vector<double> MeanTimeOffsets(90,0.0);

  ifstream offsetsfile2(fname_meantime_offsets.Data());
  
  if( offsetsfile2 ){
    int noffsets = 0;
    while( offsetsfile2 >> MeanTimeOffsets[noffsets] ){
      cout << "bar, meantime offset = " << noffsets << ", " << MeanTimeOffsets[noffsets] << endl;
      noffsets++;
    }
    fixoffsets = true;
    cout << "Number of mean-time offsets loaded = " << noffsets << endl;
    if( noffsets != 90 ){
      cout << "Error: meantime offsets file does not contain 90 entries, abort" << endl;
      fixoffsets = false;
      exit(-1);
    }
  }
	   
  vector<double> VSCINT_OLD(90,vscint_default);
  vector<double> wL_OLD(90,wL_default);
  vector<double> wR_OLD(90,wR_default);

  ifstream vscint_file(fname_vscint_old.Data());
  if( vscint_file ){
    int npar=0;
    while( vscint_file >> VSCINT_OLD[npar] ){
      cout << "bar, vscint = " << npar << ", " << VSCINT_OLD[npar] << endl;
      npar++;
    }

    if( npar != 90 ){
      cout << "Error: vscint file does not contain 90 entries, abort" << endl;
      exit(-1);
    }
  }

  ifstream walkcorr_file(fname_walkcorr_old);

  if( walkcorr_file ){
    int npar=0;
    double par = 0.0;
    while( walkcorr_file >> par ){
      if( par < 0 ) par *= -1.0;
      if( npar<90 ) wL_OLD[npar] = par;
      if( npar >= 90 ) wR_OLD[npar-90] = par;
      npar++;
    }

    if( npar != 180 ){
      cout << "Error: walk correction file given does not contain 180 entries, abort" << endl;
      exit(-1);
    }
  }
  
  thetaHCAL *= TMath::Pi()/180.0;
  
  TVector3 zaxis_HCAL(-sin(thetaHCAL),0,cos(thetaHCAL));
  TVector3 xaxis_HCAL(0,-1,0);
  TVector3 yaxis_HCAL = (zaxis_HCAL.Cross(xaxis_HCAL)).Unit();

  TVector3 HCALorigin = dHCAL*zaxis_HCAL;

  double pp_ptheta_central = 2.0*Mp*Ebeam*(Mp+Ebeam)*cos(thetaHCAL)/(pow(Mp,2)+2.0*Mp*Ebeam + pow(Ebeam*sin(thetaHCAL),2));
  double Tp_ptheta_central = sqrt(pow(pp_ptheta_central,2)+pow(Mp,2))-Mp;

  double beta_ptheta_central = pp_ptheta_central/(Tp_ptheta_central+Mp);

  double ptof_central_default = dHCAL/(beta_ptheta_central*0.299792458);

  if( nfieldsettings > 0 ){
    for( int iset=0; iset<nfieldsettings; iset++ ){
      double BdL = sbsfield[iset]*sbsmaxfield*Dgap;
      double Radius = pp_ptheta_central / 0.299792458 / (sbsfield[iset]*sbsmaxfield);
      double x_in = 0.0;
      double z_in = dSBS;
      double x_center = x_in - Radius;
      double z_center = z_in;
      double z_out = dSBS + Dgap;
      double A = 1.0;
      double B = -2.0*x_center;
      double C = pow(x_center,2) + pow(z_out - z_center,2) - pow(Radius,2);
      double x_out = (-B + sqrt(pow(B,2)-4.0*A*C))/(2.0*A);
      double theta_in = 0.0;
      double theta_out = asin( (z_center-z_out)/Radius );
      //need to calculate path length:
      double pathlength = dSBS + Radius * (theta_in - theta_out) + (dHCAL-z_out )/cos(theta_out);
      ptof_central.push_back( pathlength/(beta_ptheta_central*0.299792458)); //ns 
    }
  }
  
  
  //int nparams_hodo = 3*180; //offset, walk correction, propagation speed for left and right PMTs for 90 bars minus offset fixed at zero for reference PMT
  //actually, vscint should be determined per paddle;
  int nparams_hodo = 180+2*90; //We actually only really want one zero offset (paddle-specific) that applies to the mean time, so we need 90 zero offsets, 90 vscint values, and 180 walk correction slopes.
  int nparams_HCAL = 288*2; //offset and two walk correction parameters for 288 modules.

  int nparams_hodo_t0fixed= 3*90; // wL, wR, vscint
  
  //Ordering of params:
  // t0, walkileft, walkiright, 1/vscint
  
  TMatrixD Mhodo(nparams_hodo,nparams_hodo);
  TMatrixD Mhcal(nparams_HCAL,nparams_HCAL);
  TVectorD bhodo(nparams_hodo);
  TVectorD bhcal(nparams_HCAL);

  //This set of equations will determine propagation speed and walk correction slope with fixed t0 offsets:
  TMatrixD Mhodo_t0fixed( nparams_hodo_t0fixed, nparams_hodo_t0fixed );
  TVectorD bhodo_t0fixed( nparams_hodo_t0fixed );
  
  for( int ipar=0; ipar<nparams_hodo; ipar++ ){
    for( int jpar=0; jpar<nparams_hodo; jpar++ ){
      Mhodo(ipar,jpar) = 0.0;
    }
    bhodo(ipar) = 0.0;
  }

  for( int ipar=0; ipar<nparams_hodo_t0fixed; ipar++ ){
    for( int jpar=0; jpar<nparams_hodo_t0fixed; jpar++ ){
      Mhodo_t0fixed(ipar,jpar) = 0.0;
    }
    bhodo_t0fixed(ipar) = 0.0;
  }

  for( int ipar=0; ipar<nparams_HCAL; ipar++ ){
    for( int jpar=0; jpar<nparams_HCAL; jpar++ ){
      Mhcal(ipar,jpar) = 0.0;
    }
    bhcal(ipar) = 0.0;
  }
  
  //HCAL and hodo TDCs both use BigBite trigger as reference time, so both have reference time subtracted:

  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut",globalcut,C);

  gmn_tree *T = new gmn_tree(C);

  C->SetBranchStatus("*",0);
  C->SetBranchStatus("g.*",1);
  C->SetBranchStatus("bb.hodotdc.*",1);
  C->SetBranchStatus("sbs.hcal.tdctimeblk",1);
  C->SetBranchStatus("sbs.hcal.eblk",1);
  C->SetBranchStatus("sbs.hcal.idblk",1);
  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus("sbs.hcal.clus_blk.*",1);
  C->SetBranchStatus("bb.tr.*",1);
  C->SetBranchStatus("bb.gem.track.*",1);
  C->SetBranchStatus("e.kine.*",1);
  C->SetBranchStatus("bb.tdctrig.*",1);
  C->SetBranchStatus("sbs.tdctrig.*",1);
  C->SetBranchStatus("bb.ps.e",1);
  C->SetBranchStatus("bb.etot_over_p",1);
  C->SetBranchStatus("bb.gem.trigtime",1);
  
  int treenum=-1,oldtreenum=-1;

  long nevent=0;

  vector<double> nevent_vs_bar(90,0.0);

  TH2D *htmeancorr_vs_ID_old = new TH2D("htmeancorr_vs_ID_old", "OLD; bar ID; corrected t_{mean} (ns)", 90, -0.5, 89.5, 300,-30,30);
  
  
  while( T->GetEntry(nevent) ){

    if( nevent % 10000 == 0 ){
      cout << "nevent = " << nevent << ", run number = " << T->g_runnum << endl;
    }
    
    treenum = C->GetTreeNumber();
    if( nevent == 0 || treenum != oldtreenum ){
      oldtreenum = treenum;
      GlobalCut->UpdateFormulaLeaves();
    }

    bool passed_global_cut = GlobalCut->EvalInstance(0) != 0;

    if( passed_global_cut ){ //do the things:

      //Grab reference time and RF time for calibration:
      bool goodRF = false, goodTRIG = false;

      //We'll use bb.tdctrig branches for now; later we'll use the new convenience variables:
      double RFtime, TrigTime;
      // for( int ihit=0; ihit<T->; ihit++ ){
      if( fabs( T->bb_tdctrig_rftime <= 1.e5 ) ){ // RF Time:
	goodRF = true;
	RFtime = T->bb_tdctrig_rftime;
      }
	
	//if( T->bb_tdctrig_tdcelemID[ihit] == 5 ){ //Trigger/reference time:
      if( fabs( T->bb_tdctrig_trigtime <= 1.e5 ) ){
	goodTRIG = true;
	TrigTime = T->bb_tdctrig_trigtime;
      }
      
      //grab needed track parameters:
      double vz = T->bb_tr_vz[0];
      double pathl = T->bb_tr_pathl[0];
      double etof = pathl/0.299792458; //electron TOF from vertex to hodo.

      //These are "raw" leading-edge times and TOT values for the max bar in the

      if( T->bb_hodotdc_clus_trackindex[0]==0 && goodRF && goodTRIG){

	//double Tstart = fmod( RFtime - TrigTime, bunch_spacing_ns );


	//Tstart = 0.0;
	//double Tstart = RFtime - TrigTime - bunch_spacing_ns * floor( (RFtime-TrigTime)/bunch_spacing_ns );
	
	// cout << "(Rftime,TrigTime,Tstart,bunchspace) = (" << RFtime << ", "
	//      << TrigTime << ", " << Tstart << ", " << bunch_spacing_ns << ") ns" << endl;
	
	double yhodo = T->bb_tr_y[0]+zhodo*T->bb_tr_ph[0];
	double xhodo = T->bb_tr_x[0]+zhodo*T->bb_tr_th[0];

	double dLEFT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 - yhodo));
	double dRIGHT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 + yhodo));

	double RFcorr = RFtime + vz/0.299792458;
	//	RFcorr -= bunch_spacing_ns * (floor(RFcorr/bunch_spacing_ns)+0.5);

	double tmean0 = 0.0;
	
	for( int ibar=0; ibar<int(T->bb_hodotdc_clus_size[0]); ibar++ ){
	    
	  // Now that we have a prescription for aligning all the RF-corrected
	  // hodo mean times into a common "bucket", the thing we want to minimize
	  // is to find the meantime offsets that align the following quantity at zero:
	  
	  
	  //tleft and tright are raw leading-edge times (ref time and TDC offset subtracted)
	  double tleft = T->bb_hodotdc_clus_bar_tdc_tleft[ibar];
	  double tright = T->bb_hodotdc_clus_bar_tdc_tright[ibar];
	  double totleft = T->bb_hodotdc_clus_bar_tdc_totleft[ibar];
	  double totright = T->bb_hodotdc_clus_bar_tdc_totright[ibar];
	  
	  //Some additional variables for hit quality checks:
	  double tmean = T->bb_hodotdc_clus_bar_tdc_meantime[ibar];
	  double ypos = T->bb_hodotdc_clus_bar_tdc_timehitpos[ibar];
	  double xpos = T->bb_hodotdc_clus_bar_tdc_vpos[ibar];
	  
	  if( ibar == 0 ) {
	    tmean0 = tmean;
	  }
	  
	  int ID = int(T->bb_hodotdc_clus_bar_tdc_id[ibar]);

	  RFcorr = RFtime + vz/0.299792458 + RFoffsets[ID];
	  
	  bool goodhit = totleft >= totmin && totleft <= totmax &&
	    totright >= totmin && totright <= totmax &&
	    fabs( totleft - totright ) <= totdiff_max &&
	    fabs( ypos - yhodo ) <= ydiff_max &&
	    fabs( xpos - xhodo ) <= xdiff_max &&
	    fabs( tmean - tmean0 ) <= tdiff_max;

	  if( goodhit ){
	    nevent_vs_bar[ID] += 1.;
	  
	    // how to define our chi2 statistic? We want to correct all hodo PMT times to tvertex = 0:
	    // How will the algebra change if we want to align the hodoscope mean times to the RF signal instead of the BigBite trigger signal?
	    // Basically we want to add back in the reference time to get the "raw" times of the PMTs, and then subtract the RF time modulo bunch spacing (or some
	    // integer multiple thereof:
	    // Corrected raw PMT time is Tcorr = tPMT + ref time - RF time
	    // 
	    // tPMT = etof + t0 - walk * TOT + d/vscint
	    // chi2 = sum_{i=1}^N \sum_{j=1}^Nhit (tPMT-etof-t0+walk*TOT - d/vscint)^2
	    // dchi2/dt0_k = 2*\sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*-1
	    // dchi2/dwalk_k = 2*sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*TOT
	    // dchi2/d(1/v)_k = 2*sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*-d
	
	  
	
	    // cout << "Best hodo hit (tL,tR, totL, totR, yhodo, dL, dR, pathl, etof)=("
	    // 	   << tleft << ", " << tright << ", " << totleft << ", " << totright << ", " << yhodo << ", "
	    // 	   << dLEFT << ", " << dRIGHT << ", " << pathl << ", " << etof << ")" << endl;
	
	    // we want to minimize
	    // chi2 = sum_{i=1}^Nevent (tleft_corr^2 + tright_corr^2)
	    // dchi2/dt0_j = sum_i (2*tleft_corr_i*dtleft_corr_i/dt0_j + 2*tright_corr_i*dtright_corr_i/dt0_j) = 0
	    //we COULD just use TMinuit but that's lazy.
	    //let parameters 0-89 be the t0L, 90-179 be t0R, 180-269 be wL, 270-359 be wR, and 360-449 be 1/vscint:
	    int ipar_t0 = ID;
	    int ipar_vinv = ID+90;
	    //int ipar_t0R = ID+90;
	    int ipar_wL = ID+180;
	    int ipar_wR = ID+270;

	    // double coeff_t0L = 
	
	    //sum_i (tL - etof - t0 + w*TOTL - dL/v) = Tstart  
	    // tL-etof = t0L -w*TOTL + dL/v

	    // tL = TRFcorr + etof + RFoffset + t0 - w*TOT - D/v

	    //Correct for TOF:
	    
	    double tleft_corr = tleft - (etof-etof0) + wL_OLD[ID]*totleft - dLEFT/VSCINT_OLD[ID] - MeanTimeOffsets[ID];
	    //Correct everything to a common beam bucket:
	    //This is the thing we want to align to zero!
	    //tleft_corr -= bunch_spacing_ns * (floor( tleft_corr/bunch_spacing_ns ) + 0.5);
	    double tright_corr = tright - (etof-etof0) + wR_OLD[ID]*totright - dRIGHT/VSCINT_OLD[ID] - MeanTimeOffsets[ID];

	    double tmean_corr = 0.5*(tleft_corr+tright_corr);

	    //tmean_corr contains all the old corrections

	    if( useRFcorr ){
	      tleft_corr -= RFcorr;
	      tright_corr -= RFcorr;
	      tmean_corr -= RFcorr;
	      
	      double RF_bucket = floor( tmean_corr/bunch_spacing_ns ); //largest integer smaller than tmeancorr/bunch spacing

	      // This RF bucket subtraction will only work if all other offsets/corrections/etc are reasonably well-known to begin with!
	      tmean_corr -= bunch_spacing_ns * (RF_bucket + 0.5);
	      tleft_corr -= bunch_spacing_ns * (RF_bucket + 0.5);
	      tright_corr -= bunch_spacing_ns * (RF_bucket + 0.5);
	    }
	    // if( fixoffsets ){
	    //   tleft_corr -= MeanTimeOffsets[ID];
	    //   tright_corr -= MeanTimeOffsets[ID];
	    // }
	    //tright_corr -= bunch_spacing_ns * (floor(tright_corr/bunch_spacing_ns) + 0.5);
	    //Correcting L and R PMT times separately to a common "beam bucket" does not work; one needs to correct the mean time:
	    // if( useRFcorr ){
	    //   tleft_corr -= RFcorr + RFoffsets[ID];
	    //   tright_corr -= RFcorr + RFoffsets[ID];
	    // }
	    
	    //  if( useRFcorr ){
	    //Commenting this until we can get a fit to converge with RF corrections included and resolve the "bucket ambiguity"

	    //Now that we've aligned all corrected times to a common RF bucket, let's REMOVE the walk and propagation delay corrections
	    // (but not the t0 offset). But we retain the TOF correction and the (old) mean time and RF offset:
	    double tleft_raw = tleft_corr - wL_OLD[ID]*totleft + dLEFT/VSCINT_OLD[ID];
	    double tright_raw = tright_corr - wR_OLD[ID]*totright + dRIGHT/VSCINT_OLD[ID];
	  
	    //there is a more elegant way to do this:
	    //Left PMT:
	    Mhodo(ipar_t0,ipar_t0) += 1.0;
	    Mhodo(ipar_t0,ipar_wL) += -totleft;
	    Mhodo(ipar_t0,ipar_vinv) += dLEFT;
	
	    Mhodo(ipar_wL, ipar_t0) += -totleft;
	    Mhodo(ipar_wL, ipar_wL) += pow(totleft,2);
	    Mhodo(ipar_wL, ipar_vinv) += -totleft*dLEFT;
	
	    Mhodo(ipar_vinv, ipar_t0) += dLEFT;
	    Mhodo(ipar_vinv, ipar_wL) += -totleft*dLEFT;
	    Mhodo(ipar_vinv, ipar_vinv) += pow(dLEFT,2);
	
	    bhodo(ipar_t0) += (tleft_raw)*1.0;
	    bhodo(ipar_wL) += (tleft_raw)*(-totleft);
	    bhodo(ipar_vinv) += (tleft_raw)*dLEFT;
	
	
	    //Right PMT:
	    Mhodo(ipar_t0,ipar_t0) += 1.0;
	    Mhodo(ipar_t0,ipar_wR) += -totright;
	    Mhodo(ipar_t0,ipar_vinv) += dRIGHT;
	
	    Mhodo(ipar_wR, ipar_t0) += -totright;
	    Mhodo(ipar_wR, ipar_wR) += pow(totright,2);
	    Mhodo(ipar_wR, ipar_vinv) += -totright*dRIGHT;
	    
	    Mhodo(ipar_vinv, ipar_t0) += dRIGHT;
	    Mhodo(ipar_vinv, ipar_wR) += -totright*dRIGHT;
	    Mhodo(ipar_vinv, ipar_vinv) += pow(dRIGHT,2);
	    
	    bhodo(ipar_t0) += (tright_raw)*1.0;
	    bhodo(ipar_wR) += (tright_raw)*(-totright);
	    bhodo(ipar_vinv) += (tright_raw)*dRIGHT;
	  }
	}
      }
      
    }
    nevent++;
  }

  //Before we fit, we want to fix all parameters for paddles with less than, say, 50 events (we'll replace the offset with zero and
  // the vscint and walk correction slopes with the average of all the other paddles):
  for( int i=0; i<90; i++ ){
    if( nevent_vs_bar[i] < double(mineventsperbar) ){
      int ipar_t0 = i;
      int ipar_vinv = i+90;
      int ipar_wL = i+180;
      int ipar_wR = i+270;

      Mhodo(ipar_t0, ipar_t0) = 1.;
      Mhodo(ipar_vinv, ipar_vinv) = 1.;
      Mhodo(ipar_wL, ipar_wL) = 1.;
      Mhodo(ipar_wR, ipar_wR) = 1.;
      for( int j=0; j<nparams_hodo; j++ ){
	if( j != ipar_t0 ) {
	  Mhodo(ipar_t0,j) = Mhodo(j,ipar_t0) = 0.0;
	}
	if( j != ipar_vinv ){
	  Mhodo(ipar_vinv, j) = Mhodo(j,ipar_vinv) = 0.0;
	}
	if( j != ipar_wL ){
	  Mhodo( ipar_wL, j ) = Mhodo(j,ipar_wL) = 0.0;
	}
	if( j != ipar_wR ){
	  Mhodo( ipar_wR, j ) = Mhodo(j,ipar_wL) = 0.0;
	}
      }
      bhodo( ipar_t0 ) = 0.0;
      bhodo( ipar_vinv ) = 0.0;
      bhodo( ipar_wL ) = 0.0;
      bhodo( ipar_wR ) = 0.0;
      
    }
  }
  
  TDecompSVD Ahodo(Mhodo);
  Ahodo.Solve(bhodo);

  bhodo.Print();

  double HODObarID[90];
  double HODOt0[90], HODOwL[90], HODOwR[90], vscint[90];

  for( int i=0; i<90; i++ ){
    HODObarID[i] = i;
    HODOt0[i] = bhodo[i];
    //HODOt0R[i] = bhodo[i+90];
    HODOwL[i] = bhodo[i+180];
    HODOwR[i] = bhodo[i+270];
    vscint[i] = 1.0/bhodo[i+90];
  }

  nevent = 0;

  treenum=-1,oldtreenum=-1;

  cout << "2nd loop, testing new hodoscope calibrations..." << endl;

  fout->cd();
  
  TH2D *htmean_vs_ID_old = new TH2D("htmean_vs_ID_old","OLD ;bar ID; bar t_{mean} (ns)", 90,-0.5,89.5, 300,-30,30);
  TH2D *htmean_vs_ID_new = new TH2D("htmean_vs_ID_new","NEW (with RF corrections);bar ID; bar t_{mean} (ns)", 90,-0.5,89.5, 300,-30,30);
  TH2D *htmean_vs_ID_new_noRFcorr = new TH2D("htmean_vs_ID_new_noRFcorr","NEW (no RF corrections);bar ID; bar t_{mean} (ns)",90,-0.5,89.5,300,-30,30);
  
  TH2D *htdiff_yhodo_old = new TH2D("htdiff_yhodo_old", "OLD; y_{hodo} (m); bar t_{L}-t_{R} (ns)", 150,-0.3,0.3,300,-15,15);
  TH2D *htdiff_yhodo_new = new TH2D("htdiff_yhodo_new", "NEW; y_{hodo} (m); bar t_{L}-t_{R} (ns)", 150,-0.3,0.3,300,-15,15);

  TH2D *htleft_corr_vs_ID_old = new TH2D("htleft_corr_vs_ID_old", "OLD; bar ID; t_{left}^{corr} (ns)",90,-0.5,89.5,300,-30,30);
  TH2D *htright_corr_vs_ID_old = new TH2D("htright_corr_vs_ID_old", "OLD; bar ID; t_{right}^{corr} (ns)",90,-0.5,89.5,300,-30,30);

  //Let's make some additional diagnostic plots on walk corrections (with all other corrections applied:
  TH2D *hwalkcorr_allpmts = new TH2D("hwalkcorr_allpmts", "All corrections except timewalk; Time-over-threshold (ns); T_{PMT}",200,0,60,300,-30,30);

  TClonesArray *hwalkcorr_by_PMT = new TClonesArray("TH2D",180);
  TClonesArray *htdiff_vs_y_by_bar = new TClonesArray("TH2D",90);
  
  for( int i=0; i<90; i++ ){
    TString histname;
    TString title;
    new( (*hwalkcorr_by_PMT)[i] ) TH2D(histname.Format("hwalkcorr_bar%d_left",i), title.Format("Bar %d Left: all corrections except timewalk; time-over-threshold (ns); corrected time (ns)",i), 200,0,60,300,-30,30);
    new( (*hwalkcorr_by_PMT)[i+90] ) TH2D(histname.Format("hwalkcorr_bar%d_right",i), title.Format("Bar %d Right: all corrections except timewalk; time-over-threshold (ns); corrected time (ns)",i), 200,0,60,300,-30,30);
    new( (*htdiff_vs_y_by_bar)[i] ) TH2D(histname.Format("htdiff_vs_y_bar%d",i), title.Format("Bar %d: all corrections except propagation delay; y_{hodo} (m); T_{L}-T_{R} (ns)",i), 150,-0.3,0.3,300,-15,15);
  }
  
  //Really want a large quantity of zero-field LH2 data for this purpose:
  TH2D *htdiffHCAL_vs_HCALID_new = new TH2D("htdiffHCAL_vs_HCALID_new", "NEW hodo, old HCAL; HCAL ID; t_{HCAL}-t_{HODO,corr} (ns)",289,-0.5,288.5,250,-25,25); 

  TH2D *htdiffHCAL_vs_eblk_all = new TH2D("htdiffHCAL_vs_eblk_all", "NEW hodo, old HCAL; energy deposit (GeV); #Deltat(ns)",250,0.0,1.0,250,-25,25);
  
  TH2D *hdxdy = new TH2D("hdxdy",";#Deltay (m); #Deltax (m)",250,-1.25,1.25,250,-1.25,1.25);

  vector<double> nevent_blk_HCAL(288,0.0);

  int bar3counter=0;
  
  while( T->GetEntry(nevent) ){

    if( nevent % 10000 == 0 ){
      cout << "nevent = " << nevent << ", run number = " << T->g_runnum << endl;
    }
    
    treenum = C->GetTreeNumber();
    if( nevent == 0 || treenum != oldtreenum ){
      oldtreenum = treenum;
      GlobalCut->UpdateFormulaLeaves();
    }

    bool passed_global_cut = GlobalCut->EvalInstance(0) != 0;

    if( passed_global_cut ){ //do the things:

      //Grab reference time and RF time for calibration:
      bool goodRF = false, goodTRIG = false;

      //We'll use bb.tdctrig branches for now; later we'll use the new convenience variables:
      double RFtime, TrigTime;
      // for( int ihit=0; ihit<T->; ihit++ ){
      if( fabs( T->bb_tdctrig_rftime <= 1.e5 ) ){ // RF Time:
	goodRF = true;
	RFtime = T->bb_tdctrig_rftime;
      }
	
	//if( T->bb_tdctrig_tdcelemID[ihit] == 5 ){ //Trigger/reference time:
      if( fabs( T->bb_tdctrig_trigtime <= 1.e5 ) ){
	goodTRIG = true;
	TrigTime = T->bb_tdctrig_trigtime;
      }

      double vz = T->bb_tr_vz[0];
      double pathl = T->bb_tr_pathl[0];
      double etof = pathl/0.299792458;

      double yhodo = T->bb_tr_y[0]+zhodo*T->bb_tr_ph[0];
	  
      // how to define our chi2 statistic? We want to correct all hodo PMT times to tvertex = 0:
      // tPMT = etof + t0 - walk * TOT + d/vscint
      // chi2 = sum_{i=1}^N \sum_{j=1}^Nhit (tPMT-etof-t0+walk*TOT - d/vscint)^2
      // dchi2/dt0_k = 2*\sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*-1
      // dchi2/dwalk_k = 2*sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*TOT
      // dchi2/d(1/v)_k = 2*sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*-d
      
      double dLEFT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 - yhodo));
      double dRIGHT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 + yhodo));
      
      if( T->bb_hodotdc_clus_trackindex[0]==0 && goodRF && goodTRIG ){

	double RFcorr = RFtime + vz/0.299792458;
	//RFcorr -= bunch_spacing_ns * (floor(RFcorr/bunch_spacing_ns)+0.5);
	
	//	double Tstart = fmod( RFtime - TrigTime, bunch_spacing_ns );

	double Tstart = 0.0;
	
	double Tmean_corr = -1000.;
	
	for( int ibar=0; ibar<int(T->bb_hodotdc_clus_size[0]); ibar++ ){
	  double tleft = T->bb_hodotdc_clus_bar_tdc_tleft[ibar];
	  double tright = T->bb_hodotdc_clus_bar_tdc_tright[ibar];
	  double totleft = T->bb_hodotdc_clus_bar_tdc_totleft[ibar];
	  double totright = T->bb_hodotdc_clus_bar_tdc_totright[ibar];
	  
	  int ID = int(T->bb_hodotdc_clus_bar_tdc_id[ibar]);

	  // if( ID == 3 ){
	  //   bar3counter++;
	  //   std::cout << "(size, ibar, ID, int(ID), nhits)=(" << int(T->bb_hodotdc_clus_size[0]) << ", " << ibar << ", " << T->bb_hodotdc_clus_bar_tdc_id[ibar]
	  // 	      << ", " << ID << ", " << bar3counter << ")" << std::endl;
	  // }
	  
	  double tleft_CORR = tleft - (etof-etof0) - HODOt0[ID] + HODOwL[ID]*totleft - dLEFT/vscint[ID] - RFcorr - RFoffsets[ID] - MeanTimeOffsets[ID];
	  double tright_CORR = tright - (etof-etof0) - HODOt0[ID] + HODOwR[ID]*totright - dRIGHT/vscint[ID] - RFcorr - RFoffsets[ID] - MeanTimeOffsets[ID];

	  // if( fixoffsets ){
	  //   tleft_CORR -= MeanTimeOffsets[ID];
	  //   tright_CORR -= MeanTimeOffsets[ID];
	  // }
	  
	  double tmean_CORR = 0.5*(tleft_CORR+tright_CORR);

	  double RF_bucket = floor( tmean_CORR/bunch_spacing_ns );

	  tmean_CORR -= bunch_spacing_ns * (RF_bucket + 0.5);
	  tleft_CORR -= bunch_spacing_ns * (RF_bucket + 0.5);
	  tright_CORR -= bunch_spacing_ns * (RF_bucket + 0.5);
	  
	  htmeancorr_vs_ID_old->Fill( ID, tmean_CORR + HODOt0[ID] );
	  
	  // double tmean_CORR_nowalk = 0.5*(tleft_CORR - HODOwL[ID]*totleft +
	  // 				  tright_CORR - HODOwR[ID]*totright);
	  
	  //tmean_CORR -= bunch_spacing_ns * (floor(tmean_CORR/bunch_spacing_ns)+0.5);
	  //	  tleft_CORR -= bunch_spacing_ns * (floor( tleft_CORR/bunch_spacing_ns )+0.5);
	  //tright_CORR -= bunch_spacing_ns * (floor( tright_CORR/bunch_spacing_ns )+0.5);
		  
	  //  tleft_CORR -= bunch_spacing_ns * ( floor(tmean_CORR/bunch_spacing_ns)+0.5 );
	  //tright_CORR -= bunch_spacing_ns * ( floor(tmean_CORR/bunch_spacing_ns)+0.5 );
	  
	  
	  double tmean_old = T->bb_hodotdc_clus_bar_tdc_meantime[ibar];
	  double tdiff_old = T->bb_hodotdc_clus_bar_tdc_timediff[ibar];
	  htleft_corr_vs_ID_old->Fill( ID, tleft_CORR );
	  htright_corr_vs_ID_old->Fill( ID, tright_CORR );
	  
	  hwalkcorr_allpmts->Fill( totleft, tleft_CORR - HODOwL[ID] * totleft );
	  hwalkcorr_allpmts->Fill( totright, tright_CORR - HODOwR[ID] * totright );

	  ( (TH2D*) (*hwalkcorr_by_PMT)[ID] )->Fill( totleft, tleft_CORR - HODOwL[ID]*totleft );
	  ( (TH2D*) (*hwalkcorr_by_PMT)[ID+90] )->Fill( totright, tright_CORR - HODOwR[ID]*totright );
	  
	  //calculate a corrected time difference without the propagation correction:
	  //This only includes walk correction and offset
	  //double tdiff_CORR = (tleft - HODOt0[ID]+ HODOwL[ID]*totleft)-
	  //  (tright-HODOt0[ID]+HODOwR[ID]*totright);

	  double tdiff_CORR = (tleft_CORR + dLEFT/vscint[ID]) -
	    (tright_CORR + dRIGHT/vscint[ID]);

	  ( (TH2D*) (*htdiff_vs_y_by_bar)[ID] )->Fill( yhodo, tdiff_CORR);
	  
	  
	  //double tmean_CORR = fmod(0.5*(tleft_CORR + tright_CORR),bunch_spacing_ns);
	  //tmean_CORR -= bunch_spacing_ns*(floor(tmean_CORR/bunch_spacing_ns)+0.5);
	  
	  if( ibar== 0 ) Tmean_corr = tmean_CORR;
	  
	  htmean_vs_ID_old->Fill( ID, tmean_old );
	  htmean_vs_ID_new->Fill( ID, tmean_CORR );
	  htmean_vs_ID_new_noRFcorr->Fill( ID, tmean_CORR + RFcorr + RFoffsets[ID] );
	  htdiff_yhodo_old->Fill( yhodo, tdiff_old );
	  htdiff_yhodo_new->Fill( yhodo, tdiff_CORR );
	}

	// HCAL TDC times have already been reference time subtracted;
	// But for comparison to the "RF-corrected" hodoscope mean time we want to remove the influence of the reference time subtraction from HCAL as well. 
	
	double RFtime_HCAL = T->sbs_tdctrig_rftime;
	double TrigTime_HCAL = T->sbs_tdctrig_trigtime;

	bool goodRF_HCAL = abs(RFtime_HCAL)<1.e5;
	bool goodTrig_HCAL = abs(TrigTime_HCAL)<1.e5;

	
	//if( goodRF_HCAL ){
	//HCAL RF time has already been reference subtracted
	double RFcorr_HCAL = RFtime_HCAL + vz/0.299792458;
	
	if( !goodRF_HCAL ){ //Replace HCAL RF time with BB RF time:
	  // Up to some integer multiple of RF spacing (160.321 ns), we have:
	  // RFtime HCAL - RFtime BB = Trigger Time HCAL - Trigger Time BB
	  // HOWEVER since RF time HCAL has already been reference subtracted, we
	  // don't need to add it back in
	  RFcorr_HCAL = RFtime - TrigTime;
	}
	RFcorr_HCAL -= bunch_spacing_ns*(floor(RFcorr_HCAL/bunch_spacing_ns)+0.5);
	
	double tHCAL = T->sbs_hcal_tdctimeblk - RFcorr_HCAL;
	double eblkHCAL = T->sbs_hcal_eblk;
	int idblkHCAL = int(T->sbs_hcal_idblk);
	double xHCAL = T->sbs_hcal_x;
	double yHCAL = T->sbs_hcal_y;
	
	//Eventually we'll need to develop parametrized TOF for protons/neutrons for HCAL for different field settings.
	
	TVector3 ep(T->bb_tr_px[0], T->bb_tr_py[0], T->bb_tr_pz[0]);
	TLorentzVector k(0,0,Ebeam,Ebeam);
	TLorentzVector P(0,0,0,Mp);
	
	TLorentzVector kprime(ep,ep.Mag());
	TLorentzVector q = k-kprime;
	TLorentzVector Pprime = P + q;
	
	double Q2 = q.M2();
	
	double etheta = ep.Theta();
	double ephi = ep.Phi();

	//We are ASSUMING that we're dealing with PROTONs here:
	double ptof_central_temp = ptof_central_default;
	double sbsfield_temp = sbsfield_default * sbsmaxfield;
	if( nfieldsettings>0 ){
	  for( int iset=0; iset<nfieldsettings; iset++ ){
	    if( T->g_runnum >= rmin[iset] && T->g_runnum <= rmax[iset] ){
	      sbsfield_temp = sbsmaxfield * sbsfield[iset];
	      ptof_central_temp = ptof_central[iset];
	    }
	  }
	}
	
	double Eprime_etheta = Ebeam/(1.0+Ebeam/Mp*(1.0-cos(etheta)));
	double Tp_etheta = Ebeam-Eprime_etheta;
	double pp_etheta = sqrt(pow(Tp_etheta,2)+2.0*Mp*Tp_etheta);
	double ptheta_etheta = acos( (Ebeam-Eprime_etheta*cos(etheta))/(pp_etheta));
	double pphi_ephi = ephi + TMath::Pi();
	
	
	double ptheta_4vect = q.Vect().Theta();
	double pphi_4vect = q.Vect().Phi();
	
	TVector3 pnhat_expect( sin(ptheta_etheta)*cos(pphi_ephi),
			       sin(ptheta_etheta)*sin(pphi_ephi),
			       cos(ptheta_etheta) );
	
	TVector3 vertex(0,0,vz);
	
	double sintersect = (HCALorigin-vertex).Dot(zaxis_HCAL)/(pnhat_expect.Dot(zaxis_HCAL));
	
	TVector3 HCAL_intersect = vertex + sintersect * pnhat_expect;
	
	double xHCAL_expect = (HCAL_intersect - HCALorigin).Dot( xaxis_HCAL );
	double yHCAL_expect = (HCAL_intersect - HCALorigin).Dot( yaxis_HCAL );
	
	double Lpath_HCAL_expect = (HCAL_intersect - vertex).Mag();
	
	double beta_proton = pp_etheta/sqrt(pow(pp_etheta,2)+pow(Mp,2));
	
	double TOF_HCAL_expect = Lpath_HCAL_expect/(beta_proton*0.299792458);

	double protondeflection = 0.0;
	
	if( sbsfield_temp != 0. ){ //Calculated expected TOF with correction for magnetic deflection:
	  double BdL = sbsfield_temp * Dgap;
	  double Radius = pp_etheta/0.299792458/sbsfield_temp;
	  double thtar_expect = atan( pnhat_expect.Dot( xaxis_HCAL ) / pnhat_expect.Dot( zaxis_HCAL ) );

	  double x_in = (dSBS - vertex.Dot(zaxis_HCAL))*tan(thtar_expect);
	  double z_in = dSBS;
	  double x_center = x_in - Radius * cos(thtar_expect);
	  double z_center = z_in + Radius * sin(thtar_expect);
	  double z_out = dSBS+Dgap;
	  double A = 1.0;
	  double B = -2.0*x_center;
	  double C = pow(x_center,2)+pow( z_out-z_center,2) - pow(Radius,2);

	  double x_out = (-B + sqrt(pow(B,2)-4.0*A*C))/(2.0*A);
	  double theta_out = asin( (z_center-z_out)/Radius );
	  double phtar_expect = atan( pnhat_expect.Dot( yaxis_HCAL )/pnhat_expect.Dot( zaxis_HCAL ) );
	  //note that this neglects
	  double pathL = (dSBS-vertex.Dot(zaxis_HCAL))*sqrt(1.0+pow(tan(thtar_expect),2)+pow(tan(phtar_expect),2)) + Radius * (thtar_expect - theta_out) + (dHCAL-z_out)*sqrt(1.0+pow(tan(theta_out),2)+pow(tan(phtar_expect),2));
	  TOF_HCAL_expect = pathL/(beta_proton*0.299792458);

	  double xproj_HCAL = x_out + tan(theta_out) * (dHCAL-z_out);
	  protondeflection = xHCAL_expect - xproj_HCAL;
	}

	//	cout << "proton deflection = " << protondeflection << endl;
	
	double W2 = T->e_kine_W2;
	
	// thetabend = 0.3/p * BdL ;
	
	//double BdL = Dgap * sbsmaxfield * sbsfield;
	//double thetabend = 0.3/pp_etheta * BdL;
	
	//	double protondeflection = tan(thetabend)*(dHCAL-(dSBS+Dgap/2.0)); 
	
	//int idblkHCAL = T->sbs_hcal_idblk;
	
	double deltax = xHCAL-(xHCAL_expect-protondeflection);
	double deltay = yHCAL-(yHCAL_expect);

	if( eblkHCAL > 0.02 && W2>W2min && W2<W2max ){
	  hdxdy->Fill( deltay,deltax );
	}
	
	//if( W2min < W2 && W2 < W2max && sqrt(pow(deltax,2)+pow(deltay,2))<=0.24 && eblkHCAL>0.02 ){
	if( W2>W2min && W2<W2max && sqrt(pow(deltax,2)+pow(deltay,2))<=0.24 ){

	  for( int iblk=0; iblk<T->Ndata_sbs_hcal_clus_blk_tdctime; iblk++ ){
	    if( T->sbs_hcal_clus_blk_tdctime[iblk] != -1000. && T->sbs_hcal_clus_blk_e[iblk] >= 0.025 ){
	      double eblk = T->sbs_hcal_clus_blk_e[iblk];
	      double tblk = T->sbs_hcal_clus_blk_tdctime[iblk] - RFcorr_HCAL;
	      double idblk = T->sbs_hcal_clus_blk_id[iblk];

	      double tblk_walkcorr = tblk - 1.25/sqrt(eblk) + 0.0285/eblk;
	      
	      htdiffHCAL_vs_eblk_all->Fill( eblk, tblk-Tmean_corr-(TOF_HCAL_expect - ptof_central_temp));
	  
	      htdiffHCAL_vs_HCALID_new->Fill( idblk, tblk_walkcorr-Tmean_corr-(TOF_HCAL_expect - ptof_central_temp));
	      int ipar_hcalt0 = idblk-1;
	      int ipar_hcalw0 = idblk+288-1;
	      //	      int ipar_hcalw1 = idblk+2*288-1;

	      double tblk_corr = tblk-Tmean_corr-(TOF_HCAL_expect-ptof_central_temp);
	      
	      // tHCAL = tstart + t0 + w0/sqrt(E)
	      Mhcal(ipar_hcalt0,ipar_hcalt0) += 1.0;
	      Mhcal(ipar_hcalt0,ipar_hcalw0) += 1.0/sqrt(eblk);
	      //Mhcal(ipar_hcalt0,ipar_hcalw1) += -1.0/eblk;

	      Mhcal(ipar_hcalw0,ipar_hcalt0) += 1.0/sqrt(eblk);
	      Mhcal(ipar_hcalw0,ipar_hcalw0) += 1.0/eblk;
	      //Mhcal(ipar_hcalw0,ipar_hcalw1) += pow(eblk,-1.5);

	      // Mhcal(ipar_hcalw1,ipar_hcalt0) += -1.0/eblk;
	      // Mhcal(ipar_hcalw1,ipar_hcalw0) += pow(eblk,-1.5);
	      // Mhcal(ipar_hcalw1,ipar_hcalw1) += pow(eblk,-2);

	      bhcal(ipar_hcalt0) += (tblk_corr)*1.0;
	      bhcal(ipar_hcalw0) += (tblk_corr)/sqrt(eblk);
	      //    bhcal(ipar_hcalw1) += -(tblk_corr)/eblk;

	      nevent_blk_HCAL[idblk-1]+=1.0;
	      
	    }
	  }
	  //}
	}
      }
    }
    nevent++;
  }

  double ngoodbars=0;
  double vscint_avg_good = 0.0;
  double wL_avg_good = 0.0;
  double wR_avg_good = 0.0;
  
  for( int i=0; i<90; i++ ){
    if( nevent_vs_bar[i] >= double(mineventsperbar) ){
      ngoodbars += nevent_vs_bar[i];
      vscint_avg_good += vscint[i]*nevent_vs_bar[i];
      wL_avg_good += HODOwL[i]*nevent_vs_bar[i];
      wR_avg_good += HODOwR[i]*nevent_vs_bar[i];
    }
  }

  vscint_avg_good /= ngoodbars;
  wL_avg_good /= ngoodbars;
  wR_avg_good /= ngoodbars;
  
  for( int i=0; i<90; i++ ){
    if( nevent_vs_bar[i] < mineventsperbar ){
      HODOt0[i] = 0.0;
      vscint[i] = vscint_avg_good;
      HODOwL[i] = wL_avg_good;
      HODOwR[i] = wR_avg_good;
    }
  }
  
  //bhodo.Print();

  ofstream wLfile("hodowL_temp.dat");
  ofstream wRfile("hodowR_temp.dat");
  ofstream vscintfile("hodovscint_temp.dat");
  ofstream t0file("hodot0_temp.dat");
  
  TString dbfilename(outputfilename);
  dbfilename.ReplaceAll(".root",".dat");
  ofstream dbfile(dbfilename.Data());

  dbfile << "bb.hodotdc.timewalk0map = " << endl;
  for( int i=0; i<90; i++ ){
    if( i>0 && i%9 == 0 ) {
      dbfile << endl;
      wLfile << endl;
    }
    TString entry;
    dbfile << entry.Format("%15.6g", -HODOwL[i]) << " ";
    wLfile << entry.Format("%15.6g", HODOwL[i]) << " ";
  }

  dbfile << endl;

  //dbfile << "bb.hodotdc.timewalk1map = " << endl;
  for( int i=0; i<90; i++ ){
    if( i>0 && i%9 == 0 ) {
      dbfile << endl;
      wRfile << endl;
    }
    TString entry;
    dbfile << entry.Format("%15.6g", -HODOwR[i]) << " ";
    wRfile << entry.Format("%15.6g", HODOwR[i]) << " ";
  }

  dbfile << endl << endl;
  dbfile << "bb.hodotdc.timewalk1map = " << endl;

  for( int i=0; i<180; i++ ){
    if( i>0 && i%9 == 0 ) dbfile << endl;
    TString entry;
    dbfile << entry.Format("%15.6g", 0.0) << " ";
  }
  dbfile << endl << endl;

  dbfile << "bb.hodotdc.vscint = " << endl;
  for( int i=0; i<90; i++ ){
    if( i>0 && i%9 == 0 ){
      dbfile << endl;
      vscintfile << endl;
    }
    TString entry;
    dbfile << entry.Format("%15.6g", vscint[i] ) << " ";
    vscintfile << entry.Format("%15.6g", vscint[i] ) << " ";
  }
  dbfile << endl << endl;

  // In the calibration procedure, we have aligned the following quantity at zero:
  // t corrected = 0 = t - (etof-etof0) - t0 + w*TOT - d/v - (RFtime+vz/c + RFoffset) - MeanTimeOffset 
  dbfile << "bb.hodotdc.meantime_offset = " << endl;
  for( int i=0; i<90; i++ ){
    if( i>0 && i%9 == 0 ){
      dbfile << endl;
      t0file << endl;
    }
    TString entry;
    dbfile << entry.Format("%15.6g", HODOt0[i]+MeanTimeOffsets[i]) << " ";
    t0file << entry.Format("%15.6g", HODOt0[i]+MeanTimeOffsets[i]) << " ";
  }
  dbfile << endl << endl;
  dbfile << "bb.hodotdc.rf_offset = " << endl;
  for( int i=0; i<90; i++ ){
    if( i>0 && i%9 == 0 ){
      dbfile << endl;
    }
    TString entry;
    dbfile << entry.Format("%15.6g", RFoffsets[i]) << " ";
  }
  dbfile << endl << endl;

  dbfile << "bb.hodotdc.tdiffoffset = " << endl;
  for( int i=0; i<90; i++ ){
    if( i>0 && i%9 == 0 ){
      dbfile << endl;
    }
    TString entry;
    dbfile << entry.Format("%15.6g", 0.0) << " ";
  }
  dbfile << endl;
  
  // Print out some diagnostic plots:
  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  

  TString pdffilename = outputfilename;
  pdffilename.ReplaceAll(".root",".pdf");

  htmean_vs_ID_new->Draw("colz");
  gPad->Modified();
  c1->Update();
  gSystem->ProcessEvents();
  TString tempname = pdffilename;
  tempname+= "(";
  c1->Print(tempname.Data());

  
  
  //  double lmargin=0.2, bmargin=0.12, rmargin=0.03, tmargin=0.04;

  c1->Clear();
  c1->Divide(3,3,.001,.001);

  for( int i=0; i<90; i++ ){
    c1->cd(i%9+1);
    TString histname,histtitle;
    TH1D *htemp = htmean_vs_ID_new->ProjectionY(histname.Format("htmean_corr_bar%d",i),i+1,i+1);
    htemp->SetTitle(histtitle.Format("Bar %d, all corrections; t_{mean} (ns);",i));
    htemp->Draw();
    gPad->Modified();
    c1->Update();
    gSystem->ProcessEvents();
    tempname = pdffilename;
    if( i%9 == 8 ) c1->Print(tempname.Data());
  }
  
  for( int i=0; i<90; i++ ){
    c1->cd(i%9+1);
    ( (TH2D*) (*htdiff_vs_y_by_bar)[i] )->Draw("colz");
    gPad->Modified();
    c1->Update();
    gSystem->ProcessEvents();

    if( i%9 == 8 ){ //Print this page:
      tempname=pdffilename;

      c1->Print(tempname.Data());
    }
  }

  //Walk for left pmts:
  for( int i=0; i<90; i++ ){
    c1->cd(i%9+1);
    ( (TH2D*) (*hwalkcorr_by_PMT)[i] )->Draw("colz");
    gPad->Modified();
    c1->Update();
    gSystem->ProcessEvents();

    if( i%9 == 8 ){ //Print this page:
      tempname=pdffilename;
      //if( i == 8 ) tempname += "(";

      c1->Print(tempname.Data());
    }
  }
  //walk for right pmts:
  for( int i=0; i<90; i++ ){
    c1->cd(i%9+1);
    ( (TH2D*) (*hwalkcorr_by_PMT)[i+90] )->Draw("colz");
    gPad->Modified();
    c1->Update();
    gSystem->ProcessEvents();

    if( i%9 == 8 ){ //Print this page:
      tempname=pdffilename;
      //if( i == 8 ) tempname += "(";
      //if( i == 89 ) tempname += ")";
      c1->Print(tempname.Data());
    }
  }

  //Do we now want to solve HCAL and do a THIRD loop over the data?

  for( int i=0; i<288; i++ ){
    if( nevent_blk_HCAL[i] < 100. ){
      int ipar_t0 = i;
      int ipar_w0 = i+288;
      //      int ipar_w1 = i+288;
      Mhcal(ipar_t0,ipar_t0)=1.;
      Mhcal(ipar_w0,ipar_w0)=1.;
      //Mhcal(ipar_w1,ipar_w1)=1.;
      bhcal(ipar_t0) = 0.;
      bhcal(ipar_w0) = 0.;
      //bhcal(ipar_w1) = 0.;
      for( int j=0; j<nparams_HCAL; j++){ //zero out all off-diagonal elements:
	if( j != ipar_t0 ){
	  Mhcal(ipar_t0,j) = Mhcal(j,ipar_t0) = 0.0;
	}
	if( j != ipar_w0 ){
	  Mhcal(ipar_w0,j) = Mhcal(j,ipar_w0) = 0.0;
	}
	// if( j != ipar_w1 ){
	//   Mhcal(ipar_w1,j) = Mhcal(j,ipar_w1) = 0.0;
	// }
      }
    }
  }

  TDecompSVD Ahcal(Mhcal);
  Ahcal.Solve(bhcal);
  bhcal.Print();

  double t0avg_hcal = 0.0;
  double walkcorr_avg_hcal = 0.0;
  double nblocks_hcal = 0.0;

  double t0hcal[288];
  double w0hcal[288];
  double blkIDhcal[288];
  for( int i=0; i<288; i++ ){
    if( nevent_blk_HCAL[i] >= 100. ){
      walkcorr_avg_hcal += nevent_blk_HCAL[i]*bhcal[i+288];
      nblocks_hcal += nevent_blk_HCAL[i];
      t0avg_hcal += nevent_blk_HCAL[i]*bhcal[i];
      t0hcal[i] = bhcal[i];
      w0hcal[i] = bhcal[i+288];
      
    }
    blkIDhcal[i]=i+1;
  }

  walkcorr_avg_hcal /= nblocks_hcal;
  t0avg_hcal /= nblocks_hcal;
  
  for( int i=0; i<288; i++ ){
    if( nevent_blk_HCAL[i] < 100. ){
      t0hcal[i] = t0avg_hcal;
      w0hcal[i] = walkcorr_avg_hcal;
    }
  }
  
  nevent=0;

  TGraph *gHODOt0 = new TGraph(90,HODObarID,HODOt0);
  TGraph *gHODOwL = new TGraph(90,HODObarID,HODOwL);
  TGraph *gHODOwR = new TGraph(90,HODObarID,HODOwR);
  TGraph *gHODOvscint = new TGraph(90,HODObarID,vscint);

  TGraph *gHCALt0 = new TGraph(288,blkIDhcal,t0hcal);
  TGraph *gHCALw0 = new TGraph(288,blkIDhcal,w0hcal);
  
  gHODOt0->SetMarkerStyle(20);
  gHODOwL->SetMarkerStyle(20);
  gHODOwR->SetMarkerStyle(21);
  gHODOvscint->SetMarkerStyle(20);
  gHCALt0->SetMarkerStyle(20);
  gHCALw0->SetMarkerStyle(20);
  
  c1->Clear();
  c1->Divide(2,2,.001,.001);
  c1->cd(1);
  gHODOt0->Draw("AP");
  gHODOt0->SetTitle("Hodoscope fit results;bar ID;mean time offset (ns)");
  c1->cd(2);
  gHODOvscint->Draw("AP");
  gHODOvscint->SetTitle("Hodoscope fit results;bar ID;bar propagation speed (m/ns)");
  c1->cd(3);
  gHODOwL->Draw("AP");
  gHODOwL->SetTitle("Hodoscope fit results; bar ID; Left PMT walk correction slope");

  c1->cd(4);
  gHODOwR->Draw("AP");
  gHODOwR->SetTitle("Hodoscope fit results; bar ID; right PMT walk correction slope");

  c1->Update();
  gPad->Modified();
  gSystem->ProcessEvents();
  c1->Print(pdffilename);

  c1->Clear();
  c1->Divide(2,1,0.001,0.001);
  c1->cd(1);
  gHCALt0->Draw("AP");
  gHCALt0->SetTitle("HCAL fit results; Block ID; time offset (ns)");
  c1->cd(2);
  gHCALw0->Draw("AP");
  gHCALw0->SetTitle("HCAL fit results; Block ID; walk correction (1/sqrt(e)) slope");

  c1->Update();
  gPad->Modified();
  gSystem->ProcessEvents();
  pdffilename += ")";
  c1->Print(pdffilename);

  
  
  


  gHODOt0->Write("gHODOt0");
  gHODOwL->Write("gHODOwL");
  gHODOwR->Write("gHODOwR");
  gHODOvscint->Write("gHODOvscint");
  gHCALt0->Write("gHCALt0");
  gHCALw0->Write("gHCALw0");
  
  fout->Write();
  
}
