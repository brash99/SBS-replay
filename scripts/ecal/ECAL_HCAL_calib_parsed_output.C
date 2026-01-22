#include "TTree.h"
#include "TChain.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "gep_tree_parsed.C"
#include <vector>
#include <set>
#include <map>
#include <vector>
#include "TCut.h"
#include "TTreeFormula.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TSystem.h"
#include "TLine.h"

using namespace std;

const int nchan_ECAL = 1656;
const int nchan_HCAL = 288;

const int nrows_HCAL = 24;
const int ncols_HCAL = 12;

const int nrows = 69;
const int ncols[nrows] =
  { 12, 12, 12,
    18, 18, 18,
    21, 21, 21, 21, 21, 21,
    24, 24, 24, 24, 24, 24,
    27, 27, 27, 27, 27, 27,
    27, 27, 27, 27, 27, 27,
    27, 27, 27, 27, 27, 27,
    27, 27, 27, 27, 27, 27,
    27, 27, 27, 27, 27, 27,
    27, 27, 27, 27, 27, 27,
    24, 24, 24, 24, 24, 24,
    21, 21, 21, 21, 21, 21,
    18, 18, 18 };
    
    


void ECAL_HCAL_calib_parsed_output(const char *configfilename, const char *outfilename="GEP_ECAL_HCAL_calib_temp.root"){
  // What inputs do we need?
  //Main thing is we want a list of "bad/dead" channels. To start with we'll simply ignore the effects of dead channels until we have a proof-of-concept:
  // We need to define the list of files, the list of global cut and/or cut(s), minimum number of events per channel, etc.   

  TFile *fout = new TFile(outfilename,"RECREATE");
  
  ifstream configfile(configfilename);

  if( !configfile ) return;

  TString currentline;

  TChain *C = new TChain("Tout");
  
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

  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  c1->Divide(2,1,.001,.001);
  
  int minhits_ECAL = 100;
  int minhits_HCAL = 500; 

  int niter = 3;
  double minMoverB = 0.01;
  double mincoeff_ECAL = 0.001;
  double maxcoeff_ECAL = 0.1;

  double mincoeff_HCAL = 0.003;
  double maxcoeff_HCAL = 0.3;
    
  //We'll also need to load the "xpos" and "ypos" values of the
  //bad channels (might as well define for all the channels) to exclude them effectively

  vector<double> xposECAL(nchan_ECAL);
  vector<double> yposECAL(nchan_ECAL);
  vector<double> defaultGainECAL(nchan_ECAL,0.01); //to be loaded from a file
  vector<double> defaultGainHCAL(nchan_HCAL,0.02);
  
  bool gotoldgain_ECAL = false;

  bool gotoldgain_HCAL = false; 
  
  TString oldgainfname="";
  TString fnamexpos="";
  TString fnameypos="";

  TString oldgainfname_HCAL="";
  
  set<int> badchanlist_ECAL;
  set<int> badchanlist_HCAL; 

  double dtcut_ECAL = 10.0; //Max time difference from primary
  double dtcut_HCAL = 10.0; //Max time difference from primary

  double Rmax_ECAL = 0.15; //cm (max distance from seed to include in the fit)
  double Rmax_HCAL = 0.5; //cm (max distance from seed to include in the fit)

  double hit_thresh_ECAL = 0.01; //20 MeV
  double seed_thresh_ECAL = 0.1; // 100 MeV
  double clust_thresh_ECAL = 0.5; // 500 MeV

  // These are in GeV
  double hit_thresh_HCAL = 0.02;
  double seed_thresh_HCAL = 0.1;
  double clust_thresh_HCAL = 0.22; 

  double thetaFPPmaxHCAL = 0.015;
  
  while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endconfig") ){
    if( !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" ");

      if( tokens->GetEntries() >= 2 ){
	TString skey = ( (TObjString*) (*tokens)[0] )->GetString();

	TString sval = ( (TObjString*) (*tokens)[1] )->GetString();

	if( skey == "niter" ){
	  niter = sval.Atoi();
	}

 	if( skey == "mincoeff_ECAL" ){
	  mincoeff_ECAL = sval.Atof();
	}

	if( skey == "maxcoeff_ECAL" ){
	  maxcoeff_ECAL = sval.Atof();
	}

	if( skey == "mincoeff_HCAL" ){
	  mincoeff_HCAL = sval.Atof();
	}

	if( skey == "maxcoeff_HCAL" ){
	  maxcoeff_HCAL = sval.Atof();
	}
	
	if( skey == "minhits_ECAL" ){
	  minhits_ECAL = sval.Atoi();
	}

	if( skey == "minhits_HCAL" ){
	  minhits_HCAL = sval.Atoi();
	}

	if( skey == "dtcut_ECAL" ){
	  dtcut_ECAL = sval.Atof();
	}
	
	if( skey == "dtcut_HCAL" ){
	  dtcut_HCAL = sval.Atof();
	}
	if( skey == "Rmax_ECAL" ){
	  Rmax_ECAL = sval.Atof();
	}
	if( skey == "Rmax_HCAL" ){
	  Rmax_HCAL = sval.Atof();
	}

	if( skey == "badchanlist_ECAL" ){
	  for( int i=1; i<tokens->GetEntries(); i++ ){
	    TString schan = ( (TObjString*) (*tokens)[i] )->GetString(); 
	    badchanlist_ECAL.insert( schan.Atoi() );
	  }
	}

	if( skey == "badchanlist_HCAL" ){
	  for( int i=1; i<tokens->GetEntries(); i++ ){
	    TString schan = ( (TObjString*) (*tokens)[i] )->GetString(); 
	    badchanlist_HCAL.insert( schan.Atoi() );
	  }
	}
	//TODO: make hit/seed/cluster thresholds configurable
	if( skey == "oldgaincoeffs_ECAL" ){
	  oldgainfname = sval;

	  ifstream oldgainfile(oldgainfname.Data());

	  if( oldgainfile ){
	    int ncoeffsread = 0;
	    double coefftemp;
	    while( oldgainfile >> coefftemp && ncoeffsread < nchan_ECAL ){
	      defaultGainECAL[ncoeffsread] = coefftemp;
	      ncoeffsread++;
	    }
	    if( ncoeffsread == nchan_ECAL ) {
	      gotoldgain_ECAL = true;
	      cout << "successfully loaded old gain coeffs from file " << oldgainfname << endl;
	    }
	  }
	  
	}

	if( skey == "oldgaincoeffs_HCAL" ){
	  oldgainfname_HCAL = sval;
	  ifstream oldgainfile_HCAL( oldgainfname_HCAL.Data() );
	  if( oldgainfile_HCAL ){
	    int ncoeffsread = 0;
	    double coefftemp;
	    while( oldgainfile_HCAL >> coefftemp && ncoeffsread < nchan_HCAL ){
	      defaultGainHCAL[ncoeffsread] = coefftemp;
	      ncoeffsread++;
	    }

	    if( ncoeffsread == nchan_HCAL ){
	      gotoldgain_HCAL = true;
	      cout << "successfully loaded HCAL gain coeffs from file " << oldgainfname_HCAL << endl;
	    }
	  }
	}

	if( skey == "thetaFPPmaxHCAL" ){
	  thetaFPPmaxHCAL = sval.Atof();
	}
	
      }
    }
  }

  vector<int> Nhit_chan_ECAL(nchan_ECAL,0);
  vector<int> Nhit_chan_HCAL(nchan_HCAL,0);
  
  
  TMatrixD MECAL(nchan_ECAL,nchan_ECAL);
  TVectorD bECAL(nchan_ECAL);

  vector<bool> goodcalibECAL(nchan_ECAL,false);
  vector<double> PreviousSolutionECAL(nchan_ECAL);
  vector<double> sumEmiss_ratio(nchan_ECAL,0.0);
  vector<double> neventsEmiss_ratio(nchan_ECAL,0.0);

  set<int> list_uncalibrated; 
  
  for( int i=0; i<nchan_ECAL; i++ ){
    bECAL(i) = 0.0;
    //goodcalibECAL[i] = false;
    for( int j=0; j<nchan_ECAL; j++ ){
      MECAL(i,j) = 0.0;
    }
  }
  
  TMatrixD MHCAL(nchan_HCAL,nchan_HCAL);
  TVectorD bHCAL(nchan_HCAL);

  vector<bool> goodcalibHCAL(nchan_HCAL,false);
  vector<double> PreviousSolutionHCAL(nchan_HCAL);
  vector<double> sumEmiss_ratio_HCAL(nchan_HCAL,0.0);
  vector<double> neventsEmiss_ratio_HCAL(nchan_HCAL,0.0);

  set<int> list_uncalibrated_HCAL;
  
  for( int i=0; i<nchan_HCAL; i++ ){
    bHCAL(i) = 0.0;
    for( int j=0; j<nchan_HCAL; j++ ){
      MHCAL(i,j) = 0.0;
    }
  }
  
  
  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut, C );

  //What branches do we need to activate? Let's copy-paste from the "ElasticParse" macro:

  gep_tree_parsed *T = new gep_tree_parsed(C);

  long nevent=0;

  int treenum=-1, oldtreenum = -1;

  fout->cd();
  
  //Plots to make: OLD E/p
  TH1D *hEoverP_old_ECAL = new TH1D("hEoverP_old_ECAL","OLD;ECAL E/p; ", 200,0,2.0);
  TH2D *hEoverP_x_old_ECAL = new TH2D("hEoverP_x_old_ECAL", "OLD;xECAL (m);E/p",160,-1.6,1.6,200,0,2);
  TH2D *hEoverP_y_old_ECAL = new TH2D("hEoverP_y_old_ECAL", "OLD;yECAL (m);E/p",160,-0.8,0.8,200,0,2);

  TH1D *hEoverP_temp_ECAL = new TH1D("hEoverP_temp_ECAL","E/p",500,0,5);

  TH1D *hsampfrac_old_HCAL = new TH1D("hsampfrac_old_HCAL", "OLD; HCAL E/T_{p};", 500,0.0,5.0);
  TH2D *hsampfrac_x_old_HCAL = new TH2D("hsampfrac_x_old_HCAL", "OLD; xHCAL (m); HCAL E/T_{p}",200,-2.75,1.25,125,0,5);
  TH2D *hsampfrac_y_old_HCAL = new TH2D("hsampfrac_y_old_HCAL", "OLD; yHCAL (m); HCAL E/T_{p}",180,-0.9,0.9,125,0,5);

  TH1D *hSampFrac_temp_HCAL = new TH1D("hSampFrac_temp_HCAL", "E_{HCAL}/T_{p}", 250,0,5);
  
  double SampFracMeanHCAL = 1.0;
  double SampFracSigmaHCAL = 0.4;

  double SampFracMinHCAL = 0.0;
  double SampFracMaxHCAL = 3.0;
  
  double meanEoverP = 1.0;
  double sigmaEoverP = 0.13;
  
  for( int iter=0; iter<niter; iter++ ){
    cout << "ECAL/HCAL calib iteration " << iter << ": " <<  endl;

    nevent = 0;
    treenum = -1;
    oldtreenum = -1;

    hEoverP_temp_ECAL->SetTitle(Form("ECAL iteration %d; E/p;",iter));
    hEoverP_temp_ECAL->Reset();

    hSampFrac_temp_HCAL->SetTitle(Form("HCAL iteration %d; E/T_{p};",iter));
    hSampFrac_temp_HCAL->Reset();
    
    
    for( int i=0; i<nchan_ECAL; i++ ){
      bECAL(i) = 0.0;
      Nhit_chan_ECAL[i] = 0;
      for( int j=0; j<nchan_ECAL; j++ ){
	MECAL(i,j) = 0.0;
      }
      sumEmiss_ratio[i] = 0.0;
      neventsEmiss_ratio[i] = 0.0;
    }

    for( int i=0; i<nchan_HCAL; i++ ){
      bHCAL(i) = 0.0;
      Nhit_chan_HCAL[i] = 0;
      for( int j=0; j<nchan_HCAL; j++ ){
	MHCAL(i,j) = 0.0;
      }
      sumEmiss_ratio_HCAL[i] = 0.0;
      neventsEmiss_ratio_HCAL[i] = 0.0;
    }
    
    while( C->GetEntry(nevent++) ){
      treenum = C->GetTreeNumber();
      
      if( treenum != oldtreenum ){
	oldtreenum = treenum;
	GlobalCut->UpdateFormulaLeaves();
      }

      bool passedcut = GlobalCut->EvalInstance(0) != 0;

      //chi^2 = sum_i [ (E_i - sum_j c_j A_j)^2 ]
      // dchi2/dc_k = -2 sum_i [ E_i - sum_j c_j A_j ] A_k = 0
      // coeff multiplying c_j in dchi2/dc_k = sum_j A_j A_k
      // const term = E_i A_k
      
      if( passedcut ){
	//loop on ECAL and HCAL hits, increment matrices:
      
	double Eprime = T->eprime_eth;
	double Tp = T->Q2_eth/(2.0*0.938272);
	
	double xseedECAL, yseedECAL, eseedECAL, tseedECAL; int idseedECAL;
	
	double Esum_good_old = 0.0;
	// double xclust_old = T->xECAL;
	// double yclust_old = T->yECAL;
	double Esum_good_new = 0.0;
	
	bool goodcluster = true;

	int nbadcalib_thiscluster = 0;
	
	if( iter > 0 ){ //on calibration iterations after the first one (if applicable), calculate new EoverP value using the result from previous iteration and check it 
	  for( int i=0; i<T->nblkECAL; i++ ){
	    double A_i = T->ECAL_clusblk_e[i]/T->ECAL_clusblk_again[i];
	    double t_i = T->ECAL_clusblk_atime[i];
	    double x_i = T->ECAL_clusblk_x[i];
	    double y_i = T->ECAL_clusblk_y[i];
	    double e_i = A_i * PreviousSolutionECAL[T->ECAL_clusblk_id[i]];
	    if( i == 0 ){
	      tseedECAL = t_i;
	      xseedECAL = x_i;
	      yseedECAL = y_i;
	      eseedECAL = e_i;
	    }

	    if( fabs( t_i - tseedECAL ) <= dtcut_ECAL &&
		sqrt(pow(x_i - xseedECAL,2) + pow( y_i-yseedECAL,2)) <= Rmax_ECAL &&
		e_i > hit_thresh_ECAL ){
	      Esum_good_new += e_i;
	    }

	    if( !goodcalibECAL[T->ECAL_clusblk_id[i]] ){
	      nbadcalib_thiscluster++;
	    }
	  }

	  if( eseedECAL < seed_thresh_ECAL ) goodcluster = false;
	  
	  if ( fabs( Esum_good_new/Eprime - meanEoverP ) > 3.5*sigmaEoverP ){
	    goodcluster = false;
	  }
	}
	
	if( goodcluster ){
	  for( int i=0; i<T->nblkECAL; i++ ){
	    int id_i = T->ECAL_clusblk_id[i];
	    
	    double x_i = T->ECAL_clusblk_x[i];
	    double y_i = T->ECAL_clusblk_y[i];
	    double e_i = T->ECAL_clusblk_e[i];
	    double t_i = T->ECAL_clusblk_atime[i];
	    double oldgain_i = T->ECAL_clusblk_again[i];
	    double A_i = e_i / oldgain_i;
	    
	    if( iter > 0 ){
	      e_i = A_i * PreviousSolutionECAL[id_i];
	    }
	    
	    if( i == 0 ){
	      xseedECAL = T->ECAL_clusblk_x[i];
	      yseedECAL = T->ECAL_clusblk_y[i];
	      eseedECAL = e_i;
	      tseedECAL = T->ECAL_clusblk_atime[i];
	      idseedECAL = T->ECAL_clusblk_id[i];
	    }	
	    
	    if( eseedECAL < seed_thresh_ECAL ) continue;
	    
	    bool goodblock_i = true;
	    
	    //	if( i > 0 ){
	    double distfromseed = sqrt(pow(x_i-xseedECAL,2)+pow(y_i-yseedECAL,2));
	    if( fabs( t_i - tseedECAL ) > dtcut_ECAL ||
		distfromseed > Rmax_ECAL || e_i < hit_thresh_ECAL ){
	      goodblock_i = false;
	    }
	    
	    if( !goodblock_i ) continue;

	    // If this block DOESN'T have a good calibration and it is the only one
	    // in an otherwise good cluster that doesn't, calculate the average of the
	    // ratio of the calibrated energy of this channel to the "missing" energy from the cluster:
	    
	    if( iter > 0 && !goodcalibECAL[id_i] && nbadcalib_thiscluster == 1 ){
	      // Esum_good_new - e_i is the sum of all the other energies in the cluster
	      // e_i / (Eprime - Eother) is the ratio of the calibrated energy of this block
	      // to the "missing energy" from the cluster
	      sumEmiss_ratio[id_i] += e_i/(Eprime - (Esum_good_new-e_i) );
	      neventsEmiss_ratio[id_i] += 1.0;
	    }
	      
	    if( iter == 0 ) Esum_good_old += e_i;
	    
	    for( int j=0; j<T->nblkECAL; j++ ){
	      int id_j = T->ECAL_clusblk_id[j];
	      
	      double x_j = T->ECAL_clusblk_x[j];
	      double y_j = T->ECAL_clusblk_y[j];
	      double e_j = T->ECAL_clusblk_e[j];
	      double t_j = T->ECAL_clusblk_atime[j];
	      double oldgain_j = T->ECAL_clusblk_again[j];
	      double A_j = e_j / oldgain_j;

	      if( iter > 0 ) e_j = A_j * PreviousSolutionECAL[id_j];
	      
	      bool goodblock_j = true;
	      
	      distfromseed = sqrt(pow(x_j-xseedECAL,2)+pow(y_j-yseedECAL,2));
	      
	      if( fabs( t_j - tseedECAL ) > dtcut_ECAL ||
		  distfromseed > Rmax_ECAL || e_j < hit_thresh_ECAL ){
		goodblock_j = false;
	      }
	      
	      if( !goodblock_j ) continue;
	      
	      MECAL(id_i, id_j) += A_i * A_j;
	      
	    }
	  
	    bECAL(id_i) += A_i * Eprime;
	  
	    Nhit_chan_ECAL[id_i] += 1;
	  
	  }
	}

	if( iter == 0 ){
	  hEoverP_old_ECAL->Fill( Esum_good_old / Eprime );
	  hEoverP_x_old_ECAL->Fill( T->xECAL, Esum_good_old / Eprime );
	  hEoverP_y_old_ECAL->Fill( T->yECAL, Esum_good_old / Eprime );
	  Esum_good_new = Esum_good_old;
	}

	hEoverP_temp_ECAL->Fill( Esum_good_new/Eprime );

	///////////////////// Start HCAL stuff HERE: //////////////////////////
	// EVERYTHING WITH HCAL HAS TO ACCOUNT FOR THE FACT THAT THE "ID" OF //
	// HCAL BLOCKS STARTS AT 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! //
	///////////////////////////////////////////////////////////////////////
	
	// In addition to the global cut, we also require small-angle scattering
	// to use this event for HCAL calibration:
	if( T->thetaFPP <= thetaFPPmaxHCAL ){
	  bool goodhcalclust = true;
	  int nbadcalib_thishcalclust = 0;

	  double xseedHCAL, yseedHCAL, eseedHCAL, tseedHCAL; int idseedHCAL;

	  double Esum_good_hcal_old = 0.0;
	  double Esum_good_hcal_new = 0.0;
	  if( iter > 0 ){ //on calibration iterations after the first (if applicable), calculate new
	    //HCAL energy sum using the result of the previous iteration, check it:
	    for( int i=0; i<T->nblkHCAL; i++ ){
	      double A_i = T->HCAL_clusblk_e[i]/T->HCAL_clusblk_again[i];
	      double t_i = T->HCAL_clusblk_atime[i];
	      double x_i = T->HCAL_clusblk_x[i];
	      double y_i = T->HCAL_clusblk_y[i];
	      int id_i = T->HCAL_clusblk_id[i]-1;
	      double e_i = A_i * PreviousSolutionHCAL[id_i];
	      if( i == 0 ){
		tseedHCAL = t_i;
		xseedHCAL = x_i;
		yseedHCAL = y_i;
		eseedHCAL = e_i;
	      }
	      if( fabs( t_i - tseedHCAL ) <= dtcut_HCAL &&
		  sqrt(pow(x_i-xseedHCAL,2)+pow(y_i-yseedHCAL,2))<=Rmax_HCAL &&
		  e_i > hit_thresh_HCAL ){
		Esum_good_hcal_new += e_i;
	      }
	      if( !goodcalibHCAL[id_i] ){
		nbadcalib_thishcalclust++;
	      }
	    }
	    if( eseedHCAL  < seed_thresh_HCAL ) goodhcalclust = false;
	    
	    //	    if( fabs( Esum_good_hcal_new/Tp - SampFracMeanHCAL ) > 3.5*SampFracSigmaHCAL ){
	    //  goodhcalclust = false;
	    //}
	    
	    //if( Esum_good_hcal_new/Tp < SampFracMeanHCAL - 2.0*SampFracSigmaHCAL || Esum_good_hcal_new/Tp > SampFracMeanHCAL + 2.0*SampFracSigmaHCAL ) goodhcalclust = false;

	    // if( Esum_good_hcal_new/Tp < SampFracMinHCAL || Esum_good_hcal_new/Tp > SampFracMaxHCAL ){
	    //   goodhcalclust = false;
	    // }
	  }

	  if( goodhcalclust ){

	    bool useevent_in_calib = false;
	    if( iter == 0 || (Esum_good_hcal_new/Tp >= SampFracMinHCAL && Esum_good_hcal_new/Tp <= SampFracMaxHCAL ) ){
	      useevent_in_calib = true;
	    }

	    if( useevent_in_calib ){
	      for( int i=0; i<T->nblkHCAL; i++ ){
		//	      double A_i = T->HCAL_clusblk_e[i]/T->HCAL_clusblk_again[i];
		double t_i = T->HCAL_clusblk_atime[i];
		double x_i = T->HCAL_clusblk_x[i];
		double y_i = T->HCAL_clusblk_y[i];
		int id_i = T->HCAL_clusblk_id[i]-1;
		double oldgain_i = T->HCAL_clusblk_again[i];
		double e_i = T->HCAL_clusblk_e[i];
		double A_i = e_i/oldgain_i;
		
		if( iter > 0 ){
		  e_i = A_i * PreviousSolutionHCAL[id_i];
		}
		
		if( i == 0 ){
		  xseedHCAL = x_i;
		  yseedHCAL = y_i;
		  eseedHCAL = e_i;
		  tseedHCAL = t_i;
		  idseedHCAL = id_i;
		}

		if( eseedHCAL < seed_thresh_HCAL ) continue;

		bool goodblock_i = true;
		
		double distfromseed = sqrt(pow(x_i-xseedHCAL,2)+pow(y_i-yseedHCAL,2));
		if( fabs( t_i - tseedHCAL ) > dtcut_HCAL ||
		    distfromseed > Rmax_HCAL || e_i < hit_thresh_HCAL ){
		  goodblock_i = false;
		}
		if( !goodblock_i ) continue;
		
		if( iter > 0 && !goodcalibHCAL[id_i] && nbadcalib_thishcalclust == 1 ){
		  sumEmiss_ratio[id_i] += e_i/(Tp - (Esum_good_hcal_new-e_i) );
		  neventsEmiss_ratio[id_i] += 1.0;
		}
		
		if( iter == 0 ) Esum_good_hcal_old += e_i;
		
		
		
		for( int j=0; j<T->nblkHCAL; j++ ){
		  int id_j = T->HCAL_clusblk_id[j]-1;
		  double x_j = T->HCAL_clusblk_x[j];
		  double y_j = T->HCAL_clusblk_y[j];
		  double e_j = T->HCAL_clusblk_e[j];
		  double t_j = T->HCAL_clusblk_atime[j];
		  double oldgain_j = T->HCAL_clusblk_again[j];
		  double A_j = e_j/oldgain_j;
		  
		  if( iter > 0 ) e_j = A_j * PreviousSolutionHCAL[id_j];
		  
		  bool goodblock_j = true;
		  
		  distfromseed = sqrt(pow(x_j-xseedHCAL,2)+pow(y_j-yseedHCAL,2));
		  if( fabs( t_j - tseedHCAL ) > dtcut_HCAL ||
		      distfromseed > Rmax_HCAL || e_j < hit_thresh_HCAL ){
		    goodblock_j = false;
		  }
		  if( !goodblock_j ) continue;
		  
		  MHCAL(id_i, id_j) += A_i * A_j;
		}
		
		bHCAL(id_i) += A_i * Tp;
		
		Nhit_chan_HCAL[id_i] += 1;
	      } // loop over HCAL blocks
	    } //check if event is to be used in calibration (first iteration or sampling fraction within limits
	    
	    if( iter == 0 ){
	      hsampfrac_old_HCAL->Fill( Esum_good_hcal_old / Tp );
	      hsampfrac_x_old_HCAL->Fill( T->xHCAL, Esum_good_hcal_old / Tp );
	      hsampfrac_y_old_HCAL->Fill( T->yHCAL, Esum_good_hcal_old / Tp );
	      Esum_good_hcal_new = Esum_good_hcal_old;
	    }
	    
	    hSampFrac_temp_HCAL->Fill( Esum_good_hcal_new / Tp );
	  } //check if cluster is "good"
	
	} // small-angle cut on FPP track;
      }
    
      if( nevent % 1000 == 0 ) cout << "Event " << nevent << endl;
    }

    c1->cd(1);
    hEoverP_temp_ECAL->Draw();
    c1->cd(2);
    hSampFrac_temp_HCAL->Draw();
    
    list_uncalibrated.clear();
    list_uncalibrated_HCAL.clear();
    //list of uncalibrated channels:

    for( int i=0; i<nchan_ECAL; i++ ){
      if( Nhit_chan_ECAL[i] < minhits_ECAL ){
	bECAL(i) = 1.0;
	MECAL(i,i) = 1.0;
	
	list_uncalibrated.insert( i );
	for( int j=0; j<nchan_ECAL; j++ ){
	  if( j != i ) {
	    MECAL(i,j) = 0.0;
	    MECAL(j,i) = 0.0; //Not necessary, but also not harmful;;
	  }
	}
      }
    }


    for( int i=0; i<nchan_HCAL; i++ ){
      if( Nhit_chan_HCAL[i] < minhits_HCAL ){
	bHCAL(i) = 1.0;
	MHCAL(i,i) = 1.0;
	list_uncalibrated_HCAL.insert(i);

	for( int j=0; j<nchan_HCAL; j++ ){
	  if( j != i ){
	    MHCAL(i,j) = 0.0;
	    MHCAL(j,i) = 0.0;
	  }
	}
      }
    }
    
    cout << "constructing SVD matrix... ";
    
    TDecompSVD AECAL(MECAL);
    
    cout << "done. Solving for coefficients... ";
    AECAL.Solve( bECAL );
    cout << "done" << endl;

    cout << "constructing HCAL SVD matrix... ";
    TDecompSVD AHCAL(MHCAL);
    cout << "done. Solving for HCAL coefficients... ";
    AHCAL.Solve( bHCAL );
    cout << "done" << endl;
    
    if( iter > 0 ){

      c1->cd(1);
      hEoverP_temp_ECAL->Fit("gaus","SQ","",meanEoverP-1.5*sigmaEoverP,meanEoverP+1.5*sigmaEoverP);

      meanEoverP = ( (TF1*) (hEoverP_temp_ECAL->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
      sigmaEoverP = ( (TF1*) (hEoverP_temp_ECAL->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");


      // hSampFrac_temp_HCAL->Fit("gaus", "SQ", "", SampFracMeanHCAL-0.7*SampFracSigmaHCAL,
      // 			       SampFracMeanHCAL + 2.0*SampFracSigmaHCAL );

      //     SampFracMeanHCAL = ( (TF1*) (hSampFrac_temp_HCAL->GetListOfFunctions()->FindObject("gaus")))->GetParameter("Mean");
      //SampFracSigmaHCAL = ( (TF1*) (hSampFrac_temp_HCAL->GetListOfFunctions()->FindObject("gaus")))->GetParameter("Sigma");

      c1->cd(2);
      int maxbinhcal = hSampFrac_temp_HCAL->GetMaximumBin();

      double HCALmax = hSampFrac_temp_HCAL->GetBinContent(maxbinhcal);
      
      int binlow = maxbinhcal;
      int binhigh = maxbinhcal;

      int binlow7 = maxbinhcal;
      int binhigh7 = maxbinhcal;
      
      while( binlow >= 1 && hSampFrac_temp_HCAL->Integral(1,binlow) > 0.001*hSampFrac_temp_HCAL->Integral() ){
	binlow--;
      }

      int nbinsxhcal = hSampFrac_temp_HCAL->GetNbinsX();
      while( binhigh <= nbinsxhcal && hSampFrac_temp_HCAL->Integral(binhigh,nbinsxhcal) > 0.02*hSampFrac_temp_HCAL->Integral() ){
	binhigh++;
      }

      while( binlow7 >= 1 && hSampFrac_temp_HCAL->GetBinContent(binlow7) >= 0.7*HCALmax ){binlow7--;}
      while( binhigh7 <= nbinsxhcal && hSampFrac_temp_HCAL->GetBinContent(binhigh7) >= 0.7*HCALmax){binhigh7++;}
      
      SampFracMeanHCAL = hSampFrac_temp_HCAL->GetMean();
      SampFracSigmaHCAL = hSampFrac_temp_HCAL->GetRMS();

      SampFracMinHCAL = hSampFrac_temp_HCAL->GetBinLowEdge(binlow+1);
      SampFracMaxHCAL = hSampFrac_temp_HCAL->GetBinLowEdge(binhigh);

      hSampFrac_temp_HCAL->Fit("gaus","SQ","",hSampFrac_temp_HCAL->GetBinLowEdge(binlow7),hSampFrac_temp_HCAL->GetBinLowEdge(binhigh7+1));

      SampFracMeanHCAL = ( (TF1*) (hSampFrac_temp_HCAL->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Mean");
      SampFracSigmaHCAL = ( (TF1*) (hSampFrac_temp_HCAL->GetListOfFunctions()->FindObject("gaus") ) )->GetParameter("Sigma");
      c1->Update();
      gPad->Modified();
      gSystem->ProcessEvents();
    }

    double ngoodcalib = 0;
    double sum_calib = 0.0;
    for( int i=0; i<nchan_ECAL; i++ ){
      if( list_uncalibrated.find( i ) == list_uncalibrated.end() && bECAL(i) >= mincoeff_ECAL && bECAL(i) <= maxcoeff_ECAL ){
	ngoodcalib += 1.0;
	sum_calib += bECAL(i);
	PreviousSolutionECAL[i] = bECAL(i);
	goodcalibECAL[i] = true;
      } else {
	goodcalibECAL[i] = false;
	if( iter > 0 ){
	  if( neventsEmiss_ratio[i] >= 10. ){
	    double avgEmiss_ratio = sumEmiss_ratio[i]/neventsEmiss_ratio[i]; 
	    // ratio of e_i / (Etrue - Eother) ~= e_i / e_i(true)
	    // NOTE: we will have used "PreviousSolution" to do the energy reconstruction in the preceding loop.
	    double newcoeff = PreviousSolutionECAL[i] / avgEmiss_ratio; // old coeff / (e_i/e_i(true)) = old coeff * e_i("true")/e_i (calib)
	    if( newcoeff >= mincoeff_ECAL && newcoeff <= maxcoeff_ECAL ){
	      cout << "NOTE: ECAL channel " << i << " gain being calibrated based on missing energy ratio, (ratio, newcoeff)=(" << avgEmiss_ratio << ", " << newcoeff << ")" << endl;
	      goodcalibECAL[i] = true;
	      PreviousSolutionECAL[i] = newcoeff; //unclear how many additional channels we might be able to calibrate this way.
	      bECAL(i) = newcoeff; //Just replace bECAL with new coeff here. This should fix the possibility to have factors outside the user-defined range
	    }
	  }
	}
      }
    }

    double avgcoeff = sum_calib/ngoodcalib;

    for( int i=0; i<nchan_ECAL; i++ ){
      if( !goodcalibECAL[i] ){ //replace this channel's gain with the average gain of all "good" channels or with the default gain from the file if provided by the user
	if( gotoldgain_ECAL ){
	  bECAL(i) = defaultGainECAL[i];
	  PreviousSolutionECAL[i] = bECAL(i);
	} else { // if the user didn't provide anything, use the average "good" calibration result:
	  bECAL(i) = avgcoeff;
	  PreviousSolutionECAL[i] = avgcoeff;
	}
      }
    }

    cout << "After iteration " << iter << " of ECAL calib: (ngood, meanEoverP, sigmaEoverP)=("
	 << ngoodcalib << ", " << meanEoverP << ", " << sigmaEoverP << ")" << endl;
    
    //Now do HCAL stuff:

    //Can't hurt to re-use these local variables:
    ngoodcalib = 0.0;
    sum_calib = 0.0;

    for( int i=0; i<nchan_HCAL; i++ ){
      if( list_uncalibrated_HCAL.find( i ) == list_uncalibrated_HCAL.end() && bHCAL(i) >= mincoeff_HCAL && bHCAL(i) <= maxcoeff_HCAL ){
	ngoodcalib += 1.0;
	sum_calib += bHCAL(i);
	PreviousSolutionHCAL[i] = bHCAL(i);
	goodcalibHCAL[i] = true;
      } else {
	goodcalibHCAL[i] = false;
	if( iter > 0 ){
	  if( neventsEmiss_ratio_HCAL[i] >= 10.0 ){
	    double avgEmiss_ratio = sumEmiss_ratio_HCAL[i]/neventsEmiss_ratio_HCAL[i];
	    double newcoeff = PreviousSolutionHCAL[i]/avgEmiss_ratio;
	    if( newcoeff >= mincoeff_HCAL && newcoeff <= maxcoeff_HCAL ){
	      cout << "NOTE: HCAL channel " << i << " gain being calibrated based on missing energy ratio, (ratio, newcoeff)=(" << avgEmiss_ratio << ", " << newcoeff << ")" << endl;
	      goodcalibHCAL[i] = true;
	      PreviousSolutionHCAL[i] = newcoeff;
	      bHCAL(i) = newcoeff;
	    }
	  }
	}
      }
    }

    avgcoeff = sum_calib/ngoodcalib;

    for( int i=0; i<nchan_HCAL; i++ ){
      if( !goodcalibHCAL[i] ){
	if( gotoldgain_HCAL ){
	  bHCAL(i) = defaultGainHCAL[i];
	  PreviousSolutionHCAL[i] = bHCAL(i);
	} else {
	  bHCAL(i) = avgcoeff;
	  PreviousSolutionHCAL[i] = avgcoeff;
	}
      }
    }

    cout << "after iteration " << iter << " of HCAL calib: (ngood, meanSampFrac, sigmaSampFrac)=("
	 << ngoodcalib << ", " << SampFracMeanHCAL << ", " << SampFracSigmaHCAL << ")" << endl;
    
  } //end loop over iterations:


  //Before we write the DB file, calculate "fudge factors" to apply to both ECAL and HCAL gains:
  //int binmax_ECAL

  double fudgefactor_ECAL = 1.0/meanEoverP;
  double fudgefactor_HCAL = 1.0/SampFracMeanHCAL;

  for( int i=0; i<nchan_ECAL; i++ ){
    bECAL(i) *= fudgefactor_ECAL;
  }
  for( int i=0; i<nchan_HCAL; i++ ){
    bHCAL(i) *= fudgefactor_HCAL;
  }
  
  //Now loop over the data again and see how we did:

  nevent = 0;

  oldtreenum = -1;

  TH1D *hEoverP_new_ECAL = new TH1D("hEoverP_new_ECAL","NEW;ECAL E/p; ", 250,0,5.0);
  TH2D *hEoverP_x_new_ECAL = new TH2D("hEoverP_x_new_ECAL", "NEW;xECAL (m);E/p",160,-1.6,1.6,250,0,5);
  TH2D *hEoverP_y_new_ECAL = new TH2D("hEoverP_y_new_ECAL", "NEW;yECAL (m);E/p",160,-0.8,0.8,250,0,5);

  TH1D *hGainCoeffNew = new TH1D("hGainCoeffNew", "NEW ECAL Gain Coeffs; Gain (GeV/pC);", 200, 0.0, 0.2 ); 

  double chan[nchan_ECAL], newcoeff[nchan_ECAL], rationewold[nchan_ECAL];

  double sumcoeff = 0.0, ngoodcalib = 0.0;

  //  cout << "(mincoeff_ECAL,maxcoeff_ECAL)=(" << mincoeff_ECAL << ", " << maxcoeff_ECAL << ")" << endl;
  //New HCAL diagnostics:

  TH1D *hSampFracNew_HCAL = new TH1D("hSampFracNew_HCAL", "NEW; HCAL E/T_{p};", 250,0,5);
  TH2D *hSampFracNew_x_HCAL = new TH2D("hSampFracNew_x_HCAL", "NEW; x_{HCAL} (m); E/T_{p}", 200, -2.75, 1.25, 250, 0, 5 );
  TH2D *hSampFracNew_y_HCAL = new TH2D("hSampFracNew_y_HCAL", "NEW; y_{HCAL} (m); E/T_{p}", 180, -0.9, 0.9, 250, 0, 5 );

  TH1D *hGainCoeffNew_HCAL = new TH1D("hGainCoeffNew_HCAL", "NEW HCAL Gain Coeffs; GeV/pC;", 200, 0.0, 0.5 );

  double hcalchan[nchan_HCAL], newhcalcoeff[nchan_HCAL], rationewoldhcal[nchan_HCAL];
  double sumcoeff_hcal = 0.0, ngoodcalibhcal = 0.0;
  
  for( int i=0; i<nchan_ECAL; i++ ){
    //if( list_uncalibrated.find( i ) == list_uncalibrated.end() ){
    if( goodcalibECAL[i] && list_uncalibrated.find(i) == list_uncalibrated.end() ){
      hGainCoeffNew->Fill( bECAL(i) );
      sumcoeff += bECAL(i);
      ngoodcalib += 1.0;
      newcoeff[i] = bECAL(i);
    
      
      if( gotoldgain_ECAL ){
	rationewold[i] = newcoeff[i]/defaultGainECAL[i]; 
      }
    }
    
    chan[i] = i;
  }
  
  double avgcoeff = sumcoeff/ngoodcalib;

  
  for( int i=0; i<nchan_ECAL; i++ ){
    if( list_uncalibrated.find(i) != list_uncalibrated.end() || !goodcalibECAL[i] ){
      if( gotoldgain_ECAL ){
	newcoeff[i] = defaultGainECAL[i]*fudgefactor_ECAL;
	rationewold[i] = newcoeff[i] / defaultGainECAL[i]; //1
	//newcoeff[i] *= fudgefactor_ECAL; // for database:
      } else {
	newcoeff[i] = avgcoeff;
      }
      
    }
  }

  for( int i=0; i<nchan_HCAL; i++ ){
    if( goodcalibHCAL[i] && list_uncalibrated_HCAL.find(i) == list_uncalibrated_HCAL.end() ){
      hGainCoeffNew_HCAL->Fill( bHCAL(i) );
      sumcoeff_hcal += bHCAL(i);
      ngoodcalibhcal += 1.0;
      newhcalcoeff[i] = bHCAL(i);

      if( gotoldgain_HCAL ){
	rationewoldhcal[i] = newhcalcoeff[i]/defaultGainHCAL[i];
      }
    }
    hcalchan[i] = i+1;
  }

  double avgcoeffhcal = sumcoeff_hcal / ngoodcalibhcal;

  for( int i=0; i<nchan_HCAL; i++ ){
    if( list_uncalibrated_HCAL.find(i) != list_uncalibrated_HCAL.end() || !goodcalibHCAL[i] ){
      if( gotoldgain_HCAL ){
	newhcalcoeff[i] = defaultGainHCAL[i]*fudgefactor_HCAL;
	rationewoldhcal[i] = newhcalcoeff[i]/defaultGainHCAL[i];
      } else {
	newcoeff[i] = avgcoeffhcal;
      }
    }
  }
  
  TString dbfilename = outfilename;

  dbfilename.ReplaceAll(".root",".dat");
  dbfilename.Prepend("db_");
  
  ofstream dbfile( dbfilename.Data() );

  dbfile << "earm.ecal.adc.gain = " << endl;
  int rowtemp=0, coltemp=0;
  for( int i=0; i<nchan_ECAL; i++ ){
    if( coltemp == ncols[rowtemp] ){
      rowtemp++;
      coltemp = 0;
      dbfile << endl;
    }
    dbfile << Form( "%10.6g ", newcoeff[i] );
    coltemp++;
  }

  TGraph *gNewCoeff = new TGraph( nchan_ECAL, chan, newcoeff );
  gNewCoeff->SetMarkerStyle(20);
  fout->cd();
  gNewCoeff->Write("gNewCoeff");

  if( gotoldgain_ECAL ){
    TGraph *gNewOldRatio = new TGraph( nchan_ECAL, chan, rationewold );
    gNewOldRatio->SetMarkerStyle(20);

    TH1D *hNewOldRatio = new TH1D("hNewOldRatio", "ECAL NEW/OLD; Ratio;", 500,0,5);

    for( int i=0; i<nchan_ECAL; i++ ){
      hNewOldRatio->Fill( rationewold[i] );
    }
    
    fout->cd();
    gNewOldRatio->Write("gNewOldRatioECAL");
  }

  /// HCAL stuff:
  dbfile << endl << endl << "########################## HCAL ##############################" << endl
	 << endl;
  dbfile << "sbs.hcal.adc.gain = " << endl;
  for( int i=0; i<nchan_HCAL; i++ ){
    if( i%12 == 0 ) dbfile << endl;
    dbfile << Form( "%10.6g ", newhcalcoeff[i] );
  }

  TGraph *gNewHCALcoeff = new TGraph(nchan_HCAL, hcalchan, newhcalcoeff );
  gNewHCALcoeff->SetMarkerStyle(21);
  gNewHCALcoeff->SetMarkerColor(2);
  fout->cd();
  gNewHCALcoeff->Write("gNewHCALcoeff");

  if( gotoldgain_HCAL ){
    TGraph *gNewOldRatioHCAL = new TGraph(nchan_HCAL, hcalchan, rationewoldhcal );

    gNewOldRatioHCAL->SetMarkerStyle(21);
    gNewOldRatioHCAL->SetMarkerColor(2);
    
    TH1D *hNewOldRatioHCAL = new TH1D( "hNewOldRatioHCAL", "HCAL NEW/OLD; Ratio;", 125, 0, 5 );
    for( int i=0; i<nchan_HCAL; i++){
      hNewOldRatioHCAL->Fill( rationewoldhcal[i] );
    }
    fout->cd();
    gNewOldRatioHCAL->Write( "gNewOldRatioHCAL" );
  }
  
  while( C->GetEntry(nevent++) ){
    treenum = C->GetTreeNumber();
    
    if( treenum != oldtreenum ){
      oldtreenum = treenum;
      GlobalCut->UpdateFormulaLeaves();
    }
    
    bool passedcut = GlobalCut->EvalInstance(0) != 0;
    
    //chi^2 = sum_i [ (E_i - sum_j c_j A_j)^2 ]
    // dchi2/dc_k = -2 sum_i [ E_i - sum_j c_j A_j ] A_k = 0
    // coeff multiplying c_j in dchi2/dc_k = sum_j A_j A_k
    // const term = E_i A_k

    if( passedcut ){
      //loop on ECAL and HCAL hits, increment matrices:
      

      double Eprime = T->eprime_eth;
      double Tp = T->Q2_eth/(2.0*0.938272);
      
      double xseedECAL, yseedECAL, eseedECAL, tseedECAL; int idseedECAL;
      double xseedHCAL, yseedHCAL, eseedHCAL, tseedHCAL; int idseedHCAL;
      
      double Esum_good_new = 0.0;
      double xclust_old = T->xECAL;
      double yclust_old = T->yECAL;
      
      for( int i=0; i<T->nblkECAL; i++ ){
	int id_i = T->ECAL_clusblk_id[i];
	double x_i = T->ECAL_clusblk_x[i];
	double y_i = T->ECAL_clusblk_y[i];
	double e_i = T->ECAL_clusblk_e[i];
	double t_i = T->ECAL_clusblk_atime[i];
	double oldgain_i = T->ECAL_clusblk_again[i];
	double A_i = e_i / oldgain_i;

	e_i = A_i * bECAL(id_i);
	
	if( i == 0 ){
	  // xseedECAL = T->ECAL_clusblk_x[i];
	  // yseedECAL = T->ECAL_clusblk_y[i];
	  // eseedECAL = T->ECAL_clusblk_e[i];
	  // tseedECAL = T->ECAL_clusblk_atime[i];
	  // idseedECAL = T->ECAL_clusblk_id[i];
	  xseedECAL = x_i;
	  yseedECAL = y_i;
	  eseedECAL = e_i;
	  tseedECAL = t_i;
	  idseedECAL = id_i;
	}

	if( eseedECAL < seed_thresh_ECAL ) continue;

	bool goodblock_i = true;

	//	if( i > 0 ){
	double distfromseed = sqrt(pow(x_i-xseedECAL,2)+pow(y_i-yseedECAL,2));
	if( fabs( t_i - tseedECAL ) > dtcut_ECAL ||
	    distfromseed > Rmax_ECAL || e_i < hit_thresh_ECAL ){
	  goodblock_i = false;
	}
	
	if( !goodblock_i ) continue;
	
	Esum_good_new += e_i;
      }

      hEoverP_new_ECAL->Fill( Esum_good_new/ Eprime );
      hEoverP_x_new_ECAL->Fill( T->xECAL, Esum_good_new / Eprime );
      hEoverP_y_new_ECAL->Fill( T->yECAL, Esum_good_new / Eprime );


      if( T->thetaFPP < thetaFPPmaxHCAL ){

	double Esum_good_new_hcal = 0.0;
	
	for( int i=0; i<T->nblkHCAL; i++ ){
	  int id_i = T->HCAL_clusblk_id[i]-1;
	  double x_i = T->HCAL_clusblk_x[i];
	  double y_i = T->HCAL_clusblk_y[i];
	  double e_i = T->HCAL_clusblk_e[i];
	  double t_i = T->HCAL_clusblk_atime[i];
	  double oldgain_i = T->HCAL_clusblk_again[i];
	  double A_i = e_i / oldgain_i;

	  e_i = A_i * bHCAL(id_i);

	  if( i == 0 ){
	    xseedHCAL = x_i;
	    yseedHCAL = y_i;
	    eseedHCAL = e_i;
	    tseedHCAL = t_i;
	    idseedHCAL = id_i;
	  }

	  if( eseedHCAL < seed_thresh_HCAL ) continue;
	  
	  bool goodblock_i = true;
	  double distfromseed = sqrt(pow(x_i-xseedHCAL,2)+pow(y_i-yseedHCAL,2));

	  if( fabs( t_i - tseedHCAL ) > dtcut_HCAL ||
	      distfromseed > Rmax_HCAL || e_i < hit_thresh_HCAL ){
	    goodblock_i = false;
	  }
	  if( !goodblock_i ) continue;

	  Esum_good_new_hcal += e_i;
	  
	}

	hSampFracNew_HCAL->Fill( Esum_good_new_hcal/Tp );
	hSampFracNew_x_HCAL->Fill( T->xHCAL, Esum_good_new_hcal/Tp );
	hSampFracNew_y_HCAL->Fill( T->yHCAL, Esum_good_new_hcal/Tp );

	
							
					    

      } // check small-angle cut for FPP tracks:
    }
    if( nevent % 1000 == 0 ) cout << "Final loop, checking results, event " << nevent << endl;
  }
	
  
  fout->Write();
  

}
