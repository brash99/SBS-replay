//#include "GEM_cosmic_tracks.C"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TEventList.h"
#include "TCut.h"
#include <iostream>
#include <fstream>
#include "TMinuit.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TEventList.h"

double PI = TMath::Pi();

double Mp = 0.938272;
double Mn = 0.939565;

void MomentumCalibration( const char *configfilename, const char *outputfilename="NewMomentumFit.txt" ){
  
  ifstream configfile(configfilename);

  //so as usual the main things to define are:
  // 1. List of files
  // 2. Global cuts
  // 3. Nominal beam energy (perhaps corrected for average energy loss along half target thickness
  // 4. File name for old momentum coefficients (optional?)
  // 5. BigBite central angle
  // 6. SBS central angle
  // 7. HCAL distance

  TChain *C = new TChain("T");
  
  TString currentline;

  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline);
    }
  }

  TCut globalcut="";
  
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
    }
  }

  //int fitorder = 2;
  
  double ebeam=1.916;
  double bbtheta = 51.0*TMath::DegToRad();
  double sbstheta = 34.5*TMath::DegToRad();
  double hcaldist = 13.5; //meters
  double Ltgt = 15.0; //cm
  double rho_tgt = 0.0723; //g/cc

  double celldiameter = 5.0; //cm, right now this is a guess
  
  double Ztgt = 1.0;
  double Atgt = 1.0;
  double Mmol_tgt = 1.008; //g/mol
  
  //For energy-loss correction to beam energy:
  double dEdx_tgt=0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy:

  //Mean energy loss of the beam prior to the scattering:
  double MeanEloss = Ltgt/2.0 * rho_tgt * dEdx_tgt; //approximately 3 MeV:

  double MeanEloss_outgoing = celldiameter/2.0/sin(bbtheta) * rho_tgt * dEdx_tgt; //Approximately 1 MeV
  
  int usehcalcut = 0; //default to false:

  double dx0=0.0, dy0=0.0; //zero offsets for deltax and deltay cuts:

  double dxsigma = 0.06, dysigma = 0.1; //sigmas for deltax and deltay cuts:

  double GEMpitch = 10.0*TMath::DegToRad();

  double dpelmin_fit = -0.03;
  double dpelmax_fit = 0.03;
  double Wmin_fit = 0.86;
  double Wmax_fit = 1.02; 
  
  int order = 2; 
  
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endconfig") ){
    if( !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" ");

      int ntokens = tokens->GetEntries();

      if( ntokens >= 2 ){
	TString skey = ( (TObjString*) (*tokens)[0] )->GetString();

	if ( skey == "order" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  order = stemp.Atoi();
	}

	if( skey == "dEdx" ){ //assumed to be given in GeV/(g/cm^2)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dEdx_tgt = stemp.Atof();
	}

	if( skey == "ebeam" ){ //primary beam energy (GeV):
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  ebeam = stemp.Atof();
	}

	if( skey == "bbtheta" ){ //BigBite central angle (deg) (assumed to be on beam left)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  bbtheta = stemp.Atof() * TMath::DegToRad();
	}

	if( skey == "sbstheta" ){ //SBS/HCAL central angle (deg) (assumed to be on beam right)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  sbstheta = stemp.Atof() * TMath::DegToRad();
	}

	if( skey == "hcaldist" ){ //HCAL distance from target (m)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  hcaldist = stemp.Atof();
	}

	if( skey == "Ltgt" ){ //target length (cm):
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Ltgt = stemp.Atof();
	}

	if( skey == "rho_tgt" ){ //target density (g/cm^3)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  rho_tgt = stemp.Atof();
	}

	if( skey == "celldiameter" ){ //cell diameter (cm)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  celldiameter = stemp.Atof();
	}
	
	
	if( skey == "usehcalcut" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  usehcalcut = stemp.Atoi();
	}

	if( skey == "dx0" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dx0 = stemp.Atof();
	}

	if( skey == "dy0" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dy0 = stemp.Atof();
	}

	if( skey == "dxsigma" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dxsigma = stemp.Atof();
	}

	if( skey == "dysigma" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dysigma = stemp.Atof();
	}

	if( skey == "GEMpitch" ){ //calculate dxHCAL, dyHCAL and cut on them (2.5 sigma)
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  GEMpitch = stemp.Atof()*TMath::DegToRad();
	}

	if( skey == "dpel_min" ){ //p/pel - 1 minimum to include in fit
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dpelmin_fit = stemp.Atof();
	}

	if( skey == "dpel_max" ){ //p/pel - 1 maximum to include in fit
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  dpelmax_fit = stemp.Atof();
	}

	if( skey == "Wmin" ){ //Wmin to include in fit
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Wmin_fit = stemp.Atof();
	}

	if( skey == "Wmax" ){ //Wmax to include in fit
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  Wmax_fit = stemp.Atof();
	}
	
      }
      
      tokens->Delete();
    }
  }
  
  //Note that both of these calculations neglect the Aluminum end windows and cell walls:
  
  MeanEloss = Ltgt/2.0 * rho_tgt * dEdx_tgt; //hopefully the user has provided these quantities with correct units
  
  MeanEloss_outgoing = celldiameter/2.0/sin(bbtheta) * rho_tgt * dEdx_tgt; 

  
  cout << "Mean E loss = " << MeanEloss << endl;
  cout << "BigBite theta = " << bbtheta * TMath::RadToDeg() << endl;
  cout << "Mean E loss (outgoing electron) = " << MeanEloss_outgoing << endl;
  cout << "Ebeam (GeV) = " << ebeam - MeanEloss << endl;

						   
					   
					   
  
  // In order to get 

  TEventList *elist = new TEventList("elist","Event list for BigBite momentum calibration");

  C->Draw(">>elist",globalcut);

  int MAXNTRACKS=1000;

  //Declare variables to hold branch addresses:
  
  double ntrack;
  double px[MAXNTRACKS], py[MAXNTRACKS], pz[MAXNTRACKS], p[MAXNTRACKS];
  double vx[MAXNTRACKS], vy[MAXNTRACKS], vz[MAXNTRACKS];

  //Use "rotated" versions of the focal plane variables:
  double xfp[MAXNTRACKS], yfp[MAXNTRACKS], thfp[MAXNTRACKS], phfp[MAXNTRACKS];
  //Also need target variables:
  double xtgt[MAXNTRACKS], ytgt[MAXNTRACKS], thtgt[MAXNTRACKS], phtgt[MAXNTRACKS];
  
  double xHCAL, yHCAL, EHCAL;
  double EPS, ESH, xSH, ySH, xPS;

  //Ignore hodo variables for meow:
  // int maxHODOclusters=100;

  // double nHODOclusters;
  
  // double xHODO[maxHODOclusters],yHODO[maxHODOclusters];
  // double tmeanHODO[maxHODOclusters], tdiffHODO[maxHODOclusters];
  
  C->SetBranchStatus("*",0);
  //Minimal set of HCAL variables:
  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);

  //BigBite track variables:
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.p",1);
  C->SetBranchStatus("bb.tr.vx",1);
  C->SetBranchStatus("bb.tr.vy",1);
  C->SetBranchStatus("bb.tr.vz",1);
  
  C->SetBranchStatus("bb.tr.r_x",1);
  C->SetBranchStatus("bb.tr.r_y",1);
  C->SetBranchStatus("bb.tr.r_th",1);
  C->SetBranchStatus("bb.tr.r_ph",1);
  
  C->SetBranchStatus("bb.tr.tg_x",1);
  C->SetBranchStatus("bb.tr.tg_y",1);
  C->SetBranchStatus("bb.tr.tg_th",1);
  C->SetBranchStatus("bb.tr.tg_ph",1);
  
  //Shower and preshower variables:
  C->SetBranchStatus("bb.ps.e",1);
  C->SetBranchStatus("bb.ps.x",1);
  C->SetBranchStatus("bb.ps.y",1);
  C->SetBranchStatus("bb.sh.e",1);
  C->SetBranchStatus("bb.sh.x",1);
  C->SetBranchStatus("bb.sh.y",1);

  C->SetBranchAddress("bb.tr.n",&ntrack);
  C->SetBranchAddress("bb.tr.p",p);
  C->SetBranchAddress("bb.tr.px",px);
  C->SetBranchAddress("bb.tr.py",py);
  C->SetBranchAddress("bb.tr.pz",pz);
  C->SetBranchAddress("bb.tr.vx",vx);
  C->SetBranchAddress("bb.tr.vy",vy);
  C->SetBranchAddress("bb.tr.vz",vz);

  //Focal-plane track variables (use "rotated" versions):
  C->SetBranchAddress("bb.tr.r_x",xfp);
  C->SetBranchAddress("bb.tr.r_y",yfp);
  C->SetBranchAddress("bb.tr.r_th",thfp);
  C->SetBranchAddress("bb.tr.r_ph",phfp);

  //Target track variables (other than momentum):
  C->SetBranchAddress("bb.tr.tg_x",xtgt);
  C->SetBranchAddress("bb.tr.tg_y",ytgt);
  C->SetBranchAddress("bb.tr.tg_th",thtgt);
  C->SetBranchAddress("bb.tr.tg_ph",phtgt);

  C->SetBranchAddress("sbs.hcal.x",&xHCAL);
  C->SetBranchAddress("sbs.hcal.y",&yHCAL);
  C->SetBranchAddress("sbs.hcal.e",&EHCAL);

  C->SetBranchAddress("bb.ps.e",&EPS);
  C->SetBranchAddress("bb.sh.e",&ESH);
  C->SetBranchAddress("bb.ps.x",&xPS);
  C->SetBranchAddress("bb.sh.x",&xSH);
  C->SetBranchAddress("bb.sh.y",&ySH);
  
  int nparams = 0;

  vector<int> xtar_expon;
  vector<int> xfp_expon;
  vector<int> xpfp_expon;
  vector<int> yfp_expon;
  vector<int> ypfp_expon;
  
  
  for(int i=0; i<=order; i++){
    for(int j=0; j<=order-i; j++){
      for(int k=0; k<=order-i-j; k++){
	for(int l=0; l<=order-i-j-k; l++){
	  for(int m=0; m<=order-i-j-k-l; m++){
	    nparams++;
	    xtar_expon.push_back( i );
	    xfp_expon.push_back( m );
	    yfp_expon.push_back( l );
	    xpfp_expon.push_back( k );
	    ypfp_expon.push_back( j );
	  }
	}
      }
    }
  }

  //For our purposes here we only need to set up the augmented matrix for the p*thetabend matrix elements:
  TMatrixD M(nparams,nparams);
  TVectorD b_pth(nparams);

  for( int i=0; i<nparams; i++ ){
    for( int j=0; j<nparams; j++ ){
      M(i,j) = 0.0;
    }
    b_pth(i) = 0.0;
  }


  double pcentral = ebeam/(1.+ebeam/Mp*(1.-cos(bbtheta)));

  TFile *fout = new TFile("MomentumCalibration_temp.root","RECREATE");

  //After we get this working, we'll also want correlation plots of dpel and W vs. focal plane and target variables, new and old:
  TH1D *hdpel_old = new TH1D("hdpel_old","Old MEs;p/p_{elastic}(#theta)-1;",250,-0.125,0.125);
  TH1D *hW_old = new TH1D("hW_old","Old MEs;W (GeV);",250,Mp-0.2,Mp+0.2);
  TH1D *hp_old = new TH1D("hp_old","Old MEs;p (GeV);",250,0.5*pcentral,1.5*pcentral);
  TH1D *hpthetabend_old = new TH1D("hpthetabend_old","Old MEs; p#theta_{bend};",250,0.2,0.35);

  TH2D *hxpfp_vs_xfp_old = new TH2D("hxpfp_vs_xfp_old",";xfp (m); xpfp",250,-0.75,0.75,250,-0.4,0.4);
  
  TH2D *hdpel_xfp_old = new TH2D("hdpel_xfp_old", "Old MEs; xfp (m); p/p_{elastic}(#theta)-1",250,-0.75,0.75,250,-0.125,0.125);
  TH2D *hdpel_yfp_old = new TH2D("hdpel_yfp_old", "Old MEs; yfp (m); p/p_{elastic}(#theta)-1",250,-0.25,0.25,250,-0.125,0.125);
  TH2D *hdpel_xpfp_old = new TH2D("hdpel_xpfp_old", "Old MEs; xpfp; p/p_{elastic}(#theta)-1",250,-0.4,0.4,250,-0.125,0.125);
  TH2D *hdpel_ypfp_old = new TH2D("hdpel_ypfp_old", "Old MEs; ypfp; p/p_{elastic}(#theta)-1",250,-0.15,0.15,250,-0.125,0.125);

  TH2D *hpthetabend_xfp_old = new TH2D("hpthetabend_xfp_old", "Old MEs; xfp (m); p#theta_{bend}",250,-0.75,0.75,250,0.2,0.35);
  TH2D *hpthetabend_yfp_old = new TH2D("hpthetabend_yfp_old", "Old MEs; yfp (m); p#theta_{bend}",250,-0.25,0.25,250,0.2,0.35);
  TH2D *hpthetabend_xpfp_old = new TH2D("hpthetabend_xpfp_old", "Old MEs; xpfp; p#theta_{bend}",250,-0.4,0.4,250,0.2,0.35);
  TH2D *hpthetabend_ypfp_old = new TH2D("hpthetabend_ypfp_old", "Old MEs; ypfp; p#theta_{bend}",250,-0.15,0.15,250,0.2,0.35);

  TH2D *hpthetabend_xptar_old = new TH2D("hpthetabend_xptar_old", "Old MEs; xptar; p#theta_{bend}",250,-0.4,0.4,250,0.2,0.35);
  TH2D *hpthetabend_yptar_old = new TH2D("hpthetabend_yptar_old", "Old MEs; yptar; p#theta_{bend}",250,-0.15,0.15,250,0.2,0.35);
  TH2D *hpthetabend_ytar_old = new TH2D("hpthetabend_ytar_old", "Old MEs; ytar; p#theta_{bend}",250,-0.15,0.15,250,0.2,0.35);

  TH2D *hdpel_xptar_old = new TH2D("hdpel_xptar_old", "Old MEs; xptar; p/p_{elastic}(#theta)-1",250,-0.4,0.4,250,-0.125,0.125);
  TH2D *hdpel_yptar_old = new TH2D("hdpel_yptar_old", "Old MEs; yptar; p/p_{elastic}(#theta)-1",250,-0.15,0.15,250,-0.125,0.125);
  TH2D *hdpel_ytar_old = new TH2D("hdpel_ytar_old", "Old MEs; ytar; p/p_{elastic}(#theta)-1",250,-0.15,0.15,250,-0.125,0.125);

  
  TH1D *hW_new = new TH1D("hW_new","New MEs;W (GeV);",250,Mp-0.2,Mp+0.2);
  TH1D *hdpel_new = new TH1D("hdpel_new","New MEs;p/p_{elastic}(#theta)-1;",250,-0.125,0.125);
  
  long nevent=0;

  //First pass: accumulate sums required for the fit: 
  
  while( C->GetEntry(elist->GetEntry(nevent++)) ){
    if( nevent % 1000 == 0 ) cout << nevent << endl;

    if( int(ntrack) == 1 ){
      //The first thing we want to do is to calculate the "true" electron momentum incident on BigBite:
      double Ebeam_corrected = ebeam - MeanEloss;
      
      double etheta = acos(pz[0]/p[0]);

      // Calculate the expected momentum of an elastically scattered electron at the reconstructed scattering angle and then correct it for the mean energy loss of the
      // electron on its way out of the target:
      double pelastic = Ebeam_corrected/(1.+Ebeam_corrected/Mp*(1.0-cos(etheta))); 
      
      double precon = p[0] + MeanEloss_outgoing; //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)

      double nu_recon = Ebeam_corrected - precon;
      
      double Q2recon = 2.0*Ebeam_corrected*precon*(1.0-cos(etheta));
      double W2recon = pow(Mp,2) + 2.0*Mp*nu_recon - Q2recon;

      double Wrecon = sqrt( std::max(0.0,W2recon) );
      
      double dpel = precon/pelastic-1.0;
      
      bool passed_HCAL_cut = true;
      
      if( usehcalcut != 0 ){
	passed_HCAL_cut = false; //go through the machinery to calculate dx and dy
	
      }

      if( passed_HCAL_cut && W2recon > 0.0 ){

	hdpel_old->Fill( dpel );
	hW_old->Fill( sqrt(W2recon) );

	hdpel_xfp_old->Fill( xfp[0], dpel );
	hdpel_yfp_old->Fill( yfp[0], dpel );
	hdpel_xpfp_old->Fill( thfp[0], dpel );
	hdpel_ypfp_old->Fill( phfp[0], dpel );

	
	hdpel_xptar_old->Fill( thtgt[0], dpel );
	hdpel_yptar_old->Fill( phtgt[0], dpel );
	hdpel_ytar_old->Fill( ytgt[0], dpel );
	
	
	double pincident = pelastic - MeanEloss_outgoing;
	
	//Now we need to calculate the "true" trajectory bend angle for the electron from the reconstructed angles:
	TVector3 enhat_tgt( thtgt[0], phtgt[0], 1.0 );
	enhat_tgt = enhat_tgt.Unit();
	
	TVector3 enhat_fp( thfp[0], phfp[0], 1.0 );
	enhat_fp = enhat_fp.Unit();

        TVector3 GEMzaxis(-sin(GEMpitch),0,cos(GEMpitch));
	TVector3 GEMyaxis(0,1,0);
	TVector3 GEMxaxis = (GEMyaxis.Cross(GEMzaxis)).Unit();
	
	TVector3 enhat_fp_rot = enhat_fp.X() * GEMxaxis + enhat_fp.Y() * GEMyaxis + enhat_fp.Z() * GEMzaxis;

	double thetabend = acos( enhat_fp_rot.Dot( enhat_tgt ) );
	
	// cout << "Reconstructed bend angle = " << thetabend * TMath::RadToDeg() << endl;
	// cout << "p*thetabend = " << pincident * thetabend << endl;

	//Increment sums only if dpel and W are within limits:

	if( dpel >= dpelmin_fit && dpel <= dpelmax_fit &&
	    Wrecon >= Wmin_fit && Wrecon <= Wmax_fit ){
	  //cout << "dpel, W = " << dpel << ", " << Wrecon << endl;

	  hp_old->Fill( precon );
	  hpthetabend_old->Fill( pincident * thetabend );

	  hpthetabend_xfp_old->Fill( xfp[0], pincident*thetabend );
	  hpthetabend_yfp_old->Fill( yfp[0], pincident*thetabend );
	  hpthetabend_xpfp_old->Fill( thfp[0], pincident*thetabend );
	  hpthetabend_ypfp_old->Fill( phfp[0], pincident*thetabend );

	  hpthetabend_xptar_old->Fill( thtgt[0], pincident*thetabend );
	  hpthetabend_yptar_old->Fill( phtgt[0], pincident*thetabend );
	  hpthetabend_ytar_old->Fill( ytgt[0], pincident*thetabend );
	  
	  hxpfp_vs_xfp_old->Fill( xfp[0], thfp[0] );
	  
	  vector<double> term(nparams);

	  int ipar=0;
	  
	  double xtar = 0.0;
	  
	  for(int i=0; i<=order; i++){
	    for(int j=0; j<=order-i; j++){
	      for(int k=0; k<=order-i-j; k++){
		for(int l=0; l<=order-i-j-k; l++){
		  for(int m=0; m<=order-i-j-k-l; m++){
		    term[ipar] = pow(xfp[0],m)*pow(yfp[0],l)*pow(thfp[0],k)*pow(phfp[0],j)*pow(xtar,i);
		    
		    b_pth(ipar) += term[ipar]*pincident*thetabend; //use the "true" momentum for an elastically scattered electron 
		    
		    ipar++;
		  }
		}
	      }
	    }
	  }
	  
	  for( ipar=0; ipar<nparams; ipar++ ){
	    for( int jpar=0; jpar<nparams; jpar++ ){
	      M(ipar,jpar) += term[ipar]*term[jpar];
	    }
	  }
	  
	}
      }
      
    }
  }

  //For now, we set all matrix elements involving non-zero powers of x target to zero:
  if( true ){
    for( int ipar=0; ipar<nparams; ipar++ ){
      if( xtar_expon[ipar] != 0 ){
	M(ipar,ipar) = 1.0;
	b_pth(ipar) = 0.0;
	for( int jpar=0; jpar<nparams; jpar++ ){
	  if( jpar != ipar ) M(ipar,jpar) = 0.0;
	}
      }
    }
  }

  TDecompSVD A_pth(M);

  bool goodfit = A_pth.Solve( b_pth );
  cout << "Momentum fit done, success = " << goodfit << endl;

  cout << "New p*thetabend coefficients = " << endl;
  
  for( int ipar=0; ipar<nparams; ipar++ ){
    cout << "Coefficient, exponents = " << b_pth(ipar) << " " << xfp_expon[ipar] << "," << yfp_expon[ipar] << "," << xpfp_expon[ipar] << "," << ypfp_expon[ipar]
	 << "," << xtar_expon[ipar] << endl;
  }


  //Do a second loop over the data to fill dp new and W new histograms plus correlation plots:

  cout << "Doing a second loop to plot results with new coefficients:" << endl;
  
  nevent = 0;
  while( C->GetEntry( elist->GetEntry( nevent++ ) ) ){
    if( nevent % 1000 == 0 ) cout << nevent << endl;

    if( int(ntrack) == 1 ){
      double pthetabend_sum = 0.0;
      double xtar = 0.0;

      int ipar=0;
      
      for(int i=0; i<=order; i++){
	for(int j=0; j<=order-i; j++){
	  for(int k=0; k<=order-i-j; k++){
	    for(int l=0; l<=order-i-j-k; l++){
	      for(int m=0; m<=order-i-j-k-l; m++){
		double term = pow(xfp[0],m)*pow(yfp[0],l)*pow(thfp[0],k)*pow(phfp[0],j)*pow(xtar,i);

		pthetabend_sum += term * b_pth( ipar );

		ipar++;
	      }
	    }
	  }
	}
      }

      //Now we need to calculate the "true" trajectory bend angle for the electron from the reconstructed angles:
      TVector3 enhat_tgt( thtgt[0], phtgt[0], 1.0 );
      enhat_tgt = enhat_tgt.Unit();
	
      TVector3 enhat_fp( thfp[0], phfp[0], 1.0 );
      enhat_fp = enhat_fp.Unit();

      TVector3 GEMzaxis(-sin(GEMpitch),0,cos(GEMpitch));
      TVector3 GEMyaxis(0,1,0);
      TVector3 GEMxaxis = (GEMyaxis.Cross(GEMzaxis)).Unit();
	
      TVector3 enhat_fp_rot = enhat_fp.X() * GEMxaxis + enhat_fp.Y() * GEMyaxis + enhat_fp.Z() * GEMzaxis;

      //Calculate trajectory bend angle the same way as we did above:
      double thetabend = acos( enhat_fp_rot.Dot( enhat_tgt ) );

      double pnew = pthetabend_sum/thetabend;

      double Ebeam_corrected = ebeam - MeanEloss;
      
      double etheta = acos(pz[0]/p[0]);

      // Calculate the expected momentum of an elastically scattered electron at the reconstructed scattering angle and then correct it for the mean energy loss of the
      // electron on its way out of the target:
      double pelastic = Ebeam_corrected/(1.+Ebeam_corrected/Mp*(1.0-cos(etheta)));
      double pnew_corr = pnew + MeanEloss_outgoing;

      double nu = Ebeam_corrected - pnew_corr;
      double Q2 = 2.0*Ebeam_corrected * pnew_corr*(1.0-cos(etheta));
      double W2 = pow(Mp,2) + 2.0*Mp*nu - Q2;

      double Wnew = sqrt(std::max(0.0,W2));

      double dpelnew = pnew_corr/pelastic - 1.0;

      hdpel_new->Fill( dpelnew );
      hW_new->Fill( Wnew );
      
    }
    
  }
  
  elist->Delete();

  fout->Write();
}
