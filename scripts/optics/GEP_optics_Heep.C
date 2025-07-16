#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TCut.h"
#include "TTreeFormula.h"
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>
#include "TMinuit.h"
#include "gep_tree_data.C"
#include <vector>
#include <set>
#include <map>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRotation.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TMath.h"
//Now the theory behind this alignment is to take a clean sample of elastic ep events and, ASSUMING that the ECAL position and the SBS optics model are correct (including vertex reconstruction), determine the best set of translational and rotational offsets of the front tracker in order to optimize the proton's kinematic reconstruction. Offsets will include
// x,y,z offsets of FT
// ax,ay,az angles of FT
// beam energy offset
// In theory this procedure should not require enormous statistics

//parameters of old optics model:

//offsets will be considered RELATIVE to existing database values:
int NMAX = 100000;

double EBEAM0 = 10.688;

double Z0targ = -0.09; //meters

TVector3 ECALxaxis,ECALyaxis,ECALzaxis;
TVector3 SBSxaxis,SBSyaxis,SBSzaxis;
TVector3 ECALPOS;

double ECALtheta, SBStheta;
double ECALdist;
const double PI = TMath::Pi();
const double Mp = 0.93827208816;

const double rho_LH2 = 0.072; //g/cc
const double dEdx_LH2 = 0.0059; // dEdx in GeV * cm^2/g (it is approximately constant in the energy range relevant to GEP):

vector<double> xptar_coeff, yptar_coeff, ytar_coeff, pth_coeff;
vector<int> xfp_expon,yfp_expon,xpfp_expon,ypfp_expon,xtar_expon;

//Declare arrays to hold relevant event variables here:
vector<double> XFP,YFP,THFP,PHFP,RXFP,RYFP,RTHFP,RPHFP;
vector<double> THTAR,PHTAR,XTAR,YTAR,P;
vector<double> ZVERTEX;
vector<double> xECAL, yECAL, ECALO, TECAL;

void CHI2_FCN( int &npar, double *gin, double &f, double *par, int flag ){
  double chi2 = 0.0;

  //In this script, the GEM position is NOT defined globally, but relative to the
  //ideal optics origin, which is why we use the "r_*" tracking variables
  double GEMX0 = par[0];
  double GEMY0 = par[1];
  double GEMZ0 = par[2];
  double GEMax = par[3];
  double GEMay = par[4];
  double GEMaz = par[5];

  double dE = par[6];

  double dxECAL = par[7];
  double dyECAL = par[8];
  double dzECAL = par[9];
  
  TVector3 GEMPOS(GEMX0,GEMY0,GEMZ0);
  
  for( int i=0; i<XFP.size(); i++ ){
    TVector3 trackpos_local( RXFP[i], RYFP[i], 0.0 );
    TVector3 trackdir_local( RTHFP[i], RPHFP[i], 1.0 );
    trackdir_local = trackdir_local.Unit();

    TRotation R;
    R.RotateX(GEMax);
    R.RotateY(GEMay);
    R.RotateZ(GEMaz);
    
    TVector3 trackdir_global = R * trackdir_local;
    TVector3 trackpos_global = GEMPOS + R * trackpos_local;

    double xpfprot = trackdir_global.X()/trackdir_global.Z();
    double ypfprot = trackdir_global.Y()/trackdir_global.Z();

    double xfprot = trackpos_global.X() - xpfprot * trackpos_global.Z();
    double yfprot = trackpos_global.Y() - ypfprot * trackpos_global.Z();

    double xtar = 0.0;
    double xptar = 0.0, yptar = 0.0, ytar = 0.0, pthetabend = 0.0;
    for( int iterm=0; iterm<xptar_coeff.size(); iterm++ ){
      double term = pow(xfprot, xfp_expon[iterm])*pow(yfprot, yfp_expon[iterm])*pow(xpfprot, xpfp_expon[iterm])*pow(ypfprot,ypfp_expon[iterm])*pow(xtar,xtar_expon[iterm]);

      xptar += xptar_coeff[iterm]*term;
      yptar += yptar_coeff[iterm]*term;
      ytar += ytar_coeff[iterm]*term;
      pthetabend += pth_coeff[iterm]*term;
    }

    // ytar = vz * sin(theta) - yptar * vz * cos(theta)
    // vz = ytar/(sin(theta)-yptar*Cos(theta))
    
    double vz = ytar / (sin(SBStheta)-yptar*cos(SBStheta));
    
    TVector3 targdir(xptar,yptar,1.0);
    targdir = targdir.Unit();

    //For SBS, we know that the z axes of the ideal transport system at both target and FP coincide, so the bend angle calculation is simple:
    double thetabend = acos( targdir.Dot( trackdir_global ) );
    double p = pthetabend / thetabend;
    
    TVector3 Pproton_transport = p * targdir;
    TVector3 Pproton_hall = Pproton_transport.X() * SBSxaxis +
      Pproton_transport.Y() * SBSyaxis + Pproton_transport.Z() * SBSzaxis;

    double ptheta = Pproton_hall.Theta();
    double pphi = Pproton_hall.Phi();
    
    TVector3 epos_hall = ECALPOS + (xECAL[i]+dxECAL) * ECALxaxis + (yECAL[i]+dyECAL) * ECALyaxis + dzECAL * ECALzaxis;

    TVector3 vertex(0,0,vz);

    TVector3 edir_hall = (epos_hall - vertex).Unit();

    double etheta = edir_hall.Theta();
    double ephi = edir_hall.Phi();
    
    double ebeam = EBEAM0 + dE;

    //let's correct for beam energy loss in LH2 prior to the scattering and proton energy loss exiting the target:
    double thick_LH2_gcm2 = rho_LH2 * std::max(0.0,std::min(100.0*(vz-Z0targ),30.0));

    //Mean collisional energy loss before the scattering:
    double eloss_mean = dEdx_LH2 * thick_LH2_gcm2; 
    double ebeam_corr = ebeam - eloss_mean;
    
    double Eprime_eth = ebeam_corr/(1.0 + ebeam_corr/Mp * (1.0-cos(etheta)));
    double nu_eth = ebeam_corr - Eprime_eth;
    double Q2_eth = 2.0*Mp*nu_eth;
    double tau_eth = Q2_eth/(4.0*pow(Mp,2));
    double pp_eth = sqrt(Q2_eth*(1.0+tau_eth));

    double ptheta_eth = acos( (ebeam_corr-Eprime_eth*cos(etheta))/pp_eth );
    double pphi_eph = ephi + PI;

    TVector3 pnhat_hall_eangles(sin(ptheta_eth)*cos(pphi_eph),sin(ptheta_eth)*sin(pphi_eph),cos(ptheta_eth));

    TVector3 pnhat_transport_eangles( pnhat_hall_eangles.Dot(SBSxaxis),
				      pnhat_hall_eangles.Dot(SBSyaxis),
				      pnhat_hall_eangles.Dot(SBSzaxis) );

    double thtar_e = pnhat_transport_eangles.X()/pnhat_transport_eangles.Z();
    double phtar_e = pnhat_transport_eangles.Y()/pnhat_transport_eangles.Z();

    double pp_ptheta = 2.0*Mp*ebeam_corr*(Mp+ebeam_corr)*cos(ptheta)/(pow(Mp,2)+2.0*Mp*ebeam_corr + pow(ebeam_corr*sin(ptheta),2));

    // Shall we try to correct the proton momentum for eloss in
    // the target? We don't do this when we generate the optics model, and the
    // correction is small and roughly constant. So let's ignore.
    
    double dpp = p/pp_ptheta - 1.;
    
    chi2 += pow((xptar-thtar_e)/0.001,2) + pow( (yptar-phtar_e)/0.001, 2 ) +
      pow( dpp/0.008, 2 );
  }

  f = chi2;
}

void GEP_optics_Heep( const char *configfilename, const char *outfilename="GEPheep.root" ){
  ifstream infile(configfilename);
  if( !infile ) return;

 
  
  TString currentline;

  TChain *C = new TChain("T");
  
  while( currentline.ReadLine(infile) && !(currentline.BeginsWith("endlist")) ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline.Data());
    }
  }

  TCut globalcut = "";
  while( currentline.ReadLine(infile) && !(currentline.BeginsWith("endcut")) ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
    }
  }

  TString oldopticsfname="";
  
  //  double ebea = 10.688;  

  double ECALX0=0,ECALY0=0;
  
  double dpp0 = 0, dppsigma = 0.01;
  double dx0 = 0, dxsigma = 0.015;
  double dy0 = 0, dysigma = 0.021;
  double dt0 = 15, dtsigma = 5.0;
  double EoverP_cent = 1.0, EoverP_sigma = 0.1;
  
  while( currentline.ReadLine(infile) && !(currentline.BeginsWith("endconfig")) ){
    if( !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      if( tokens->GetEntries() >= 2 ){
	TString skey = ( (TObjString*) (*tokens)[0] )->GetString();
	TString sval = ( (TObjString*) (*tokens)[1] )->GetString();

	if( skey == "EoverPcent" ){
	  EoverP_cent = sval.Atof();
	}

	if( skey == "EoverPsigma" ){
	  EoverP_sigma = sval.Atof();
	}
	
	if( skey == "dpp0" ){
	  dpp0 = sval.Atof();
	}
	if( skey == "dx0" ){
	  dx0 = sval.Atof();
	}

	if( skey == "dy0" ){
	  dy0 = sval.Atof();
	}

	if( skey == "dt0" ){
	  dt0 = sval.Atof();
	}

	if( skey == "dppsigma" ){
	  dppsigma = sval.Atof();
	}
	if( skey == "dxsigma" ){
	  dxsigma = sval.Atof();
	}

	if( skey == "dysigma" ){
	  dysigma = sval.Atof();
	}

	if( skey == "dtsigma" ){
	  dtsigma = sval.Atof();
	}
	
	
	if( skey == "oldopticsfile" ){
	  oldopticsfname = sval;
	}

	if( skey == "ebeam" ){
	  EBEAM0 = sval.Atof();
	}

	if( skey == "ECALtheta" ){
	  ECALtheta = sval.Atof() * PI/180.0;
	}

	if( skey == "SBStheta" ){
	  SBStheta = sval.Atof() * PI/180.0;
	}

	if( skey == "ECALdist" ){
	  ECALdist = sval.Atof();
	}

	if( skey == "ECALX0" ){
	  ECALX0 = sval.Atof();
	}

	if( skey == "ECALY0" ){
	  ECALY0 = sval.Atof();
	}
      }
    }
  }

  ECALzaxis.SetXYZ(sin(ECALtheta),0,cos(ECALtheta));
  ECALxaxis.SetXYZ(0,-1,0);
  ECALyaxis = ECALzaxis.Cross(ECALxaxis).Unit();

  ECALPOS = ECALzaxis * ECALdist + ECALxaxis * ECALX0 + ECALyaxis * ECALY0;
  
  SBSzaxis.SetXYZ(-sin(SBStheta),0,cos(SBStheta));
  SBSxaxis.SetXYZ(0,-1,0);
  SBSyaxis = SBSzaxis.Cross(SBSxaxis).Unit();
  
  C->SetBranchStatus("*",0);
  C->SetBranchStatus("heep.*",1);
  C->SetBranchStatus("sbs.tr.*",1);
  C->SetBranchStatus("sbs.gemFT.track.*",1);
  C->SetBranchStatus("sbs.gemFPP.track.*",1);
  C->SetBranchStatus("earm.ecal.x",1);
  C->SetBranchStatus("earm.ecal.y",1);
  C->SetBranchStatus("earm.ecal.e",1);
  C->SetBranchStatus("earm.ecal.atimeblk",1);

  //we may add others as we go. 
  
  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut, C );

  gep_tree_data *T = new gep_tree_data(C);

  if( oldopticsfname == "" ) return;

  ifstream oldopticsfile(oldopticsfname);

  int nterms;

  oldopticsfile >> nterms;

  xptar_coeff.resize(nterms);
  yptar_coeff.resize(nterms);
  ytar_coeff.resize(nterms);
  pth_coeff.resize(nterms);

  xfp_expon.resize(nterms);
  yfp_expon.resize(nterms);
  xpfp_expon.resize(nterms);
  ypfp_expon.resize(nterms);
  xtar_expon.resize(nterms);
  
  for( int i=0; i<nterms; i++ ){
    oldopticsfile >> xptar_coeff[i] >> yptar_coeff[i] >> ytar_coeff[i] >> pth_coeff[i]
		  >> xfp_expon[i] >> yfp_expon[i] >> xpfp_expon[i] >> ypfp_expon[i] >> xtar_expon[i];
  }

  long nevent = 0;

  int treenum=-1, oldtreenum=treenum;

  while( C->GetEntry(nevent++) ){
    treenum = C->GetTreeNumber();
    if( treenum != oldtreenum ){
      oldtreenum = treenum;
      GlobalCut->UpdateFormulaLeaves();
    }
    if( nevent % 1000 == 0 ) cout << "event " << nevent << ", file = " << C->GetFile()->GetName() << endl;

    bool passedcut = GlobalCut->EvalInstance(0) != 0;

    if( passedcut && XFP.size() < NMAX ){
      XFP.push_back( T->sbs_tr_x[0] );
      YFP.push_back( T->sbs_tr_y[0] );
      THFP.push_back( T->sbs_tr_th[0] );
      PHFP.push_back( T->sbs_tr_ph[0] );

      RXFP.push_back( T->sbs_tr_r_x[0] );
      RYFP.push_back( T->sbs_tr_r_y[0] );
      RTHFP.push_back( T->sbs_tr_r_th[0] );
      RPHFP.push_back( T->sbs_tr_r_ph[0] );

      THTAR.push_back( T->sbs_tr_tg_th[0] );
      PHTAR.push_back( T->sbs_tr_tg_ph[0] );
      XTAR.push_back( 0.0 );
      YTAR.push_back( T->sbs_tr_tg_y[0] );
      P.push_back( T->sbs_tr_p[0] );

      ZVERTEX.push_back( T->sbs_tr_vz[0] );
      
      xECAL.push_back( T->earm_ecal_x );
      yECAL.push_back( T->earm_ecal_y );
      ECALO.push_back( T->earm_ecal_e );
      TECAL.push_back( T->earm_ecal_atimeblk );
      
    }

    if( XFP.size() == NMAX ) break;
   
  }

  cout << "N tracks passing cut = " << XFP.size() << endl;
  
  //number of parameters = 6 + ??
  // FT GEM x,y,z, ax,ay,az,
  // (up to) 10 parameters: FT GEM x,y,z, ax,ay,az, ebeam, ECALx,y,z
  TMinuit *FitFunc = new TMinuit( 10 ); 

  FitFunc->SetFCN( CHI2_FCN );

  int ierflg = 0;
  
  FitFunc->mnparm( 0, "GEMX0", 0, 0.01,0,0,ierflg);
  FitFunc->mnparm( 1, "GEMY0", 0, 0.001,0,0,ierflg);
  FitFunc->mnparm( 2, "GEMZ0", 0, 0.015,0,0,ierflg);
  FitFunc->mnparm( 3, "GEMax", 0, 0.002,0,0,ierflg);
  FitFunc->mnparm( 4, "GEMay", 0, 0.002,0,0,ierflg);
  FitFunc->mnparm( 5, "GEMaz", 0, 0.002,0,0,ierflg);
  FitFunc->mnparm( 6, "dEbeam", 0, 0.0005*EBEAM0,0,0,ierflg);
  FitFunc->mnparm( 7, "ECALX0", 0, 0.005, 0,0,ierflg);
  FitFunc->mnparm( 8, "ECALY0", 0, 0.005, 0,0,ierflg);
  FitFunc->mnparm( 9, "ECALZ0", 0, 0.01, 0,0,ierflg);

  FitFunc->FixParameter(1);
  //  FitFunc->FixParameter(2);
  FitFunc->FixParameter(3);
  FitFunc->FixParameter(6);
  FitFunc->FixParameter(7);
  FitFunc->FixParameter(8);
  FitFunc->FixParameter(9);
  
  double arglist[10];
  arglist[0] = 1;
  FitFunc->mnexcm( "SET ERR", arglist, 1, ierflg );
  arglist[0] = 5000;
  arglist[1] = 1.;
  
  FitFunc->mnexcm("MIGRAD",arglist,2,ierflg);

  double GEMX0,GEMY0,GEMZ0,dGEMX0,dGEMY0,dGEMZ0;
  double GEMax,GEMay,GEMaz,dGEMax,dGEMay,dGEMaz;
  
  FitFunc->GetParameter(0, GEMX0, dGEMX0);
  FitFunc->GetParameter(1, GEMY0, dGEMY0);
  FitFunc->GetParameter(2, GEMZ0, dGEMZ0);
  FitFunc->GetParameter(3, GEMax, dGEMax);
  FitFunc->GetParameter(4, GEMay, dGEMay);
  FitFunc->GetParameter(5, GEMaz, dGEMaz);

  TVector3 GEMPOS(GEMX0,GEMY0,GEMZ0);
  
  nevent = 0;

  TFile *fout = new TFile(outfilename,"RECREATE");
  
  TH1D *hdppold = new TH1D("hdppold","OLD;p_{p}/p_{el}(#theta_{p})-1;",250,-0.25,0.25);
  TH1D *hdppnew = new TH1D("hdppnew","NEW;p_{p}/p_{el}(#theta_{p})-1;",250,-0.25,0.25);

  TH2D *hdxdyECAL_old = new TH2D("hdxdyECAL_old","OLD;#Deltay (m);#Deltax (m)",250,-0.25,0.25,250,-0.25,0.25);
  TH2D *hdxdyECAL_new = new TH2D("hdxdyECAL_new","NEW;#Deltay (m);#Deltax (m)",250,-0.25,0.25,250,-0.25,0.25);

  TH1D *hdthtar_old = new TH1D("hdthtar_old","OLD ;#Delta#theta_{tgt} (rad);",150,-0.03,0.03);
  TH1D *hdphtar_old = new TH1D("hdphtar_old","OLD ;#Delta#phi_{tgt} (rad);",150,-0.03,0.03);

  TH1D *hdthtar_new = new TH1D("hdthtar_new","NEW ;#Delta#theta_{tgt} (rad);",150,-0.03,0.03);
  TH1D *hdphtar_new = new TH1D("hdphtar_new","NEW ;#Delta#phi_{tgt} (rad);",150,-0.03,0.03);

  TTree *Tout = new TTree("Tout","GEP heep optics");
  //Basic kinematics, old and new:
  double T_xfp_old, T_yfp_old, T_thfp_old, T_phfp_old;
  double T_xfp_new, T_yfp_new, T_thfp_new, T_phfp_new;
  double T_thtar_old, T_phtar_old, T_ytar_old, T_p_old;
  double T_thtar_new, T_phtar_new, T_ytar_new, T_p_new;
  double T_vz_old, T_vz_new, T_xECAL, T_yECAL;

  //Physics angles and calorimeter energy:
  double T_ptheta_old, T_pphi_old, T_Ecalo;
  double T_ptheta_new, T_pphi_new;

  //Exclusivity cut variables:
  double T_dpp_old, T_dpp_new;
  double T_dthtar_old, T_dthtar_new;
  double T_dphtar_old, T_dphtar_new;
  double T_dxECAL_old, T_dxECAL_new;
  double T_dyECAL_old, T_dyECAL_new;

  Tout->Branch("xfpold",&T_xfp_old);
  Tout->Branch("yfpold",&T_yfp_old);
  Tout->Branch("thfpold",&T_thfp_old);
  Tout->Branch("phfpold",&T_phfp_old);
  Tout->Branch("xfpnew",&T_xfp_new);
  Tout->Branch("yfpnew",&T_yfp_new);
  Tout->Branch("thfpnew",&T_thfp_new);
  Tout->Branch("phfpnew",&T_phfp_new);

  Tout->Branch("thtarold",&T_thtar_old);
  Tout->Branch("phtarold",&T_phtar_old);
  Tout->Branch("ytarold",&T_ytar_old);
  Tout->Branch("pold",&T_p_old);
  
  Tout->Branch("thtarnew",&T_thtar_new);
  Tout->Branch("phtarnew",&T_phtar_new);
  Tout->Branch("ytarnew",&T_ytar_new);
  Tout->Branch("pnew",&T_p_new);

  Tout->Branch("vzold",&T_vz_old);
  Tout->Branch("vznew",&T_vz_new);
  Tout->Branch("xECAL",&T_xECAL);
  Tout->Branch("yECAL",&T_yECAL);

  Tout->Branch("ECALO",&T_Ecalo);
  Tout->Branch("pthetaold",&T_ptheta_old);
  Tout->Branch("pphiold",&T_pphi_old);
  Tout->Branch("pthetanew",&T_ptheta_new);
  Tout->Branch("pphinew",&T_pphi_new);
  
  
  
  treenum = 0; oldtreenum = -1;
  while( C->GetEntry(nevent++) ){
    treenum = C->GetTreeNumber();
    if( treenum != oldtreenum ){
      oldtreenum = treenum;
      GlobalCut->UpdateFormulaLeaves();
    }
    if( nevent % 1000 == 0 ) cout << "event " << nevent << ", file = " << C->GetFile()->GetName() << endl;
    
    if( T->sbs_tr_n > 0 && fabs(T->heep_dt_ADC-dt0)<3.0*dtsigma &&
	fabs( T->sbs_tr_vz[0] + Z0targ ) <= 0.18 &&
	T->heep_ecalo/T->heep_eprime_eth > EoverP_cent - 3.0*EoverP_sigma &&
	T->sbs_gemFPP_track_ntrack > 0 ){

      if( fabs( T->heep_dpp - dpp0 )<=3.0*dppsigma ){
	hdxdyECAL_old->Fill( T->heep_dyECAL, T->heep_dxECAL );
      }
      
      if( sqrt(pow((T->heep_dxECAL-dx0)/dxsigma, 2)+pow((T->heep_dyECAL-dy0)/dysigma,2))<=3.0 ){
	hdppold->Fill( T->heep_dpp );
      }

      //Calculate "old" delta thetatarget and delta phi target:
      TVector3 epos = ECALPOS + T->earm_ecal_x * ECALxaxis + T->earm_ecal_y * ECALyaxis;

      //      TVector3 vertex_old(0,0,T->sbs_tr_vz[0]);

      //TVector3 edir_old = (epos - vertex_old).Unit();

      double ptheta_eth_old = T->heep_pth_eth;
      double pphi_eph_old = T->heep_ephi + PI;

      TVector3 pdir_hall_old(sin(ptheta_eth_old)*cos(pphi_eph_old),
			     sin(ptheta_eth_old)*sin(pphi_eph_old),
			     cos(ptheta_eth_old) );

      TVector3 pdir_transport_old( pdir_hall_old.Dot(SBSxaxis),
				   pdir_hall_old.Dot(SBSyaxis),
				   pdir_hall_old.Dot(SBSzaxis) );

      double pthtar_e_old = pdir_transport_old.X()/pdir_transport_old.Z();
      double pphtar_e_old = pdir_transport_old.Y()/pdir_transport_old.Z();

      hdthtar_old->Fill( T->sbs_tr_tg_th[0]-pthtar_e_old );
      hdphtar_old->Fill( T->sbs_tr_tg_ph[0]-pphtar_e_old );


      TVector3 trackpos_local( T->sbs_tr_r_x[0], T->sbs_tr_r_y[0], 0.0 );
      TVector3 trackdir_local( T->sbs_tr_r_th[0], T->sbs_tr_r_ph[0], 1.0 );
      trackdir_local = trackdir_local.Unit();
      
      TRotation R;
      R.RotateX(GEMax);
      R.RotateY(GEMay);
      R.RotateZ(GEMaz);
      
      TVector3 trackpos_global = GEMPOS + R * trackpos_local;
      TVector3 trackdir_global = R * trackdir_local;
      
      double xpfprot = trackdir_global.X()/trackdir_global.Z();
      double ypfprot = trackdir_global.Y()/trackdir_global.Z();
      
      double xfprot = trackpos_global.X() - xpfprot * trackpos_global.Z();
      double yfprot = trackpos_global.Y() - ypfprot * trackpos_global.Z();

      //Now calculate "new" exclusivity variables:
      double xptar_new = 0.0;
      double yptar_new = 0.0;
      double ytar_new = 0.0;
      double pthetabend_new = 0.0;
      double xtar = 0.0;
      
      for( int iterm=0; iterm<xptar_coeff.size(); iterm++ ){
	double term = pow(xfprot, xfp_expon[iterm])*pow(yfprot, yfp_expon[iterm])*pow(xpfprot, xpfp_expon[iterm])*pow(ypfprot,ypfp_expon[iterm])*pow(xtar,xtar_expon[iterm]);
	
	xptar_new += xptar_coeff[iterm]*term;
	yptar_new += yptar_coeff[iterm]*term;
	ytar_new += ytar_coeff[iterm]*term;
	pthetabend_new += pth_coeff[iterm]*term;
      }

      double vz_new = ytar_new/(sin(SBStheta)-yptar_new*cos(SBStheta));
      TVector3 targdir_new(xptar_new,yptar_new,1.0);
      targdir_new = targdir_new.Unit();
      double thetabend_new = acos( targdir_new.Dot( trackdir_global ) );
      double p_new = pthetabend_new/thetabend_new;

      TVector3 Pproton_transport_new = p_new * targdir_new;
      TVector3 Pproton_hall_new = Pproton_transport_new.X() * SBSxaxis +
	Pproton_transport_new.Y() * SBSyaxis +
	Pproton_transport_new.Z() * SBSzaxis;
	
      double ptheta_new = Pproton_hall_new.Theta();
      double pphi_new = Pproton_hall_new.Phi();

      double thick_LH2_gcm2 = rho_LH2 * std::max(0.0,std::min(100.0*(vz_new-Z0targ),30.0));
      double eloss_new = dEdx_LH2 * thick_LH2_gcm2;
      double ebeam_corr_new = EBEAM0 - eloss_new;

      double pp_ptheta_new = 2.0*Mp*ebeam_corr_new*(Mp+ebeam_corr_new)*cos(ptheta_new)/(pow(Mp,2)+2.0*Mp*ebeam_corr_new + pow(ebeam_corr_new*sin(ptheta_new),2));
      double dpp_new = p_new/pp_ptheta_new - 1.0;
      hdppnew->Fill( dpp_new );
      //      hdthtar_new->Fill( xptar_new - 

      TVector3 vertex_new(0,0,vz_new);
      TVector3 edir_new = (epos-vertex_new).Unit();

      double etheta_new = edir_new.Theta();
      double ephi_new = edir_new.Phi();
      
      double pphi_eph_new = ephi_new + PI;

      double Eprime_eth_new = ebeam_corr_new / (1.0 + ebeam_corr_new/Mp*(1.0-cos(etheta_new)));
      double nu_new = ebeam_corr_new - Eprime_eth_new;
      double Q2_eth_new = 2.0*Mp*nu_new;

      double tau_eth_new = Q2_eth_new/(4.0*pow(Mp,2));

      double pp_eth_new = sqrt(Q2_eth_new*(1.0 + tau_eth_new));

      double ptheta_eth_new = acos( (ebeam_corr_new-Eprime_eth_new*cos(etheta_new))/pp_eth_new);

      TVector3 pdir_hall_new( sin(ptheta_eth_new)*cos(pphi_eph_new),
			      sin(ptheta_eth_new)*sin(pphi_eph_new),
			      cos(ptheta_eth_new) );

      TVector3 pdir_transport_new( pdir_hall_new.Dot( SBSxaxis ),
				   pdir_hall_new.Dot( SBSyaxis ),
				   pdir_hall_new.Dot( SBSzaxis ) );

      double pthtar_e_new = pdir_transport_new.X()/pdir_transport_new.Z();
      double pphtar_e_new = pdir_transport_new.Y()/pdir_transport_new.Z();

      hdthtar_new->Fill( xptar_new - pthtar_e_new );
      hdphtar_new->Fill( yptar_new - pphtar_e_new );
      
      
      
    }

  }
  fout->Write();
}
