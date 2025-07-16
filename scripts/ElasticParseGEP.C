#include "gep_tree_data.C"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include "TCut.h"
#include "TTreeFormula.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TLorentzVector.h"

//Generic elastic parsing script; baseline for physics analysis.
//Initial purpose is quick-and-dirty elastic event selection and optics calibration; read root trees, create output tree with basic, loos cuts applied:

const double PI = TMath::Pi();
const double Mp = 0.938272;
const double mu_p = 2.793;

void ElasticParseGEP(const char *configfilename, const char *outfilename="temp.root", double sbsthetadeg=25.7){

  double THSBS = sbsthetadeg*PI/180.0;

  TVector3 SBSzaxis( -sin(THSBS),0,cos(THSBS) );
  TVector3 SBSxaxis(0,-1,0);
  TVector3 SBSyaxis = SBSzaxis.Cross(SBSxaxis).Unit();
  
  ifstream infile(configfilename);

  if( !infile ){
    return;
  } else {
    TString currentline;

    TChain *C = new TChain("T");
    TFile *fout = new TFile(outfilename,"RECREATE");
    
    while( currentline.ReadLine( infile ) && !currentline.BeginsWith("endlist")){
      if( !currentline.BeginsWith("#") ){
	C->Add(currentline.Data());
      }
    }

    TCut globalcut = "";
    
    while( currentline.ReadLine( infile ) && !currentline.BeginsWith("endcut")){
      if( !currentline.BeginsWith("#") ){
	globalcut += currentline.Data();
      }
    }

    fout->cd();

    TTree *Tout = new TTree("Tout", "GEP elastic parsing");

    int runnum;
    double etheta,ephi,ptheta,pphi,ecalo,pp,eprime_eth,pp_eth,pth_eth,pp_pth,eprime_pth,eth_pth,eprime_pp,eth_pp,Q2_pp,Q2_eth,Q2_pth,Q2_p4vect,Q2_e4vect,eps_pp, eps_eth,eps_pth,eps_4vect,K_eth,K_pth,K_pp,K_p4vect,dpp,dpe,dphi,acopl,inel,dxECAL,dyECAL,dxECAL4vect,dyECAL4vect,dtADC;

    Tout->Branch("runnum",&runnum);
    Tout->Branch("etheta",&etheta);
    Tout->Branch("ephi",&ephi);
    Tout->Branch("ptheta",&ptheta);
    Tout->Branch("pphi",&pphi);
    Tout->Branch("ecalo",&ecalo);
    Tout->Branch("pp",&pp);
    Tout->Branch("eprime_eth",&eprime_eth);
    Tout->Branch("pp_eth",&pp_eth);
    Tout->Branch("pp_pth",&pp_pth);
    Tout->Branch("eprime_pth",&eprime_pth);
    Tout->Branch("eprime_pp",&eprime_pp);
    Tout->Branch("eth_pp",&eth_pp);
    Tout->Branch("Q2_pp",&Q2_pp);
    Tout->Branch("Q2_eth",&Q2_eth);
    Tout->Branch("Q2_pth",&Q2_pth);
    Tout->Branch("Q2_e4vect",&Q2_e4vect);
    Tout->Branch("eps_pp",&eps_pp);
    Tout->Branch("eps_eth",&eps_eth);
    Tout->Branch("eps_pth",&eps_pth);
    Tout->Branch("eps_4vect",&eps_4vect);
    Tout->Branch("K_eth",&K_eth);
    Tout->Branch("K_pth",&K_pth);
    Tout->Branch("K_pp",&K_pp);
    Tout->Branch("K_p4vect",&K_p4vect);
    Tout->Branch("dpp",&dpp);
    Tout->Branch("dpe",&dpe);
    Tout->Branch("dphi",&dphi);
    Tout->Branch("acopl",&acopl);
    Tout->Branch("inel",&inel);
    Tout->Branch("dxECAL",&dxECAL);
    Tout->Branch("dyECAL",&dyECAL);
    Tout->Branch("dxECAL4vect",&dxECAL4vect);
    Tout->Branch("dyECAL4vect",&dyECAL4vect);
    Tout->Branch("dtADC",&dtADC);
    
    double pthtar_e,pphtar_e; //SBS thetatarget,phitarget from electron scattering angles

    Tout->Branch("pthtar_e",&pthtar_e);
    Tout->Branch("pphtar_e",&pphtar_e);
    
    //SBS reconstructed track info:
    double thtar,phtar,ytar,vz;

    Tout->Branch("thtar",&thtar);
    Tout->Branch("phtar",&phtar);
    Tout->Branch("ytar",&ytar);
    Tout->Branch("vz",&vz);
    
    double xfp,yfp,thfp,phfp;

    Tout->Branch("xfp",&xfp);
    Tout->Branch("yfp",&yfp);
    Tout->Branch("thfp",&thfp);
    Tout->Branch("phfp",&phfp);

    double rxfp,ryfp,rthfp,rphfp;
    Tout->Branch("rxfp",&rxfp);
    Tout->Branch("ryfp",&ryfp);
    Tout->Branch("rthfp",&rthfp);
    Tout->Branch("rphfp",&rphfp);

    double dxfp,dyfp,dthfp,dphfp;
    Tout->Branch("dxfp",&dxfp);
    Tout->Branch("dyfp",&dyfp);
    Tout->Branch("dthfp",&dthfp);
    Tout->Branch("dphfp",&dphfp);
    
    int ntrackFPP,ntrackFT;
    double thetaFPP,phiFPP,scloseFPP,zcloseFPP;
    double xFPP,yFPP,xpFPP,ypFPP;

    Tout->Branch("ntrackFT",&ntrackFT);
    Tout->Branch("ntrackFPP",&ntrackFPP);
    Tout->Branch("thetaFPP",&thetaFPP);
    Tout->Branch("phiFPP",&phiFPP);
    Tout->Branch("scloseFPP",&scloseFPP);
    Tout->Branch("zcloseFPP",&zcloseFPP);
    Tout->Branch("xFPP",&xFPP);
    Tout->Branch("yFPP",&yFPP);
    Tout->Branch("xpFPP",&xpFPP);
    Tout->Branch("ypFPP",&ypFPP);
    
    double xHCAL,yHCAL,EHCAL,tHCAL_ADC;
    double xECAL,yECAL,EECAL,tECAL_ADC;
    int nblkECAL,nblkHCAL;

    Tout->Branch("xHCAL",&xHCAL);
    Tout->Branch("yHCAL",&yHCAL);
    Tout->Branch("EHCAL",&EHCAL);
    Tout->Branch("tHCAL_ADC",&tHCAL_ADC);
    Tout->Branch("nblkHCAL",&nblkHCAL);
    Tout->Branch("xECAL",&xECAL);
    Tout->Branch("yECAL",&yECAL);
    Tout->Branch("tECAL_ADC",&tECAL_ADC);
    Tout->Branch("nblkECAL",&nblkECAL);

    int ECAL_MAXBLOCKS = 30;
    
    double ECAL_clusblk_e[ECAL_MAXBLOCKS];
    double ECAL_clusblk_atime[ECAL_MAXBLOCKS];
    double ECAL_clusblk_again[ECAL_MAXBLOCKS];
    double ECAL_clusblk_id[ECAL_MAXBLOCKS];
    double ECAL_clusblk_x[ECAL_MAXBLOCKS];
    double ECAL_clusblk_y[ECAL_MAXBLOCKS];

    Tout->Branch("ECAL_clusblk_e", ECAL_clusblk_e);
    Tout->Branch("ECAL_clusblk_atime", ECAL_clusblk_atime);
    Tout->Branch("ECAL_clusblk_again", ECAL_clusblk_again);
    Tout->Branch("ECAL_clusblk_id", ECAL_clusblk_id);
    Tout->Branch("ECAL_clusblk_x", ECAL_clusblk_x);
    Tout->Branch("ECAL_clusblk_y", ECAL_clusblk_y);
    
    int HCAL_MAXBLOCKS = 30;
    double HCAL_clusblk_e[HCAL_MAXBLOCKS];
    double HCAL_clusblk_atime[HCAL_MAXBLOCKS];
    double HCAL_clusblk_again[HCAL_MAXBLOCKS];
    double HCAL_clusblk_id[HCAL_MAXBLOCKS];
    double HCAL_clusblk_x[HCAL_MAXBLOCKS];
    double HCAL_clusblk_y[HCAL_MAXBLOCKS];

    Tout->Branch("HCAL_clusblk_e", HCAL_clusblk_e);
    Tout->Branch("HCAL_clusblk_atime", HCAL_clusblk_atime);
    Tout->Branch("HCAL_clusblk_again", HCAL_clusblk_again);
    Tout->Branch("HCAL_clusblk_id", HCAL_clusblk_id);
    Tout->Branch("HCAL_clusblk_x", HCAL_clusblk_x);
    Tout->Branch("HCAL_clusblk_y", HCAL_clusblk_y);

    
    
    double helicity;
    double IHWP;

    Tout->Branch("helicity",&helicity);
    Tout->Branch("IHWP",&IHWP);
    
    //branches: do all H(e,e'p) branches and a few others:
    
    gep_tree_data *T = new gep_tree_data(C);

    C->SetBranchStatus("*",0);
    C->SetBranchStatus("sbs.tr.*",1);
    C->SetBranchStatus("heep.*",1);
    C->SetBranchStatus("sbs.gemFT.track.*",1);
    C->SetBranchStatus("sbs.gemFPP.track.*",1);
    C->SetBranchStatus("sbs.hcal.x",1);
    C->SetBranchStatus("sbs.hcal.y",1);
    C->SetBranchStatus("sbs.hcal.e",1);
    C->SetBranchStatus("sbs.hcal.atimeblk",1);
    C->SetBranchStatus("sbs.hcal.nblk",1);
    C->SetBranchStatus("earm.ecal.e",1);
    C->SetBranchStatus("earm.ecal.x",1);
    C->SetBranchStatus("earm.ecal.y",1);
    C->SetBranchStatus("earm.ecal.atimeblk",1);
    C->SetBranchStatus("earm.ecal.nblk",1);
    C->SetBranchStatus("sbs.x_fcp",1);
    C->SetBranchStatus("sbs.y_fcp",1);
    C->SetBranchStatus("sbs.y_fcp",1);
    C->SetBranchStatus("sbs.x_bcp",1);
    C->SetBranchStatus("sbs.y_bcp",1);
    C->SetBranchStatus("sbs.y_bcp",1);
    C->SetBranchStatus("sbs.ecal.clus_blk.*",1);
    C->SetBranchStatus("sbs.hcal.clus_blk.*",1);
    C->SetBranchStatus("g.*",1);
    C->SetBranchStatus("scalhel.*",1);
    C->SetBranchStatus("IGL1I00OD16_16",1);
    //    C->SetBranchStatus("sbs.x_bcp",1);
    
    TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);

    int treenum = -1, oldtreenum = -1;

    long nevent = 0;
    while( C->GetEntry(nevent++) ){
      if( nevent % 10000 == 0 ) cout << "event " << nevent << ", run number = " << int(T->g_runnum) << endl;

      treenum = C->GetTreeNumber();
      if( treenum != oldtreenum ){
	std::cout << "Switching trees, run number, file num, file name, nevent = " << runnum << ", " << treenum << ", " 
		  << C->GetFile()->GetName() << ", " << nevent << std::endl;
	GlobalCut->UpdateFormulaLeaves();
	oldtreenum = treenum;
      }

      bool passedcut = GlobalCut->EvalInstance(0) != 0;
      if( passedcut ){

	runnum = int(T->g_runnum);
	etheta = T->heep_etheta;
	ephi = T->heep_ephi;
	ptheta = T->heep_ptheta;
	pphi = T->heep_pphi;
	ecalo = T->heep_ecalo;
	pp = T->heep_pp;
	eprime_eth = T->heep_eprime_eth;
	pp_eth = T->heep_pp_eth;
	pth_eth = T->heep_pth_eth;
	pp_pth = T->heep_pp_pth;
	eprime_pth = T->heep_eprime_pth;
	eth_pth = T->heep_eth_pth;
	eprime_pp = T->heep_eprime_pp;
	eth_pp = T->heep_eth_pp;
	Q2_pp = T->heep_Q2_pp;
	Q2_eth = T->heep_Q2_eth;
	Q2_pth = T->heep_Q2_pth;
	Q2_p4vect = T->heep_Q2_p4vect;
	Q2_e4vect = T->heep_Q2_e4vect;
	eps_pp = T->heep_eps_pp;
	eps_eth = T->heep_eps_eth;
	eps_pth = T->heep_eps_pth;
	eps_4vect = T->heep_eps_4vect;
	K_eth = T->heep_K_eth;
	K_pth = T->heep_K_pth;
	K_pp = T->heep_K_pp;
	K_p4vect = T->heep_K_p4vect;
	dpp = T->heep_dpp;
	dpe = T->heep_dpe;
	dphi = T->heep_dphi;
	acopl = T->heep_acoplanarity;
	inel = T->heep_inelasticity;
	dxECAL = T->heep_dxECAL;
	dyECAL = T->heep_dyECAL;
	dxECAL4vect = T->heep_dxECAL4vect;
	dyECAL4vect = T->heep_dyECAL4vect;
	dtADC = T->heep_dt_ADC;

	TVector3 pnhat_e_global(sin(pth_eth)*cos(ephi+PI),sin(pth_eth)*sin(ephi+PI),cos(pth_eth));

	pthtar_e = pnhat_e_global.Dot(SBSxaxis)/pnhat_e_global.Dot(SBSzaxis);
	pphtar_e = pnhat_e_global.Dot(SBSyaxis)/pnhat_e_global.Dot(SBSzaxis);

	thtar = T->sbs_tr_tg_th[0];
	phtar = T->sbs_tr_tg_ph[0];
	ytar = T->sbs_tr_tg_y[0];
	vz = T->sbs_tr_vz[0];

	xfp = T->sbs_tr_x[0];
	yfp = T->sbs_tr_y[0];
	thfp = T->sbs_tr_th[0];
	phfp = T->sbs_tr_ph[0];

	rxfp = T->sbs_tr_r_x[0];
	ryfp = T->sbs_tr_r_y[0];
	rthfp = T->sbs_tr_r_th[0];
	rphfp = T->sbs_tr_r_ph[0];

	dxfp = T->sbs_tr_d_x[0];
	dyfp = T->sbs_tr_d_y[0];
	dthfp = T->sbs_tr_d_th[0];
	dphfp = T->sbs_tr_d_ph[0];
	
	ntrackFT = int( T->sbs_tr_n );
	ntrackFPP = int( T->sbs_gemFPP_track_ntrack );

	thetaFPP = -10000;
	phiFPP = -10000.;
	scloseFPP = -10000.;
	zcloseFPP = -10000.;

	xFPP = -10000.;
	yFPP = xpFPP = ypFPP = -10000;

	if( ntrackFPP > 0 ){
	  thetaFPP = T->sbs_gemFPP_track_theta[0];
	  phiFPP = T->sbs_gemFPP_track_phi[0];
	  scloseFPP = T->sbs_gemFPP_track_sclose[0];
	  zcloseFPP = T->sbs_gemFPP_track_zclose[0];
	  xFPP = T->sbs_gemFPP_track_x[0];
	  yFPP = T->sbs_gemFPP_track_y[0];
	  xpFPP = T->sbs_gemFPP_track_xp[0];
	  ypFPP = T->sbs_gemFPP_track_yp[0];
	}

	xHCAL = T->sbs_hcal_x;
	yHCAL = T->sbs_hcal_y;
	EHCAL = T->sbs_hcal_e;
	tHCAL_ADC = T->sbs_hcal_atimeblk;
	nblkHCAL = int(T->sbs_hcal_nblk);

	xECAL = T->earm_ecal_x;
	yECAL = T->earm_ecal_y;
	EECAL = T->earm_ecal_e;
	tECAL_ADC = T->earm_ecal_atimeblk;

	nblkECAL = int(T->earm_ecal_nblk);

	helicity = T->scalhel_true_hel;
	IHWP = T->IGL1I00OD16_16; //half-wave-plate state
	
	//eventually we want to add the ECAL and HCAL cluster individual block stuff:
	
	
	Tout->Fill();
      }
      
    }

    fout->Write();
  }

 
  
}
