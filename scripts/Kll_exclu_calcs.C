#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>

#include "genrp_tree.C"
#include "TVector3.h"
#include "TF1.h"
#include "TLorentzVector.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <iomanip>
#include <cmath>
#include <chrono>
#include "TTreeFormula.h"
#include "TCut.h"

using namespace std;

// Running the script:
//
// Arguments
// 1. The name of the input rootfile(s), with the file extension .root
// 2. Name you want for the output rootfile, with the .root file extension

double TriangleFunc(double a, double b, double c){
  return pow(a,2)+pow(b,2)+pow(c,2) - 2.*a*b - 2.*a*c - 2.*b*c;
}

void Kll_exclu_calcs( const char *rootfile, const char *outputFile, const char *cutstring="" ){

  TCut globalcut = cutstring;

  
  
  // Timer
  auto start = chrono::high_resolution_clock::now(); 
  
  // Double check that outputFileName argument is written in .root format
  string outFileName = outputFile;
  if ( outFileName.find(".root") == string::npos) {
    cerr << "Warning: Output file name does not have a .root extension. Adding it automatically..." << endl;
    outFileName += ".root";
  }

  // Set-up a TChain to handle multiple input root files
  TChain *tree = new TChain("T"); 
  tree->Add(rootfile);

  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut",globalcut,tree);
  
  genrp_tree *T = new genrp_tree(tree); 

  // Check tree entries
  // if ( tree->GetEntries() == 0 ) {
  //   cerr << "ERROR: No valid entries found in input root files!" << endl;
  //   return;
  // }

  // Speed up processing, disables all unnecessary branches
  tree->SetBranchStatus("*", 0);

  // Enable necessary branches
  // tree->SetBranchStatus("bb.tr.n", 1);
  // tree->SetBranchStatus("sbs.tr.n", 1);
  // tree->SetBranchStatus("bb.tr.p", 1);
  // tree->SetBranchStatus("sbs.tr.p", 1);
  // tree->SetBranchStatus("bb.tr.px", 1);
  // tree->SetBranchStatus("bb.tr.py", 1);
  // tree->SetBranchStatus("bb.tr.pz", 1);
  tree->SetBranchStatus("bb.tr.*",1);
  tree->SetBranchStatus("sbs.tr.*",1);
  tree->SetBranchStatus("scalhel.hel",1);
  tree->SetBranchStatus("sbs.gemCeR.track.*",1);
  tree->SetBranchStatus("sbs.hcal.e",1);
  tree->SetBranchStatus("sbs.hcal.atimeblk",1);
  tree->SetBranchStatus("sbs.hcal.x",1);
  tree->SetBranchStatus("sbs.hcal.y",1);

  tree->SetBranchStatus("bb.ps.e",1);
  tree->SetBranchStatus("bb.ps.atimeblk",1);
  tree->SetBranchStatus("bb.sh.atimeblk",1);
  tree->SetBranchStatus("bb.sh.e",1);
  
  // tree->SetBranchStatus("bb.x_fcp",1);
  // tree->SetBranchStatus("bb.y_fcp",1);
  // tree->SetBranchStatus("bb.z_fcp",1);

  // tree->SetBranchStatus("bb.x_bcp",1);
  // tree->SetBranchStatus("bb.y_bcp",1);
  // tree->SetBranchStatus("bb.z_bcp",1);

  // tree->SetBranchStatus("sbs.x_fcp",1);
  // tree->SetBranchStatus("sbs.y_fcp",1);
  // tree->SetBranchStatus("sbs.z_fcp",1);

  // tree->SetBranchStatus("sbs.x_bcp",1);
  // tree->SetBranchStatus("sbs.y_bcp",1);
  // tree->SetBranchStatus("sbs.z_bcp",1);
  
  // tree->SetBranchStatus("sbs.tr.px", 1);
  // tree->SetBranchStatus("sbs.tr.py", 1);
  // tree->SetBranchStatus("sbs.tr.pz", 1);

  // Define variables
  Double_t Epi, Ep, Ppi, Pp, ThetaPi, ThetaP, PhiPi, PhiP;
  const Double_t mn = 0.939565; // neutron mass (GeV)
  const Double_t mp = 0.938272; // protron mass (GeV)
  const Double_t mpi = 0.13957; // pion mass (GeVC)

  // Define varibales to be computed
  Double_t Eg_pi, Eg_p, Eg_comb, s_pi_p;
  Double_t Emiss, Pmiss_par, Pmiss_perp;
  Double_t Ptrans_sq, Mmiss_sq;
  Double_t Emiss_pi;
  
  // Three-vector components for Ppion and Pproton
  Double_t Ppion_x, Ppion_y, Ppion_z;
  Double_t Pproton_x, Pproton_y, Pproton_z;
  
  // Four-vector components for Ppion4vect, Pproton4vect and total 4-vector
  Double_t Ppion4vect_px, Ppion4vect_py, Ppion4vect_pz, Ppion4vect_E;
  Double_t Pproton4vect_px, Pproton4vect_py, Pproton4vect_pz, Pproton4vect_E;
  Double_t Ptotal4vect_px, Ptotal4vect_py, Ptotal4vect_pz, Ptotal4vect_E;
  Double_t acoplanarity;

  // Create an output root file
  TFile *outFile = new TFile(outputFile, "RECREATE");
  if ( !outFile || outFile->IsZombie() ) {
    cerr << "ERROR: could not create output root file!!" << endl;
    return;
  }

  outFile->cd();  // IDK if this is needed 

  // enables statistical error handling for histograms
  TH1D::SetDefaultSumw2();
  
  // Create new tress to store computed results
  TTree *outTree = new TTree("KLL", "Tree with computed variables");

  // Store computed variables in TTree
  outTree->Branch("Eg_pi", &Eg_pi, "Eg_pi/D");     // Eqn (10)
  outTree->Branch("Eg_p", &Eg_p, "Eg_p/D");        // Eqn (11)
  outTree->Branch("Eg_comb", &Eg_comb, "Eg_comb/D");    // Eqn (12)
  outTree->Branch("s_pi_p", &s_pi_p, "s_pi_p/D");     // Eqn (13)
  outTree->Branch("Emiss", &Emiss, "Emiss/D");  // Missing Energy
  outTree->Branch("Emiss_pi",&Emiss_pi, "Emiss_pi/D");
  outTree->Branch("Pmiss_par", &Pmiss_par, "Pmiss_par/D");   // Parallel missing momentum 
  outTree->Branch("Pmiss_perp", &Pmiss_perp, "Pmiss_perp/D");  // Perpendicular missing momentum 
  outTree->Branch("Ptrans_sq", &Ptrans_sq, "Ptrans_sq/D");  // Squared transverse momentum 
  outTree->Branch("Mmiss_sq", &Mmiss_sq, "Mmiss_sq/D");     // Missing mass squared
  outTree->Branch("ThetaPi", &ThetaPi, "ThetaPi/D");
  outTree->Branch("ThetaP", &ThetaP, "ThetaP/D");
  outTree->Branch("PhiPi", &PhiPi, "PhiPi/D");
  outTree->Branch("PhiP", &PhiP, "PhiP/D");

  outTree->Branch("acoplanarity", &acoplanarity, "acoplanarity/D");
  
  // Energies
  outTree->Branch("Epi", &Epi, "Epi/D");
  outTree->Branch("Ep", &Ep, "Ep/D");

  outTree->Branch("Ppi", &Ppi, "Ppi/D");
  outTree->Branch("Pp", &Pp, "Pp/D");
  
  // 3-vector components for Ppion
  outTree->Branch("Ppion_x", &Ppion_x, "Ppion_x/D");
  outTree->Branch("Ppion_y", &Ppion_y, "Ppion_y/D");
  outTree->Branch("Ppion_z", &Ppion_z, "Ppion_z/D");
  
  // 3-vector components for Pproton
  outTree->Branch("Pproton_x", &Pproton_x, "Pproton_x/D");
  outTree->Branch("Pproton_y", &Pproton_y, "Pproton_y/D");
  outTree->Branch("Pproton_z", &Pproton_z, "Pproton_z/D");
  
  // Four-vector components for Ppion4vect
  outTree->Branch("Ppion4vect_px", &Ppion4vect_px, "Ppion4vect_px/D");
  outTree->Branch("Ppion4vect_py", &Ppion4vect_py, "Ppion4vect_py/D");
  outTree->Branch("Ppion4vect_pz", &Ppion4vect_pz, "Ppion4vect_pz/D");
  outTree->Branch("Ppion4vect_E", &Ppion4vect_E, "Ppion4vect_E/D");
  
  // Four-vector components for Pproton4vect
  outTree->Branch("Pproton4vect_px", &Pproton4vect_px, "Pproton4vect_px/D");
  outTree->Branch("Pproton4vect_py", &Pproton4vect_py, "Pproton4vect_py/D");
  outTree->Branch("Pproton4vect_pz", &Pproton4vect_pz, "Pproton4vect_pz/D");
  outTree->Branch("Pproton4vect_E", &Pproton4vect_E, "Pproton4vect_E/D");
  
  // Four-vector components for total 4-vector
  outTree->Branch("Ptotal4vect_px", &Ptotal4vect_px, "Ptotal4vect_px/D");
  outTree->Branch("Ptotal4vect_py", &Ptotal4vect_py, "Ptotal4vect_py/D");
  outTree->Branch("Ptotal4vect_pz", &Ptotal4vect_pz, "Ptotal4vect_pz/D");
  outTree->Branch("Ptotal4vect_E", &Ptotal4vect_E, "Ptotal4vect_E/D");

  double thetaFPP, phiFPP, scloseFPP, zcloseFPP;
  double vzBB, vzSBS;
  double HCALatime, BBSHatime, BBPSatime;
  double EHCAL, xHCAL, yHCAL;
  double ESH, EPS;
  
  outTree->Branch("vzBB", &vzBB);
  outTree->Branch("vzSBS", &vzSBS);
  outTree->Branch("THCAL",&HCALatime );
  outTree->Branch("EHCAL",&EHCAL );
  outTree->Branch("xHCAL",&xHCAL );
  outTree->Branch("yHCAL",&yHCAL );
  outTree->Branch("ESH",&ESH);
  outTree->Branch("EPS",&EPS);
  outTree->Branch("TSH", &BBSHatime);
  outTree->Branch("TPS", &BBPSatime);
  

  outTree->Branch( "thetaFPP", &thetaFPP );
  outTree->Branch( "phiFPP", &phiFPP );
  outTree->Branch( "scloseFPP", &scloseFPP );
  outTree->Branch( "zcloseFPP", &zcloseFPP );

  
  
  double helicity;
  outTree->Branch( "helicity", &helicity );
  outTree->Branch( "thetaFPP", &thetaFPP );
  outTree->Branch( "phiFPP", &phiFPP );
  outTree->Branch( "scloseFPP", &scloseFPP );
  outTree->Branch( "zcloseFPP", &zcloseFPP );

  double s_Mandelstam, t_Mandelstam, u_Mandelstam;
  double costhCM;
  outTree->Branch( "s", &s_Mandelstam );
  outTree->Branch( "t", &t_Mandelstam );
  outTree->Branch( "u", &u_Mandelstam );
  outTree->Branch( "costhCM", &costhCM );
  
  // Create histograms
  // Mandelstam variable histograms
  TH1D *hs = new TH1D("hs", "s (GeV^2)", 150, 6.0, 18.0);
  TH1D *ht = new TH1D("ht", "-t (GeV^2)", 150, 3.0, 18.0);
  TH1D *hu = new TH1D("hu", "-u (GeV^2)", 150, 6.0, 18.0);
  
  TH1D *hEg_pi = new TH1D("hEg_pi", "Reconstructed Photon Energy (Pion); E_{#gamma} (GeV);rate/bin (Hz)", 100, 0, 4);
  TH1D *hEg_p = new TH1D("hEg_p", "Reconstructed Photon Energy (Photon); E_{#gamma} (GeV);rate/bin (Hz)", 100, 0, 4);
  TH1D *hEg_comb= new TH1D("hEg_comb", "Combined Photon Energy; E_{#gamma} (GeV);rate/bin (Hz)", 100, 0, 4);
  TH1D *hEmiss = new TH1D("hEmiss", "Missing Energy;E_{miss} (GeV);Rate (Hz/bin)", 100, -1, 1);
  TH1D *hPmiss_par = new TH1D("hPmiss_par", "Parallel Missing Momentum;P_{miss,||} (GeV);Rate (Hz/bin)", 100, -0.5, 0.5);
  TH1D *hPmiss_perp = new TH1D("hPmiss_perp", "Perpendicular Missing Momentum;P_{miss,⊥} (GeV);Rate (Hz/bin)", 100, 0, 1);
  TH1D *hPtrans_sq = new TH1D("hPtrans_sq", "Squared Transverse Momentum;P_{T}^2 (GeV^2);Rate (Hz/bin)", 100, 0, 0.1);
  TH1D *hMmiss_sq = new TH1D("hMmiss_sq", "Missing Mass Squared;M_{miss}^2 (GeV^2);Rate (Hz/bin)", 100, 0, 0.5);

  // Energies
  TH1D *hEpi = new TH1D("hEpi", "Pion Energy; E_{pi} (GeV);Counts", 100, 0, 3);
  TH1D *hEp  = new TH1D("hEp", "Proton Energy; E_{p} (GeV);Counts", 100, 0, 3);

  // 3-vector components for Ppion
  TH1D *hPpion_x = new TH1D("hPpion_x", "Pion momentum X; P_{pi,x} (GeV/c);Counts", 100, -3, 3);
  TH1D *hPpion_y = new TH1D("hPpion_y", "Pion momentum Y; P_{pi,y} (GeV/c);Counts", 100, -3, 3);
  TH1D *hPpion_z = new TH1D("hPpion_z", "Pion momentum Z; P_{pi,z} (GeV/c);Counts", 100, -3, 3);

  // 3-vector components for Pproton
  TH1D *hPproton_x = new TH1D("hPproton_x", "Proton momentum X; P_{p,x} (GeV/c);Counts", 100, -3, 3);
  TH1D *hPproton_y = new TH1D("hPproton_y", "Proton momentum Y; P_{p,y} (GeV/c);Counts", 100, -3, 3);
  TH1D *hPproton_z = new TH1D("hPproton_z", "Proton momentum Z; P_{p,z} (GeV/c);Counts", 100, -3, 3);

  // Four-vector components for Ppion4vect
  TH1D *hPpion4vect_px = new TH1D("hPpion4vect_px", "Pion 4-vector px; P_{pi,4,x} (GeV/c);Counts", 100, -3, 3);
  TH1D *hPpion4vect_py = new TH1D("hPpion4vect_py", "Pion 4-vector py; P_{pi,4,y} (GeV/c);Counts", 100, -3, 3);
  TH1D *hPpion4vect_pz = new TH1D("hPpion4vect_pz", "Pion 4-vector pz; P_{pi,4,z} (GeV/c);Counts", 100, -3, 3);
  TH1D *hPpion4vect_E  = new TH1D("hPpion4vect_E", "Pion 4-vector Energy; E_{pi,4} (GeV);Counts", 100, 0, 3);

  // Four-vector components for Pproton4vect
  TH1D *hPproton4vect_px = new TH1D("hPproton4vect_px", "Proton 4-vector px; P_{p,4,x} (GeV/c);Counts", 100, -3, 3);
  TH1D *hPproton4vect_py = new TH1D("hPproton4vect_py", "Proton 4-vector py; P_{p,4,y} (GeV/c);Counts", 100, -3, 3);
  TH1D *hPproton4vect_pz = new TH1D("hPproton4vect_pz", "Proton 4-vector pz; P_{p,4,z} (GeV/c);Counts", 100, -3, 3);
  TH1D *hPproton4vect_E  = new TH1D("hPproton4vect_E", "Proton 4-vector Energy; E_{p,4} (GeV);Counts", 100, 0, 3);

  // Four-vector components for total 4-vector
  TH1D *hPtotal4vect_px = new TH1D("hPtotal4vect_px", "Total 4-vector px; P_{total,4,x} (GeV/c);Counts", 100, -3, 3);
  TH1D *hPtotal4vect_py = new TH1D("hPtotal4vect_py", "Total 4-vector py; P_{total,4,y} (GeV/c);Counts", 100, -3, 3);
  TH1D *hPtotal4vect_pz = new TH1D("hPtotal4vect_pz", "Total 4-vector pz; P_{total,4,z} (GeV/c);Counts", 100, 0, 4);
  TH1D *hPtotal4vect_E  = new TH1D("hPtotal4vect_E", "Total 4-vector Energy; E_{total,4} (GeV);Counts", 100, 1, 5);

  // Histograms for angular variables
  TH1D *hThetaPi = new TH1D("hThetaPi", "Theta (Pion); #theta_{pi} (rad);Counts", 100, 0, TMath::Pi());
  TH1D *hThetaP  = new TH1D("hThetaP", "Theta (Proton); #theta_{p} (rad);Counts", 100, 0, TMath::Pi());
  TH1D *hPhiPi   = new TH1D("hPhiPi", "Phi (Pion); #phi_{pi} (rad);Counts", 100, -TMath::Pi(), TMath::Pi());
  TH1D *hPhiP    = new TH1D("hPhiP", "Phi (Proton); #phi_{p} (rad);Counts", 100, -TMath::Pi(), TMath::Pi());

  

  // Process entries
  //  Long64_t nEntries = tree->GetEntries();
  //  cout << "Processing " << nEntries << " entries..." << endl;

  // loop through entries
  //for ( Long64_t i = 0; i < nEntries; i++) {

  int treenum=-1, oldtreenum=-1;

  long i=0;
  
  while( tree->GetEntry(i++) ){
    
    if( i%10000 == 0 ) cout << i << endl;
    //if (i % (nEntries / 10) == 0) cout << "Processed " << (i * 100.0 / nEntries) << "% entries." << endl;
    
    treenum = tree->GetTreeNumber();
    if( treenum != oldtreenum ){
      oldtreenum = treenum;
      GlobalCut->UpdateFormulaLeaves();
    }

    bool passed_globalcut = GlobalCut->EvalInstance(0) != 0;
    
    
    // Apply cut
    if ( T->bb_tr_n >= 1. && T->sbs_tr_n >= 1. && passed_globalcut && T->sbs_gemCeR_track_ntrack>=1. ){

      // Energy for pion and Proton
      Epi = sqrt(pow(T->bb_tr_p[0],2) + pow(mpi, 2)); 
      Ep = sqrt(pow(T->sbs_tr_p[0],2) + pow(mp, 2)); 

      // Constructing 3D and 4-vector momentum 
      TVector3 Ppion( T->bb_tr_px[0], T->bb_tr_py[0], T->bb_tr_pz[0] ); 
      TVector3 Pproton( T->sbs_tr_px[0], T->sbs_tr_py[0], T->sbs_tr_pz[0] ); 
      TLorentzVector Ppion4vect( Ppion, Epi ); 
      TLorentzVector Pproton4vect( Pproton, Ep );
      TLorentzVector Ptotal4vect = Ppion4vect+Pproton4vect;
      
      Ppi = T->bb_tr_p[0]; 
      Pp = T->sbs_tr_p[0]; 
      ThetaPi = acos(T->bb_tr_pz[0]/T->bb_tr_p[0]);
      ThetaP = acos(T->sbs_tr_pz[0]/T->sbs_tr_p[0]);
      PhiPi = atan2(T->bb_tr_py[0], T->bb_tr_px[0]);
      PhiP = atan2(T->sbs_tr_py[0], T->sbs_tr_px[0]);

      // Compute the photon energies (Eqn 10, 11, 12, 13)
      Eg_pi = ( 2.*mn*Epi + mp*mp - mpi*mpi - mn*mn ) / ( 2.*(mn + Ppi*TMath::Cos(ThetaPi) - Epi));
      Eg_p = ( 2.*mn*Ep + mpi*mpi - mp*mp - mn*mn ) / ( 2.*(mn + Pp*TMath::Cos(ThetaP) - Ep));

      //     double s_ppi = pow(Ep+Epi,2) - (Ppion+Pproton).Mag2()
      double s_ppi = Ptotal4vect.Mag2(); 
      Eg_comb = ( s_ppi - mn*mn ) / ( 2.*(Epi + Ep - Ppi*TMath::Cos(ThetaPi) - Pp*TMath::Cos(ThetaP)));

      TLorentzVector Pgamma4vect(0,0,Eg_comb, Eg_comb);
      TLorentzVector Ptarg4vect(0,0,0,mn);

      s_Mandelstam = (Pgamma4vect+Ptarg4vect).M2();

      // t = (p1 - p3)^2 
      t_Mandelstam = (Pgamma4vect - Ppion4vect).M2();
      // s + t + u = mn^2 + mp^2 + mpi^2
      u_Mandelstam = pow(mn,2)+pow(mp,2)+pow(mpi,2)-s_Mandelstam - t_Mandelstam;

      costhCM = (s_Mandelstam * ( t_Mandelstam - u_Mandelstam ) - pow(mn,2)*(pow(mpi,2)-pow(mp,2))) / (sqrt(TriangleFunc( s_Mandelstam,0.0,pow(mn,2)) ) * sqrt( TriangleFunc( s_Mandelstam, pow(mpi,2), pow(mp,2))));

      
      s_pi_p = s_ppi;
      
      // Compute missing kinematics
      Emiss = Eg_p + mn - Ep - Epi;

      Emiss_pi = Eg_pi + mn - Ep - Epi;
      //TODO: fix formulas below using TVector3, TLorentzVector, etc
      // Use "WAPP_FOM_quick_and_dirty.C" from g4sbs as example. 
      // Compute missing momentum
      TLorentzVector Pmiss4vect = TLorentzVector(0, 0, Eg_comb, Eg_comb) + TLorentzVector(0, 0, 0, mn) - Ppion4vect - Pproton4vect;
      TVector3 Pmiss = Pmiss4vect.Vect();  // Extract three-momentum part

      // Compute four-momentum transfer Q4vect
      TLorentzVector Q4vect = TLorentzVector(0, 0, Eg_comb, Eg_comb) - Ppion4vect;
      TVector3 q3vect = Q4vect.Vect();  // Extract three-momentum component

      //Emiss = Pmiss4vect.E();

      TVector3 nhat_pi = Pgamma4vect.Vect().Cross( Ppion4vect.Vect() ).Unit();
      TVector3 nhat_proton = Pgamma4vect.Vect().Cross( Pproton4vect.Vect() ).Unit();
      
      acoplanarity = acos(-nhat_pi.Dot(nhat_proton));
      
      // Compute Pmiss components
      Pmiss_par = Pmiss.Dot(q3vect.Unit());  // Projection onto q direction
      Pmiss_perp = (Pmiss - q3vect.Unit() * Pmiss.Dot(q3vect.Unit())).Mag();  // Perpendicular component
      //Ptrans_sq = (Pmiss - Pmiss_par * q3vect.Unit()).Mag2();
      //Ptrans_sq = (Ppi + Pp)*(Ppi + Pp) - (Epi + Ep)*(Epi + Ep);
      Ptrans_sq = Ptotal4vect.Perp2();

      Mmiss_sq = Pmiss4vect.M2();
      //      Mmiss_sq = pow(Eg_comb + mn - Ep - Epi, 2);

      // --- Fill new variables ---
      // 3-vector components
      Ppion_x = T->bb_tr_px[0];
      Ppion_y = T->bb_tr_py[0];
      Ppion_z = T->bb_tr_pz[0];
      
      Pproton_x = T->sbs_tr_px[0];
      Pproton_y = T->sbs_tr_py[0];
      Pproton_z = T->sbs_tr_pz[0];
      
      // 4-vector components for pion
      Ppion4vect_px = Ppion4vect.Px();
      Ppion4vect_py = Ppion4vect.Py();
      Ppion4vect_pz = Ppion4vect.Pz();
      Ppion4vect_E  = Ppion4vect.E();
      
      // 4-vector components for proton
      Pproton4vect_px = Pproton4vect.Px();
      Pproton4vect_py = Pproton4vect.Py();
      Pproton4vect_pz = Pproton4vect.Pz();
      Pproton4vect_E  = Pproton4vect.E();
      
      // 4-vector components for total
      Ptotal4vect_px = Ptotal4vect.Px();
      Ptotal4vect_py = Ptotal4vect.Py();
      Ptotal4vect_pz = Ptotal4vect.Pz();
      Ptotal4vect_E  = Ptotal4vect.E();
      
      // Fill Histograms
      hEg_pi->Fill(Eg_pi);
      hEg_p->Fill(Eg_p);
      hEg_comb->Fill(Eg_comb);
      hEmiss->Fill(Emiss);
      hPmiss_par->Fill(Pmiss_par);
      hPmiss_perp->Fill(Pmiss_perp);
      hPtrans_sq->Fill(Ptrans_sq);
      hMmiss_sq->Fill(Mmiss_sq);

      // New histograms filling
      hEpi->Fill(Epi);
      hEp->Fill(Ep);
      
      hPpion_x->Fill(Ppion_x);
      hPpion_y->Fill(Ppion_y);
      hPpion_z->Fill(Ppion_z);
      
      hPproton_x->Fill(Pproton_x);
      hPproton_y->Fill(Pproton_y);
      hPproton_z->Fill(Pproton_z);
      
      hPpion4vect_px->Fill(Ppion4vect_px);
      hPpion4vect_py->Fill(Ppion4vect_py);
      hPpion4vect_pz->Fill(Ppion4vect_pz);
      hPpion4vect_E->Fill(Ppion4vect_E);
      
      hPproton4vect_px->Fill(Pproton4vect_px);
      hPproton4vect_py->Fill(Pproton4vect_py);
      hPproton4vect_pz->Fill(Pproton4vect_pz);
      hPproton4vect_E->Fill(Pproton4vect_E);
      
      hPtotal4vect_px->Fill(Ptotal4vect_px);
      hPtotal4vect_py->Fill(Ptotal4vect_py);
      hPtotal4vect_pz->Fill(Ptotal4vect_pz);
      hPtotal4vect_E->Fill(Ptotal4vect_E);
      
      hThetaPi->Fill(ThetaPi);
      hThetaP->Fill(ThetaP);
      hPhiPi->Fill(PhiPi);
      hPhiP->Fill(PhiP);

      vzBB = T->bb_tr_vz[0];
      vzSBS = T->sbs_tr_vz[0];
      HCALatime = T->sbs_hcal_atimeblk;
      BBSHatime = T->bb_sh_atimeblk;
      BBPSatime = T->bb_ps_atimeblk;
      EHCAL = T->sbs_hcal_e;
      xHCAL = T->sbs_hcal_x;
      yHCAL = T->sbs_hcal_y;
      ESH = T->bb_sh_e;
      EPS = T->bb_ps_e;

      helicity = T->scalhel_hel;
      thetaFPP = T->sbs_gemCeR_track_theta[0];
      phiFPP = T->sbs_gemCeR_track_phi[0];
      scloseFPP = T->sbs_gemCeR_track_sclose[0];
      zcloseFPP = T->sbs_gemCeR_track_zclose[0];

      // Store event in TTree
      outTree->Fill(); 
    }
  }


  // Write tree and histograms to the output file
  outTree->Write();
  hs->Write();
  ht->Write();
  hu->Write();
  hEg_pi->Write();
  hEg_p->Write();
  hEg_comb->Write();
  hEmiss->Write();
  hPmiss_par->Write();
  hPmiss_perp->Write();
  hPtrans_sq->Write();
  hMmiss_sq->Write();
  
  hEpi->Write();
  hEp->Write();
  
  hPpion_x->Write();
  hPpion_y->Write();
  hPpion_z->Write();
  
  hPproton_x->Write();
  hPproton_y->Write();
  hPproton_z->Write();
  
  hPpion4vect_px->Write();
  hPpion4vect_py->Write();
  hPpion4vect_pz->Write();
  hPpion4vect_E->Write();
  
  hPproton4vect_px->Write();
  hPproton4vect_py->Write();
  hPproton4vect_pz->Write();
  hPproton4vect_E->Write();
  
  hPtotal4vect_px->Write();
  hPtotal4vect_py->Write();
  hPtotal4vect_pz->Write();
  hPtotal4vect_E->Write();
  
  hThetaPi->Write();
  hThetaP->Write();
  hPhiPi->Write();
  hPhiP->Write();

  outFile->Close();

  // end timer
  auto end = chrono::high_resolution_clock::now();
  chrono::duration<double> elapsed = end - start;
  double elapsedTime = elapsed.count();
  double eventRate = i / elapsedTime; 
  
  cout << "=========================================" << endl;
  cout << "Processing complete! Results saved to: " << outputFile << endl;
  cout << "Total execution time: " << elapsedTime << " seconds" << endl;
  cout << "Processing speed: " << eventRate << " Hz" << endl; 

  
}
