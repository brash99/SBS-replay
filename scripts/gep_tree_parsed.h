//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Aug  8 08:32:20 2025 by ROOT version 6.30/04
// from TTree Tout/GEP elastic parsing
// found on file: /volatile/halla/sbs/puckett/gep3parsed_temp_aug6_2025.root
//////////////////////////////////////////////////////////

#ifndef gep_tree_parsed_h
#define gep_tree_parsed_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class gep_tree_parsed {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           runnum;
   Double_t        etheta;
   Double_t        ephi;
   Double_t        ptheta;
   Double_t        pphi;
   Double_t        ECALO;
   Double_t        pp;
   Double_t        eprime_eth;
   Double_t        pp_eth;
   Double_t        pp_pth;
   Double_t        eprime_pth;
   Double_t        eprime_pp;
   Double_t        eth_pp;
   Double_t        Q2_pp;
   Double_t        Q2_eth;
   Double_t        Q2_pth;
   Double_t        Q2_e4vect;
   Double_t        eps_pp;
   Double_t        eps_eth;
   Double_t        eps_pth;
   Double_t        eps_4vect;
   Double_t        K_eth;
   Double_t        K_pth;
   Double_t        K_pp;
   Double_t        K_p4vect;
   Double_t        dpp;
   Double_t        dpe;
   Double_t        dphi;
   Double_t        acopl;
   Double_t        inel;
   Double_t        dxECAL;
   Double_t        dyECAL;
   Double_t        dxECAL4vect;
   Double_t        dyECAL4vect;
   Double_t        dtADC;
   Double_t        pthtar_e;
   Double_t        pphtar_e;
   Double_t        thtar;
   Double_t        phtar;
   Double_t        ytar;
   Double_t        vz;
   Double_t        xfp;
   Double_t        yfp;
   Double_t        thfp;
   Double_t        phfp;
   Double_t        rxfp;
   Double_t        ryfp;
   Double_t        rthfp;
   Double_t        rphfp;
   Double_t        dxfp;
   Double_t        dyfp;
   Double_t        dthfp;
   Double_t        dphfp;
   Int_t           ntrackFT;
   Int_t           ntrackFPP;
   Double_t        thetaFPP;
   Double_t        phiFPP;
   Double_t        scloseFPP;
   Double_t        zcloseFPP;
   Double_t        xFPP;
   Double_t        yFPP;
   Double_t        xpFPP;
   Double_t        ypFPP;
   Int_t           nhitsFPP;
   Int_t           nhitsFT;
   Int_t           ngoodhitsFPP;
   Int_t           ngoodhitsFT;
   Double_t        chi2ndfFPP;
   Double_t        chi2ndfFT;
   Double_t        xHCAL;
   Double_t        yHCAL;
   Double_t        EHCAL;
   Double_t        tHCAL_ADC;
   Int_t           nblkHCAL;
   Double_t        xECAL;
   Double_t        yECAL;
   Double_t        tECAL_ADC;
   Int_t           nblkECAL;
   Double_t        ECAL_clusblk_e[33];   //[nblkECAL]
   Double_t        ECAL_clusblk_atime[33];   //[nblkECAL]
   Double_t        ECAL_clusblk_again[33];   //[nblkECAL]
   Int_t           ECAL_clusblk_id[33];   //[nblkECAL]
   Double_t        ECAL_clusblk_x[33];   //[nblkECAL]
   Double_t        ECAL_clusblk_y[33];   //[nblkECAL]
   Double_t        HCAL_clusblk_e[26];   //[nblkHCAL]
   Double_t        HCAL_clusblk_atime[26];   //[nblkHCAL]
   Double_t        HCAL_clusblk_again[26];   //[nblkHCAL]
   Int_t           HCAL_clusblk_id[26];   //[nblkHCAL]
   Double_t        HCAL_clusblk_x[26];   //[nblkHCAL]
   Double_t        HCAL_clusblk_y[26];   //[nblkHCAL]
   Double_t        helicity;
   Double_t        IHWP;

   // List of branches
   TBranch        *b_runnum;   //!
   TBranch        *b_etheta;   //!
   TBranch        *b_ephi;   //!
   TBranch        *b_ptheta;   //!
   TBranch        *b_pphi;   //!
   TBranch        *b_ECALO;   //!
   TBranch        *b_pp;   //!
   TBranch        *b_eprime_eth;   //!
   TBranch        *b_pp_eth;   //!
   TBranch        *b_pp_pth;   //!
   TBranch        *b_eprime_pth;   //!
   TBranch        *b_eprime_pp;   //!
   TBranch        *b_eth_pp;   //!
   TBranch        *b_Q2_pp;   //!
   TBranch        *b_Q2_eth;   //!
   TBranch        *b_Q2_pth;   //!
   TBranch        *b_Q2_e4vect;   //!
   TBranch        *b_eps_pp;   //!
   TBranch        *b_eps_eth;   //!
   TBranch        *b_eps_pth;   //!
   TBranch        *b_eps_4vect;   //!
   TBranch        *b_K_eth;   //!
   TBranch        *b_K_pth;   //!
   TBranch        *b_K_pp;   //!
   TBranch        *b_K_p4vect;   //!
   TBranch        *b_dpp;   //!
   TBranch        *b_dpe;   //!
   TBranch        *b_dphi;   //!
   TBranch        *b_acopl;   //!
   TBranch        *b_inel;   //!
   TBranch        *b_dxECAL;   //!
   TBranch        *b_dyECAL;   //!
   TBranch        *b_dxECAL4vect;   //!
   TBranch        *b_dyECAL4vect;   //!
   TBranch        *b_dtADC;   //!
   TBranch        *b_pthtar_e;   //!
   TBranch        *b_pphtar_e;   //!
   TBranch        *b_thtar;   //!
   TBranch        *b_phtar;   //!
   TBranch        *b_ytar;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_xfp;   //!
   TBranch        *b_yfp;   //!
   TBranch        *b_thfp;   //!
   TBranch        *b_phfp;   //!
   TBranch        *b_rxfp;   //!
   TBranch        *b_ryfp;   //!
   TBranch        *b_rthfp;   //!
   TBranch        *b_rphfp;   //!
   TBranch        *b_dxfp;   //!
   TBranch        *b_dyfp;   //!
   TBranch        *b_dthfp;   //!
   TBranch        *b_dphfp;   //!
   TBranch        *b_ntrackFT;   //!
   TBranch        *b_ntrackFPP;   //!
   TBranch        *b_thetaFPP;   //!
   TBranch        *b_phiFPP;   //!
   TBranch        *b_scloseFPP;   //!
   TBranch        *b_zcloseFPP;   //!
   TBranch        *b_xFPP;   //!
   TBranch        *b_yFPP;   //!
   TBranch        *b_xpFPP;   //!
   TBranch        *b_ypFPP;   //!
   TBranch        *b_nhitsFPP;   //!
   TBranch        *b_nhitsFT;   //!
   TBranch        *b_ngoodhitsFPP;   //!
   TBranch        *b_ngoodhitsFT;   //!
   TBranch        *b_chi2ndfFPP;   //!
   TBranch        *b_chi2ndfFT;   //!
   TBranch        *b_xHCAL;   //!
   TBranch        *b_yHCAL;   //!
   TBranch        *b_EHCAL;   //!
   TBranch        *b_tHCAL_ADC;   //!
   TBranch        *b_nblkHCAL;   //!
   TBranch        *b_xECAL;   //!
   TBranch        *b_yECAL;   //!
   TBranch        *b_tECAL_ADC;   //!
   TBranch        *b_nblkECAL;   //!
   TBranch        *b_ECAL_clusblk_e;   //!
   TBranch        *b_ECAL_clusblk_atime;   //!
   TBranch        *b_ECAL_clusblk_again;   //!
   TBranch        *b_ECAL_clusblk_id;   //!
   TBranch        *b_ECAL_clusblk_x;   //!
   TBranch        *b_ECAL_clusblk_y;   //!
   TBranch        *b_HCAL_clusblk_e;   //!
   TBranch        *b_HCAL_clusblk_atime;   //!
   TBranch        *b_HCAL_clusblk_again;   //!
   TBranch        *b_HCAL_clusblk_id;   //!
   TBranch        *b_HCAL_clusblk_x;   //!
   TBranch        *b_HCAL_clusblk_y;   //!
   TBranch        *b_helicity;   //!
   TBranch        *b_IHWP;   //!

   gep_tree_parsed(TTree *tree=0);
   virtual ~gep_tree_parsed();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef gep_tree_parsed_cxx
gep_tree_parsed::gep_tree_parsed(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/volatile/halla/sbs/puckett/gep3parsed_temp_aug6_2025.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/volatile/halla/sbs/puckett/gep3parsed_temp_aug6_2025.root");
      }
      f->GetObject("Tout",tree);

   }
   Init(tree);
}

gep_tree_parsed::~gep_tree_parsed()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gep_tree_parsed::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gep_tree_parsed::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void gep_tree_parsed::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runnum", &runnum, &b_runnum);
   fChain->SetBranchAddress("etheta", &etheta, &b_etheta);
   fChain->SetBranchAddress("ephi", &ephi, &b_ephi);
   fChain->SetBranchAddress("ptheta", &ptheta, &b_ptheta);
   fChain->SetBranchAddress("pphi", &pphi, &b_pphi);
   fChain->SetBranchAddress("ECALO", &ECALO, &b_ECALO);
   fChain->SetBranchAddress("pp", &pp, &b_pp);
   fChain->SetBranchAddress("eprime_eth", &eprime_eth, &b_eprime_eth);
   fChain->SetBranchAddress("pp_eth", &pp_eth, &b_pp_eth);
   fChain->SetBranchAddress("pp_pth", &pp_pth, &b_pp_pth);
   fChain->SetBranchAddress("eprime_pth", &eprime_pth, &b_eprime_pth);
   fChain->SetBranchAddress("eprime_pp", &eprime_pp, &b_eprime_pp);
   fChain->SetBranchAddress("eth_pp", &eth_pp, &b_eth_pp);
   fChain->SetBranchAddress("Q2_pp", &Q2_pp, &b_Q2_pp);
   fChain->SetBranchAddress("Q2_eth", &Q2_eth, &b_Q2_eth);
   fChain->SetBranchAddress("Q2_pth", &Q2_pth, &b_Q2_pth);
   fChain->SetBranchAddress("Q2_e4vect", &Q2_e4vect, &b_Q2_e4vect);
   fChain->SetBranchAddress("eps_pp", &eps_pp, &b_eps_pp);
   fChain->SetBranchAddress("eps_eth", &eps_eth, &b_eps_eth);
   fChain->SetBranchAddress("eps_pth", &eps_pth, &b_eps_pth);
   fChain->SetBranchAddress("eps_4vect", &eps_4vect, &b_eps_4vect);
   fChain->SetBranchAddress("K_eth", &K_eth, &b_K_eth);
   fChain->SetBranchAddress("K_pth", &K_pth, &b_K_pth);
   fChain->SetBranchAddress("K_pp", &K_pp, &b_K_pp);
   fChain->SetBranchAddress("K_p4vect", &K_p4vect, &b_K_p4vect);
   fChain->SetBranchAddress("dpp", &dpp, &b_dpp);
   fChain->SetBranchAddress("dpe", &dpe, &b_dpe);
   fChain->SetBranchAddress("dphi", &dphi, &b_dphi);
   fChain->SetBranchAddress("acopl", &acopl, &b_acopl);
   fChain->SetBranchAddress("inel", &inel, &b_inel);
   fChain->SetBranchAddress("dxECAL", &dxECAL, &b_dxECAL);
   fChain->SetBranchAddress("dyECAL", &dyECAL, &b_dyECAL);
   fChain->SetBranchAddress("dxECAL4vect", &dxECAL4vect, &b_dxECAL4vect);
   fChain->SetBranchAddress("dyECAL4vect", &dyECAL4vect, &b_dyECAL4vect);
   fChain->SetBranchAddress("dtADC", &dtADC, &b_dtADC);
   fChain->SetBranchAddress("pthtar_e", &pthtar_e, &b_pthtar_e);
   fChain->SetBranchAddress("pphtar_e", &pphtar_e, &b_pphtar_e);
   fChain->SetBranchAddress("thtar", &thtar, &b_thtar);
   fChain->SetBranchAddress("phtar", &phtar, &b_phtar);
   fChain->SetBranchAddress("ytar", &ytar, &b_ytar);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("xfp", &xfp, &b_xfp);
   fChain->SetBranchAddress("yfp", &yfp, &b_yfp);
   fChain->SetBranchAddress("thfp", &thfp, &b_thfp);
   fChain->SetBranchAddress("phfp", &phfp, &b_phfp);
   fChain->SetBranchAddress("rxfp", &rxfp, &b_rxfp);
   fChain->SetBranchAddress("ryfp", &ryfp, &b_ryfp);
   fChain->SetBranchAddress("rthfp", &rthfp, &b_rthfp);
   fChain->SetBranchAddress("rphfp", &rphfp, &b_rphfp);
   fChain->SetBranchAddress("dxfp", &dxfp, &b_dxfp);
   fChain->SetBranchAddress("dyfp", &dyfp, &b_dyfp);
   fChain->SetBranchAddress("dthfp", &dthfp, &b_dthfp);
   fChain->SetBranchAddress("dphfp", &dphfp, &b_dphfp);
   fChain->SetBranchAddress("ntrackFT", &ntrackFT, &b_ntrackFT);
   fChain->SetBranchAddress("ntrackFPP", &ntrackFPP, &b_ntrackFPP);
   fChain->SetBranchAddress("thetaFPP", &thetaFPP, &b_thetaFPP);
   fChain->SetBranchAddress("phiFPP", &phiFPP, &b_phiFPP);
   fChain->SetBranchAddress("scloseFPP", &scloseFPP, &b_scloseFPP);
   fChain->SetBranchAddress("zcloseFPP", &zcloseFPP, &b_zcloseFPP);
   fChain->SetBranchAddress("xFPP", &xFPP, &b_xFPP);
   fChain->SetBranchAddress("yFPP", &yFPP, &b_yFPP);
   fChain->SetBranchAddress("xpFPP", &xpFPP, &b_xpFPP);
   fChain->SetBranchAddress("ypFPP", &ypFPP, &b_ypFPP);
   fChain->SetBranchAddress("nhitsFPP", &nhitsFPP, &b_nhitsFPP);
   fChain->SetBranchAddress("nhitsFT", &nhitsFT, &b_nhitsFT);
   fChain->SetBranchAddress("ngoodhitsFPP", &ngoodhitsFPP, &b_ngoodhitsFPP);
   fChain->SetBranchAddress("ngoodhitsFT", &ngoodhitsFT, &b_ngoodhitsFT);
   fChain->SetBranchAddress("chi2ndfFPP", &chi2ndfFPP, &b_chi2ndfFPP);
   fChain->SetBranchAddress("chi2ndfFT", &chi2ndfFT, &b_chi2ndfFT);
   fChain->SetBranchAddress("xHCAL", &xHCAL, &b_xHCAL);
   fChain->SetBranchAddress("yHCAL", &yHCAL, &b_yHCAL);
   fChain->SetBranchAddress("EHCAL", &EHCAL, &b_EHCAL);
   fChain->SetBranchAddress("tHCAL_ADC", &tHCAL_ADC, &b_tHCAL_ADC);
   fChain->SetBranchAddress("nblkHCAL", &nblkHCAL, &b_nblkHCAL);
   fChain->SetBranchAddress("xECAL", &xECAL, &b_xECAL);
   fChain->SetBranchAddress("yECAL", &yECAL, &b_yECAL);
   fChain->SetBranchAddress("tECAL_ADC", &tECAL_ADC, &b_tECAL_ADC);
   fChain->SetBranchAddress("nblkECAL", &nblkECAL, &b_nblkECAL);
   fChain->SetBranchAddress("ECAL_clusblk_e", ECAL_clusblk_e, &b_ECAL_clusblk_e);
   fChain->SetBranchAddress("ECAL_clusblk_atime", ECAL_clusblk_atime, &b_ECAL_clusblk_atime);
   fChain->SetBranchAddress("ECAL_clusblk_again", ECAL_clusblk_again, &b_ECAL_clusblk_again);
   fChain->SetBranchAddress("ECAL_clusblk_id", ECAL_clusblk_id, &b_ECAL_clusblk_id);
   fChain->SetBranchAddress("ECAL_clusblk_x", ECAL_clusblk_x, &b_ECAL_clusblk_x);
   fChain->SetBranchAddress("ECAL_clusblk_y", ECAL_clusblk_y, &b_ECAL_clusblk_y);
   fChain->SetBranchAddress("HCAL_clusblk_e", HCAL_clusblk_e, &b_HCAL_clusblk_e);
   fChain->SetBranchAddress("HCAL_clusblk_atime", HCAL_clusblk_atime, &b_HCAL_clusblk_atime);
   fChain->SetBranchAddress("HCAL_clusblk_again", HCAL_clusblk_again, &b_HCAL_clusblk_again);
   fChain->SetBranchAddress("HCAL_clusblk_id", HCAL_clusblk_id, &b_HCAL_clusblk_id);
   fChain->SetBranchAddress("HCAL_clusblk_x", HCAL_clusblk_x, &b_HCAL_clusblk_x);
   fChain->SetBranchAddress("HCAL_clusblk_y", HCAL_clusblk_y, &b_HCAL_clusblk_y);
   fChain->SetBranchAddress("helicity", &helicity, &b_helicity);
   fChain->SetBranchAddress("IHWP", &IHWP, &b_IHWP);
   Notify();
}

Bool_t gep_tree_parsed::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gep_tree_parsed::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gep_tree_parsed::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gep_tree_parsed_cxx
