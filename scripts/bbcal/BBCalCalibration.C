#include <iostream>
#include "TChain.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "gmn_rec_tree.C"

void BBCalCalibration(const char* filename = "replayed_simdigtest_2.root")
{
  TFile *_file0 = TFile::Open(filename);
  TChain *C = (TChain*)_file0->Get("T");
  gmn_rec_tree *T = new gmn_rec_tree(C);
  
  TH1D* h1_DeltaE = new TH1D("h1_DeltaE", ";1-E_{Rec}/p_{MC}", 100, -1, 1);
  
  TH2D* h1_EclusViClus_ps = new TH2D("h1_EclusViClus_ps", "PS;i_{clus};E_clus (GeV)", 15, 0, 15, 140, 0, 7);
  TH2D* h1_EclusViClus_sh = new TH2D("h1_EclusViClus_sh", "SH;i_{clus};E_clus (GeV)", 15, 0, 15, 140, 0, 7);
  
  Long64_t NEvts = C->GetEntries();
  
  //convention: 0-188: SH; 189-240: PS.
  TMatrixD M_jk(241, 241);
  TMatrixD B_j(241, 1);
    
  double A_j[241];
  
  double E_i = 0;
  double A_i = 0;
  
  for(Long64_t nevent = 0; nevent<NEvts; nevent++){
    if( nevent%1000 == 0){
      cout << nevent << "/" << NEvts << endl;
    }
    T->GetEntry(nevent);
    
    E_i = A_i = 0;
    memset(A_j, 0, 241*sizeof(double));
    
    //select good track first:
    //assumption: the best track has the highest momentum
    double pmax = -10;
    int imax = -1;
    for(int i = 0; i<T->bb_tr_n; i++){
      if(T->bb_tr_p[i]>pmax){
	pmax = T->bb_tr_p[i];
	imax = i;
      }
    }
    E_i = pmax;
    // for(int i = 0; i<T->MC_nbbtracks; i++){
    //   if(T->MC_bbtrack_mid[i]==0 && T->MC_bbtrack_pid[i]==11){
    // 	E_i = T->MC_bbtrack_p[i];
    //   }
    // }
    
    //if no clean reference, continue;
    if(E_i==0 || imax<0)continue;

    /*
    for(int i = 0; i<T->bb_sh_nclus; i++){
      h1_EclusViClus_sh->Fill(i, T->bb_sh_clus_e_c[i]);
    }
    for(int i = 0; i<T->bb_sh_nclus; i++){
      h1_EclusViClus_ps->Fill(i, T->bb_ps_clus_e_c[i]);
    }
    A_i=T->bb_ps_clus_e_c[0]+T->bb_sh_clus_e_c[0];
    h1_DeltaE->Fill(1.-A_i/E_i);
    */
    
    double x_tr = T->bb_tr_x[imax]+1.95715*T->bb_tr_th[imax];
    double y_tr = T->bb_tr_y[imax]+1.95715*T->bb_tr_ph[imax];
    
    //cout << T->bb_tr_x[imax] << " " << T->bb_tr_th[imax] << " " << T->bb_tr_p[imax] << endl;
    
    for(int i = 0; i<T->bb_sh_ngoodADChits; i++){
      if(fabs(T->bb_sh_goodblock_x[i]-x_tr)<0.085*2 || fabs(T->bb_sh_goodblock_y[i]-y_tr)<0.085*2)
	A_j[int(T->bb_sh_goodblock_id[i])]+=T->bb_sh_goodblock_e[i];
    }

    x_tr = T->bb_tr_x[imax]+1.75045*T->bb_tr_th[imax];
    y_tr = T->bb_tr_y[imax]+1.75045*T->bb_tr_ph[imax];
    
    for(int i = 0; i<T->bb_ps_ngoodADChits; i++){
      if(fabs(T->bb_ps_goodblock_x[i]-x_tr)<0.09*2 || fabs(T->bb_sh_goodblock_y[i]-y_tr)<0.225)
	 A_j[189+int(T->bb_ps_goodblock_id[i])]+=T->bb_ps_goodblock_e[i];
    }
    
    for(int j = 0; j<241; j++){
      B_j(j, 0)+= A_j[j]*E_i;
      for(int k = 0; k<241; k++){
	M_jk(j, k)+= A_j[j]*A_j[k];
      } 
    }
    
  }
  
  //B_j.Print();
  
  //M_jk.Print();
  
  TMatrixD M_jk_inv = M_jk.Invert();
  TMatrixD C_j = M_jk_inv*B_j;
  
  cout << endl;
  int k = 0;
  //C_j.Print();
  cout << " *** add the following lines in your bigbite shower database *** " << endl;
  cout << endl;
  cout << "bb.sh.adc.gain = " << endl;
  for(int i = 0; i<27; i++){
    for(int j = 0; j<7; j++){
      cout << abs(C_j(k, 0)) << "  ";
      k++;
    }
    cout << endl;
  }
  cout << endl;
  
  cout << " *** add the following lines in your bigbite preshower database *** " << endl;
  cout << endl;
  cout << "bb.ps.adc.gain = " << endl;
  for(int i = 0; i<26; i++){
    for(int j = 0; j<2; j++){
      cout << abs(C_j(k, 0)) << "  ";
      k++;
    }
    cout << endl;
  }
  cout << endl;
  
  /*
  TCanvas* C0 = new TCanvas();
  h1_DeltaE->Draw();
  TCanvas* C1 = new TCanvas("C1", "C1", 800, 400);
  C1->Divide(2, 1);
  C1->cd(1);
  h1_EclusViClus_ps->Draw("colz");
  C1->cd(2);
  h1_EclusViClus_sh->Draw("colz");
  */
}
