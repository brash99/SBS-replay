#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TGraph.h"

void CheckECALHCALGainFactors(const char *rootfilename, const char *gainfilename_ecal, const char *gainfilename_hcal, int nev=500000){

  TChain *C = new TChain("T");
  C->Add(rootfilename);

  TProfile *hgain_chan_ECAL = new TProfile("hgain_chan_ECAL","ECAL;block ID;gain/DB gain",1656,-0.5,1655.5);

  TProfile *hgain_chan_HCAL = new TProfile("hgain_chan_HCAL","HCAL;block ID;gain/DB gain",288,0.5,288.5);

  //  TProfile *hgainratio_chan_ECAL = new TProfile("hgainratio_chan_ECAL","ECAL;block ID;gain/DB gain",1656,-0.5,1655.5,200,0,2.0);

  //TProfile *hgainratio_chan_HCAL = new TH2D("hgainratio_chan_HCAL","HCAL;block ID;gain/DB gain",288,0.5,288.5,200,0,2.0);

  C->Project("hgain_chan_ECAL","earm.ecal.clus_blk.again:earm.ecal.clus_blk.id","earm.ecal.clus_blk.id[0]==earm.ecal.idblk","",nev);
  C->Project("hgain_chan_HCAL","sbs.hcal.clus_blk.again:sbs.hcal.clus_blk.id","sbs.hcal.clus_blk.id[0]==sbs.hcal.idblk","",nev);

  ifstream gainfile_ecal(gainfilename_ecal);
  //if( gainfile_ecal )
  double dbgain_ecal[1656];
  double dbgain_hcal[288];
  
  int nread=0;
  double gaintemp;
  while( gainfile_ecal >> gaintemp && nread < 1656 ){
    dbgain_ecal[nread++] = gaintemp;
  }

  nread = 0;
  ifstream gainfile_hcal(gainfilename_hcal);

  while( gainfile_hcal >> gaintemp && nread < 288 ){
    dbgain_hcal[nread++] = gaintemp;
  }

  double chan_ecal[1656];
  double ratio_ecal[1656];
  double ratio_hcal[288];
  double chan_hcal[288];
  
  for( int i=0; i<1656; i++ ){
    //hgainratio_chan_ECAL->SetBinContent(i+1, hgain_chan_ECAL->GetBinContent(i+1)/dbgain_ecal[i]);
    chan_ecal[i] = i;
    ratio_ecal[i] = hgain_chan_ECAL->GetBinContent(i+1)/dbgain_ecal[i];
    cout << "(ecal chan, gain, db gain, ratio)=(" << i << ", " << hgain_chan_ECAL->GetBinContent(i+1) << ", " << dbgain_ecal[i] << ", " << ratio_ecal[i] << ")" << endl;
  }

  for( int i=0; i<288; i++ ){
    chan_hcal[i] = i+1;
    //    hgainratio_chan_HCAL->SetBinContent(i+1, hgain_chan_HCAL->GetBinContent(i+1)/dbgain_hcal[i]);
    ratio_hcal[i] = hgain_chan_HCAL->GetBinContent(i+1)/dbgain_hcal[i];
    cout << "(hcal chan, gain, db gain, ratio)=(" << i+1 << ", " << hgain_chan_HCAL->GetBinContent(i+1) << ", " << dbgain_hcal[i] << ", " << ratio_hcal[i] << ")" << endl;
  }

  TGraph *gRatio_ECAL = new TGraph(1656,chan_ecal,ratio_ecal);
  TGraph *gRatio_HCAL = new TGraph(288,chan_hcal,ratio_hcal);

  gRatio_ECAL->SetMarkerStyle(20);
  gRatio_HCAL->SetMarkerStyle(20);
  
  TCanvas *c1 = new TCanvas("c1","c1",1600,1200);

  c1->Divide(2,2,0.001,0.001);

  c1->cd(1);
  hgain_chan_ECAL->Draw("colz");
  c1->cd(2);
  //  hgainratio_chan_HCAL->Draw("colz");
  gRatio_ECAL->Draw("AP");
  
  c1->cd(3);
  hgain_chan_HCAL->Draw("colz");
  c1->cd(4);
  //  hgain_chan_HCAL->Draw("colz");
  gRatio_HCAL->Draw("AP");
}
