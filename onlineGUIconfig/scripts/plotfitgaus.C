#include "TFile.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TPad.h"

void plotfitgaus(TString histname, double fracmax=0.5, int logy=1){
  TH1D *htemp;
  //gDirectory->GetObject(histname,htemp);
  htemp = (TH1D*) (gFile->Get(histname));

  gPad->SetLogy(logy);
  
  //htemp->Print();
  
  if( htemp ){
  
    htemp->Draw();
  
    int binmax = htemp->GetMaximumBin();

    int binlow=binmax, binhigh=binmax;

    double hmax = htemp->GetBinContent(binmax);

    while( htemp->GetBinContent(binlow) >= fracmax * hmax && binlow > 1){binlow--;}
    while( htemp->GetBinContent(binhigh) >= fracmax * hmax && binhigh < htemp->GetNbinsX() ){ binhigh++; }

    gStyle->SetStatW(0.35);
    gStyle->SetStatH(0.3);
  
    gStyle->SetOptFit();

    htemp->Fit("gaus","","",htemp->GetBinCenter(binlow),htemp->GetBinCenter(binhigh));

    
    //htemp->SetStatH(0.3);
  }
}
