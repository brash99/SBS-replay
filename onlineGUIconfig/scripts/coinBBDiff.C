void coinBBDiff(){

  cout << "Processing coinBBDiff macro.." << endl;

  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);

  TPad *pads = (TPad*) gPad;
  pads->Divide(2,2);

  Int_t Ndata_bb_tdctrig_tdcelemID = 0;
  Double_t bb_sh_nclus = 0., sbs_hcal_nclus = 0.;
  Double_t bb_tdctrig_tdc[6] = {0.}, bb_tdctrig_tdcelemID[6] = {0.};

  TH1D *h1_coinBBDiff = new TH1D("h1_coinBBDiff","Coincidence/BBCalTrigger Difference; ns", 1000, -1000, 1000);
  TH1D *h1_coinBBDiff_slow = new TH1D("h1_coinBBDiff_slow","Coincidence/BBCal Difference; ns", 1000, -1000, 1000);
  TH1D *h1_calDiff = new TH1D("h1_calDiff","HCal/BBCal Difference; ns", 250, 400, 650);
  TH1D *h1_RF_time = new TH1D("h1_RF_time","Accelerator RF Time; ns", 1000, -1000, 1000);

  // Declare branches
  TTree *T = (TTree*) gDirectory->Get("T");
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("bb.tdctrig.tdc",1);
  T->SetBranchStatus("bb.sh.nclus",1);
  T->SetBranchStatus("sbs.hcal.nclus",1);
  T->SetBranchStatus("Ndata.bb.tdctrig.tdcelemID",1);
  T->SetBranchStatus("bb.tdctrig.tdcelemID",1);

  T->SetBranchAddress("bb.tdctrig.tdc", &bb_tdctrig_tdc);
  T->SetBranchAddress("bb.tdctrig.tdcelemID", &bb_tdctrig_tdcelemID);
  T->SetBranchAddress("Ndata.bb.tdctrig.tdcelemID", &Ndata_bb_tdctrig_tdcelemID);
  T->SetBranchAddress("bb.sh.nclus", &bb_sh_nclus);
  T->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus);

  // Loop through events
  Long64_t nevents = T->GetEntries();
  for(Long64_t nevent=0; nevent<nevents; nevent++){

    T->GetEntry(nevent);

    //if( bb_sh_nclus==0 || sbs_hcal_nclus==0 ) continue;

    Double_t bb_time=0.;
    Double_t coin_time=0., hcal_time=0, bbslow_time=0, rf_time=0;
    for(Int_t ihit=0; ihit<Ndata_bb_tdctrig_tdcelemID; ihit++){
      if(bb_tdctrig_tdcelemID[ihit]==0) hcal_time=bb_tdctrig_tdc[ihit];
      if(bb_tdctrig_tdcelemID[ihit]==1) coin_time=bb_tdctrig_tdc[ihit];
      if(bb_tdctrig_tdcelemID[ihit]==2) bbslow_time=bb_tdctrig_tdc[ihit];
      if(bb_tdctrig_tdcelemID[ihit]==4) rf_time=bb_tdctrig_tdc[ihit];
      if(bb_tdctrig_tdcelemID[ihit]==5) bb_time=bb_tdctrig_tdc[ihit];
    }
    //Take difference to eliminate 25ns jitter
    h1_coinBBDiff->Fill(coin_time-bb_time);
    h1_coinBBDiff_slow->Fill(coin_time-bbslow_time);
    h1_calDiff->Fill(hcal_time-bb_time);
    h1_RF_time->Fill(rf_time-bb_time);
  }
  
  pads->cd(1);
  h1_coinBBDiff->SetStats(0);
  h1_coinBBDiff->Draw();

  pads->cd(2);
  h1_coinBBDiff_slow->SetStats(0);
  h1_coinBBDiff_slow->Draw();

  pads->cd(3);
  h1_calDiff->SetStats(0);
  h1_calDiff->Draw();

  pads->cd(4);
  h1_RF_time->SetStats(0);
  h1_RF_time->Draw();

  cout << "Processed macro with " << nevents << " entries." << endl;
  cout << "Note: a large spread indicates that HCal is setting the trigger time." << endl;

  st->Stop();
  cout << "CPU time = " << st->CpuTime() << " s " << " Real time = " << st->RealTime() << " s " << endl;
}
