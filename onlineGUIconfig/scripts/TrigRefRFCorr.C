void TrigRefRFCorr(){

  cout << "Processing TrigRefRFCorr.." << endl;

  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);

  Int_t Ndata_bb_tdctrig_tdcelemID = 0;
  Double_t bb_sh_nclus = 0., sbs_hcal_nclus = 0.;
  Double_t bb_tdctrig_tdc[6] = {0.}, bb_tdctrig_tdcelemID[6] = {0.};

  TH1D *h1_TrigRefRFCorr = new TH1D("h1_TrigRefRFCorr","BBCal Trigger RF Corrected; ns", 1000, 499, 999);

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

    Double_t bb_time=0.; //L1A
    Double_t rf_time=0.; //RF time
    for(Int_t ihit=0; ihit<Ndata_bb_tdctrig_tdcelemID; ihit++){
      if(bb_tdctrig_tdcelemID[ihit]==4) rf_time=bb_tdctrig_tdc[ihit];
      if(bb_tdctrig_tdcelemID[ihit]==5) bb_time=bb_tdctrig_tdc[ihit];
    }
    //Double_t rfmod = std::fmod( rf_time, 4. );
    Double_t rfmod = std::fmod( bb_time-rf_time, 4. );

    //cout << rfmod << endl;

    //h1_TrigRefRFCorr->Fill(rfmod);
    h1_TrigRefRFCorr->Fill(bb_time-rfmod);
    //h1_TrigRefRFCorr->Fill(bb_time-rf_time);
  }
  
  h1_TrigRefRFCorr->SetStats(0);
  h1_TrigRefRFCorr->Draw();

  cout << "Processed macro with " << nevents << " entries." << endl;
  cout << "Note: This gives the difference in time between the raw BBCal and HCal triggers." << endl;

  st->Stop();
  cout << "CPU time = " << st->CpuTime() << " s " << " Real time = " << st->RealTime() << " s " << endl;
}
