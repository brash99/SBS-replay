void hcal_trigAtime(){

  cout << "Processing hcal_trigAtime.." << endl;

  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);

  Int_t Ndata_sbs_trig_adcelemID = 0;
  Double_t bb_sh_nclus = 0., sbs_hcal_nclus = 0.;
  Double_t sbs_trig_a_time[11] = {0.}, sbs_trig_adcelemID[11] = {0.};

  TH1D *h1_hcal_trigAtime = new TH1D("h1_hcal_trigAtime","HCal Trigger ADC Time; ns", 150, 100, 250);

  // Declare branches
  TTree *T = (TTree*) gDirectory->Get("T");
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("sbs.trig.a_time",1);
  T->SetBranchStatus("bb.sh.nclus",1);
  T->SetBranchStatus("sbs.hcal.nclus",1);
  T->SetBranchStatus("Ndata.sbs.trig.adcelemID",1);
  T->SetBranchStatus("sbs.trig.adcelemID",1);

  T->SetBranchAddress("sbs.trig.a_time", &sbs_trig_a_time);
  T->SetBranchAddress("sbs.trig.adcelemID", &sbs_trig_adcelemID);
  T->SetBranchAddress("Ndata.sbs.trig.adcelemID", &Ndata_sbs_trig_adcelemID);
  T->SetBranchAddress("bb.sh.nclus", &bb_sh_nclus);
  T->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus);

  // Loop through events
  Long64_t nevents = T->GetEntries();
  for(Long64_t nevent=0; nevent<nevents; nevent++){

    T->GetEntry(nevent);

    //if( bb_sh_nclus==0 || sbs_hcal_nclus==0 ) continue;

    Double_t hcal_trig_time=0.;
    for(Int_t ihit=0; ihit<Ndata_sbs_trig_adcelemID; ihit++){
      //Fill histo with only the trigger, ignoring the individual sums
      if(sbs_trig_adcelemID[ihit]==0) hcal_trig_time=sbs_trig_a_time[ihit]; 
    }

    h1_hcal_trigAtime->Fill(hcal_trig_time);
  }
  
  h1_hcal_trigAtime->SetStats(0);
  h1_hcal_trigAtime->Draw();

  cout << "Processed macro with " << nevents << " entries." << endl;

  st->Stop();
  cout << "CPU time = " << st->CpuTime() << " s " << " Real time = " << st->RealTime() << " s " << endl;
}
