void coinTRIGdiff(){

  cout << "Processing coinTRIGdiff macro.." << endl;

  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);

  TPad *pads = (TPad*) gPad;
  pads->Divide(3,2);

  Double_t sbs_hcal_nclus = 0.;
  Double_t hcal_reftdc[7];
  Double_t hcal_reftdc_elemID[7];

  Int_t Ndata_hcal_reftdc=0;
  
  TH1D *h1_ECAL = new TH1D("h1_ECAL","ECAL Trig time; ns", 400, 0., 400.);
  TH1D *h1_HCAL = new TH1D("h1_HCAL","HCAL Trig time; ns", 400, 0., 400.);
  TH1D *h1_COIN = new TH1D("h1_COIN","COIN Trig time; ns", 400, 0., 400.);
  TH1D *h1_COINmatch = new TH1D("h1_COINmatch","COIN-Match Trig time; ns", 400, 0., 400.);
  TH1D *h1_calDiff = new TH1D("h1_calDiff","HCal/ECAL Difference; ns", 500, -250, 250);
  TH1D *h1_RF_time = new TH1D("h1_RF_time","Accelerator RF Time; ns;",320, 0.,320.);

  double trigbits;
  
  // Declare branches
  TTree *T = (TTree*) gDirectory->Get("T");
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("sbs.hcal.Ref.tdc",1);
  T->SetBranchStatus("sbs.hcal.Ref.tdcelemID",1);
  T->SetBranchStatus("Ndata.sbs.hcal.Ref.tdcelemID",1);
  T->SetBranchStatus("g.trigbits",1);
  
  T->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus);
  T->SetBranchAddress("sbs.hcal.Ref.tdc",hcal_reftdc);
  T->SetBranchAddress("sbs.hcal.Ref.tdcelemID",hcal_reftdc_elemID);
  T->SetBranchAddress("Ndata.sbs.hcal.Ref.tdcelemID",&Ndata_hcal_reftdc);
  T->SetBranchAddress("g.trigbits", &trigbits);
  
  // Loop through events
  Long64_t nevents = T->GetEntries();

  int nbadRF_HCAL = 0;

  for(Long64_t nevent=0; nevent<nevents; nevent++){

    T->GetEntry(nevent);

    //if( bb_sh_nclus==0 || sbs_hcal_nclus==0 ) continue;


    Double_t coin_time=-1., coinMatch_time=-1, hcal_time=-1., ecal_time=-1., rf_time=-1.;

    bool anyHCALRF = false;
    
    for( int ihit=0; ihit<Ndata_hcal_reftdc; ihit++ ){
      if( hcal_reftdc_elemID[ihit]==3 ){
	rf_time = hcal_reftdc[ihit];
	anyHCALRF = true;
      }
      if( hcal_reftdc_elemID[ihit]==4 ){ //coin-Match trigger time:
	coinMatch_time = hcal_reftdc[ihit];
      }
      if( hcal_reftdc_elemID[ihit]==5 ){ //coin trigger time:
	coin_time = hcal_reftdc[ihit];
      }
      if( hcal_reftdc_elemID[ihit]==6 ){ //ecal trigger time:
	ecal_time = hcal_reftdc[ihit];
      }
      if( hcal_reftdc_elemID[ihit]==7 ){ //hcal trigger time:
	hcal_time = hcal_reftdc[ihit];
      }
    }
    
    if( !anyHCALRF ) nbadRF_HCAL++;
    
    //Take difference to eliminate 25ns jitter

    //if( trigbits == 4 ){
      if(ecal_time>=0) h1_ECAL->Fill(ecal_time);
      if(hcal_time>=0) h1_HCAL->Fill(hcal_time);
      if(ecal_time>=0 && hcal_time>=0) h1_calDiff->Fill(ecal_time-hcal_time);
      if(coin_time>=0) h1_COIN->Fill(coin_time);
      if(coinMatch_time>=0) h1_COINmatch->Fill(coinMatch_time);
      
      h1_RF_time->Fill(rf_time);
      
    //}
  }
  
  pads->cd(1);//->SetLogy();
  //h1_ECAL->SetStats(0);
  h1_ECAL->Draw();

  pads->cd(2);//->SetLogy();
  //h1_HCAL->SetStats(0);
  h1_HCAL->Draw();

  pads->cd(3);//->SetLogy();
  //h1_calDiff->SetStats(0);
  h1_calDiff->Draw();

  TString stitle;
  
  pads->cd(4);//->SetLogy();
  //h1_COIN->SetStats(0);
  h1_COIN->Draw();

  pads->cd(5);//->SetLogy();
  //h1_COINmatch->SetStats(0);
  h1_COINmatch->Draw();

  pads->cd(6);//->SetLogy();
  h1_RF_time->SetTitle(stitle.Format("RF time HCAL (N bad = %d);ns;events", nbadRF_HCAL));
  //h1_RF_time->SetStats(0);
  h1_RF_time->Draw();

 
 
  cout << "Processed macro with " << nevents << " entries." << endl;
  //cout << "Note: a large spread indicates that HCal is setting the trigger time." << endl;

  st->Stop();
  cout << "CPU time = " << st->CpuTime() << " s " << " Real time = " << st->RealTime() << " s " << endl;
}
