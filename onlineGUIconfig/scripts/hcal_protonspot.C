void hcal_protonspot(){

  cout << "Processing hcal_protonspot.." << endl;

  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);
  
  Double_t sbs_hcal_clus_blk_e[30] = {0.}, sbs_hcal_clus_blk_atime[30] = {0.}, sbs_hcal_nclus = 0., sbs_hcal_ledbit = -2., sbs_hcal_colblk=-1., sbs_hcal_rowblk=-1., e_kine_W2=0., sbs_hcal_y=-10., sbs_hcal_x=-10., bb_tr_n=-1.; 
  UInt_t fEvtHdr_fTrigBits = -1;

  TH2D *h2_hcal_N_XvY_Cuts = new TH2D("h2_hcal_N_XvY_Cuts","HCal X vs Y with Elastic Cuts; Y (col N); X (row N)", 12, -0.86, 0.86, 24, -2.2, 1.5 );

  // Declare branches
  TTree *T = (TTree*) gDirectory->Get("T");
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("sbs.hcal.clus_blk.atime",1);
  T->SetBranchStatus("sbs.hcal.clus_blk.e",1);
  T->SetBranchStatus("sbs.hcal.nclus",1);
  T->SetBranchStatus("sbs.hcal.ledbit",1);
  T->SetBranchStatus("sbs.hcal.colblk",1);
  T->SetBranchStatus("sbs.hcal.rowblk",1);
  T->SetBranchStatus("sbs.hcal.y",1);
  T->SetBranchStatus("sbs.hcal.x",1);
  T->SetBranchStatus("e.kine.W2",1);
  T->SetBranchStatus("bb.tr.n", 1);
  T->SetBranchStatus("fEvtHdr.fTrigBits",1);

  T->SetBranchAddress("sbs.hcal.clus_blk.atime", &sbs_hcal_clus_blk_atime);
  T->SetBranchAddress("sbs.hcal.clus_blk.e", &sbs_hcal_clus_blk_e);
  T->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus);
  T->SetBranchAddress("sbs.hcal.ledbit", &sbs_hcal_ledbit);
  T->SetBranchAddress("sbs.hcal.colblk", &sbs_hcal_colblk);
  T->SetBranchAddress("sbs.hcal.rowblk", &sbs_hcal_rowblk);
  T->SetBranchAddress("sbs.hcal.y", &sbs_hcal_y);
  T->SetBranchAddress("sbs.hcal.x", &sbs_hcal_x);
  T->SetBranchAddress("e.kine.W2", &e_kine_W2);
  T->SetBranchAddress("bb.tr.n", &bb_tr_n);
  T->SetBranchAddress("fEvtHdr.fTrigBits", &fEvtHdr_fTrigBits);

  // Acquire the number of entries
  Long64_t nevents = T->GetEntries();
  
  for(Long64_t nevent=0; nevent<nevents; nevent++){

    T->GetEntry(nevent);

    //Use if BB optics are trustworthy
    //if( e_kine_W2<0.8 || e_kine_W2>1.2 || sbs_hcal_clus_blk_atime[0]<40 || sbs_hcal_clus_blk_atime[0]>80 || fEvtHdr_fTrigBits==32 || sbs_hcal_nclus==0 || bb_tr_n!=1 ) continue;
    //Use if BB optics are NOT trustworthy
    if( sbs_hcal_clus_blk_atime[0]<40 || sbs_hcal_clus_blk_atime[0]>80 || fEvtHdr_fTrigBits==32 || sbs_hcal_nclus==0 ) continue;

    h2_hcal_N_XvY_Cuts->Fill( sbs_hcal_y, sbs_hcal_x );

  }
  
  //h2_hcal_N_XvY_Cuts->SetMinimum(-0.1);
  h2_hcal_N_XvY_Cuts->SetStats(0);
  h2_hcal_N_XvY_Cuts->Draw("colz");

  cout << "Processed macro with " << nevents << " entries." << endl;
  
  st->Stop();
  cout << "CPU time = " << st->CpuTime() << " s " << " Real time = " << st->RealTime() << " s " << endl;
}
