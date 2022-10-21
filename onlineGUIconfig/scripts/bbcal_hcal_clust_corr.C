/*
  Script to plot the difference between HCAL and BBCAL trigger times. 
*/

void bbcal_hcal_clust_corr() {

  TH1F *h_bbcal_hcal_trigT_diff = new TH1F("h_bbh_trigT_diff","HCal-BBCal Trigger Time Diff. (ns)",200,0,990);
  // TH2F *h2_bbcal_hcal_corr = new TH2F("h2_bbh_corr","BBCal-HCal Cluster Correlation; BB Shower Rows; HCal Rows",27,1,28,24,1,25);

  // Declare trees
  TTree *T = (TTree*) gDirectory->Get("T");
  T->SetBranchStatus("*",0);
  Double_t nclusSH;      T->SetBranchStatus("bb.sh.nclus",1);     T->SetBranchAddress("bb.sh.nclus", &nclusSH);
  Double_t rblkSH;       T->SetBranchStatus("bb.sh.rowblk",1);    T->SetBranchAddress("bb.sh.rowblk", &rblkSH);
  Double_t nclusHCAL;    T->SetBranchStatus("sbs.hcal.nclus",1);  T->SetBranchAddress("sbs.hcal.nclus", &nclusHCAL);
  Double_t rblkHCAL;     T->SetBranchStatus("sbs.hcal.rowblk",1); T->SetBranchAddress("sbs.hcal.rowblk", &rblkHCAL);
  Double_t tdcBBtrig[6]; T->SetBranchStatus("bb.tdctrig.tdc",1);  T->SetBranchAddress("bb.tdctrig.tdc", &tdcBBtrig);
  Int_t NtdcelemBBtrig;  T->SetBranchStatus("Ndata.bb.tdctrig.tdcelemID",1); T->SetBranchAddress("Ndata.bb.tdctrig.tdcelemID", &NtdcelemBBtrig);
  Double_t tdcelemBBtrig[6]; T->SetBranchStatus("bb.tdctrig.tdcelemID",1);   T->SetBranchAddress("bb.tdctrig.tdcelemID", &tdcelemBBtrig);

  // Loop through events
  Long64_t nevents = T->GetEntries();
  for(Long64_t nevent=0; nevent<nevents; nevent++){
    T->GetEntry(nevent);

    if(nclusSH==0 || nclusHCAL==0) continue;

    Double_t bbcal_time=0., hcal_time=0.;
    for(Int_t ihit=0; ihit<NtdcelemBBtrig; ihit++){
      if(tdcelemBBtrig[ihit]==5) bbcal_time=tdcBBtrig[ihit];
      if(tdcelemBBtrig[ihit]==0) hcal_time=tdcBBtrig[ihit];
    }

    Double_t diff = hcal_time - bbcal_time; 
    h_bbcal_hcal_trigT_diff->Fill(diff);
    // if(fabs(diff-510.)<20.){
    //   h2_bbcal_hcal_corr->Fill(rblkSH+1, rblkHCAL+1);
    // }
  }
  
  // h2_bbcal_hcal_corr->SetMinimum(-0.1);
  // h2_bbcal_hcal_corr->SetStats(0);
  // h2_bbcal_hcal_corr->Draw("colz");

  gStyle->SetOptStat(1100);
  h_bbcal_hcal_trigT_diff->Draw();
}
