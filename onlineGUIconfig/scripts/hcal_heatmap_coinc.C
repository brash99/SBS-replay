//void hcal_heatmap_coinc(Int_t run, Int_t evt){
void hcal_heatmap_coinc(){

  cout<<"Processing hcal_heatmap_coinc"<<endl;
  TStopwatch *st = new TStopwatch();
  st-> Start(kTRUE);

   TH2D *hist_heatmap_cut = new TH2D("hcal_heatmap_coinc","Position of best Cluster: HCAL coincidence cut; Column; Row",12,0,12,24,0,24);

  //TTree *tree = (TTree*) gFile->Get("T");
  TTree *tree =  (TTree*) gDirectory->Get("T");
  tree->SetBranchStatus("*",0);
  Double_t trigbits; tree->SetBranchStatus("g.trigbits",1); tree->SetBranchAddress("g.trigbits",&trigbits);
  Double_t nclus; tree->SetBranchStatus("sbs.hcal.nclus",1); tree->SetBranchAddress("sbs.hcal.nclus",&nclus);
  Double_t rowblk; tree->SetBranchStatus("sbs.hcal.rowblk",1); tree->SetBranchAddress("sbs.hcal.rowblk",&rowblk);
  Double_t colblk; tree->SetBranchStatus("sbs.hcal.colblk",1); tree->SetBranchAddress("sbs.hcal.colblk",&colblk);
  Double_t tdctrig[6]; tree->SetBranchStatus("bb.tdctrig.tdc",1); tree->SetBranchAddress("bb.tdctrig.tdc",&tdctrig);
  Double_t tdctrigID[6]; tree->SetBranchStatus("bb.tdctrig.tdcelemID",1); tree->SetBranchAddress("bb.tdctrig.tdcelemID",&tdctrigID);
  
  Int_t Nentries = tree->GetEntries();
  bool trig_bits=false;
  bool coin_time=false;

  for(int entry=0;entry<Nentries;++entry){
    tree->GetEntry(entry);

    if(trigbits==4 || trigbits==5){
      trig_bits=true;
    }else{ continue;}
    
    if(nclus>0&&trig_bits==true&&abs(tdctrig[2]-tdctrig[0]-50)<20) hist_heatmap_cut->Fill(colblk,rowblk);

  }
  


  /*
  TFile *file = gDirectory->GetFile();
  TTreeReader tree("T",file);

  TTreeReaderValue<double> g_trigbits(tree,"g.trigbits");
  TTreeReaderValue<double> sbs_hcal_nclus(tree,"sbs.hcal.nclus");
  TTreeReaderValue<double> sbs_hcal_rowblk(tree,"sbs.hcal.rowblk");
  TTreeReaderValue<double> sbs_hcal_colblk(tree,"sbs.hcal.colblk");
  TTreeReaderArray<double> bb_tdctrig_tdc(tree,"bb.tdctrig.tdc");
  
  while(tr.Next()){
    if(*g_trigbits==4.0||*g_trigbits==5.0) hist_heatmap->Fill(*sbs_hcal_colblk,*sbs_hcal_rowblk);
    //cout<<*sbs_hcal_colblk<<"\t"<<*sbs_hcal_rowblk<<endl;
  }
  */  
    
    gStyle->SetPalette(53);
    hist_heatmap_cut->Draw("colz");
    hist_heatmap_cut->SetStats(0);
    gPad->SetLogz();
    gPad->SetMargin(0.15,0.15,0.08,0.12);

    st->Stop();

    cout<<"time of processing "<<st->CpuTime() << "s"<<endl;
 }
