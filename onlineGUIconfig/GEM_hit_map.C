void GEM_hit_map(int layer){

  double width = .60;
  double height = 2.0;
  double edge = 0.05;

  if(layer < 4){
    width = 0.4;
    height = 1.5;
  }

  TTree *t = (TTree*) gFile->Get("T");

  TH2F *hit2d = new TH2F("hit2d",Form("Layer %i;y (m); x (m)",layer),120,-width/2 - edge,width/2 + edge,120,-height/2 - edge,height/2 + edge);

  
  t->Draw("bb.gem.hit.xglobal:bb.gem.hit.yglobal>>hit2d",Form("bb.gem.hit.layer == %i",layer),"colz");

  gStyle->SetOptStat(0);  
 
}
