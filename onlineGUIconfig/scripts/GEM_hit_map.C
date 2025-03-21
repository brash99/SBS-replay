void GEM_hit_map(int layer, TString det = "bb.gem", int what_gep_tracker = 0 /* Input 1 for FT and 2 for FPP */){

  gStyle->SetPalette(1);

  double width = .60;
  double height = 2.0;
  double edge = 0.05;

  if(layer < 4){
    width = 0.4;
    height = 1.5;
  }

  if ( what_gep_tracker == 1 )
  {
    if ( layer < 6 )
    {
      width = 0.4;
      height = 1.5;
    }
    if ( layer == 6 || layer == 7 )
    {
      width = 0.6;
      height = 2.0;
    }
  }
  else if ( what_gep_tracker == 2 )
  {
    width = 0.6;
    height = 2.0;
  }

  TTree *t = (TTree*) gFile->Get("T");

  TH2F *hit2d = new TH2F(Form("hit2d_layer%d",layer),Form("Layer %i;y (m); x (m)",layer),120,-width/2 - edge,width/2 + edge,120,-height/2 - edge,height/2 + edge);

  
  t->Draw( Form("" + det + ".hit.xglobal:" + det + ".hit.yglobal>>hit2d_layer%d",layer),Form("" + det + ".hit.layer == %d&& " + det + ".track.ntrack==1&&" + det + ".track.nhits[" + det + ".hit.trackindex]>=3",layer),"colz" );

  gStyle->SetOptStat(0);  
 
}
