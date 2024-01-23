void plotHCALeff( TH1D *heff, double xmin=-3, double xmax=1.1 ){
  //heff->GetXaxis()->SetRangeUser( xmin - 0.2*(xmax-xmin), 
  //				  xmax + 0.2*(xmax-xmin) );

  gStyle->SetOptFit(0); 
  gStyle->SetOptStat(0);
  
  TFitResultPtr fitresult = heff->Fit( "pol0", "SQ", "", xmin, xmax );
  
  double eff = fitresult->Parameter(0);
  double efferr = (fitresult->GetErrors())[0] * sqrt( fitresult->Chi2()/fitresult->Ndf());

  gPad->SetGrid();

  heff->Draw();

  TString newtitle; 

  newtitle.Form( "Avg. Efficiency = %6.4g #pm %6.4g", eff, efferr );
  heff->SetTitle(newtitle.Data());

  heff->GetYaxis()->SetTitle( "HCAL proton efficiency" );


}
