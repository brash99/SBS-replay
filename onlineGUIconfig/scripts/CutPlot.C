void CutPlot(const char *varexp, const char *cutexp, const char *treename="T"){

  TTree *t = (TTree*) gFile->Get(treename);

  if( t ){
    t->Draw(varexp,cutexp);
  }
}
