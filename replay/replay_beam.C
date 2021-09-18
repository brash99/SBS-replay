R__ADD_INCLUDE_PATH($SBS/include)
R__ADD_LIBRARY_PATH($SBS/lib64)
R__ADD_LIBRARY_PATH($SBS/lib)
R__LOAD_LIBRARY(libsbs.so)

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <cstdlib>
#include <iostream>
#include <string> 

#include "TStyle.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSystem.h"

#include "THaRun.h"
#include "THaEvent.h"
#include "THaAnalyzer.h"
#include "THaApparatus.h"
#include "THaVarList.h"
#include "THaInterface.h"
#include "THaHRS.h"
#include "THaDecData.h"

#include "SBSRasteredBeam.h"
#endif 

const std::string SCRIPT = "[replay_beam]: "; 

//______________________________________________________________________________
void replay_beam(const char *filePath,int runNum,unsigned int nev){

   // output and cut definition files 
   std::string odef_path = "./def/output_rastersize.def"; 
   std::string cdef_path = "./def/cuts_rastersize.def"; 

   // TString rootFileName = Form("replayed_beam_%d.root",run);

   // input file path 
   TString filePrefix = Form("lhrs_trigtest");
   TString fullPath   = Form("%s/%s_%d.evio.0",filePath,filePrefix.Data(),runNum);

   // output ROOT file destination and name 
   TString out_dir = gSystem->Getenv("OUT_DIR");
   if( out_dir.IsNull() ) out_dir = ".";
   TString out_file = out_dir + "/" + Form("replayed_%s_%d.root", filePrefix.Data(),runNum);

   THaEvent* event = new THaEvent;

   // Set up the analyzer
   THaAnalyzer* analyzer = THaAnalyzer::GetInstance();
   if(analyzer){
      analyzer->Close();
   }else{
      analyzer = new THaAnalyzer;
   }
 
  // add the LHRS (where the beam signals are)  
  THaHRS* HRSL = new THaHRS("L","Left arm HRS");
  HRSL->AutoStandardDetectors(kFALSE);
  gHaApps->Add( HRSL );

  // add decoder
  THaApparatus* decL = new THaDecData("DL","Misc. Decoder Data");
  gHaApps->Add( decL );

  // add *rastered* beam
  THaApparatus* Lrb = new SBSRasteredBeam("Lrb","Raster Beamline for FADC");
  gHaApps->Add(Lrb);

  // FIXME: add scalers (BCMs...). Need a db_LeftScalevt.dat file!  
  // THaScalerEvtHandler *lScaler = new THaScalerEvtHandler("Left","HA scaler event type 140");
  // gHaEvtHandlers->Add(lScaler);

  analyzer->SetEvent(event);

  analyzer->SetCompressionLevel(1);
  analyzer->SetOdefFile(odef_path.c_str());
  analyzer->SetCutFile(cdef_path.c_str());

  std::cout << SCRIPT << "Opening file: " << fullPath << std::endl;
  THaRun* run = new THaRun(fullPath);
  // run->SetFirstEvent(1);
  if(nev>0) run->SetLastEvent(nev);
  
  analyzer->SetOutFile( out_file.Data() );
  std::cout << SCRIPT << "Output file " << out_file.Data() << " set up " << std::endl; 

  // File to record cuts accounting information
  analyzer->SetSummaryFile("sbs_beam_test.log"); // optional
 
  // run analysis 
  analyzer->SetVerbosity(10);   
  analyzer->Process(*run);
 
  // clean up  
  delete analyzer;
  gHaVars->Clear();
  gHaPhysics->Delete();
  gHaApps->Delete();
}
