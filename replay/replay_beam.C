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

void replay_beam(const char *codaFilePath,int runNum,unsigned int firstEv,unsigned int lastEv,const char *outfileName){

   // set up file paths

   // output and cut definition files
   // TString replayDir = gSystem->Getenv("SBS-replay");
   // TString odef_path = Form("%s/replay/output_beam_raster.def",replayDir.Data());  
   // TString cdef_path = Form("%s/replay/cuts_beam_raster.def"  ,replayDir.Data());  
   std::string odef_path = "/adaqfs/home/a-onl/sbs/sbs_devel/SBS-replay/replay/output_beam_raster.def"; 
   std::string cdef_path = "/adaqfs/home/a-onl/sbs/sbs_devel/SBS-replay/replay/cuts_beam_raster.def"; 
  
   // output ROOT file destination and name
   // TString out_dir = gSystem->Getenv("OUT_DIR");
   // if( out_dir.IsNull() ) out_dir = ".";
   // TString out_file = out_dir + "/" + Form("replayed_%s_%d.root", filePrefix.Data(),runNum);
   TString out_file = Form("%s",outfileName);
 
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

  // FIXME: Check db_LeftScalevt.dat file!  
  THaScalerEvtHandler *lScaler = new THaScalerEvtHandler("Left","HA scaler event type 140");
  gHaEvtHandlers->Add(lScaler);

  analyzer->SetEvent(event);

  analyzer->SetCompressionLevel(1);
  analyzer->SetOdefFile(odef_path.Data());
  analyzer->SetCutFile(cdef_path.Data());

  std::cout << SCRIPT << "Opening file: " << codaFilePath << std::endl;
  THaRun* run = new THaRun(codaFilePath);

  // FIXME: to address prescale factor issue (temporary) 
  TDatime now; 
  run->SetDate(now); 
  run->SetDataRequired(0); 

  if(firstEv>0) run->SetFirstEvent(firstEv);
  if(lastEv>0)  run->SetLastEvent(lastEv);
  
  analyzer->SetOutFile( out_file.Data() );
  std::cout << SCRIPT << "Output file " << out_file.Data() << " set up " << std::endl; 

  // File to record cuts accounting information
  analyzer->SetSummaryFile("sbs_beam_test.log"); // optional
 
  // run analysis 
  analyzer->SetVerbosity(10);   
  analyzer->Process(*run);

  std::cout << SCRIPT << "Replay of run " << runNum << " COMPLETE!" << std::endl;
 
  // clean up  
  delete analyzer;
  gHaVars->Clear();
  gHaPhysics->Delete();
  gHaApps->Delete();
}
