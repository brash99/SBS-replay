// R__ADD_INCLUDE_PATH($SBS/include)
// R__ADD_LIBRARY_PATH($SBS/lib64)
// R__ADD_LIBRARY_PATH($SBS/lib)
// R__LOAD_LIBRARY(libsbs.so)

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
#include "SBSScalerEvtHandler.h"
#include "LHRSScalerEvtHandler.h"

const std::string SCRIPT = "[replay_beam]: "; 

// void replay_beam(const char *codaFilePath,int runNum,unsigned int firstEv,unsigned int lastEv,const char *outfileName){
void replay_beam(int runNum,Long_t firstevent=0,Long_t nevents=-1,int maxsegments=1){

   std::cout << "=============================" << std::endl;
   std::cout << SCRIPT << "Input parameters: "  << std::endl;
   std::cout << "Run:          " << runNum      << std::endl;
   std::cout << "First event:  " << firstevent  << std::endl;
   std::cout << "N events:     " << nevents     << std::endl;
   std::cout << "max segments: " << maxsegments << std::endl;
   std::cout << "=============================" << std::endl;

   // Set up the analyzer
   THaAnalyzer* analyzer = THaAnalyzer::GetInstance();
   if(analyzer){
      analyzer->Close();
   }else{
      analyzer = new THaAnalyzer;
   }

   // unsigned int firstEv = 1; 
   // unsigned int lastEv  = 50000;
  
   TString prefix = gSystem->Getenv("DATA_DIR");

   // find all files associated with the run 
   char fname_prefix[200]; 
   sprintf(fname_prefix,"e1209019"); 

   int stream=0;

   TClonesArray *filelist = new TClonesArray("THaRun",10);

   // search paths to find the raw data
   std::vector<TString> pathlist;
   pathlist.push_back( prefix );

   if(prefix!="/adaqeb1/data1" )  pathlist.push_back( "/adaqeb1/data1" );
   if(prefix!="/adaq1/data1/sbs") pathlist.push_back( "/adaq1/data1/sbs" );
   if(prefix!="/cache/mss/halla/sbs/raw") pathlist.push_back( "/cache/mss/halla/sbs/raw" );

   for( int i=0; i<pathlist.size(); i++ ){
      std::cout << "search paths = " << pathlist[i] << std::endl;
   }

   // int maxsegments=1;
   int max1 = maxsegments;
   int segcounter=0,segment=0,firstsegment=0,lastsegment=0; 

   bool segmentexists;

   //This loop adds all file segments found to the list of THaRuns to process:
   while( segcounter < max1 && segment - firstsegment < maxsegments ){
      TString codafilename;
      //codafilename.Form( "%s/bbgem_%d.evio.%d", prefix.Data(), runNum, segment );
      codafilename.Form("%s_%d.evio.%d.%d", fname_prefix, runNum, stream, segment );
      
      TString ftest(fname_prefix);
      if( ftest == "bbgem" || ftest == "e1209019_trigtest" ){
              codafilename.Form("%s_%d.evio.%d", fname_prefix, runNum, segment );
      }
      
      segmentexists = false;
      
      cout << "codafilename = " << codafilename << endl;
      
      for( int ipath=0; ipath<pathlist.size(); ipath++ ){
         TString searchname;
         searchname.Form( "%s/%s", pathlist[ipath].Data(), codafilename.Data() );
         if( !gSystem->AccessPathName( searchname.Data() ) ){
            segmentexists = true;
            break;
         }
      }
      
      if( segmentexists ){
         new( (*filelist)[segcounter] ) THaRun( pathlist, codafilename.Data(), "GMN run" );
         cout << "Added segment " << segment << ", CODA file name = " << codafilename << endl;
      }
      if( segmentexists ){
         segcounter++;
         lastsegment = segment;
      }
      segment++;
   }

   std::cout << "n segments to analyze = " << segcounter << endl;
   
   // set up file paths
   prefix = gSystem->Getenv("OUT_DIR");
   
   TString outfilename = Form("%s/gmn_replayed-beam_%d_stream%d_seg%d_%d.root",prefix.Data(),runNum,
		              stream,firstsegment,lastsegment); 

   prefix = gSystem->Getenv("LOG_DIR");
   analyzer->SetSummaryFile(Form("%s/replay_beam_gmn.log", prefix.Data()));

   prefix = gSystem->Getenv("SBS_REPLAY");
   prefix += "/replay/";

   TString odef_filename = Form("%s/output_beam_both.def",prefix.Data());
   TString cdef_filename = Form("%s/cuts_beam_raster.def",prefix.Data());

   odef_filename.Prepend( prefix );

   analyzer->SetOdefFile(odef_filename);
   analyzer->SetCutFile(cdef_filename);
   
   THaEvent* event = new THaEvent;

   // add the LHRS (where the beam signals are)  
   THaHRS* HRSL = new THaHRS("L","Left arm HRS");
   HRSL->AutoStandardDetectors(kFALSE);
   gHaApps->Add( HRSL );

   // add decoder
   THaApparatus* decL = new THaDecData("DL","Misc. Decoder Data");
   gHaApps->Add( decL );

   // add *rastered* beam (LHRS)
   THaApparatus* Lrb = new SBSRasteredBeam("Lrb","Raster Beamline for FADC");
   gHaApps->Add(Lrb);

   // add *rastered* beam (SBS)
   THaApparatus* sbs = new SBSRasteredBeam("SBSrb","Raster Beamline for FADC");
   gHaApps->Add(sbs);

   std::ofstream debugFile; 
   debugFile.open("lhrs-scaler-dump.txt");

   // LHRS scaler data (ROC10)  
   LHRSScalerEvtHandler *lScaler = new LHRSScalerEvtHandler("Left","HA scaler event type 140");
   lScaler->SetDebugFile(&debugFile);
   gHaEvtHandlers->Add(lScaler);

   // SBS scaler data (sbsvme29) 
   SBSScalerEvtHandler *sbsScaler = new SBSScalerEvtHandler("sbs","SBS Scaler Bank event type 1");
   // sbsScaler->AddEvType(1);             // Repeat for each event type with scaler banks
   sbsScaler->SetUseFirstEvent(kTRUE);
   gHaEvtHandlers->Add(sbsScaler);

   analyzer->SetEvent(event);

   analyzer->SetCompressionLevel(1);
   // analyzer->SetOdefFile(odef_path.Data());
   // analyzer->SetCutFile(cdef_path.Data());
   
   analyzer->SetOutFile( outfilename.Data() );
   std::cout << SCRIPT << "Output file " << outfilename.Data() << " set up " << std::endl; 

   // File to record cuts accounting information
   analyzer->SetSummaryFile("sbs_beam_test.log"); // optional

   analyzer->SetVerbosity(10);   

   filelist->Compress();

   for(int iseg=0; iseg<filelist->GetEntries(); iseg++ ){
      std::cout << SCRIPT << "Trying file segment " << iseg << std::endl;
      THaRun *run = ( (THaRun*) (*filelist)[iseg] );
      if(nevents>0) run->SetLastEvent(nevents); //not sure if this will work as we want it to for multiple file segments chained together
      run->SetFirstEvent( firstevent );
      run->SetDataRequired(THaRunBase::kDate|THaRunBase::kRunNumber);
      // FIXME: to address prescale factor issue (temporary) 
      // TDatime now; 
      // run->SetDate(now); 
      // run->SetDataRequired(0); 
      if( run->GetSegment() >= firstsegment && run->GetSegment() - firstsegment < maxsegments ){
         analyzer->Process(run);     // start the actual analysis
         std::cout << SCRIPT << "--> Done!" << std::endl;
      }
   }

   std::cout << SCRIPT << "Replay of run " << runNum << " COMPLETE!" << std::endl;

   // clean up  
   delete analyzer;
   gHaVars->Clear();
   gHaPhysics->Delete();
   gHaApps->Delete();
}
