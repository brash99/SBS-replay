R__ADD_INCLUDE_PATH($SBS/include)
R__ADD_LIBRARY_PATH($SBS/lib64)
R__ADD_LIBRARY_PATH($SBS/lib)
R__LOAD_LIBRARY(libsbs.so)

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <iostream>

#include "TSystem.h"
#include "THaGlobals.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TObject.h"

#include "THaEvData.h"
#include "THaRun.h"
#include "THaAnalyzer.h"
#include "THaVarList.h"
#include "THaInterface.h"

#include "SBSBigBite.h"
#include "SBSBBShower.h"
#include "SBSBBTotalShower.h"
#include "SBSGRINCH.h"
#include "SBSEArm.h"
#include "SBSHCal.h"
#include "SBSGEMStand.h"
#include "SBSTimingHodoscope.h"
#endif

void replay_gmn(int runnum=220, int firstsegment=0, int maxsegments=1, const char *fname_prefix="e1209019_trigtest", long firstevent=0, long nevents=-1, int pedestalmode=0)
{
  SBSBigBite* bigbite = new SBSBigBite("bb", "BigBite spectrometer" );
  //bigbite->AddDetector( new SBSBBShower("ps", "BigBite preshower") );
  //bigbite->AddDetector( new SBSBBShower("sh", "BigBite shower") );
  SBSBBTotalShower* ts= new SBSBBTotalShower("ts", "sh", "ps", "BigBite shower");
  ts->SetDataOutputLevel(0);
  bigbite->AddDetector( ts );
  SBSGenericDetector* bbtrig= new SBSGenericDetector("trig","BigBite shower trig");
  bbtrig->SetModeADC(SBSModeADC::kWaveform);
  bigbite->AddDetector( bbtrig );
  
  bigbite->AddDetector( new SBSGRINCH("grinch", "GRINCH PID") );
  bigbite->AddDetector( new SBSTimingHodoscope("hodotdc", "timing hodo") );
  //bigbite->AddDetector( new SBSGEMSpectrometerTracker("gem", "GEM tracker") );
  SBSGEMSpectrometerTracker *bbgem = new SBSGEMSpectrometerTracker("gem", "BigBite Hall A GEM data");
  bool pm =  ( pedestalmode != 0 );
  //this will override the database setting:
  ( static_cast<SBSGEMTrackerBase *> (bbgem) )->SetPedestalMode( pm );
  bigbite->AddDetector(bbgem);
  gHaApps->Add(bigbite);
    
  SBSEArm *harm = new SBSEArm("sbs","Hadron Arm with HCal");
  SBSHCal* hcal =  new SBSHCal("hcal","HCAL");
  hcal->SetStoreRawHits(kTRUE);
  harm->AddDetector(hcal);

  SBSGenericDetector* sbstrig= new SBSGenericDetector("trig","HCal trigs");
  sbstrig->SetModeADC(SBSModeADC::kWaveform);
  sbstrig->SetStoreRawHits(kTRUE);
  harm->AddDetector( sbstrig );  
  
  gHaApps->Add(harm);
  
  //gHaEvtHandlers->Add (new THaScalerEvtHandler("Left","HA scaler event type 140"));
  //gHaEvtHandlers->Add (new THaScalerEvtHandler("SBS","HA scaler event type 141"));
  
  //bigbite->SetDebug(2);
  //harm->SetDebug(2);

  THaAnalyzer* analyzer = new THaAnalyzer;
  
  //THaInterface::SetDecoder( SBSSimDecoder::Class() );
  THaEvent* event = new THaEvent;

  TString prefix = gSystem->Getenv("DATA_DIR");
  
  bool segmentexists = true;
  int segment=firstsegment; 

  int lastsegment=firstsegment;
  
  TClonesArray *filelist = new TClonesArray("THaRun",10);

  TDatime now = TDatime();
  
  //EPAF: copied the following from replay_BBGEM.C, as this script seems to be thought to handle splits properly.
  int segcounter=0;
  //This loop adds all file segments found to the list of THaRuns to process:
  while( segcounter < maxsegments && segment - firstsegment < maxsegments ){

    TString codafilename;
    //codafilename.Form( "%s/bbgem_%d.evio.%d", prefix.Data(), runnum, segment );
    codafilename.Form("%s/%s_%d.evio.%d", prefix.Data(), fname_prefix, runnum, segment );

    segmentexists = true;
    
    if( gSystem->AccessPathName( codafilename.Data() ) ){
      segmentexists = false;
    } else if( segcounter == 0 ){
      new( (*filelist)[segcounter] ) THaRun( codafilename.Data() );
      cout << "Added segment " << segcounter << ", CODA file name = " << codafilename << endl;
      //( (THaRun*) (*filelist)[segcounter] )->SetDate( now );
      ( (THaRun*) (*filelist)[segcounter] )->SetDataRequired(1);
      ( (THaRun*) (*filelist)[segcounter] )->SetNumber( runnum );
      ( (THaRun*) (*filelist)[segcounter] )->Init();
    } else {
      THaRun *rtemp = ( (THaRun*) (*filelist)[segcounter-1] ); //make otherwise identical copy of previous run in all respects except coda file name:
      new( (*filelist)[segcounter] ) THaRun( *rtemp );
      ( (THaRun*) (*filelist)[segcounter] )->SetFilename( codafilename.Data() );
      ( (THaRun*) (*filelist)[segcounter] )->SetNumber( runnum );
      cout << "Added segment " << segcounter << ", CODA file name = " << codafilename << endl;
    }
    if( segmentexists ){
      segcounter++;
      lastsegment = segment;
    }
    segment++;
  }

  cout << "n segments to analyze = " << segcounter << endl;

  prefix = gSystem->Getenv("OUT_DIR");

  TString outfilename;
  outfilename.Form( "%s/gmn_replayed_%d_seg%d_%d.root", prefix.Data(), runnum,
		    firstsegment, lastsegment );

  // Define the run(s) that we want to analyze.
  // We just set up one, but this could be many.
  //  THaRun* run = new THaRun( "prod12_4100V_TrigRate25_4.dat" );
  //THaRun* run = new THaRun( "5GEM_sample.dat" );
  //THaRun* run = new THaRun( "/Users/puckett/WORK/GEM_ALIGNMENT/RAWDATA/gem_cleanroom_2811.evio.31" );
  //THaRun* run = new THaRun( codafilename.Data() );
  //THaRun* run = new THaRun( "/Users/puckett/WORK/GEM_ALIGNMENT/RAWDATA/gem_cleanroom_2805.evio.0" );

  

  analyzer->SetVerbosity(2);
  analyzer->SetMarkInterval(100);

  analyzer->EnableBenchmarks();
  
  // Define the analysis parameters
  analyzer->SetEvent( event );
  analyzer->SetOutFile( outfilename.Data() );
  // File to record cuts accounting information
  analyzer->SetSummaryFile("replay_gmn.log"); // optional

  prefix = gSystem->Getenv("SBS_REPLAY");
  prefix += "/replay/";

  TString odef_filename = "replay_gmn.odef";
  
  odef_filename.Prepend( prefix );

  analyzer->SetOdefFile( odef_filename );
  
  //analyzer->SetCompressionLevel(0); // turn off compression

  filelist->Compress();

  for( int iseg=0; iseg<filelist->GetEntries(); iseg++ ){
    THaRun *run = ( (THaRun*) (*filelist)[iseg] );
    if( nevents > 0 ) run->SetLastEvent(nevents); //not sure if this will work as we want it to for multiple file segments chained together

    run->SetFirstEvent( firstevent );
    
    run->SetDataRequired(1);
    

    analyzer->Process(run);     // start the actual analysis
  }
  
  /*
  TString run_file = Form("%s.root", filebase);
  if(std::getenv("DATA_DIR")){
    run_file = Form("%s/%s.root", std::string(std::getenv("DATA_DIR")).c_str(), filebase);
  }
  
  if( gSystem->AccessPathName(run_file) ) {
    Error("replay.C", "Input file does not exist: %s", run_file.Data() );
    exit(-1);
  }

  THaRunBase *run = new SBSSimFile(run_file.Data(), "gmn", "");
  run->SetFirstEvent(0);

  cout << "Number of events to replay (-1=all)? ";
  //if( nev > 0 )
  //run->SetFirstEvent(110);
  run->SetLastEvent(nev);
  
  run->SetDataRequired(0);
  run->SetDate(TDatime());
  
  TString out_dir = gSystem->Getenv("OUT_DIR");
  if( out_dir.IsNull() )
    out_dir = ".";
  TString out_file = out_dir + "/" + Form("replayed_%s.root", filebase);
  
  analyzer->SetOutFile( out_file.Data() );
  cout << "output file " << out_file.Data() << " set up " << endl; 
  // File to record cuts accounting information
  analyzer->SetSummaryFile("sbs_hcal_test.log"); // optional

  // Change the cratemap to point to the sim one
  analyzer->SetCrateMapFileName("sbssim_cratemap");

  cout << "sim crate map setup " << endl;
  
  analyzer->SetCutFile( "replay_gmn.cdef" );
  analyzer->SetOdefFile( "replay_gmn.odef" );
  
  cout << "cut file and out file processed " << endl;
  
  analyzer->SetVerbosity(2);  // write cut summary to stdout
  analyzer->EnableBenchmarks();
  
  run->Print();
  
  cout << "about to process " << endl;
  analyzer->Process(run);

  // Clean up

  analyzer->Close();
  delete analyzer;
  //gHaCuts->Clear();
  gHaVars->Clear();
  gHaPhysics->Delete();
  gHaApps->Delete();
  */
  
}


int main(int argc, char *argv[])
{
  new THaInterface( "The Hall A analyzer", &argc, argv, 0, 0, 1 );
  uint runnum = 220;
  uint firstsegment = 0;
  uint maxsegments = 0;
  string name_prefix = "e1209019_trigtest";
  ulong firstevent=0;
  ulong nevents=-1; 
  int pedestalmode=0;
  uint nev = -1;
  if(argc<2 || argc>8){
    cout << "Usage: replay_gmn runnum(uint) firstsegment(uint) maxsegments(uint)" 
	 << endl 
	 << "                  name_prefix(string) firstevent(ulong) nevents(ulong)"
	 << endl 
	 << "                  pedestalmode(int)"
	 << endl;
    return -1;
  }
  if(argc>=2) runnum = atoi(argv[1]);
  if(argc>=3) firstsegment = atoi(argv[2]);
  if(argc>=4) maxsegments = atoi(argv[3]);
  if(argc>=5) name_prefix = argv[4];
  if(argc>=6) firstevent = atoi(argv[5]);
  if(argc>=7) nevents = atoi(argv[6]);
  if(argc>=8) pedestalmode = atoi(argv[7]);

  replay_gmn(runnum, firstsegment, maxsegments, 
	     name_prefix.c_str(), firstevent, nevents, 
	     pedestalmode);
  return 0;
}
