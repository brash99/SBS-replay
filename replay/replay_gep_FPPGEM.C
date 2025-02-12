#include "TSystem.h"
#include "TList.h"
#include "THaRun.h"
#include "THaEvent.h"
#include "THaAnalyzer.h"
#include "THaApparatus.h"
#include "TString.h"
#include "TClonesArray.h"
#include "SBSGenericDetector.h"
#include "SBSBBTotalShower.h"
#include "SBSBBShower.h"


#include "SBSGEMSpectrometerTracker.h"
#include "SBSBigBite.h"


// This replay script will only replay data from the GEp front trakcer i.e. first eight layers.


void replay_gep_FPPGEM( UInt_t runnum=1234, Long_t nevents=5000, Long_t firstevent=0, const char *fname_prefix="gep5", UInt_t firstsegment=0, UInt_t maxsegments=1, Int_t pedestalmode=0, Int_t cmplots=1, Int_t parallelemode=0, Int_t singlemode=0, Int_t file_seg=0 )
{	
	//  gSystem->Load("libsbs.so");

 	int stream = 0;

 	if(pedestalmode == 1)
 	{
	    cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
	    cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;
	    cout << "-----Running in pedestal mode-----" << endl;
	    cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;
	    cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
  }

  SBSBigBite *sbs = new SBSBigBite("sbs", "Generic apparatus");

  SBSGEMSpectrometerTracker *sbsgem = new SBSGEMSpectrometerTracker("gemFPP", "SBS GEp focal plane polarimeter GEM data"); 

  sbs->AddDetector(sbsgem);

  //this will override the database setting:
  ( static_cast<SBSGEMTrackerBase *> (sbsgem) )->SetPedestalMode( pedestalmode );
  ( static_cast<SBSGEMTrackerBase *> (sbsgem) )->SetMakeCommonModePlots( cmplots );

  // Set up the analyzer - we use the standard one,
  // but this could be an experiment-specific one as well.
  // The Analyzer controls the reading of the data, executes
  // tests/cuts, loops over Apparatus's and PhysicsModules,
  // and executes the output routines.
  THaAnalyzer* analyzer = new THaAnalyzer;

  gHaApps->Add(sbs);


  // A simple event class to be output to the resulting tree.
  // Creating your own descendant of THaEvent is one way of
  // defining and controlling the output.
  THaEvent* event = new THaEvent;

  TString prefix = gSystem->Getenv("DATA_DIR");

  bool segmentexists = true;
  int segment = firstsegment;
  int lastsegment = firstsegment;

  TClonesArray *filelist = new TClonesArray("THaRun",10);

  TDatime now = TDatime();

  vector<TString> pathlist;

  pathlist.push_back(prefix);

  // 
  // 
  //
  if( prefix != "/adaqeb2/data1" )
   pathlist.push_back( "/adaqeb2/data1" );

  if( prefix != "/data/raw" )
   pathlist.push_back( "/data/raw" );

  if( prefix != "/adaq1/data1/sbs" )
  pathlist.push_back( "/adaq1/data1/sbs" );

  if( prefix != "/cache/halla/sbs/raw" )
  pathlist.push_back( "/cache/halla/sbs/raw" );

	if( prefix != "/cache/halla/sbs/GEnRP/raw" )
	pathlist.push_back("/cache/halla/sbs/GEnRP/raw");

  if( prefix != "/cache/halla/sbs/GEp/raw" )
	pathlist.push_back("/cache/halla/sbs/GEp/raw");

  for( int i=0; i<pathlist.size(); i++ ){
   	cout << "search paths = " << pathlist[i] << endl;
  }

  TDatime RunDate = TDatime(); 

 	int max1 = maxsegments;

  int segcounter=0;
    	
  //This loop adds all file segments found to the list of THaRuns to process:
  while( segcounter < max1 && segment - firstsegment < maxsegments )
  {
		TString codafilename;
    
	  TString ftest(fname_prefix);

	  if( ftest == "gep5" )
	  {
		 	codafilename.Form("%s_%d.evio.0.%d", fname_prefix, runnum, segment);
	  }
    else if ( ftest == "gem" )
    {
      codafilename.Form("%s_%d.dat.0.%d", fname_prefix, runnum, segment);

    }
    	
    segmentexists = false;

	  cout << "codafilename = " << codafilename << endl;

	  for( int ipath=0; ipath<pathlist.size(); ipath++ )
	  {
	    TString searchname;
	    searchname.Form( "%s/%s", pathlist[ipath].Data(), codafilename.Data() );
	    cout << "---------********-----------********-------------" << endl;
	    cout << " Search name: " << searchname << endl;
	    cout << "---------********-----------********-------------" << endl;
	      
	  	if( !gSystem->AccessPathName( searchname.Data() ) )
	  	{
				segmentexists = true;
				break;
	    }
	  }
   
	  if( segmentexists )
	  {
	    new( (*filelist)[segcounter] ) THaRun( pathlist, codafilename.Data(), "GEp run" );
	    cout << "Added segment " << segment << ", CODA file name = " << codafilename << endl;

	      //( (THaRun*) (*filelist)[segcounter] )->SetDate( now );

	      //if( stream == 0 && segment == 0 ){
	      //( (THaRun*) (*filelist)[segcounter] )->SetDataRequired(THaRunBase::kDate|THaRunBase::kRunNumber);
	      //( (THaRun*) (*filelist)[segcounter] )->Init();

	      //RunDate = ( (THaRun*) (*filelist)[segcounter] )->GetDate();
	      //} else {
	      //	( (THaRun*) (*filelist)[segcounter] )->SetDataRequired(0);

	      //cout << "Warning: setting date to " << RunDate.AsString() << " for stream " << stream << " segment " << segment 
	      //     << endl; 

	      //	( (THaRun*) (*filelist)[segcounter] )->SetDate(RunDate);
	      //	( (THaRun*) (*filelist)[segcounter] )->SetNumber(UInt_t(runnum));
	      //}
	      //( (THaRun*) (*filelist)[segcounter] )->SetNumber( runnum );
	      //( (THaRun*) (*filelist)[segcounter] )->Init();
	    } // else {
	    //   THaRun *rtemp = ( (THaRun*) (*filelist)[segcounter-1] ); //make otherwise identical copy of previous run in all respects except coda file name:
	    //   new( (*filelist)[segcounter] ) THaRun( *rtemp );
	    //   ( (THaRun*) (*filelist)[segcounter] )->SetFilename( codafilename.Data() );
	    //   ( (THaRun*) (*filelist)[segcounter] )->SetNumber( runnum );
	    //   cout << "Added segment " << segcounter << ", CODA file name = " << codafilename << endl;
	    // }
	    if( segmentexists )
	    {
	      segcounter++;
	      lastsegment = segment;
	    }

	    segment++;
  	}

  cout << "n segments to analyze = " << segcounter << endl;

  prefix = gSystem->Getenv("OUT_DIR");

  TString outfilebase, outfilename;

  if( nevents > 0 ){ 

    outfilebase.Form( "%s_replayed_gepFPPGEMs_%d_seg%d_%d_firstevent%d_nevent%d", fname_prefix, runnum, firstsegment, lastsegment, firstevent, nevents );
    outfilename.Form( "%s/%s.root", prefix.Data(), outfilebase.Data());
  }
  else {
    outfilebase.Form( "%s_fullreplay_gepFPPGEMs_%d_seg%d_%d", fname_prefix, runnum, firstsegment, lastsegment );
    outfilename.Form( "%s/%s.root", prefix.Data(), outfilebase.Data());
  }

  // TString outfilename;
  // outfilename.Form( "%s/%s_SBSGEMs_%d.evio.0.%d_%d.root", prefix.Data(), fname_prefix, runnum, firstsegment, lastsegment );
  /*TString fntest(fname_prefix);
  if( fntest == "eel125"){
    outfilename.Form( "%s/%s_replayed_%d_newDB.root", prefix.Data(), fname_prefix, runnum);
  }*/
  
 	//TString fntest(fname_prefix);
  // 	if( fntest == "eel125"){
  //   	outfilename.Form( "%s/%s_replayed_%d_seg%d_%d.root", prefix.Data(), fname_prefix, runnum, firstsegment, lastsegment );
  //   } 
  // 	if( parallelmode == 1 && fntest == "eel125"){

	 //    if( singlemode == 1 ){
	 //      outfilename.Form("%s/%s_replayed_%d_seg_%d.root", prefix.Data(), fname_prefix, runnum, file_seg);
	 //    }
	 //    else{
		//     outfilename.Form("%s/%s_replayed_%d_seg_%d.root", prefix.Data(), fname_prefix, runnum, firstsegment);
	 //    }	
  // 	}

  // if( fntest == "e1209016"){
  //   outfilename.Form( "%s/%s_SBSGEMs_%d.evio.%d_%d.root", prefix.Data(), fname_prefix, runnum, firstsegment, lastsegment );
  //   } 
  // if( parallelmode == 1 && fntest == "e1209016"){

  //   if( singlemode == 1 ){
  //     outfilename.Form("%s/%s_SBSGEMs_%d.evio.0.%d.root", prefix.Data(), fname_prefix, runnum, file_seg);
  //   }
  //   else{
  //     outfilename.Form("%s/%s_SBSGEMs_%d.evio.0.%d.root", prefix.Data(), fname_prefix, runnum, firstsegment);
  //   } 
  // }


  // if( fntest == "e1217004" )
  // {
  // 	outfilename.Form("%s/%s_replayed_SBSGEMs_%d_seg_%d_%d.root", prefix.Data(), fname_prefix, runnum, firstsegment, lastsegment);
  // }

  analyzer->SetVerbosity(2);
  analyzer->SetMarkInterval(100);

  analyzer->EnableBenchmarks();
  
  // Define the analysis parameters
  analyzer->SetEvent( event );
  analyzer->SetOutFile( outfilename.Data() );
  // File to record cuts accounting information

  prefix = gSystem->Getenv("LOG_DIR");
  analyzer->SetSummaryFile(Form("%s/%s.log", prefix.Data(), outfilebase.Data())); // optional

  prefix = gSystem->Getenv("SBS_REPLAY");
  prefix += "/replay/";

  TString odef_filename = "replay_FPPGEM_gep.odef";
  
  odef_filename.Prepend( prefix );

  analyzer->SetOdefFile( odef_filename );
  
  //analyzer->SetCompressionLevel(0); // turn off compression

  filelist->Compress();

  for( int iseg=0; iseg<filelist->GetEntries(); iseg++ ){
    THaRun *run = ( (THaRun*) (*filelist)[iseg] );

    run->Init();

    if( nevents > 0 ) run->SetLastEvent(nevents); //not sure if this will work as we want it to for multiple file segments chained together

    run->SetFirstEvent( firstevent );
    
    run->SetDataRequired(THaRunBase::kDate|THaRunBase::kRunNumber);

    // run->SetDataRequired(0);
    // run->SetDate(now);
    
    if( run->GetSegment() >= firstsegment && run->GetSegment() - firstsegment < maxsegments ){
      analyzer->Process(run);     // start the actual analysis
    }
  }

}