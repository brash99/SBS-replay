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

#include "SBSEArm.h"
#include "SBSGEMSpectrometerTracker.h"
#include "SBSBigBite.h"

using namespace std;

// This replay script will do tracking analysis for both the GEp FT and FPP trackers.

void replay_gep_allGEM( UInt_t runnum=1234, Long_t nevents=5000, Long_t firstevent=0, const char *fname_prefix="gep5", UInt_t firstsegment=0, UInt_t maxsegments=1, Int_t pedestalmode=0, Int_t cmplots=1, Int_t parallelemode=0, Int_t singlemode=0, Int_t file_seg=0 )
{	
    int stream = 0;

 	if(pedestalmode == 1)
 	{
	    cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
	    cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;
	    cout << "-----Running in pedestal mode-----" << endl;
	    cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;
	    cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
    }

    SBSBigBite *harm = new SBSBigBite("sbs", "GEp Hadron Arm with FT FPP and HCal");

    SBSGEMSpectrometerTracker *gemFT = new SBSGEMSpectrometerTracker("gemFT", "SBS GEp front tracker GEM data"); 
    ( static_cast<SBSGEMTrackerBase *> (gemFT) )->SetPedestalMode( pedestalmode );
    ( static_cast<SBSGEMTrackerBase *> (gemFT) )->SetMakeCommonModePlots( cmplots );
    harm->AddDetector(gemFT);

    SBSGEMSpectrometerTracker *gemFPP = new SBSGEMSpectrometerTracker("gemFPP", "SBS GEp polarimeter tracker GEM data"); 
    ( static_cast<SBSGEMTrackerBase *> (gemFPP) )->SetPedestalMode( pedestalmode );
    ( static_cast<SBSGEMTrackerBase *> (gemFPP) )->SetMakeCommonModePlots( cmplots );
    harm->AddDetector(gemFPP);
        
    THaAnalyzer* analyzer = new THaAnalyzer;

    gHaApps->Add(harm);

    THaEvent* event = new THaEvent;

    TString prefix = gSystem->Getenv("DATA_DIR");

    bool segmentexists = true;
    int segment = firstsegment;
    int lastsegment = firstsegment;

    TClonesArray *filelist = new TClonesArray("THaRun",10);

    TDatime now = TDatime();

    vector<TString> pathlist;

    pathlist.push_back(prefix);

    if( prefix != "/adaqeb2/data1" ) pathlist.push_back( "/adaqeb2/data1" );

    if( prefix != "/data/raw" ) pathlist.push_back( "/data/raw" );

    if( prefix != "/adaq1/data1" ) pathlist.push_back( "/adaq1/data1" );

    if( prefix != "/cache/halla/sbs/raw" ) pathlist.push_back( "/cache/halla/sbs/raw" );

    if( prefix != "/cache/halla/sbs/GEnRP/raw" ) pathlist.push_back("/cache/halla/sbs/GEnRP/raw");

    if( prefix != "/cache/halla/sbs/GEp/raw" ) pathlist.push_back("/cache/halla/sbs/GEp/raw");

    for( int i=0; i<pathlist.size(); i++ )
    {
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

        } 
        
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

    if( nevents > 0 )
    { 
        outfilebase.Form( "%s_replayed_allGEMs_%d_seg%d_%d_firstevent%d_nevent%d", fname_prefix, runnum, firstsegment, lastsegment, firstevent, nevents );
        outfilename.Form( "%s/%s.root", prefix.Data(), outfilebase.Data());
    }
    else {
        outfilebase.Form( "%s_fullreplay_allGEMs_%d_seg%d_%d", fname_prefix, runnum, firstsegment, lastsegment );
        outfilename.Form( "%s/%s.root", prefix.Data(), outfilebase.Data());
    }

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

    TString odef_filename = "replay_gep_allGEM.odef";
    
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