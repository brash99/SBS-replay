#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <iostream>
#include <fstream>

void MakeOdef_and_cfg_files( const char *configfilename ){ 
  ifstream configfile(configfilename);

  TString detname="bb.gem";
  TString detname_nodots;

  int nmodules=12;
  int nlayers=5;
  
  double maxlayersizeX = 2.2; //meters
  double maxlayersizeY = 0.75; //meters 

  double maxTrackXp = 1.2;
  double maxTrackYp = 0.4; 
  
  vector<int> mod_nstripu,mod_nstripv;
  vector<double> mod_xsize, mod_ysize;
  vector<double> mod_x0, mod_y0, mod_z0;
  vector<TString> modname;
  vector<TString> modname_nodots;
  
  TString currentline;

  while( currentline.ReadLine( configfile ) ){
    if( !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize( " " );

      if( tokens->GetEntries() >= 2 ){
	TString skey = ( (TObjString*) (*tokens)[0] )->GetString();

	if( skey == "name" ){
	  detname = ( (TObjString*) (*tokens)[1] )->GetString();
	}

	if( skey == "nmodules" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();

	  nmodules = stemp.Atoi();
	  mod_nstripu.resize( nmodules );
	  mod_nstripv.resize( nmodules );
	  mod_xsize.resize( nmodules );
	  mod_ysize.resize( nmodules );
	  mod_x0.resize( nmodules );
	  mod_y0.resize( nmodules );
	  mod_z0.resize( nmodules );

	  for( int i=0; i<nmodules; i++ ){
	    TString nametemp = detname;

	    modname.push_back( nametemp.Format( "%s.m%d", nametemp.Data(), i ) );
	    
	    nametemp.ReplaceAll(".","_");

	    modname_nodots.push_back( nametemp.Format("%s_m%d", nametemp.Data(), i ) );
	    
	    
	  }
	}

	if( skey == "nlayers" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();

	  nlayers = stemp.Atoi();
	}

	if( skey == "layersizeX" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();

	  maxlayersizeX = stemp.Atof();
	}

	if( skey == "layersizeY" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();

	  maxlayersizeY = stemp.Atof();
	}

	if( skey == "maxtrackxp" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();

	  maxTrackXp = stemp.Atof();
	}

	if( skey == "maxtrackyp" ){
	  TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();

	  maxTrackYp = stemp.Atof();
	}
	
	if( skey == "mod_nstripu" && tokens->GetEntries() >= nmodules+1 ){
	  for( int i=1; i<=nmodules; i++ ){
	    TString sn = ( (TObjString*) (*tokens)[i] )->GetString();

	    int mod = i-1;
	    
	    mod_nstripu[mod] = sn.Atoi();
	    
	  }
	}

	if( skey == "mod_nstripv" && tokens->GetEntries() >= nmodules + 1){
	  for( int i=1; i<=nmodules; i++ ){
	    TString sn = ( (TObjString*) (*tokens)[i] )->GetString();

	    int mod = i-1;
	    
	    mod_nstripv[mod] = sn.Atoi();
	    
	  }
	}

	if( skey == "mod_xsize" && tokens->GetEntries() >= nmodules + 1 ){
	  for( int i=1; i<=nmodules; i++ ){
	    TString sn = ( (TObjString*) (*tokens)[i] )->GetString();

	    int mod = i-1;
	    
	    mod_xsize[mod] = sn.Atof();
	    
	  }
	}

	if( skey == "mod_ysize" && tokens->GetEntries() >= nmodules + 1 ){
	  for( int i=1; i<=nmodules; i++ ){
	    TString sn = ( (TObjString*) (*tokens)[i] )->GetString();

	    int mod = i-1;
	    
	    mod_ysize[mod] = sn.Atof();
	    
	  }
	}

	if( skey == "mod_x0" && tokens->GetEntries() >= nmodules + 1 ){
	  for( int i=1; i<=nmodules; i++ ){
	    TString sn = ( (TObjString*) (*tokens)[i] )->GetString();

	    int mod = i-1;
	    
	    mod_x0[mod] = sn.Atof();
	    
	  }
	}


	if( skey == "mod_y0" && tokens->GetEntries() >= nmodules + 1 ){
	  for( int i=1; i<=nmodules; i++ ){
	    TString sn = ( (TObjString*) (*tokens)[i] )->GetString();

	    int mod = i-1;
	    
	    mod_y0[mod] = sn.Atof();
	    
	  }
	}

	if( skey == "mod_z0" && tokens->GetEntries() >= nmodules + 1 ){
	  for( int i=1; i<=nmodules; i++ ){
	    TString sn = ( (TObjString*) (*tokens)[i] )->GetString();

	    int mod = i-1;
	    
	    mod_z0[mod] = sn.Atof();
	    
	  }
	}	
      }
    } 
  }

  detname_nodots = detname;
  detname_nodots.ReplaceAll(".","_");
  
  TString odef_fname;
  odef_fname.Form( "replay_%s.odef", detname_nodots.Data() );

  ofstream odef_file(odef_fname.Data());

  TString odefline;

  vector<TString> treevarnames;
  treevarnames.push_back( "track.*" );
  treevarnames.push_back( "hit.*" );
  treevarnames.push_back( "nlayershit" );
  treevarnames.push_back( "nlayershitu" );
  treevarnames.push_back( "nlayershitv" );
  treevarnames.push_back( "nlayershituv" );
  treevarnames.push_back( "nstripsu_layer" );
  treevarnames.push_back( "nstripsv_layer" );
  treevarnames.push_back( "nclustu_layer" );
  treevarnames.push_back( "nclustv_layer" );
  treevarnames.push_back( "n2Dhit_layer" );
  treevarnames.push_back( "*.time.T0_by_APV" );
  treevarnames.push_back( "*.time.Tref_coarse" );
  treevarnames.push_back( "*.time.Tcoarse_by_APV" );
  treevarnames.push_back( "*.time.Tfine_by_APV" );
  treevarnames.push_back( "*.time.EventCount_by_APV" );
  treevarnames.push_back( "*.time.T_ns_by_APV" );

  odef_file << "#tree variables: " << endl;
  
  for( int i=0; i<treevarnames.size(); i++ ){
    if( treevarnames[i].Contains("*") ){ //assume a block of variables:
      odefline.Form( "block %s.%s", detname.Data(), treevarnames[i].Data() );
    } else {
      odefline.Form( "variable %s.%s", detname.Data(), treevarnames[i].Data() );
    }

    odef_file << odefline << endl;
  }

  odef_file << endl;
  //don't worry about the 
  // vector<TString> cutnames just yet, we can add those manually

  odef_file << "#Cuts: " << endl;
  TString trackcut;
  trackcut.Form( "cut singletrack %s.track.ntrack == 1", detname.Data() );
  
  odef_file << trackcut << endl;
  odef_file << trackcut.Format( "cut anytrack %s.track.ntrack > 0", detname.Data() ) << endl << endl;
  
  odef_file << "# Histogram definitions: " << endl;
  odef_file << "# layers fired, strip and cluster multiplicities" << endl;

  TString histdef;

  
  odef_file << histdef.Format( "th1d h%s_Nlayers_hit 'Number of layers with fired strips (u/x OR v/y)' %s.nlayershit %d %g %g",
			       detname_nodots.Data(), detname.Data(), nlayers+1, -0.5, nlayers+0.5 ) << endl;
  odef_file << histdef.Format( "th1d h%s_Nlayers_hitu 'Number of layers with fired U/X strips' %s.nlayershitu %d %g %g",
			       detname_nodots.Data(), detname.Data(), nlayers+1, -0.5, nlayers+0.5 ) << endl;

  odef_file << histdef.Format( "th1d h%s_Nlayers_hitv 'Number of layers with fired V/Y strips' %s.nlayershitv %d %g %g",
			       detname_nodots.Data(), detname.Data(), nlayers+1, -0.5, nlayers+0.5 ) << endl;
  odef_file << histdef.Format( "th1d h%s_Nlayers_hituv 'Number of layers with 2D hits' %s.nlayershituv %d %g %g",
			       detname_nodots.Data(), detname.Data(), nlayers+1, -0.5, nlayers+0.5 ) << endl << endl;


  
  odef_file << histdef.Format( "th2d h%s_NstripsU_layer 'Number of U/X strips fired by layer' [I] %s.nstripsu_layer[I] %d %g %g 1001 -0.5 1000.5",
			       detname_nodots.Data(), detname.Data(), nlayers, -0.5, nlayers-.5 ) << endl;
  odef_file << histdef.Format( "th2d h%s_NstripsV_layer 'Number of V/Y strips fired by layer' [I] %s.nstripsv_layer[I] %d %g %g 1001 -0.5 1000.5",
			       detname_nodots.Data(), detname.Data(), nlayers, -0.5, nlayers-0.5 ) << endl;


  odef_file << histdef.Format( "th2d h%s_NclustU_layer 'Number of U/X clusters by layer' [I] %s.nclustu_layer[I] %d %g %g 1001 -0.5 1000.5",
			       detname_nodots.Data(), detname.Data(), nlayers, -0.5, nlayers-0.5 ) << endl;
  odef_file << histdef.Format( "th2d h%s_NclustV_layer 'Number of V/Y clusters by layer' [I] %s.nclustv_layer[I] %d %g %g 1001 -0.5 1000.5",
			       detname_nodots.Data(), detname.Data(), nlayers, -0.5, nlayers-0.5 ) << endl;
  odef_file << histdef.Format( "th2d h%s_Nclust2D_layer 'Number of 2D hits by layer' [I] %s.n2Dhit_layer[I] %d %g %g 1001 -0.5 1000.5",
			       detname_nodots.Data(), detname.Data(), nlayers, -0.5, nlayers-0.5 ) << endl << endl;


  odef_file << histdef.Format( "th1d h%s_clustwidthU 'U/X cluster size in strips (clusters on tracks)' %s.hit.nstripu 10 0.5 10.5 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;
  odef_file << histdef.Format( "th1d h%s_clustwidthV 'V/Y cluster size in strips (clusters on tracks)' %s.hit.nstripv 10 0.5 10.5 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;
  odef_file << histdef.Format( "th2d h%s_clustwidthVvsU 'V/Y cluster size vs U/X cluster size' %s.hit.nstripu %s.hit.nstripv 10 0.5 10.5 10 0.5 10.5 singletrack",
			       detname_nodots.Data(), detname.Data(), detname.Data() ) << endl;


  odef_file << histdef.Format( "th2d h%s_clustwidthU_vs_module 'U/X cluster size by module' %s.hit.module %s.hit.nstripu %d -0.5 %g 10 0.5 10.5 singletrack",
			       detname_nodots.Data(), detname.Data(), detname.Data(), nmodules, nmodules-0.5 ) << endl;
  odef_file << histdef.Format( "th2d h%s_clustwidthV_vs_module 'V/Y cluster size by module' %s.hit.module %s.hit.nstripv %d -0.5 %g 10 0.5 10.5 singletrack",
			       detname_nodots.Data(), detname.Data(), detname.Data(), nmodules, nmodules-0.5 ) << endl;
  
  
  odef_file << histdef.Format( "th1d h%s_clust_Utime 'U/X cluster time (ns)' %s.hit.Utime 150 0 150 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;
  odef_file << histdef.Format( "th1d h%s_clust_Vtime 'V/Y cluster time (ns)' %s.hit.Vtime 150 0 150 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;
  odef_file << histdef.Format( "th2d h%s_clust_UVtime 'V/Y vs U/X cluster time (ns)' %s.hit.Utime %s.hit.Vtime 150 0 150 150 0 150 singletrack",
			       detname_nodots.Data(), detname.Data(), detname.Data() ) << endl;

  odef_file << histdef.Format( "th1d h%s_clust_deltat 't_{U/X} - t_{V/Y} (ns)' %s.hit.deltat 200 -50 50 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;
  odef_file << histdef.Format( "th2d h%s_clust_deltat_vs_tavg 't_{diff} vs t_{avg} (ns)' %s.hit.Tavg %s.hit.deltat 150 0 150 200 -50 50 singletrack",
			       detname_nodots.Data(), detname.Data(), detname.Data() ) << endl;
  
  odef_file << histdef.Format( "th1d h%s_iSampMaxU 'peak time sample, cluster-summed U/X strips' %s.hit.isampmaxUclust 6 -0.5 5.5",
			       detname_nodots.Data(), detname.Data() ) << endl;
  odef_file << histdef.Format( "th1d h%s_iSampMaxV 'peak time sample, cluster-summed V/Y strips' %s.hit.isampmaxVclust 6 -0.5 5.5",
			       detname_nodots.Data(), detname.Data() ) << endl;
  
  odef_file << histdef.Format( "th1d h%s_CorrCoeff_clust 'U/V (X/Y) correlation coeff., cluster-summed ADC samples' %s.hit.ccor_clust 202 -1.01 1.01 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;
  odef_file << histdef.Format( "th1d h%s_CorrCoeff_maxstrip 'U/V (X/Y) correlation coeff., max strips in cluster' %s.hit.ccor_strip 202 -1.01 1.01 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;

  odef_file << histdef.Format( "th1d h%s_ADCasym 'ADC asymmetry (ADCU-ADCV)/(ADCU+ADCV), hits on tracks' %s.hit.ADCasym 202 -1.01 1.01 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl << endl;

  odef_file << histdef.Format( "th1d h%s_ADCavg '(ADCU+ADCV)/2, hits on tracks' %s.hit.ADCavg 1500 0 30000 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;
  odef_file << histdef.Format( "th2d h%s_ADCasym_vs_ADCavg 'ADC asym. vs ADC avg' %s.hit.ADCavg %s.hit.ADCasym 250 0 25000 250 -1.01 1.01",
			       detname_nodots.Data(), detname.Data(), detname.Data() ) << endl;
  
  odef_file << histdef.Format( "th1d h%s_ADCU_clust 'ADC cluster sum (U/X strips)' %s.hit.ADCU 1500 0 30000 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;
  odef_file << histdef.Format( "th1d h%s_ADCV_clust 'ADC cluster sum (V/Y strips)' %s.hit.ADCV 1500 0 30000 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;
  odef_file << histdef.Format( "th2d h%s_ADCVvsU_clust 'ADC cluster sum, V/Y vs U/X' %s.hit.ADCU %s.hit.ADCV 200 0 30000 200 0 30000 singletrack",
			       detname_nodots.Data(), detname.Data(), detname.Data() ) << endl;

  odef_file << histdef.Format( "th2d h%s_ADCavg_vs_module 'Cluster ADC avg vs module' %s.hit.module %s.hit.ADCavg %d -0.5 %g 250 0 25000",
			       detname_nodots.Data(), detname.Data(), detname.Data(), nmodules, nmodules-0.5 ) << endl;
  
  odef_file << histdef.Format( "th1d h%s_ADCU_maxstrip 'Max strip ADC sum (U/X strips)' %s.hit.ADCmaxstripU 1500 0 15000 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;
  odef_file << histdef.Format( "th1d h%s_ADCV_maxstrip 'Max strip ADC sum (V/Y strips)' %s.hit.ADCmaxstripV 1500 0 15000 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;
  odef_file << histdef.Format( "th2d h%s_ADCVvsU_maxstrip 'Max strip ADC sum, V/Y vs U/X' %s.hit.ADCmaxstripU %s.hit.ADCmaxstripV 200 0 12000 200 0 12000 singletrack",
			       detname_nodots.Data(), detname.Data(), detname.Data() ) << endl;



  odef_file << histdef.Format( "th1d h%s_ADCU_maxsample 'Max ADC sample in max strip (U/X strips)' %s.hit.ADCmaxsampU 500 0 4000 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;
  odef_file << histdef.Format( "th1d h%s_ADCV_maxsample 'Max ADC sample in max strip (V/Y strips)' %s.hit.ADCmaxsampV 500 0 4000 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;
  odef_file << histdef.Format( "th2d h%s_ADCVvsU_maxsample 'Max ADC sample in max strip, V/Y vs U/X' %s.hit.ADCmaxsampU %s.hit.ADCmaxsampV 200 0 4000 200 0 4000 singletrack",
			       detname_nodots.Data(), detname.Data(), detname.Data() ) << endl;


  odef_file << histdef.Format( "th1d h%s_ADCU_maxclustsample 'Max cluster-summed time sample (U/X strips)' %s.hit.ADCmaxsampUclust 500 0 8000 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;
  odef_file << histdef.Format( "th1d h%s_ADCV_maxclustsample 'Max cluster-summed time sample (V/Y strips)' %s.hit.ADCmaxsampVclust 500 0 8000 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;
  odef_file << histdef.Format( "th2d h%s_ADCVvsU_maxclustsample 'Max cluster-summed time sample, V/Y vs U/X' %s.hit.ADCmaxsampUclust %s.hit.ADCmaxsampVclust 200 0 8000 200 0 8000 singletrack",
			       detname_nodots.Data(), detname.Data(), detname.Data() ) << endl << endl;
  
  
  //Now put together the track stuff:

  odef_file << histdef.Format( "th1d h%s_Ntracks 'Number of tracks found' %s.track.ntrack 11 -0.5 10.5",
			       detname_nodots.Data(), detname.Data() ) << endl;
  
  odef_file << histdef.Format( "th1d h%s_TrackNhits 'Number of hits on first/best track' %s.track.nhits[I] %d -0.5 %g singletrack",
			       detname_nodots.Data(), detname.Data(), nlayers+1, nlayers+0.5 ) << endl;
  odef_file << histdef.Format( "th1d h%s_TrackX 'Best track X(z=0), m' %s.track.x[I] 200 %g %g singletrack",
			       detname_nodots.Data(), detname.Data(), -maxlayersizeX/2.0, maxlayersizeX/2.0 ) << endl;
  odef_file << histdef.Format( "th1d h%s_TrackY 'Best track Y(z=0), m' %s.track.y[I] 200 %g %g singletrack",
			       detname_nodots.Data(), detname.Data(), -maxlayersizeY/2.0, maxlayersizeY/2.0 ) << endl;

  odef_file << histdef.Format( "th1d h%s_TrackXp 'Best track dx/dz' %s.track.xp[I] 200 %g %g singletrack",
			       detname_nodots.Data(), detname.Data(), -maxTrackXp, maxTrackXp ) << endl;
  odef_file << histdef.Format( "th1d h%s_TrackYp 'Best track dy/dz' %s.track.yp[I] 200 %g %g singletrack",
			       detname_nodots.Data(), detname.Data(), -maxTrackYp, maxTrackYp) << endl;
  
  odef_file << histdef.Format( "th2d h%s_TrackXY 'Best track x vs y, m' %s.track.y[I] %s.track.x[I] 150 %g %g 150 %g %g singletrack",
			       detname_nodots.Data(), detname.Data(), detname.Data(), -maxlayersizeY/2.0, maxlayersizeY/2.0, -maxlayersizeX/2.0, maxlayersizeX/2.0 ) << endl;

  odef_file << histdef.Format( "th2d h%s_TrackXpYp 'Best track dx/dz vs dy/dz' %s.track.yp[I] %s.track.xp[I] 150 %g %g 150 %g %g singletrack",
			       detname_nodots.Data(), detname.Data(), detname.Data(), -maxTrackYp, maxTrackYp, -maxTrackXp, maxTrackXp ) << endl;

  odef_file << histdef.Format( "th1d h%s_TrackChi2NDF 'Best track chi2/ndf' %s.track.chi2ndf[I] 200 0 200 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl << endl;


  //what is next: tracking residuals inclusive and exclusive, all hits:
  odef_file << histdef.Format( "th1d h%s_residu_allhits 'Track u/x incl. residuals, all hits in all layers' %s.hit.residu 500 -0.0025 0.0025 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;
  odef_file << histdef.Format( "th1d h%s_residv_allhits 'Track v/y incl. residuals, all hits in all layers' %s.hit.residv 500 -0.0025 0.0025 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;

  odef_file << histdef.Format( "th1d h%s_eresidu_allhits 'Track u/x excl. residuals, all hits in all layers' %s.hit.eresidu 500 -0.0025 0.0025 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;
  odef_file << histdef.Format( "th1d h%s_eresidv_allhits 'Track v/y excl. residuals, all hits in all layers' %s.hit.eresidv 500 -0.0025 0.0025 singletrack",
			       detname_nodots.Data(), detname.Data() ) << endl;

  //inclusive tracking residuals by layer and module:
  odef_file << histdef.Format( "th2d h%s_residu_vs_layer 'Track u/x incl. residuals vs layer' %s.hit.layer %s.hit.residu %d -0.5 %g 250 -0.0025 0.0025 singletrack",
			       detname_nodots.Data(), detname.Data(), detname.Data(), nlayers, nlayers - 0.5 ) << endl;
  odef_file << histdef.Format( "th2d h%s_residv_vs_layer 'Track v/y incl. residuals vs layer' %s.hit.layer %s.hit.residv %d -0.5 %g 250 -0.0025 0.0025 singletrack",
			       detname_nodots.Data(), detname.Data(), detname.Data(), nlayers, nlayers - 0.5 ) << endl;

  odef_file << histdef.Format( "th2d h%s_residu_vs_module 'Track u/x incl. residuals vs module' %s.hit.module %s.hit.residu %d -0.5 %g 250 -0.0025 0.0025 singletrack",
			       detname_nodots.Data(), detname.Data(), detname.Data(), nmodules, nmodules - 0.5 ) << endl;
  odef_file << histdef.Format( "th2d h%s_residv_vs_module 'Track v/y incl. residuals vs module' %s.hit.module %s.hit.residv %d -0.5 %g 250 -0.0025 0.0025 singletrack",
			       detname_nodots.Data(), detname.Data(), detname.Data(), nmodules, nmodules - 0.5 ) << endl;


  //exclusive tracking residuals by layer and module:
  odef_file << histdef.Format( "th2d h%s_eresidu_vs_layer 'Track u/x excl. residuals vs layer' %s.hit.layer %s.hit.eresidu %d -0.5 %g 250 -0.0025 0.0025 singletrack",
			       detname_nodots.Data(), detname.Data(), detname.Data(), nlayers, nlayers - 0.5 ) << endl;
  odef_file << histdef.Format( "th2d h%s_eresidv_vs_layer 'Track v/y excl. residuals vs layer' %s.hit.layer %s.hit.eresidv %d -0.5 %g 250 -0.0025 0.0025 singletrack",
			       detname_nodots.Data(), detname.Data(), detname.Data(), nlayers, nlayers - 0.5 ) << endl;

  odef_file << histdef.Format( "th2d h%s_eresidu_vs_module 'Track u/x excl. residuals vs module' %s.hit.module %s.hit.eresidu %d -0.5 %g 250 -0.0025 0.0025 singletrack",
			       detname_nodots.Data(), detname.Data(), detname.Data(), nmodules, nmodules - 0.5 ) << endl;
  odef_file << histdef.Format( "th2d h%s_eresidv_vs_module 'Track v/y excl. residuals vs module' %s.hit.module %s.hit.eresidv %d -0.5 %g 250 -0.0025 0.0025 singletrack",
			       detname_nodots.Data(), detname.Data(), detname.Data(), nmodules, nmodules - 0.5 ) << endl;

 
  //Next up: hit maps by layer/module and strip heat maps/ADC distributions by module:
  odef_file << histdef.Format( "th2d h%s_hit_xglobal_vs_layer 'Hit x global (m) by layer' %s.hit.layer %s.hit.xglobal %d -0.5 %g 250 %g %g",
			       detname_nodots.Data(), detname.Data(), detname.Data(), nlayers, nlayers-0.5, -maxlayersizeX/2.0, maxlayersizeX/2.0 ) << endl;
  odef_file << histdef.Format( "th2d h%s_hit_yglobal_vs_layer 'Hit y global (m) by layer' %s.hit.layer %s.hit.yglobal %d -0.5 %g 250 %g %g",
			       detname_nodots.Data(), detname.Data(), detname.Data(), nlayers, nlayers-0.5, -maxlayersizeY/2.0, maxlayersizeY/2.0 ) << endl << endl;


  odef_file << "# Module-specific histograms: " << endl;
    
  
  //All module-specific histograms defined in one loop:
  for( int i=0; i<nmodules; i++ ){
    odef_file << "# Module " << i << " strip histograms, all strips: " << endl;
    //Strip counts, require track cut:
    odef_file << histdef.Format( "th1d h%s_nstripstot 'Num. fired strips (U/X) + (V/Y), Module %d' %s.strip.nstripsfired %d -0.5 %g %s.ontrack",
				 modname_nodots[i].Data(), i, modname[i].Data(), mod_nstripu[i], mod_nstripu[i]-0.5, modname[i].Data() ) << endl;

    odef_file << histdef.Format( "th1d h%s_nclustU 'num. U clusters, Module %d' %s.clust.nclustu 251 -0.5 250.5 %s.ontrack",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data() ) << endl;

    odef_file << histdef.Format( "th1d h%s_nclustV 'num. V clusters, Module %d' %s.clust.nclustv 251 -0.5 250.5 %s.ontrack",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data() ) << endl;
    
    //No cuts for strip heat map:
    odef_file << histdef.Format( "th1d h%s_stripU_all 'U/X strip index, Module %d (all strips)' %s.strip.istrip[I] %d -0.5 %g %s.strip.IsU[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), mod_nstripu[i], mod_nstripu[i]-0.5, modname[i].Data() ) << endl;
    //No cuts for strip heat map:
    odef_file << histdef.Format( "th1d h%s_stripV_all 'V/Y strip index, Module %d (all strips)' %s.strip.istrip[I] %d -0.5 %g %s.strip.IsV[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), mod_nstripv[i], mod_nstripv[i]-0.5, modname[i].Data() ) << endl;

    

    //1D histograms for strip T values and ADC values:
    odef_file << histdef.Format( "th1d h%s_ADCsumU_all 'Strip ADC sum, U/X strips, Module %d' %s.strip.ADCsum[I] 1000 0 15000 %s.strip.IsU[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data() ) << endl;
    odef_file << histdef.Format( "th1d h%s_ADCsumV_all 'Strip ADC sum, V/Y strips, Module %d' %s.strip.ADCsum[I] 1000 0 15000 %s.strip.IsV[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data() ) << endl;
    odef_file << histdef.Format( "th1d h%s_ADCmaxU_all 'Strip ADC max, U/X strips, Module %d' %s.strip.ADCmax[I] 500 0 4000 %s.strip.IsU[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data() ) << endl;
    odef_file << histdef.Format( "th1d h%s_ADCmaxV_all 'Strip ADC max, V/Y strips, Module %d' %s.strip.ADCmax[I] 500 0 4000 %s.strip.IsV[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data() ) << endl;
    odef_file << histdef.Format( "th1d h%s_TimeU_all 'Strip Time, U/X strips, Module %d' %s.strip.Tmean[I] 150 0 150 %s.strip.IsU[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data() ) << endl;
    odef_file << histdef.Format( "th1d h%s_TimeV_all 'Strip Time, V/Y strips, Module %d' %s.strip.Tmean[I] 150 0 150 %s.strip.IsV[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data() ) << endl;
    odef_file << histdef.Format( "th1d h%s_iSampMaxU_all 'Peak time sample, U/X strips, Module %d' %s.strip.isampmax[I] 6 -0.5 5.5 %s.strip.IsU[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data() ) << endl;
    odef_file << histdef.Format( "th1d h%s_iSampMaxV_all 'Peak time sample, V/Y strips, Module %d' %s.strip.isampmax[I] 6 -0.5 5.5 %s.strip.IsV[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data() ) << endl;


    odef_file << endl << "# Module " << i << " strip histograms, strips on tracks: " << endl;
    //Strips on tracks:
    odef_file << histdef.Format( "th1d h%s_stripU_good 'U/X strip index, Module %d (strips on tracks)' %s.strip.istrip[I] %d -0.5 %g %s.strip.ontrackU[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), mod_nstripu[i], mod_nstripu[i]-0.5, modname[i].Data() ) << endl;
    
    odef_file << histdef.Format( "th1d h%s_stripV_good 'V/Y strip index, Module %d (strips on tracks)' %s.strip.istrip[I] %d -0.5 %g %s.strip.ontrackV[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), mod_nstripv[i], mod_nstripv[i]-0.5, modname[i].Data() ) << endl;

    //1D histograms for strip T values and ADC values:
    odef_file << histdef.Format( "th1d h%s_ADCsumU_good 'Strip ADC sum, U/X strips, Module %d' %s.strip.ADCsum[I] 1000 0 15000 %s.strip.ontrackU[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data() ) << endl;
    odef_file << histdef.Format( "th1d h%s_ADCsumV_good 'Strip ADC sum, V/Y strips, Module %d' %s.strip.ADCsum[I] 1000 0 15000 %s.strip.ontrackV[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data() ) << endl;
    odef_file << histdef.Format( "th1d h%s_ADCmaxU_good 'Strip ADC max, U/X strips, Module %d' %s.strip.ADCmax[I] 500 0 4000 %s.strip.ontrackU[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data() ) << endl;
    odef_file << histdef.Format( "th1d h%s_ADCmaxV_good 'Strip ADC max, V/Y strips, Module %d' %s.strip.ADCmax[I] 500 0 4000 %s.strip.ontrackV[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data() ) << endl;
    odef_file << histdef.Format( "th1d h%s_TimeU_good 'Strip Time, U/X strips, Module %d' %s.strip.Tmean[I] 150 0 150 %s.strip.ontrackU[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data() ) << endl;
    odef_file << histdef.Format( "th1d h%s_TimeV_good 'Strip Time, V/Y strips, Module %d' %s.strip.Tmean[I] 150 0 150 %s.strip.ontrackV[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data() ) << endl;
    odef_file << histdef.Format( "th1d h%s_iSampMaxU_good 'Peak time sample, U/X strips, Module %d' %s.strip.isampmax[I] 6 -0.5 5.5 %s.strip.ontrackU[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data() ) << endl;
    odef_file << histdef.Format( "th1d h%s_iSampMaxV_good 'Peak time sample, V/Y strips, Module %d' %s.strip.isampmax[I] 6 -0.5 5.5 %s.strip.ontrackV[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data() ) << endl;

    //ADC sum dist. vs strip:
    odef_file << histdef.Format( "th2d h%s_ADCsum_vs_Ustrip_good 'U strip ADC sum vs strip, module %d' %s.strip.istrip[I] %s.strip.ADCsum[I] %d -0.5 %g 250 0 12000 %s.strip.ontrackU[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data(), mod_nstripu[i], mod_nstripu[i]-0.5, modname[i].Data() ) << endl;
    odef_file << histdef.Format( "th2d h%s_ADCsum_vs_Vstrip_good 'V strip ADC sum vs strip, module %d' %s.strip.istrip[I] %s.strip.ADCsum[I] %d -0.5 %g 250 0 12000 %s.strip.ontrackV[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data(), mod_nstripv[i], mod_nstripv[i]-0.5, modname[i].Data() ) << endl;

    //ADC max dist. vs strip:
    odef_file << histdef.Format( "th2d h%s_ADCmax_vs_Ustrip_good 'U strip ADC max vs strip, module %d' %s.strip.istrip[I] %s.strip.ADCmax[I] %d -0.5 %g 250 0 4000 %s.strip.ontrackU[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data(), mod_nstripu[i], mod_nstripu[i]-0.5, modname[i].Data() ) << endl;
    odef_file << histdef.Format( "th2d h%s_ADCmax_vs_Vstrip_good 'V strip ADC max vs strip, module %d' %s.strip.istrip[I] %s.strip.ADCmax[I] %d -0.5 %g 250 0 4000 %s.strip.ontrackV[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data(), mod_nstripv[i], mod_nstripv[i]-0.5, modname[i].Data() ) << endl;
    
    
    odef_file << endl << "# Module " << i << " hit maps " << endl;

    odef_file << histdef.Format( "th1d h%s_hit_xlocal 'Hit local X position (m), Module %d' %s.hit.hitx[I] 250 %g %g %s.hit.ontrack[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), -mod_xsize[i]/2.0, mod_xsize[i]/2.0, modname[i].Data() ) << endl;
    odef_file << histdef.Format( "th1d h%s_hit_ylocal 'Hit local Y position (m), Module %d' %s.hit.hity[I] 250 %g %g %s.hit.ontrack[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), -mod_ysize[i]/2.0, mod_ysize[i]/2.0, modname[i].Data() ) << endl;
    odef_file << histdef.Format( "th2d h%s_hit_xylocal 'Hit local X vs Y (m), Module %d' %s.hit.hity[I] %s.hit.hitx[I] 200 %g %g 200 %g %g %s.hit.ontrack[I]",
				 modname_nodots[i].Data(), i, modname[i].Data(), modname[i].Data(), -mod_ysize[i]/2.0, mod_ysize[i]/2.0, -mod_xsize[i]/2.0, mod_xsize[i]/2.0,
				 modname[i].Data() ) << endl;

    
    odef_file << endl;
  }
  

  
  
}
