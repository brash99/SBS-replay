//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct  9 12:05:13 2024 by ROOT version 6.30/04
// from TTree T/Hall A Analyzer Output DST
// found on file: /volatile/halla/sbs/puckett/GENRP_REPLAYS/test/rootfiles/e1209019_replayed_13486_stream0_seg0_0_firstevent0_nevent100.root
//////////////////////////////////////////////////////////

#ifndef gmn_tree_allstripoutputs_h
#define gmn_tree_allstripoutputs_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "THaEvent.h"
#include "THaEvent.h"

class gmn_tree_allstripoutputs {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        ax_new;
   Double_t        bx_new;
   Double_t        ay_new;
   Double_t        by_new;
   Double_t        targx;
   Double_t        targy;
   Int_t           Ndata_bb_bbtrig_a_amp_p;
   Double_t        bb_bbtrig_a_amp_p[25];   //[Ndata.bb.bbtrig.a_amp_p]
   Int_t           Ndata_bb_bbtrig_a_time;
   Double_t        bb_bbtrig_a_time[25];   //[Ndata.bb.bbtrig.a_time]
   Int_t           Ndata_bb_bbtrig_adcelemID;
   Double_t        bb_bbtrig_adcelemID[25];   //[Ndata.bb.bbtrig.adcelemID]
   Int_t           Ndata_bb_bbtrig_tdc;
   Double_t        bb_bbtrig_tdc[3];   //[Ndata.bb.bbtrig.tdc]
   Int_t           Ndata_bb_bbtrig_tdcelemID;
   Double_t        bb_bbtrig_tdcelemID[3];   //[Ndata.bb.bbtrig.tdcelemID]
   Int_t           Ndata_bb_eps_over_etot;
   Double_t        bb_eps_over_etot[1];   //[Ndata.bb.eps_over_etot]
   Int_t           Ndata_bb_etot_over_p;
   Double_t        bb_etot_over_p[1];   //[Ndata.bb.etot_over_p]
   Int_t           Ndata_bb_gem_hit_ADCU;
   Double_t        bb_gem_hit_ADCU[5];   //[Ndata.bb.gem.hit.ADCU]
   Int_t           Ndata_bb_gem_hit_ADCU_deconv;
   Double_t        bb_gem_hit_ADCU_deconv[5];   //[Ndata.bb.gem.hit.ADCU_deconv]
   Int_t           Ndata_bb_gem_hit_ADCV;
   Double_t        bb_gem_hit_ADCV[5];   //[Ndata.bb.gem.hit.ADCV]
   Int_t           Ndata_bb_gem_hit_ADCV_deconv;
   Double_t        bb_gem_hit_ADCV_deconv[5];   //[Ndata.bb.gem.hit.ADCV_deconv]
   Int_t           Ndata_bb_gem_hit_ADCasym;
   Double_t        bb_gem_hit_ADCasym[5];   //[Ndata.bb.gem.hit.ADCasym]
   Int_t           Ndata_bb_gem_hit_ADCasym_deconv;
   Double_t        bb_gem_hit_ADCasym_deconv[5];   //[Ndata.bb.gem.hit.ADCasym_deconv]
   Int_t           Ndata_bb_gem_hit_ADCavg;
   Double_t        bb_gem_hit_ADCavg[5];   //[Ndata.bb.gem.hit.ADCavg]
   Int_t           Ndata_bb_gem_hit_ADCavg_deconv;
   Double_t        bb_gem_hit_ADCavg_deconv[5];   //[Ndata.bb.gem.hit.ADCavg_deconv]
   Int_t           Ndata_bb_gem_hit_ADCfrac0_Umax;
   Double_t        bb_gem_hit_ADCfrac0_Umax[5];   //[Ndata.bb.gem.hit.ADCfrac0_Umax]
   Int_t           Ndata_bb_gem_hit_ADCfrac0_Vmax;
   Double_t        bb_gem_hit_ADCfrac0_Vmax[5];   //[Ndata.bb.gem.hit.ADCfrac0_Vmax]
   Int_t           Ndata_bb_gem_hit_ADCfrac1_Umax;
   Double_t        bb_gem_hit_ADCfrac1_Umax[5];   //[Ndata.bb.gem.hit.ADCfrac1_Umax]
   Int_t           Ndata_bb_gem_hit_ADCfrac1_Vmax;
   Double_t        bb_gem_hit_ADCfrac1_Vmax[5];   //[Ndata.bb.gem.hit.ADCfrac1_Vmax]
   Int_t           Ndata_bb_gem_hit_ADCfrac2_Umax;
   Double_t        bb_gem_hit_ADCfrac2_Umax[5];   //[Ndata.bb.gem.hit.ADCfrac2_Umax]
   Int_t           Ndata_bb_gem_hit_ADCfrac2_Vmax;
   Double_t        bb_gem_hit_ADCfrac2_Vmax[5];   //[Ndata.bb.gem.hit.ADCfrac2_Vmax]
   Int_t           Ndata_bb_gem_hit_ADCfrac3_Umax;
   Double_t        bb_gem_hit_ADCfrac3_Umax[5];   //[Ndata.bb.gem.hit.ADCfrac3_Umax]
   Int_t           Ndata_bb_gem_hit_ADCfrac3_Vmax;
   Double_t        bb_gem_hit_ADCfrac3_Vmax[5];   //[Ndata.bb.gem.hit.ADCfrac3_Vmax]
   Int_t           Ndata_bb_gem_hit_ADCfrac4_Umax;
   Double_t        bb_gem_hit_ADCfrac4_Umax[5];   //[Ndata.bb.gem.hit.ADCfrac4_Umax]
   Int_t           Ndata_bb_gem_hit_ADCfrac4_Vmax;
   Double_t        bb_gem_hit_ADCfrac4_Vmax[5];   //[Ndata.bb.gem.hit.ADCfrac4_Vmax]
   Int_t           Ndata_bb_gem_hit_ADCfrac5_Umax;
   Double_t        bb_gem_hit_ADCfrac5_Umax[5];   //[Ndata.bb.gem.hit.ADCfrac5_Umax]
   Int_t           Ndata_bb_gem_hit_ADCfrac5_Vmax;
   Double_t        bb_gem_hit_ADCfrac5_Vmax[5];   //[Ndata.bb.gem.hit.ADCfrac5_Vmax]
   Int_t           Ndata_bb_gem_hit_ADCmaxcomboUclust_deconv;
   Double_t        bb_gem_hit_ADCmaxcomboUclust_deconv[5];   //[Ndata.bb.gem.hit.ADCmaxcomboUclust_deconv]
   Int_t           Ndata_bb_gem_hit_ADCmaxcomboVclust_deconv;
   Double_t        bb_gem_hit_ADCmaxcomboVclust_deconv[5];   //[Ndata.bb.gem.hit.ADCmaxcomboVclust_deconv]
   Int_t           Ndata_bb_gem_hit_ADCmaxsampU;
   Double_t        bb_gem_hit_ADCmaxsampU[5];   //[Ndata.bb.gem.hit.ADCmaxsampU]
   Int_t           Ndata_bb_gem_hit_ADCmaxsampUclust;
   Double_t        bb_gem_hit_ADCmaxsampUclust[5];   //[Ndata.bb.gem.hit.ADCmaxsampUclust]
   Int_t           Ndata_bb_gem_hit_ADCmaxsampUclust_deconv;
   Double_t        bb_gem_hit_ADCmaxsampUclust_deconv[5];   //[Ndata.bb.gem.hit.ADCmaxsampUclust_deconv]
   Int_t           Ndata_bb_gem_hit_ADCmaxsampV;
   Double_t        bb_gem_hit_ADCmaxsampV[5];   //[Ndata.bb.gem.hit.ADCmaxsampV]
   Int_t           Ndata_bb_gem_hit_ADCmaxsampVclust;
   Double_t        bb_gem_hit_ADCmaxsampVclust[5];   //[Ndata.bb.gem.hit.ADCmaxsampVclust]
   Int_t           Ndata_bb_gem_hit_ADCmaxsampVclust_deconv;
   Double_t        bb_gem_hit_ADCmaxsampVclust_deconv[5];   //[Ndata.bb.gem.hit.ADCmaxsampVclust_deconv]
   Int_t           Ndata_bb_gem_hit_ADCmaxstripU;
   Double_t        bb_gem_hit_ADCmaxstripU[5];   //[Ndata.bb.gem.hit.ADCmaxstripU]
   Int_t           Ndata_bb_gem_hit_ADCmaxstripV;
   Double_t        bb_gem_hit_ADCmaxstripV[5];   //[Ndata.bb.gem.hit.ADCmaxstripV]
   Int_t           Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_U;
   Double_t        bb_gem_hit_BUILD_ALL_SAMPLES_U[5];   //[Ndata.bb.gem.hit.BUILD_ALL_SAMPLES_U]
   Int_t           Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_V;
   Double_t        bb_gem_hit_BUILD_ALL_SAMPLES_V[5];   //[Ndata.bb.gem.hit.BUILD_ALL_SAMPLES_V]
   Int_t           Ndata_bb_gem_hit_CM_GOOD_U;
   Double_t        bb_gem_hit_CM_GOOD_U[5];   //[Ndata.bb.gem.hit.CM_GOOD_U]
   Int_t           Ndata_bb_gem_hit_CM_GOOD_V;
   Double_t        bb_gem_hit_CM_GOOD_V[5];   //[Ndata.bb.gem.hit.CM_GOOD_V]
   Int_t           Ndata_bb_gem_hit_DeconvADC0_Umax;
   Double_t        bb_gem_hit_DeconvADC0_Umax[5];   //[Ndata.bb.gem.hit.DeconvADC0_Umax]
   Int_t           Ndata_bb_gem_hit_DeconvADC0_Vmax;
   Double_t        bb_gem_hit_DeconvADC0_Vmax[5];   //[Ndata.bb.gem.hit.DeconvADC0_Vmax]
   Int_t           Ndata_bb_gem_hit_DeconvADC1_Umax;
   Double_t        bb_gem_hit_DeconvADC1_Umax[5];   //[Ndata.bb.gem.hit.DeconvADC1_Umax]
   Int_t           Ndata_bb_gem_hit_DeconvADC1_Vmax;
   Double_t        bb_gem_hit_DeconvADC1_Vmax[5];   //[Ndata.bb.gem.hit.DeconvADC1_Vmax]
   Int_t           Ndata_bb_gem_hit_DeconvADC2_Umax;
   Double_t        bb_gem_hit_DeconvADC2_Umax[5];   //[Ndata.bb.gem.hit.DeconvADC2_Umax]
   Int_t           Ndata_bb_gem_hit_DeconvADC2_Vmax;
   Double_t        bb_gem_hit_DeconvADC2_Vmax[5];   //[Ndata.bb.gem.hit.DeconvADC2_Vmax]
   Int_t           Ndata_bb_gem_hit_DeconvADC3_Umax;
   Double_t        bb_gem_hit_DeconvADC3_Umax[5];   //[Ndata.bb.gem.hit.DeconvADC3_Umax]
   Int_t           Ndata_bb_gem_hit_DeconvADC3_Vmax;
   Double_t        bb_gem_hit_DeconvADC3_Vmax[5];   //[Ndata.bb.gem.hit.DeconvADC3_Vmax]
   Int_t           Ndata_bb_gem_hit_DeconvADC4_Umax;
   Double_t        bb_gem_hit_DeconvADC4_Umax[5];   //[Ndata.bb.gem.hit.DeconvADC4_Umax]
   Int_t           Ndata_bb_gem_hit_DeconvADC4_Vmax;
   Double_t        bb_gem_hit_DeconvADC4_Vmax[5];   //[Ndata.bb.gem.hit.DeconvADC4_Vmax]
   Int_t           Ndata_bb_gem_hit_DeconvADC5_Umax;
   Double_t        bb_gem_hit_DeconvADC5_Umax[5];   //[Ndata.bb.gem.hit.DeconvADC5_Umax]
   Int_t           Ndata_bb_gem_hit_DeconvADC5_Vmax;
   Double_t        bb_gem_hit_DeconvADC5_Vmax[5];   //[Ndata.bb.gem.hit.DeconvADC5_Vmax]
   Int_t           Ndata_bb_gem_hit_DeconvADCmaxcomboU;
   Double_t        bb_gem_hit_DeconvADCmaxcomboU[5];   //[Ndata.bb.gem.hit.DeconvADCmaxcomboU]
   Int_t           Ndata_bb_gem_hit_DeconvADCmaxcomboV;
   Double_t        bb_gem_hit_DeconvADCmaxcomboV[5];   //[Ndata.bb.gem.hit.DeconvADCmaxcomboV]
   Int_t           Ndata_bb_gem_hit_DeconvADCmaxsampU;
   Double_t        bb_gem_hit_DeconvADCmaxsampU[5];   //[Ndata.bb.gem.hit.DeconvADCmaxsampU]
   Int_t           Ndata_bb_gem_hit_DeconvADCmaxsampV;
   Double_t        bb_gem_hit_DeconvADCmaxsampV[5];   //[Ndata.bb.gem.hit.DeconvADCmaxsampV]
   Int_t           Ndata_bb_gem_hit_DeconvADCmaxstripU;
   Double_t        bb_gem_hit_DeconvADCmaxstripU[5];   //[Ndata.bb.gem.hit.DeconvADCmaxstripU]
   Int_t           Ndata_bb_gem_hit_DeconvADCmaxstripV;
   Double_t        bb_gem_hit_DeconvADCmaxstripV[5];   //[Ndata.bb.gem.hit.DeconvADCmaxstripV]
   Int_t           Ndata_bb_gem_hit_ENABLE_CM_U;
   Double_t        bb_gem_hit_ENABLE_CM_U[5];   //[Ndata.bb.gem.hit.ENABLE_CM_U]
   Int_t           Ndata_bb_gem_hit_ENABLE_CM_V;
   Double_t        bb_gem_hit_ENABLE_CM_V[5];   //[Ndata.bb.gem.hit.ENABLE_CM_V]
   Int_t           Ndata_bb_gem_hit_TSchi2_Umax;
   Double_t        bb_gem_hit_TSchi2_Umax[5];   //[Ndata.bb.gem.hit.TSchi2_Umax]
   Int_t           Ndata_bb_gem_hit_TSchi2_Vmax;
   Double_t        bb_gem_hit_TSchi2_Vmax[5];   //[Ndata.bb.gem.hit.TSchi2_Vmax]
   Int_t           Ndata_bb_gem_hit_TSprob_Umax;
   Double_t        bb_gem_hit_TSprob_Umax[5];   //[Ndata.bb.gem.hit.TSprob_Umax]
   Int_t           Ndata_bb_gem_hit_TSprob_Vmax;
   Double_t        bb_gem_hit_TSprob_Vmax[5];   //[Ndata.bb.gem.hit.TSprob_Vmax]
   Int_t           Ndata_bb_gem_hit_Tavg;
   Double_t        bb_gem_hit_Tavg[5];   //[Ndata.bb.gem.hit.Tavg]
   Int_t           Ndata_bb_gem_hit_Tavg_corr;
   Double_t        bb_gem_hit_Tavg_corr[5];   //[Ndata.bb.gem.hit.Tavg_corr]
   Int_t           Ndata_bb_gem_hit_Tavg_deconv;
   Double_t        bb_gem_hit_Tavg_deconv[5];   //[Ndata.bb.gem.hit.Tavg_deconv]
   Int_t           Ndata_bb_gem_hit_Tavg_fit;
   Double_t        bb_gem_hit_Tavg_fit[5];   //[Ndata.bb.gem.hit.Tavg_fit]
   Int_t           Ndata_bb_gem_hit_Ugain;
   Double_t        bb_gem_hit_Ugain[5];   //[Ndata.bb.gem.hit.Ugain]
   Int_t           Ndata_bb_gem_hit_Utime;
   Double_t        bb_gem_hit_Utime[5];   //[Ndata.bb.gem.hit.Utime]
   Int_t           Ndata_bb_gem_hit_UtimeDeconv;
   Double_t        bb_gem_hit_UtimeDeconv[5];   //[Ndata.bb.gem.hit.UtimeDeconv]
   Int_t           Ndata_bb_gem_hit_UtimeFit;
   Double_t        bb_gem_hit_UtimeFit[5];   //[Ndata.bb.gem.hit.UtimeFit]
   Int_t           Ndata_bb_gem_hit_UtimeMaxStrip;
   Double_t        bb_gem_hit_UtimeMaxStrip[5];   //[Ndata.bb.gem.hit.UtimeMaxStrip]
   Int_t           Ndata_bb_gem_hit_UtimeMaxStripDeconv;
   Double_t        bb_gem_hit_UtimeMaxStripDeconv[5];   //[Ndata.bb.gem.hit.UtimeMaxStripDeconv]
   Int_t           Ndata_bb_gem_hit_UtimeMaxStripFit;
   Double_t        bb_gem_hit_UtimeMaxStripFit[5];   //[Ndata.bb.gem.hit.UtimeMaxStripFit]
   Int_t           Ndata_bb_gem_hit_Vgain;
   Double_t        bb_gem_hit_Vgain[5];   //[Ndata.bb.gem.hit.Vgain]
   Int_t           Ndata_bb_gem_hit_Vtime;
   Double_t        bb_gem_hit_Vtime[5];   //[Ndata.bb.gem.hit.Vtime]
   Int_t           Ndata_bb_gem_hit_VtimeDeconv;
   Double_t        bb_gem_hit_VtimeDeconv[5];   //[Ndata.bb.gem.hit.VtimeDeconv]
   Int_t           Ndata_bb_gem_hit_VtimeFit;
   Double_t        bb_gem_hit_VtimeFit[5];   //[Ndata.bb.gem.hit.VtimeFit]
   Int_t           Ndata_bb_gem_hit_VtimeMaxStrip;
   Double_t        bb_gem_hit_VtimeMaxStrip[5];   //[Ndata.bb.gem.hit.VtimeMaxStrip]
   Int_t           Ndata_bb_gem_hit_VtimeMaxStripDeconv;
   Double_t        bb_gem_hit_VtimeMaxStripDeconv[5];   //[Ndata.bb.gem.hit.VtimeMaxStripDeconv]
   Int_t           Ndata_bb_gem_hit_VtimeMaxStripFit;
   Double_t        bb_gem_hit_VtimeMaxStripFit[5];   //[Ndata.bb.gem.hit.VtimeMaxStripFit]
   Int_t           Ndata_bb_gem_hit_ccor_clust;
   Double_t        bb_gem_hit_ccor_clust[5];   //[Ndata.bb.gem.hit.ccor_clust]
   Int_t           Ndata_bb_gem_hit_ccor_clust_deconv;
   Double_t        bb_gem_hit_ccor_clust_deconv[5];   //[Ndata.bb.gem.hit.ccor_clust_deconv]
   Int_t           Ndata_bb_gem_hit_ccor_strip;
   Double_t        bb_gem_hit_ccor_strip[5];   //[Ndata.bb.gem.hit.ccor_strip]
   Int_t           Ndata_bb_gem_hit_ccor_strip_deconv;
   Double_t        bb_gem_hit_ccor_strip_deconv[5];   //[Ndata.bb.gem.hit.ccor_strip_deconv]
   Int_t           Ndata_bb_gem_hit_deltat;
   Double_t        bb_gem_hit_deltat[5];   //[Ndata.bb.gem.hit.deltat]
   Int_t           Ndata_bb_gem_hit_deltat_deconv;
   Double_t        bb_gem_hit_deltat_deconv[5];   //[Ndata.bb.gem.hit.deltat_deconv]
   Int_t           Ndata_bb_gem_hit_deltat_fit;
   Double_t        bb_gem_hit_deltat_fit[5];   //[Ndata.bb.gem.hit.deltat_fit]
   Int_t           Ndata_bb_gem_hit_eresidu;
   Double_t        bb_gem_hit_eresidu[5];   //[Ndata.bb.gem.hit.eresidu]
   Int_t           Ndata_bb_gem_hit_eresidv;
   Double_t        bb_gem_hit_eresidv[5];   //[Ndata.bb.gem.hit.eresidv]
   Int_t           Ndata_bb_gem_hit_icombomaxUclustDeconv;
   Double_t        bb_gem_hit_icombomaxUclustDeconv[5];   //[Ndata.bb.gem.hit.icombomaxUclustDeconv]
   Int_t           Ndata_bb_gem_hit_icombomaxUstripDeconv;
   Double_t        bb_gem_hit_icombomaxUstripDeconv[5];   //[Ndata.bb.gem.hit.icombomaxUstripDeconv]
   Int_t           Ndata_bb_gem_hit_icombomaxVclustDeconv;
   Double_t        bb_gem_hit_icombomaxVclustDeconv[5];   //[Ndata.bb.gem.hit.icombomaxVclustDeconv]
   Int_t           Ndata_bb_gem_hit_icombomaxVstripDeconv;
   Double_t        bb_gem_hit_icombomaxVstripDeconv[5];   //[Ndata.bb.gem.hit.icombomaxVstripDeconv]
   Int_t           Ndata_bb_gem_hit_isampmaxUclust;
   Double_t        bb_gem_hit_isampmaxUclust[5];   //[Ndata.bb.gem.hit.isampmaxUclust]
   Int_t           Ndata_bb_gem_hit_isampmaxUclustDeconv;
   Double_t        bb_gem_hit_isampmaxUclustDeconv[5];   //[Ndata.bb.gem.hit.isampmaxUclustDeconv]
   Int_t           Ndata_bb_gem_hit_isampmaxUstrip;
   Double_t        bb_gem_hit_isampmaxUstrip[5];   //[Ndata.bb.gem.hit.isampmaxUstrip]
   Int_t           Ndata_bb_gem_hit_isampmaxUstripDeconv;
   Double_t        bb_gem_hit_isampmaxUstripDeconv[5];   //[Ndata.bb.gem.hit.isampmaxUstripDeconv]
   Int_t           Ndata_bb_gem_hit_isampmaxVclust;
   Double_t        bb_gem_hit_isampmaxVclust[5];   //[Ndata.bb.gem.hit.isampmaxVclust]
   Int_t           Ndata_bb_gem_hit_isampmaxVclustDeconv;
   Double_t        bb_gem_hit_isampmaxVclustDeconv[5];   //[Ndata.bb.gem.hit.isampmaxVclustDeconv]
   Int_t           Ndata_bb_gem_hit_isampmaxVstrip;
   Double_t        bb_gem_hit_isampmaxVstrip[5];   //[Ndata.bb.gem.hit.isampmaxVstrip]
   Int_t           Ndata_bb_gem_hit_isampmaxVstripDeconv;
   Double_t        bb_gem_hit_isampmaxVstripDeconv[5];   //[Ndata.bb.gem.hit.isampmaxVstripDeconv]
   Int_t           Ndata_bb_gem_hit_layer;
   Double_t        bb_gem_hit_layer[5];   //[Ndata.bb.gem.hit.layer]
   Int_t           Ndata_bb_gem_hit_module;
   Double_t        bb_gem_hit_module[5];   //[Ndata.bb.gem.hit.module]
   Int_t           Ndata_bb_gem_hit_nstripu;
   Double_t        bb_gem_hit_nstripu[5];   //[Ndata.bb.gem.hit.nstripu]
   Int_t           Ndata_bb_gem_hit_nstripv;
   Double_t        bb_gem_hit_nstripv[5];   //[Ndata.bb.gem.hit.nstripv]
   Int_t           Ndata_bb_gem_hit_residu;
   Double_t        bb_gem_hit_residu[5];   //[Ndata.bb.gem.hit.residu]
   Int_t           Ndata_bb_gem_hit_residv;
   Double_t        bb_gem_hit_residv[5];   //[Ndata.bb.gem.hit.residv]
   Int_t           Ndata_bb_gem_hit_trackindex;
   Double_t        bb_gem_hit_trackindex[5];   //[Ndata.bb.gem.hit.trackindex]
   Int_t           Ndata_bb_gem_hit_u;
   Double_t        bb_gem_hit_u[5];   //[Ndata.bb.gem.hit.u]
   Int_t           Ndata_bb_gem_hit_umoment;
   Double_t        bb_gem_hit_umoment[5];   //[Ndata.bb.gem.hit.umoment]
   Int_t           Ndata_bb_gem_hit_usigma;
   Double_t        bb_gem_hit_usigma[5];   //[Ndata.bb.gem.hit.usigma]
   Int_t           Ndata_bb_gem_hit_ustriphi;
   Double_t        bb_gem_hit_ustriphi[5];   //[Ndata.bb.gem.hit.ustriphi]
   Int_t           Ndata_bb_gem_hit_ustriplo;
   Double_t        bb_gem_hit_ustriplo[5];   //[Ndata.bb.gem.hit.ustriplo]
   Int_t           Ndata_bb_gem_hit_ustripmax;
   Double_t        bb_gem_hit_ustripmax[5];   //[Ndata.bb.gem.hit.ustripmax]
   Int_t           Ndata_bb_gem_hit_v;
   Double_t        bb_gem_hit_v[5];   //[Ndata.bb.gem.hit.v]
   Int_t           Ndata_bb_gem_hit_vmoment;
   Double_t        bb_gem_hit_vmoment[5];   //[Ndata.bb.gem.hit.vmoment]
   Int_t           Ndata_bb_gem_hit_vsigma;
   Double_t        bb_gem_hit_vsigma[5];   //[Ndata.bb.gem.hit.vsigma]
   Int_t           Ndata_bb_gem_hit_vstriphi;
   Double_t        bb_gem_hit_vstriphi[5];   //[Ndata.bb.gem.hit.vstriphi]
   Int_t           Ndata_bb_gem_hit_vstriplo;
   Double_t        bb_gem_hit_vstriplo[5];   //[Ndata.bb.gem.hit.vstriplo]
   Int_t           Ndata_bb_gem_hit_vstripmax;
   Double_t        bb_gem_hit_vstripmax[5];   //[Ndata.bb.gem.hit.vstripmax]
   Int_t           Ndata_bb_gem_hit_xglobal;
   Double_t        bb_gem_hit_xglobal[5];   //[Ndata.bb.gem.hit.xglobal]
   Int_t           Ndata_bb_gem_hit_xlocal;
   Double_t        bb_gem_hit_xlocal[5];   //[Ndata.bb.gem.hit.xlocal]
   Int_t           Ndata_bb_gem_hit_yglobal;
   Double_t        bb_gem_hit_yglobal[5];   //[Ndata.bb.gem.hit.yglobal]
   Int_t           Ndata_bb_gem_hit_ylocal;
   Double_t        bb_gem_hit_ylocal[5];   //[Ndata.bb.gem.hit.ylocal]
   Int_t           Ndata_bb_gem_hit_zglobal;
   Double_t        bb_gem_hit_zglobal[5];   //[Ndata.bb.gem.hit.zglobal]
   Int_t           Ndata_bb_gem_m0_clust_clustu_adc;
   Double_t        bb_gem_m0_clust_clustu_adc[19];   //[Ndata.bb.gem.m0.clust.clustu_adc]
   Int_t           Ndata_bb_gem_m0_clust_clustu_pos;
   Double_t        bb_gem_m0_clust_clustu_pos[19];   //[Ndata.bb.gem.m0.clust.clustu_pos]
   Int_t           Ndata_bb_gem_m0_clust_clustu_strips;
   Double_t        bb_gem_m0_clust_clustu_strips[19];   //[Ndata.bb.gem.m0.clust.clustu_strips]
   Int_t           Ndata_bb_gem_m0_clust_clustu_time;
   Double_t        bb_gem_m0_clust_clustu_time[19];   //[Ndata.bb.gem.m0.clust.clustu_time]
   Int_t           Ndata_bb_gem_m0_clust_clustv_adc;
   Double_t        bb_gem_m0_clust_clustv_adc[18];   //[Ndata.bb.gem.m0.clust.clustv_adc]
   Int_t           Ndata_bb_gem_m0_clust_clustv_pos;
   Double_t        bb_gem_m0_clust_clustv_pos[18];   //[Ndata.bb.gem.m0.clust.clustv_pos]
   Int_t           Ndata_bb_gem_m0_clust_clustv_strips;
   Double_t        bb_gem_m0_clust_clustv_strips[18];   //[Ndata.bb.gem.m0.clust.clustv_strips]
   Int_t           Ndata_bb_gem_m0_clust_clustv_time;
   Double_t        bb_gem_m0_clust_clustv_time[18];   //[Ndata.bb.gem.m0.clust.clustv_time]
   Int_t           Ndata_bb_gem_m0_clust_isnegativeU;
   Double_t        bb_gem_m0_clust_isnegativeU[19];   //[Ndata.bb.gem.m0.clust.isnegativeU]
   Int_t           Ndata_bb_gem_m0_clust_isnegativeV;
   Double_t        bb_gem_m0_clust_isnegativeV[18];   //[Ndata.bb.gem.m0.clust.isnegativeV]
   Int_t           Ndata_bb_gem_m0_clust_isnegontrackU;
   Double_t        bb_gem_m0_clust_isnegontrackU[19];   //[Ndata.bb.gem.m0.clust.isnegontrackU]
   Int_t           Ndata_bb_gem_m0_clust_isnegontrackV;
   Double_t        bb_gem_m0_clust_isnegontrackV[18];   //[Ndata.bb.gem.m0.clust.isnegontrackV]
   Int_t           Ndata_bb_gem_m0_strip_ADCavg;
   Double_t        bb_gem_m0_strip_ADCavg[1138];   //[Ndata.bb.gem.m0.strip.ADCavg]
   Int_t           Ndata_bb_gem_m0_strip_ADCmax;
   Double_t        bb_gem_m0_strip_ADCmax[1138];   //[Ndata.bb.gem.m0.strip.ADCmax]
   Int_t           Ndata_bb_gem_m0_strip_ADCsamples;
   Double_t        bb_gem_m0_strip_ADCsamples[6828];   //[Ndata.bb.gem.m0.strip.ADCsamples]
   Int_t           Ndata_bb_gem_m0_strip_ADCsum;
   Double_t        bb_gem_m0_strip_ADCsum[1138];   //[Ndata.bb.gem.m0.strip.ADCsum]
   Int_t           Ndata_bb_gem_m0_strip_BUILD_ALL_SAMPLES;
   Double_t        bb_gem_m0_strip_BUILD_ALL_SAMPLES[1138];   //[Ndata.bb.gem.m0.strip.BUILD_ALL_SAMPLES]
   Int_t           Ndata_bb_gem_m0_strip_CM_GOOD;
   Double_t        bb_gem_m0_strip_CM_GOOD[1138];   //[Ndata.bb.gem.m0.strip.CM_GOOD]
   Int_t           Ndata_bb_gem_m0_strip_CorrCoeff;
   Double_t        bb_gem_m0_strip_CorrCoeff[1138];   //[Ndata.bb.gem.m0.strip.CorrCoeff]
   Int_t           Ndata_bb_gem_m0_strip_DeconvADCmax;
   Double_t        bb_gem_m0_strip_DeconvADCmax[1138];   //[Ndata.bb.gem.m0.strip.DeconvADCmax]
   Int_t           Ndata_bb_gem_m0_strip_DeconvADCmaxCombo;
   Double_t        bb_gem_m0_strip_DeconvADCmaxCombo[1138];   //[Ndata.bb.gem.m0.strip.DeconvADCmaxCombo]
   Int_t           Ndata_bb_gem_m0_strip_DeconvADCsamples;
   Double_t        bb_gem_m0_strip_DeconvADCsamples[6828];   //[Ndata.bb.gem.m0.strip.DeconvADCsamples]
   Int_t           Ndata_bb_gem_m0_strip_DeconvADCsum;
   Double_t        bb_gem_m0_strip_DeconvADCsum[1138];   //[Ndata.bb.gem.m0.strip.DeconvADCsum]
   Int_t           Ndata_bb_gem_m0_strip_ENABLE_CM;
   Double_t        bb_gem_m0_strip_ENABLE_CM[1138];   //[Ndata.bb.gem.m0.strip.ENABLE_CM]
   Int_t           Ndata_bb_gem_m0_strip_IsU;
   Double_t        bb_gem_m0_strip_IsU[1138];   //[Ndata.bb.gem.m0.strip.IsU]
   Int_t           Ndata_bb_gem_m0_strip_IsV;
   Double_t        bb_gem_m0_strip_IsV[1138];   //[Ndata.bb.gem.m0.strip.IsV]
   Int_t           Ndata_bb_gem_m0_strip_TSchi2;
   Double_t        bb_gem_m0_strip_TSchi2[1138];   //[Ndata.bb.gem.m0.strip.TSchi2]
   Int_t           Ndata_bb_gem_m0_strip_TSprob;
   Double_t        bb_gem_m0_strip_TSprob[1138];   //[Ndata.bb.gem.m0.strip.TSprob]
   Int_t           Ndata_bb_gem_m0_strip_Tcorr;
   Double_t        bb_gem_m0_strip_Tcorr[1138];   //[Ndata.bb.gem.m0.strip.Tcorr]
   Int_t           Ndata_bb_gem_m0_strip_Tdiff;
   Double_t        bb_gem_m0_strip_Tdiff[1138];   //[Ndata.bb.gem.m0.strip.Tdiff]
   Int_t           Ndata_bb_gem_m0_strip_Tfit;
   Double_t        bb_gem_m0_strip_Tfit[1138];   //[Ndata.bb.gem.m0.strip.Tfit]
   Int_t           Ndata_bb_gem_m0_strip_Tmean;
   Double_t        bb_gem_m0_strip_Tmean[1138];   //[Ndata.bb.gem.m0.strip.Tmean]
   Int_t           Ndata_bb_gem_m0_strip_TmeanDeconv;
   Double_t        bb_gem_m0_strip_TmeanDeconv[1138];   //[Ndata.bb.gem.m0.strip.TmeanDeconv]
   Int_t           Ndata_bb_gem_m0_strip_Tsigma;
   Double_t        bb_gem_m0_strip_Tsigma[1138];   //[Ndata.bb.gem.m0.strip.Tsigma]
   Int_t           Ndata_bb_gem_m0_strip_adc_id;
   Double_t        bb_gem_m0_strip_adc_id[1138];   //[Ndata.bb.gem.m0.strip.adc_id]
   Int_t           Ndata_bb_gem_m0_strip_crate;
   Double_t        bb_gem_m0_strip_crate[1138];   //[Ndata.bb.gem.m0.strip.crate]
   Int_t           Ndata_bb_gem_m0_strip_isampmax;
   Double_t        bb_gem_m0_strip_isampmax[1138];   //[Ndata.bb.gem.m0.strip.isampmax]
   Int_t           Ndata_bb_gem_m0_strip_isampmaxDeconv;
   Double_t        bb_gem_m0_strip_isampmaxDeconv[1138];   //[Ndata.bb.gem.m0.strip.isampmaxDeconv]
   Int_t           Ndata_bb_gem_m0_strip_isampmaxDeconvCombo;
   Double_t        bb_gem_m0_strip_isampmaxDeconvCombo[1138];   //[Ndata.bb.gem.m0.strip.isampmaxDeconvCombo]
   Int_t           Ndata_bb_gem_m0_strip_isnegative;
   Double_t        bb_gem_m0_strip_isnegative[1138];   //[Ndata.bb.gem.m0.strip.isnegative]
   Int_t           Ndata_bb_gem_m0_strip_isnegativeU;
   Double_t        bb_gem_m0_strip_isnegativeU[1138];   //[Ndata.bb.gem.m0.strip.isnegativeU]
   Int_t           Ndata_bb_gem_m0_strip_isnegativeV;
   Double_t        bb_gem_m0_strip_isnegativeV[1138];   //[Ndata.bb.gem.m0.strip.isnegativeV]
   Int_t           Ndata_bb_gem_m0_strip_isnegontrack;
   Double_t        bb_gem_m0_strip_isnegontrack[1138];   //[Ndata.bb.gem.m0.strip.isnegontrack]
   Int_t           Ndata_bb_gem_m0_strip_isnegontrackU;
   Double_t        bb_gem_m0_strip_isnegontrackU[1138];   //[Ndata.bb.gem.m0.strip.isnegontrackU]
   Int_t           Ndata_bb_gem_m0_strip_isnegontrackV;
   Double_t        bb_gem_m0_strip_isnegontrackV[1138];   //[Ndata.bb.gem.m0.strip.isnegontrackV]
   Int_t           Ndata_bb_gem_m0_strip_istrip;
   Double_t        bb_gem_m0_strip_istrip[1138];   //[Ndata.bb.gem.m0.strip.istrip]
   Int_t           Ndata_bb_gem_m0_strip_itrack;
   Double_t        bb_gem_m0_strip_itrack[1138];   //[Ndata.bb.gem.m0.strip.itrack]
   Int_t           Ndata_bb_gem_m0_strip_mpd;
   Double_t        bb_gem_m0_strip_mpd[1138];   //[Ndata.bb.gem.m0.strip.mpd]
   Int_t           Ndata_bb_gem_m0_strip_ontrack;
   Double_t        bb_gem_m0_strip_ontrack[1138];   //[Ndata.bb.gem.m0.strip.ontrack]
   Int_t           Ndata_bb_gem_m0_strip_ontrackU;
   Double_t        bb_gem_m0_strip_ontrackU[1138];   //[Ndata.bb.gem.m0.strip.ontrackU]
   Int_t           Ndata_bb_gem_m0_strip_ontrackV;
   Double_t        bb_gem_m0_strip_ontrackV[1138];   //[Ndata.bb.gem.m0.strip.ontrackV]
   Int_t           Ndata_bb_gem_m0_strip_rawADCsamples;
   Double_t        bb_gem_m0_strip_rawADCsamples[6828];   //[Ndata.bb.gem.m0.strip.rawADCsamples]
   Int_t           Ndata_bb_gem_m1_clust_clustu_adc;
   Double_t        bb_gem_m1_clust_clustu_adc[28];   //[Ndata.bb.gem.m1.clust.clustu_adc]
   Int_t           Ndata_bb_gem_m1_clust_clustu_pos;
   Double_t        bb_gem_m1_clust_clustu_pos[28];   //[Ndata.bb.gem.m1.clust.clustu_pos]
   Int_t           Ndata_bb_gem_m1_clust_clustu_strips;
   Double_t        bb_gem_m1_clust_clustu_strips[28];   //[Ndata.bb.gem.m1.clust.clustu_strips]
   Int_t           Ndata_bb_gem_m1_clust_clustu_time;
   Double_t        bb_gem_m1_clust_clustu_time[28];   //[Ndata.bb.gem.m1.clust.clustu_time]
   Int_t           Ndata_bb_gem_m1_clust_clustv_adc;
   Double_t        bb_gem_m1_clust_clustv_adc[37];   //[Ndata.bb.gem.m1.clust.clustv_adc]
   Int_t           Ndata_bb_gem_m1_clust_clustv_pos;
   Double_t        bb_gem_m1_clust_clustv_pos[37];   //[Ndata.bb.gem.m1.clust.clustv_pos]
   Int_t           Ndata_bb_gem_m1_clust_clustv_strips;
   Double_t        bb_gem_m1_clust_clustv_strips[37];   //[Ndata.bb.gem.m1.clust.clustv_strips]
   Int_t           Ndata_bb_gem_m1_clust_clustv_time;
   Double_t        bb_gem_m1_clust_clustv_time[37];   //[Ndata.bb.gem.m1.clust.clustv_time]
   Int_t           Ndata_bb_gem_m1_clust_isnegativeU;
   Double_t        bb_gem_m1_clust_isnegativeU[28];   //[Ndata.bb.gem.m1.clust.isnegativeU]
   Int_t           Ndata_bb_gem_m1_clust_isnegativeV;
   Double_t        bb_gem_m1_clust_isnegativeV[37];   //[Ndata.bb.gem.m1.clust.isnegativeV]
   Int_t           Ndata_bb_gem_m1_clust_isnegontrackU;
   Double_t        bb_gem_m1_clust_isnegontrackU[28];   //[Ndata.bb.gem.m1.clust.isnegontrackU]
   Int_t           Ndata_bb_gem_m1_clust_isnegontrackV;
   Double_t        bb_gem_m1_clust_isnegontrackV[37];   //[Ndata.bb.gem.m1.clust.isnegontrackV]
   Int_t           Ndata_bb_gem_m1_strip_ADCavg;
   Double_t        bb_gem_m1_strip_ADCavg[1143];   //[Ndata.bb.gem.m1.strip.ADCavg]
   Int_t           Ndata_bb_gem_m1_strip_ADCmax;
   Double_t        bb_gem_m1_strip_ADCmax[1143];   //[Ndata.bb.gem.m1.strip.ADCmax]
   Int_t           Ndata_bb_gem_m1_strip_ADCsamples;
   Double_t        bb_gem_m1_strip_ADCsamples[6858];   //[Ndata.bb.gem.m1.strip.ADCsamples]
   Int_t           Ndata_bb_gem_m1_strip_ADCsum;
   Double_t        bb_gem_m1_strip_ADCsum[1143];   //[Ndata.bb.gem.m1.strip.ADCsum]
   Int_t           Ndata_bb_gem_m1_strip_BUILD_ALL_SAMPLES;
   Double_t        bb_gem_m1_strip_BUILD_ALL_SAMPLES[1143];   //[Ndata.bb.gem.m1.strip.BUILD_ALL_SAMPLES]
   Int_t           Ndata_bb_gem_m1_strip_CM_GOOD;
   Double_t        bb_gem_m1_strip_CM_GOOD[1143];   //[Ndata.bb.gem.m1.strip.CM_GOOD]
   Int_t           Ndata_bb_gem_m1_strip_CorrCoeff;
   Double_t        bb_gem_m1_strip_CorrCoeff[1143];   //[Ndata.bb.gem.m1.strip.CorrCoeff]
   Int_t           Ndata_bb_gem_m1_strip_DeconvADCmax;
   Double_t        bb_gem_m1_strip_DeconvADCmax[1143];   //[Ndata.bb.gem.m1.strip.DeconvADCmax]
   Int_t           Ndata_bb_gem_m1_strip_DeconvADCmaxCombo;
   Double_t        bb_gem_m1_strip_DeconvADCmaxCombo[1143];   //[Ndata.bb.gem.m1.strip.DeconvADCmaxCombo]
   Int_t           Ndata_bb_gem_m1_strip_DeconvADCsamples;
   Double_t        bb_gem_m1_strip_DeconvADCsamples[6858];   //[Ndata.bb.gem.m1.strip.DeconvADCsamples]
   Int_t           Ndata_bb_gem_m1_strip_DeconvADCsum;
   Double_t        bb_gem_m1_strip_DeconvADCsum[1143];   //[Ndata.bb.gem.m1.strip.DeconvADCsum]
   Int_t           Ndata_bb_gem_m1_strip_ENABLE_CM;
   Double_t        bb_gem_m1_strip_ENABLE_CM[1143];   //[Ndata.bb.gem.m1.strip.ENABLE_CM]
   Int_t           Ndata_bb_gem_m1_strip_IsU;
   Double_t        bb_gem_m1_strip_IsU[1143];   //[Ndata.bb.gem.m1.strip.IsU]
   Int_t           Ndata_bb_gem_m1_strip_IsV;
   Double_t        bb_gem_m1_strip_IsV[1143];   //[Ndata.bb.gem.m1.strip.IsV]
   Int_t           Ndata_bb_gem_m1_strip_TSchi2;
   Double_t        bb_gem_m1_strip_TSchi2[1143];   //[Ndata.bb.gem.m1.strip.TSchi2]
   Int_t           Ndata_bb_gem_m1_strip_TSprob;
   Double_t        bb_gem_m1_strip_TSprob[1143];   //[Ndata.bb.gem.m1.strip.TSprob]
   Int_t           Ndata_bb_gem_m1_strip_Tcorr;
   Double_t        bb_gem_m1_strip_Tcorr[1143];   //[Ndata.bb.gem.m1.strip.Tcorr]
   Int_t           Ndata_bb_gem_m1_strip_Tdiff;
   Double_t        bb_gem_m1_strip_Tdiff[1143];   //[Ndata.bb.gem.m1.strip.Tdiff]
   Int_t           Ndata_bb_gem_m1_strip_Tfit;
   Double_t        bb_gem_m1_strip_Tfit[1143];   //[Ndata.bb.gem.m1.strip.Tfit]
   Int_t           Ndata_bb_gem_m1_strip_Tmean;
   Double_t        bb_gem_m1_strip_Tmean[1143];   //[Ndata.bb.gem.m1.strip.Tmean]
   Int_t           Ndata_bb_gem_m1_strip_TmeanDeconv;
   Double_t        bb_gem_m1_strip_TmeanDeconv[1143];   //[Ndata.bb.gem.m1.strip.TmeanDeconv]
   Int_t           Ndata_bb_gem_m1_strip_Tsigma;
   Double_t        bb_gem_m1_strip_Tsigma[1143];   //[Ndata.bb.gem.m1.strip.Tsigma]
   Int_t           Ndata_bb_gem_m1_strip_adc_id;
   Double_t        bb_gem_m1_strip_adc_id[1143];   //[Ndata.bb.gem.m1.strip.adc_id]
   Int_t           Ndata_bb_gem_m1_strip_crate;
   Double_t        bb_gem_m1_strip_crate[1143];   //[Ndata.bb.gem.m1.strip.crate]
   Int_t           Ndata_bb_gem_m1_strip_isampmax;
   Double_t        bb_gem_m1_strip_isampmax[1143];   //[Ndata.bb.gem.m1.strip.isampmax]
   Int_t           Ndata_bb_gem_m1_strip_isampmaxDeconv;
   Double_t        bb_gem_m1_strip_isampmaxDeconv[1143];   //[Ndata.bb.gem.m1.strip.isampmaxDeconv]
   Int_t           Ndata_bb_gem_m1_strip_isampmaxDeconvCombo;
   Double_t        bb_gem_m1_strip_isampmaxDeconvCombo[1143];   //[Ndata.bb.gem.m1.strip.isampmaxDeconvCombo]
   Int_t           Ndata_bb_gem_m1_strip_isnegative;
   Double_t        bb_gem_m1_strip_isnegative[1143];   //[Ndata.bb.gem.m1.strip.isnegative]
   Int_t           Ndata_bb_gem_m1_strip_isnegativeU;
   Double_t        bb_gem_m1_strip_isnegativeU[1143];   //[Ndata.bb.gem.m1.strip.isnegativeU]
   Int_t           Ndata_bb_gem_m1_strip_isnegativeV;
   Double_t        bb_gem_m1_strip_isnegativeV[1143];   //[Ndata.bb.gem.m1.strip.isnegativeV]
   Int_t           Ndata_bb_gem_m1_strip_isnegontrack;
   Double_t        bb_gem_m1_strip_isnegontrack[1143];   //[Ndata.bb.gem.m1.strip.isnegontrack]
   Int_t           Ndata_bb_gem_m1_strip_isnegontrackU;
   Double_t        bb_gem_m1_strip_isnegontrackU[1143];   //[Ndata.bb.gem.m1.strip.isnegontrackU]
   Int_t           Ndata_bb_gem_m1_strip_isnegontrackV;
   Double_t        bb_gem_m1_strip_isnegontrackV[1143];   //[Ndata.bb.gem.m1.strip.isnegontrackV]
   Int_t           Ndata_bb_gem_m1_strip_istrip;
   Double_t        bb_gem_m1_strip_istrip[1143];   //[Ndata.bb.gem.m1.strip.istrip]
   Int_t           Ndata_bb_gem_m1_strip_itrack;
   Double_t        bb_gem_m1_strip_itrack[1143];   //[Ndata.bb.gem.m1.strip.itrack]
   Int_t           Ndata_bb_gem_m1_strip_mpd;
   Double_t        bb_gem_m1_strip_mpd[1143];   //[Ndata.bb.gem.m1.strip.mpd]
   Int_t           Ndata_bb_gem_m1_strip_ontrack;
   Double_t        bb_gem_m1_strip_ontrack[1143];   //[Ndata.bb.gem.m1.strip.ontrack]
   Int_t           Ndata_bb_gem_m1_strip_ontrackU;
   Double_t        bb_gem_m1_strip_ontrackU[1143];   //[Ndata.bb.gem.m1.strip.ontrackU]
   Int_t           Ndata_bb_gem_m1_strip_ontrackV;
   Double_t        bb_gem_m1_strip_ontrackV[1143];   //[Ndata.bb.gem.m1.strip.ontrackV]
   Int_t           Ndata_bb_gem_m1_strip_rawADCsamples;
   Double_t        bb_gem_m1_strip_rawADCsamples[6858];   //[Ndata.bb.gem.m1.strip.rawADCsamples]
   Int_t           Ndata_bb_gem_m2_clust_clustu_adc;
   Double_t        bb_gem_m2_clust_clustu_adc[15];   //[Ndata.bb.gem.m2.clust.clustu_adc]
   Int_t           Ndata_bb_gem_m2_clust_clustu_pos;
   Double_t        bb_gem_m2_clust_clustu_pos[15];   //[Ndata.bb.gem.m2.clust.clustu_pos]
   Int_t           Ndata_bb_gem_m2_clust_clustu_strips;
   Double_t        bb_gem_m2_clust_clustu_strips[15];   //[Ndata.bb.gem.m2.clust.clustu_strips]
   Int_t           Ndata_bb_gem_m2_clust_clustu_time;
   Double_t        bb_gem_m2_clust_clustu_time[15];   //[Ndata.bb.gem.m2.clust.clustu_time]
   Int_t           Ndata_bb_gem_m2_clust_clustv_adc;
   Double_t        bb_gem_m2_clust_clustv_adc[24];   //[Ndata.bb.gem.m2.clust.clustv_adc]
   Int_t           Ndata_bb_gem_m2_clust_clustv_pos;
   Double_t        bb_gem_m2_clust_clustv_pos[24];   //[Ndata.bb.gem.m2.clust.clustv_pos]
   Int_t           Ndata_bb_gem_m2_clust_clustv_strips;
   Double_t        bb_gem_m2_clust_clustv_strips[24];   //[Ndata.bb.gem.m2.clust.clustv_strips]
   Int_t           Ndata_bb_gem_m2_clust_clustv_time;
   Double_t        bb_gem_m2_clust_clustv_time[24];   //[Ndata.bb.gem.m2.clust.clustv_time]
   Int_t           Ndata_bb_gem_m2_clust_isnegativeU;
   Double_t        bb_gem_m2_clust_isnegativeU[15];   //[Ndata.bb.gem.m2.clust.isnegativeU]
   Int_t           Ndata_bb_gem_m2_clust_isnegativeV;
   Double_t        bb_gem_m2_clust_isnegativeV[24];   //[Ndata.bb.gem.m2.clust.isnegativeV]
   Int_t           Ndata_bb_gem_m2_clust_isnegontrackU;
   Double_t        bb_gem_m2_clust_isnegontrackU[15];   //[Ndata.bb.gem.m2.clust.isnegontrackU]
   Int_t           Ndata_bb_gem_m2_clust_isnegontrackV;
   Double_t        bb_gem_m2_clust_isnegontrackV[24];   //[Ndata.bb.gem.m2.clust.isnegontrackV]
   Int_t           Ndata_bb_gem_m2_strip_ADCavg;
   Double_t        bb_gem_m2_strip_ADCavg[1068];   //[Ndata.bb.gem.m2.strip.ADCavg]
   Int_t           Ndata_bb_gem_m2_strip_ADCmax;
   Double_t        bb_gem_m2_strip_ADCmax[1068];   //[Ndata.bb.gem.m2.strip.ADCmax]
   Int_t           Ndata_bb_gem_m2_strip_ADCsamples;
   Double_t        bb_gem_m2_strip_ADCsamples[6408];   //[Ndata.bb.gem.m2.strip.ADCsamples]
   Int_t           Ndata_bb_gem_m2_strip_ADCsum;
   Double_t        bb_gem_m2_strip_ADCsum[1068];   //[Ndata.bb.gem.m2.strip.ADCsum]
   Int_t           Ndata_bb_gem_m2_strip_BUILD_ALL_SAMPLES;
   Double_t        bb_gem_m2_strip_BUILD_ALL_SAMPLES[1068];   //[Ndata.bb.gem.m2.strip.BUILD_ALL_SAMPLES]
   Int_t           Ndata_bb_gem_m2_strip_CM_GOOD;
   Double_t        bb_gem_m2_strip_CM_GOOD[1068];   //[Ndata.bb.gem.m2.strip.CM_GOOD]
   Int_t           Ndata_bb_gem_m2_strip_CorrCoeff;
   Double_t        bb_gem_m2_strip_CorrCoeff[1068];   //[Ndata.bb.gem.m2.strip.CorrCoeff]
   Int_t           Ndata_bb_gem_m2_strip_DeconvADCmax;
   Double_t        bb_gem_m2_strip_DeconvADCmax[1068];   //[Ndata.bb.gem.m2.strip.DeconvADCmax]
   Int_t           Ndata_bb_gem_m2_strip_DeconvADCmaxCombo;
   Double_t        bb_gem_m2_strip_DeconvADCmaxCombo[1068];   //[Ndata.bb.gem.m2.strip.DeconvADCmaxCombo]
   Int_t           Ndata_bb_gem_m2_strip_DeconvADCsamples;
   Double_t        bb_gem_m2_strip_DeconvADCsamples[6408];   //[Ndata.bb.gem.m2.strip.DeconvADCsamples]
   Int_t           Ndata_bb_gem_m2_strip_DeconvADCsum;
   Double_t        bb_gem_m2_strip_DeconvADCsum[1068];   //[Ndata.bb.gem.m2.strip.DeconvADCsum]
   Int_t           Ndata_bb_gem_m2_strip_ENABLE_CM;
   Double_t        bb_gem_m2_strip_ENABLE_CM[1068];   //[Ndata.bb.gem.m2.strip.ENABLE_CM]
   Int_t           Ndata_bb_gem_m2_strip_IsU;
   Double_t        bb_gem_m2_strip_IsU[1068];   //[Ndata.bb.gem.m2.strip.IsU]
   Int_t           Ndata_bb_gem_m2_strip_IsV;
   Double_t        bb_gem_m2_strip_IsV[1068];   //[Ndata.bb.gem.m2.strip.IsV]
   Int_t           Ndata_bb_gem_m2_strip_TSchi2;
   Double_t        bb_gem_m2_strip_TSchi2[1068];   //[Ndata.bb.gem.m2.strip.TSchi2]
   Int_t           Ndata_bb_gem_m2_strip_TSprob;
   Double_t        bb_gem_m2_strip_TSprob[1068];   //[Ndata.bb.gem.m2.strip.TSprob]
   Int_t           Ndata_bb_gem_m2_strip_Tcorr;
   Double_t        bb_gem_m2_strip_Tcorr[1068];   //[Ndata.bb.gem.m2.strip.Tcorr]
   Int_t           Ndata_bb_gem_m2_strip_Tdiff;
   Double_t        bb_gem_m2_strip_Tdiff[1068];   //[Ndata.bb.gem.m2.strip.Tdiff]
   Int_t           Ndata_bb_gem_m2_strip_Tfit;
   Double_t        bb_gem_m2_strip_Tfit[1068];   //[Ndata.bb.gem.m2.strip.Tfit]
   Int_t           Ndata_bb_gem_m2_strip_Tmean;
   Double_t        bb_gem_m2_strip_Tmean[1068];   //[Ndata.bb.gem.m2.strip.Tmean]
   Int_t           Ndata_bb_gem_m2_strip_TmeanDeconv;
   Double_t        bb_gem_m2_strip_TmeanDeconv[1068];   //[Ndata.bb.gem.m2.strip.TmeanDeconv]
   Int_t           Ndata_bb_gem_m2_strip_Tsigma;
   Double_t        bb_gem_m2_strip_Tsigma[1068];   //[Ndata.bb.gem.m2.strip.Tsigma]
   Int_t           Ndata_bb_gem_m2_strip_adc_id;
   Double_t        bb_gem_m2_strip_adc_id[1068];   //[Ndata.bb.gem.m2.strip.adc_id]
   Int_t           Ndata_bb_gem_m2_strip_crate;
   Double_t        bb_gem_m2_strip_crate[1068];   //[Ndata.bb.gem.m2.strip.crate]
   Int_t           Ndata_bb_gem_m2_strip_isampmax;
   Double_t        bb_gem_m2_strip_isampmax[1068];   //[Ndata.bb.gem.m2.strip.isampmax]
   Int_t           Ndata_bb_gem_m2_strip_isampmaxDeconv;
   Double_t        bb_gem_m2_strip_isampmaxDeconv[1068];   //[Ndata.bb.gem.m2.strip.isampmaxDeconv]
   Int_t           Ndata_bb_gem_m2_strip_isampmaxDeconvCombo;
   Double_t        bb_gem_m2_strip_isampmaxDeconvCombo[1068];   //[Ndata.bb.gem.m2.strip.isampmaxDeconvCombo]
   Int_t           Ndata_bb_gem_m2_strip_isnegative;
   Double_t        bb_gem_m2_strip_isnegative[1068];   //[Ndata.bb.gem.m2.strip.isnegative]
   Int_t           Ndata_bb_gem_m2_strip_isnegativeU;
   Double_t        bb_gem_m2_strip_isnegativeU[1068];   //[Ndata.bb.gem.m2.strip.isnegativeU]
   Int_t           Ndata_bb_gem_m2_strip_isnegativeV;
   Double_t        bb_gem_m2_strip_isnegativeV[1068];   //[Ndata.bb.gem.m2.strip.isnegativeV]
   Int_t           Ndata_bb_gem_m2_strip_isnegontrack;
   Double_t        bb_gem_m2_strip_isnegontrack[1068];   //[Ndata.bb.gem.m2.strip.isnegontrack]
   Int_t           Ndata_bb_gem_m2_strip_isnegontrackU;
   Double_t        bb_gem_m2_strip_isnegontrackU[1068];   //[Ndata.bb.gem.m2.strip.isnegontrackU]
   Int_t           Ndata_bb_gem_m2_strip_isnegontrackV;
   Double_t        bb_gem_m2_strip_isnegontrackV[1068];   //[Ndata.bb.gem.m2.strip.isnegontrackV]
   Int_t           Ndata_bb_gem_m2_strip_istrip;
   Double_t        bb_gem_m2_strip_istrip[1068];   //[Ndata.bb.gem.m2.strip.istrip]
   Int_t           Ndata_bb_gem_m2_strip_itrack;
   Double_t        bb_gem_m2_strip_itrack[1068];   //[Ndata.bb.gem.m2.strip.itrack]
   Int_t           Ndata_bb_gem_m2_strip_mpd;
   Double_t        bb_gem_m2_strip_mpd[1068];   //[Ndata.bb.gem.m2.strip.mpd]
   Int_t           Ndata_bb_gem_m2_strip_ontrack;
   Double_t        bb_gem_m2_strip_ontrack[1068];   //[Ndata.bb.gem.m2.strip.ontrack]
   Int_t           Ndata_bb_gem_m2_strip_ontrackU;
   Double_t        bb_gem_m2_strip_ontrackU[1068];   //[Ndata.bb.gem.m2.strip.ontrackU]
   Int_t           Ndata_bb_gem_m2_strip_ontrackV;
   Double_t        bb_gem_m2_strip_ontrackV[1068];   //[Ndata.bb.gem.m2.strip.ontrackV]
   Int_t           Ndata_bb_gem_m2_strip_rawADCsamples;
   Double_t        bb_gem_m2_strip_rawADCsamples[6408];   //[Ndata.bb.gem.m2.strip.rawADCsamples]
   Int_t           Ndata_bb_gem_m3_clust_clustu_adc;
   Double_t        bb_gem_m3_clust_clustu_adc[23];   //[Ndata.bb.gem.m3.clust.clustu_adc]
   Int_t           Ndata_bb_gem_m3_clust_clustu_pos;
   Double_t        bb_gem_m3_clust_clustu_pos[23];   //[Ndata.bb.gem.m3.clust.clustu_pos]
   Int_t           Ndata_bb_gem_m3_clust_clustu_strips;
   Double_t        bb_gem_m3_clust_clustu_strips[23];   //[Ndata.bb.gem.m3.clust.clustu_strips]
   Int_t           Ndata_bb_gem_m3_clust_clustu_time;
   Double_t        bb_gem_m3_clust_clustu_time[23];   //[Ndata.bb.gem.m3.clust.clustu_time]
   Int_t           Ndata_bb_gem_m3_clust_clustv_adc;
   Double_t        bb_gem_m3_clust_clustv_adc[36];   //[Ndata.bb.gem.m3.clust.clustv_adc]
   Int_t           Ndata_bb_gem_m3_clust_clustv_pos;
   Double_t        bb_gem_m3_clust_clustv_pos[36];   //[Ndata.bb.gem.m3.clust.clustv_pos]
   Int_t           Ndata_bb_gem_m3_clust_clustv_strips;
   Double_t        bb_gem_m3_clust_clustv_strips[36];   //[Ndata.bb.gem.m3.clust.clustv_strips]
   Int_t           Ndata_bb_gem_m3_clust_clustv_time;
   Double_t        bb_gem_m3_clust_clustv_time[36];   //[Ndata.bb.gem.m3.clust.clustv_time]
   Int_t           Ndata_bb_gem_m3_clust_isnegativeU;
   Double_t        bb_gem_m3_clust_isnegativeU[23];   //[Ndata.bb.gem.m3.clust.isnegativeU]
   Int_t           Ndata_bb_gem_m3_clust_isnegativeV;
   Double_t        bb_gem_m3_clust_isnegativeV[36];   //[Ndata.bb.gem.m3.clust.isnegativeV]
   Int_t           Ndata_bb_gem_m3_clust_isnegontrackU;
   Double_t        bb_gem_m3_clust_isnegontrackU[23];   //[Ndata.bb.gem.m3.clust.isnegontrackU]
   Int_t           Ndata_bb_gem_m3_clust_isnegontrackV;
   Double_t        bb_gem_m3_clust_isnegontrackV[36];   //[Ndata.bb.gem.m3.clust.isnegontrackV]
   Int_t           Ndata_bb_gem_m3_strip_ADCavg;
   Double_t        bb_gem_m3_strip_ADCavg[1128];   //[Ndata.bb.gem.m3.strip.ADCavg]
   Int_t           Ndata_bb_gem_m3_strip_ADCmax;
   Double_t        bb_gem_m3_strip_ADCmax[1128];   //[Ndata.bb.gem.m3.strip.ADCmax]
   Int_t           Ndata_bb_gem_m3_strip_ADCsamples;
   Double_t        bb_gem_m3_strip_ADCsamples[6768];   //[Ndata.bb.gem.m3.strip.ADCsamples]
   Int_t           Ndata_bb_gem_m3_strip_ADCsum;
   Double_t        bb_gem_m3_strip_ADCsum[1128];   //[Ndata.bb.gem.m3.strip.ADCsum]
   Int_t           Ndata_bb_gem_m3_strip_BUILD_ALL_SAMPLES;
   Double_t        bb_gem_m3_strip_BUILD_ALL_SAMPLES[1128];   //[Ndata.bb.gem.m3.strip.BUILD_ALL_SAMPLES]
   Int_t           Ndata_bb_gem_m3_strip_CM_GOOD;
   Double_t        bb_gem_m3_strip_CM_GOOD[1128];   //[Ndata.bb.gem.m3.strip.CM_GOOD]
   Int_t           Ndata_bb_gem_m3_strip_CorrCoeff;
   Double_t        bb_gem_m3_strip_CorrCoeff[1128];   //[Ndata.bb.gem.m3.strip.CorrCoeff]
   Int_t           Ndata_bb_gem_m3_strip_DeconvADCmax;
   Double_t        bb_gem_m3_strip_DeconvADCmax[1128];   //[Ndata.bb.gem.m3.strip.DeconvADCmax]
   Int_t           Ndata_bb_gem_m3_strip_DeconvADCmaxCombo;
   Double_t        bb_gem_m3_strip_DeconvADCmaxCombo[1128];   //[Ndata.bb.gem.m3.strip.DeconvADCmaxCombo]
   Int_t           Ndata_bb_gem_m3_strip_DeconvADCsamples;
   Double_t        bb_gem_m3_strip_DeconvADCsamples[6768];   //[Ndata.bb.gem.m3.strip.DeconvADCsamples]
   Int_t           Ndata_bb_gem_m3_strip_DeconvADCsum;
   Double_t        bb_gem_m3_strip_DeconvADCsum[1128];   //[Ndata.bb.gem.m3.strip.DeconvADCsum]
   Int_t           Ndata_bb_gem_m3_strip_ENABLE_CM;
   Double_t        bb_gem_m3_strip_ENABLE_CM[1128];   //[Ndata.bb.gem.m3.strip.ENABLE_CM]
   Int_t           Ndata_bb_gem_m3_strip_IsU;
   Double_t        bb_gem_m3_strip_IsU[1128];   //[Ndata.bb.gem.m3.strip.IsU]
   Int_t           Ndata_bb_gem_m3_strip_IsV;
   Double_t        bb_gem_m3_strip_IsV[1128];   //[Ndata.bb.gem.m3.strip.IsV]
   Int_t           Ndata_bb_gem_m3_strip_TSchi2;
   Double_t        bb_gem_m3_strip_TSchi2[1128];   //[Ndata.bb.gem.m3.strip.TSchi2]
   Int_t           Ndata_bb_gem_m3_strip_TSprob;
   Double_t        bb_gem_m3_strip_TSprob[1128];   //[Ndata.bb.gem.m3.strip.TSprob]
   Int_t           Ndata_bb_gem_m3_strip_Tcorr;
   Double_t        bb_gem_m3_strip_Tcorr[1128];   //[Ndata.bb.gem.m3.strip.Tcorr]
   Int_t           Ndata_bb_gem_m3_strip_Tdiff;
   Double_t        bb_gem_m3_strip_Tdiff[1128];   //[Ndata.bb.gem.m3.strip.Tdiff]
   Int_t           Ndata_bb_gem_m3_strip_Tfit;
   Double_t        bb_gem_m3_strip_Tfit[1128];   //[Ndata.bb.gem.m3.strip.Tfit]
   Int_t           Ndata_bb_gem_m3_strip_Tmean;
   Double_t        bb_gem_m3_strip_Tmean[1128];   //[Ndata.bb.gem.m3.strip.Tmean]
   Int_t           Ndata_bb_gem_m3_strip_TmeanDeconv;
   Double_t        bb_gem_m3_strip_TmeanDeconv[1128];   //[Ndata.bb.gem.m3.strip.TmeanDeconv]
   Int_t           Ndata_bb_gem_m3_strip_Tsigma;
   Double_t        bb_gem_m3_strip_Tsigma[1128];   //[Ndata.bb.gem.m3.strip.Tsigma]
   Int_t           Ndata_bb_gem_m3_strip_adc_id;
   Double_t        bb_gem_m3_strip_adc_id[1128];   //[Ndata.bb.gem.m3.strip.adc_id]
   Int_t           Ndata_bb_gem_m3_strip_crate;
   Double_t        bb_gem_m3_strip_crate[1128];   //[Ndata.bb.gem.m3.strip.crate]
   Int_t           Ndata_bb_gem_m3_strip_isampmax;
   Double_t        bb_gem_m3_strip_isampmax[1128];   //[Ndata.bb.gem.m3.strip.isampmax]
   Int_t           Ndata_bb_gem_m3_strip_isampmaxDeconv;
   Double_t        bb_gem_m3_strip_isampmaxDeconv[1128];   //[Ndata.bb.gem.m3.strip.isampmaxDeconv]
   Int_t           Ndata_bb_gem_m3_strip_isampmaxDeconvCombo;
   Double_t        bb_gem_m3_strip_isampmaxDeconvCombo[1128];   //[Ndata.bb.gem.m3.strip.isampmaxDeconvCombo]
   Int_t           Ndata_bb_gem_m3_strip_isnegative;
   Double_t        bb_gem_m3_strip_isnegative[1128];   //[Ndata.bb.gem.m3.strip.isnegative]
   Int_t           Ndata_bb_gem_m3_strip_isnegativeU;
   Double_t        bb_gem_m3_strip_isnegativeU[1128];   //[Ndata.bb.gem.m3.strip.isnegativeU]
   Int_t           Ndata_bb_gem_m3_strip_isnegativeV;
   Double_t        bb_gem_m3_strip_isnegativeV[1128];   //[Ndata.bb.gem.m3.strip.isnegativeV]
   Int_t           Ndata_bb_gem_m3_strip_isnegontrack;
   Double_t        bb_gem_m3_strip_isnegontrack[1128];   //[Ndata.bb.gem.m3.strip.isnegontrack]
   Int_t           Ndata_bb_gem_m3_strip_isnegontrackU;
   Double_t        bb_gem_m3_strip_isnegontrackU[1128];   //[Ndata.bb.gem.m3.strip.isnegontrackU]
   Int_t           Ndata_bb_gem_m3_strip_isnegontrackV;
   Double_t        bb_gem_m3_strip_isnegontrackV[1128];   //[Ndata.bb.gem.m3.strip.isnegontrackV]
   Int_t           Ndata_bb_gem_m3_strip_istrip;
   Double_t        bb_gem_m3_strip_istrip[1128];   //[Ndata.bb.gem.m3.strip.istrip]
   Int_t           Ndata_bb_gem_m3_strip_itrack;
   Double_t        bb_gem_m3_strip_itrack[1128];   //[Ndata.bb.gem.m3.strip.itrack]
   Int_t           Ndata_bb_gem_m3_strip_mpd;
   Double_t        bb_gem_m3_strip_mpd[1128];   //[Ndata.bb.gem.m3.strip.mpd]
   Int_t           Ndata_bb_gem_m3_strip_ontrack;
   Double_t        bb_gem_m3_strip_ontrack[1128];   //[Ndata.bb.gem.m3.strip.ontrack]
   Int_t           Ndata_bb_gem_m3_strip_ontrackU;
   Double_t        bb_gem_m3_strip_ontrackU[1128];   //[Ndata.bb.gem.m3.strip.ontrackU]
   Int_t           Ndata_bb_gem_m3_strip_ontrackV;
   Double_t        bb_gem_m3_strip_ontrackV[1128];   //[Ndata.bb.gem.m3.strip.ontrackV]
   Int_t           Ndata_bb_gem_m3_strip_rawADCsamples;
   Double_t        bb_gem_m3_strip_rawADCsamples[6768];   //[Ndata.bb.gem.m3.strip.rawADCsamples]
   Int_t           Ndata_bb_gem_m4_clust_clustu_adc;
   Double_t        bb_gem_m4_clust_clustu_adc[3];   //[Ndata.bb.gem.m4.clust.clustu_adc]
   Int_t           Ndata_bb_gem_m4_clust_clustu_pos;
   Double_t        bb_gem_m4_clust_clustu_pos[3];   //[Ndata.bb.gem.m4.clust.clustu_pos]
   Int_t           Ndata_bb_gem_m4_clust_clustu_strips;
   Double_t        bb_gem_m4_clust_clustu_strips[3];   //[Ndata.bb.gem.m4.clust.clustu_strips]
   Int_t           Ndata_bb_gem_m4_clust_clustu_time;
   Double_t        bb_gem_m4_clust_clustu_time[3];   //[Ndata.bb.gem.m4.clust.clustu_time]
   Int_t           Ndata_bb_gem_m4_clust_clustv_adc;
   Double_t        bb_gem_m4_clust_clustv_adc[18];   //[Ndata.bb.gem.m4.clust.clustv_adc]
   Int_t           Ndata_bb_gem_m4_clust_clustv_pos;
   Double_t        bb_gem_m4_clust_clustv_pos[18];   //[Ndata.bb.gem.m4.clust.clustv_pos]
   Int_t           Ndata_bb_gem_m4_clust_clustv_strips;
   Double_t        bb_gem_m4_clust_clustv_strips[18];   //[Ndata.bb.gem.m4.clust.clustv_strips]
   Int_t           Ndata_bb_gem_m4_clust_clustv_time;
   Double_t        bb_gem_m4_clust_clustv_time[18];   //[Ndata.bb.gem.m4.clust.clustv_time]
   Int_t           Ndata_bb_gem_m4_clust_isnegativeU;
   Double_t        bb_gem_m4_clust_isnegativeU[3];   //[Ndata.bb.gem.m4.clust.isnegativeU]
   Int_t           Ndata_bb_gem_m4_clust_isnegativeV;
   Double_t        bb_gem_m4_clust_isnegativeV[18];   //[Ndata.bb.gem.m4.clust.isnegativeV]
   Int_t           Ndata_bb_gem_m4_clust_isnegontrackU;
   Double_t        bb_gem_m4_clust_isnegontrackU[3];   //[Ndata.bb.gem.m4.clust.isnegontrackU]
   Int_t           Ndata_bb_gem_m4_clust_isnegontrackV;
   Double_t        bb_gem_m4_clust_isnegontrackV[18];   //[Ndata.bb.gem.m4.clust.isnegontrackV]
   Int_t           Ndata_bb_gem_m4_strip_ADCavg;
   Double_t        bb_gem_m4_strip_ADCavg[756];   //[Ndata.bb.gem.m4.strip.ADCavg]
   Int_t           Ndata_bb_gem_m4_strip_ADCmax;
   Double_t        bb_gem_m4_strip_ADCmax[756];   //[Ndata.bb.gem.m4.strip.ADCmax]
   Int_t           Ndata_bb_gem_m4_strip_ADCsamples;
   Double_t        bb_gem_m4_strip_ADCsamples[4536];   //[Ndata.bb.gem.m4.strip.ADCsamples]
   Int_t           Ndata_bb_gem_m4_strip_ADCsum;
   Double_t        bb_gem_m4_strip_ADCsum[756];   //[Ndata.bb.gem.m4.strip.ADCsum]
   Int_t           Ndata_bb_gem_m4_strip_BUILD_ALL_SAMPLES;
   Double_t        bb_gem_m4_strip_BUILD_ALL_SAMPLES[756];   //[Ndata.bb.gem.m4.strip.BUILD_ALL_SAMPLES]
   Int_t           Ndata_bb_gem_m4_strip_CM_GOOD;
   Double_t        bb_gem_m4_strip_CM_GOOD[756];   //[Ndata.bb.gem.m4.strip.CM_GOOD]
   Int_t           Ndata_bb_gem_m4_strip_CorrCoeff;
   Double_t        bb_gem_m4_strip_CorrCoeff[756];   //[Ndata.bb.gem.m4.strip.CorrCoeff]
   Int_t           Ndata_bb_gem_m4_strip_DeconvADCmax;
   Double_t        bb_gem_m4_strip_DeconvADCmax[756];   //[Ndata.bb.gem.m4.strip.DeconvADCmax]
   Int_t           Ndata_bb_gem_m4_strip_DeconvADCmaxCombo;
   Double_t        bb_gem_m4_strip_DeconvADCmaxCombo[756];   //[Ndata.bb.gem.m4.strip.DeconvADCmaxCombo]
   Int_t           Ndata_bb_gem_m4_strip_DeconvADCsamples;
   Double_t        bb_gem_m4_strip_DeconvADCsamples[4536];   //[Ndata.bb.gem.m4.strip.DeconvADCsamples]
   Int_t           Ndata_bb_gem_m4_strip_DeconvADCsum;
   Double_t        bb_gem_m4_strip_DeconvADCsum[756];   //[Ndata.bb.gem.m4.strip.DeconvADCsum]
   Int_t           Ndata_bb_gem_m4_strip_ENABLE_CM;
   Double_t        bb_gem_m4_strip_ENABLE_CM[756];   //[Ndata.bb.gem.m4.strip.ENABLE_CM]
   Int_t           Ndata_bb_gem_m4_strip_IsU;
   Double_t        bb_gem_m4_strip_IsU[756];   //[Ndata.bb.gem.m4.strip.IsU]
   Int_t           Ndata_bb_gem_m4_strip_IsV;
   Double_t        bb_gem_m4_strip_IsV[756];   //[Ndata.bb.gem.m4.strip.IsV]
   Int_t           Ndata_bb_gem_m4_strip_TSchi2;
   Double_t        bb_gem_m4_strip_TSchi2[756];   //[Ndata.bb.gem.m4.strip.TSchi2]
   Int_t           Ndata_bb_gem_m4_strip_TSprob;
   Double_t        bb_gem_m4_strip_TSprob[756];   //[Ndata.bb.gem.m4.strip.TSprob]
   Int_t           Ndata_bb_gem_m4_strip_Tcorr;
   Double_t        bb_gem_m4_strip_Tcorr[756];   //[Ndata.bb.gem.m4.strip.Tcorr]
   Int_t           Ndata_bb_gem_m4_strip_Tdiff;
   Double_t        bb_gem_m4_strip_Tdiff[756];   //[Ndata.bb.gem.m4.strip.Tdiff]
   Int_t           Ndata_bb_gem_m4_strip_Tfit;
   Double_t        bb_gem_m4_strip_Tfit[756];   //[Ndata.bb.gem.m4.strip.Tfit]
   Int_t           Ndata_bb_gem_m4_strip_Tmean;
   Double_t        bb_gem_m4_strip_Tmean[756];   //[Ndata.bb.gem.m4.strip.Tmean]
   Int_t           Ndata_bb_gem_m4_strip_TmeanDeconv;
   Double_t        bb_gem_m4_strip_TmeanDeconv[756];   //[Ndata.bb.gem.m4.strip.TmeanDeconv]
   Int_t           Ndata_bb_gem_m4_strip_Tsigma;
   Double_t        bb_gem_m4_strip_Tsigma[756];   //[Ndata.bb.gem.m4.strip.Tsigma]
   Int_t           Ndata_bb_gem_m4_strip_adc_id;
   Double_t        bb_gem_m4_strip_adc_id[756];   //[Ndata.bb.gem.m4.strip.adc_id]
   Int_t           Ndata_bb_gem_m4_strip_crate;
   Double_t        bb_gem_m4_strip_crate[756];   //[Ndata.bb.gem.m4.strip.crate]
   Int_t           Ndata_bb_gem_m4_strip_isampmax;
   Double_t        bb_gem_m4_strip_isampmax[756];   //[Ndata.bb.gem.m4.strip.isampmax]
   Int_t           Ndata_bb_gem_m4_strip_isampmaxDeconv;
   Double_t        bb_gem_m4_strip_isampmaxDeconv[756];   //[Ndata.bb.gem.m4.strip.isampmaxDeconv]
   Int_t           Ndata_bb_gem_m4_strip_isampmaxDeconvCombo;
   Double_t        bb_gem_m4_strip_isampmaxDeconvCombo[756];   //[Ndata.bb.gem.m4.strip.isampmaxDeconvCombo]
   Int_t           Ndata_bb_gem_m4_strip_isnegative;
   Double_t        bb_gem_m4_strip_isnegative[756];   //[Ndata.bb.gem.m4.strip.isnegative]
   Int_t           Ndata_bb_gem_m4_strip_isnegativeU;
   Double_t        bb_gem_m4_strip_isnegativeU[756];   //[Ndata.bb.gem.m4.strip.isnegativeU]
   Int_t           Ndata_bb_gem_m4_strip_isnegativeV;
   Double_t        bb_gem_m4_strip_isnegativeV[756];   //[Ndata.bb.gem.m4.strip.isnegativeV]
   Int_t           Ndata_bb_gem_m4_strip_isnegontrack;
   Double_t        bb_gem_m4_strip_isnegontrack[756];   //[Ndata.bb.gem.m4.strip.isnegontrack]
   Int_t           Ndata_bb_gem_m4_strip_isnegontrackU;
   Double_t        bb_gem_m4_strip_isnegontrackU[756];   //[Ndata.bb.gem.m4.strip.isnegontrackU]
   Int_t           Ndata_bb_gem_m4_strip_isnegontrackV;
   Double_t        bb_gem_m4_strip_isnegontrackV[756];   //[Ndata.bb.gem.m4.strip.isnegontrackV]
   Int_t           Ndata_bb_gem_m4_strip_istrip;
   Double_t        bb_gem_m4_strip_istrip[756];   //[Ndata.bb.gem.m4.strip.istrip]
   Int_t           Ndata_bb_gem_m4_strip_itrack;
   Double_t        bb_gem_m4_strip_itrack[756];   //[Ndata.bb.gem.m4.strip.itrack]
   Int_t           Ndata_bb_gem_m4_strip_mpd;
   Double_t        bb_gem_m4_strip_mpd[756];   //[Ndata.bb.gem.m4.strip.mpd]
   Int_t           Ndata_bb_gem_m4_strip_ontrack;
   Double_t        bb_gem_m4_strip_ontrack[756];   //[Ndata.bb.gem.m4.strip.ontrack]
   Int_t           Ndata_bb_gem_m4_strip_ontrackU;
   Double_t        bb_gem_m4_strip_ontrackU[756];   //[Ndata.bb.gem.m4.strip.ontrackU]
   Int_t           Ndata_bb_gem_m4_strip_ontrackV;
   Double_t        bb_gem_m4_strip_ontrackV[756];   //[Ndata.bb.gem.m4.strip.ontrackV]
   Int_t           Ndata_bb_gem_m4_strip_rawADCsamples;
   Double_t        bb_gem_m4_strip_rawADCsamples[4536];   //[Ndata.bb.gem.m4.strip.rawADCsamples]
   Int_t           Ndata_bb_gem_m5_clust_clustu_adc;
   Double_t        bb_gem_m5_clust_clustu_adc[17];   //[Ndata.bb.gem.m5.clust.clustu_adc]
   Int_t           Ndata_bb_gem_m5_clust_clustu_pos;
   Double_t        bb_gem_m5_clust_clustu_pos[17];   //[Ndata.bb.gem.m5.clust.clustu_pos]
   Int_t           Ndata_bb_gem_m5_clust_clustu_strips;
   Double_t        bb_gem_m5_clust_clustu_strips[17];   //[Ndata.bb.gem.m5.clust.clustu_strips]
   Int_t           Ndata_bb_gem_m5_clust_clustu_time;
   Double_t        bb_gem_m5_clust_clustu_time[17];   //[Ndata.bb.gem.m5.clust.clustu_time]
   Int_t           Ndata_bb_gem_m5_clust_clustv_adc;
   Double_t        bb_gem_m5_clust_clustv_adc[15];   //[Ndata.bb.gem.m5.clust.clustv_adc]
   Int_t           Ndata_bb_gem_m5_clust_clustv_pos;
   Double_t        bb_gem_m5_clust_clustv_pos[15];   //[Ndata.bb.gem.m5.clust.clustv_pos]
   Int_t           Ndata_bb_gem_m5_clust_clustv_strips;
   Double_t        bb_gem_m5_clust_clustv_strips[15];   //[Ndata.bb.gem.m5.clust.clustv_strips]
   Int_t           Ndata_bb_gem_m5_clust_clustv_time;
   Double_t        bb_gem_m5_clust_clustv_time[15];   //[Ndata.bb.gem.m5.clust.clustv_time]
   Int_t           Ndata_bb_gem_m5_clust_isnegativeU;
   Double_t        bb_gem_m5_clust_isnegativeU[17];   //[Ndata.bb.gem.m5.clust.isnegativeU]
   Int_t           Ndata_bb_gem_m5_clust_isnegativeV;
   Double_t        bb_gem_m5_clust_isnegativeV[15];   //[Ndata.bb.gem.m5.clust.isnegativeV]
   Int_t           Ndata_bb_gem_m5_clust_isnegontrackU;
   Double_t        bb_gem_m5_clust_isnegontrackU[17];   //[Ndata.bb.gem.m5.clust.isnegontrackU]
   Int_t           Ndata_bb_gem_m5_clust_isnegontrackV;
   Double_t        bb_gem_m5_clust_isnegontrackV[15];   //[Ndata.bb.gem.m5.clust.isnegontrackV]
   Int_t           Ndata_bb_gem_m5_strip_ADCavg;
   Double_t        bb_gem_m5_strip_ADCavg[645];   //[Ndata.bb.gem.m5.strip.ADCavg]
   Int_t           Ndata_bb_gem_m5_strip_ADCmax;
   Double_t        bb_gem_m5_strip_ADCmax[645];   //[Ndata.bb.gem.m5.strip.ADCmax]
   Int_t           Ndata_bb_gem_m5_strip_ADCsamples;
   Double_t        bb_gem_m5_strip_ADCsamples[3870];   //[Ndata.bb.gem.m5.strip.ADCsamples]
   Int_t           Ndata_bb_gem_m5_strip_ADCsum;
   Double_t        bb_gem_m5_strip_ADCsum[645];   //[Ndata.bb.gem.m5.strip.ADCsum]
   Int_t           Ndata_bb_gem_m5_strip_BUILD_ALL_SAMPLES;
   Double_t        bb_gem_m5_strip_BUILD_ALL_SAMPLES[645];   //[Ndata.bb.gem.m5.strip.BUILD_ALL_SAMPLES]
   Int_t           Ndata_bb_gem_m5_strip_CM_GOOD;
   Double_t        bb_gem_m5_strip_CM_GOOD[645];   //[Ndata.bb.gem.m5.strip.CM_GOOD]
   Int_t           Ndata_bb_gem_m5_strip_CorrCoeff;
   Double_t        bb_gem_m5_strip_CorrCoeff[645];   //[Ndata.bb.gem.m5.strip.CorrCoeff]
   Int_t           Ndata_bb_gem_m5_strip_DeconvADCmax;
   Double_t        bb_gem_m5_strip_DeconvADCmax[645];   //[Ndata.bb.gem.m5.strip.DeconvADCmax]
   Int_t           Ndata_bb_gem_m5_strip_DeconvADCmaxCombo;
   Double_t        bb_gem_m5_strip_DeconvADCmaxCombo[645];   //[Ndata.bb.gem.m5.strip.DeconvADCmaxCombo]
   Int_t           Ndata_bb_gem_m5_strip_DeconvADCsamples;
   Double_t        bb_gem_m5_strip_DeconvADCsamples[3870];   //[Ndata.bb.gem.m5.strip.DeconvADCsamples]
   Int_t           Ndata_bb_gem_m5_strip_DeconvADCsum;
   Double_t        bb_gem_m5_strip_DeconvADCsum[645];   //[Ndata.bb.gem.m5.strip.DeconvADCsum]
   Int_t           Ndata_bb_gem_m5_strip_ENABLE_CM;
   Double_t        bb_gem_m5_strip_ENABLE_CM[645];   //[Ndata.bb.gem.m5.strip.ENABLE_CM]
   Int_t           Ndata_bb_gem_m5_strip_IsU;
   Double_t        bb_gem_m5_strip_IsU[645];   //[Ndata.bb.gem.m5.strip.IsU]
   Int_t           Ndata_bb_gem_m5_strip_IsV;
   Double_t        bb_gem_m5_strip_IsV[645];   //[Ndata.bb.gem.m5.strip.IsV]
   Int_t           Ndata_bb_gem_m5_strip_TSchi2;
   Double_t        bb_gem_m5_strip_TSchi2[645];   //[Ndata.bb.gem.m5.strip.TSchi2]
   Int_t           Ndata_bb_gem_m5_strip_TSprob;
   Double_t        bb_gem_m5_strip_TSprob[645];   //[Ndata.bb.gem.m5.strip.TSprob]
   Int_t           Ndata_bb_gem_m5_strip_Tcorr;
   Double_t        bb_gem_m5_strip_Tcorr[645];   //[Ndata.bb.gem.m5.strip.Tcorr]
   Int_t           Ndata_bb_gem_m5_strip_Tdiff;
   Double_t        bb_gem_m5_strip_Tdiff[645];   //[Ndata.bb.gem.m5.strip.Tdiff]
   Int_t           Ndata_bb_gem_m5_strip_Tfit;
   Double_t        bb_gem_m5_strip_Tfit[645];   //[Ndata.bb.gem.m5.strip.Tfit]
   Int_t           Ndata_bb_gem_m5_strip_Tmean;
   Double_t        bb_gem_m5_strip_Tmean[645];   //[Ndata.bb.gem.m5.strip.Tmean]
   Int_t           Ndata_bb_gem_m5_strip_TmeanDeconv;
   Double_t        bb_gem_m5_strip_TmeanDeconv[645];   //[Ndata.bb.gem.m5.strip.TmeanDeconv]
   Int_t           Ndata_bb_gem_m5_strip_Tsigma;
   Double_t        bb_gem_m5_strip_Tsigma[645];   //[Ndata.bb.gem.m5.strip.Tsigma]
   Int_t           Ndata_bb_gem_m5_strip_adc_id;
   Double_t        bb_gem_m5_strip_adc_id[645];   //[Ndata.bb.gem.m5.strip.adc_id]
   Int_t           Ndata_bb_gem_m5_strip_crate;
   Double_t        bb_gem_m5_strip_crate[645];   //[Ndata.bb.gem.m5.strip.crate]
   Int_t           Ndata_bb_gem_m5_strip_isampmax;
   Double_t        bb_gem_m5_strip_isampmax[645];   //[Ndata.bb.gem.m5.strip.isampmax]
   Int_t           Ndata_bb_gem_m5_strip_isampmaxDeconv;
   Double_t        bb_gem_m5_strip_isampmaxDeconv[645];   //[Ndata.bb.gem.m5.strip.isampmaxDeconv]
   Int_t           Ndata_bb_gem_m5_strip_isampmaxDeconvCombo;
   Double_t        bb_gem_m5_strip_isampmaxDeconvCombo[645];   //[Ndata.bb.gem.m5.strip.isampmaxDeconvCombo]
   Int_t           Ndata_bb_gem_m5_strip_isnegative;
   Double_t        bb_gem_m5_strip_isnegative[645];   //[Ndata.bb.gem.m5.strip.isnegative]
   Int_t           Ndata_bb_gem_m5_strip_isnegativeU;
   Double_t        bb_gem_m5_strip_isnegativeU[645];   //[Ndata.bb.gem.m5.strip.isnegativeU]
   Int_t           Ndata_bb_gem_m5_strip_isnegativeV;
   Double_t        bb_gem_m5_strip_isnegativeV[645];   //[Ndata.bb.gem.m5.strip.isnegativeV]
   Int_t           Ndata_bb_gem_m5_strip_isnegontrack;
   Double_t        bb_gem_m5_strip_isnegontrack[645];   //[Ndata.bb.gem.m5.strip.isnegontrack]
   Int_t           Ndata_bb_gem_m5_strip_isnegontrackU;
   Double_t        bb_gem_m5_strip_isnegontrackU[645];   //[Ndata.bb.gem.m5.strip.isnegontrackU]
   Int_t           Ndata_bb_gem_m5_strip_isnegontrackV;
   Double_t        bb_gem_m5_strip_isnegontrackV[645];   //[Ndata.bb.gem.m5.strip.isnegontrackV]
   Int_t           Ndata_bb_gem_m5_strip_istrip;
   Double_t        bb_gem_m5_strip_istrip[645];   //[Ndata.bb.gem.m5.strip.istrip]
   Int_t           Ndata_bb_gem_m5_strip_itrack;
   Double_t        bb_gem_m5_strip_itrack[645];   //[Ndata.bb.gem.m5.strip.itrack]
   Int_t           Ndata_bb_gem_m5_strip_mpd;
   Double_t        bb_gem_m5_strip_mpd[645];   //[Ndata.bb.gem.m5.strip.mpd]
   Int_t           Ndata_bb_gem_m5_strip_ontrack;
   Double_t        bb_gem_m5_strip_ontrack[645];   //[Ndata.bb.gem.m5.strip.ontrack]
   Int_t           Ndata_bb_gem_m5_strip_ontrackU;
   Double_t        bb_gem_m5_strip_ontrackU[645];   //[Ndata.bb.gem.m5.strip.ontrackU]
   Int_t           Ndata_bb_gem_m5_strip_ontrackV;
   Double_t        bb_gem_m5_strip_ontrackV[645];   //[Ndata.bb.gem.m5.strip.ontrackV]
   Int_t           Ndata_bb_gem_m5_strip_rawADCsamples;
   Double_t        bb_gem_m5_strip_rawADCsamples[3870];   //[Ndata.bb.gem.m5.strip.rawADCsamples]
   Int_t           Ndata_bb_gem_m6_clust_clustu_adc;
   Double_t        bb_gem_m6_clust_clustu_adc[6];   //[Ndata.bb.gem.m6.clust.clustu_adc]
   Int_t           Ndata_bb_gem_m6_clust_clustu_pos;
   Double_t        bb_gem_m6_clust_clustu_pos[6];   //[Ndata.bb.gem.m6.clust.clustu_pos]
   Int_t           Ndata_bb_gem_m6_clust_clustu_strips;
   Double_t        bb_gem_m6_clust_clustu_strips[6];   //[Ndata.bb.gem.m6.clust.clustu_strips]
   Int_t           Ndata_bb_gem_m6_clust_clustu_time;
   Double_t        bb_gem_m6_clust_clustu_time[6];   //[Ndata.bb.gem.m6.clust.clustu_time]
   Int_t           Ndata_bb_gem_m6_clust_clustv_adc;
   Double_t        bb_gem_m6_clust_clustv_adc[15];   //[Ndata.bb.gem.m6.clust.clustv_adc]
   Int_t           Ndata_bb_gem_m6_clust_clustv_pos;
   Double_t        bb_gem_m6_clust_clustv_pos[15];   //[Ndata.bb.gem.m6.clust.clustv_pos]
   Int_t           Ndata_bb_gem_m6_clust_clustv_strips;
   Double_t        bb_gem_m6_clust_clustv_strips[15];   //[Ndata.bb.gem.m6.clust.clustv_strips]
   Int_t           Ndata_bb_gem_m6_clust_clustv_time;
   Double_t        bb_gem_m6_clust_clustv_time[15];   //[Ndata.bb.gem.m6.clust.clustv_time]
   Int_t           Ndata_bb_gem_m6_clust_isnegativeU;
   Double_t        bb_gem_m6_clust_isnegativeU[6];   //[Ndata.bb.gem.m6.clust.isnegativeU]
   Int_t           Ndata_bb_gem_m6_clust_isnegativeV;
   Double_t        bb_gem_m6_clust_isnegativeV[15];   //[Ndata.bb.gem.m6.clust.isnegativeV]
   Int_t           Ndata_bb_gem_m6_clust_isnegontrackU;
   Double_t        bb_gem_m6_clust_isnegontrackU[6];   //[Ndata.bb.gem.m6.clust.isnegontrackU]
   Int_t           Ndata_bb_gem_m6_clust_isnegontrackV;
   Double_t        bb_gem_m6_clust_isnegontrackV[15];   //[Ndata.bb.gem.m6.clust.isnegontrackV]
   Int_t           Ndata_bb_gem_m6_strip_ADCavg;
   Double_t        bb_gem_m6_strip_ADCavg[657];   //[Ndata.bb.gem.m6.strip.ADCavg]
   Int_t           Ndata_bb_gem_m6_strip_ADCmax;
   Double_t        bb_gem_m6_strip_ADCmax[657];   //[Ndata.bb.gem.m6.strip.ADCmax]
   Int_t           Ndata_bb_gem_m6_strip_ADCsamples;
   Double_t        bb_gem_m6_strip_ADCsamples[3942];   //[Ndata.bb.gem.m6.strip.ADCsamples]
   Int_t           Ndata_bb_gem_m6_strip_ADCsum;
   Double_t        bb_gem_m6_strip_ADCsum[657];   //[Ndata.bb.gem.m6.strip.ADCsum]
   Int_t           Ndata_bb_gem_m6_strip_BUILD_ALL_SAMPLES;
   Double_t        bb_gem_m6_strip_BUILD_ALL_SAMPLES[657];   //[Ndata.bb.gem.m6.strip.BUILD_ALL_SAMPLES]
   Int_t           Ndata_bb_gem_m6_strip_CM_GOOD;
   Double_t        bb_gem_m6_strip_CM_GOOD[657];   //[Ndata.bb.gem.m6.strip.CM_GOOD]
   Int_t           Ndata_bb_gem_m6_strip_CorrCoeff;
   Double_t        bb_gem_m6_strip_CorrCoeff[657];   //[Ndata.bb.gem.m6.strip.CorrCoeff]
   Int_t           Ndata_bb_gem_m6_strip_DeconvADCmax;
   Double_t        bb_gem_m6_strip_DeconvADCmax[657];   //[Ndata.bb.gem.m6.strip.DeconvADCmax]
   Int_t           Ndata_bb_gem_m6_strip_DeconvADCmaxCombo;
   Double_t        bb_gem_m6_strip_DeconvADCmaxCombo[657];   //[Ndata.bb.gem.m6.strip.DeconvADCmaxCombo]
   Int_t           Ndata_bb_gem_m6_strip_DeconvADCsamples;
   Double_t        bb_gem_m6_strip_DeconvADCsamples[3942];   //[Ndata.bb.gem.m6.strip.DeconvADCsamples]
   Int_t           Ndata_bb_gem_m6_strip_DeconvADCsum;
   Double_t        bb_gem_m6_strip_DeconvADCsum[657];   //[Ndata.bb.gem.m6.strip.DeconvADCsum]
   Int_t           Ndata_bb_gem_m6_strip_ENABLE_CM;
   Double_t        bb_gem_m6_strip_ENABLE_CM[657];   //[Ndata.bb.gem.m6.strip.ENABLE_CM]
   Int_t           Ndata_bb_gem_m6_strip_IsU;
   Double_t        bb_gem_m6_strip_IsU[657];   //[Ndata.bb.gem.m6.strip.IsU]
   Int_t           Ndata_bb_gem_m6_strip_IsV;
   Double_t        bb_gem_m6_strip_IsV[657];   //[Ndata.bb.gem.m6.strip.IsV]
   Int_t           Ndata_bb_gem_m6_strip_TSchi2;
   Double_t        bb_gem_m6_strip_TSchi2[657];   //[Ndata.bb.gem.m6.strip.TSchi2]
   Int_t           Ndata_bb_gem_m6_strip_TSprob;
   Double_t        bb_gem_m6_strip_TSprob[657];   //[Ndata.bb.gem.m6.strip.TSprob]
   Int_t           Ndata_bb_gem_m6_strip_Tcorr;
   Double_t        bb_gem_m6_strip_Tcorr[657];   //[Ndata.bb.gem.m6.strip.Tcorr]
   Int_t           Ndata_bb_gem_m6_strip_Tdiff;
   Double_t        bb_gem_m6_strip_Tdiff[657];   //[Ndata.bb.gem.m6.strip.Tdiff]
   Int_t           Ndata_bb_gem_m6_strip_Tfit;
   Double_t        bb_gem_m6_strip_Tfit[657];   //[Ndata.bb.gem.m6.strip.Tfit]
   Int_t           Ndata_bb_gem_m6_strip_Tmean;
   Double_t        bb_gem_m6_strip_Tmean[657];   //[Ndata.bb.gem.m6.strip.Tmean]
   Int_t           Ndata_bb_gem_m6_strip_TmeanDeconv;
   Double_t        bb_gem_m6_strip_TmeanDeconv[657];   //[Ndata.bb.gem.m6.strip.TmeanDeconv]
   Int_t           Ndata_bb_gem_m6_strip_Tsigma;
   Double_t        bb_gem_m6_strip_Tsigma[657];   //[Ndata.bb.gem.m6.strip.Tsigma]
   Int_t           Ndata_bb_gem_m6_strip_adc_id;
   Double_t        bb_gem_m6_strip_adc_id[657];   //[Ndata.bb.gem.m6.strip.adc_id]
   Int_t           Ndata_bb_gem_m6_strip_crate;
   Double_t        bb_gem_m6_strip_crate[657];   //[Ndata.bb.gem.m6.strip.crate]
   Int_t           Ndata_bb_gem_m6_strip_isampmax;
   Double_t        bb_gem_m6_strip_isampmax[657];   //[Ndata.bb.gem.m6.strip.isampmax]
   Int_t           Ndata_bb_gem_m6_strip_isampmaxDeconv;
   Double_t        bb_gem_m6_strip_isampmaxDeconv[657];   //[Ndata.bb.gem.m6.strip.isampmaxDeconv]
   Int_t           Ndata_bb_gem_m6_strip_isampmaxDeconvCombo;
   Double_t        bb_gem_m6_strip_isampmaxDeconvCombo[657];   //[Ndata.bb.gem.m6.strip.isampmaxDeconvCombo]
   Int_t           Ndata_bb_gem_m6_strip_isnegative;
   Double_t        bb_gem_m6_strip_isnegative[657];   //[Ndata.bb.gem.m6.strip.isnegative]
   Int_t           Ndata_bb_gem_m6_strip_isnegativeU;
   Double_t        bb_gem_m6_strip_isnegativeU[657];   //[Ndata.bb.gem.m6.strip.isnegativeU]
   Int_t           Ndata_bb_gem_m6_strip_isnegativeV;
   Double_t        bb_gem_m6_strip_isnegativeV[657];   //[Ndata.bb.gem.m6.strip.isnegativeV]
   Int_t           Ndata_bb_gem_m6_strip_isnegontrack;
   Double_t        bb_gem_m6_strip_isnegontrack[657];   //[Ndata.bb.gem.m6.strip.isnegontrack]
   Int_t           Ndata_bb_gem_m6_strip_isnegontrackU;
   Double_t        bb_gem_m6_strip_isnegontrackU[657];   //[Ndata.bb.gem.m6.strip.isnegontrackU]
   Int_t           Ndata_bb_gem_m6_strip_isnegontrackV;
   Double_t        bb_gem_m6_strip_isnegontrackV[657];   //[Ndata.bb.gem.m6.strip.isnegontrackV]
   Int_t           Ndata_bb_gem_m6_strip_istrip;
   Double_t        bb_gem_m6_strip_istrip[657];   //[Ndata.bb.gem.m6.strip.istrip]
   Int_t           Ndata_bb_gem_m6_strip_itrack;
   Double_t        bb_gem_m6_strip_itrack[657];   //[Ndata.bb.gem.m6.strip.itrack]
   Int_t           Ndata_bb_gem_m6_strip_mpd;
   Double_t        bb_gem_m6_strip_mpd[657];   //[Ndata.bb.gem.m6.strip.mpd]
   Int_t           Ndata_bb_gem_m6_strip_ontrack;
   Double_t        bb_gem_m6_strip_ontrack[657];   //[Ndata.bb.gem.m6.strip.ontrack]
   Int_t           Ndata_bb_gem_m6_strip_ontrackU;
   Double_t        bb_gem_m6_strip_ontrackU[657];   //[Ndata.bb.gem.m6.strip.ontrackU]
   Int_t           Ndata_bb_gem_m6_strip_ontrackV;
   Double_t        bb_gem_m6_strip_ontrackV[657];   //[Ndata.bb.gem.m6.strip.ontrackV]
   Int_t           Ndata_bb_gem_m6_strip_rawADCsamples;
   Double_t        bb_gem_m6_strip_rawADCsamples[3942];   //[Ndata.bb.gem.m6.strip.rawADCsamples]
   Int_t           Ndata_bb_gem_m7_clust_clustu_adc;
   Double_t        bb_gem_m7_clust_clustu_adc[11];   //[Ndata.bb.gem.m7.clust.clustu_adc]
   Int_t           Ndata_bb_gem_m7_clust_clustu_pos;
   Double_t        bb_gem_m7_clust_clustu_pos[11];   //[Ndata.bb.gem.m7.clust.clustu_pos]
   Int_t           Ndata_bb_gem_m7_clust_clustu_strips;
   Double_t        bb_gem_m7_clust_clustu_strips[11];   //[Ndata.bb.gem.m7.clust.clustu_strips]
   Int_t           Ndata_bb_gem_m7_clust_clustu_time;
   Double_t        bb_gem_m7_clust_clustu_time[11];   //[Ndata.bb.gem.m7.clust.clustu_time]
   Int_t           Ndata_bb_gem_m7_clust_clustv_adc;
   Double_t        bb_gem_m7_clust_clustv_adc[16];   //[Ndata.bb.gem.m7.clust.clustv_adc]
   Int_t           Ndata_bb_gem_m7_clust_clustv_pos;
   Double_t        bb_gem_m7_clust_clustv_pos[16];   //[Ndata.bb.gem.m7.clust.clustv_pos]
   Int_t           Ndata_bb_gem_m7_clust_clustv_strips;
   Double_t        bb_gem_m7_clust_clustv_strips[16];   //[Ndata.bb.gem.m7.clust.clustv_strips]
   Int_t           Ndata_bb_gem_m7_clust_clustv_time;
   Double_t        bb_gem_m7_clust_clustv_time[16];   //[Ndata.bb.gem.m7.clust.clustv_time]
   Int_t           Ndata_bb_gem_m7_clust_isnegativeU;
   Double_t        bb_gem_m7_clust_isnegativeU[11];   //[Ndata.bb.gem.m7.clust.isnegativeU]
   Int_t           Ndata_bb_gem_m7_clust_isnegativeV;
   Double_t        bb_gem_m7_clust_isnegativeV[16];   //[Ndata.bb.gem.m7.clust.isnegativeV]
   Int_t           Ndata_bb_gem_m7_clust_isnegontrackU;
   Double_t        bb_gem_m7_clust_isnegontrackU[11];   //[Ndata.bb.gem.m7.clust.isnegontrackU]
   Int_t           Ndata_bb_gem_m7_clust_isnegontrackV;
   Double_t        bb_gem_m7_clust_isnegontrackV[16];   //[Ndata.bb.gem.m7.clust.isnegontrackV]
   Int_t           Ndata_bb_gem_m7_strip_ADCavg;
   Double_t        bb_gem_m7_strip_ADCavg[632];   //[Ndata.bb.gem.m7.strip.ADCavg]
   Int_t           Ndata_bb_gem_m7_strip_ADCmax;
   Double_t        bb_gem_m7_strip_ADCmax[632];   //[Ndata.bb.gem.m7.strip.ADCmax]
   Int_t           Ndata_bb_gem_m7_strip_ADCsamples;
   Double_t        bb_gem_m7_strip_ADCsamples[3792];   //[Ndata.bb.gem.m7.strip.ADCsamples]
   Int_t           Ndata_bb_gem_m7_strip_ADCsum;
   Double_t        bb_gem_m7_strip_ADCsum[632];   //[Ndata.bb.gem.m7.strip.ADCsum]
   Int_t           Ndata_bb_gem_m7_strip_BUILD_ALL_SAMPLES;
   Double_t        bb_gem_m7_strip_BUILD_ALL_SAMPLES[632];   //[Ndata.bb.gem.m7.strip.BUILD_ALL_SAMPLES]
   Int_t           Ndata_bb_gem_m7_strip_CM_GOOD;
   Double_t        bb_gem_m7_strip_CM_GOOD[632];   //[Ndata.bb.gem.m7.strip.CM_GOOD]
   Int_t           Ndata_bb_gem_m7_strip_CorrCoeff;
   Double_t        bb_gem_m7_strip_CorrCoeff[632];   //[Ndata.bb.gem.m7.strip.CorrCoeff]
   Int_t           Ndata_bb_gem_m7_strip_DeconvADCmax;
   Double_t        bb_gem_m7_strip_DeconvADCmax[632];   //[Ndata.bb.gem.m7.strip.DeconvADCmax]
   Int_t           Ndata_bb_gem_m7_strip_DeconvADCmaxCombo;
   Double_t        bb_gem_m7_strip_DeconvADCmaxCombo[632];   //[Ndata.bb.gem.m7.strip.DeconvADCmaxCombo]
   Int_t           Ndata_bb_gem_m7_strip_DeconvADCsamples;
   Double_t        bb_gem_m7_strip_DeconvADCsamples[3792];   //[Ndata.bb.gem.m7.strip.DeconvADCsamples]
   Int_t           Ndata_bb_gem_m7_strip_DeconvADCsum;
   Double_t        bb_gem_m7_strip_DeconvADCsum[632];   //[Ndata.bb.gem.m7.strip.DeconvADCsum]
   Int_t           Ndata_bb_gem_m7_strip_ENABLE_CM;
   Double_t        bb_gem_m7_strip_ENABLE_CM[632];   //[Ndata.bb.gem.m7.strip.ENABLE_CM]
   Int_t           Ndata_bb_gem_m7_strip_IsU;
   Double_t        bb_gem_m7_strip_IsU[632];   //[Ndata.bb.gem.m7.strip.IsU]
   Int_t           Ndata_bb_gem_m7_strip_IsV;
   Double_t        bb_gem_m7_strip_IsV[632];   //[Ndata.bb.gem.m7.strip.IsV]
   Int_t           Ndata_bb_gem_m7_strip_TSchi2;
   Double_t        bb_gem_m7_strip_TSchi2[632];   //[Ndata.bb.gem.m7.strip.TSchi2]
   Int_t           Ndata_bb_gem_m7_strip_TSprob;
   Double_t        bb_gem_m7_strip_TSprob[632];   //[Ndata.bb.gem.m7.strip.TSprob]
   Int_t           Ndata_bb_gem_m7_strip_Tcorr;
   Double_t        bb_gem_m7_strip_Tcorr[632];   //[Ndata.bb.gem.m7.strip.Tcorr]
   Int_t           Ndata_bb_gem_m7_strip_Tdiff;
   Double_t        bb_gem_m7_strip_Tdiff[632];   //[Ndata.bb.gem.m7.strip.Tdiff]
   Int_t           Ndata_bb_gem_m7_strip_Tfit;
   Double_t        bb_gem_m7_strip_Tfit[632];   //[Ndata.bb.gem.m7.strip.Tfit]
   Int_t           Ndata_bb_gem_m7_strip_Tmean;
   Double_t        bb_gem_m7_strip_Tmean[632];   //[Ndata.bb.gem.m7.strip.Tmean]
   Int_t           Ndata_bb_gem_m7_strip_TmeanDeconv;
   Double_t        bb_gem_m7_strip_TmeanDeconv[632];   //[Ndata.bb.gem.m7.strip.TmeanDeconv]
   Int_t           Ndata_bb_gem_m7_strip_Tsigma;
   Double_t        bb_gem_m7_strip_Tsigma[632];   //[Ndata.bb.gem.m7.strip.Tsigma]
   Int_t           Ndata_bb_gem_m7_strip_adc_id;
   Double_t        bb_gem_m7_strip_adc_id[632];   //[Ndata.bb.gem.m7.strip.adc_id]
   Int_t           Ndata_bb_gem_m7_strip_crate;
   Double_t        bb_gem_m7_strip_crate[632];   //[Ndata.bb.gem.m7.strip.crate]
   Int_t           Ndata_bb_gem_m7_strip_isampmax;
   Double_t        bb_gem_m7_strip_isampmax[632];   //[Ndata.bb.gem.m7.strip.isampmax]
   Int_t           Ndata_bb_gem_m7_strip_isampmaxDeconv;
   Double_t        bb_gem_m7_strip_isampmaxDeconv[632];   //[Ndata.bb.gem.m7.strip.isampmaxDeconv]
   Int_t           Ndata_bb_gem_m7_strip_isampmaxDeconvCombo;
   Double_t        bb_gem_m7_strip_isampmaxDeconvCombo[632];   //[Ndata.bb.gem.m7.strip.isampmaxDeconvCombo]
   Int_t           Ndata_bb_gem_m7_strip_isnegative;
   Double_t        bb_gem_m7_strip_isnegative[632];   //[Ndata.bb.gem.m7.strip.isnegative]
   Int_t           Ndata_bb_gem_m7_strip_isnegativeU;
   Double_t        bb_gem_m7_strip_isnegativeU[632];   //[Ndata.bb.gem.m7.strip.isnegativeU]
   Int_t           Ndata_bb_gem_m7_strip_isnegativeV;
   Double_t        bb_gem_m7_strip_isnegativeV[632];   //[Ndata.bb.gem.m7.strip.isnegativeV]
   Int_t           Ndata_bb_gem_m7_strip_isnegontrack;
   Double_t        bb_gem_m7_strip_isnegontrack[632];   //[Ndata.bb.gem.m7.strip.isnegontrack]
   Int_t           Ndata_bb_gem_m7_strip_isnegontrackU;
   Double_t        bb_gem_m7_strip_isnegontrackU[632];   //[Ndata.bb.gem.m7.strip.isnegontrackU]
   Int_t           Ndata_bb_gem_m7_strip_isnegontrackV;
   Double_t        bb_gem_m7_strip_isnegontrackV[632];   //[Ndata.bb.gem.m7.strip.isnegontrackV]
   Int_t           Ndata_bb_gem_m7_strip_istrip;
   Double_t        bb_gem_m7_strip_istrip[632];   //[Ndata.bb.gem.m7.strip.istrip]
   Int_t           Ndata_bb_gem_m7_strip_itrack;
   Double_t        bb_gem_m7_strip_itrack[632];   //[Ndata.bb.gem.m7.strip.itrack]
   Int_t           Ndata_bb_gem_m7_strip_mpd;
   Double_t        bb_gem_m7_strip_mpd[632];   //[Ndata.bb.gem.m7.strip.mpd]
   Int_t           Ndata_bb_gem_m7_strip_ontrack;
   Double_t        bb_gem_m7_strip_ontrack[632];   //[Ndata.bb.gem.m7.strip.ontrack]
   Int_t           Ndata_bb_gem_m7_strip_ontrackU;
   Double_t        bb_gem_m7_strip_ontrackU[632];   //[Ndata.bb.gem.m7.strip.ontrackU]
   Int_t           Ndata_bb_gem_m7_strip_ontrackV;
   Double_t        bb_gem_m7_strip_ontrackV[632];   //[Ndata.bb.gem.m7.strip.ontrackV]
   Int_t           Ndata_bb_gem_m7_strip_rawADCsamples;
   Double_t        bb_gem_m7_strip_rawADCsamples[3792];   //[Ndata.bb.gem.m7.strip.rawADCsamples]
   Int_t           Ndata_bb_gem_n2Dhit_layer;
   Double_t        bb_gem_n2Dhit_layer[5];   //[Ndata.bb.gem.n2Dhit_layer]
   Int_t           Ndata_bb_gem_nclustu_layer;
   Double_t        bb_gem_nclustu_layer[5];   //[Ndata.bb.gem.nclustu_layer]
   Int_t           Ndata_bb_gem_nclustv_layer;
   Double_t        bb_gem_nclustv_layer[5];   //[Ndata.bb.gem.nclustv_layer]
   Int_t           Ndata_bb_gem_nstripsu_layer;
   Double_t        bb_gem_nstripsu_layer[5];   //[Ndata.bb.gem.nstripsu_layer]
   Int_t           Ndata_bb_gem_nstripsv_layer;
   Double_t        bb_gem_nstripsv_layer[5];   //[Ndata.bb.gem.nstripsv_layer]
   Int_t           Ndata_bb_gem_track_chi2ndf;
   Double_t        bb_gem_track_chi2ndf[1];   //[Ndata.bb.gem.track.chi2ndf]
   Int_t           Ndata_bb_gem_track_chi2ndf_hitquality;
   Double_t        bb_gem_track_chi2ndf_hitquality[1];   //[Ndata.bb.gem.track.chi2ndf_hitquality]
   Int_t           Ndata_bb_gem_track_ngoodhits;
   Double_t        bb_gem_track_ngoodhits[1];   //[Ndata.bb.gem.track.ngoodhits]
   Int_t           Ndata_bb_gem_track_nhits;
   Double_t        bb_gem_track_nhits[1];   //[Ndata.bb.gem.track.nhits]
   Int_t           Ndata_bb_gem_track_t0;
   Double_t        bb_gem_track_t0[1];   //[Ndata.bb.gem.track.t0]
   Int_t           Ndata_bb_gem_track_x;
   Double_t        bb_gem_track_x[1];   //[Ndata.bb.gem.track.x]
   Int_t           Ndata_bb_gem_track_xp;
   Double_t        bb_gem_track_xp[1];   //[Ndata.bb.gem.track.xp]
   Int_t           Ndata_bb_gem_track_y;
   Double_t        bb_gem_track_y[1];   //[Ndata.bb.gem.track.y]
   Int_t           Ndata_bb_gem_track_yp;
   Double_t        bb_gem_track_yp[1];   //[Ndata.bb.gem.track.yp]
   Int_t           Ndata_bb_grinch_tdc_allclus_adc;
   Double_t        bb_grinch_tdc_allclus_adc[14];   //[Ndata.bb.grinch_tdc.allclus.adc]
   Int_t           Ndata_bb_grinch_tdc_allclus_dx;
   Double_t        bb_grinch_tdc_allclus_dx[14];   //[Ndata.bb.grinch_tdc.allclus.dx]
   Int_t           Ndata_bb_grinch_tdc_allclus_dy;
   Double_t        bb_grinch_tdc_allclus_dy[14];   //[Ndata.bb.grinch_tdc.allclus.dy]
   Int_t           Ndata_bb_grinch_tdc_allclus_mirrorindex;
   Double_t        bb_grinch_tdc_allclus_mirrorindex[14];   //[Ndata.bb.grinch_tdc.allclus.mirrorindex]
   Int_t           Ndata_bb_grinch_tdc_allclus_size;
   Double_t        bb_grinch_tdc_allclus_size[14];   //[Ndata.bb.grinch_tdc.allclus.size]
   Int_t           Ndata_bb_grinch_tdc_allclus_t_mean;
   Double_t        bb_grinch_tdc_allclus_t_mean[14];   //[Ndata.bb.grinch_tdc.allclus.t_mean]
   Int_t           Ndata_bb_grinch_tdc_allclus_t_rms;
   Double_t        bb_grinch_tdc_allclus_t_rms[14];   //[Ndata.bb.grinch_tdc.allclus.t_rms]
   Int_t           Ndata_bb_grinch_tdc_allclus_tot_mean;
   Double_t        bb_grinch_tdc_allclus_tot_mean[14];   //[Ndata.bb.grinch_tdc.allclus.tot_mean]
   Int_t           Ndata_bb_grinch_tdc_allclus_trackindex;
   Double_t        bb_grinch_tdc_allclus_trackindex[14];   //[Ndata.bb.grinch_tdc.allclus.trackindex]
   Int_t           Ndata_bb_grinch_tdc_allclus_x_mean;
   Double_t        bb_grinch_tdc_allclus_x_mean[14];   //[Ndata.bb.grinch_tdc.allclus.x_mean]
   Int_t           Ndata_bb_grinch_tdc_allclus_y_mean;
   Double_t        bb_grinch_tdc_allclus_y_mean[14];   //[Ndata.bb.grinch_tdc.allclus.y_mean]
   Int_t           Ndata_bb_grinch_tdc_hit_amp;
   Double_t        bb_grinch_tdc_hit_amp[28];   //[Ndata.bb.grinch_tdc.hit.amp]
   Int_t           Ndata_bb_grinch_tdc_hit_clustindex;
   Double_t        bb_grinch_tdc_hit_clustindex[28];   //[Ndata.bb.grinch_tdc.hit.clustindex]
   Int_t           Ndata_bb_grinch_tdc_hit_col;
   Double_t        bb_grinch_tdc_hit_col[28];   //[Ndata.bb.grinch_tdc.hit.col]
   Int_t           Ndata_bb_grinch_tdc_hit_pmtnum;
   Double_t        bb_grinch_tdc_hit_pmtnum[28];   //[Ndata.bb.grinch_tdc.hit.pmtnum]
   Int_t           Ndata_bb_grinch_tdc_hit_row;
   Double_t        bb_grinch_tdc_hit_row[28];   //[Ndata.bb.grinch_tdc.hit.row]
   Int_t           Ndata_bb_grinch_tdc_hit_time;
   Double_t        bb_grinch_tdc_hit_time[28];   //[Ndata.bb.grinch_tdc.hit.time]
   Int_t           Ndata_bb_grinch_tdc_hit_trackindex;
   Double_t        bb_grinch_tdc_hit_trackindex[28];   //[Ndata.bb.grinch_tdc.hit.trackindex]
   Int_t           Ndata_bb_grinch_tdc_hit_xhit;
   Double_t        bb_grinch_tdc_hit_xhit[28];   //[Ndata.bb.grinch_tdc.hit.xhit]
   Int_t           Ndata_bb_grinch_tdc_hit_yhit;
   Double_t        bb_grinch_tdc_hit_yhit[28];   //[Ndata.bb.grinch_tdc.hit.yhit]
   Int_t           Ndata_bb_hodotdc_Ref_tdc;
   Double_t        bb_hodotdc_Ref_tdc[2];   //[Ndata.bb.hodotdc.Ref.tdc]
   Int_t           Ndata_bb_hodotdc_Ref_tdc_mult;
   Double_t        bb_hodotdc_Ref_tdc_mult[2];   //[Ndata.bb.hodotdc.Ref.tdc_mult]
   Int_t           Ndata_bb_hodotdc_Ref_tdc_te;
   Double_t        bb_hodotdc_Ref_tdc_te[2];   //[Ndata.bb.hodotdc.Ref.tdc_te]
   Int_t           Ndata_bb_hodotdc_Ref_tdc_tot;
   Double_t        bb_hodotdc_Ref_tdc_tot[2];   //[Ndata.bb.hodotdc.Ref.tdc_tot]
   Int_t           Ndata_bb_hodotdc_Ref_tdcelemID;
   Double_t        bb_hodotdc_Ref_tdcelemID[2];   //[Ndata.bb.hodotdc.Ref.tdcelemID]
   Int_t           Ndata_bb_hodotdc_clus_bar_tdc_id;
   Double_t        bb_hodotdc_clus_bar_tdc_id[4];   //[Ndata.bb.hodotdc.clus.bar.tdc.id]
   Int_t           Ndata_bb_hodotdc_clus_bar_tdc_itrack;
   Double_t        bb_hodotdc_clus_bar_tdc_itrack[4];   //[Ndata.bb.hodotdc.clus.bar.tdc.itrack]
   Int_t           Ndata_bb_hodotdc_clus_bar_tdc_meantime;
   Double_t        bb_hodotdc_clus_bar_tdc_meantime[4];   //[Ndata.bb.hodotdc.clus.bar.tdc.meantime]
   Int_t           Ndata_bb_hodotdc_clus_bar_tdc_meantot;
   Double_t        bb_hodotdc_clus_bar_tdc_meantot[4];   //[Ndata.bb.hodotdc.clus.bar.tdc.meantot]
   Int_t           Ndata_bb_hodotdc_clus_bar_tdc_timediff;
   Double_t        bb_hodotdc_clus_bar_tdc_timediff[4];   //[Ndata.bb.hodotdc.clus.bar.tdc.timediff]
   Int_t           Ndata_bb_hodotdc_clus_bar_tdc_timehitpos;
   Double_t        bb_hodotdc_clus_bar_tdc_timehitpos[4];   //[Ndata.bb.hodotdc.clus.bar.tdc.timehitpos]
   Int_t           Ndata_bb_hodotdc_clus_bar_tdc_tleft;
   Double_t        bb_hodotdc_clus_bar_tdc_tleft[4];   //[Ndata.bb.hodotdc.clus.bar.tdc.tleft]
   Int_t           Ndata_bb_hodotdc_clus_bar_tdc_totleft;
   Double_t        bb_hodotdc_clus_bar_tdc_totleft[4];   //[Ndata.bb.hodotdc.clus.bar.tdc.totleft]
   Int_t           Ndata_bb_hodotdc_clus_bar_tdc_totright;
   Double_t        bb_hodotdc_clus_bar_tdc_totright[4];   //[Ndata.bb.hodotdc.clus.bar.tdc.totright]
   Int_t           Ndata_bb_hodotdc_clus_bar_tdc_tright;
   Double_t        bb_hodotdc_clus_bar_tdc_tright[4];   //[Ndata.bb.hodotdc.clus.bar.tdc.tright]
   Int_t           Ndata_bb_hodotdc_clus_bar_tdc_vpos;
   Double_t        bb_hodotdc_clus_bar_tdc_vpos[4];   //[Ndata.bb.hodotdc.clus.bar.tdc.vpos]
   Int_t           Ndata_bb_hodotdc_clus_id;
   Double_t        bb_hodotdc_clus_id[1];   //[Ndata.bb.hodotdc.clus.id]
   Int_t           Ndata_bb_hodotdc_clus_size;
   Double_t        bb_hodotdc_clus_size[1];   //[Ndata.bb.hodotdc.clus.size]
   Int_t           Ndata_bb_hodotdc_clus_tdiff;
   Double_t        bb_hodotdc_clus_tdiff[1];   //[Ndata.bb.hodotdc.clus.tdiff]
   Int_t           Ndata_bb_hodotdc_clus_tleft;
   Double_t        bb_hodotdc_clus_tleft[1];   //[Ndata.bb.hodotdc.clus.tleft]
   Int_t           Ndata_bb_hodotdc_clus_tmean;
   Double_t        bb_hodotdc_clus_tmean[1];   //[Ndata.bb.hodotdc.clus.tmean]
   Int_t           Ndata_bb_hodotdc_clus_totleft;
   Double_t        bb_hodotdc_clus_totleft[1];   //[Ndata.bb.hodotdc.clus.totleft]
   Int_t           Ndata_bb_hodotdc_clus_totmean;
   Double_t        bb_hodotdc_clus_totmean[1];   //[Ndata.bb.hodotdc.clus.totmean]
   Int_t           Ndata_bb_hodotdc_clus_totright;
   Double_t        bb_hodotdc_clus_totright[1];   //[Ndata.bb.hodotdc.clus.totright]
   Int_t           Ndata_bb_hodotdc_clus_trackindex;
   Double_t        bb_hodotdc_clus_trackindex[1];   //[Ndata.bb.hodotdc.clus.trackindex]
   Int_t           Ndata_bb_hodotdc_clus_tright;
   Double_t        bb_hodotdc_clus_tright[1];   //[Ndata.bb.hodotdc.clus.tright]
   Int_t           Ndata_bb_hodotdc_clus_xmean;
   Double_t        bb_hodotdc_clus_xmean[1];   //[Ndata.bb.hodotdc.clus.xmean]
   Int_t           Ndata_bb_hodotdc_clus_ymean;
   Double_t        bb_hodotdc_clus_ymean[1];   //[Ndata.bb.hodotdc.clus.ymean]
   Int_t           Ndata_bb_ps_clus_again;
   Double_t        bb_ps_clus_again[3];   //[Ndata.bb.ps.clus.again]
   Int_t           Ndata_bb_ps_clus_atime;
   Double_t        bb_ps_clus_atime[3];   //[Ndata.bb.ps.clus.atime]
   Int_t           Ndata_bb_ps_clus_col;
   Double_t        bb_ps_clus_col[3];   //[Ndata.bb.ps.clus.col]
   Int_t           Ndata_bb_ps_clus_e;
   Double_t        bb_ps_clus_e[3];   //[Ndata.bb.ps.clus.e]
   Int_t           Ndata_bb_ps_clus_eblk;
   Double_t        bb_ps_clus_eblk[3];   //[Ndata.bb.ps.clus.eblk]
   Int_t           Ndata_bb_ps_clus_id;
   Double_t        bb_ps_clus_id[3];   //[Ndata.bb.ps.clus.id]
   Int_t           Ndata_bb_ps_clus_nblk;
   Double_t        bb_ps_clus_nblk[3];   //[Ndata.bb.ps.clus.nblk]
   Int_t           Ndata_bb_ps_clus_row;
   Double_t        bb_ps_clus_row[3];   //[Ndata.bb.ps.clus.row]
   Int_t           Ndata_bb_ps_clus_tdctime;
   Double_t        bb_ps_clus_tdctime[3];   //[Ndata.bb.ps.clus.tdctime]
   Int_t           Ndata_bb_ps_clus_x;
   Double_t        bb_ps_clus_x[3];   //[Ndata.bb.ps.clus.x]
   Int_t           Ndata_bb_ps_clus_y;
   Double_t        bb_ps_clus_y[3];   //[Ndata.bb.ps.clus.y]
   Int_t           Ndata_bb_ps_clus_blk_again;
   Double_t        bb_ps_clus_blk_again[3];   //[Ndata.bb.ps.clus_blk.again]
   Int_t           Ndata_bb_ps_clus_blk_atime;
   Double_t        bb_ps_clus_blk_atime[3];   //[Ndata.bb.ps.clus_blk.atime]
   Int_t           Ndata_bb_ps_clus_blk_col;
   Double_t        bb_ps_clus_blk_col[3];   //[Ndata.bb.ps.clus_blk.col]
   Int_t           Ndata_bb_ps_clus_blk_e;
   Double_t        bb_ps_clus_blk_e[3];   //[Ndata.bb.ps.clus_blk.e]
   Int_t           Ndata_bb_ps_clus_blk_id;
   Double_t        bb_ps_clus_blk_id[3];   //[Ndata.bb.ps.clus_blk.id]
   Int_t           Ndata_bb_ps_clus_blk_row;
   Double_t        bb_ps_clus_blk_row[3];   //[Ndata.bb.ps.clus_blk.row]
   Int_t           Ndata_bb_ps_clus_blk_tdctime;
   Double_t        bb_ps_clus_blk_tdctime[3];   //[Ndata.bb.ps.clus_blk.tdctime]
   Int_t           Ndata_bb_ps_clus_blk_x;
   Double_t        bb_ps_clus_blk_x[3];   //[Ndata.bb.ps.clus_blk.x]
   Int_t           Ndata_bb_ps_clus_blk_y;
   Double_t        bb_ps_clus_blk_y[3];   //[Ndata.bb.ps.clus_blk.y]
   Int_t           Ndata_bb_sh_clus_again;
   Double_t        bb_sh_clus_again[3];   //[Ndata.bb.sh.clus.again]
   Int_t           Ndata_bb_sh_clus_atime;
   Double_t        bb_sh_clus_atime[3];   //[Ndata.bb.sh.clus.atime]
   Int_t           Ndata_bb_sh_clus_col;
   Double_t        bb_sh_clus_col[3];   //[Ndata.bb.sh.clus.col]
   Int_t           Ndata_bb_sh_clus_e;
   Double_t        bb_sh_clus_e[3];   //[Ndata.bb.sh.clus.e]
   Int_t           Ndata_bb_sh_clus_eblk;
   Double_t        bb_sh_clus_eblk[3];   //[Ndata.bb.sh.clus.eblk]
   Int_t           Ndata_bb_sh_clus_id;
   Double_t        bb_sh_clus_id[3];   //[Ndata.bb.sh.clus.id]
   Int_t           Ndata_bb_sh_clus_nblk;
   Double_t        bb_sh_clus_nblk[3];   //[Ndata.bb.sh.clus.nblk]
   Int_t           Ndata_bb_sh_clus_row;
   Double_t        bb_sh_clus_row[3];   //[Ndata.bb.sh.clus.row]
   Int_t           Ndata_bb_sh_clus_tdctime;
   Double_t        bb_sh_clus_tdctime[3];   //[Ndata.bb.sh.clus.tdctime]
   Int_t           Ndata_bb_sh_clus_x;
   Double_t        bb_sh_clus_x[3];   //[Ndata.bb.sh.clus.x]
   Int_t           Ndata_bb_sh_clus_y;
   Double_t        bb_sh_clus_y[3];   //[Ndata.bb.sh.clus.y]
   Int_t           Ndata_bb_sh_clus_blk_again;
   Double_t        bb_sh_clus_blk_again[9];   //[Ndata.bb.sh.clus_blk.again]
   Int_t           Ndata_bb_sh_clus_blk_atime;
   Double_t        bb_sh_clus_blk_atime[9];   //[Ndata.bb.sh.clus_blk.atime]
   Int_t           Ndata_bb_sh_clus_blk_col;
   Double_t        bb_sh_clus_blk_col[9];   //[Ndata.bb.sh.clus_blk.col]
   Int_t           Ndata_bb_sh_clus_blk_e;
   Double_t        bb_sh_clus_blk_e[9];   //[Ndata.bb.sh.clus_blk.e]
   Int_t           Ndata_bb_sh_clus_blk_id;
   Double_t        bb_sh_clus_blk_id[9];   //[Ndata.bb.sh.clus_blk.id]
   Int_t           Ndata_bb_sh_clus_blk_row;
   Double_t        bb_sh_clus_blk_row[9];   //[Ndata.bb.sh.clus_blk.row]
   Int_t           Ndata_bb_sh_clus_blk_tdctime;
   Double_t        bb_sh_clus_blk_tdctime[9];   //[Ndata.bb.sh.clus_blk.tdctime]
   Int_t           Ndata_bb_sh_clus_blk_x;
   Double_t        bb_sh_clus_blk_x[9];   //[Ndata.bb.sh.clus_blk.x]
   Int_t           Ndata_bb_sh_clus_blk_y;
   Double_t        bb_sh_clus_blk_y[9];   //[Ndata.bb.sh.clus_blk.y]
   Int_t           Ndata_bb_tdctrig_tdc;
   Double_t        bb_tdctrig_tdc[6];   //[Ndata.bb.tdctrig.tdc]
   Int_t           Ndata_bb_tdctrig_tdcelemID;
   Double_t        bb_tdctrig_tdcelemID[6];   //[Ndata.bb.tdctrig.tdcelemID]
   Int_t           Ndata_bb_tr_beta;
   Double_t        bb_tr_beta[1];   //[Ndata.bb.tr.beta]
   Int_t           Ndata_bb_tr_chi2;
   Double_t        bb_tr_chi2[1];   //[Ndata.bb.tr.chi2]
   Int_t           Ndata_bb_tr_d_ph;
   Double_t        bb_tr_d_ph[1];   //[Ndata.bb.tr.d_ph]
   Int_t           Ndata_bb_tr_d_th;
   Double_t        bb_tr_d_th[1];   //[Ndata.bb.tr.d_th]
   Int_t           Ndata_bb_tr_d_x;
   Double_t        bb_tr_d_x[1];   //[Ndata.bb.tr.d_x]
   Int_t           Ndata_bb_tr_d_y;
   Double_t        bb_tr_d_y[1];   //[Ndata.bb.tr.d_y]
   Int_t           Ndata_bb_tr_dbeta;
   Double_t        bb_tr_dbeta[1];   //[Ndata.bb.tr.dbeta]
   Int_t           Ndata_bb_tr_dtime;
   Double_t        bb_tr_dtime[1];   //[Ndata.bb.tr.dtime]
   Int_t           Ndata_bb_tr_flag;
   Double_t        bb_tr_flag[1];   //[Ndata.bb.tr.flag]
   Int_t           Ndata_bb_tr_ndof;
   Double_t        bb_tr_ndof[1];   //[Ndata.bb.tr.ndof]
   Int_t           Ndata_bb_tr_p;
   Double_t        bb_tr_p[1];   //[Ndata.bb.tr.p]
   Int_t           Ndata_bb_tr_pathl;
   Double_t        bb_tr_pathl[1];   //[Ndata.bb.tr.pathl]
   Int_t           Ndata_bb_tr_ph;
   Double_t        bb_tr_ph[1];   //[Ndata.bb.tr.ph]
   Int_t           Ndata_bb_tr_px;
   Double_t        bb_tr_px[1];   //[Ndata.bb.tr.px]
   Int_t           Ndata_bb_tr_py;
   Double_t        bb_tr_py[1];   //[Ndata.bb.tr.py]
   Int_t           Ndata_bb_tr_pz;
   Double_t        bb_tr_pz[1];   //[Ndata.bb.tr.pz]
   Int_t           Ndata_bb_tr_r_ph;
   Double_t        bb_tr_r_ph[1];   //[Ndata.bb.tr.r_ph]
   Int_t           Ndata_bb_tr_r_th;
   Double_t        bb_tr_r_th[1];   //[Ndata.bb.tr.r_th]
   Int_t           Ndata_bb_tr_r_x;
   Double_t        bb_tr_r_x[1];   //[Ndata.bb.tr.r_x]
   Int_t           Ndata_bb_tr_r_y;
   Double_t        bb_tr_r_y[1];   //[Ndata.bb.tr.r_y]
   Int_t           Ndata_bb_tr_tg_dp;
   Double_t        bb_tr_tg_dp[1];   //[Ndata.bb.tr.tg_dp]
   Int_t           Ndata_bb_tr_tg_ph;
   Double_t        bb_tr_tg_ph[1];   //[Ndata.bb.tr.tg_ph]
   Int_t           Ndata_bb_tr_tg_th;
   Double_t        bb_tr_tg_th[1];   //[Ndata.bb.tr.tg_th]
   Int_t           Ndata_bb_tr_tg_x;
   Double_t        bb_tr_tg_x[1];   //[Ndata.bb.tr.tg_x]
   Int_t           Ndata_bb_tr_tg_y;
   Double_t        bb_tr_tg_y[1];   //[Ndata.bb.tr.tg_y]
   Int_t           Ndata_bb_tr_th;
   Double_t        bb_tr_th[1];   //[Ndata.bb.tr.th]
   Int_t           Ndata_bb_tr_time;
   Double_t        bb_tr_time[1];   //[Ndata.bb.tr.time]
   Int_t           Ndata_bb_tr_vx;
   Double_t        bb_tr_vx[1];   //[Ndata.bb.tr.vx]
   Int_t           Ndata_bb_tr_vy;
   Double_t        bb_tr_vy[1];   //[Ndata.bb.tr.vy]
   Int_t           Ndata_bb_tr_vz;
   Double_t        bb_tr_vz[1];   //[Ndata.bb.tr.vz]
   Int_t           Ndata_bb_tr_x;
   Double_t        bb_tr_x[1];   //[Ndata.bb.tr.x]
   Int_t           Ndata_bb_tr_y;
   Double_t        bb_tr_y[1];   //[Ndata.bb.tr.y]
   Int_t           Ndata_bb_x_bcp;
   Double_t        bb_x_bcp[1];   //[Ndata.bb.x_bcp]
   Int_t           Ndata_bb_x_fcp;
   Double_t        bb_x_fcp[1];   //[Ndata.bb.x_fcp]
   Int_t           Ndata_bb_y_bcp;
   Double_t        bb_y_bcp[1];   //[Ndata.bb.y_bcp]
   Int_t           Ndata_bb_y_fcp;
   Double_t        bb_y_fcp[1];   //[Ndata.bb.y_fcp]
   Int_t           Ndata_bb_z_bcp;
   Double_t        bb_z_bcp[1];   //[Ndata.bb.z_bcp]
   Int_t           Ndata_bb_z_fcp;
   Double_t        bb_z_fcp[1];   //[Ndata.bb.z_fcp]
   Int_t           Ndata_sbs_hcal_Ref_tdc;
   Double_t        sbs_hcal_Ref_tdc[3];   //[Ndata.sbs.hcal.Ref.tdc]
   Int_t           Ndata_sbs_hcal_Ref_tdc_mult;
   Double_t        sbs_hcal_Ref_tdc_mult[3];   //[Ndata.sbs.hcal.Ref.tdc_mult]
   Int_t           Ndata_sbs_hcal_Ref_tdcelemID;
   Double_t        sbs_hcal_Ref_tdcelemID[3];   //[Ndata.sbs.hcal.Ref.tdcelemID]
   Int_t           Ndata_sbs_hcal_clus_again;
   Double_t        sbs_hcal_clus_again[10];   //[Ndata.sbs.hcal.clus.again]
   Int_t           Ndata_sbs_hcal_clus_atime;
   Double_t        sbs_hcal_clus_atime[10];   //[Ndata.sbs.hcal.clus.atime]
   Int_t           Ndata_sbs_hcal_clus_col;
   Double_t        sbs_hcal_clus_col[10];   //[Ndata.sbs.hcal.clus.col]
   Int_t           Ndata_sbs_hcal_clus_e;
   Double_t        sbs_hcal_clus_e[10];   //[Ndata.sbs.hcal.clus.e]
   Int_t           Ndata_sbs_hcal_clus_eblk;
   Double_t        sbs_hcal_clus_eblk[10];   //[Ndata.sbs.hcal.clus.eblk]
   Int_t           Ndata_sbs_hcal_clus_id;
   Double_t        sbs_hcal_clus_id[10];   //[Ndata.sbs.hcal.clus.id]
   Int_t           Ndata_sbs_hcal_clus_nblk;
   Double_t        sbs_hcal_clus_nblk[10];   //[Ndata.sbs.hcal.clus.nblk]
   Int_t           Ndata_sbs_hcal_clus_row;
   Double_t        sbs_hcal_clus_row[10];   //[Ndata.sbs.hcal.clus.row]
   Int_t           Ndata_sbs_hcal_clus_tdctime;
   Double_t        sbs_hcal_clus_tdctime[10];   //[Ndata.sbs.hcal.clus.tdctime]
   Int_t           Ndata_sbs_hcal_clus_x;
   Double_t        sbs_hcal_clus_x[10];   //[Ndata.sbs.hcal.clus.x]
   Int_t           Ndata_sbs_hcal_clus_y;
   Double_t        sbs_hcal_clus_y[10];   //[Ndata.sbs.hcal.clus.y]
   Int_t           Ndata_sbs_hcal_clus_blk_again;
   Double_t        sbs_hcal_clus_blk_again[5];   //[Ndata.sbs.hcal.clus_blk.again]
   Int_t           Ndata_sbs_hcal_clus_blk_atime;
   Double_t        sbs_hcal_clus_blk_atime[5];   //[Ndata.sbs.hcal.clus_blk.atime]
   Int_t           Ndata_sbs_hcal_clus_blk_col;
   Double_t        sbs_hcal_clus_blk_col[5];   //[Ndata.sbs.hcal.clus_blk.col]
   Int_t           Ndata_sbs_hcal_clus_blk_e;
   Double_t        sbs_hcal_clus_blk_e[5];   //[Ndata.sbs.hcal.clus_blk.e]
   Int_t           Ndata_sbs_hcal_clus_blk_id;
   Double_t        sbs_hcal_clus_blk_id[5];   //[Ndata.sbs.hcal.clus_blk.id]
   Int_t           Ndata_sbs_hcal_clus_blk_row;
   Double_t        sbs_hcal_clus_blk_row[5];   //[Ndata.sbs.hcal.clus_blk.row]
   Int_t           Ndata_sbs_hcal_clus_blk_tdctime;
   Double_t        sbs_hcal_clus_blk_tdctime[5];   //[Ndata.sbs.hcal.clus_blk.tdctime]
   Int_t           Ndata_sbs_hcal_clus_blk_x;
   Double_t        sbs_hcal_clus_blk_x[5];   //[Ndata.sbs.hcal.clus_blk.x]
   Int_t           Ndata_sbs_hcal_clus_blk_y;
   Double_t        sbs_hcal_clus_blk_y[5];   //[Ndata.sbs.hcal.clus_blk.y]
   Int_t           Ndata_sbs_hcal_goodblock_atime;
   Double_t        sbs_hcal_goodblock_atime[22];   //[Ndata.sbs.hcal.goodblock.atime]
   Int_t           Ndata_sbs_hcal_goodblock_cid;
   Double_t        sbs_hcal_goodblock_cid[22];   //[Ndata.sbs.hcal.goodblock.cid]
   Int_t           Ndata_sbs_hcal_goodblock_col;
   Double_t        sbs_hcal_goodblock_col[22];   //[Ndata.sbs.hcal.goodblock.col]
   Int_t           Ndata_sbs_hcal_goodblock_e;
   Double_t        sbs_hcal_goodblock_e[22];   //[Ndata.sbs.hcal.goodblock.e]
   Int_t           Ndata_sbs_hcal_goodblock_id;
   Double_t        sbs_hcal_goodblock_id[22];   //[Ndata.sbs.hcal.goodblock.id]
   Int_t           Ndata_sbs_hcal_goodblock_row;
   Double_t        sbs_hcal_goodblock_row[22];   //[Ndata.sbs.hcal.goodblock.row]
   Int_t           Ndata_sbs_hcal_goodblock_tdctime;
   Double_t        sbs_hcal_goodblock_tdctime[22];   //[Ndata.sbs.hcal.goodblock.tdctime]
   Int_t           Ndata_sbs_hcal_goodblock_x;
   Double_t        sbs_hcal_goodblock_x[22];   //[Ndata.sbs.hcal.goodblock.x]
   Int_t           Ndata_sbs_hcal_goodblock_y;
   Double_t        sbs_hcal_goodblock_y[22];   //[Ndata.sbs.hcal.goodblock.y]
   Int_t           Ndata_sbs_tr_beta;
   Double_t        sbs_tr_beta[1];   //[Ndata.sbs.tr.beta]
   Int_t           Ndata_sbs_tr_chi2;
   Double_t        sbs_tr_chi2[1];   //[Ndata.sbs.tr.chi2]
   Int_t           Ndata_sbs_tr_d_ph;
   Double_t        sbs_tr_d_ph[1];   //[Ndata.sbs.tr.d_ph]
   Int_t           Ndata_sbs_tr_d_th;
   Double_t        sbs_tr_d_th[1];   //[Ndata.sbs.tr.d_th]
   Int_t           Ndata_sbs_tr_d_x;
   Double_t        sbs_tr_d_x[1];   //[Ndata.sbs.tr.d_x]
   Int_t           Ndata_sbs_tr_d_y;
   Double_t        sbs_tr_d_y[1];   //[Ndata.sbs.tr.d_y]
   Int_t           Ndata_sbs_tr_dbeta;
   Double_t        sbs_tr_dbeta[1];   //[Ndata.sbs.tr.dbeta]
   Int_t           Ndata_sbs_tr_dtime;
   Double_t        sbs_tr_dtime[1];   //[Ndata.sbs.tr.dtime]
   Int_t           Ndata_sbs_tr_flag;
   Double_t        sbs_tr_flag[1];   //[Ndata.sbs.tr.flag]
   Int_t           Ndata_sbs_tr_ndof;
   Double_t        sbs_tr_ndof[1];   //[Ndata.sbs.tr.ndof]
   Int_t           Ndata_sbs_tr_p;
   Double_t        sbs_tr_p[1];   //[Ndata.sbs.tr.p]
   Int_t           Ndata_sbs_tr_pathl;
   Double_t        sbs_tr_pathl[1];   //[Ndata.sbs.tr.pathl]
   Int_t           Ndata_sbs_tr_ph;
   Double_t        sbs_tr_ph[1];   //[Ndata.sbs.tr.ph]
   Int_t           Ndata_sbs_tr_px;
   Double_t        sbs_tr_px[1];   //[Ndata.sbs.tr.px]
   Int_t           Ndata_sbs_tr_py;
   Double_t        sbs_tr_py[1];   //[Ndata.sbs.tr.py]
   Int_t           Ndata_sbs_tr_pz;
   Double_t        sbs_tr_pz[1];   //[Ndata.sbs.tr.pz]
   Int_t           Ndata_sbs_tr_r_ph;
   Double_t        sbs_tr_r_ph[1];   //[Ndata.sbs.tr.r_ph]
   Int_t           Ndata_sbs_tr_r_th;
   Double_t        sbs_tr_r_th[1];   //[Ndata.sbs.tr.r_th]
   Int_t           Ndata_sbs_tr_r_x;
   Double_t        sbs_tr_r_x[1];   //[Ndata.sbs.tr.r_x]
   Int_t           Ndata_sbs_tr_r_y;
   Double_t        sbs_tr_r_y[1];   //[Ndata.sbs.tr.r_y]
   Int_t           Ndata_sbs_tr_tg_dp;
   Double_t        sbs_tr_tg_dp[1];   //[Ndata.sbs.tr.tg_dp]
   Int_t           Ndata_sbs_tr_tg_ph;
   Double_t        sbs_tr_tg_ph[1];   //[Ndata.sbs.tr.tg_ph]
   Int_t           Ndata_sbs_tr_tg_th;
   Double_t        sbs_tr_tg_th[1];   //[Ndata.sbs.tr.tg_th]
   Int_t           Ndata_sbs_tr_tg_y;
   Double_t        sbs_tr_tg_y[1];   //[Ndata.sbs.tr.tg_y]
   Int_t           Ndata_sbs_tr_th;
   Double_t        sbs_tr_th[1];   //[Ndata.sbs.tr.th]
   Int_t           Ndata_sbs_tr_time;
   Double_t        sbs_tr_time[1];   //[Ndata.sbs.tr.time]
   Int_t           Ndata_sbs_tr_vx;
   Double_t        sbs_tr_vx[1];   //[Ndata.sbs.tr.vx]
   Int_t           Ndata_sbs_tr_vy;
   Double_t        sbs_tr_vy[1];   //[Ndata.sbs.tr.vy]
   Int_t           Ndata_sbs_tr_vz;
   Double_t        sbs_tr_vz[1];   //[Ndata.sbs.tr.vz]
   Int_t           Ndata_sbs_tr_x;
   Double_t        sbs_tr_x[1];   //[Ndata.sbs.tr.x]
   Int_t           Ndata_sbs_tr_y;
   Double_t        sbs_tr_y[1];   //[Ndata.sbs.tr.y]
   Int_t           Ndata_sbs_trig_a_amp_p;
   Double_t        sbs_trig_a_amp_p[6];   //[Ndata.sbs.trig.a_amp_p]
   Int_t           Ndata_sbs_trig_a_p;
   Double_t        sbs_trig_a_p[6];   //[Ndata.sbs.trig.a_p]
   Int_t           Ndata_sbs_trig_a_time;
   Double_t        sbs_trig_a_time[6];   //[Ndata.sbs.trig.a_time]
   Int_t           Ndata_sbs_trig_adcelemID;
   Double_t        sbs_trig_adcelemID[6];   //[Ndata.sbs.trig.adcelemID]
   Int_t           Ndata_sbs_x_bcp;
   Double_t        sbs_x_bcp[1];   //[Ndata.sbs.x_bcp]
   Int_t           Ndata_sbs_x_fcp;
   Double_t        sbs_x_fcp[1];   //[Ndata.sbs.x_fcp]
   Int_t           Ndata_sbs_y_bcp;
   Double_t        sbs_y_bcp[1];   //[Ndata.sbs.y_bcp]
   Int_t           Ndata_sbs_y_fcp;
   Double_t        sbs_y_fcp[1];   //[Ndata.sbs.y_fcp]
   Int_t           Ndata_sbs_z_bcp;
   Double_t        sbs_z_bcp[1];   //[Ndata.sbs.z_bcp]
   Int_t           Ndata_sbs_z_fcp;
   Double_t        sbs_z_fcp[1];   //[Ndata.sbs.z_fcp]
   Double_t        BB_gold_beta;
   Double_t        BB_gold_dp;
   Double_t        BB_gold_index;
   Double_t        BB_gold_ok;
   Double_t        BB_gold_p;
   Double_t        BB_gold_ph;
   Double_t        BB_gold_px;
   Double_t        BB_gold_py;
   Double_t        BB_gold_pz;
   Double_t        BB_gold_th;
   Double_t        BB_gold_x;
   Double_t        BB_gold_y;
   Double_t        Lrb_BPMA_rawcur_1;
   Double_t        Lrb_BPMA_rawcur_2;
   Double_t        Lrb_BPMA_rawcur_3;
   Double_t        Lrb_BPMA_rawcur_4;
   Double_t        Lrb_BPMA_rotpos1;
   Double_t        Lrb_BPMA_rotpos2;
   Double_t        Lrb_BPMA_x;
   Double_t        Lrb_BPMA_y;
   Double_t        Lrb_BPMA_z;
   Double_t        Lrb_BPMB_rawcur_1;
   Double_t        Lrb_BPMB_rawcur_2;
   Double_t        Lrb_BPMB_rawcur_3;
   Double_t        Lrb_BPMB_rawcur_4;
   Double_t        Lrb_BPMB_rotpos1;
   Double_t        Lrb_BPMB_rotpos2;
   Double_t        Lrb_BPMB_x;
   Double_t        Lrb_BPMB_y;
   Double_t        Lrb_BPMB_z;
   Double_t        Lrb_Raster_bpma_x;
   Double_t        Lrb_Raster_bpma_y;
   Double_t        Lrb_Raster_bpma_z;
   Double_t        Lrb_Raster_bpmb_x;
   Double_t        Lrb_Raster_bpmb_y;
   Double_t        Lrb_Raster_bpmb_z;
   Double_t        Lrb_Raster_rawcur_x;
   Double_t        Lrb_Raster_rawcur_y;
   Double_t        Lrb_Raster_rawslope_x;
   Double_t        Lrb_Raster_rawslope_y;
   Double_t        Lrb_Raster_target_dir_x;
   Double_t        Lrb_Raster_target_dir_y;
   Double_t        Lrb_Raster_target_dir_z;
   Double_t        Lrb_Raster_target_x;
   Double_t        Lrb_Raster_target_y;
   Double_t        Lrb_Raster_target_z;
   Double_t        Lrb_Raster2_bpma_x;
   Double_t        Lrb_Raster2_bpma_y;
   Double_t        Lrb_Raster2_bpma_z;
   Double_t        Lrb_Raster2_bpmb_x;
   Double_t        Lrb_Raster2_bpmb_y;
   Double_t        Lrb_Raster2_bpmb_z;
   Double_t        Lrb_Raster2_rawcur_x;
   Double_t        Lrb_Raster2_rawcur_y;
   Double_t        Lrb_Raster2_rawslope_x;
   Double_t        Lrb_Raster2_rawslope_y;
   Double_t        Lrb_Raster2_target_dir_x;
   Double_t        Lrb_Raster2_target_dir_y;
   Double_t        Lrb_Raster2_target_dir_z;
   Double_t        Lrb_Raster2_target_x;
   Double_t        Lrb_Raster2_target_y;
   Double_t        Lrb_Raster2_target_z;
   Double_t        SBSrb_BPMA_rawcur_1;
   Double_t        SBSrb_BPMA_rawcur_2;
   Double_t        SBSrb_BPMA_rawcur_3;
   Double_t        SBSrb_BPMA_rawcur_4;
   Double_t        SBSrb_BPMA_rotpos1;
   Double_t        SBSrb_BPMA_rotpos2;
   Double_t        SBSrb_BPMA_x;
   Double_t        SBSrb_BPMA_y;
   Double_t        SBSrb_BPMA_z;
   Double_t        SBSrb_BPMB_rawcur_1;
   Double_t        SBSrb_BPMB_rawcur_2;
   Double_t        SBSrb_BPMB_rawcur_3;
   Double_t        SBSrb_BPMB_rawcur_4;
   Double_t        SBSrb_BPMB_rotpos1;
   Double_t        SBSrb_BPMB_rotpos2;
   Double_t        SBSrb_BPMB_x;
   Double_t        SBSrb_BPMB_y;
   Double_t        SBSrb_BPMB_z;
   Double_t        SBSrb_Raster_bpma_x;
   Double_t        SBSrb_Raster_bpma_y;
   Double_t        SBSrb_Raster_bpma_z;
   Double_t        SBSrb_Raster_bpmb_x;
   Double_t        SBSrb_Raster_bpmb_y;
   Double_t        SBSrb_Raster_bpmb_z;
   Double_t        SBSrb_Raster_rawcur_x;
   Double_t        SBSrb_Raster_rawcur_y;
   Double_t        SBSrb_Raster_rawslope_x;
   Double_t        SBSrb_Raster_rawslope_y;
   Double_t        SBSrb_Raster_target_dir_x;
   Double_t        SBSrb_Raster_target_dir_y;
   Double_t        SBSrb_Raster_target_dir_z;
   Double_t        SBSrb_Raster_target_x;
   Double_t        SBSrb_Raster_target_y;
   Double_t        SBSrb_Raster_target_z;
   Double_t        SBSrb_Raster2_bpma_x;
   Double_t        SBSrb_Raster2_bpma_y;
   Double_t        SBSrb_Raster2_bpma_z;
   Double_t        SBSrb_Raster2_bpmb_x;
   Double_t        SBSrb_Raster2_bpmb_y;
   Double_t        SBSrb_Raster2_bpmb_z;
   Double_t        SBSrb_Raster2_rawcur_x;
   Double_t        SBSrb_Raster2_rawcur_y;
   Double_t        SBSrb_Raster2_rawslope_x;
   Double_t        SBSrb_Raster2_rawslope_y;
   Double_t        SBSrb_Raster2_target_dir_x;
   Double_t        SBSrb_Raster2_target_dir_y;
   Double_t        SBSrb_Raster2_target_dir_z;
   Double_t        SBSrb_Raster2_target_x;
   Double_t        SBSrb_Raster2_target_y;
   Double_t        SBSrb_Raster2_target_z;
   Double_t        bb_gem_hit_ngoodhits;
   Double_t        bb_gem_m0_clust_nclustu;
   Double_t        bb_gem_m0_clust_nclustu_neg;
   Double_t        bb_gem_m0_clust_nclustu_tot;
   Double_t        bb_gem_m0_clust_nclustv;
   Double_t        bb_gem_m0_clust_nclustv_neg;
   Double_t        bb_gem_m0_clust_nclustv_tot;
   Double_t        bb_gem_m0_strip_nstrips_keep;
   Double_t        bb_gem_m0_strip_nstrips_keepU;
   Double_t        bb_gem_m0_strip_nstrips_keepV;
   Double_t        bb_gem_m0_strip_nstrips_keep_lmax;
   Double_t        bb_gem_m0_strip_nstrips_keep_lmaxU;
   Double_t        bb_gem_m0_strip_nstrips_keep_lmaxV;
   Double_t        bb_gem_m0_strip_nstripsfired;
   Double_t        bb_gem_m0_strip_nstripsfired_neg;
   Double_t        bb_gem_m1_clust_nclustu;
   Double_t        bb_gem_m1_clust_nclustu_neg;
   Double_t        bb_gem_m1_clust_nclustu_tot;
   Double_t        bb_gem_m1_clust_nclustv;
   Double_t        bb_gem_m1_clust_nclustv_neg;
   Double_t        bb_gem_m1_clust_nclustv_tot;
   Double_t        bb_gem_m1_strip_nstrips_keep;
   Double_t        bb_gem_m1_strip_nstrips_keepU;
   Double_t        bb_gem_m1_strip_nstrips_keepV;
   Double_t        bb_gem_m1_strip_nstrips_keep_lmax;
   Double_t        bb_gem_m1_strip_nstrips_keep_lmaxU;
   Double_t        bb_gem_m1_strip_nstrips_keep_lmaxV;
   Double_t        bb_gem_m1_strip_nstripsfired;
   Double_t        bb_gem_m1_strip_nstripsfired_neg;
   Double_t        bb_gem_m2_clust_nclustu;
   Double_t        bb_gem_m2_clust_nclustu_neg;
   Double_t        bb_gem_m2_clust_nclustu_tot;
   Double_t        bb_gem_m2_clust_nclustv;
   Double_t        bb_gem_m2_clust_nclustv_neg;
   Double_t        bb_gem_m2_clust_nclustv_tot;
   Double_t        bb_gem_m2_strip_nstrips_keep;
   Double_t        bb_gem_m2_strip_nstrips_keepU;
   Double_t        bb_gem_m2_strip_nstrips_keepV;
   Double_t        bb_gem_m2_strip_nstrips_keep_lmax;
   Double_t        bb_gem_m2_strip_nstrips_keep_lmaxU;
   Double_t        bb_gem_m2_strip_nstrips_keep_lmaxV;
   Double_t        bb_gem_m2_strip_nstripsfired;
   Double_t        bb_gem_m2_strip_nstripsfired_neg;
   Double_t        bb_gem_m3_clust_nclustu;
   Double_t        bb_gem_m3_clust_nclustu_neg;
   Double_t        bb_gem_m3_clust_nclustu_tot;
   Double_t        bb_gem_m3_clust_nclustv;
   Double_t        bb_gem_m3_clust_nclustv_neg;
   Double_t        bb_gem_m3_clust_nclustv_tot;
   Double_t        bb_gem_m3_strip_nstrips_keep;
   Double_t        bb_gem_m3_strip_nstrips_keepU;
   Double_t        bb_gem_m3_strip_nstrips_keepV;
   Double_t        bb_gem_m3_strip_nstrips_keep_lmax;
   Double_t        bb_gem_m3_strip_nstrips_keep_lmaxU;
   Double_t        bb_gem_m3_strip_nstrips_keep_lmaxV;
   Double_t        bb_gem_m3_strip_nstripsfired;
   Double_t        bb_gem_m3_strip_nstripsfired_neg;
   Double_t        bb_gem_m4_clust_nclustu;
   Double_t        bb_gem_m4_clust_nclustu_neg;
   Double_t        bb_gem_m4_clust_nclustu_tot;
   Double_t        bb_gem_m4_clust_nclustv;
   Double_t        bb_gem_m4_clust_nclustv_neg;
   Double_t        bb_gem_m4_clust_nclustv_tot;
   Double_t        bb_gem_m4_strip_nstrips_keep;
   Double_t        bb_gem_m4_strip_nstrips_keepU;
   Double_t        bb_gem_m4_strip_nstrips_keepV;
   Double_t        bb_gem_m4_strip_nstrips_keep_lmax;
   Double_t        bb_gem_m4_strip_nstrips_keep_lmaxU;
   Double_t        bb_gem_m4_strip_nstrips_keep_lmaxV;
   Double_t        bb_gem_m4_strip_nstripsfired;
   Double_t        bb_gem_m4_strip_nstripsfired_neg;
   Double_t        bb_gem_m5_clust_nclustu;
   Double_t        bb_gem_m5_clust_nclustu_neg;
   Double_t        bb_gem_m5_clust_nclustu_tot;
   Double_t        bb_gem_m5_clust_nclustv;
   Double_t        bb_gem_m5_clust_nclustv_neg;
   Double_t        bb_gem_m5_clust_nclustv_tot;
   Double_t        bb_gem_m5_strip_nstrips_keep;
   Double_t        bb_gem_m5_strip_nstrips_keepU;
   Double_t        bb_gem_m5_strip_nstrips_keepV;
   Double_t        bb_gem_m5_strip_nstrips_keep_lmax;
   Double_t        bb_gem_m5_strip_nstrips_keep_lmaxU;
   Double_t        bb_gem_m5_strip_nstrips_keep_lmaxV;
   Double_t        bb_gem_m5_strip_nstripsfired;
   Double_t        bb_gem_m5_strip_nstripsfired_neg;
   Double_t        bb_gem_m6_clust_nclustu;
   Double_t        bb_gem_m6_clust_nclustu_neg;
   Double_t        bb_gem_m6_clust_nclustu_tot;
   Double_t        bb_gem_m6_clust_nclustv;
   Double_t        bb_gem_m6_clust_nclustv_neg;
   Double_t        bb_gem_m6_clust_nclustv_tot;
   Double_t        bb_gem_m6_strip_nstrips_keep;
   Double_t        bb_gem_m6_strip_nstrips_keepU;
   Double_t        bb_gem_m6_strip_nstrips_keepV;
   Double_t        bb_gem_m6_strip_nstrips_keep_lmax;
   Double_t        bb_gem_m6_strip_nstrips_keep_lmaxU;
   Double_t        bb_gem_m6_strip_nstrips_keep_lmaxV;
   Double_t        bb_gem_m6_strip_nstripsfired;
   Double_t        bb_gem_m6_strip_nstripsfired_neg;
   Double_t        bb_gem_m7_clust_nclustu;
   Double_t        bb_gem_m7_clust_nclustu_neg;
   Double_t        bb_gem_m7_clust_nclustu_tot;
   Double_t        bb_gem_m7_clust_nclustv;
   Double_t        bb_gem_m7_clust_nclustv_neg;
   Double_t        bb_gem_m7_clust_nclustv_tot;
   Double_t        bb_gem_m7_strip_nstrips_keep;
   Double_t        bb_gem_m7_strip_nstrips_keepU;
   Double_t        bb_gem_m7_strip_nstrips_keepV;
   Double_t        bb_gem_m7_strip_nstrips_keep_lmax;
   Double_t        bb_gem_m7_strip_nstrips_keep_lmaxU;
   Double_t        bb_gem_m7_strip_nstrips_keep_lmaxV;
   Double_t        bb_gem_m7_strip_nstripsfired;
   Double_t        bb_gem_m7_strip_nstripsfired_neg;
   Double_t        bb_gem_nlayershit;
   Double_t        bb_gem_nlayershitu;
   Double_t        bb_gem_nlayershituv;
   Double_t        bb_gem_nlayershitv;
   Double_t        bb_gem_track_besttrack;
   Double_t        bb_gem_track_ntrack;
   Double_t        bb_gem_trigtime;
   Double_t        bb_grinch_tdc_bestcluster;
   Double_t        bb_grinch_tdc_clus_adc;
   Double_t        bb_grinch_tdc_clus_dx;
   Double_t        bb_grinch_tdc_clus_dy;
   Double_t        bb_grinch_tdc_clus_mirrorindex;
   Double_t        bb_grinch_tdc_clus_size;
   Double_t        bb_grinch_tdc_clus_t_mean;
   Double_t        bb_grinch_tdc_clus_t_rms;
   Double_t        bb_grinch_tdc_clus_tot_mean;
   Double_t        bb_grinch_tdc_clus_trackindex;
   Double_t        bb_grinch_tdc_clus_x_mean;
   Double_t        bb_grinch_tdc_clus_y_mean;
   Double_t        bb_grinch_tdc_nclus;
   Double_t        bb_grinch_tdc_ngoodhits;
   Double_t        bb_grinch_tdc_ntrackmatch;
   Double_t        bb_hodotdc_nclus;
   Double_t        bb_ps_againblk;
   Double_t        bb_ps_atimeblk;
   Double_t        bb_ps_colblk;
   Double_t        bb_ps_e;
   Double_t        bb_ps_eblk;
   Double_t        bb_ps_idblk;
   Double_t        bb_ps_index;
   Double_t        bb_ps_nblk;
   Double_t        bb_ps_nclus;
   Double_t        bb_ps_ngoodADChits;
   Double_t        bb_ps_rowblk;
   Double_t        bb_ps_x;
   Double_t        bb_ps_y;
   Double_t        bb_sh_againblk;
   Double_t        bb_sh_atimeblk;
   Double_t        bb_sh_colblk;
   Double_t        bb_sh_e;
   Double_t        bb_sh_eblk;
   Double_t        bb_sh_idblk;
   Double_t        bb_sh_index;
   Double_t        bb_sh_nblk;
   Double_t        bb_sh_nclus;
   Double_t        bb_sh_ngoodADChits;
   Double_t        bb_sh_rowblk;
   Double_t        bb_sh_x;
   Double_t        bb_sh_y;
   Double_t        bb_tr_n;
   Double_t        bb_ts_over_threshold;
   Double_t        e_kine_Q2;
   Double_t        e_kine_W2;
   Double_t        e_kine_angle;
   Double_t        e_kine_epsilon;
   Double_t        e_kine_nu;
   Double_t        e_kine_omega;
   Double_t        e_kine_ph_q;
   Double_t        e_kine_q3m;
   Double_t        e_kine_q_x;
   Double_t        e_kine_q_y;
   Double_t        e_kine_q_z;
   Double_t        e_kine_th_q;
   Double_t        e_kine_x_bj;
   Double_t        g_datatype;
   Double_t        g_evlen;
   Double_t        g_evnum;
   Double_t        g_evtime;
   Double_t        g_evtyp;
   Double_t        g_runnum;
   Double_t        g_runtime;
   Double_t        g_runtype;
   Double_t        g_trigbits;
   Double_t        sbs_HCALdir_x;
   Double_t        sbs_HCALdir_y;
   Double_t        sbs_HCALdir_z;
   Double_t        sbs_HCALph_n;
   Double_t        sbs_HCALth_n;
   Double_t        sbs_hcal_againblk;
   Double_t        sbs_hcal_atimeblk;
   Double_t        sbs_hcal_colblk;
   Double_t        sbs_hcal_e;
   Double_t        sbs_hcal_eblk;
   Double_t        sbs_hcal_idblk;
   Double_t        sbs_hcal_index;
   Double_t        sbs_hcal_nblk;
   Double_t        sbs_hcal_nclus;
   Double_t        sbs_hcal_ngoodADChits;
   Double_t        sbs_hcal_rowblk;
   Double_t        sbs_hcal_tdctimeblk;
   Double_t        sbs_hcal_x;
   Double_t        sbs_hcal_y;
   Double_t        sbs_tr_n;
   Double_t        singletrack;
   Double_t        anytrack;
   Double_t        HALLA_p;
 //THaEvent        *Event_Branch;
   ULong64_t       fEvtHdr_fEvtTime;
   UInt_t          fEvtHdr_fEvtNum;
   UInt_t          fEvtHdr_fEvtType;
   UInt_t          fEvtHdr_fEvtLen;
   Int_t           fEvtHdr_fHelicity;
   UInt_t          fEvtHdr_fTrigBits;
   UInt_t          fEvtHdr_fRun;

   // List of branches
   TBranch        *b_ax_new;   //!
   TBranch        *b_bx_new;   //!
   TBranch        *b_ay_new;   //!
   TBranch        *b_by_new;   //!
   TBranch        *b_targx;   //!
   TBranch        *b_targy;   //!
   TBranch        *b_Ndata_bb_bbtrig_a_amp_p;   //!
   TBranch        *b_bb_bbtrig_a_amp_p;   //!
   TBranch        *b_Ndata_bb_bbtrig_a_time;   //!
   TBranch        *b_bb_bbtrig_a_time;   //!
   TBranch        *b_Ndata_bb_bbtrig_adcelemID;   //!
   TBranch        *b_bb_bbtrig_adcelemID;   //!
   TBranch        *b_Ndata_bb_bbtrig_tdc;   //!
   TBranch        *b_bb_bbtrig_tdc;   //!
   TBranch        *b_Ndata_bb_bbtrig_tdcelemID;   //!
   TBranch        *b_bb_bbtrig_tdcelemID;   //!
   TBranch        *b_Ndata_bb_eps_over_etot;   //!
   TBranch        *b_bb_eps_over_etot;   //!
   TBranch        *b_Ndata_bb_etot_over_p;   //!
   TBranch        *b_bb_etot_over_p;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCU;   //!
   TBranch        *b_bb_gem_hit_ADCU;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCU_deconv;   //!
   TBranch        *b_bb_gem_hit_ADCU_deconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCV;   //!
   TBranch        *b_bb_gem_hit_ADCV;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCV_deconv;   //!
   TBranch        *b_bb_gem_hit_ADCV_deconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCasym;   //!
   TBranch        *b_bb_gem_hit_ADCasym;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCasym_deconv;   //!
   TBranch        *b_bb_gem_hit_ADCasym_deconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCavg;   //!
   TBranch        *b_bb_gem_hit_ADCavg;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCavg_deconv;   //!
   TBranch        *b_bb_gem_hit_ADCavg_deconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCfrac0_Umax;   //!
   TBranch        *b_bb_gem_hit_ADCfrac0_Umax;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCfrac0_Vmax;   //!
   TBranch        *b_bb_gem_hit_ADCfrac0_Vmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCfrac1_Umax;   //!
   TBranch        *b_bb_gem_hit_ADCfrac1_Umax;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCfrac1_Vmax;   //!
   TBranch        *b_bb_gem_hit_ADCfrac1_Vmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCfrac2_Umax;   //!
   TBranch        *b_bb_gem_hit_ADCfrac2_Umax;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCfrac2_Vmax;   //!
   TBranch        *b_bb_gem_hit_ADCfrac2_Vmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCfrac3_Umax;   //!
   TBranch        *b_bb_gem_hit_ADCfrac3_Umax;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCfrac3_Vmax;   //!
   TBranch        *b_bb_gem_hit_ADCfrac3_Vmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCfrac4_Umax;   //!
   TBranch        *b_bb_gem_hit_ADCfrac4_Umax;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCfrac4_Vmax;   //!
   TBranch        *b_bb_gem_hit_ADCfrac4_Vmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCfrac5_Umax;   //!
   TBranch        *b_bb_gem_hit_ADCfrac5_Umax;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCfrac5_Vmax;   //!
   TBranch        *b_bb_gem_hit_ADCfrac5_Vmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxcomboUclust_deconv;   //!
   TBranch        *b_bb_gem_hit_ADCmaxcomboUclust_deconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxcomboVclust_deconv;   //!
   TBranch        *b_bb_gem_hit_ADCmaxcomboVclust_deconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxsampU;   //!
   TBranch        *b_bb_gem_hit_ADCmaxsampU;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxsampUclust;   //!
   TBranch        *b_bb_gem_hit_ADCmaxsampUclust;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxsampUclust_deconv;   //!
   TBranch        *b_bb_gem_hit_ADCmaxsampUclust_deconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxsampV;   //!
   TBranch        *b_bb_gem_hit_ADCmaxsampV;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxsampVclust;   //!
   TBranch        *b_bb_gem_hit_ADCmaxsampVclust;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxsampVclust_deconv;   //!
   TBranch        *b_bb_gem_hit_ADCmaxsampVclust_deconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxstripU;   //!
   TBranch        *b_bb_gem_hit_ADCmaxstripU;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxstripV;   //!
   TBranch        *b_bb_gem_hit_ADCmaxstripV;   //!
   TBranch        *b_Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_U;   //!
   TBranch        *b_bb_gem_hit_BUILD_ALL_SAMPLES_U;   //!
   TBranch        *b_Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_V;   //!
   TBranch        *b_bb_gem_hit_BUILD_ALL_SAMPLES_V;   //!
   TBranch        *b_Ndata_bb_gem_hit_CM_GOOD_U;   //!
   TBranch        *b_bb_gem_hit_CM_GOOD_U;   //!
   TBranch        *b_Ndata_bb_gem_hit_CM_GOOD_V;   //!
   TBranch        *b_bb_gem_hit_CM_GOOD_V;   //!
   TBranch        *b_Ndata_bb_gem_hit_DeconvADC0_Umax;   //!
   TBranch        *b_bb_gem_hit_DeconvADC0_Umax;   //!
   TBranch        *b_Ndata_bb_gem_hit_DeconvADC0_Vmax;   //!
   TBranch        *b_bb_gem_hit_DeconvADC0_Vmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_DeconvADC1_Umax;   //!
   TBranch        *b_bb_gem_hit_DeconvADC1_Umax;   //!
   TBranch        *b_Ndata_bb_gem_hit_DeconvADC1_Vmax;   //!
   TBranch        *b_bb_gem_hit_DeconvADC1_Vmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_DeconvADC2_Umax;   //!
   TBranch        *b_bb_gem_hit_DeconvADC2_Umax;   //!
   TBranch        *b_Ndata_bb_gem_hit_DeconvADC2_Vmax;   //!
   TBranch        *b_bb_gem_hit_DeconvADC2_Vmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_DeconvADC3_Umax;   //!
   TBranch        *b_bb_gem_hit_DeconvADC3_Umax;   //!
   TBranch        *b_Ndata_bb_gem_hit_DeconvADC3_Vmax;   //!
   TBranch        *b_bb_gem_hit_DeconvADC3_Vmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_DeconvADC4_Umax;   //!
   TBranch        *b_bb_gem_hit_DeconvADC4_Umax;   //!
   TBranch        *b_Ndata_bb_gem_hit_DeconvADC4_Vmax;   //!
   TBranch        *b_bb_gem_hit_DeconvADC4_Vmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_DeconvADC5_Umax;   //!
   TBranch        *b_bb_gem_hit_DeconvADC5_Umax;   //!
   TBranch        *b_Ndata_bb_gem_hit_DeconvADC5_Vmax;   //!
   TBranch        *b_bb_gem_hit_DeconvADC5_Vmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_DeconvADCmaxcomboU;   //!
   TBranch        *b_bb_gem_hit_DeconvADCmaxcomboU;   //!
   TBranch        *b_Ndata_bb_gem_hit_DeconvADCmaxcomboV;   //!
   TBranch        *b_bb_gem_hit_DeconvADCmaxcomboV;   //!
   TBranch        *b_Ndata_bb_gem_hit_DeconvADCmaxsampU;   //!
   TBranch        *b_bb_gem_hit_DeconvADCmaxsampU;   //!
   TBranch        *b_Ndata_bb_gem_hit_DeconvADCmaxsampV;   //!
   TBranch        *b_bb_gem_hit_DeconvADCmaxsampV;   //!
   TBranch        *b_Ndata_bb_gem_hit_DeconvADCmaxstripU;   //!
   TBranch        *b_bb_gem_hit_DeconvADCmaxstripU;   //!
   TBranch        *b_Ndata_bb_gem_hit_DeconvADCmaxstripV;   //!
   TBranch        *b_bb_gem_hit_DeconvADCmaxstripV;   //!
   TBranch        *b_Ndata_bb_gem_hit_ENABLE_CM_U;   //!
   TBranch        *b_bb_gem_hit_ENABLE_CM_U;   //!
   TBranch        *b_Ndata_bb_gem_hit_ENABLE_CM_V;   //!
   TBranch        *b_bb_gem_hit_ENABLE_CM_V;   //!
   TBranch        *b_Ndata_bb_gem_hit_TSchi2_Umax;   //!
   TBranch        *b_bb_gem_hit_TSchi2_Umax;   //!
   TBranch        *b_Ndata_bb_gem_hit_TSchi2_Vmax;   //!
   TBranch        *b_bb_gem_hit_TSchi2_Vmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_TSprob_Umax;   //!
   TBranch        *b_bb_gem_hit_TSprob_Umax;   //!
   TBranch        *b_Ndata_bb_gem_hit_TSprob_Vmax;   //!
   TBranch        *b_bb_gem_hit_TSprob_Vmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_Tavg;   //!
   TBranch        *b_bb_gem_hit_Tavg;   //!
   TBranch        *b_Ndata_bb_gem_hit_Tavg_corr;   //!
   TBranch        *b_bb_gem_hit_Tavg_corr;   //!
   TBranch        *b_Ndata_bb_gem_hit_Tavg_deconv;   //!
   TBranch        *b_bb_gem_hit_Tavg_deconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_Tavg_fit;   //!
   TBranch        *b_bb_gem_hit_Tavg_fit;   //!
   TBranch        *b_Ndata_bb_gem_hit_Ugain;   //!
   TBranch        *b_bb_gem_hit_Ugain;   //!
   TBranch        *b_Ndata_bb_gem_hit_Utime;   //!
   TBranch        *b_bb_gem_hit_Utime;   //!
   TBranch        *b_Ndata_bb_gem_hit_UtimeDeconv;   //!
   TBranch        *b_bb_gem_hit_UtimeDeconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_UtimeFit;   //!
   TBranch        *b_bb_gem_hit_UtimeFit;   //!
   TBranch        *b_Ndata_bb_gem_hit_UtimeMaxStrip;   //!
   TBranch        *b_bb_gem_hit_UtimeMaxStrip;   //!
   TBranch        *b_Ndata_bb_gem_hit_UtimeMaxStripDeconv;   //!
   TBranch        *b_bb_gem_hit_UtimeMaxStripDeconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_UtimeMaxStripFit;   //!
   TBranch        *b_bb_gem_hit_UtimeMaxStripFit;   //!
   TBranch        *b_Ndata_bb_gem_hit_Vgain;   //!
   TBranch        *b_bb_gem_hit_Vgain;   //!
   TBranch        *b_Ndata_bb_gem_hit_Vtime;   //!
   TBranch        *b_bb_gem_hit_Vtime;   //!
   TBranch        *b_Ndata_bb_gem_hit_VtimeDeconv;   //!
   TBranch        *b_bb_gem_hit_VtimeDeconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_VtimeFit;   //!
   TBranch        *b_bb_gem_hit_VtimeFit;   //!
   TBranch        *b_Ndata_bb_gem_hit_VtimeMaxStrip;   //!
   TBranch        *b_bb_gem_hit_VtimeMaxStrip;   //!
   TBranch        *b_Ndata_bb_gem_hit_VtimeMaxStripDeconv;   //!
   TBranch        *b_bb_gem_hit_VtimeMaxStripDeconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_VtimeMaxStripFit;   //!
   TBranch        *b_bb_gem_hit_VtimeMaxStripFit;   //!
   TBranch        *b_Ndata_bb_gem_hit_ccor_clust;   //!
   TBranch        *b_bb_gem_hit_ccor_clust;   //!
   TBranch        *b_Ndata_bb_gem_hit_ccor_clust_deconv;   //!
   TBranch        *b_bb_gem_hit_ccor_clust_deconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_ccor_strip;   //!
   TBranch        *b_bb_gem_hit_ccor_strip;   //!
   TBranch        *b_Ndata_bb_gem_hit_ccor_strip_deconv;   //!
   TBranch        *b_bb_gem_hit_ccor_strip_deconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_deltat;   //!
   TBranch        *b_bb_gem_hit_deltat;   //!
   TBranch        *b_Ndata_bb_gem_hit_deltat_deconv;   //!
   TBranch        *b_bb_gem_hit_deltat_deconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_deltat_fit;   //!
   TBranch        *b_bb_gem_hit_deltat_fit;   //!
   TBranch        *b_Ndata_bb_gem_hit_eresidu;   //!
   TBranch        *b_bb_gem_hit_eresidu;   //!
   TBranch        *b_Ndata_bb_gem_hit_eresidv;   //!
   TBranch        *b_bb_gem_hit_eresidv;   //!
   TBranch        *b_Ndata_bb_gem_hit_icombomaxUclustDeconv;   //!
   TBranch        *b_bb_gem_hit_icombomaxUclustDeconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_icombomaxUstripDeconv;   //!
   TBranch        *b_bb_gem_hit_icombomaxUstripDeconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_icombomaxVclustDeconv;   //!
   TBranch        *b_bb_gem_hit_icombomaxVclustDeconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_icombomaxVstripDeconv;   //!
   TBranch        *b_bb_gem_hit_icombomaxVstripDeconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_isampmaxUclust;   //!
   TBranch        *b_bb_gem_hit_isampmaxUclust;   //!
   TBranch        *b_Ndata_bb_gem_hit_isampmaxUclustDeconv;   //!
   TBranch        *b_bb_gem_hit_isampmaxUclustDeconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_isampmaxUstrip;   //!
   TBranch        *b_bb_gem_hit_isampmaxUstrip;   //!
   TBranch        *b_Ndata_bb_gem_hit_isampmaxUstripDeconv;   //!
   TBranch        *b_bb_gem_hit_isampmaxUstripDeconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_isampmaxVclust;   //!
   TBranch        *b_bb_gem_hit_isampmaxVclust;   //!
   TBranch        *b_Ndata_bb_gem_hit_isampmaxVclustDeconv;   //!
   TBranch        *b_bb_gem_hit_isampmaxVclustDeconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_isampmaxVstrip;   //!
   TBranch        *b_bb_gem_hit_isampmaxVstrip;   //!
   TBranch        *b_Ndata_bb_gem_hit_isampmaxVstripDeconv;   //!
   TBranch        *b_bb_gem_hit_isampmaxVstripDeconv;   //!
   TBranch        *b_Ndata_bb_gem_hit_layer;   //!
   TBranch        *b_bb_gem_hit_layer;   //!
   TBranch        *b_Ndata_bb_gem_hit_module;   //!
   TBranch        *b_bb_gem_hit_module;   //!
   TBranch        *b_Ndata_bb_gem_hit_nstripu;   //!
   TBranch        *b_bb_gem_hit_nstripu;   //!
   TBranch        *b_Ndata_bb_gem_hit_nstripv;   //!
   TBranch        *b_bb_gem_hit_nstripv;   //!
   TBranch        *b_Ndata_bb_gem_hit_residu;   //!
   TBranch        *b_bb_gem_hit_residu;   //!
   TBranch        *b_Ndata_bb_gem_hit_residv;   //!
   TBranch        *b_bb_gem_hit_residv;   //!
   TBranch        *b_Ndata_bb_gem_hit_trackindex;   //!
   TBranch        *b_bb_gem_hit_trackindex;   //!
   TBranch        *b_Ndata_bb_gem_hit_u;   //!
   TBranch        *b_bb_gem_hit_u;   //!
   TBranch        *b_Ndata_bb_gem_hit_umoment;   //!
   TBranch        *b_bb_gem_hit_umoment;   //!
   TBranch        *b_Ndata_bb_gem_hit_usigma;   //!
   TBranch        *b_bb_gem_hit_usigma;   //!
   TBranch        *b_Ndata_bb_gem_hit_ustriphi;   //!
   TBranch        *b_bb_gem_hit_ustriphi;   //!
   TBranch        *b_Ndata_bb_gem_hit_ustriplo;   //!
   TBranch        *b_bb_gem_hit_ustriplo;   //!
   TBranch        *b_Ndata_bb_gem_hit_ustripmax;   //!
   TBranch        *b_bb_gem_hit_ustripmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_v;   //!
   TBranch        *b_bb_gem_hit_v;   //!
   TBranch        *b_Ndata_bb_gem_hit_vmoment;   //!
   TBranch        *b_bb_gem_hit_vmoment;   //!
   TBranch        *b_Ndata_bb_gem_hit_vsigma;   //!
   TBranch        *b_bb_gem_hit_vsigma;   //!
   TBranch        *b_Ndata_bb_gem_hit_vstriphi;   //!
   TBranch        *b_bb_gem_hit_vstriphi;   //!
   TBranch        *b_Ndata_bb_gem_hit_vstriplo;   //!
   TBranch        *b_bb_gem_hit_vstriplo;   //!
   TBranch        *b_Ndata_bb_gem_hit_vstripmax;   //!
   TBranch        *b_bb_gem_hit_vstripmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_xglobal;   //!
   TBranch        *b_bb_gem_hit_xglobal;   //!
   TBranch        *b_Ndata_bb_gem_hit_xlocal;   //!
   TBranch        *b_bb_gem_hit_xlocal;   //!
   TBranch        *b_Ndata_bb_gem_hit_yglobal;   //!
   TBranch        *b_bb_gem_hit_yglobal;   //!
   TBranch        *b_Ndata_bb_gem_hit_ylocal;   //!
   TBranch        *b_bb_gem_hit_ylocal;   //!
   TBranch        *b_Ndata_bb_gem_hit_zglobal;   //!
   TBranch        *b_bb_gem_hit_zglobal;   //!
   TBranch        *b_Ndata_bb_gem_m0_clust_clustu_adc;   //!
   TBranch        *b_bb_gem_m0_clust_clustu_adc;   //!
   TBranch        *b_Ndata_bb_gem_m0_clust_clustu_pos;   //!
   TBranch        *b_bb_gem_m0_clust_clustu_pos;   //!
   TBranch        *b_Ndata_bb_gem_m0_clust_clustu_strips;   //!
   TBranch        *b_bb_gem_m0_clust_clustu_strips;   //!
   TBranch        *b_Ndata_bb_gem_m0_clust_clustu_time;   //!
   TBranch        *b_bb_gem_m0_clust_clustu_time;   //!
   TBranch        *b_Ndata_bb_gem_m0_clust_clustv_adc;   //!
   TBranch        *b_bb_gem_m0_clust_clustv_adc;   //!
   TBranch        *b_Ndata_bb_gem_m0_clust_clustv_pos;   //!
   TBranch        *b_bb_gem_m0_clust_clustv_pos;   //!
   TBranch        *b_Ndata_bb_gem_m0_clust_clustv_strips;   //!
   TBranch        *b_bb_gem_m0_clust_clustv_strips;   //!
   TBranch        *b_Ndata_bb_gem_m0_clust_clustv_time;   //!
   TBranch        *b_bb_gem_m0_clust_clustv_time;   //!
   TBranch        *b_Ndata_bb_gem_m0_clust_isnegativeU;   //!
   TBranch        *b_bb_gem_m0_clust_isnegativeU;   //!
   TBranch        *b_Ndata_bb_gem_m0_clust_isnegativeV;   //!
   TBranch        *b_bb_gem_m0_clust_isnegativeV;   //!
   TBranch        *b_Ndata_bb_gem_m0_clust_isnegontrackU;   //!
   TBranch        *b_bb_gem_m0_clust_isnegontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m0_clust_isnegontrackV;   //!
   TBranch        *b_bb_gem_m0_clust_isnegontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_ADCavg;   //!
   TBranch        *b_bb_gem_m0_strip_ADCavg;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_ADCmax;   //!
   TBranch        *b_bb_gem_m0_strip_ADCmax;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_ADCsamples;   //!
   TBranch        *b_bb_gem_m0_strip_ADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_ADCsum;   //!
   TBranch        *b_bb_gem_m0_strip_ADCsum;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_BUILD_ALL_SAMPLES;   //!
   TBranch        *b_bb_gem_m0_strip_BUILD_ALL_SAMPLES;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_CM_GOOD;   //!
   TBranch        *b_bb_gem_m0_strip_CM_GOOD;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_CorrCoeff;   //!
   TBranch        *b_bb_gem_m0_strip_CorrCoeff;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_DeconvADCmax;   //!
   TBranch        *b_bb_gem_m0_strip_DeconvADCmax;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_DeconvADCmaxCombo;   //!
   TBranch        *b_bb_gem_m0_strip_DeconvADCmaxCombo;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_DeconvADCsamples;   //!
   TBranch        *b_bb_gem_m0_strip_DeconvADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_DeconvADCsum;   //!
   TBranch        *b_bb_gem_m0_strip_DeconvADCsum;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_ENABLE_CM;   //!
   TBranch        *b_bb_gem_m0_strip_ENABLE_CM;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_IsU;   //!
   TBranch        *b_bb_gem_m0_strip_IsU;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_IsV;   //!
   TBranch        *b_bb_gem_m0_strip_IsV;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_TSchi2;   //!
   TBranch        *b_bb_gem_m0_strip_TSchi2;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_TSprob;   //!
   TBranch        *b_bb_gem_m0_strip_TSprob;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_Tcorr;   //!
   TBranch        *b_bb_gem_m0_strip_Tcorr;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_Tdiff;   //!
   TBranch        *b_bb_gem_m0_strip_Tdiff;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_Tfit;   //!
   TBranch        *b_bb_gem_m0_strip_Tfit;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_Tmean;   //!
   TBranch        *b_bb_gem_m0_strip_Tmean;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_TmeanDeconv;   //!
   TBranch        *b_bb_gem_m0_strip_TmeanDeconv;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_Tsigma;   //!
   TBranch        *b_bb_gem_m0_strip_Tsigma;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_adc_id;   //!
   TBranch        *b_bb_gem_m0_strip_adc_id;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_crate;   //!
   TBranch        *b_bb_gem_m0_strip_crate;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_isampmax;   //!
   TBranch        *b_bb_gem_m0_strip_isampmax;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_isampmaxDeconv;   //!
   TBranch        *b_bb_gem_m0_strip_isampmaxDeconv;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_isampmaxDeconvCombo;   //!
   TBranch        *b_bb_gem_m0_strip_isampmaxDeconvCombo;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_isnegative;   //!
   TBranch        *b_bb_gem_m0_strip_isnegative;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_isnegativeU;   //!
   TBranch        *b_bb_gem_m0_strip_isnegativeU;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_isnegativeV;   //!
   TBranch        *b_bb_gem_m0_strip_isnegativeV;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_isnegontrack;   //!
   TBranch        *b_bb_gem_m0_strip_isnegontrack;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_isnegontrackU;   //!
   TBranch        *b_bb_gem_m0_strip_isnegontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_isnegontrackV;   //!
   TBranch        *b_bb_gem_m0_strip_isnegontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_istrip;   //!
   TBranch        *b_bb_gem_m0_strip_istrip;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_itrack;   //!
   TBranch        *b_bb_gem_m0_strip_itrack;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_mpd;   //!
   TBranch        *b_bb_gem_m0_strip_mpd;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_ontrack;   //!
   TBranch        *b_bb_gem_m0_strip_ontrack;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_ontrackU;   //!
   TBranch        *b_bb_gem_m0_strip_ontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_ontrackV;   //!
   TBranch        *b_bb_gem_m0_strip_ontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m0_strip_rawADCsamples;   //!
   TBranch        *b_bb_gem_m0_strip_rawADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m1_clust_clustu_adc;   //!
   TBranch        *b_bb_gem_m1_clust_clustu_adc;   //!
   TBranch        *b_Ndata_bb_gem_m1_clust_clustu_pos;   //!
   TBranch        *b_bb_gem_m1_clust_clustu_pos;   //!
   TBranch        *b_Ndata_bb_gem_m1_clust_clustu_strips;   //!
   TBranch        *b_bb_gem_m1_clust_clustu_strips;   //!
   TBranch        *b_Ndata_bb_gem_m1_clust_clustu_time;   //!
   TBranch        *b_bb_gem_m1_clust_clustu_time;   //!
   TBranch        *b_Ndata_bb_gem_m1_clust_clustv_adc;   //!
   TBranch        *b_bb_gem_m1_clust_clustv_adc;   //!
   TBranch        *b_Ndata_bb_gem_m1_clust_clustv_pos;   //!
   TBranch        *b_bb_gem_m1_clust_clustv_pos;   //!
   TBranch        *b_Ndata_bb_gem_m1_clust_clustv_strips;   //!
   TBranch        *b_bb_gem_m1_clust_clustv_strips;   //!
   TBranch        *b_Ndata_bb_gem_m1_clust_clustv_time;   //!
   TBranch        *b_bb_gem_m1_clust_clustv_time;   //!
   TBranch        *b_Ndata_bb_gem_m1_clust_isnegativeU;   //!
   TBranch        *b_bb_gem_m1_clust_isnegativeU;   //!
   TBranch        *b_Ndata_bb_gem_m1_clust_isnegativeV;   //!
   TBranch        *b_bb_gem_m1_clust_isnegativeV;   //!
   TBranch        *b_Ndata_bb_gem_m1_clust_isnegontrackU;   //!
   TBranch        *b_bb_gem_m1_clust_isnegontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m1_clust_isnegontrackV;   //!
   TBranch        *b_bb_gem_m1_clust_isnegontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_ADCavg;   //!
   TBranch        *b_bb_gem_m1_strip_ADCavg;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_ADCmax;   //!
   TBranch        *b_bb_gem_m1_strip_ADCmax;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_ADCsamples;   //!
   TBranch        *b_bb_gem_m1_strip_ADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_ADCsum;   //!
   TBranch        *b_bb_gem_m1_strip_ADCsum;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_BUILD_ALL_SAMPLES;   //!
   TBranch        *b_bb_gem_m1_strip_BUILD_ALL_SAMPLES;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_CM_GOOD;   //!
   TBranch        *b_bb_gem_m1_strip_CM_GOOD;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_CorrCoeff;   //!
   TBranch        *b_bb_gem_m1_strip_CorrCoeff;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_DeconvADCmax;   //!
   TBranch        *b_bb_gem_m1_strip_DeconvADCmax;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_DeconvADCmaxCombo;   //!
   TBranch        *b_bb_gem_m1_strip_DeconvADCmaxCombo;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_DeconvADCsamples;   //!
   TBranch        *b_bb_gem_m1_strip_DeconvADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_DeconvADCsum;   //!
   TBranch        *b_bb_gem_m1_strip_DeconvADCsum;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_ENABLE_CM;   //!
   TBranch        *b_bb_gem_m1_strip_ENABLE_CM;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_IsU;   //!
   TBranch        *b_bb_gem_m1_strip_IsU;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_IsV;   //!
   TBranch        *b_bb_gem_m1_strip_IsV;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_TSchi2;   //!
   TBranch        *b_bb_gem_m1_strip_TSchi2;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_TSprob;   //!
   TBranch        *b_bb_gem_m1_strip_TSprob;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_Tcorr;   //!
   TBranch        *b_bb_gem_m1_strip_Tcorr;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_Tdiff;   //!
   TBranch        *b_bb_gem_m1_strip_Tdiff;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_Tfit;   //!
   TBranch        *b_bb_gem_m1_strip_Tfit;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_Tmean;   //!
   TBranch        *b_bb_gem_m1_strip_Tmean;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_TmeanDeconv;   //!
   TBranch        *b_bb_gem_m1_strip_TmeanDeconv;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_Tsigma;   //!
   TBranch        *b_bb_gem_m1_strip_Tsigma;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_adc_id;   //!
   TBranch        *b_bb_gem_m1_strip_adc_id;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_crate;   //!
   TBranch        *b_bb_gem_m1_strip_crate;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_isampmax;   //!
   TBranch        *b_bb_gem_m1_strip_isampmax;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_isampmaxDeconv;   //!
   TBranch        *b_bb_gem_m1_strip_isampmaxDeconv;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_isampmaxDeconvCombo;   //!
   TBranch        *b_bb_gem_m1_strip_isampmaxDeconvCombo;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_isnegative;   //!
   TBranch        *b_bb_gem_m1_strip_isnegative;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_isnegativeU;   //!
   TBranch        *b_bb_gem_m1_strip_isnegativeU;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_isnegativeV;   //!
   TBranch        *b_bb_gem_m1_strip_isnegativeV;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_isnegontrack;   //!
   TBranch        *b_bb_gem_m1_strip_isnegontrack;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_isnegontrackU;   //!
   TBranch        *b_bb_gem_m1_strip_isnegontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_isnegontrackV;   //!
   TBranch        *b_bb_gem_m1_strip_isnegontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_istrip;   //!
   TBranch        *b_bb_gem_m1_strip_istrip;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_itrack;   //!
   TBranch        *b_bb_gem_m1_strip_itrack;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_mpd;   //!
   TBranch        *b_bb_gem_m1_strip_mpd;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_ontrack;   //!
   TBranch        *b_bb_gem_m1_strip_ontrack;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_ontrackU;   //!
   TBranch        *b_bb_gem_m1_strip_ontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_ontrackV;   //!
   TBranch        *b_bb_gem_m1_strip_ontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m1_strip_rawADCsamples;   //!
   TBranch        *b_bb_gem_m1_strip_rawADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m2_clust_clustu_adc;   //!
   TBranch        *b_bb_gem_m2_clust_clustu_adc;   //!
   TBranch        *b_Ndata_bb_gem_m2_clust_clustu_pos;   //!
   TBranch        *b_bb_gem_m2_clust_clustu_pos;   //!
   TBranch        *b_Ndata_bb_gem_m2_clust_clustu_strips;   //!
   TBranch        *b_bb_gem_m2_clust_clustu_strips;   //!
   TBranch        *b_Ndata_bb_gem_m2_clust_clustu_time;   //!
   TBranch        *b_bb_gem_m2_clust_clustu_time;   //!
   TBranch        *b_Ndata_bb_gem_m2_clust_clustv_adc;   //!
   TBranch        *b_bb_gem_m2_clust_clustv_adc;   //!
   TBranch        *b_Ndata_bb_gem_m2_clust_clustv_pos;   //!
   TBranch        *b_bb_gem_m2_clust_clustv_pos;   //!
   TBranch        *b_Ndata_bb_gem_m2_clust_clustv_strips;   //!
   TBranch        *b_bb_gem_m2_clust_clustv_strips;   //!
   TBranch        *b_Ndata_bb_gem_m2_clust_clustv_time;   //!
   TBranch        *b_bb_gem_m2_clust_clustv_time;   //!
   TBranch        *b_Ndata_bb_gem_m2_clust_isnegativeU;   //!
   TBranch        *b_bb_gem_m2_clust_isnegativeU;   //!
   TBranch        *b_Ndata_bb_gem_m2_clust_isnegativeV;   //!
   TBranch        *b_bb_gem_m2_clust_isnegativeV;   //!
   TBranch        *b_Ndata_bb_gem_m2_clust_isnegontrackU;   //!
   TBranch        *b_bb_gem_m2_clust_isnegontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m2_clust_isnegontrackV;   //!
   TBranch        *b_bb_gem_m2_clust_isnegontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_ADCavg;   //!
   TBranch        *b_bb_gem_m2_strip_ADCavg;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_ADCmax;   //!
   TBranch        *b_bb_gem_m2_strip_ADCmax;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_ADCsamples;   //!
   TBranch        *b_bb_gem_m2_strip_ADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_ADCsum;   //!
   TBranch        *b_bb_gem_m2_strip_ADCsum;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_BUILD_ALL_SAMPLES;   //!
   TBranch        *b_bb_gem_m2_strip_BUILD_ALL_SAMPLES;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_CM_GOOD;   //!
   TBranch        *b_bb_gem_m2_strip_CM_GOOD;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_CorrCoeff;   //!
   TBranch        *b_bb_gem_m2_strip_CorrCoeff;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_DeconvADCmax;   //!
   TBranch        *b_bb_gem_m2_strip_DeconvADCmax;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_DeconvADCmaxCombo;   //!
   TBranch        *b_bb_gem_m2_strip_DeconvADCmaxCombo;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_DeconvADCsamples;   //!
   TBranch        *b_bb_gem_m2_strip_DeconvADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_DeconvADCsum;   //!
   TBranch        *b_bb_gem_m2_strip_DeconvADCsum;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_ENABLE_CM;   //!
   TBranch        *b_bb_gem_m2_strip_ENABLE_CM;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_IsU;   //!
   TBranch        *b_bb_gem_m2_strip_IsU;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_IsV;   //!
   TBranch        *b_bb_gem_m2_strip_IsV;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_TSchi2;   //!
   TBranch        *b_bb_gem_m2_strip_TSchi2;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_TSprob;   //!
   TBranch        *b_bb_gem_m2_strip_TSprob;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_Tcorr;   //!
   TBranch        *b_bb_gem_m2_strip_Tcorr;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_Tdiff;   //!
   TBranch        *b_bb_gem_m2_strip_Tdiff;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_Tfit;   //!
   TBranch        *b_bb_gem_m2_strip_Tfit;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_Tmean;   //!
   TBranch        *b_bb_gem_m2_strip_Tmean;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_TmeanDeconv;   //!
   TBranch        *b_bb_gem_m2_strip_TmeanDeconv;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_Tsigma;   //!
   TBranch        *b_bb_gem_m2_strip_Tsigma;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_adc_id;   //!
   TBranch        *b_bb_gem_m2_strip_adc_id;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_crate;   //!
   TBranch        *b_bb_gem_m2_strip_crate;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_isampmax;   //!
   TBranch        *b_bb_gem_m2_strip_isampmax;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_isampmaxDeconv;   //!
   TBranch        *b_bb_gem_m2_strip_isampmaxDeconv;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_isampmaxDeconvCombo;   //!
   TBranch        *b_bb_gem_m2_strip_isampmaxDeconvCombo;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_isnegative;   //!
   TBranch        *b_bb_gem_m2_strip_isnegative;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_isnegativeU;   //!
   TBranch        *b_bb_gem_m2_strip_isnegativeU;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_isnegativeV;   //!
   TBranch        *b_bb_gem_m2_strip_isnegativeV;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_isnegontrack;   //!
   TBranch        *b_bb_gem_m2_strip_isnegontrack;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_isnegontrackU;   //!
   TBranch        *b_bb_gem_m2_strip_isnegontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_isnegontrackV;   //!
   TBranch        *b_bb_gem_m2_strip_isnegontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_istrip;   //!
   TBranch        *b_bb_gem_m2_strip_istrip;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_itrack;   //!
   TBranch        *b_bb_gem_m2_strip_itrack;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_mpd;   //!
   TBranch        *b_bb_gem_m2_strip_mpd;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_ontrack;   //!
   TBranch        *b_bb_gem_m2_strip_ontrack;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_ontrackU;   //!
   TBranch        *b_bb_gem_m2_strip_ontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_ontrackV;   //!
   TBranch        *b_bb_gem_m2_strip_ontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m2_strip_rawADCsamples;   //!
   TBranch        *b_bb_gem_m2_strip_rawADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m3_clust_clustu_adc;   //!
   TBranch        *b_bb_gem_m3_clust_clustu_adc;   //!
   TBranch        *b_Ndata_bb_gem_m3_clust_clustu_pos;   //!
   TBranch        *b_bb_gem_m3_clust_clustu_pos;   //!
   TBranch        *b_Ndata_bb_gem_m3_clust_clustu_strips;   //!
   TBranch        *b_bb_gem_m3_clust_clustu_strips;   //!
   TBranch        *b_Ndata_bb_gem_m3_clust_clustu_time;   //!
   TBranch        *b_bb_gem_m3_clust_clustu_time;   //!
   TBranch        *b_Ndata_bb_gem_m3_clust_clustv_adc;   //!
   TBranch        *b_bb_gem_m3_clust_clustv_adc;   //!
   TBranch        *b_Ndata_bb_gem_m3_clust_clustv_pos;   //!
   TBranch        *b_bb_gem_m3_clust_clustv_pos;   //!
   TBranch        *b_Ndata_bb_gem_m3_clust_clustv_strips;   //!
   TBranch        *b_bb_gem_m3_clust_clustv_strips;   //!
   TBranch        *b_Ndata_bb_gem_m3_clust_clustv_time;   //!
   TBranch        *b_bb_gem_m3_clust_clustv_time;   //!
   TBranch        *b_Ndata_bb_gem_m3_clust_isnegativeU;   //!
   TBranch        *b_bb_gem_m3_clust_isnegativeU;   //!
   TBranch        *b_Ndata_bb_gem_m3_clust_isnegativeV;   //!
   TBranch        *b_bb_gem_m3_clust_isnegativeV;   //!
   TBranch        *b_Ndata_bb_gem_m3_clust_isnegontrackU;   //!
   TBranch        *b_bb_gem_m3_clust_isnegontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m3_clust_isnegontrackV;   //!
   TBranch        *b_bb_gem_m3_clust_isnegontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_ADCavg;   //!
   TBranch        *b_bb_gem_m3_strip_ADCavg;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_ADCmax;   //!
   TBranch        *b_bb_gem_m3_strip_ADCmax;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_ADCsamples;   //!
   TBranch        *b_bb_gem_m3_strip_ADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_ADCsum;   //!
   TBranch        *b_bb_gem_m3_strip_ADCsum;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_BUILD_ALL_SAMPLES;   //!
   TBranch        *b_bb_gem_m3_strip_BUILD_ALL_SAMPLES;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_CM_GOOD;   //!
   TBranch        *b_bb_gem_m3_strip_CM_GOOD;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_CorrCoeff;   //!
   TBranch        *b_bb_gem_m3_strip_CorrCoeff;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_DeconvADCmax;   //!
   TBranch        *b_bb_gem_m3_strip_DeconvADCmax;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_DeconvADCmaxCombo;   //!
   TBranch        *b_bb_gem_m3_strip_DeconvADCmaxCombo;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_DeconvADCsamples;   //!
   TBranch        *b_bb_gem_m3_strip_DeconvADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_DeconvADCsum;   //!
   TBranch        *b_bb_gem_m3_strip_DeconvADCsum;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_ENABLE_CM;   //!
   TBranch        *b_bb_gem_m3_strip_ENABLE_CM;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_IsU;   //!
   TBranch        *b_bb_gem_m3_strip_IsU;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_IsV;   //!
   TBranch        *b_bb_gem_m3_strip_IsV;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_TSchi2;   //!
   TBranch        *b_bb_gem_m3_strip_TSchi2;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_TSprob;   //!
   TBranch        *b_bb_gem_m3_strip_TSprob;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_Tcorr;   //!
   TBranch        *b_bb_gem_m3_strip_Tcorr;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_Tdiff;   //!
   TBranch        *b_bb_gem_m3_strip_Tdiff;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_Tfit;   //!
   TBranch        *b_bb_gem_m3_strip_Tfit;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_Tmean;   //!
   TBranch        *b_bb_gem_m3_strip_Tmean;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_TmeanDeconv;   //!
   TBranch        *b_bb_gem_m3_strip_TmeanDeconv;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_Tsigma;   //!
   TBranch        *b_bb_gem_m3_strip_Tsigma;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_adc_id;   //!
   TBranch        *b_bb_gem_m3_strip_adc_id;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_crate;   //!
   TBranch        *b_bb_gem_m3_strip_crate;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_isampmax;   //!
   TBranch        *b_bb_gem_m3_strip_isampmax;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_isampmaxDeconv;   //!
   TBranch        *b_bb_gem_m3_strip_isampmaxDeconv;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_isampmaxDeconvCombo;   //!
   TBranch        *b_bb_gem_m3_strip_isampmaxDeconvCombo;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_isnegative;   //!
   TBranch        *b_bb_gem_m3_strip_isnegative;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_isnegativeU;   //!
   TBranch        *b_bb_gem_m3_strip_isnegativeU;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_isnegativeV;   //!
   TBranch        *b_bb_gem_m3_strip_isnegativeV;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_isnegontrack;   //!
   TBranch        *b_bb_gem_m3_strip_isnegontrack;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_isnegontrackU;   //!
   TBranch        *b_bb_gem_m3_strip_isnegontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_isnegontrackV;   //!
   TBranch        *b_bb_gem_m3_strip_isnegontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_istrip;   //!
   TBranch        *b_bb_gem_m3_strip_istrip;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_itrack;   //!
   TBranch        *b_bb_gem_m3_strip_itrack;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_mpd;   //!
   TBranch        *b_bb_gem_m3_strip_mpd;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_ontrack;   //!
   TBranch        *b_bb_gem_m3_strip_ontrack;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_ontrackU;   //!
   TBranch        *b_bb_gem_m3_strip_ontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_ontrackV;   //!
   TBranch        *b_bb_gem_m3_strip_ontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m3_strip_rawADCsamples;   //!
   TBranch        *b_bb_gem_m3_strip_rawADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m4_clust_clustu_adc;   //!
   TBranch        *b_bb_gem_m4_clust_clustu_adc;   //!
   TBranch        *b_Ndata_bb_gem_m4_clust_clustu_pos;   //!
   TBranch        *b_bb_gem_m4_clust_clustu_pos;   //!
   TBranch        *b_Ndata_bb_gem_m4_clust_clustu_strips;   //!
   TBranch        *b_bb_gem_m4_clust_clustu_strips;   //!
   TBranch        *b_Ndata_bb_gem_m4_clust_clustu_time;   //!
   TBranch        *b_bb_gem_m4_clust_clustu_time;   //!
   TBranch        *b_Ndata_bb_gem_m4_clust_clustv_adc;   //!
   TBranch        *b_bb_gem_m4_clust_clustv_adc;   //!
   TBranch        *b_Ndata_bb_gem_m4_clust_clustv_pos;   //!
   TBranch        *b_bb_gem_m4_clust_clustv_pos;   //!
   TBranch        *b_Ndata_bb_gem_m4_clust_clustv_strips;   //!
   TBranch        *b_bb_gem_m4_clust_clustv_strips;   //!
   TBranch        *b_Ndata_bb_gem_m4_clust_clustv_time;   //!
   TBranch        *b_bb_gem_m4_clust_clustv_time;   //!
   TBranch        *b_Ndata_bb_gem_m4_clust_isnegativeU;   //!
   TBranch        *b_bb_gem_m4_clust_isnegativeU;   //!
   TBranch        *b_Ndata_bb_gem_m4_clust_isnegativeV;   //!
   TBranch        *b_bb_gem_m4_clust_isnegativeV;   //!
   TBranch        *b_Ndata_bb_gem_m4_clust_isnegontrackU;   //!
   TBranch        *b_bb_gem_m4_clust_isnegontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m4_clust_isnegontrackV;   //!
   TBranch        *b_bb_gem_m4_clust_isnegontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_ADCavg;   //!
   TBranch        *b_bb_gem_m4_strip_ADCavg;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_ADCmax;   //!
   TBranch        *b_bb_gem_m4_strip_ADCmax;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_ADCsamples;   //!
   TBranch        *b_bb_gem_m4_strip_ADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_ADCsum;   //!
   TBranch        *b_bb_gem_m4_strip_ADCsum;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_BUILD_ALL_SAMPLES;   //!
   TBranch        *b_bb_gem_m4_strip_BUILD_ALL_SAMPLES;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_CM_GOOD;   //!
   TBranch        *b_bb_gem_m4_strip_CM_GOOD;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_CorrCoeff;   //!
   TBranch        *b_bb_gem_m4_strip_CorrCoeff;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_DeconvADCmax;   //!
   TBranch        *b_bb_gem_m4_strip_DeconvADCmax;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_DeconvADCmaxCombo;   //!
   TBranch        *b_bb_gem_m4_strip_DeconvADCmaxCombo;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_DeconvADCsamples;   //!
   TBranch        *b_bb_gem_m4_strip_DeconvADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_DeconvADCsum;   //!
   TBranch        *b_bb_gem_m4_strip_DeconvADCsum;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_ENABLE_CM;   //!
   TBranch        *b_bb_gem_m4_strip_ENABLE_CM;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_IsU;   //!
   TBranch        *b_bb_gem_m4_strip_IsU;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_IsV;   //!
   TBranch        *b_bb_gem_m4_strip_IsV;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_TSchi2;   //!
   TBranch        *b_bb_gem_m4_strip_TSchi2;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_TSprob;   //!
   TBranch        *b_bb_gem_m4_strip_TSprob;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_Tcorr;   //!
   TBranch        *b_bb_gem_m4_strip_Tcorr;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_Tdiff;   //!
   TBranch        *b_bb_gem_m4_strip_Tdiff;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_Tfit;   //!
   TBranch        *b_bb_gem_m4_strip_Tfit;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_Tmean;   //!
   TBranch        *b_bb_gem_m4_strip_Tmean;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_TmeanDeconv;   //!
   TBranch        *b_bb_gem_m4_strip_TmeanDeconv;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_Tsigma;   //!
   TBranch        *b_bb_gem_m4_strip_Tsigma;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_adc_id;   //!
   TBranch        *b_bb_gem_m4_strip_adc_id;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_crate;   //!
   TBranch        *b_bb_gem_m4_strip_crate;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_isampmax;   //!
   TBranch        *b_bb_gem_m4_strip_isampmax;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_isampmaxDeconv;   //!
   TBranch        *b_bb_gem_m4_strip_isampmaxDeconv;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_isampmaxDeconvCombo;   //!
   TBranch        *b_bb_gem_m4_strip_isampmaxDeconvCombo;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_isnegative;   //!
   TBranch        *b_bb_gem_m4_strip_isnegative;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_isnegativeU;   //!
   TBranch        *b_bb_gem_m4_strip_isnegativeU;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_isnegativeV;   //!
   TBranch        *b_bb_gem_m4_strip_isnegativeV;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_isnegontrack;   //!
   TBranch        *b_bb_gem_m4_strip_isnegontrack;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_isnegontrackU;   //!
   TBranch        *b_bb_gem_m4_strip_isnegontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_isnegontrackV;   //!
   TBranch        *b_bb_gem_m4_strip_isnegontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_istrip;   //!
   TBranch        *b_bb_gem_m4_strip_istrip;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_itrack;   //!
   TBranch        *b_bb_gem_m4_strip_itrack;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_mpd;   //!
   TBranch        *b_bb_gem_m4_strip_mpd;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_ontrack;   //!
   TBranch        *b_bb_gem_m4_strip_ontrack;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_ontrackU;   //!
   TBranch        *b_bb_gem_m4_strip_ontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_ontrackV;   //!
   TBranch        *b_bb_gem_m4_strip_ontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m4_strip_rawADCsamples;   //!
   TBranch        *b_bb_gem_m4_strip_rawADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m5_clust_clustu_adc;   //!
   TBranch        *b_bb_gem_m5_clust_clustu_adc;   //!
   TBranch        *b_Ndata_bb_gem_m5_clust_clustu_pos;   //!
   TBranch        *b_bb_gem_m5_clust_clustu_pos;   //!
   TBranch        *b_Ndata_bb_gem_m5_clust_clustu_strips;   //!
   TBranch        *b_bb_gem_m5_clust_clustu_strips;   //!
   TBranch        *b_Ndata_bb_gem_m5_clust_clustu_time;   //!
   TBranch        *b_bb_gem_m5_clust_clustu_time;   //!
   TBranch        *b_Ndata_bb_gem_m5_clust_clustv_adc;   //!
   TBranch        *b_bb_gem_m5_clust_clustv_adc;   //!
   TBranch        *b_Ndata_bb_gem_m5_clust_clustv_pos;   //!
   TBranch        *b_bb_gem_m5_clust_clustv_pos;   //!
   TBranch        *b_Ndata_bb_gem_m5_clust_clustv_strips;   //!
   TBranch        *b_bb_gem_m5_clust_clustv_strips;   //!
   TBranch        *b_Ndata_bb_gem_m5_clust_clustv_time;   //!
   TBranch        *b_bb_gem_m5_clust_clustv_time;   //!
   TBranch        *b_Ndata_bb_gem_m5_clust_isnegativeU;   //!
   TBranch        *b_bb_gem_m5_clust_isnegativeU;   //!
   TBranch        *b_Ndata_bb_gem_m5_clust_isnegativeV;   //!
   TBranch        *b_bb_gem_m5_clust_isnegativeV;   //!
   TBranch        *b_Ndata_bb_gem_m5_clust_isnegontrackU;   //!
   TBranch        *b_bb_gem_m5_clust_isnegontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m5_clust_isnegontrackV;   //!
   TBranch        *b_bb_gem_m5_clust_isnegontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_ADCavg;   //!
   TBranch        *b_bb_gem_m5_strip_ADCavg;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_ADCmax;   //!
   TBranch        *b_bb_gem_m5_strip_ADCmax;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_ADCsamples;   //!
   TBranch        *b_bb_gem_m5_strip_ADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_ADCsum;   //!
   TBranch        *b_bb_gem_m5_strip_ADCsum;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_BUILD_ALL_SAMPLES;   //!
   TBranch        *b_bb_gem_m5_strip_BUILD_ALL_SAMPLES;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_CM_GOOD;   //!
   TBranch        *b_bb_gem_m5_strip_CM_GOOD;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_CorrCoeff;   //!
   TBranch        *b_bb_gem_m5_strip_CorrCoeff;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_DeconvADCmax;   //!
   TBranch        *b_bb_gem_m5_strip_DeconvADCmax;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_DeconvADCmaxCombo;   //!
   TBranch        *b_bb_gem_m5_strip_DeconvADCmaxCombo;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_DeconvADCsamples;   //!
   TBranch        *b_bb_gem_m5_strip_DeconvADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_DeconvADCsum;   //!
   TBranch        *b_bb_gem_m5_strip_DeconvADCsum;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_ENABLE_CM;   //!
   TBranch        *b_bb_gem_m5_strip_ENABLE_CM;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_IsU;   //!
   TBranch        *b_bb_gem_m5_strip_IsU;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_IsV;   //!
   TBranch        *b_bb_gem_m5_strip_IsV;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_TSchi2;   //!
   TBranch        *b_bb_gem_m5_strip_TSchi2;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_TSprob;   //!
   TBranch        *b_bb_gem_m5_strip_TSprob;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_Tcorr;   //!
   TBranch        *b_bb_gem_m5_strip_Tcorr;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_Tdiff;   //!
   TBranch        *b_bb_gem_m5_strip_Tdiff;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_Tfit;   //!
   TBranch        *b_bb_gem_m5_strip_Tfit;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_Tmean;   //!
   TBranch        *b_bb_gem_m5_strip_Tmean;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_TmeanDeconv;   //!
   TBranch        *b_bb_gem_m5_strip_TmeanDeconv;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_Tsigma;   //!
   TBranch        *b_bb_gem_m5_strip_Tsigma;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_adc_id;   //!
   TBranch        *b_bb_gem_m5_strip_adc_id;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_crate;   //!
   TBranch        *b_bb_gem_m5_strip_crate;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_isampmax;   //!
   TBranch        *b_bb_gem_m5_strip_isampmax;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_isampmaxDeconv;   //!
   TBranch        *b_bb_gem_m5_strip_isampmaxDeconv;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_isampmaxDeconvCombo;   //!
   TBranch        *b_bb_gem_m5_strip_isampmaxDeconvCombo;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_isnegative;   //!
   TBranch        *b_bb_gem_m5_strip_isnegative;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_isnegativeU;   //!
   TBranch        *b_bb_gem_m5_strip_isnegativeU;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_isnegativeV;   //!
   TBranch        *b_bb_gem_m5_strip_isnegativeV;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_isnegontrack;   //!
   TBranch        *b_bb_gem_m5_strip_isnegontrack;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_isnegontrackU;   //!
   TBranch        *b_bb_gem_m5_strip_isnegontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_isnegontrackV;   //!
   TBranch        *b_bb_gem_m5_strip_isnegontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_istrip;   //!
   TBranch        *b_bb_gem_m5_strip_istrip;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_itrack;   //!
   TBranch        *b_bb_gem_m5_strip_itrack;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_mpd;   //!
   TBranch        *b_bb_gem_m5_strip_mpd;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_ontrack;   //!
   TBranch        *b_bb_gem_m5_strip_ontrack;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_ontrackU;   //!
   TBranch        *b_bb_gem_m5_strip_ontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_ontrackV;   //!
   TBranch        *b_bb_gem_m5_strip_ontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m5_strip_rawADCsamples;   //!
   TBranch        *b_bb_gem_m5_strip_rawADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m6_clust_clustu_adc;   //!
   TBranch        *b_bb_gem_m6_clust_clustu_adc;   //!
   TBranch        *b_Ndata_bb_gem_m6_clust_clustu_pos;   //!
   TBranch        *b_bb_gem_m6_clust_clustu_pos;   //!
   TBranch        *b_Ndata_bb_gem_m6_clust_clustu_strips;   //!
   TBranch        *b_bb_gem_m6_clust_clustu_strips;   //!
   TBranch        *b_Ndata_bb_gem_m6_clust_clustu_time;   //!
   TBranch        *b_bb_gem_m6_clust_clustu_time;   //!
   TBranch        *b_Ndata_bb_gem_m6_clust_clustv_adc;   //!
   TBranch        *b_bb_gem_m6_clust_clustv_adc;   //!
   TBranch        *b_Ndata_bb_gem_m6_clust_clustv_pos;   //!
   TBranch        *b_bb_gem_m6_clust_clustv_pos;   //!
   TBranch        *b_Ndata_bb_gem_m6_clust_clustv_strips;   //!
   TBranch        *b_bb_gem_m6_clust_clustv_strips;   //!
   TBranch        *b_Ndata_bb_gem_m6_clust_clustv_time;   //!
   TBranch        *b_bb_gem_m6_clust_clustv_time;   //!
   TBranch        *b_Ndata_bb_gem_m6_clust_isnegativeU;   //!
   TBranch        *b_bb_gem_m6_clust_isnegativeU;   //!
   TBranch        *b_Ndata_bb_gem_m6_clust_isnegativeV;   //!
   TBranch        *b_bb_gem_m6_clust_isnegativeV;   //!
   TBranch        *b_Ndata_bb_gem_m6_clust_isnegontrackU;   //!
   TBranch        *b_bb_gem_m6_clust_isnegontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m6_clust_isnegontrackV;   //!
   TBranch        *b_bb_gem_m6_clust_isnegontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_ADCavg;   //!
   TBranch        *b_bb_gem_m6_strip_ADCavg;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_ADCmax;   //!
   TBranch        *b_bb_gem_m6_strip_ADCmax;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_ADCsamples;   //!
   TBranch        *b_bb_gem_m6_strip_ADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_ADCsum;   //!
   TBranch        *b_bb_gem_m6_strip_ADCsum;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_BUILD_ALL_SAMPLES;   //!
   TBranch        *b_bb_gem_m6_strip_BUILD_ALL_SAMPLES;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_CM_GOOD;   //!
   TBranch        *b_bb_gem_m6_strip_CM_GOOD;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_CorrCoeff;   //!
   TBranch        *b_bb_gem_m6_strip_CorrCoeff;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_DeconvADCmax;   //!
   TBranch        *b_bb_gem_m6_strip_DeconvADCmax;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_DeconvADCmaxCombo;   //!
   TBranch        *b_bb_gem_m6_strip_DeconvADCmaxCombo;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_DeconvADCsamples;   //!
   TBranch        *b_bb_gem_m6_strip_DeconvADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_DeconvADCsum;   //!
   TBranch        *b_bb_gem_m6_strip_DeconvADCsum;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_ENABLE_CM;   //!
   TBranch        *b_bb_gem_m6_strip_ENABLE_CM;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_IsU;   //!
   TBranch        *b_bb_gem_m6_strip_IsU;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_IsV;   //!
   TBranch        *b_bb_gem_m6_strip_IsV;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_TSchi2;   //!
   TBranch        *b_bb_gem_m6_strip_TSchi2;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_TSprob;   //!
   TBranch        *b_bb_gem_m6_strip_TSprob;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_Tcorr;   //!
   TBranch        *b_bb_gem_m6_strip_Tcorr;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_Tdiff;   //!
   TBranch        *b_bb_gem_m6_strip_Tdiff;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_Tfit;   //!
   TBranch        *b_bb_gem_m6_strip_Tfit;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_Tmean;   //!
   TBranch        *b_bb_gem_m6_strip_Tmean;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_TmeanDeconv;   //!
   TBranch        *b_bb_gem_m6_strip_TmeanDeconv;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_Tsigma;   //!
   TBranch        *b_bb_gem_m6_strip_Tsigma;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_adc_id;   //!
   TBranch        *b_bb_gem_m6_strip_adc_id;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_crate;   //!
   TBranch        *b_bb_gem_m6_strip_crate;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_isampmax;   //!
   TBranch        *b_bb_gem_m6_strip_isampmax;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_isampmaxDeconv;   //!
   TBranch        *b_bb_gem_m6_strip_isampmaxDeconv;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_isampmaxDeconvCombo;   //!
   TBranch        *b_bb_gem_m6_strip_isampmaxDeconvCombo;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_isnegative;   //!
   TBranch        *b_bb_gem_m6_strip_isnegative;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_isnegativeU;   //!
   TBranch        *b_bb_gem_m6_strip_isnegativeU;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_isnegativeV;   //!
   TBranch        *b_bb_gem_m6_strip_isnegativeV;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_isnegontrack;   //!
   TBranch        *b_bb_gem_m6_strip_isnegontrack;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_isnegontrackU;   //!
   TBranch        *b_bb_gem_m6_strip_isnegontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_isnegontrackV;   //!
   TBranch        *b_bb_gem_m6_strip_isnegontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_istrip;   //!
   TBranch        *b_bb_gem_m6_strip_istrip;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_itrack;   //!
   TBranch        *b_bb_gem_m6_strip_itrack;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_mpd;   //!
   TBranch        *b_bb_gem_m6_strip_mpd;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_ontrack;   //!
   TBranch        *b_bb_gem_m6_strip_ontrack;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_ontrackU;   //!
   TBranch        *b_bb_gem_m6_strip_ontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_ontrackV;   //!
   TBranch        *b_bb_gem_m6_strip_ontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m6_strip_rawADCsamples;   //!
   TBranch        *b_bb_gem_m6_strip_rawADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m7_clust_clustu_adc;   //!
   TBranch        *b_bb_gem_m7_clust_clustu_adc;   //!
   TBranch        *b_Ndata_bb_gem_m7_clust_clustu_pos;   //!
   TBranch        *b_bb_gem_m7_clust_clustu_pos;   //!
   TBranch        *b_Ndata_bb_gem_m7_clust_clustu_strips;   //!
   TBranch        *b_bb_gem_m7_clust_clustu_strips;   //!
   TBranch        *b_Ndata_bb_gem_m7_clust_clustu_time;   //!
   TBranch        *b_bb_gem_m7_clust_clustu_time;   //!
   TBranch        *b_Ndata_bb_gem_m7_clust_clustv_adc;   //!
   TBranch        *b_bb_gem_m7_clust_clustv_adc;   //!
   TBranch        *b_Ndata_bb_gem_m7_clust_clustv_pos;   //!
   TBranch        *b_bb_gem_m7_clust_clustv_pos;   //!
   TBranch        *b_Ndata_bb_gem_m7_clust_clustv_strips;   //!
   TBranch        *b_bb_gem_m7_clust_clustv_strips;   //!
   TBranch        *b_Ndata_bb_gem_m7_clust_clustv_time;   //!
   TBranch        *b_bb_gem_m7_clust_clustv_time;   //!
   TBranch        *b_Ndata_bb_gem_m7_clust_isnegativeU;   //!
   TBranch        *b_bb_gem_m7_clust_isnegativeU;   //!
   TBranch        *b_Ndata_bb_gem_m7_clust_isnegativeV;   //!
   TBranch        *b_bb_gem_m7_clust_isnegativeV;   //!
   TBranch        *b_Ndata_bb_gem_m7_clust_isnegontrackU;   //!
   TBranch        *b_bb_gem_m7_clust_isnegontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m7_clust_isnegontrackV;   //!
   TBranch        *b_bb_gem_m7_clust_isnegontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_ADCavg;   //!
   TBranch        *b_bb_gem_m7_strip_ADCavg;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_ADCmax;   //!
   TBranch        *b_bb_gem_m7_strip_ADCmax;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_ADCsamples;   //!
   TBranch        *b_bb_gem_m7_strip_ADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_ADCsum;   //!
   TBranch        *b_bb_gem_m7_strip_ADCsum;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_BUILD_ALL_SAMPLES;   //!
   TBranch        *b_bb_gem_m7_strip_BUILD_ALL_SAMPLES;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_CM_GOOD;   //!
   TBranch        *b_bb_gem_m7_strip_CM_GOOD;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_CorrCoeff;   //!
   TBranch        *b_bb_gem_m7_strip_CorrCoeff;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_DeconvADCmax;   //!
   TBranch        *b_bb_gem_m7_strip_DeconvADCmax;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_DeconvADCmaxCombo;   //!
   TBranch        *b_bb_gem_m7_strip_DeconvADCmaxCombo;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_DeconvADCsamples;   //!
   TBranch        *b_bb_gem_m7_strip_DeconvADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_DeconvADCsum;   //!
   TBranch        *b_bb_gem_m7_strip_DeconvADCsum;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_ENABLE_CM;   //!
   TBranch        *b_bb_gem_m7_strip_ENABLE_CM;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_IsU;   //!
   TBranch        *b_bb_gem_m7_strip_IsU;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_IsV;   //!
   TBranch        *b_bb_gem_m7_strip_IsV;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_TSchi2;   //!
   TBranch        *b_bb_gem_m7_strip_TSchi2;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_TSprob;   //!
   TBranch        *b_bb_gem_m7_strip_TSprob;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_Tcorr;   //!
   TBranch        *b_bb_gem_m7_strip_Tcorr;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_Tdiff;   //!
   TBranch        *b_bb_gem_m7_strip_Tdiff;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_Tfit;   //!
   TBranch        *b_bb_gem_m7_strip_Tfit;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_Tmean;   //!
   TBranch        *b_bb_gem_m7_strip_Tmean;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_TmeanDeconv;   //!
   TBranch        *b_bb_gem_m7_strip_TmeanDeconv;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_Tsigma;   //!
   TBranch        *b_bb_gem_m7_strip_Tsigma;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_adc_id;   //!
   TBranch        *b_bb_gem_m7_strip_adc_id;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_crate;   //!
   TBranch        *b_bb_gem_m7_strip_crate;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_isampmax;   //!
   TBranch        *b_bb_gem_m7_strip_isampmax;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_isampmaxDeconv;   //!
   TBranch        *b_bb_gem_m7_strip_isampmaxDeconv;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_isampmaxDeconvCombo;   //!
   TBranch        *b_bb_gem_m7_strip_isampmaxDeconvCombo;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_isnegative;   //!
   TBranch        *b_bb_gem_m7_strip_isnegative;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_isnegativeU;   //!
   TBranch        *b_bb_gem_m7_strip_isnegativeU;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_isnegativeV;   //!
   TBranch        *b_bb_gem_m7_strip_isnegativeV;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_isnegontrack;   //!
   TBranch        *b_bb_gem_m7_strip_isnegontrack;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_isnegontrackU;   //!
   TBranch        *b_bb_gem_m7_strip_isnegontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_isnegontrackV;   //!
   TBranch        *b_bb_gem_m7_strip_isnegontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_istrip;   //!
   TBranch        *b_bb_gem_m7_strip_istrip;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_itrack;   //!
   TBranch        *b_bb_gem_m7_strip_itrack;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_mpd;   //!
   TBranch        *b_bb_gem_m7_strip_mpd;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_ontrack;   //!
   TBranch        *b_bb_gem_m7_strip_ontrack;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_ontrackU;   //!
   TBranch        *b_bb_gem_m7_strip_ontrackU;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_ontrackV;   //!
   TBranch        *b_bb_gem_m7_strip_ontrackV;   //!
   TBranch        *b_Ndata_bb_gem_m7_strip_rawADCsamples;   //!
   TBranch        *b_bb_gem_m7_strip_rawADCsamples;   //!
   TBranch        *b_Ndata_bb_gem_n2Dhit_layer;   //!
   TBranch        *b_bb_gem_n2Dhit_layer;   //!
   TBranch        *b_Ndata_bb_gem_nclustu_layer;   //!
   TBranch        *b_bb_gem_nclustu_layer;   //!
   TBranch        *b_Ndata_bb_gem_nclustv_layer;   //!
   TBranch        *b_bb_gem_nclustv_layer;   //!
   TBranch        *b_Ndata_bb_gem_nstripsu_layer;   //!
   TBranch        *b_bb_gem_nstripsu_layer;   //!
   TBranch        *b_Ndata_bb_gem_nstripsv_layer;   //!
   TBranch        *b_bb_gem_nstripsv_layer;   //!
   TBranch        *b_Ndata_bb_gem_track_chi2ndf;   //!
   TBranch        *b_bb_gem_track_chi2ndf;   //!
   TBranch        *b_Ndata_bb_gem_track_chi2ndf_hitquality;   //!
   TBranch        *b_bb_gem_track_chi2ndf_hitquality;   //!
   TBranch        *b_Ndata_bb_gem_track_ngoodhits;   //!
   TBranch        *b_bb_gem_track_ngoodhits;   //!
   TBranch        *b_Ndata_bb_gem_track_nhits;   //!
   TBranch        *b_bb_gem_track_nhits;   //!
   TBranch        *b_Ndata_bb_gem_track_t0;   //!
   TBranch        *b_bb_gem_track_t0;   //!
   TBranch        *b_Ndata_bb_gem_track_x;   //!
   TBranch        *b_bb_gem_track_x;   //!
   TBranch        *b_Ndata_bb_gem_track_xp;   //!
   TBranch        *b_bb_gem_track_xp;   //!
   TBranch        *b_Ndata_bb_gem_track_y;   //!
   TBranch        *b_bb_gem_track_y;   //!
   TBranch        *b_Ndata_bb_gem_track_yp;   //!
   TBranch        *b_bb_gem_track_yp;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_allclus_adc;   //!
   TBranch        *b_bb_grinch_tdc_allclus_adc;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_allclus_dx;   //!
   TBranch        *b_bb_grinch_tdc_allclus_dx;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_allclus_dy;   //!
   TBranch        *b_bb_grinch_tdc_allclus_dy;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_allclus_mirrorindex;   //!
   TBranch        *b_bb_grinch_tdc_allclus_mirrorindex;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_allclus_size;   //!
   TBranch        *b_bb_grinch_tdc_allclus_size;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_allclus_t_mean;   //!
   TBranch        *b_bb_grinch_tdc_allclus_t_mean;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_allclus_t_rms;   //!
   TBranch        *b_bb_grinch_tdc_allclus_t_rms;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_allclus_tot_mean;   //!
   TBranch        *b_bb_grinch_tdc_allclus_tot_mean;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_allclus_trackindex;   //!
   TBranch        *b_bb_grinch_tdc_allclus_trackindex;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_allclus_x_mean;   //!
   TBranch        *b_bb_grinch_tdc_allclus_x_mean;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_allclus_y_mean;   //!
   TBranch        *b_bb_grinch_tdc_allclus_y_mean;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_hit_amp;   //!
   TBranch        *b_bb_grinch_tdc_hit_amp;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_hit_clustindex;   //!
   TBranch        *b_bb_grinch_tdc_hit_clustindex;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_hit_col;   //!
   TBranch        *b_bb_grinch_tdc_hit_col;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_hit_pmtnum;   //!
   TBranch        *b_bb_grinch_tdc_hit_pmtnum;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_hit_row;   //!
   TBranch        *b_bb_grinch_tdc_hit_row;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_hit_time;   //!
   TBranch        *b_bb_grinch_tdc_hit_time;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_hit_trackindex;   //!
   TBranch        *b_bb_grinch_tdc_hit_trackindex;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_hit_xhit;   //!
   TBranch        *b_bb_grinch_tdc_hit_xhit;   //!
   TBranch        *b_Ndata_bb_grinch_tdc_hit_yhit;   //!
   TBranch        *b_bb_grinch_tdc_hit_yhit;   //!
   TBranch        *b_Ndata_bb_hodotdc_Ref_tdc;   //!
   TBranch        *b_bb_hodotdc_Ref_tdc;   //!
   TBranch        *b_Ndata_bb_hodotdc_Ref_tdc_mult;   //!
   TBranch        *b_bb_hodotdc_Ref_tdc_mult;   //!
   TBranch        *b_Ndata_bb_hodotdc_Ref_tdc_te;   //!
   TBranch        *b_bb_hodotdc_Ref_tdc_te;   //!
   TBranch        *b_Ndata_bb_hodotdc_Ref_tdc_tot;   //!
   TBranch        *b_bb_hodotdc_Ref_tdc_tot;   //!
   TBranch        *b_Ndata_bb_hodotdc_Ref_tdcelemID;   //!
   TBranch        *b_bb_hodotdc_Ref_tdcelemID;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_bar_tdc_id;   //!
   TBranch        *b_bb_hodotdc_clus_bar_tdc_id;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_bar_tdc_itrack;   //!
   TBranch        *b_bb_hodotdc_clus_bar_tdc_itrack;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_bar_tdc_meantime;   //!
   TBranch        *b_bb_hodotdc_clus_bar_tdc_meantime;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_bar_tdc_meantot;   //!
   TBranch        *b_bb_hodotdc_clus_bar_tdc_meantot;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_bar_tdc_timediff;   //!
   TBranch        *b_bb_hodotdc_clus_bar_tdc_timediff;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_bar_tdc_timehitpos;   //!
   TBranch        *b_bb_hodotdc_clus_bar_tdc_timehitpos;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_bar_tdc_tleft;   //!
   TBranch        *b_bb_hodotdc_clus_bar_tdc_tleft;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_bar_tdc_totleft;   //!
   TBranch        *b_bb_hodotdc_clus_bar_tdc_totleft;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_bar_tdc_totright;   //!
   TBranch        *b_bb_hodotdc_clus_bar_tdc_totright;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_bar_tdc_tright;   //!
   TBranch        *b_bb_hodotdc_clus_bar_tdc_tright;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_bar_tdc_vpos;   //!
   TBranch        *b_bb_hodotdc_clus_bar_tdc_vpos;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_id;   //!
   TBranch        *b_bb_hodotdc_clus_id;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_size;   //!
   TBranch        *b_bb_hodotdc_clus_size;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_tdiff;   //!
   TBranch        *b_bb_hodotdc_clus_tdiff;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_tleft;   //!
   TBranch        *b_bb_hodotdc_clus_tleft;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_tmean;   //!
   TBranch        *b_bb_hodotdc_clus_tmean;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_totleft;   //!
   TBranch        *b_bb_hodotdc_clus_totleft;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_totmean;   //!
   TBranch        *b_bb_hodotdc_clus_totmean;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_totright;   //!
   TBranch        *b_bb_hodotdc_clus_totright;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_trackindex;   //!
   TBranch        *b_bb_hodotdc_clus_trackindex;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_tright;   //!
   TBranch        *b_bb_hodotdc_clus_tright;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_xmean;   //!
   TBranch        *b_bb_hodotdc_clus_xmean;   //!
   TBranch        *b_Ndata_bb_hodotdc_clus_ymean;   //!
   TBranch        *b_bb_hodotdc_clus_ymean;   //!
   TBranch        *b_Ndata_bb_ps_clus_again;   //!
   TBranch        *b_bb_ps_clus_again;   //!
   TBranch        *b_Ndata_bb_ps_clus_atime;   //!
   TBranch        *b_bb_ps_clus_atime;   //!
   TBranch        *b_Ndata_bb_ps_clus_col;   //!
   TBranch        *b_bb_ps_clus_col;   //!
   TBranch        *b_Ndata_bb_ps_clus_e;   //!
   TBranch        *b_bb_ps_clus_e;   //!
   TBranch        *b_Ndata_bb_ps_clus_eblk;   //!
   TBranch        *b_bb_ps_clus_eblk;   //!
   TBranch        *b_Ndata_bb_ps_clus_id;   //!
   TBranch        *b_bb_ps_clus_id;   //!
   TBranch        *b_Ndata_bb_ps_clus_nblk;   //!
   TBranch        *b_bb_ps_clus_nblk;   //!
   TBranch        *b_Ndata_bb_ps_clus_row;   //!
   TBranch        *b_bb_ps_clus_row;   //!
   TBranch        *b_Ndata_bb_ps_clus_tdctime;   //!
   TBranch        *b_bb_ps_clus_tdctime;   //!
   TBranch        *b_Ndata_bb_ps_clus_x;   //!
   TBranch        *b_bb_ps_clus_x;   //!
   TBranch        *b_Ndata_bb_ps_clus_y;   //!
   TBranch        *b_bb_ps_clus_y;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_again;   //!
   TBranch        *b_bb_ps_clus_blk_again;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_atime;   //!
   TBranch        *b_bb_ps_clus_blk_atime;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_col;   //!
   TBranch        *b_bb_ps_clus_blk_col;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_e;   //!
   TBranch        *b_bb_ps_clus_blk_e;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_id;   //!
   TBranch        *b_bb_ps_clus_blk_id;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_row;   //!
   TBranch        *b_bb_ps_clus_blk_row;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_tdctime;   //!
   TBranch        *b_bb_ps_clus_blk_tdctime;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_x;   //!
   TBranch        *b_bb_ps_clus_blk_x;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_y;   //!
   TBranch        *b_bb_ps_clus_blk_y;   //!
   TBranch        *b_Ndata_bb_sh_clus_again;   //!
   TBranch        *b_bb_sh_clus_again;   //!
   TBranch        *b_Ndata_bb_sh_clus_atime;   //!
   TBranch        *b_bb_sh_clus_atime;   //!
   TBranch        *b_Ndata_bb_sh_clus_col;   //!
   TBranch        *b_bb_sh_clus_col;   //!
   TBranch        *b_Ndata_bb_sh_clus_e;   //!
   TBranch        *b_bb_sh_clus_e;   //!
   TBranch        *b_Ndata_bb_sh_clus_eblk;   //!
   TBranch        *b_bb_sh_clus_eblk;   //!
   TBranch        *b_Ndata_bb_sh_clus_id;   //!
   TBranch        *b_bb_sh_clus_id;   //!
   TBranch        *b_Ndata_bb_sh_clus_nblk;   //!
   TBranch        *b_bb_sh_clus_nblk;   //!
   TBranch        *b_Ndata_bb_sh_clus_row;   //!
   TBranch        *b_bb_sh_clus_row;   //!
   TBranch        *b_Ndata_bb_sh_clus_tdctime;   //!
   TBranch        *b_bb_sh_clus_tdctime;   //!
   TBranch        *b_Ndata_bb_sh_clus_x;   //!
   TBranch        *b_bb_sh_clus_x;   //!
   TBranch        *b_Ndata_bb_sh_clus_y;   //!
   TBranch        *b_bb_sh_clus_y;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_again;   //!
   TBranch        *b_bb_sh_clus_blk_again;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_atime;   //!
   TBranch        *b_bb_sh_clus_blk_atime;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_col;   //!
   TBranch        *b_bb_sh_clus_blk_col;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_e;   //!
   TBranch        *b_bb_sh_clus_blk_e;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_id;   //!
   TBranch        *b_bb_sh_clus_blk_id;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_row;   //!
   TBranch        *b_bb_sh_clus_blk_row;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_tdctime;   //!
   TBranch        *b_bb_sh_clus_blk_tdctime;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_x;   //!
   TBranch        *b_bb_sh_clus_blk_x;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_y;   //!
   TBranch        *b_bb_sh_clus_blk_y;   //!
   TBranch        *b_Ndata_bb_tdctrig_tdc;   //!
   TBranch        *b_bb_tdctrig_tdc;   //!
   TBranch        *b_Ndata_bb_tdctrig_tdcelemID;   //!
   TBranch        *b_bb_tdctrig_tdcelemID;   //!
   TBranch        *b_Ndata_bb_tr_beta;   //!
   TBranch        *b_bb_tr_beta;   //!
   TBranch        *b_Ndata_bb_tr_chi2;   //!
   TBranch        *b_bb_tr_chi2;   //!
   TBranch        *b_Ndata_bb_tr_d_ph;   //!
   TBranch        *b_bb_tr_d_ph;   //!
   TBranch        *b_Ndata_bb_tr_d_th;   //!
   TBranch        *b_bb_tr_d_th;   //!
   TBranch        *b_Ndata_bb_tr_d_x;   //!
   TBranch        *b_bb_tr_d_x;   //!
   TBranch        *b_Ndata_bb_tr_d_y;   //!
   TBranch        *b_bb_tr_d_y;   //!
   TBranch        *b_Ndata_bb_tr_dbeta;   //!
   TBranch        *b_bb_tr_dbeta;   //!
   TBranch        *b_Ndata_bb_tr_dtime;   //!
   TBranch        *b_bb_tr_dtime;   //!
   TBranch        *b_Ndata_bb_tr_flag;   //!
   TBranch        *b_bb_tr_flag;   //!
   TBranch        *b_Ndata_bb_tr_ndof;   //!
   TBranch        *b_bb_tr_ndof;   //!
   TBranch        *b_Ndata_bb_tr_p;   //!
   TBranch        *b_bb_tr_p;   //!
   TBranch        *b_Ndata_bb_tr_pathl;   //!
   TBranch        *b_bb_tr_pathl;   //!
   TBranch        *b_Ndata_bb_tr_ph;   //!
   TBranch        *b_bb_tr_ph;   //!
   TBranch        *b_Ndata_bb_tr_px;   //!
   TBranch        *b_bb_tr_px;   //!
   TBranch        *b_Ndata_bb_tr_py;   //!
   TBranch        *b_bb_tr_py;   //!
   TBranch        *b_Ndata_bb_tr_pz;   //!
   TBranch        *b_bb_tr_pz;   //!
   TBranch        *b_Ndata_bb_tr_r_ph;   //!
   TBranch        *b_bb_tr_r_ph;   //!
   TBranch        *b_Ndata_bb_tr_r_th;   //!
   TBranch        *b_bb_tr_r_th;   //!
   TBranch        *b_Ndata_bb_tr_r_x;   //!
   TBranch        *b_bb_tr_r_x;   //!
   TBranch        *b_Ndata_bb_tr_r_y;   //!
   TBranch        *b_bb_tr_r_y;   //!
   TBranch        *b_Ndata_bb_tr_tg_dp;   //!
   TBranch        *b_bb_tr_tg_dp;   //!
   TBranch        *b_Ndata_bb_tr_tg_ph;   //!
   TBranch        *b_bb_tr_tg_ph;   //!
   TBranch        *b_Ndata_bb_tr_tg_th;   //!
   TBranch        *b_bb_tr_tg_th;   //!
   TBranch        *b_Ndata_bb_tr_tg_x;   //!
   TBranch        *b_bb_tr_tg_x;   //!
   TBranch        *b_Ndata_bb_tr_tg_y;   //!
   TBranch        *b_bb_tr_tg_y;   //!
   TBranch        *b_Ndata_bb_tr_th;   //!
   TBranch        *b_bb_tr_th;   //!
   TBranch        *b_Ndata_bb_tr_time;   //!
   TBranch        *b_bb_tr_time;   //!
   TBranch        *b_Ndata_bb_tr_vx;   //!
   TBranch        *b_bb_tr_vx;   //!
   TBranch        *b_Ndata_bb_tr_vy;   //!
   TBranch        *b_bb_tr_vy;   //!
   TBranch        *b_Ndata_bb_tr_vz;   //!
   TBranch        *b_bb_tr_vz;   //!
   TBranch        *b_Ndata_bb_tr_x;   //!
   TBranch        *b_bb_tr_x;   //!
   TBranch        *b_Ndata_bb_tr_y;   //!
   TBranch        *b_bb_tr_y;   //!
   TBranch        *b_Ndata_bb_x_bcp;   //!
   TBranch        *b_bb_x_bcp;   //!
   TBranch        *b_Ndata_bb_x_fcp;   //!
   TBranch        *b_bb_x_fcp;   //!
   TBranch        *b_Ndata_bb_y_bcp;   //!
   TBranch        *b_bb_y_bcp;   //!
   TBranch        *b_Ndata_bb_y_fcp;   //!
   TBranch        *b_bb_y_fcp;   //!
   TBranch        *b_Ndata_bb_z_bcp;   //!
   TBranch        *b_bb_z_bcp;   //!
   TBranch        *b_Ndata_bb_z_fcp;   //!
   TBranch        *b_bb_z_fcp;   //!
   TBranch        *b_Ndata_sbs_hcal_Ref_tdc;   //!
   TBranch        *b_sbs_hcal_Ref_tdc;   //!
   TBranch        *b_Ndata_sbs_hcal_Ref_tdc_mult;   //!
   TBranch        *b_sbs_hcal_Ref_tdc_mult;   //!
   TBranch        *b_Ndata_sbs_hcal_Ref_tdcelemID;   //!
   TBranch        *b_sbs_hcal_Ref_tdcelemID;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_again;   //!
   TBranch        *b_sbs_hcal_clus_again;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_atime;   //!
   TBranch        *b_sbs_hcal_clus_atime;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_col;   //!
   TBranch        *b_sbs_hcal_clus_col;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_e;   //!
   TBranch        *b_sbs_hcal_clus_e;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_eblk;   //!
   TBranch        *b_sbs_hcal_clus_eblk;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_id;   //!
   TBranch        *b_sbs_hcal_clus_id;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_nblk;   //!
   TBranch        *b_sbs_hcal_clus_nblk;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_row;   //!
   TBranch        *b_sbs_hcal_clus_row;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_tdctime;   //!
   TBranch        *b_sbs_hcal_clus_tdctime;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_x;   //!
   TBranch        *b_sbs_hcal_clus_x;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_y;   //!
   TBranch        *b_sbs_hcal_clus_y;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_again;   //!
   TBranch        *b_sbs_hcal_clus_blk_again;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_atime;   //!
   TBranch        *b_sbs_hcal_clus_blk_atime;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_col;   //!
   TBranch        *b_sbs_hcal_clus_blk_col;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_e;   //!
   TBranch        *b_sbs_hcal_clus_blk_e;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_id;   //!
   TBranch        *b_sbs_hcal_clus_blk_id;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_row;   //!
   TBranch        *b_sbs_hcal_clus_blk_row;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_tdctime;   //!
   TBranch        *b_sbs_hcal_clus_blk_tdctime;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_x;   //!
   TBranch        *b_sbs_hcal_clus_blk_x;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_y;   //!
   TBranch        *b_sbs_hcal_clus_blk_y;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_atime;   //!
   TBranch        *b_sbs_hcal_goodblock_atime;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_cid;   //!
   TBranch        *b_sbs_hcal_goodblock_cid;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_col;   //!
   TBranch        *b_sbs_hcal_goodblock_col;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_e;   //!
   TBranch        *b_sbs_hcal_goodblock_e;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_id;   //!
   TBranch        *b_sbs_hcal_goodblock_id;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_row;   //!
   TBranch        *b_sbs_hcal_goodblock_row;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_tdctime;   //!
   TBranch        *b_sbs_hcal_goodblock_tdctime;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_x;   //!
   TBranch        *b_sbs_hcal_goodblock_x;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_y;   //!
   TBranch        *b_sbs_hcal_goodblock_y;   //!
   TBranch        *b_Ndata_sbs_tr_beta;   //!
   TBranch        *b_sbs_tr_beta;   //!
   TBranch        *b_Ndata_sbs_tr_chi2;   //!
   TBranch        *b_sbs_tr_chi2;   //!
   TBranch        *b_Ndata_sbs_tr_d_ph;   //!
   TBranch        *b_sbs_tr_d_ph;   //!
   TBranch        *b_Ndata_sbs_tr_d_th;   //!
   TBranch        *b_sbs_tr_d_th;   //!
   TBranch        *b_Ndata_sbs_tr_d_x;   //!
   TBranch        *b_sbs_tr_d_x;   //!
   TBranch        *b_Ndata_sbs_tr_d_y;   //!
   TBranch        *b_sbs_tr_d_y;   //!
   TBranch        *b_Ndata_sbs_tr_dbeta;   //!
   TBranch        *b_sbs_tr_dbeta;   //!
   TBranch        *b_Ndata_sbs_tr_dtime;   //!
   TBranch        *b_sbs_tr_dtime;   //!
   TBranch        *b_Ndata_sbs_tr_flag;   //!
   TBranch        *b_sbs_tr_flag;   //!
   TBranch        *b_Ndata_sbs_tr_ndof;   //!
   TBranch        *b_sbs_tr_ndof;   //!
   TBranch        *b_Ndata_sbs_tr_p;   //!
   TBranch        *b_sbs_tr_p;   //!
   TBranch        *b_Ndata_sbs_tr_pathl;   //!
   TBranch        *b_sbs_tr_pathl;   //!
   TBranch        *b_Ndata_sbs_tr_ph;   //!
   TBranch        *b_sbs_tr_ph;   //!
   TBranch        *b_Ndata_sbs_tr_px;   //!
   TBranch        *b_sbs_tr_px;   //!
   TBranch        *b_Ndata_sbs_tr_py;   //!
   TBranch        *b_sbs_tr_py;   //!
   TBranch        *b_Ndata_sbs_tr_pz;   //!
   TBranch        *b_sbs_tr_pz;   //!
   TBranch        *b_Ndata_sbs_tr_r_ph;   //!
   TBranch        *b_sbs_tr_r_ph;   //!
   TBranch        *b_Ndata_sbs_tr_r_th;   //!
   TBranch        *b_sbs_tr_r_th;   //!
   TBranch        *b_Ndata_sbs_tr_r_x;   //!
   TBranch        *b_sbs_tr_r_x;   //!
   TBranch        *b_Ndata_sbs_tr_r_y;   //!
   TBranch        *b_sbs_tr_r_y;   //!
   TBranch        *b_Ndata_sbs_tr_tg_dp;   //!
   TBranch        *b_sbs_tr_tg_dp;   //!
   TBranch        *b_Ndata_sbs_tr_tg_ph;   //!
   TBranch        *b_sbs_tr_tg_ph;   //!
   TBranch        *b_Ndata_sbs_tr_tg_th;   //!
   TBranch        *b_sbs_tr_tg_th;   //!
   TBranch        *b_Ndata_sbs_tr_tg_y;   //!
   TBranch        *b_sbs_tr_tg_y;   //!
   TBranch        *b_Ndata_sbs_tr_th;   //!
   TBranch        *b_sbs_tr_th;   //!
   TBranch        *b_Ndata_sbs_tr_time;   //!
   TBranch        *b_sbs_tr_time;   //!
   TBranch        *b_Ndata_sbs_tr_vx;   //!
   TBranch        *b_sbs_tr_vx;   //!
   TBranch        *b_Ndata_sbs_tr_vy;   //!
   TBranch        *b_sbs_tr_vy;   //!
   TBranch        *b_Ndata_sbs_tr_vz;   //!
   TBranch        *b_sbs_tr_vz;   //!
   TBranch        *b_Ndata_sbs_tr_x;   //!
   TBranch        *b_sbs_tr_x;   //!
   TBranch        *b_Ndata_sbs_tr_y;   //!
   TBranch        *b_sbs_tr_y;   //!
   TBranch        *b_Ndata_sbs_trig_a_amp_p;   //!
   TBranch        *b_sbs_trig_a_amp_p;   //!
   TBranch        *b_Ndata_sbs_trig_a_p;   //!
   TBranch        *b_sbs_trig_a_p;   //!
   TBranch        *b_Ndata_sbs_trig_a_time;   //!
   TBranch        *b_sbs_trig_a_time;   //!
   TBranch        *b_Ndata_sbs_trig_adcelemID;   //!
   TBranch        *b_sbs_trig_adcelemID;   //!
   TBranch        *b_Ndata_sbs_x_bcp;   //!
   TBranch        *b_sbs_x_bcp;   //!
   TBranch        *b_Ndata_sbs_x_fcp;   //!
   TBranch        *b_sbs_x_fcp;   //!
   TBranch        *b_Ndata_sbs_y_bcp;   //!
   TBranch        *b_sbs_y_bcp;   //!
   TBranch        *b_Ndata_sbs_y_fcp;   //!
   TBranch        *b_sbs_y_fcp;   //!
   TBranch        *b_Ndata_sbs_z_bcp;   //!
   TBranch        *b_sbs_z_bcp;   //!
   TBranch        *b_Ndata_sbs_z_fcp;   //!
   TBranch        *b_sbs_z_fcp;   //!
   TBranch        *b_BB_gold_beta;   //!
   TBranch        *b_BB_gold_dp;   //!
   TBranch        *b_BB_gold_index;   //!
   TBranch        *b_BB_gold_ok;   //!
   TBranch        *b_BB_gold_p;   //!
   TBranch        *b_BB_gold_ph;   //!
   TBranch        *b_BB_gold_px;   //!
   TBranch        *b_BB_gold_py;   //!
   TBranch        *b_BB_gold_pz;   //!
   TBranch        *b_BB_gold_th;   //!
   TBranch        *b_BB_gold_x;   //!
   TBranch        *b_BB_gold_y;   //!
   TBranch        *b_Lrb_BPMA_rawcur_1;   //!
   TBranch        *b_Lrb_BPMA_rawcur_2;   //!
   TBranch        *b_Lrb_BPMA_rawcur_3;   //!
   TBranch        *b_Lrb_BPMA_rawcur_4;   //!
   TBranch        *b_Lrb_BPMA_rotpos1;   //!
   TBranch        *b_Lrb_BPMA_rotpos2;   //!
   TBranch        *b_Lrb_BPMA_x;   //!
   TBranch        *b_Lrb_BPMA_y;   //!
   TBranch        *b_Lrb_BPMA_z;   //!
   TBranch        *b_Lrb_BPMB_rawcur_1;   //!
   TBranch        *b_Lrb_BPMB_rawcur_2;   //!
   TBranch        *b_Lrb_BPMB_rawcur_3;   //!
   TBranch        *b_Lrb_BPMB_rawcur_4;   //!
   TBranch        *b_Lrb_BPMB_rotpos1;   //!
   TBranch        *b_Lrb_BPMB_rotpos2;   //!
   TBranch        *b_Lrb_BPMB_x;   //!
   TBranch        *b_Lrb_BPMB_y;   //!
   TBranch        *b_Lrb_BPMB_z;   //!
   TBranch        *b_Lrb_Raster_bpma_x;   //!
   TBranch        *b_Lrb_Raster_bpma_y;   //!
   TBranch        *b_Lrb_Raster_bpma_z;   //!
   TBranch        *b_Lrb_Raster_bpmb_x;   //!
   TBranch        *b_Lrb_Raster_bpmb_y;   //!
   TBranch        *b_Lrb_Raster_bpmb_z;   //!
   TBranch        *b_Lrb_Raster_rawcur_x;   //!
   TBranch        *b_Lrb_Raster_rawcur_y;   //!
   TBranch        *b_Lrb_Raster_rawslope_x;   //!
   TBranch        *b_Lrb_Raster_rawslope_y;   //!
   TBranch        *b_Lrb_Raster_target_dir_x;   //!
   TBranch        *b_Lrb_Raster_target_dir_y;   //!
   TBranch        *b_Lrb_Raster_target_dir_z;   //!
   TBranch        *b_Lrb_Raster_target_x;   //!
   TBranch        *b_Lrb_Raster_target_y;   //!
   TBranch        *b_Lrb_Raster_target_z;   //!
   TBranch        *b_Lrb_Raster2_bpma_x;   //!
   TBranch        *b_Lrb_Raster2_bpma_y;   //!
   TBranch        *b_Lrb_Raster2_bpma_z;   //!
   TBranch        *b_Lrb_Raster2_bpmb_x;   //!
   TBranch        *b_Lrb_Raster2_bpmb_y;   //!
   TBranch        *b_Lrb_Raster2_bpmb_z;   //!
   TBranch        *b_Lrb_Raster2_rawcur_x;   //!
   TBranch        *b_Lrb_Raster2_rawcur_y;   //!
   TBranch        *b_Lrb_Raster2_rawslope_x;   //!
   TBranch        *b_Lrb_Raster2_rawslope_y;   //!
   TBranch        *b_Lrb_Raster2_target_dir_x;   //!
   TBranch        *b_Lrb_Raster2_target_dir_y;   //!
   TBranch        *b_Lrb_Raster2_target_dir_z;   //!
   TBranch        *b_Lrb_Raster2_target_x;   //!
   TBranch        *b_Lrb_Raster2_target_y;   //!
   TBranch        *b_Lrb_Raster2_target_z;   //!
   TBranch        *b_SBSrb_BPMA_rawcur_1;   //!
   TBranch        *b_SBSrb_BPMA_rawcur_2;   //!
   TBranch        *b_SBSrb_BPMA_rawcur_3;   //!
   TBranch        *b_SBSrb_BPMA_rawcur_4;   //!
   TBranch        *b_SBSrb_BPMA_rotpos1;   //!
   TBranch        *b_SBSrb_BPMA_rotpos2;   //!
   TBranch        *b_SBSrb_BPMA_x;   //!
   TBranch        *b_SBSrb_BPMA_y;   //!
   TBranch        *b_SBSrb_BPMA_z;   //!
   TBranch        *b_SBSrb_BPMB_rawcur_1;   //!
   TBranch        *b_SBSrb_BPMB_rawcur_2;   //!
   TBranch        *b_SBSrb_BPMB_rawcur_3;   //!
   TBranch        *b_SBSrb_BPMB_rawcur_4;   //!
   TBranch        *b_SBSrb_BPMB_rotpos1;   //!
   TBranch        *b_SBSrb_BPMB_rotpos2;   //!
   TBranch        *b_SBSrb_BPMB_x;   //!
   TBranch        *b_SBSrb_BPMB_y;   //!
   TBranch        *b_SBSrb_BPMB_z;   //!
   TBranch        *b_SBSrb_Raster_bpma_x;   //!
   TBranch        *b_SBSrb_Raster_bpma_y;   //!
   TBranch        *b_SBSrb_Raster_bpma_z;   //!
   TBranch        *b_SBSrb_Raster_bpmb_x;   //!
   TBranch        *b_SBSrb_Raster_bpmb_y;   //!
   TBranch        *b_SBSrb_Raster_bpmb_z;   //!
   TBranch        *b_SBSrb_Raster_rawcur_x;   //!
   TBranch        *b_SBSrb_Raster_rawcur_y;   //!
   TBranch        *b_SBSrb_Raster_rawslope_x;   //!
   TBranch        *b_SBSrb_Raster_rawslope_y;   //!
   TBranch        *b_SBSrb_Raster_target_dir_x;   //!
   TBranch        *b_SBSrb_Raster_target_dir_y;   //!
   TBranch        *b_SBSrb_Raster_target_dir_z;   //!
   TBranch        *b_SBSrb_Raster_target_x;   //!
   TBranch        *b_SBSrb_Raster_target_y;   //!
   TBranch        *b_SBSrb_Raster_target_z;   //!
   TBranch        *b_SBSrb_Raster2_bpma_x;   //!
   TBranch        *b_SBSrb_Raster2_bpma_y;   //!
   TBranch        *b_SBSrb_Raster2_bpma_z;   //!
   TBranch        *b_SBSrb_Raster2_bpmb_x;   //!
   TBranch        *b_SBSrb_Raster2_bpmb_y;   //!
   TBranch        *b_SBSrb_Raster2_bpmb_z;   //!
   TBranch        *b_SBSrb_Raster2_rawcur_x;   //!
   TBranch        *b_SBSrb_Raster2_rawcur_y;   //!
   TBranch        *b_SBSrb_Raster2_rawslope_x;   //!
   TBranch        *b_SBSrb_Raster2_rawslope_y;   //!
   TBranch        *b_SBSrb_Raster2_target_dir_x;   //!
   TBranch        *b_SBSrb_Raster2_target_dir_y;   //!
   TBranch        *b_SBSrb_Raster2_target_dir_z;   //!
   TBranch        *b_SBSrb_Raster2_target_x;   //!
   TBranch        *b_SBSrb_Raster2_target_y;   //!
   TBranch        *b_SBSrb_Raster2_target_z;   //!
   TBranch        *b_bb_gem_hit_ngoodhits;   //!
   TBranch        *b_bb_gem_m0_clust_nclustu;   //!
   TBranch        *b_bb_gem_m0_clust_nclustu_neg;   //!
   TBranch        *b_bb_gem_m0_clust_nclustu_tot;   //!
   TBranch        *b_bb_gem_m0_clust_nclustv;   //!
   TBranch        *b_bb_gem_m0_clust_nclustv_neg;   //!
   TBranch        *b_bb_gem_m0_clust_nclustv_tot;   //!
   TBranch        *b_bb_gem_m0_strip_nstrips_keep;   //!
   TBranch        *b_bb_gem_m0_strip_nstrips_keepU;   //!
   TBranch        *b_bb_gem_m0_strip_nstrips_keepV;   //!
   TBranch        *b_bb_gem_m0_strip_nstrips_keep_lmax;   //!
   TBranch        *b_bb_gem_m0_strip_nstrips_keep_lmaxU;   //!
   TBranch        *b_bb_gem_m0_strip_nstrips_keep_lmaxV;   //!
   TBranch        *b_bb_gem_m0_strip_nstripsfired;   //!
   TBranch        *b_bb_gem_m0_strip_nstripsfired_neg;   //!
   TBranch        *b_bb_gem_m1_clust_nclustu;   //!
   TBranch        *b_bb_gem_m1_clust_nclustu_neg;   //!
   TBranch        *b_bb_gem_m1_clust_nclustu_tot;   //!
   TBranch        *b_bb_gem_m1_clust_nclustv;   //!
   TBranch        *b_bb_gem_m1_clust_nclustv_neg;   //!
   TBranch        *b_bb_gem_m1_clust_nclustv_tot;   //!
   TBranch        *b_bb_gem_m1_strip_nstrips_keep;   //!
   TBranch        *b_bb_gem_m1_strip_nstrips_keepU;   //!
   TBranch        *b_bb_gem_m1_strip_nstrips_keepV;   //!
   TBranch        *b_bb_gem_m1_strip_nstrips_keep_lmax;   //!
   TBranch        *b_bb_gem_m1_strip_nstrips_keep_lmaxU;   //!
   TBranch        *b_bb_gem_m1_strip_nstrips_keep_lmaxV;   //!
   TBranch        *b_bb_gem_m1_strip_nstripsfired;   //!
   TBranch        *b_bb_gem_m1_strip_nstripsfired_neg;   //!
   TBranch        *b_bb_gem_m2_clust_nclustu;   //!
   TBranch        *b_bb_gem_m2_clust_nclustu_neg;   //!
   TBranch        *b_bb_gem_m2_clust_nclustu_tot;   //!
   TBranch        *b_bb_gem_m2_clust_nclustv;   //!
   TBranch        *b_bb_gem_m2_clust_nclustv_neg;   //!
   TBranch        *b_bb_gem_m2_clust_nclustv_tot;   //!
   TBranch        *b_bb_gem_m2_strip_nstrips_keep;   //!
   TBranch        *b_bb_gem_m2_strip_nstrips_keepU;   //!
   TBranch        *b_bb_gem_m2_strip_nstrips_keepV;   //!
   TBranch        *b_bb_gem_m2_strip_nstrips_keep_lmax;   //!
   TBranch        *b_bb_gem_m2_strip_nstrips_keep_lmaxU;   //!
   TBranch        *b_bb_gem_m2_strip_nstrips_keep_lmaxV;   //!
   TBranch        *b_bb_gem_m2_strip_nstripsfired;   //!
   TBranch        *b_bb_gem_m2_strip_nstripsfired_neg;   //!
   TBranch        *b_bb_gem_m3_clust_nclustu;   //!
   TBranch        *b_bb_gem_m3_clust_nclustu_neg;   //!
   TBranch        *b_bb_gem_m3_clust_nclustu_tot;   //!
   TBranch        *b_bb_gem_m3_clust_nclustv;   //!
   TBranch        *b_bb_gem_m3_clust_nclustv_neg;   //!
   TBranch        *b_bb_gem_m3_clust_nclustv_tot;   //!
   TBranch        *b_bb_gem_m3_strip_nstrips_keep;   //!
   TBranch        *b_bb_gem_m3_strip_nstrips_keepU;   //!
   TBranch        *b_bb_gem_m3_strip_nstrips_keepV;   //!
   TBranch        *b_bb_gem_m3_strip_nstrips_keep_lmax;   //!
   TBranch        *b_bb_gem_m3_strip_nstrips_keep_lmaxU;   //!
   TBranch        *b_bb_gem_m3_strip_nstrips_keep_lmaxV;   //!
   TBranch        *b_bb_gem_m3_strip_nstripsfired;   //!
   TBranch        *b_bb_gem_m3_strip_nstripsfired_neg;   //!
   TBranch        *b_bb_gem_m4_clust_nclustu;   //!
   TBranch        *b_bb_gem_m4_clust_nclustu_neg;   //!
   TBranch        *b_bb_gem_m4_clust_nclustu_tot;   //!
   TBranch        *b_bb_gem_m4_clust_nclustv;   //!
   TBranch        *b_bb_gem_m4_clust_nclustv_neg;   //!
   TBranch        *b_bb_gem_m4_clust_nclustv_tot;   //!
   TBranch        *b_bb_gem_m4_strip_nstrips_keep;   //!
   TBranch        *b_bb_gem_m4_strip_nstrips_keepU;   //!
   TBranch        *b_bb_gem_m4_strip_nstrips_keepV;   //!
   TBranch        *b_bb_gem_m4_strip_nstrips_keep_lmax;   //!
   TBranch        *b_bb_gem_m4_strip_nstrips_keep_lmaxU;   //!
   TBranch        *b_bb_gem_m4_strip_nstrips_keep_lmaxV;   //!
   TBranch        *b_bb_gem_m4_strip_nstripsfired;   //!
   TBranch        *b_bb_gem_m4_strip_nstripsfired_neg;   //!
   TBranch        *b_bb_gem_m5_clust_nclustu;   //!
   TBranch        *b_bb_gem_m5_clust_nclustu_neg;   //!
   TBranch        *b_bb_gem_m5_clust_nclustu_tot;   //!
   TBranch        *b_bb_gem_m5_clust_nclustv;   //!
   TBranch        *b_bb_gem_m5_clust_nclustv_neg;   //!
   TBranch        *b_bb_gem_m5_clust_nclustv_tot;   //!
   TBranch        *b_bb_gem_m5_strip_nstrips_keep;   //!
   TBranch        *b_bb_gem_m5_strip_nstrips_keepU;   //!
   TBranch        *b_bb_gem_m5_strip_nstrips_keepV;   //!
   TBranch        *b_bb_gem_m5_strip_nstrips_keep_lmax;   //!
   TBranch        *b_bb_gem_m5_strip_nstrips_keep_lmaxU;   //!
   TBranch        *b_bb_gem_m5_strip_nstrips_keep_lmaxV;   //!
   TBranch        *b_bb_gem_m5_strip_nstripsfired;   //!
   TBranch        *b_bb_gem_m5_strip_nstripsfired_neg;   //!
   TBranch        *b_bb_gem_m6_clust_nclustu;   //!
   TBranch        *b_bb_gem_m6_clust_nclustu_neg;   //!
   TBranch        *b_bb_gem_m6_clust_nclustu_tot;   //!
   TBranch        *b_bb_gem_m6_clust_nclustv;   //!
   TBranch        *b_bb_gem_m6_clust_nclustv_neg;   //!
   TBranch        *b_bb_gem_m6_clust_nclustv_tot;   //!
   TBranch        *b_bb_gem_m6_strip_nstrips_keep;   //!
   TBranch        *b_bb_gem_m6_strip_nstrips_keepU;   //!
   TBranch        *b_bb_gem_m6_strip_nstrips_keepV;   //!
   TBranch        *b_bb_gem_m6_strip_nstrips_keep_lmax;   //!
   TBranch        *b_bb_gem_m6_strip_nstrips_keep_lmaxU;   //!
   TBranch        *b_bb_gem_m6_strip_nstrips_keep_lmaxV;   //!
   TBranch        *b_bb_gem_m6_strip_nstripsfired;   //!
   TBranch        *b_bb_gem_m6_strip_nstripsfired_neg;   //!
   TBranch        *b_bb_gem_m7_clust_nclustu;   //!
   TBranch        *b_bb_gem_m7_clust_nclustu_neg;   //!
   TBranch        *b_bb_gem_m7_clust_nclustu_tot;   //!
   TBranch        *b_bb_gem_m7_clust_nclustv;   //!
   TBranch        *b_bb_gem_m7_clust_nclustv_neg;   //!
   TBranch        *b_bb_gem_m7_clust_nclustv_tot;   //!
   TBranch        *b_bb_gem_m7_strip_nstrips_keep;   //!
   TBranch        *b_bb_gem_m7_strip_nstrips_keepU;   //!
   TBranch        *b_bb_gem_m7_strip_nstrips_keepV;   //!
   TBranch        *b_bb_gem_m7_strip_nstrips_keep_lmax;   //!
   TBranch        *b_bb_gem_m7_strip_nstrips_keep_lmaxU;   //!
   TBranch        *b_bb_gem_m7_strip_nstrips_keep_lmaxV;   //!
   TBranch        *b_bb_gem_m7_strip_nstripsfired;   //!
   TBranch        *b_bb_gem_m7_strip_nstripsfired_neg;   //!
   TBranch        *b_bb_gem_nlayershit;   //!
   TBranch        *b_bb_gem_nlayershitu;   //!
   TBranch        *b_bb_gem_nlayershituv;   //!
   TBranch        *b_bb_gem_nlayershitv;   //!
   TBranch        *b_bb_gem_track_besttrack;   //!
   TBranch        *b_bb_gem_track_ntrack;   //!
   TBranch        *b_bb_gem_trigtime;   //!
   TBranch        *b_bb_grinch_tdc_bestcluster;   //!
   TBranch        *b_bb_grinch_tdc_clus_adc;   //!
   TBranch        *b_bb_grinch_tdc_clus_dx;   //!
   TBranch        *b_bb_grinch_tdc_clus_dy;   //!
   TBranch        *b_bb_grinch_tdc_clus_mirrorindex;   //!
   TBranch        *b_bb_grinch_tdc_clus_size;   //!
   TBranch        *b_bb_grinch_tdc_clus_t_mean;   //!
   TBranch        *b_bb_grinch_tdc_clus_t_rms;   //!
   TBranch        *b_bb_grinch_tdc_clus_tot_mean;   //!
   TBranch        *b_bb_grinch_tdc_clus_trackindex;   //!
   TBranch        *b_bb_grinch_tdc_clus_x_mean;   //!
   TBranch        *b_bb_grinch_tdc_clus_y_mean;   //!
   TBranch        *b_bb_grinch_tdc_nclus;   //!
   TBranch        *b_bb_grinch_tdc_ngoodhits;   //!
   TBranch        *b_bb_grinch_tdc_ntrackmatch;   //!
   TBranch        *b_bb_hodotdc_nclus;   //!
   TBranch        *b_bb_ps_againblk;   //!
   TBranch        *b_bb_ps_atimeblk;   //!
   TBranch        *b_bb_ps_colblk;   //!
   TBranch        *b_bb_ps_e;   //!
   TBranch        *b_bb_ps_eblk;   //!
   TBranch        *b_bb_ps_idblk;   //!
   TBranch        *b_bb_ps_index;   //!
   TBranch        *b_bb_ps_nblk;   //!
   TBranch        *b_bb_ps_nclus;   //!
   TBranch        *b_bb_ps_ngoodADChits;   //!
   TBranch        *b_bb_ps_rowblk;   //!
   TBranch        *b_bb_ps_x;   //!
   TBranch        *b_bb_ps_y;   //!
   TBranch        *b_bb_sh_againblk;   //!
   TBranch        *b_bb_sh_atimeblk;   //!
   TBranch        *b_bb_sh_colblk;   //!
   TBranch        *b_bb_sh_e;   //!
   TBranch        *b_bb_sh_eblk;   //!
   TBranch        *b_bb_sh_idblk;   //!
   TBranch        *b_bb_sh_index;   //!
   TBranch        *b_bb_sh_nblk;   //!
   TBranch        *b_bb_sh_nclus;   //!
   TBranch        *b_bb_sh_ngoodADChits;   //!
   TBranch        *b_bb_sh_rowblk;   //!
   TBranch        *b_bb_sh_x;   //!
   TBranch        *b_bb_sh_y;   //!
   TBranch        *b_bb_tr_n;   //!
   TBranch        *b_bb_ts_over_threshold;   //!
   TBranch        *b_e_kine_Q2;   //!
   TBranch        *b_e_kine_W2;   //!
   TBranch        *b_e_kine_angle;   //!
   TBranch        *b_e_kine_epsilon;   //!
   TBranch        *b_e_kine_nu;   //!
   TBranch        *b_e_kine_omega;   //!
   TBranch        *b_e_kine_ph_q;   //!
   TBranch        *b_e_kine_q3m;   //!
   TBranch        *b_e_kine_q_x;   //!
   TBranch        *b_e_kine_q_y;   //!
   TBranch        *b_e_kine_q_z;   //!
   TBranch        *b_e_kine_th_q;   //!
   TBranch        *b_e_kine_x_bj;   //!
   TBranch        *b_g_datatype;   //!
   TBranch        *b_g_evlen;   //!
   TBranch        *b_g_evnum;   //!
   TBranch        *b_g_evtime;   //!
   TBranch        *b_g_evtyp;   //!
   TBranch        *b_g_runnum;   //!
   TBranch        *b_g_runtime;   //!
   TBranch        *b_g_runtype;   //!
   TBranch        *b_g_trigbits;   //!
   TBranch        *b_sbs_HCALdir_x;   //!
   TBranch        *b_sbs_HCALdir_y;   //!
   TBranch        *b_sbs_HCALdir_z;   //!
   TBranch        *b_sbs_HCALph_n;   //!
   TBranch        *b_sbs_HCALth_n;   //!
   TBranch        *b_sbs_hcal_againblk;   //!
   TBranch        *b_sbs_hcal_atimeblk;   //!
   TBranch        *b_sbs_hcal_colblk;   //!
   TBranch        *b_sbs_hcal_e;   //!
   TBranch        *b_sbs_hcal_eblk;   //!
   TBranch        *b_sbs_hcal_idblk;   //!
   TBranch        *b_sbs_hcal_index;   //!
   TBranch        *b_sbs_hcal_nblk;   //!
   TBranch        *b_sbs_hcal_nclus;   //!
   TBranch        *b_sbs_hcal_ngoodADChits;   //!
   TBranch        *b_sbs_hcal_rowblk;   //!
   TBranch        *b_sbs_hcal_tdctimeblk;   //!
   TBranch        *b_sbs_hcal_x;   //!
   TBranch        *b_sbs_hcal_y;   //!
   TBranch        *b_sbs_tr_n;   //!
   TBranch        *b_singletrack;   //!
   TBranch        *b_anytrack;   //!
   TBranch        *b_HALLA_p;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtTime;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtNum;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtType;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtLen;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fHelicity;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fTrigBits;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fRun;   //!

   gmn_tree_allstripoutputs(TTree *tree=0);
   virtual ~gmn_tree_allstripoutputs();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef gmn_tree_allstripoutputs_cxx
gmn_tree_allstripoutputs::gmn_tree_allstripoutputs(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/volatile/halla/sbs/puckett/GENRP_REPLAYS/test/rootfiles/e1209019_replayed_13486_stream0_seg0_0_firstevent0_nevent100.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/volatile/halla/sbs/puckett/GENRP_REPLAYS/test/rootfiles/e1209019_replayed_13486_stream0_seg0_0_firstevent0_nevent100.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

gmn_tree_allstripoutputs::~gmn_tree_allstripoutputs()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gmn_tree_allstripoutputs::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gmn_tree_allstripoutputs::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void gmn_tree_allstripoutputs::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ax_new", &ax_new, &b_ax_new);
   fChain->SetBranchAddress("bx_new", &bx_new, &b_bx_new);
   fChain->SetBranchAddress("ay_new", &ay_new, &b_ay_new);
   fChain->SetBranchAddress("by_new", &by_new, &b_by_new);
   fChain->SetBranchAddress("targx", &targx, &b_targx);
   fChain->SetBranchAddress("targy", &targy, &b_targy);
   fChain->SetBranchAddress("Ndata.bb.bbtrig.a_amp_p", &Ndata_bb_bbtrig_a_amp_p, &b_Ndata_bb_bbtrig_a_amp_p);
   fChain->SetBranchAddress("bb.bbtrig.a_amp_p", bb_bbtrig_a_amp_p, &b_bb_bbtrig_a_amp_p);
   fChain->SetBranchAddress("Ndata.bb.bbtrig.a_time", &Ndata_bb_bbtrig_a_time, &b_Ndata_bb_bbtrig_a_time);
   fChain->SetBranchAddress("bb.bbtrig.a_time", bb_bbtrig_a_time, &b_bb_bbtrig_a_time);
   fChain->SetBranchAddress("Ndata.bb.bbtrig.adcelemID", &Ndata_bb_bbtrig_adcelemID, &b_Ndata_bb_bbtrig_adcelemID);
   fChain->SetBranchAddress("bb.bbtrig.adcelemID", bb_bbtrig_adcelemID, &b_bb_bbtrig_adcelemID);
   fChain->SetBranchAddress("Ndata.bb.bbtrig.tdc", &Ndata_bb_bbtrig_tdc, &b_Ndata_bb_bbtrig_tdc);
   fChain->SetBranchAddress("bb.bbtrig.tdc", bb_bbtrig_tdc, &b_bb_bbtrig_tdc);
   fChain->SetBranchAddress("Ndata.bb.bbtrig.tdcelemID", &Ndata_bb_bbtrig_tdcelemID, &b_Ndata_bb_bbtrig_tdcelemID);
   fChain->SetBranchAddress("bb.bbtrig.tdcelemID", bb_bbtrig_tdcelemID, &b_bb_bbtrig_tdcelemID);
   fChain->SetBranchAddress("Ndata.bb.eps_over_etot", &Ndata_bb_eps_over_etot, &b_Ndata_bb_eps_over_etot);
   fChain->SetBranchAddress("bb.eps_over_etot", bb_eps_over_etot, &b_bb_eps_over_etot);
   fChain->SetBranchAddress("Ndata.bb.etot_over_p", &Ndata_bb_etot_over_p, &b_Ndata_bb_etot_over_p);
   fChain->SetBranchAddress("bb.etot_over_p", bb_etot_over_p, &b_bb_etot_over_p);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCU", &Ndata_bb_gem_hit_ADCU, &b_Ndata_bb_gem_hit_ADCU);
   fChain->SetBranchAddress("bb.gem.hit.ADCU", bb_gem_hit_ADCU, &b_bb_gem_hit_ADCU);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCU_deconv", &Ndata_bb_gem_hit_ADCU_deconv, &b_Ndata_bb_gem_hit_ADCU_deconv);
   fChain->SetBranchAddress("bb.gem.hit.ADCU_deconv", bb_gem_hit_ADCU_deconv, &b_bb_gem_hit_ADCU_deconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCV", &Ndata_bb_gem_hit_ADCV, &b_Ndata_bb_gem_hit_ADCV);
   fChain->SetBranchAddress("bb.gem.hit.ADCV", bb_gem_hit_ADCV, &b_bb_gem_hit_ADCV);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCV_deconv", &Ndata_bb_gem_hit_ADCV_deconv, &b_Ndata_bb_gem_hit_ADCV_deconv);
   fChain->SetBranchAddress("bb.gem.hit.ADCV_deconv", bb_gem_hit_ADCV_deconv, &b_bb_gem_hit_ADCV_deconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCasym", &Ndata_bb_gem_hit_ADCasym, &b_Ndata_bb_gem_hit_ADCasym);
   fChain->SetBranchAddress("bb.gem.hit.ADCasym", bb_gem_hit_ADCasym, &b_bb_gem_hit_ADCasym);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCasym_deconv", &Ndata_bb_gem_hit_ADCasym_deconv, &b_Ndata_bb_gem_hit_ADCasym_deconv);
   fChain->SetBranchAddress("bb.gem.hit.ADCasym_deconv", bb_gem_hit_ADCasym_deconv, &b_bb_gem_hit_ADCasym_deconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCavg", &Ndata_bb_gem_hit_ADCavg, &b_Ndata_bb_gem_hit_ADCavg);
   fChain->SetBranchAddress("bb.gem.hit.ADCavg", bb_gem_hit_ADCavg, &b_bb_gem_hit_ADCavg);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCavg_deconv", &Ndata_bb_gem_hit_ADCavg_deconv, &b_Ndata_bb_gem_hit_ADCavg_deconv);
   fChain->SetBranchAddress("bb.gem.hit.ADCavg_deconv", bb_gem_hit_ADCavg_deconv, &b_bb_gem_hit_ADCavg_deconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCfrac0_Umax", &Ndata_bb_gem_hit_ADCfrac0_Umax, &b_Ndata_bb_gem_hit_ADCfrac0_Umax);
   fChain->SetBranchAddress("bb.gem.hit.ADCfrac0_Umax", bb_gem_hit_ADCfrac0_Umax, &b_bb_gem_hit_ADCfrac0_Umax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCfrac0_Vmax", &Ndata_bb_gem_hit_ADCfrac0_Vmax, &b_Ndata_bb_gem_hit_ADCfrac0_Vmax);
   fChain->SetBranchAddress("bb.gem.hit.ADCfrac0_Vmax", bb_gem_hit_ADCfrac0_Vmax, &b_bb_gem_hit_ADCfrac0_Vmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCfrac1_Umax", &Ndata_bb_gem_hit_ADCfrac1_Umax, &b_Ndata_bb_gem_hit_ADCfrac1_Umax);
   fChain->SetBranchAddress("bb.gem.hit.ADCfrac1_Umax", bb_gem_hit_ADCfrac1_Umax, &b_bb_gem_hit_ADCfrac1_Umax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCfrac1_Vmax", &Ndata_bb_gem_hit_ADCfrac1_Vmax, &b_Ndata_bb_gem_hit_ADCfrac1_Vmax);
   fChain->SetBranchAddress("bb.gem.hit.ADCfrac1_Vmax", bb_gem_hit_ADCfrac1_Vmax, &b_bb_gem_hit_ADCfrac1_Vmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCfrac2_Umax", &Ndata_bb_gem_hit_ADCfrac2_Umax, &b_Ndata_bb_gem_hit_ADCfrac2_Umax);
   fChain->SetBranchAddress("bb.gem.hit.ADCfrac2_Umax", bb_gem_hit_ADCfrac2_Umax, &b_bb_gem_hit_ADCfrac2_Umax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCfrac2_Vmax", &Ndata_bb_gem_hit_ADCfrac2_Vmax, &b_Ndata_bb_gem_hit_ADCfrac2_Vmax);
   fChain->SetBranchAddress("bb.gem.hit.ADCfrac2_Vmax", bb_gem_hit_ADCfrac2_Vmax, &b_bb_gem_hit_ADCfrac2_Vmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCfrac3_Umax", &Ndata_bb_gem_hit_ADCfrac3_Umax, &b_Ndata_bb_gem_hit_ADCfrac3_Umax);
   fChain->SetBranchAddress("bb.gem.hit.ADCfrac3_Umax", bb_gem_hit_ADCfrac3_Umax, &b_bb_gem_hit_ADCfrac3_Umax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCfrac3_Vmax", &Ndata_bb_gem_hit_ADCfrac3_Vmax, &b_Ndata_bb_gem_hit_ADCfrac3_Vmax);
   fChain->SetBranchAddress("bb.gem.hit.ADCfrac3_Vmax", bb_gem_hit_ADCfrac3_Vmax, &b_bb_gem_hit_ADCfrac3_Vmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCfrac4_Umax", &Ndata_bb_gem_hit_ADCfrac4_Umax, &b_Ndata_bb_gem_hit_ADCfrac4_Umax);
   fChain->SetBranchAddress("bb.gem.hit.ADCfrac4_Umax", bb_gem_hit_ADCfrac4_Umax, &b_bb_gem_hit_ADCfrac4_Umax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCfrac4_Vmax", &Ndata_bb_gem_hit_ADCfrac4_Vmax, &b_Ndata_bb_gem_hit_ADCfrac4_Vmax);
   fChain->SetBranchAddress("bb.gem.hit.ADCfrac4_Vmax", bb_gem_hit_ADCfrac4_Vmax, &b_bb_gem_hit_ADCfrac4_Vmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCfrac5_Umax", &Ndata_bb_gem_hit_ADCfrac5_Umax, &b_Ndata_bb_gem_hit_ADCfrac5_Umax);
   fChain->SetBranchAddress("bb.gem.hit.ADCfrac5_Umax", bb_gem_hit_ADCfrac5_Umax, &b_bb_gem_hit_ADCfrac5_Umax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCfrac5_Vmax", &Ndata_bb_gem_hit_ADCfrac5_Vmax, &b_Ndata_bb_gem_hit_ADCfrac5_Vmax);
   fChain->SetBranchAddress("bb.gem.hit.ADCfrac5_Vmax", bb_gem_hit_ADCfrac5_Vmax, &b_bb_gem_hit_ADCfrac5_Vmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxcomboUclust_deconv", &Ndata_bb_gem_hit_ADCmaxcomboUclust_deconv, &b_Ndata_bb_gem_hit_ADCmaxcomboUclust_deconv);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxcomboUclust_deconv", bb_gem_hit_ADCmaxcomboUclust_deconv, &b_bb_gem_hit_ADCmaxcomboUclust_deconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxcomboVclust_deconv", &Ndata_bb_gem_hit_ADCmaxcomboVclust_deconv, &b_Ndata_bb_gem_hit_ADCmaxcomboVclust_deconv);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxcomboVclust_deconv", bb_gem_hit_ADCmaxcomboVclust_deconv, &b_bb_gem_hit_ADCmaxcomboVclust_deconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxsampU", &Ndata_bb_gem_hit_ADCmaxsampU, &b_Ndata_bb_gem_hit_ADCmaxsampU);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxsampU", bb_gem_hit_ADCmaxsampU, &b_bb_gem_hit_ADCmaxsampU);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxsampUclust", &Ndata_bb_gem_hit_ADCmaxsampUclust, &b_Ndata_bb_gem_hit_ADCmaxsampUclust);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxsampUclust", bb_gem_hit_ADCmaxsampUclust, &b_bb_gem_hit_ADCmaxsampUclust);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxsampUclust_deconv", &Ndata_bb_gem_hit_ADCmaxsampUclust_deconv, &b_Ndata_bb_gem_hit_ADCmaxsampUclust_deconv);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxsampUclust_deconv", bb_gem_hit_ADCmaxsampUclust_deconv, &b_bb_gem_hit_ADCmaxsampUclust_deconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxsampV", &Ndata_bb_gem_hit_ADCmaxsampV, &b_Ndata_bb_gem_hit_ADCmaxsampV);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxsampV", bb_gem_hit_ADCmaxsampV, &b_bb_gem_hit_ADCmaxsampV);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxsampVclust", &Ndata_bb_gem_hit_ADCmaxsampVclust, &b_Ndata_bb_gem_hit_ADCmaxsampVclust);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxsampVclust", bb_gem_hit_ADCmaxsampVclust, &b_bb_gem_hit_ADCmaxsampVclust);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxsampVclust_deconv", &Ndata_bb_gem_hit_ADCmaxsampVclust_deconv, &b_Ndata_bb_gem_hit_ADCmaxsampVclust_deconv);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxsampVclust_deconv", bb_gem_hit_ADCmaxsampVclust_deconv, &b_bb_gem_hit_ADCmaxsampVclust_deconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxstripU", &Ndata_bb_gem_hit_ADCmaxstripU, &b_Ndata_bb_gem_hit_ADCmaxstripU);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxstripU", bb_gem_hit_ADCmaxstripU, &b_bb_gem_hit_ADCmaxstripU);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxstripV", &Ndata_bb_gem_hit_ADCmaxstripV, &b_Ndata_bb_gem_hit_ADCmaxstripV);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxstripV", bb_gem_hit_ADCmaxstripV, &b_bb_gem_hit_ADCmaxstripV);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.BUILD_ALL_SAMPLES_U", &Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_U, &b_Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_U);
   fChain->SetBranchAddress("bb.gem.hit.BUILD_ALL_SAMPLES_U", bb_gem_hit_BUILD_ALL_SAMPLES_U, &b_bb_gem_hit_BUILD_ALL_SAMPLES_U);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.BUILD_ALL_SAMPLES_V", &Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_V, &b_Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_V);
   fChain->SetBranchAddress("bb.gem.hit.BUILD_ALL_SAMPLES_V", bb_gem_hit_BUILD_ALL_SAMPLES_V, &b_bb_gem_hit_BUILD_ALL_SAMPLES_V);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.CM_GOOD_U", &Ndata_bb_gem_hit_CM_GOOD_U, &b_Ndata_bb_gem_hit_CM_GOOD_U);
   fChain->SetBranchAddress("bb.gem.hit.CM_GOOD_U", bb_gem_hit_CM_GOOD_U, &b_bb_gem_hit_CM_GOOD_U);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.CM_GOOD_V", &Ndata_bb_gem_hit_CM_GOOD_V, &b_Ndata_bb_gem_hit_CM_GOOD_V);
   fChain->SetBranchAddress("bb.gem.hit.CM_GOOD_V", bb_gem_hit_CM_GOOD_V, &b_bb_gem_hit_CM_GOOD_V);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.DeconvADC0_Umax", &Ndata_bb_gem_hit_DeconvADC0_Umax, &b_Ndata_bb_gem_hit_DeconvADC0_Umax);
   fChain->SetBranchAddress("bb.gem.hit.DeconvADC0_Umax", bb_gem_hit_DeconvADC0_Umax, &b_bb_gem_hit_DeconvADC0_Umax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.DeconvADC0_Vmax", &Ndata_bb_gem_hit_DeconvADC0_Vmax, &b_Ndata_bb_gem_hit_DeconvADC0_Vmax);
   fChain->SetBranchAddress("bb.gem.hit.DeconvADC0_Vmax", bb_gem_hit_DeconvADC0_Vmax, &b_bb_gem_hit_DeconvADC0_Vmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.DeconvADC1_Umax", &Ndata_bb_gem_hit_DeconvADC1_Umax, &b_Ndata_bb_gem_hit_DeconvADC1_Umax);
   fChain->SetBranchAddress("bb.gem.hit.DeconvADC1_Umax", bb_gem_hit_DeconvADC1_Umax, &b_bb_gem_hit_DeconvADC1_Umax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.DeconvADC1_Vmax", &Ndata_bb_gem_hit_DeconvADC1_Vmax, &b_Ndata_bb_gem_hit_DeconvADC1_Vmax);
   fChain->SetBranchAddress("bb.gem.hit.DeconvADC1_Vmax", bb_gem_hit_DeconvADC1_Vmax, &b_bb_gem_hit_DeconvADC1_Vmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.DeconvADC2_Umax", &Ndata_bb_gem_hit_DeconvADC2_Umax, &b_Ndata_bb_gem_hit_DeconvADC2_Umax);
   fChain->SetBranchAddress("bb.gem.hit.DeconvADC2_Umax", bb_gem_hit_DeconvADC2_Umax, &b_bb_gem_hit_DeconvADC2_Umax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.DeconvADC2_Vmax", &Ndata_bb_gem_hit_DeconvADC2_Vmax, &b_Ndata_bb_gem_hit_DeconvADC2_Vmax);
   fChain->SetBranchAddress("bb.gem.hit.DeconvADC2_Vmax", bb_gem_hit_DeconvADC2_Vmax, &b_bb_gem_hit_DeconvADC2_Vmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.DeconvADC3_Umax", &Ndata_bb_gem_hit_DeconvADC3_Umax, &b_Ndata_bb_gem_hit_DeconvADC3_Umax);
   fChain->SetBranchAddress("bb.gem.hit.DeconvADC3_Umax", bb_gem_hit_DeconvADC3_Umax, &b_bb_gem_hit_DeconvADC3_Umax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.DeconvADC3_Vmax", &Ndata_bb_gem_hit_DeconvADC3_Vmax, &b_Ndata_bb_gem_hit_DeconvADC3_Vmax);
   fChain->SetBranchAddress("bb.gem.hit.DeconvADC3_Vmax", bb_gem_hit_DeconvADC3_Vmax, &b_bb_gem_hit_DeconvADC3_Vmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.DeconvADC4_Umax", &Ndata_bb_gem_hit_DeconvADC4_Umax, &b_Ndata_bb_gem_hit_DeconvADC4_Umax);
   fChain->SetBranchAddress("bb.gem.hit.DeconvADC4_Umax", bb_gem_hit_DeconvADC4_Umax, &b_bb_gem_hit_DeconvADC4_Umax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.DeconvADC4_Vmax", &Ndata_bb_gem_hit_DeconvADC4_Vmax, &b_Ndata_bb_gem_hit_DeconvADC4_Vmax);
   fChain->SetBranchAddress("bb.gem.hit.DeconvADC4_Vmax", bb_gem_hit_DeconvADC4_Vmax, &b_bb_gem_hit_DeconvADC4_Vmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.DeconvADC5_Umax", &Ndata_bb_gem_hit_DeconvADC5_Umax, &b_Ndata_bb_gem_hit_DeconvADC5_Umax);
   fChain->SetBranchAddress("bb.gem.hit.DeconvADC5_Umax", bb_gem_hit_DeconvADC5_Umax, &b_bb_gem_hit_DeconvADC5_Umax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.DeconvADC5_Vmax", &Ndata_bb_gem_hit_DeconvADC5_Vmax, &b_Ndata_bb_gem_hit_DeconvADC5_Vmax);
   fChain->SetBranchAddress("bb.gem.hit.DeconvADC5_Vmax", bb_gem_hit_DeconvADC5_Vmax, &b_bb_gem_hit_DeconvADC5_Vmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.DeconvADCmaxcomboU", &Ndata_bb_gem_hit_DeconvADCmaxcomboU, &b_Ndata_bb_gem_hit_DeconvADCmaxcomboU);
   fChain->SetBranchAddress("bb.gem.hit.DeconvADCmaxcomboU", bb_gem_hit_DeconvADCmaxcomboU, &b_bb_gem_hit_DeconvADCmaxcomboU);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.DeconvADCmaxcomboV", &Ndata_bb_gem_hit_DeconvADCmaxcomboV, &b_Ndata_bb_gem_hit_DeconvADCmaxcomboV);
   fChain->SetBranchAddress("bb.gem.hit.DeconvADCmaxcomboV", bb_gem_hit_DeconvADCmaxcomboV, &b_bb_gem_hit_DeconvADCmaxcomboV);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.DeconvADCmaxsampU", &Ndata_bb_gem_hit_DeconvADCmaxsampU, &b_Ndata_bb_gem_hit_DeconvADCmaxsampU);
   fChain->SetBranchAddress("bb.gem.hit.DeconvADCmaxsampU", bb_gem_hit_DeconvADCmaxsampU, &b_bb_gem_hit_DeconvADCmaxsampU);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.DeconvADCmaxsampV", &Ndata_bb_gem_hit_DeconvADCmaxsampV, &b_Ndata_bb_gem_hit_DeconvADCmaxsampV);
   fChain->SetBranchAddress("bb.gem.hit.DeconvADCmaxsampV", bb_gem_hit_DeconvADCmaxsampV, &b_bb_gem_hit_DeconvADCmaxsampV);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.DeconvADCmaxstripU", &Ndata_bb_gem_hit_DeconvADCmaxstripU, &b_Ndata_bb_gem_hit_DeconvADCmaxstripU);
   fChain->SetBranchAddress("bb.gem.hit.DeconvADCmaxstripU", bb_gem_hit_DeconvADCmaxstripU, &b_bb_gem_hit_DeconvADCmaxstripU);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.DeconvADCmaxstripV", &Ndata_bb_gem_hit_DeconvADCmaxstripV, &b_Ndata_bb_gem_hit_DeconvADCmaxstripV);
   fChain->SetBranchAddress("bb.gem.hit.DeconvADCmaxstripV", bb_gem_hit_DeconvADCmaxstripV, &b_bb_gem_hit_DeconvADCmaxstripV);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ENABLE_CM_U", &Ndata_bb_gem_hit_ENABLE_CM_U, &b_Ndata_bb_gem_hit_ENABLE_CM_U);
   fChain->SetBranchAddress("bb.gem.hit.ENABLE_CM_U", bb_gem_hit_ENABLE_CM_U, &b_bb_gem_hit_ENABLE_CM_U);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ENABLE_CM_V", &Ndata_bb_gem_hit_ENABLE_CM_V, &b_Ndata_bb_gem_hit_ENABLE_CM_V);
   fChain->SetBranchAddress("bb.gem.hit.ENABLE_CM_V", bb_gem_hit_ENABLE_CM_V, &b_bb_gem_hit_ENABLE_CM_V);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.TSchi2_Umax", &Ndata_bb_gem_hit_TSchi2_Umax, &b_Ndata_bb_gem_hit_TSchi2_Umax);
   fChain->SetBranchAddress("bb.gem.hit.TSchi2_Umax", bb_gem_hit_TSchi2_Umax, &b_bb_gem_hit_TSchi2_Umax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.TSchi2_Vmax", &Ndata_bb_gem_hit_TSchi2_Vmax, &b_Ndata_bb_gem_hit_TSchi2_Vmax);
   fChain->SetBranchAddress("bb.gem.hit.TSchi2_Vmax", bb_gem_hit_TSchi2_Vmax, &b_bb_gem_hit_TSchi2_Vmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.TSprob_Umax", &Ndata_bb_gem_hit_TSprob_Umax, &b_Ndata_bb_gem_hit_TSprob_Umax);
   fChain->SetBranchAddress("bb.gem.hit.TSprob_Umax", bb_gem_hit_TSprob_Umax, &b_bb_gem_hit_TSprob_Umax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.TSprob_Vmax", &Ndata_bb_gem_hit_TSprob_Vmax, &b_Ndata_bb_gem_hit_TSprob_Vmax);
   fChain->SetBranchAddress("bb.gem.hit.TSprob_Vmax", bb_gem_hit_TSprob_Vmax, &b_bb_gem_hit_TSprob_Vmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.Tavg", &Ndata_bb_gem_hit_Tavg, &b_Ndata_bb_gem_hit_Tavg);
   fChain->SetBranchAddress("bb.gem.hit.Tavg", bb_gem_hit_Tavg, &b_bb_gem_hit_Tavg);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.Tavg_corr", &Ndata_bb_gem_hit_Tavg_corr, &b_Ndata_bb_gem_hit_Tavg_corr);
   fChain->SetBranchAddress("bb.gem.hit.Tavg_corr", bb_gem_hit_Tavg_corr, &b_bb_gem_hit_Tavg_corr);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.Tavg_deconv", &Ndata_bb_gem_hit_Tavg_deconv, &b_Ndata_bb_gem_hit_Tavg_deconv);
   fChain->SetBranchAddress("bb.gem.hit.Tavg_deconv", bb_gem_hit_Tavg_deconv, &b_bb_gem_hit_Tavg_deconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.Tavg_fit", &Ndata_bb_gem_hit_Tavg_fit, &b_Ndata_bb_gem_hit_Tavg_fit);
   fChain->SetBranchAddress("bb.gem.hit.Tavg_fit", bb_gem_hit_Tavg_fit, &b_bb_gem_hit_Tavg_fit);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.Ugain", &Ndata_bb_gem_hit_Ugain, &b_Ndata_bb_gem_hit_Ugain);
   fChain->SetBranchAddress("bb.gem.hit.Ugain", bb_gem_hit_Ugain, &b_bb_gem_hit_Ugain);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.Utime", &Ndata_bb_gem_hit_Utime, &b_Ndata_bb_gem_hit_Utime);
   fChain->SetBranchAddress("bb.gem.hit.Utime", bb_gem_hit_Utime, &b_bb_gem_hit_Utime);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.UtimeDeconv", &Ndata_bb_gem_hit_UtimeDeconv, &b_Ndata_bb_gem_hit_UtimeDeconv);
   fChain->SetBranchAddress("bb.gem.hit.UtimeDeconv", bb_gem_hit_UtimeDeconv, &b_bb_gem_hit_UtimeDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.UtimeFit", &Ndata_bb_gem_hit_UtimeFit, &b_Ndata_bb_gem_hit_UtimeFit);
   fChain->SetBranchAddress("bb.gem.hit.UtimeFit", bb_gem_hit_UtimeFit, &b_bb_gem_hit_UtimeFit);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.UtimeMaxStrip", &Ndata_bb_gem_hit_UtimeMaxStrip, &b_Ndata_bb_gem_hit_UtimeMaxStrip);
   fChain->SetBranchAddress("bb.gem.hit.UtimeMaxStrip", bb_gem_hit_UtimeMaxStrip, &b_bb_gem_hit_UtimeMaxStrip);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.UtimeMaxStripDeconv", &Ndata_bb_gem_hit_UtimeMaxStripDeconv, &b_Ndata_bb_gem_hit_UtimeMaxStripDeconv);
   fChain->SetBranchAddress("bb.gem.hit.UtimeMaxStripDeconv", bb_gem_hit_UtimeMaxStripDeconv, &b_bb_gem_hit_UtimeMaxStripDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.UtimeMaxStripFit", &Ndata_bb_gem_hit_UtimeMaxStripFit, &b_Ndata_bb_gem_hit_UtimeMaxStripFit);
   fChain->SetBranchAddress("bb.gem.hit.UtimeMaxStripFit", bb_gem_hit_UtimeMaxStripFit, &b_bb_gem_hit_UtimeMaxStripFit);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.Vgain", &Ndata_bb_gem_hit_Vgain, &b_Ndata_bb_gem_hit_Vgain);
   fChain->SetBranchAddress("bb.gem.hit.Vgain", bb_gem_hit_Vgain, &b_bb_gem_hit_Vgain);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.Vtime", &Ndata_bb_gem_hit_Vtime, &b_Ndata_bb_gem_hit_Vtime);
   fChain->SetBranchAddress("bb.gem.hit.Vtime", bb_gem_hit_Vtime, &b_bb_gem_hit_Vtime);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.VtimeDeconv", &Ndata_bb_gem_hit_VtimeDeconv, &b_Ndata_bb_gem_hit_VtimeDeconv);
   fChain->SetBranchAddress("bb.gem.hit.VtimeDeconv", bb_gem_hit_VtimeDeconv, &b_bb_gem_hit_VtimeDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.VtimeFit", &Ndata_bb_gem_hit_VtimeFit, &b_Ndata_bb_gem_hit_VtimeFit);
   fChain->SetBranchAddress("bb.gem.hit.VtimeFit", bb_gem_hit_VtimeFit, &b_bb_gem_hit_VtimeFit);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.VtimeMaxStrip", &Ndata_bb_gem_hit_VtimeMaxStrip, &b_Ndata_bb_gem_hit_VtimeMaxStrip);
   fChain->SetBranchAddress("bb.gem.hit.VtimeMaxStrip", bb_gem_hit_VtimeMaxStrip, &b_bb_gem_hit_VtimeMaxStrip);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.VtimeMaxStripDeconv", &Ndata_bb_gem_hit_VtimeMaxStripDeconv, &b_Ndata_bb_gem_hit_VtimeMaxStripDeconv);
   fChain->SetBranchAddress("bb.gem.hit.VtimeMaxStripDeconv", bb_gem_hit_VtimeMaxStripDeconv, &b_bb_gem_hit_VtimeMaxStripDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.VtimeMaxStripFit", &Ndata_bb_gem_hit_VtimeMaxStripFit, &b_Ndata_bb_gem_hit_VtimeMaxStripFit);
   fChain->SetBranchAddress("bb.gem.hit.VtimeMaxStripFit", bb_gem_hit_VtimeMaxStripFit, &b_bb_gem_hit_VtimeMaxStripFit);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ccor_clust", &Ndata_bb_gem_hit_ccor_clust, &b_Ndata_bb_gem_hit_ccor_clust);
   fChain->SetBranchAddress("bb.gem.hit.ccor_clust", bb_gem_hit_ccor_clust, &b_bb_gem_hit_ccor_clust);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ccor_clust_deconv", &Ndata_bb_gem_hit_ccor_clust_deconv, &b_Ndata_bb_gem_hit_ccor_clust_deconv);
   fChain->SetBranchAddress("bb.gem.hit.ccor_clust_deconv", bb_gem_hit_ccor_clust_deconv, &b_bb_gem_hit_ccor_clust_deconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ccor_strip", &Ndata_bb_gem_hit_ccor_strip, &b_Ndata_bb_gem_hit_ccor_strip);
   fChain->SetBranchAddress("bb.gem.hit.ccor_strip", bb_gem_hit_ccor_strip, &b_bb_gem_hit_ccor_strip);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ccor_strip_deconv", &Ndata_bb_gem_hit_ccor_strip_deconv, &b_Ndata_bb_gem_hit_ccor_strip_deconv);
   fChain->SetBranchAddress("bb.gem.hit.ccor_strip_deconv", bb_gem_hit_ccor_strip_deconv, &b_bb_gem_hit_ccor_strip_deconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.deltat", &Ndata_bb_gem_hit_deltat, &b_Ndata_bb_gem_hit_deltat);
   fChain->SetBranchAddress("bb.gem.hit.deltat", bb_gem_hit_deltat, &b_bb_gem_hit_deltat);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.deltat_deconv", &Ndata_bb_gem_hit_deltat_deconv, &b_Ndata_bb_gem_hit_deltat_deconv);
   fChain->SetBranchAddress("bb.gem.hit.deltat_deconv", bb_gem_hit_deltat_deconv, &b_bb_gem_hit_deltat_deconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.deltat_fit", &Ndata_bb_gem_hit_deltat_fit, &b_Ndata_bb_gem_hit_deltat_fit);
   fChain->SetBranchAddress("bb.gem.hit.deltat_fit", bb_gem_hit_deltat_fit, &b_bb_gem_hit_deltat_fit);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.eresidu", &Ndata_bb_gem_hit_eresidu, &b_Ndata_bb_gem_hit_eresidu);
   fChain->SetBranchAddress("bb.gem.hit.eresidu", bb_gem_hit_eresidu, &b_bb_gem_hit_eresidu);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.eresidv", &Ndata_bb_gem_hit_eresidv, &b_Ndata_bb_gem_hit_eresidv);
   fChain->SetBranchAddress("bb.gem.hit.eresidv", bb_gem_hit_eresidv, &b_bb_gem_hit_eresidv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.icombomaxUclustDeconv", &Ndata_bb_gem_hit_icombomaxUclustDeconv, &b_Ndata_bb_gem_hit_icombomaxUclustDeconv);
   fChain->SetBranchAddress("bb.gem.hit.icombomaxUclustDeconv", bb_gem_hit_icombomaxUclustDeconv, &b_bb_gem_hit_icombomaxUclustDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.icombomaxUstripDeconv", &Ndata_bb_gem_hit_icombomaxUstripDeconv, &b_Ndata_bb_gem_hit_icombomaxUstripDeconv);
   fChain->SetBranchAddress("bb.gem.hit.icombomaxUstripDeconv", bb_gem_hit_icombomaxUstripDeconv, &b_bb_gem_hit_icombomaxUstripDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.icombomaxVclustDeconv", &Ndata_bb_gem_hit_icombomaxVclustDeconv, &b_Ndata_bb_gem_hit_icombomaxVclustDeconv);
   fChain->SetBranchAddress("bb.gem.hit.icombomaxVclustDeconv", bb_gem_hit_icombomaxVclustDeconv, &b_bb_gem_hit_icombomaxVclustDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.icombomaxVstripDeconv", &Ndata_bb_gem_hit_icombomaxVstripDeconv, &b_Ndata_bb_gem_hit_icombomaxVstripDeconv);
   fChain->SetBranchAddress("bb.gem.hit.icombomaxVstripDeconv", bb_gem_hit_icombomaxVstripDeconv, &b_bb_gem_hit_icombomaxVstripDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.isampmaxUclust", &Ndata_bb_gem_hit_isampmaxUclust, &b_Ndata_bb_gem_hit_isampmaxUclust);
   fChain->SetBranchAddress("bb.gem.hit.isampmaxUclust", bb_gem_hit_isampmaxUclust, &b_bb_gem_hit_isampmaxUclust);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.isampmaxUclustDeconv", &Ndata_bb_gem_hit_isampmaxUclustDeconv, &b_Ndata_bb_gem_hit_isampmaxUclustDeconv);
   fChain->SetBranchAddress("bb.gem.hit.isampmaxUclustDeconv", bb_gem_hit_isampmaxUclustDeconv, &b_bb_gem_hit_isampmaxUclustDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.isampmaxUstrip", &Ndata_bb_gem_hit_isampmaxUstrip, &b_Ndata_bb_gem_hit_isampmaxUstrip);
   fChain->SetBranchAddress("bb.gem.hit.isampmaxUstrip", bb_gem_hit_isampmaxUstrip, &b_bb_gem_hit_isampmaxUstrip);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.isampmaxUstripDeconv", &Ndata_bb_gem_hit_isampmaxUstripDeconv, &b_Ndata_bb_gem_hit_isampmaxUstripDeconv);
   fChain->SetBranchAddress("bb.gem.hit.isampmaxUstripDeconv", bb_gem_hit_isampmaxUstripDeconv, &b_bb_gem_hit_isampmaxUstripDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.isampmaxVclust", &Ndata_bb_gem_hit_isampmaxVclust, &b_Ndata_bb_gem_hit_isampmaxVclust);
   fChain->SetBranchAddress("bb.gem.hit.isampmaxVclust", bb_gem_hit_isampmaxVclust, &b_bb_gem_hit_isampmaxVclust);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.isampmaxVclustDeconv", &Ndata_bb_gem_hit_isampmaxVclustDeconv, &b_Ndata_bb_gem_hit_isampmaxVclustDeconv);
   fChain->SetBranchAddress("bb.gem.hit.isampmaxVclustDeconv", bb_gem_hit_isampmaxVclustDeconv, &b_bb_gem_hit_isampmaxVclustDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.isampmaxVstrip", &Ndata_bb_gem_hit_isampmaxVstrip, &b_Ndata_bb_gem_hit_isampmaxVstrip);
   fChain->SetBranchAddress("bb.gem.hit.isampmaxVstrip", bb_gem_hit_isampmaxVstrip, &b_bb_gem_hit_isampmaxVstrip);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.isampmaxVstripDeconv", &Ndata_bb_gem_hit_isampmaxVstripDeconv, &b_Ndata_bb_gem_hit_isampmaxVstripDeconv);
   fChain->SetBranchAddress("bb.gem.hit.isampmaxVstripDeconv", bb_gem_hit_isampmaxVstripDeconv, &b_bb_gem_hit_isampmaxVstripDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.layer", &Ndata_bb_gem_hit_layer, &b_Ndata_bb_gem_hit_layer);
   fChain->SetBranchAddress("bb.gem.hit.layer", bb_gem_hit_layer, &b_bb_gem_hit_layer);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.module", &Ndata_bb_gem_hit_module, &b_Ndata_bb_gem_hit_module);
   fChain->SetBranchAddress("bb.gem.hit.module", bb_gem_hit_module, &b_bb_gem_hit_module);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.nstripu", &Ndata_bb_gem_hit_nstripu, &b_Ndata_bb_gem_hit_nstripu);
   fChain->SetBranchAddress("bb.gem.hit.nstripu", bb_gem_hit_nstripu, &b_bb_gem_hit_nstripu);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.nstripv", &Ndata_bb_gem_hit_nstripv, &b_Ndata_bb_gem_hit_nstripv);
   fChain->SetBranchAddress("bb.gem.hit.nstripv", bb_gem_hit_nstripv, &b_bb_gem_hit_nstripv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.residu", &Ndata_bb_gem_hit_residu, &b_Ndata_bb_gem_hit_residu);
   fChain->SetBranchAddress("bb.gem.hit.residu", bb_gem_hit_residu, &b_bb_gem_hit_residu);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.residv", &Ndata_bb_gem_hit_residv, &b_Ndata_bb_gem_hit_residv);
   fChain->SetBranchAddress("bb.gem.hit.residv", bb_gem_hit_residv, &b_bb_gem_hit_residv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.trackindex", &Ndata_bb_gem_hit_trackindex, &b_Ndata_bb_gem_hit_trackindex);
   fChain->SetBranchAddress("bb.gem.hit.trackindex", bb_gem_hit_trackindex, &b_bb_gem_hit_trackindex);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.u", &Ndata_bb_gem_hit_u, &b_Ndata_bb_gem_hit_u);
   fChain->SetBranchAddress("bb.gem.hit.u", bb_gem_hit_u, &b_bb_gem_hit_u);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.umoment", &Ndata_bb_gem_hit_umoment, &b_Ndata_bb_gem_hit_umoment);
   fChain->SetBranchAddress("bb.gem.hit.umoment", bb_gem_hit_umoment, &b_bb_gem_hit_umoment);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.usigma", &Ndata_bb_gem_hit_usigma, &b_Ndata_bb_gem_hit_usigma);
   fChain->SetBranchAddress("bb.gem.hit.usigma", bb_gem_hit_usigma, &b_bb_gem_hit_usigma);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ustriphi", &Ndata_bb_gem_hit_ustriphi, &b_Ndata_bb_gem_hit_ustriphi);
   fChain->SetBranchAddress("bb.gem.hit.ustriphi", bb_gem_hit_ustriphi, &b_bb_gem_hit_ustriphi);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ustriplo", &Ndata_bb_gem_hit_ustriplo, &b_Ndata_bb_gem_hit_ustriplo);
   fChain->SetBranchAddress("bb.gem.hit.ustriplo", bb_gem_hit_ustriplo, &b_bb_gem_hit_ustriplo);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ustripmax", &Ndata_bb_gem_hit_ustripmax, &b_Ndata_bb_gem_hit_ustripmax);
   fChain->SetBranchAddress("bb.gem.hit.ustripmax", bb_gem_hit_ustripmax, &b_bb_gem_hit_ustripmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.v", &Ndata_bb_gem_hit_v, &b_Ndata_bb_gem_hit_v);
   fChain->SetBranchAddress("bb.gem.hit.v", bb_gem_hit_v, &b_bb_gem_hit_v);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.vmoment", &Ndata_bb_gem_hit_vmoment, &b_Ndata_bb_gem_hit_vmoment);
   fChain->SetBranchAddress("bb.gem.hit.vmoment", bb_gem_hit_vmoment, &b_bb_gem_hit_vmoment);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.vsigma", &Ndata_bb_gem_hit_vsigma, &b_Ndata_bb_gem_hit_vsigma);
   fChain->SetBranchAddress("bb.gem.hit.vsigma", bb_gem_hit_vsigma, &b_bb_gem_hit_vsigma);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.vstriphi", &Ndata_bb_gem_hit_vstriphi, &b_Ndata_bb_gem_hit_vstriphi);
   fChain->SetBranchAddress("bb.gem.hit.vstriphi", bb_gem_hit_vstriphi, &b_bb_gem_hit_vstriphi);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.vstriplo", &Ndata_bb_gem_hit_vstriplo, &b_Ndata_bb_gem_hit_vstriplo);
   fChain->SetBranchAddress("bb.gem.hit.vstriplo", bb_gem_hit_vstriplo, &b_bb_gem_hit_vstriplo);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.vstripmax", &Ndata_bb_gem_hit_vstripmax, &b_Ndata_bb_gem_hit_vstripmax);
   fChain->SetBranchAddress("bb.gem.hit.vstripmax", bb_gem_hit_vstripmax, &b_bb_gem_hit_vstripmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.xglobal", &Ndata_bb_gem_hit_xglobal, &b_Ndata_bb_gem_hit_xglobal);
   fChain->SetBranchAddress("bb.gem.hit.xglobal", bb_gem_hit_xglobal, &b_bb_gem_hit_xglobal);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.xlocal", &Ndata_bb_gem_hit_xlocal, &b_Ndata_bb_gem_hit_xlocal);
   fChain->SetBranchAddress("bb.gem.hit.xlocal", bb_gem_hit_xlocal, &b_bb_gem_hit_xlocal);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.yglobal", &Ndata_bb_gem_hit_yglobal, &b_Ndata_bb_gem_hit_yglobal);
   fChain->SetBranchAddress("bb.gem.hit.yglobal", bb_gem_hit_yglobal, &b_bb_gem_hit_yglobal);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ylocal", &Ndata_bb_gem_hit_ylocal, &b_Ndata_bb_gem_hit_ylocal);
   fChain->SetBranchAddress("bb.gem.hit.ylocal", bb_gem_hit_ylocal, &b_bb_gem_hit_ylocal);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.zglobal", &Ndata_bb_gem_hit_zglobal, &b_Ndata_bb_gem_hit_zglobal);
   fChain->SetBranchAddress("bb.gem.hit.zglobal", bb_gem_hit_zglobal, &b_bb_gem_hit_zglobal);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.clust.clustu_adc", &Ndata_bb_gem_m0_clust_clustu_adc, &b_Ndata_bb_gem_m0_clust_clustu_adc);
   fChain->SetBranchAddress("bb.gem.m0.clust.clustu_adc", bb_gem_m0_clust_clustu_adc, &b_bb_gem_m0_clust_clustu_adc);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.clust.clustu_pos", &Ndata_bb_gem_m0_clust_clustu_pos, &b_Ndata_bb_gem_m0_clust_clustu_pos);
   fChain->SetBranchAddress("bb.gem.m0.clust.clustu_pos", bb_gem_m0_clust_clustu_pos, &b_bb_gem_m0_clust_clustu_pos);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.clust.clustu_strips", &Ndata_bb_gem_m0_clust_clustu_strips, &b_Ndata_bb_gem_m0_clust_clustu_strips);
   fChain->SetBranchAddress("bb.gem.m0.clust.clustu_strips", bb_gem_m0_clust_clustu_strips, &b_bb_gem_m0_clust_clustu_strips);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.clust.clustu_time", &Ndata_bb_gem_m0_clust_clustu_time, &b_Ndata_bb_gem_m0_clust_clustu_time);
   fChain->SetBranchAddress("bb.gem.m0.clust.clustu_time", bb_gem_m0_clust_clustu_time, &b_bb_gem_m0_clust_clustu_time);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.clust.clustv_adc", &Ndata_bb_gem_m0_clust_clustv_adc, &b_Ndata_bb_gem_m0_clust_clustv_adc);
   fChain->SetBranchAddress("bb.gem.m0.clust.clustv_adc", bb_gem_m0_clust_clustv_adc, &b_bb_gem_m0_clust_clustv_adc);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.clust.clustv_pos", &Ndata_bb_gem_m0_clust_clustv_pos, &b_Ndata_bb_gem_m0_clust_clustv_pos);
   fChain->SetBranchAddress("bb.gem.m0.clust.clustv_pos", bb_gem_m0_clust_clustv_pos, &b_bb_gem_m0_clust_clustv_pos);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.clust.clustv_strips", &Ndata_bb_gem_m0_clust_clustv_strips, &b_Ndata_bb_gem_m0_clust_clustv_strips);
   fChain->SetBranchAddress("bb.gem.m0.clust.clustv_strips", bb_gem_m0_clust_clustv_strips, &b_bb_gem_m0_clust_clustv_strips);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.clust.clustv_time", &Ndata_bb_gem_m0_clust_clustv_time, &b_Ndata_bb_gem_m0_clust_clustv_time);
   fChain->SetBranchAddress("bb.gem.m0.clust.clustv_time", bb_gem_m0_clust_clustv_time, &b_bb_gem_m0_clust_clustv_time);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.clust.isnegativeU", &Ndata_bb_gem_m0_clust_isnegativeU, &b_Ndata_bb_gem_m0_clust_isnegativeU);
   fChain->SetBranchAddress("bb.gem.m0.clust.isnegativeU", bb_gem_m0_clust_isnegativeU, &b_bb_gem_m0_clust_isnegativeU);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.clust.isnegativeV", &Ndata_bb_gem_m0_clust_isnegativeV, &b_Ndata_bb_gem_m0_clust_isnegativeV);
   fChain->SetBranchAddress("bb.gem.m0.clust.isnegativeV", bb_gem_m0_clust_isnegativeV, &b_bb_gem_m0_clust_isnegativeV);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.clust.isnegontrackU", &Ndata_bb_gem_m0_clust_isnegontrackU, &b_Ndata_bb_gem_m0_clust_isnegontrackU);
   fChain->SetBranchAddress("bb.gem.m0.clust.isnegontrackU", bb_gem_m0_clust_isnegontrackU, &b_bb_gem_m0_clust_isnegontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.clust.isnegontrackV", &Ndata_bb_gem_m0_clust_isnegontrackV, &b_Ndata_bb_gem_m0_clust_isnegontrackV);
   fChain->SetBranchAddress("bb.gem.m0.clust.isnegontrackV", bb_gem_m0_clust_isnegontrackV, &b_bb_gem_m0_clust_isnegontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.ADCavg", &Ndata_bb_gem_m0_strip_ADCavg, &b_Ndata_bb_gem_m0_strip_ADCavg);
   fChain->SetBranchAddress("bb.gem.m0.strip.ADCavg", bb_gem_m0_strip_ADCavg, &b_bb_gem_m0_strip_ADCavg);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.ADCmax", &Ndata_bb_gem_m0_strip_ADCmax, &b_Ndata_bb_gem_m0_strip_ADCmax);
   fChain->SetBranchAddress("bb.gem.m0.strip.ADCmax", bb_gem_m0_strip_ADCmax, &b_bb_gem_m0_strip_ADCmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.ADCsamples", &Ndata_bb_gem_m0_strip_ADCsamples, &b_Ndata_bb_gem_m0_strip_ADCsamples);
   fChain->SetBranchAddress("bb.gem.m0.strip.ADCsamples", bb_gem_m0_strip_ADCsamples, &b_bb_gem_m0_strip_ADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.ADCsum", &Ndata_bb_gem_m0_strip_ADCsum, &b_Ndata_bb_gem_m0_strip_ADCsum);
   fChain->SetBranchAddress("bb.gem.m0.strip.ADCsum", bb_gem_m0_strip_ADCsum, &b_bb_gem_m0_strip_ADCsum);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.BUILD_ALL_SAMPLES", &Ndata_bb_gem_m0_strip_BUILD_ALL_SAMPLES, &b_Ndata_bb_gem_m0_strip_BUILD_ALL_SAMPLES);
   fChain->SetBranchAddress("bb.gem.m0.strip.BUILD_ALL_SAMPLES", bb_gem_m0_strip_BUILD_ALL_SAMPLES, &b_bb_gem_m0_strip_BUILD_ALL_SAMPLES);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.CM_GOOD", &Ndata_bb_gem_m0_strip_CM_GOOD, &b_Ndata_bb_gem_m0_strip_CM_GOOD);
   fChain->SetBranchAddress("bb.gem.m0.strip.CM_GOOD", bb_gem_m0_strip_CM_GOOD, &b_bb_gem_m0_strip_CM_GOOD);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.CorrCoeff", &Ndata_bb_gem_m0_strip_CorrCoeff, &b_Ndata_bb_gem_m0_strip_CorrCoeff);
   fChain->SetBranchAddress("bb.gem.m0.strip.CorrCoeff", bb_gem_m0_strip_CorrCoeff, &b_bb_gem_m0_strip_CorrCoeff);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.DeconvADCmax", &Ndata_bb_gem_m0_strip_DeconvADCmax, &b_Ndata_bb_gem_m0_strip_DeconvADCmax);
   fChain->SetBranchAddress("bb.gem.m0.strip.DeconvADCmax", bb_gem_m0_strip_DeconvADCmax, &b_bb_gem_m0_strip_DeconvADCmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.DeconvADCmaxCombo", &Ndata_bb_gem_m0_strip_DeconvADCmaxCombo, &b_Ndata_bb_gem_m0_strip_DeconvADCmaxCombo);
   fChain->SetBranchAddress("bb.gem.m0.strip.DeconvADCmaxCombo", bb_gem_m0_strip_DeconvADCmaxCombo, &b_bb_gem_m0_strip_DeconvADCmaxCombo);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.DeconvADCsamples", &Ndata_bb_gem_m0_strip_DeconvADCsamples, &b_Ndata_bb_gem_m0_strip_DeconvADCsamples);
   fChain->SetBranchAddress("bb.gem.m0.strip.DeconvADCsamples", bb_gem_m0_strip_DeconvADCsamples, &b_bb_gem_m0_strip_DeconvADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.DeconvADCsum", &Ndata_bb_gem_m0_strip_DeconvADCsum, &b_Ndata_bb_gem_m0_strip_DeconvADCsum);
   fChain->SetBranchAddress("bb.gem.m0.strip.DeconvADCsum", bb_gem_m0_strip_DeconvADCsum, &b_bb_gem_m0_strip_DeconvADCsum);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.ENABLE_CM", &Ndata_bb_gem_m0_strip_ENABLE_CM, &b_Ndata_bb_gem_m0_strip_ENABLE_CM);
   fChain->SetBranchAddress("bb.gem.m0.strip.ENABLE_CM", bb_gem_m0_strip_ENABLE_CM, &b_bb_gem_m0_strip_ENABLE_CM);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.IsU", &Ndata_bb_gem_m0_strip_IsU, &b_Ndata_bb_gem_m0_strip_IsU);
   fChain->SetBranchAddress("bb.gem.m0.strip.IsU", bb_gem_m0_strip_IsU, &b_bb_gem_m0_strip_IsU);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.IsV", &Ndata_bb_gem_m0_strip_IsV, &b_Ndata_bb_gem_m0_strip_IsV);
   fChain->SetBranchAddress("bb.gem.m0.strip.IsV", bb_gem_m0_strip_IsV, &b_bb_gem_m0_strip_IsV);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.TSchi2", &Ndata_bb_gem_m0_strip_TSchi2, &b_Ndata_bb_gem_m0_strip_TSchi2);
   fChain->SetBranchAddress("bb.gem.m0.strip.TSchi2", bb_gem_m0_strip_TSchi2, &b_bb_gem_m0_strip_TSchi2);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.TSprob", &Ndata_bb_gem_m0_strip_TSprob, &b_Ndata_bb_gem_m0_strip_TSprob);
   fChain->SetBranchAddress("bb.gem.m0.strip.TSprob", bb_gem_m0_strip_TSprob, &b_bb_gem_m0_strip_TSprob);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.Tcorr", &Ndata_bb_gem_m0_strip_Tcorr, &b_Ndata_bb_gem_m0_strip_Tcorr);
   fChain->SetBranchAddress("bb.gem.m0.strip.Tcorr", bb_gem_m0_strip_Tcorr, &b_bb_gem_m0_strip_Tcorr);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.Tdiff", &Ndata_bb_gem_m0_strip_Tdiff, &b_Ndata_bb_gem_m0_strip_Tdiff);
   fChain->SetBranchAddress("bb.gem.m0.strip.Tdiff", bb_gem_m0_strip_Tdiff, &b_bb_gem_m0_strip_Tdiff);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.Tfit", &Ndata_bb_gem_m0_strip_Tfit, &b_Ndata_bb_gem_m0_strip_Tfit);
   fChain->SetBranchAddress("bb.gem.m0.strip.Tfit", bb_gem_m0_strip_Tfit, &b_bb_gem_m0_strip_Tfit);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.Tmean", &Ndata_bb_gem_m0_strip_Tmean, &b_Ndata_bb_gem_m0_strip_Tmean);
   fChain->SetBranchAddress("bb.gem.m0.strip.Tmean", bb_gem_m0_strip_Tmean, &b_bb_gem_m0_strip_Tmean);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.TmeanDeconv", &Ndata_bb_gem_m0_strip_TmeanDeconv, &b_Ndata_bb_gem_m0_strip_TmeanDeconv);
   fChain->SetBranchAddress("bb.gem.m0.strip.TmeanDeconv", bb_gem_m0_strip_TmeanDeconv, &b_bb_gem_m0_strip_TmeanDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.Tsigma", &Ndata_bb_gem_m0_strip_Tsigma, &b_Ndata_bb_gem_m0_strip_Tsigma);
   fChain->SetBranchAddress("bb.gem.m0.strip.Tsigma", bb_gem_m0_strip_Tsigma, &b_bb_gem_m0_strip_Tsigma);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.adc_id", &Ndata_bb_gem_m0_strip_adc_id, &b_Ndata_bb_gem_m0_strip_adc_id);
   fChain->SetBranchAddress("bb.gem.m0.strip.adc_id", bb_gem_m0_strip_adc_id, &b_bb_gem_m0_strip_adc_id);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.crate", &Ndata_bb_gem_m0_strip_crate, &b_Ndata_bb_gem_m0_strip_crate);
   fChain->SetBranchAddress("bb.gem.m0.strip.crate", bb_gem_m0_strip_crate, &b_bb_gem_m0_strip_crate);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.isampmax", &Ndata_bb_gem_m0_strip_isampmax, &b_Ndata_bb_gem_m0_strip_isampmax);
   fChain->SetBranchAddress("bb.gem.m0.strip.isampmax", bb_gem_m0_strip_isampmax, &b_bb_gem_m0_strip_isampmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.isampmaxDeconv", &Ndata_bb_gem_m0_strip_isampmaxDeconv, &b_Ndata_bb_gem_m0_strip_isampmaxDeconv);
   fChain->SetBranchAddress("bb.gem.m0.strip.isampmaxDeconv", bb_gem_m0_strip_isampmaxDeconv, &b_bb_gem_m0_strip_isampmaxDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.isampmaxDeconvCombo", &Ndata_bb_gem_m0_strip_isampmaxDeconvCombo, &b_Ndata_bb_gem_m0_strip_isampmaxDeconvCombo);
   fChain->SetBranchAddress("bb.gem.m0.strip.isampmaxDeconvCombo", bb_gem_m0_strip_isampmaxDeconvCombo, &b_bb_gem_m0_strip_isampmaxDeconvCombo);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.isnegative", &Ndata_bb_gem_m0_strip_isnegative, &b_Ndata_bb_gem_m0_strip_isnegative);
   fChain->SetBranchAddress("bb.gem.m0.strip.isnegative", bb_gem_m0_strip_isnegative, &b_bb_gem_m0_strip_isnegative);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.isnegativeU", &Ndata_bb_gem_m0_strip_isnegativeU, &b_Ndata_bb_gem_m0_strip_isnegativeU);
   fChain->SetBranchAddress("bb.gem.m0.strip.isnegativeU", bb_gem_m0_strip_isnegativeU, &b_bb_gem_m0_strip_isnegativeU);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.isnegativeV", &Ndata_bb_gem_m0_strip_isnegativeV, &b_Ndata_bb_gem_m0_strip_isnegativeV);
   fChain->SetBranchAddress("bb.gem.m0.strip.isnegativeV", bb_gem_m0_strip_isnegativeV, &b_bb_gem_m0_strip_isnegativeV);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.isnegontrack", &Ndata_bb_gem_m0_strip_isnegontrack, &b_Ndata_bb_gem_m0_strip_isnegontrack);
   fChain->SetBranchAddress("bb.gem.m0.strip.isnegontrack", bb_gem_m0_strip_isnegontrack, &b_bb_gem_m0_strip_isnegontrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.isnegontrackU", &Ndata_bb_gem_m0_strip_isnegontrackU, &b_Ndata_bb_gem_m0_strip_isnegontrackU);
   fChain->SetBranchAddress("bb.gem.m0.strip.isnegontrackU", bb_gem_m0_strip_isnegontrackU, &b_bb_gem_m0_strip_isnegontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.isnegontrackV", &Ndata_bb_gem_m0_strip_isnegontrackV, &b_Ndata_bb_gem_m0_strip_isnegontrackV);
   fChain->SetBranchAddress("bb.gem.m0.strip.isnegontrackV", bb_gem_m0_strip_isnegontrackV, &b_bb_gem_m0_strip_isnegontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.istrip", &Ndata_bb_gem_m0_strip_istrip, &b_Ndata_bb_gem_m0_strip_istrip);
   fChain->SetBranchAddress("bb.gem.m0.strip.istrip", bb_gem_m0_strip_istrip, &b_bb_gem_m0_strip_istrip);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.itrack", &Ndata_bb_gem_m0_strip_itrack, &b_Ndata_bb_gem_m0_strip_itrack);
   fChain->SetBranchAddress("bb.gem.m0.strip.itrack", bb_gem_m0_strip_itrack, &b_bb_gem_m0_strip_itrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.mpd", &Ndata_bb_gem_m0_strip_mpd, &b_Ndata_bb_gem_m0_strip_mpd);
   fChain->SetBranchAddress("bb.gem.m0.strip.mpd", bb_gem_m0_strip_mpd, &b_bb_gem_m0_strip_mpd);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.ontrack", &Ndata_bb_gem_m0_strip_ontrack, &b_Ndata_bb_gem_m0_strip_ontrack);
   fChain->SetBranchAddress("bb.gem.m0.strip.ontrack", bb_gem_m0_strip_ontrack, &b_bb_gem_m0_strip_ontrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.ontrackU", &Ndata_bb_gem_m0_strip_ontrackU, &b_Ndata_bb_gem_m0_strip_ontrackU);
   fChain->SetBranchAddress("bb.gem.m0.strip.ontrackU", bb_gem_m0_strip_ontrackU, &b_bb_gem_m0_strip_ontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.ontrackV", &Ndata_bb_gem_m0_strip_ontrackV, &b_Ndata_bb_gem_m0_strip_ontrackV);
   fChain->SetBranchAddress("bb.gem.m0.strip.ontrackV", bb_gem_m0_strip_ontrackV, &b_bb_gem_m0_strip_ontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m0.strip.rawADCsamples", &Ndata_bb_gem_m0_strip_rawADCsamples, &b_Ndata_bb_gem_m0_strip_rawADCsamples);
   fChain->SetBranchAddress("bb.gem.m0.strip.rawADCsamples", bb_gem_m0_strip_rawADCsamples, &b_bb_gem_m0_strip_rawADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.clust.clustu_adc", &Ndata_bb_gem_m1_clust_clustu_adc, &b_Ndata_bb_gem_m1_clust_clustu_adc);
   fChain->SetBranchAddress("bb.gem.m1.clust.clustu_adc", bb_gem_m1_clust_clustu_adc, &b_bb_gem_m1_clust_clustu_adc);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.clust.clustu_pos", &Ndata_bb_gem_m1_clust_clustu_pos, &b_Ndata_bb_gem_m1_clust_clustu_pos);
   fChain->SetBranchAddress("bb.gem.m1.clust.clustu_pos", bb_gem_m1_clust_clustu_pos, &b_bb_gem_m1_clust_clustu_pos);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.clust.clustu_strips", &Ndata_bb_gem_m1_clust_clustu_strips, &b_Ndata_bb_gem_m1_clust_clustu_strips);
   fChain->SetBranchAddress("bb.gem.m1.clust.clustu_strips", bb_gem_m1_clust_clustu_strips, &b_bb_gem_m1_clust_clustu_strips);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.clust.clustu_time", &Ndata_bb_gem_m1_clust_clustu_time, &b_Ndata_bb_gem_m1_clust_clustu_time);
   fChain->SetBranchAddress("bb.gem.m1.clust.clustu_time", bb_gem_m1_clust_clustu_time, &b_bb_gem_m1_clust_clustu_time);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.clust.clustv_adc", &Ndata_bb_gem_m1_clust_clustv_adc, &b_Ndata_bb_gem_m1_clust_clustv_adc);
   fChain->SetBranchAddress("bb.gem.m1.clust.clustv_adc", bb_gem_m1_clust_clustv_adc, &b_bb_gem_m1_clust_clustv_adc);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.clust.clustv_pos", &Ndata_bb_gem_m1_clust_clustv_pos, &b_Ndata_bb_gem_m1_clust_clustv_pos);
   fChain->SetBranchAddress("bb.gem.m1.clust.clustv_pos", bb_gem_m1_clust_clustv_pos, &b_bb_gem_m1_clust_clustv_pos);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.clust.clustv_strips", &Ndata_bb_gem_m1_clust_clustv_strips, &b_Ndata_bb_gem_m1_clust_clustv_strips);
   fChain->SetBranchAddress("bb.gem.m1.clust.clustv_strips", bb_gem_m1_clust_clustv_strips, &b_bb_gem_m1_clust_clustv_strips);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.clust.clustv_time", &Ndata_bb_gem_m1_clust_clustv_time, &b_Ndata_bb_gem_m1_clust_clustv_time);
   fChain->SetBranchAddress("bb.gem.m1.clust.clustv_time", bb_gem_m1_clust_clustv_time, &b_bb_gem_m1_clust_clustv_time);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.clust.isnegativeU", &Ndata_bb_gem_m1_clust_isnegativeU, &b_Ndata_bb_gem_m1_clust_isnegativeU);
   fChain->SetBranchAddress("bb.gem.m1.clust.isnegativeU", bb_gem_m1_clust_isnegativeU, &b_bb_gem_m1_clust_isnegativeU);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.clust.isnegativeV", &Ndata_bb_gem_m1_clust_isnegativeV, &b_Ndata_bb_gem_m1_clust_isnegativeV);
   fChain->SetBranchAddress("bb.gem.m1.clust.isnegativeV", bb_gem_m1_clust_isnegativeV, &b_bb_gem_m1_clust_isnegativeV);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.clust.isnegontrackU", &Ndata_bb_gem_m1_clust_isnegontrackU, &b_Ndata_bb_gem_m1_clust_isnegontrackU);
   fChain->SetBranchAddress("bb.gem.m1.clust.isnegontrackU", bb_gem_m1_clust_isnegontrackU, &b_bb_gem_m1_clust_isnegontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.clust.isnegontrackV", &Ndata_bb_gem_m1_clust_isnegontrackV, &b_Ndata_bb_gem_m1_clust_isnegontrackV);
   fChain->SetBranchAddress("bb.gem.m1.clust.isnegontrackV", bb_gem_m1_clust_isnegontrackV, &b_bb_gem_m1_clust_isnegontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.ADCavg", &Ndata_bb_gem_m1_strip_ADCavg, &b_Ndata_bb_gem_m1_strip_ADCavg);
   fChain->SetBranchAddress("bb.gem.m1.strip.ADCavg", bb_gem_m1_strip_ADCavg, &b_bb_gem_m1_strip_ADCavg);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.ADCmax", &Ndata_bb_gem_m1_strip_ADCmax, &b_Ndata_bb_gem_m1_strip_ADCmax);
   fChain->SetBranchAddress("bb.gem.m1.strip.ADCmax", bb_gem_m1_strip_ADCmax, &b_bb_gem_m1_strip_ADCmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.ADCsamples", &Ndata_bb_gem_m1_strip_ADCsamples, &b_Ndata_bb_gem_m1_strip_ADCsamples);
   fChain->SetBranchAddress("bb.gem.m1.strip.ADCsamples", bb_gem_m1_strip_ADCsamples, &b_bb_gem_m1_strip_ADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.ADCsum", &Ndata_bb_gem_m1_strip_ADCsum, &b_Ndata_bb_gem_m1_strip_ADCsum);
   fChain->SetBranchAddress("bb.gem.m1.strip.ADCsum", bb_gem_m1_strip_ADCsum, &b_bb_gem_m1_strip_ADCsum);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.BUILD_ALL_SAMPLES", &Ndata_bb_gem_m1_strip_BUILD_ALL_SAMPLES, &b_Ndata_bb_gem_m1_strip_BUILD_ALL_SAMPLES);
   fChain->SetBranchAddress("bb.gem.m1.strip.BUILD_ALL_SAMPLES", bb_gem_m1_strip_BUILD_ALL_SAMPLES, &b_bb_gem_m1_strip_BUILD_ALL_SAMPLES);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.CM_GOOD", &Ndata_bb_gem_m1_strip_CM_GOOD, &b_Ndata_bb_gem_m1_strip_CM_GOOD);
   fChain->SetBranchAddress("bb.gem.m1.strip.CM_GOOD", bb_gem_m1_strip_CM_GOOD, &b_bb_gem_m1_strip_CM_GOOD);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.CorrCoeff", &Ndata_bb_gem_m1_strip_CorrCoeff, &b_Ndata_bb_gem_m1_strip_CorrCoeff);
   fChain->SetBranchAddress("bb.gem.m1.strip.CorrCoeff", bb_gem_m1_strip_CorrCoeff, &b_bb_gem_m1_strip_CorrCoeff);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.DeconvADCmax", &Ndata_bb_gem_m1_strip_DeconvADCmax, &b_Ndata_bb_gem_m1_strip_DeconvADCmax);
   fChain->SetBranchAddress("bb.gem.m1.strip.DeconvADCmax", bb_gem_m1_strip_DeconvADCmax, &b_bb_gem_m1_strip_DeconvADCmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.DeconvADCmaxCombo", &Ndata_bb_gem_m1_strip_DeconvADCmaxCombo, &b_Ndata_bb_gem_m1_strip_DeconvADCmaxCombo);
   fChain->SetBranchAddress("bb.gem.m1.strip.DeconvADCmaxCombo", bb_gem_m1_strip_DeconvADCmaxCombo, &b_bb_gem_m1_strip_DeconvADCmaxCombo);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.DeconvADCsamples", &Ndata_bb_gem_m1_strip_DeconvADCsamples, &b_Ndata_bb_gem_m1_strip_DeconvADCsamples);
   fChain->SetBranchAddress("bb.gem.m1.strip.DeconvADCsamples", bb_gem_m1_strip_DeconvADCsamples, &b_bb_gem_m1_strip_DeconvADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.DeconvADCsum", &Ndata_bb_gem_m1_strip_DeconvADCsum, &b_Ndata_bb_gem_m1_strip_DeconvADCsum);
   fChain->SetBranchAddress("bb.gem.m1.strip.DeconvADCsum", bb_gem_m1_strip_DeconvADCsum, &b_bb_gem_m1_strip_DeconvADCsum);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.ENABLE_CM", &Ndata_bb_gem_m1_strip_ENABLE_CM, &b_Ndata_bb_gem_m1_strip_ENABLE_CM);
   fChain->SetBranchAddress("bb.gem.m1.strip.ENABLE_CM", bb_gem_m1_strip_ENABLE_CM, &b_bb_gem_m1_strip_ENABLE_CM);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.IsU", &Ndata_bb_gem_m1_strip_IsU, &b_Ndata_bb_gem_m1_strip_IsU);
   fChain->SetBranchAddress("bb.gem.m1.strip.IsU", bb_gem_m1_strip_IsU, &b_bb_gem_m1_strip_IsU);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.IsV", &Ndata_bb_gem_m1_strip_IsV, &b_Ndata_bb_gem_m1_strip_IsV);
   fChain->SetBranchAddress("bb.gem.m1.strip.IsV", bb_gem_m1_strip_IsV, &b_bb_gem_m1_strip_IsV);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.TSchi2", &Ndata_bb_gem_m1_strip_TSchi2, &b_Ndata_bb_gem_m1_strip_TSchi2);
   fChain->SetBranchAddress("bb.gem.m1.strip.TSchi2", bb_gem_m1_strip_TSchi2, &b_bb_gem_m1_strip_TSchi2);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.TSprob", &Ndata_bb_gem_m1_strip_TSprob, &b_Ndata_bb_gem_m1_strip_TSprob);
   fChain->SetBranchAddress("bb.gem.m1.strip.TSprob", bb_gem_m1_strip_TSprob, &b_bb_gem_m1_strip_TSprob);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.Tcorr", &Ndata_bb_gem_m1_strip_Tcorr, &b_Ndata_bb_gem_m1_strip_Tcorr);
   fChain->SetBranchAddress("bb.gem.m1.strip.Tcorr", bb_gem_m1_strip_Tcorr, &b_bb_gem_m1_strip_Tcorr);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.Tdiff", &Ndata_bb_gem_m1_strip_Tdiff, &b_Ndata_bb_gem_m1_strip_Tdiff);
   fChain->SetBranchAddress("bb.gem.m1.strip.Tdiff", bb_gem_m1_strip_Tdiff, &b_bb_gem_m1_strip_Tdiff);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.Tfit", &Ndata_bb_gem_m1_strip_Tfit, &b_Ndata_bb_gem_m1_strip_Tfit);
   fChain->SetBranchAddress("bb.gem.m1.strip.Tfit", bb_gem_m1_strip_Tfit, &b_bb_gem_m1_strip_Tfit);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.Tmean", &Ndata_bb_gem_m1_strip_Tmean, &b_Ndata_bb_gem_m1_strip_Tmean);
   fChain->SetBranchAddress("bb.gem.m1.strip.Tmean", bb_gem_m1_strip_Tmean, &b_bb_gem_m1_strip_Tmean);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.TmeanDeconv", &Ndata_bb_gem_m1_strip_TmeanDeconv, &b_Ndata_bb_gem_m1_strip_TmeanDeconv);
   fChain->SetBranchAddress("bb.gem.m1.strip.TmeanDeconv", bb_gem_m1_strip_TmeanDeconv, &b_bb_gem_m1_strip_TmeanDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.Tsigma", &Ndata_bb_gem_m1_strip_Tsigma, &b_Ndata_bb_gem_m1_strip_Tsigma);
   fChain->SetBranchAddress("bb.gem.m1.strip.Tsigma", bb_gem_m1_strip_Tsigma, &b_bb_gem_m1_strip_Tsigma);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.adc_id", &Ndata_bb_gem_m1_strip_adc_id, &b_Ndata_bb_gem_m1_strip_adc_id);
   fChain->SetBranchAddress("bb.gem.m1.strip.adc_id", bb_gem_m1_strip_adc_id, &b_bb_gem_m1_strip_adc_id);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.crate", &Ndata_bb_gem_m1_strip_crate, &b_Ndata_bb_gem_m1_strip_crate);
   fChain->SetBranchAddress("bb.gem.m1.strip.crate", bb_gem_m1_strip_crate, &b_bb_gem_m1_strip_crate);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.isampmax", &Ndata_bb_gem_m1_strip_isampmax, &b_Ndata_bb_gem_m1_strip_isampmax);
   fChain->SetBranchAddress("bb.gem.m1.strip.isampmax", bb_gem_m1_strip_isampmax, &b_bb_gem_m1_strip_isampmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.isampmaxDeconv", &Ndata_bb_gem_m1_strip_isampmaxDeconv, &b_Ndata_bb_gem_m1_strip_isampmaxDeconv);
   fChain->SetBranchAddress("bb.gem.m1.strip.isampmaxDeconv", bb_gem_m1_strip_isampmaxDeconv, &b_bb_gem_m1_strip_isampmaxDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.isampmaxDeconvCombo", &Ndata_bb_gem_m1_strip_isampmaxDeconvCombo, &b_Ndata_bb_gem_m1_strip_isampmaxDeconvCombo);
   fChain->SetBranchAddress("bb.gem.m1.strip.isampmaxDeconvCombo", bb_gem_m1_strip_isampmaxDeconvCombo, &b_bb_gem_m1_strip_isampmaxDeconvCombo);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.isnegative", &Ndata_bb_gem_m1_strip_isnegative, &b_Ndata_bb_gem_m1_strip_isnegative);
   fChain->SetBranchAddress("bb.gem.m1.strip.isnegative", bb_gem_m1_strip_isnegative, &b_bb_gem_m1_strip_isnegative);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.isnegativeU", &Ndata_bb_gem_m1_strip_isnegativeU, &b_Ndata_bb_gem_m1_strip_isnegativeU);
   fChain->SetBranchAddress("bb.gem.m1.strip.isnegativeU", bb_gem_m1_strip_isnegativeU, &b_bb_gem_m1_strip_isnegativeU);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.isnegativeV", &Ndata_bb_gem_m1_strip_isnegativeV, &b_Ndata_bb_gem_m1_strip_isnegativeV);
   fChain->SetBranchAddress("bb.gem.m1.strip.isnegativeV", bb_gem_m1_strip_isnegativeV, &b_bb_gem_m1_strip_isnegativeV);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.isnegontrack", &Ndata_bb_gem_m1_strip_isnegontrack, &b_Ndata_bb_gem_m1_strip_isnegontrack);
   fChain->SetBranchAddress("bb.gem.m1.strip.isnegontrack", bb_gem_m1_strip_isnegontrack, &b_bb_gem_m1_strip_isnegontrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.isnegontrackU", &Ndata_bb_gem_m1_strip_isnegontrackU, &b_Ndata_bb_gem_m1_strip_isnegontrackU);
   fChain->SetBranchAddress("bb.gem.m1.strip.isnegontrackU", bb_gem_m1_strip_isnegontrackU, &b_bb_gem_m1_strip_isnegontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.isnegontrackV", &Ndata_bb_gem_m1_strip_isnegontrackV, &b_Ndata_bb_gem_m1_strip_isnegontrackV);
   fChain->SetBranchAddress("bb.gem.m1.strip.isnegontrackV", bb_gem_m1_strip_isnegontrackV, &b_bb_gem_m1_strip_isnegontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.istrip", &Ndata_bb_gem_m1_strip_istrip, &b_Ndata_bb_gem_m1_strip_istrip);
   fChain->SetBranchAddress("bb.gem.m1.strip.istrip", bb_gem_m1_strip_istrip, &b_bb_gem_m1_strip_istrip);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.itrack", &Ndata_bb_gem_m1_strip_itrack, &b_Ndata_bb_gem_m1_strip_itrack);
   fChain->SetBranchAddress("bb.gem.m1.strip.itrack", bb_gem_m1_strip_itrack, &b_bb_gem_m1_strip_itrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.mpd", &Ndata_bb_gem_m1_strip_mpd, &b_Ndata_bb_gem_m1_strip_mpd);
   fChain->SetBranchAddress("bb.gem.m1.strip.mpd", bb_gem_m1_strip_mpd, &b_bb_gem_m1_strip_mpd);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.ontrack", &Ndata_bb_gem_m1_strip_ontrack, &b_Ndata_bb_gem_m1_strip_ontrack);
   fChain->SetBranchAddress("bb.gem.m1.strip.ontrack", bb_gem_m1_strip_ontrack, &b_bb_gem_m1_strip_ontrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.ontrackU", &Ndata_bb_gem_m1_strip_ontrackU, &b_Ndata_bb_gem_m1_strip_ontrackU);
   fChain->SetBranchAddress("bb.gem.m1.strip.ontrackU", bb_gem_m1_strip_ontrackU, &b_bb_gem_m1_strip_ontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.ontrackV", &Ndata_bb_gem_m1_strip_ontrackV, &b_Ndata_bb_gem_m1_strip_ontrackV);
   fChain->SetBranchAddress("bb.gem.m1.strip.ontrackV", bb_gem_m1_strip_ontrackV, &b_bb_gem_m1_strip_ontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m1.strip.rawADCsamples", &Ndata_bb_gem_m1_strip_rawADCsamples, &b_Ndata_bb_gem_m1_strip_rawADCsamples);
   fChain->SetBranchAddress("bb.gem.m1.strip.rawADCsamples", bb_gem_m1_strip_rawADCsamples, &b_bb_gem_m1_strip_rawADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.clust.clustu_adc", &Ndata_bb_gem_m2_clust_clustu_adc, &b_Ndata_bb_gem_m2_clust_clustu_adc);
   fChain->SetBranchAddress("bb.gem.m2.clust.clustu_adc", bb_gem_m2_clust_clustu_adc, &b_bb_gem_m2_clust_clustu_adc);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.clust.clustu_pos", &Ndata_bb_gem_m2_clust_clustu_pos, &b_Ndata_bb_gem_m2_clust_clustu_pos);
   fChain->SetBranchAddress("bb.gem.m2.clust.clustu_pos", bb_gem_m2_clust_clustu_pos, &b_bb_gem_m2_clust_clustu_pos);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.clust.clustu_strips", &Ndata_bb_gem_m2_clust_clustu_strips, &b_Ndata_bb_gem_m2_clust_clustu_strips);
   fChain->SetBranchAddress("bb.gem.m2.clust.clustu_strips", bb_gem_m2_clust_clustu_strips, &b_bb_gem_m2_clust_clustu_strips);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.clust.clustu_time", &Ndata_bb_gem_m2_clust_clustu_time, &b_Ndata_bb_gem_m2_clust_clustu_time);
   fChain->SetBranchAddress("bb.gem.m2.clust.clustu_time", bb_gem_m2_clust_clustu_time, &b_bb_gem_m2_clust_clustu_time);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.clust.clustv_adc", &Ndata_bb_gem_m2_clust_clustv_adc, &b_Ndata_bb_gem_m2_clust_clustv_adc);
   fChain->SetBranchAddress("bb.gem.m2.clust.clustv_adc", bb_gem_m2_clust_clustv_adc, &b_bb_gem_m2_clust_clustv_adc);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.clust.clustv_pos", &Ndata_bb_gem_m2_clust_clustv_pos, &b_Ndata_bb_gem_m2_clust_clustv_pos);
   fChain->SetBranchAddress("bb.gem.m2.clust.clustv_pos", bb_gem_m2_clust_clustv_pos, &b_bb_gem_m2_clust_clustv_pos);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.clust.clustv_strips", &Ndata_bb_gem_m2_clust_clustv_strips, &b_Ndata_bb_gem_m2_clust_clustv_strips);
   fChain->SetBranchAddress("bb.gem.m2.clust.clustv_strips", bb_gem_m2_clust_clustv_strips, &b_bb_gem_m2_clust_clustv_strips);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.clust.clustv_time", &Ndata_bb_gem_m2_clust_clustv_time, &b_Ndata_bb_gem_m2_clust_clustv_time);
   fChain->SetBranchAddress("bb.gem.m2.clust.clustv_time", bb_gem_m2_clust_clustv_time, &b_bb_gem_m2_clust_clustv_time);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.clust.isnegativeU", &Ndata_bb_gem_m2_clust_isnegativeU, &b_Ndata_bb_gem_m2_clust_isnegativeU);
   fChain->SetBranchAddress("bb.gem.m2.clust.isnegativeU", bb_gem_m2_clust_isnegativeU, &b_bb_gem_m2_clust_isnegativeU);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.clust.isnegativeV", &Ndata_bb_gem_m2_clust_isnegativeV, &b_Ndata_bb_gem_m2_clust_isnegativeV);
   fChain->SetBranchAddress("bb.gem.m2.clust.isnegativeV", bb_gem_m2_clust_isnegativeV, &b_bb_gem_m2_clust_isnegativeV);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.clust.isnegontrackU", &Ndata_bb_gem_m2_clust_isnegontrackU, &b_Ndata_bb_gem_m2_clust_isnegontrackU);
   fChain->SetBranchAddress("bb.gem.m2.clust.isnegontrackU", bb_gem_m2_clust_isnegontrackU, &b_bb_gem_m2_clust_isnegontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.clust.isnegontrackV", &Ndata_bb_gem_m2_clust_isnegontrackV, &b_Ndata_bb_gem_m2_clust_isnegontrackV);
   fChain->SetBranchAddress("bb.gem.m2.clust.isnegontrackV", bb_gem_m2_clust_isnegontrackV, &b_bb_gem_m2_clust_isnegontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.ADCavg", &Ndata_bb_gem_m2_strip_ADCavg, &b_Ndata_bb_gem_m2_strip_ADCavg);
   fChain->SetBranchAddress("bb.gem.m2.strip.ADCavg", bb_gem_m2_strip_ADCavg, &b_bb_gem_m2_strip_ADCavg);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.ADCmax", &Ndata_bb_gem_m2_strip_ADCmax, &b_Ndata_bb_gem_m2_strip_ADCmax);
   fChain->SetBranchAddress("bb.gem.m2.strip.ADCmax", bb_gem_m2_strip_ADCmax, &b_bb_gem_m2_strip_ADCmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.ADCsamples", &Ndata_bb_gem_m2_strip_ADCsamples, &b_Ndata_bb_gem_m2_strip_ADCsamples);
   fChain->SetBranchAddress("bb.gem.m2.strip.ADCsamples", bb_gem_m2_strip_ADCsamples, &b_bb_gem_m2_strip_ADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.ADCsum", &Ndata_bb_gem_m2_strip_ADCsum, &b_Ndata_bb_gem_m2_strip_ADCsum);
   fChain->SetBranchAddress("bb.gem.m2.strip.ADCsum", bb_gem_m2_strip_ADCsum, &b_bb_gem_m2_strip_ADCsum);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.BUILD_ALL_SAMPLES", &Ndata_bb_gem_m2_strip_BUILD_ALL_SAMPLES, &b_Ndata_bb_gem_m2_strip_BUILD_ALL_SAMPLES);
   fChain->SetBranchAddress("bb.gem.m2.strip.BUILD_ALL_SAMPLES", bb_gem_m2_strip_BUILD_ALL_SAMPLES, &b_bb_gem_m2_strip_BUILD_ALL_SAMPLES);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.CM_GOOD", &Ndata_bb_gem_m2_strip_CM_GOOD, &b_Ndata_bb_gem_m2_strip_CM_GOOD);
   fChain->SetBranchAddress("bb.gem.m2.strip.CM_GOOD", bb_gem_m2_strip_CM_GOOD, &b_bb_gem_m2_strip_CM_GOOD);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.CorrCoeff", &Ndata_bb_gem_m2_strip_CorrCoeff, &b_Ndata_bb_gem_m2_strip_CorrCoeff);
   fChain->SetBranchAddress("bb.gem.m2.strip.CorrCoeff", bb_gem_m2_strip_CorrCoeff, &b_bb_gem_m2_strip_CorrCoeff);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.DeconvADCmax", &Ndata_bb_gem_m2_strip_DeconvADCmax, &b_Ndata_bb_gem_m2_strip_DeconvADCmax);
   fChain->SetBranchAddress("bb.gem.m2.strip.DeconvADCmax", bb_gem_m2_strip_DeconvADCmax, &b_bb_gem_m2_strip_DeconvADCmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.DeconvADCmaxCombo", &Ndata_bb_gem_m2_strip_DeconvADCmaxCombo, &b_Ndata_bb_gem_m2_strip_DeconvADCmaxCombo);
   fChain->SetBranchAddress("bb.gem.m2.strip.DeconvADCmaxCombo", bb_gem_m2_strip_DeconvADCmaxCombo, &b_bb_gem_m2_strip_DeconvADCmaxCombo);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.DeconvADCsamples", &Ndata_bb_gem_m2_strip_DeconvADCsamples, &b_Ndata_bb_gem_m2_strip_DeconvADCsamples);
   fChain->SetBranchAddress("bb.gem.m2.strip.DeconvADCsamples", bb_gem_m2_strip_DeconvADCsamples, &b_bb_gem_m2_strip_DeconvADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.DeconvADCsum", &Ndata_bb_gem_m2_strip_DeconvADCsum, &b_Ndata_bb_gem_m2_strip_DeconvADCsum);
   fChain->SetBranchAddress("bb.gem.m2.strip.DeconvADCsum", bb_gem_m2_strip_DeconvADCsum, &b_bb_gem_m2_strip_DeconvADCsum);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.ENABLE_CM", &Ndata_bb_gem_m2_strip_ENABLE_CM, &b_Ndata_bb_gem_m2_strip_ENABLE_CM);
   fChain->SetBranchAddress("bb.gem.m2.strip.ENABLE_CM", bb_gem_m2_strip_ENABLE_CM, &b_bb_gem_m2_strip_ENABLE_CM);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.IsU", &Ndata_bb_gem_m2_strip_IsU, &b_Ndata_bb_gem_m2_strip_IsU);
   fChain->SetBranchAddress("bb.gem.m2.strip.IsU", bb_gem_m2_strip_IsU, &b_bb_gem_m2_strip_IsU);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.IsV", &Ndata_bb_gem_m2_strip_IsV, &b_Ndata_bb_gem_m2_strip_IsV);
   fChain->SetBranchAddress("bb.gem.m2.strip.IsV", bb_gem_m2_strip_IsV, &b_bb_gem_m2_strip_IsV);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.TSchi2", &Ndata_bb_gem_m2_strip_TSchi2, &b_Ndata_bb_gem_m2_strip_TSchi2);
   fChain->SetBranchAddress("bb.gem.m2.strip.TSchi2", bb_gem_m2_strip_TSchi2, &b_bb_gem_m2_strip_TSchi2);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.TSprob", &Ndata_bb_gem_m2_strip_TSprob, &b_Ndata_bb_gem_m2_strip_TSprob);
   fChain->SetBranchAddress("bb.gem.m2.strip.TSprob", bb_gem_m2_strip_TSprob, &b_bb_gem_m2_strip_TSprob);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.Tcorr", &Ndata_bb_gem_m2_strip_Tcorr, &b_Ndata_bb_gem_m2_strip_Tcorr);
   fChain->SetBranchAddress("bb.gem.m2.strip.Tcorr", bb_gem_m2_strip_Tcorr, &b_bb_gem_m2_strip_Tcorr);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.Tdiff", &Ndata_bb_gem_m2_strip_Tdiff, &b_Ndata_bb_gem_m2_strip_Tdiff);
   fChain->SetBranchAddress("bb.gem.m2.strip.Tdiff", bb_gem_m2_strip_Tdiff, &b_bb_gem_m2_strip_Tdiff);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.Tfit", &Ndata_bb_gem_m2_strip_Tfit, &b_Ndata_bb_gem_m2_strip_Tfit);
   fChain->SetBranchAddress("bb.gem.m2.strip.Tfit", bb_gem_m2_strip_Tfit, &b_bb_gem_m2_strip_Tfit);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.Tmean", &Ndata_bb_gem_m2_strip_Tmean, &b_Ndata_bb_gem_m2_strip_Tmean);
   fChain->SetBranchAddress("bb.gem.m2.strip.Tmean", bb_gem_m2_strip_Tmean, &b_bb_gem_m2_strip_Tmean);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.TmeanDeconv", &Ndata_bb_gem_m2_strip_TmeanDeconv, &b_Ndata_bb_gem_m2_strip_TmeanDeconv);
   fChain->SetBranchAddress("bb.gem.m2.strip.TmeanDeconv", bb_gem_m2_strip_TmeanDeconv, &b_bb_gem_m2_strip_TmeanDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.Tsigma", &Ndata_bb_gem_m2_strip_Tsigma, &b_Ndata_bb_gem_m2_strip_Tsigma);
   fChain->SetBranchAddress("bb.gem.m2.strip.Tsigma", bb_gem_m2_strip_Tsigma, &b_bb_gem_m2_strip_Tsigma);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.adc_id", &Ndata_bb_gem_m2_strip_adc_id, &b_Ndata_bb_gem_m2_strip_adc_id);
   fChain->SetBranchAddress("bb.gem.m2.strip.adc_id", bb_gem_m2_strip_adc_id, &b_bb_gem_m2_strip_adc_id);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.crate", &Ndata_bb_gem_m2_strip_crate, &b_Ndata_bb_gem_m2_strip_crate);
   fChain->SetBranchAddress("bb.gem.m2.strip.crate", bb_gem_m2_strip_crate, &b_bb_gem_m2_strip_crate);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.isampmax", &Ndata_bb_gem_m2_strip_isampmax, &b_Ndata_bb_gem_m2_strip_isampmax);
   fChain->SetBranchAddress("bb.gem.m2.strip.isampmax", bb_gem_m2_strip_isampmax, &b_bb_gem_m2_strip_isampmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.isampmaxDeconv", &Ndata_bb_gem_m2_strip_isampmaxDeconv, &b_Ndata_bb_gem_m2_strip_isampmaxDeconv);
   fChain->SetBranchAddress("bb.gem.m2.strip.isampmaxDeconv", bb_gem_m2_strip_isampmaxDeconv, &b_bb_gem_m2_strip_isampmaxDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.isampmaxDeconvCombo", &Ndata_bb_gem_m2_strip_isampmaxDeconvCombo, &b_Ndata_bb_gem_m2_strip_isampmaxDeconvCombo);
   fChain->SetBranchAddress("bb.gem.m2.strip.isampmaxDeconvCombo", bb_gem_m2_strip_isampmaxDeconvCombo, &b_bb_gem_m2_strip_isampmaxDeconvCombo);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.isnegative", &Ndata_bb_gem_m2_strip_isnegative, &b_Ndata_bb_gem_m2_strip_isnegative);
   fChain->SetBranchAddress("bb.gem.m2.strip.isnegative", bb_gem_m2_strip_isnegative, &b_bb_gem_m2_strip_isnegative);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.isnegativeU", &Ndata_bb_gem_m2_strip_isnegativeU, &b_Ndata_bb_gem_m2_strip_isnegativeU);
   fChain->SetBranchAddress("bb.gem.m2.strip.isnegativeU", bb_gem_m2_strip_isnegativeU, &b_bb_gem_m2_strip_isnegativeU);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.isnegativeV", &Ndata_bb_gem_m2_strip_isnegativeV, &b_Ndata_bb_gem_m2_strip_isnegativeV);
   fChain->SetBranchAddress("bb.gem.m2.strip.isnegativeV", bb_gem_m2_strip_isnegativeV, &b_bb_gem_m2_strip_isnegativeV);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.isnegontrack", &Ndata_bb_gem_m2_strip_isnegontrack, &b_Ndata_bb_gem_m2_strip_isnegontrack);
   fChain->SetBranchAddress("bb.gem.m2.strip.isnegontrack", bb_gem_m2_strip_isnegontrack, &b_bb_gem_m2_strip_isnegontrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.isnegontrackU", &Ndata_bb_gem_m2_strip_isnegontrackU, &b_Ndata_bb_gem_m2_strip_isnegontrackU);
   fChain->SetBranchAddress("bb.gem.m2.strip.isnegontrackU", bb_gem_m2_strip_isnegontrackU, &b_bb_gem_m2_strip_isnegontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.isnegontrackV", &Ndata_bb_gem_m2_strip_isnegontrackV, &b_Ndata_bb_gem_m2_strip_isnegontrackV);
   fChain->SetBranchAddress("bb.gem.m2.strip.isnegontrackV", bb_gem_m2_strip_isnegontrackV, &b_bb_gem_m2_strip_isnegontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.istrip", &Ndata_bb_gem_m2_strip_istrip, &b_Ndata_bb_gem_m2_strip_istrip);
   fChain->SetBranchAddress("bb.gem.m2.strip.istrip", bb_gem_m2_strip_istrip, &b_bb_gem_m2_strip_istrip);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.itrack", &Ndata_bb_gem_m2_strip_itrack, &b_Ndata_bb_gem_m2_strip_itrack);
   fChain->SetBranchAddress("bb.gem.m2.strip.itrack", bb_gem_m2_strip_itrack, &b_bb_gem_m2_strip_itrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.mpd", &Ndata_bb_gem_m2_strip_mpd, &b_Ndata_bb_gem_m2_strip_mpd);
   fChain->SetBranchAddress("bb.gem.m2.strip.mpd", bb_gem_m2_strip_mpd, &b_bb_gem_m2_strip_mpd);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.ontrack", &Ndata_bb_gem_m2_strip_ontrack, &b_Ndata_bb_gem_m2_strip_ontrack);
   fChain->SetBranchAddress("bb.gem.m2.strip.ontrack", bb_gem_m2_strip_ontrack, &b_bb_gem_m2_strip_ontrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.ontrackU", &Ndata_bb_gem_m2_strip_ontrackU, &b_Ndata_bb_gem_m2_strip_ontrackU);
   fChain->SetBranchAddress("bb.gem.m2.strip.ontrackU", bb_gem_m2_strip_ontrackU, &b_bb_gem_m2_strip_ontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.ontrackV", &Ndata_bb_gem_m2_strip_ontrackV, &b_Ndata_bb_gem_m2_strip_ontrackV);
   fChain->SetBranchAddress("bb.gem.m2.strip.ontrackV", bb_gem_m2_strip_ontrackV, &b_bb_gem_m2_strip_ontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m2.strip.rawADCsamples", &Ndata_bb_gem_m2_strip_rawADCsamples, &b_Ndata_bb_gem_m2_strip_rawADCsamples);
   fChain->SetBranchAddress("bb.gem.m2.strip.rawADCsamples", bb_gem_m2_strip_rawADCsamples, &b_bb_gem_m2_strip_rawADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.clust.clustu_adc", &Ndata_bb_gem_m3_clust_clustu_adc, &b_Ndata_bb_gem_m3_clust_clustu_adc);
   fChain->SetBranchAddress("bb.gem.m3.clust.clustu_adc", bb_gem_m3_clust_clustu_adc, &b_bb_gem_m3_clust_clustu_adc);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.clust.clustu_pos", &Ndata_bb_gem_m3_clust_clustu_pos, &b_Ndata_bb_gem_m3_clust_clustu_pos);
   fChain->SetBranchAddress("bb.gem.m3.clust.clustu_pos", bb_gem_m3_clust_clustu_pos, &b_bb_gem_m3_clust_clustu_pos);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.clust.clustu_strips", &Ndata_bb_gem_m3_clust_clustu_strips, &b_Ndata_bb_gem_m3_clust_clustu_strips);
   fChain->SetBranchAddress("bb.gem.m3.clust.clustu_strips", bb_gem_m3_clust_clustu_strips, &b_bb_gem_m3_clust_clustu_strips);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.clust.clustu_time", &Ndata_bb_gem_m3_clust_clustu_time, &b_Ndata_bb_gem_m3_clust_clustu_time);
   fChain->SetBranchAddress("bb.gem.m3.clust.clustu_time", bb_gem_m3_clust_clustu_time, &b_bb_gem_m3_clust_clustu_time);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.clust.clustv_adc", &Ndata_bb_gem_m3_clust_clustv_adc, &b_Ndata_bb_gem_m3_clust_clustv_adc);
   fChain->SetBranchAddress("bb.gem.m3.clust.clustv_adc", bb_gem_m3_clust_clustv_adc, &b_bb_gem_m3_clust_clustv_adc);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.clust.clustv_pos", &Ndata_bb_gem_m3_clust_clustv_pos, &b_Ndata_bb_gem_m3_clust_clustv_pos);
   fChain->SetBranchAddress("bb.gem.m3.clust.clustv_pos", bb_gem_m3_clust_clustv_pos, &b_bb_gem_m3_clust_clustv_pos);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.clust.clustv_strips", &Ndata_bb_gem_m3_clust_clustv_strips, &b_Ndata_bb_gem_m3_clust_clustv_strips);
   fChain->SetBranchAddress("bb.gem.m3.clust.clustv_strips", bb_gem_m3_clust_clustv_strips, &b_bb_gem_m3_clust_clustv_strips);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.clust.clustv_time", &Ndata_bb_gem_m3_clust_clustv_time, &b_Ndata_bb_gem_m3_clust_clustv_time);
   fChain->SetBranchAddress("bb.gem.m3.clust.clustv_time", bb_gem_m3_clust_clustv_time, &b_bb_gem_m3_clust_clustv_time);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.clust.isnegativeU", &Ndata_bb_gem_m3_clust_isnegativeU, &b_Ndata_bb_gem_m3_clust_isnegativeU);
   fChain->SetBranchAddress("bb.gem.m3.clust.isnegativeU", bb_gem_m3_clust_isnegativeU, &b_bb_gem_m3_clust_isnegativeU);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.clust.isnegativeV", &Ndata_bb_gem_m3_clust_isnegativeV, &b_Ndata_bb_gem_m3_clust_isnegativeV);
   fChain->SetBranchAddress("bb.gem.m3.clust.isnegativeV", bb_gem_m3_clust_isnegativeV, &b_bb_gem_m3_clust_isnegativeV);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.clust.isnegontrackU", &Ndata_bb_gem_m3_clust_isnegontrackU, &b_Ndata_bb_gem_m3_clust_isnegontrackU);
   fChain->SetBranchAddress("bb.gem.m3.clust.isnegontrackU", bb_gem_m3_clust_isnegontrackU, &b_bb_gem_m3_clust_isnegontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.clust.isnegontrackV", &Ndata_bb_gem_m3_clust_isnegontrackV, &b_Ndata_bb_gem_m3_clust_isnegontrackV);
   fChain->SetBranchAddress("bb.gem.m3.clust.isnegontrackV", bb_gem_m3_clust_isnegontrackV, &b_bb_gem_m3_clust_isnegontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.ADCavg", &Ndata_bb_gem_m3_strip_ADCavg, &b_Ndata_bb_gem_m3_strip_ADCavg);
   fChain->SetBranchAddress("bb.gem.m3.strip.ADCavg", bb_gem_m3_strip_ADCavg, &b_bb_gem_m3_strip_ADCavg);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.ADCmax", &Ndata_bb_gem_m3_strip_ADCmax, &b_Ndata_bb_gem_m3_strip_ADCmax);
   fChain->SetBranchAddress("bb.gem.m3.strip.ADCmax", bb_gem_m3_strip_ADCmax, &b_bb_gem_m3_strip_ADCmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.ADCsamples", &Ndata_bb_gem_m3_strip_ADCsamples, &b_Ndata_bb_gem_m3_strip_ADCsamples);
   fChain->SetBranchAddress("bb.gem.m3.strip.ADCsamples", bb_gem_m3_strip_ADCsamples, &b_bb_gem_m3_strip_ADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.ADCsum", &Ndata_bb_gem_m3_strip_ADCsum, &b_Ndata_bb_gem_m3_strip_ADCsum);
   fChain->SetBranchAddress("bb.gem.m3.strip.ADCsum", bb_gem_m3_strip_ADCsum, &b_bb_gem_m3_strip_ADCsum);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.BUILD_ALL_SAMPLES", &Ndata_bb_gem_m3_strip_BUILD_ALL_SAMPLES, &b_Ndata_bb_gem_m3_strip_BUILD_ALL_SAMPLES);
   fChain->SetBranchAddress("bb.gem.m3.strip.BUILD_ALL_SAMPLES", bb_gem_m3_strip_BUILD_ALL_SAMPLES, &b_bb_gem_m3_strip_BUILD_ALL_SAMPLES);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.CM_GOOD", &Ndata_bb_gem_m3_strip_CM_GOOD, &b_Ndata_bb_gem_m3_strip_CM_GOOD);
   fChain->SetBranchAddress("bb.gem.m3.strip.CM_GOOD", bb_gem_m3_strip_CM_GOOD, &b_bb_gem_m3_strip_CM_GOOD);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.CorrCoeff", &Ndata_bb_gem_m3_strip_CorrCoeff, &b_Ndata_bb_gem_m3_strip_CorrCoeff);
   fChain->SetBranchAddress("bb.gem.m3.strip.CorrCoeff", bb_gem_m3_strip_CorrCoeff, &b_bb_gem_m3_strip_CorrCoeff);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.DeconvADCmax", &Ndata_bb_gem_m3_strip_DeconvADCmax, &b_Ndata_bb_gem_m3_strip_DeconvADCmax);
   fChain->SetBranchAddress("bb.gem.m3.strip.DeconvADCmax", bb_gem_m3_strip_DeconvADCmax, &b_bb_gem_m3_strip_DeconvADCmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.DeconvADCmaxCombo", &Ndata_bb_gem_m3_strip_DeconvADCmaxCombo, &b_Ndata_bb_gem_m3_strip_DeconvADCmaxCombo);
   fChain->SetBranchAddress("bb.gem.m3.strip.DeconvADCmaxCombo", bb_gem_m3_strip_DeconvADCmaxCombo, &b_bb_gem_m3_strip_DeconvADCmaxCombo);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.DeconvADCsamples", &Ndata_bb_gem_m3_strip_DeconvADCsamples, &b_Ndata_bb_gem_m3_strip_DeconvADCsamples);
   fChain->SetBranchAddress("bb.gem.m3.strip.DeconvADCsamples", bb_gem_m3_strip_DeconvADCsamples, &b_bb_gem_m3_strip_DeconvADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.DeconvADCsum", &Ndata_bb_gem_m3_strip_DeconvADCsum, &b_Ndata_bb_gem_m3_strip_DeconvADCsum);
   fChain->SetBranchAddress("bb.gem.m3.strip.DeconvADCsum", bb_gem_m3_strip_DeconvADCsum, &b_bb_gem_m3_strip_DeconvADCsum);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.ENABLE_CM", &Ndata_bb_gem_m3_strip_ENABLE_CM, &b_Ndata_bb_gem_m3_strip_ENABLE_CM);
   fChain->SetBranchAddress("bb.gem.m3.strip.ENABLE_CM", bb_gem_m3_strip_ENABLE_CM, &b_bb_gem_m3_strip_ENABLE_CM);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.IsU", &Ndata_bb_gem_m3_strip_IsU, &b_Ndata_bb_gem_m3_strip_IsU);
   fChain->SetBranchAddress("bb.gem.m3.strip.IsU", bb_gem_m3_strip_IsU, &b_bb_gem_m3_strip_IsU);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.IsV", &Ndata_bb_gem_m3_strip_IsV, &b_Ndata_bb_gem_m3_strip_IsV);
   fChain->SetBranchAddress("bb.gem.m3.strip.IsV", bb_gem_m3_strip_IsV, &b_bb_gem_m3_strip_IsV);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.TSchi2", &Ndata_bb_gem_m3_strip_TSchi2, &b_Ndata_bb_gem_m3_strip_TSchi2);
   fChain->SetBranchAddress("bb.gem.m3.strip.TSchi2", bb_gem_m3_strip_TSchi2, &b_bb_gem_m3_strip_TSchi2);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.TSprob", &Ndata_bb_gem_m3_strip_TSprob, &b_Ndata_bb_gem_m3_strip_TSprob);
   fChain->SetBranchAddress("bb.gem.m3.strip.TSprob", bb_gem_m3_strip_TSprob, &b_bb_gem_m3_strip_TSprob);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.Tcorr", &Ndata_bb_gem_m3_strip_Tcorr, &b_Ndata_bb_gem_m3_strip_Tcorr);
   fChain->SetBranchAddress("bb.gem.m3.strip.Tcorr", bb_gem_m3_strip_Tcorr, &b_bb_gem_m3_strip_Tcorr);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.Tdiff", &Ndata_bb_gem_m3_strip_Tdiff, &b_Ndata_bb_gem_m3_strip_Tdiff);
   fChain->SetBranchAddress("bb.gem.m3.strip.Tdiff", bb_gem_m3_strip_Tdiff, &b_bb_gem_m3_strip_Tdiff);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.Tfit", &Ndata_bb_gem_m3_strip_Tfit, &b_Ndata_bb_gem_m3_strip_Tfit);
   fChain->SetBranchAddress("bb.gem.m3.strip.Tfit", bb_gem_m3_strip_Tfit, &b_bb_gem_m3_strip_Tfit);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.Tmean", &Ndata_bb_gem_m3_strip_Tmean, &b_Ndata_bb_gem_m3_strip_Tmean);
   fChain->SetBranchAddress("bb.gem.m3.strip.Tmean", bb_gem_m3_strip_Tmean, &b_bb_gem_m3_strip_Tmean);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.TmeanDeconv", &Ndata_bb_gem_m3_strip_TmeanDeconv, &b_Ndata_bb_gem_m3_strip_TmeanDeconv);
   fChain->SetBranchAddress("bb.gem.m3.strip.TmeanDeconv", bb_gem_m3_strip_TmeanDeconv, &b_bb_gem_m3_strip_TmeanDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.Tsigma", &Ndata_bb_gem_m3_strip_Tsigma, &b_Ndata_bb_gem_m3_strip_Tsigma);
   fChain->SetBranchAddress("bb.gem.m3.strip.Tsigma", bb_gem_m3_strip_Tsigma, &b_bb_gem_m3_strip_Tsigma);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.adc_id", &Ndata_bb_gem_m3_strip_adc_id, &b_Ndata_bb_gem_m3_strip_adc_id);
   fChain->SetBranchAddress("bb.gem.m3.strip.adc_id", bb_gem_m3_strip_adc_id, &b_bb_gem_m3_strip_adc_id);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.crate", &Ndata_bb_gem_m3_strip_crate, &b_Ndata_bb_gem_m3_strip_crate);
   fChain->SetBranchAddress("bb.gem.m3.strip.crate", bb_gem_m3_strip_crate, &b_bb_gem_m3_strip_crate);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.isampmax", &Ndata_bb_gem_m3_strip_isampmax, &b_Ndata_bb_gem_m3_strip_isampmax);
   fChain->SetBranchAddress("bb.gem.m3.strip.isampmax", bb_gem_m3_strip_isampmax, &b_bb_gem_m3_strip_isampmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.isampmaxDeconv", &Ndata_bb_gem_m3_strip_isampmaxDeconv, &b_Ndata_bb_gem_m3_strip_isampmaxDeconv);
   fChain->SetBranchAddress("bb.gem.m3.strip.isampmaxDeconv", bb_gem_m3_strip_isampmaxDeconv, &b_bb_gem_m3_strip_isampmaxDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.isampmaxDeconvCombo", &Ndata_bb_gem_m3_strip_isampmaxDeconvCombo, &b_Ndata_bb_gem_m3_strip_isampmaxDeconvCombo);
   fChain->SetBranchAddress("bb.gem.m3.strip.isampmaxDeconvCombo", bb_gem_m3_strip_isampmaxDeconvCombo, &b_bb_gem_m3_strip_isampmaxDeconvCombo);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.isnegative", &Ndata_bb_gem_m3_strip_isnegative, &b_Ndata_bb_gem_m3_strip_isnegative);
   fChain->SetBranchAddress("bb.gem.m3.strip.isnegative", bb_gem_m3_strip_isnegative, &b_bb_gem_m3_strip_isnegative);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.isnegativeU", &Ndata_bb_gem_m3_strip_isnegativeU, &b_Ndata_bb_gem_m3_strip_isnegativeU);
   fChain->SetBranchAddress("bb.gem.m3.strip.isnegativeU", bb_gem_m3_strip_isnegativeU, &b_bb_gem_m3_strip_isnegativeU);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.isnegativeV", &Ndata_bb_gem_m3_strip_isnegativeV, &b_Ndata_bb_gem_m3_strip_isnegativeV);
   fChain->SetBranchAddress("bb.gem.m3.strip.isnegativeV", bb_gem_m3_strip_isnegativeV, &b_bb_gem_m3_strip_isnegativeV);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.isnegontrack", &Ndata_bb_gem_m3_strip_isnegontrack, &b_Ndata_bb_gem_m3_strip_isnegontrack);
   fChain->SetBranchAddress("bb.gem.m3.strip.isnegontrack", bb_gem_m3_strip_isnegontrack, &b_bb_gem_m3_strip_isnegontrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.isnegontrackU", &Ndata_bb_gem_m3_strip_isnegontrackU, &b_Ndata_bb_gem_m3_strip_isnegontrackU);
   fChain->SetBranchAddress("bb.gem.m3.strip.isnegontrackU", bb_gem_m3_strip_isnegontrackU, &b_bb_gem_m3_strip_isnegontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.isnegontrackV", &Ndata_bb_gem_m3_strip_isnegontrackV, &b_Ndata_bb_gem_m3_strip_isnegontrackV);
   fChain->SetBranchAddress("bb.gem.m3.strip.isnegontrackV", bb_gem_m3_strip_isnegontrackV, &b_bb_gem_m3_strip_isnegontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.istrip", &Ndata_bb_gem_m3_strip_istrip, &b_Ndata_bb_gem_m3_strip_istrip);
   fChain->SetBranchAddress("bb.gem.m3.strip.istrip", bb_gem_m3_strip_istrip, &b_bb_gem_m3_strip_istrip);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.itrack", &Ndata_bb_gem_m3_strip_itrack, &b_Ndata_bb_gem_m3_strip_itrack);
   fChain->SetBranchAddress("bb.gem.m3.strip.itrack", bb_gem_m3_strip_itrack, &b_bb_gem_m3_strip_itrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.mpd", &Ndata_bb_gem_m3_strip_mpd, &b_Ndata_bb_gem_m3_strip_mpd);
   fChain->SetBranchAddress("bb.gem.m3.strip.mpd", bb_gem_m3_strip_mpd, &b_bb_gem_m3_strip_mpd);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.ontrack", &Ndata_bb_gem_m3_strip_ontrack, &b_Ndata_bb_gem_m3_strip_ontrack);
   fChain->SetBranchAddress("bb.gem.m3.strip.ontrack", bb_gem_m3_strip_ontrack, &b_bb_gem_m3_strip_ontrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.ontrackU", &Ndata_bb_gem_m3_strip_ontrackU, &b_Ndata_bb_gem_m3_strip_ontrackU);
   fChain->SetBranchAddress("bb.gem.m3.strip.ontrackU", bb_gem_m3_strip_ontrackU, &b_bb_gem_m3_strip_ontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.ontrackV", &Ndata_bb_gem_m3_strip_ontrackV, &b_Ndata_bb_gem_m3_strip_ontrackV);
   fChain->SetBranchAddress("bb.gem.m3.strip.ontrackV", bb_gem_m3_strip_ontrackV, &b_bb_gem_m3_strip_ontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m3.strip.rawADCsamples", &Ndata_bb_gem_m3_strip_rawADCsamples, &b_Ndata_bb_gem_m3_strip_rawADCsamples);
   fChain->SetBranchAddress("bb.gem.m3.strip.rawADCsamples", bb_gem_m3_strip_rawADCsamples, &b_bb_gem_m3_strip_rawADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.clust.clustu_adc", &Ndata_bb_gem_m4_clust_clustu_adc, &b_Ndata_bb_gem_m4_clust_clustu_adc);
   fChain->SetBranchAddress("bb.gem.m4.clust.clustu_adc", bb_gem_m4_clust_clustu_adc, &b_bb_gem_m4_clust_clustu_adc);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.clust.clustu_pos", &Ndata_bb_gem_m4_clust_clustu_pos, &b_Ndata_bb_gem_m4_clust_clustu_pos);
   fChain->SetBranchAddress("bb.gem.m4.clust.clustu_pos", bb_gem_m4_clust_clustu_pos, &b_bb_gem_m4_clust_clustu_pos);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.clust.clustu_strips", &Ndata_bb_gem_m4_clust_clustu_strips, &b_Ndata_bb_gem_m4_clust_clustu_strips);
   fChain->SetBranchAddress("bb.gem.m4.clust.clustu_strips", bb_gem_m4_clust_clustu_strips, &b_bb_gem_m4_clust_clustu_strips);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.clust.clustu_time", &Ndata_bb_gem_m4_clust_clustu_time, &b_Ndata_bb_gem_m4_clust_clustu_time);
   fChain->SetBranchAddress("bb.gem.m4.clust.clustu_time", bb_gem_m4_clust_clustu_time, &b_bb_gem_m4_clust_clustu_time);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.clust.clustv_adc", &Ndata_bb_gem_m4_clust_clustv_adc, &b_Ndata_bb_gem_m4_clust_clustv_adc);
   fChain->SetBranchAddress("bb.gem.m4.clust.clustv_adc", bb_gem_m4_clust_clustv_adc, &b_bb_gem_m4_clust_clustv_adc);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.clust.clustv_pos", &Ndata_bb_gem_m4_clust_clustv_pos, &b_Ndata_bb_gem_m4_clust_clustv_pos);
   fChain->SetBranchAddress("bb.gem.m4.clust.clustv_pos", bb_gem_m4_clust_clustv_pos, &b_bb_gem_m4_clust_clustv_pos);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.clust.clustv_strips", &Ndata_bb_gem_m4_clust_clustv_strips, &b_Ndata_bb_gem_m4_clust_clustv_strips);
   fChain->SetBranchAddress("bb.gem.m4.clust.clustv_strips", bb_gem_m4_clust_clustv_strips, &b_bb_gem_m4_clust_clustv_strips);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.clust.clustv_time", &Ndata_bb_gem_m4_clust_clustv_time, &b_Ndata_bb_gem_m4_clust_clustv_time);
   fChain->SetBranchAddress("bb.gem.m4.clust.clustv_time", bb_gem_m4_clust_clustv_time, &b_bb_gem_m4_clust_clustv_time);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.clust.isnegativeU", &Ndata_bb_gem_m4_clust_isnegativeU, &b_Ndata_bb_gem_m4_clust_isnegativeU);
   fChain->SetBranchAddress("bb.gem.m4.clust.isnegativeU", bb_gem_m4_clust_isnegativeU, &b_bb_gem_m4_clust_isnegativeU);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.clust.isnegativeV", &Ndata_bb_gem_m4_clust_isnegativeV, &b_Ndata_bb_gem_m4_clust_isnegativeV);
   fChain->SetBranchAddress("bb.gem.m4.clust.isnegativeV", bb_gem_m4_clust_isnegativeV, &b_bb_gem_m4_clust_isnegativeV);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.clust.isnegontrackU", &Ndata_bb_gem_m4_clust_isnegontrackU, &b_Ndata_bb_gem_m4_clust_isnegontrackU);
   fChain->SetBranchAddress("bb.gem.m4.clust.isnegontrackU", bb_gem_m4_clust_isnegontrackU, &b_bb_gem_m4_clust_isnegontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.clust.isnegontrackV", &Ndata_bb_gem_m4_clust_isnegontrackV, &b_Ndata_bb_gem_m4_clust_isnegontrackV);
   fChain->SetBranchAddress("bb.gem.m4.clust.isnegontrackV", bb_gem_m4_clust_isnegontrackV, &b_bb_gem_m4_clust_isnegontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.ADCavg", &Ndata_bb_gem_m4_strip_ADCavg, &b_Ndata_bb_gem_m4_strip_ADCavg);
   fChain->SetBranchAddress("bb.gem.m4.strip.ADCavg", bb_gem_m4_strip_ADCavg, &b_bb_gem_m4_strip_ADCavg);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.ADCmax", &Ndata_bb_gem_m4_strip_ADCmax, &b_Ndata_bb_gem_m4_strip_ADCmax);
   fChain->SetBranchAddress("bb.gem.m4.strip.ADCmax", bb_gem_m4_strip_ADCmax, &b_bb_gem_m4_strip_ADCmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.ADCsamples", &Ndata_bb_gem_m4_strip_ADCsamples, &b_Ndata_bb_gem_m4_strip_ADCsamples);
   fChain->SetBranchAddress("bb.gem.m4.strip.ADCsamples", bb_gem_m4_strip_ADCsamples, &b_bb_gem_m4_strip_ADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.ADCsum", &Ndata_bb_gem_m4_strip_ADCsum, &b_Ndata_bb_gem_m4_strip_ADCsum);
   fChain->SetBranchAddress("bb.gem.m4.strip.ADCsum", bb_gem_m4_strip_ADCsum, &b_bb_gem_m4_strip_ADCsum);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.BUILD_ALL_SAMPLES", &Ndata_bb_gem_m4_strip_BUILD_ALL_SAMPLES, &b_Ndata_bb_gem_m4_strip_BUILD_ALL_SAMPLES);
   fChain->SetBranchAddress("bb.gem.m4.strip.BUILD_ALL_SAMPLES", bb_gem_m4_strip_BUILD_ALL_SAMPLES, &b_bb_gem_m4_strip_BUILD_ALL_SAMPLES);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.CM_GOOD", &Ndata_bb_gem_m4_strip_CM_GOOD, &b_Ndata_bb_gem_m4_strip_CM_GOOD);
   fChain->SetBranchAddress("bb.gem.m4.strip.CM_GOOD", bb_gem_m4_strip_CM_GOOD, &b_bb_gem_m4_strip_CM_GOOD);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.CorrCoeff", &Ndata_bb_gem_m4_strip_CorrCoeff, &b_Ndata_bb_gem_m4_strip_CorrCoeff);
   fChain->SetBranchAddress("bb.gem.m4.strip.CorrCoeff", bb_gem_m4_strip_CorrCoeff, &b_bb_gem_m4_strip_CorrCoeff);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.DeconvADCmax", &Ndata_bb_gem_m4_strip_DeconvADCmax, &b_Ndata_bb_gem_m4_strip_DeconvADCmax);
   fChain->SetBranchAddress("bb.gem.m4.strip.DeconvADCmax", bb_gem_m4_strip_DeconvADCmax, &b_bb_gem_m4_strip_DeconvADCmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.DeconvADCmaxCombo", &Ndata_bb_gem_m4_strip_DeconvADCmaxCombo, &b_Ndata_bb_gem_m4_strip_DeconvADCmaxCombo);
   fChain->SetBranchAddress("bb.gem.m4.strip.DeconvADCmaxCombo", bb_gem_m4_strip_DeconvADCmaxCombo, &b_bb_gem_m4_strip_DeconvADCmaxCombo);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.DeconvADCsamples", &Ndata_bb_gem_m4_strip_DeconvADCsamples, &b_Ndata_bb_gem_m4_strip_DeconvADCsamples);
   fChain->SetBranchAddress("bb.gem.m4.strip.DeconvADCsamples", bb_gem_m4_strip_DeconvADCsamples, &b_bb_gem_m4_strip_DeconvADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.DeconvADCsum", &Ndata_bb_gem_m4_strip_DeconvADCsum, &b_Ndata_bb_gem_m4_strip_DeconvADCsum);
   fChain->SetBranchAddress("bb.gem.m4.strip.DeconvADCsum", bb_gem_m4_strip_DeconvADCsum, &b_bb_gem_m4_strip_DeconvADCsum);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.ENABLE_CM", &Ndata_bb_gem_m4_strip_ENABLE_CM, &b_Ndata_bb_gem_m4_strip_ENABLE_CM);
   fChain->SetBranchAddress("bb.gem.m4.strip.ENABLE_CM", bb_gem_m4_strip_ENABLE_CM, &b_bb_gem_m4_strip_ENABLE_CM);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.IsU", &Ndata_bb_gem_m4_strip_IsU, &b_Ndata_bb_gem_m4_strip_IsU);
   fChain->SetBranchAddress("bb.gem.m4.strip.IsU", bb_gem_m4_strip_IsU, &b_bb_gem_m4_strip_IsU);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.IsV", &Ndata_bb_gem_m4_strip_IsV, &b_Ndata_bb_gem_m4_strip_IsV);
   fChain->SetBranchAddress("bb.gem.m4.strip.IsV", bb_gem_m4_strip_IsV, &b_bb_gem_m4_strip_IsV);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.TSchi2", &Ndata_bb_gem_m4_strip_TSchi2, &b_Ndata_bb_gem_m4_strip_TSchi2);
   fChain->SetBranchAddress("bb.gem.m4.strip.TSchi2", bb_gem_m4_strip_TSchi2, &b_bb_gem_m4_strip_TSchi2);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.TSprob", &Ndata_bb_gem_m4_strip_TSprob, &b_Ndata_bb_gem_m4_strip_TSprob);
   fChain->SetBranchAddress("bb.gem.m4.strip.TSprob", bb_gem_m4_strip_TSprob, &b_bb_gem_m4_strip_TSprob);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.Tcorr", &Ndata_bb_gem_m4_strip_Tcorr, &b_Ndata_bb_gem_m4_strip_Tcorr);
   fChain->SetBranchAddress("bb.gem.m4.strip.Tcorr", bb_gem_m4_strip_Tcorr, &b_bb_gem_m4_strip_Tcorr);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.Tdiff", &Ndata_bb_gem_m4_strip_Tdiff, &b_Ndata_bb_gem_m4_strip_Tdiff);
   fChain->SetBranchAddress("bb.gem.m4.strip.Tdiff", bb_gem_m4_strip_Tdiff, &b_bb_gem_m4_strip_Tdiff);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.Tfit", &Ndata_bb_gem_m4_strip_Tfit, &b_Ndata_bb_gem_m4_strip_Tfit);
   fChain->SetBranchAddress("bb.gem.m4.strip.Tfit", bb_gem_m4_strip_Tfit, &b_bb_gem_m4_strip_Tfit);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.Tmean", &Ndata_bb_gem_m4_strip_Tmean, &b_Ndata_bb_gem_m4_strip_Tmean);
   fChain->SetBranchAddress("bb.gem.m4.strip.Tmean", bb_gem_m4_strip_Tmean, &b_bb_gem_m4_strip_Tmean);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.TmeanDeconv", &Ndata_bb_gem_m4_strip_TmeanDeconv, &b_Ndata_bb_gem_m4_strip_TmeanDeconv);
   fChain->SetBranchAddress("bb.gem.m4.strip.TmeanDeconv", bb_gem_m4_strip_TmeanDeconv, &b_bb_gem_m4_strip_TmeanDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.Tsigma", &Ndata_bb_gem_m4_strip_Tsigma, &b_Ndata_bb_gem_m4_strip_Tsigma);
   fChain->SetBranchAddress("bb.gem.m4.strip.Tsigma", bb_gem_m4_strip_Tsigma, &b_bb_gem_m4_strip_Tsigma);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.adc_id", &Ndata_bb_gem_m4_strip_adc_id, &b_Ndata_bb_gem_m4_strip_adc_id);
   fChain->SetBranchAddress("bb.gem.m4.strip.adc_id", bb_gem_m4_strip_adc_id, &b_bb_gem_m4_strip_adc_id);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.crate", &Ndata_bb_gem_m4_strip_crate, &b_Ndata_bb_gem_m4_strip_crate);
   fChain->SetBranchAddress("bb.gem.m4.strip.crate", bb_gem_m4_strip_crate, &b_bb_gem_m4_strip_crate);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.isampmax", &Ndata_bb_gem_m4_strip_isampmax, &b_Ndata_bb_gem_m4_strip_isampmax);
   fChain->SetBranchAddress("bb.gem.m4.strip.isampmax", bb_gem_m4_strip_isampmax, &b_bb_gem_m4_strip_isampmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.isampmaxDeconv", &Ndata_bb_gem_m4_strip_isampmaxDeconv, &b_Ndata_bb_gem_m4_strip_isampmaxDeconv);
   fChain->SetBranchAddress("bb.gem.m4.strip.isampmaxDeconv", bb_gem_m4_strip_isampmaxDeconv, &b_bb_gem_m4_strip_isampmaxDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.isampmaxDeconvCombo", &Ndata_bb_gem_m4_strip_isampmaxDeconvCombo, &b_Ndata_bb_gem_m4_strip_isampmaxDeconvCombo);
   fChain->SetBranchAddress("bb.gem.m4.strip.isampmaxDeconvCombo", bb_gem_m4_strip_isampmaxDeconvCombo, &b_bb_gem_m4_strip_isampmaxDeconvCombo);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.isnegative", &Ndata_bb_gem_m4_strip_isnegative, &b_Ndata_bb_gem_m4_strip_isnegative);
   fChain->SetBranchAddress("bb.gem.m4.strip.isnegative", bb_gem_m4_strip_isnegative, &b_bb_gem_m4_strip_isnegative);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.isnegativeU", &Ndata_bb_gem_m4_strip_isnegativeU, &b_Ndata_bb_gem_m4_strip_isnegativeU);
   fChain->SetBranchAddress("bb.gem.m4.strip.isnegativeU", bb_gem_m4_strip_isnegativeU, &b_bb_gem_m4_strip_isnegativeU);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.isnegativeV", &Ndata_bb_gem_m4_strip_isnegativeV, &b_Ndata_bb_gem_m4_strip_isnegativeV);
   fChain->SetBranchAddress("bb.gem.m4.strip.isnegativeV", bb_gem_m4_strip_isnegativeV, &b_bb_gem_m4_strip_isnegativeV);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.isnegontrack", &Ndata_bb_gem_m4_strip_isnegontrack, &b_Ndata_bb_gem_m4_strip_isnegontrack);
   fChain->SetBranchAddress("bb.gem.m4.strip.isnegontrack", bb_gem_m4_strip_isnegontrack, &b_bb_gem_m4_strip_isnegontrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.isnegontrackU", &Ndata_bb_gem_m4_strip_isnegontrackU, &b_Ndata_bb_gem_m4_strip_isnegontrackU);
   fChain->SetBranchAddress("bb.gem.m4.strip.isnegontrackU", bb_gem_m4_strip_isnegontrackU, &b_bb_gem_m4_strip_isnegontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.isnegontrackV", &Ndata_bb_gem_m4_strip_isnegontrackV, &b_Ndata_bb_gem_m4_strip_isnegontrackV);
   fChain->SetBranchAddress("bb.gem.m4.strip.isnegontrackV", bb_gem_m4_strip_isnegontrackV, &b_bb_gem_m4_strip_isnegontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.istrip", &Ndata_bb_gem_m4_strip_istrip, &b_Ndata_bb_gem_m4_strip_istrip);
   fChain->SetBranchAddress("bb.gem.m4.strip.istrip", bb_gem_m4_strip_istrip, &b_bb_gem_m4_strip_istrip);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.itrack", &Ndata_bb_gem_m4_strip_itrack, &b_Ndata_bb_gem_m4_strip_itrack);
   fChain->SetBranchAddress("bb.gem.m4.strip.itrack", bb_gem_m4_strip_itrack, &b_bb_gem_m4_strip_itrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.mpd", &Ndata_bb_gem_m4_strip_mpd, &b_Ndata_bb_gem_m4_strip_mpd);
   fChain->SetBranchAddress("bb.gem.m4.strip.mpd", bb_gem_m4_strip_mpd, &b_bb_gem_m4_strip_mpd);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.ontrack", &Ndata_bb_gem_m4_strip_ontrack, &b_Ndata_bb_gem_m4_strip_ontrack);
   fChain->SetBranchAddress("bb.gem.m4.strip.ontrack", bb_gem_m4_strip_ontrack, &b_bb_gem_m4_strip_ontrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.ontrackU", &Ndata_bb_gem_m4_strip_ontrackU, &b_Ndata_bb_gem_m4_strip_ontrackU);
   fChain->SetBranchAddress("bb.gem.m4.strip.ontrackU", bb_gem_m4_strip_ontrackU, &b_bb_gem_m4_strip_ontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.ontrackV", &Ndata_bb_gem_m4_strip_ontrackV, &b_Ndata_bb_gem_m4_strip_ontrackV);
   fChain->SetBranchAddress("bb.gem.m4.strip.ontrackV", bb_gem_m4_strip_ontrackV, &b_bb_gem_m4_strip_ontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m4.strip.rawADCsamples", &Ndata_bb_gem_m4_strip_rawADCsamples, &b_Ndata_bb_gem_m4_strip_rawADCsamples);
   fChain->SetBranchAddress("bb.gem.m4.strip.rawADCsamples", bb_gem_m4_strip_rawADCsamples, &b_bb_gem_m4_strip_rawADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.clust.clustu_adc", &Ndata_bb_gem_m5_clust_clustu_adc, &b_Ndata_bb_gem_m5_clust_clustu_adc);
   fChain->SetBranchAddress("bb.gem.m5.clust.clustu_adc", bb_gem_m5_clust_clustu_adc, &b_bb_gem_m5_clust_clustu_adc);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.clust.clustu_pos", &Ndata_bb_gem_m5_clust_clustu_pos, &b_Ndata_bb_gem_m5_clust_clustu_pos);
   fChain->SetBranchAddress("bb.gem.m5.clust.clustu_pos", bb_gem_m5_clust_clustu_pos, &b_bb_gem_m5_clust_clustu_pos);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.clust.clustu_strips", &Ndata_bb_gem_m5_clust_clustu_strips, &b_Ndata_bb_gem_m5_clust_clustu_strips);
   fChain->SetBranchAddress("bb.gem.m5.clust.clustu_strips", bb_gem_m5_clust_clustu_strips, &b_bb_gem_m5_clust_clustu_strips);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.clust.clustu_time", &Ndata_bb_gem_m5_clust_clustu_time, &b_Ndata_bb_gem_m5_clust_clustu_time);
   fChain->SetBranchAddress("bb.gem.m5.clust.clustu_time", bb_gem_m5_clust_clustu_time, &b_bb_gem_m5_clust_clustu_time);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.clust.clustv_adc", &Ndata_bb_gem_m5_clust_clustv_adc, &b_Ndata_bb_gem_m5_clust_clustv_adc);
   fChain->SetBranchAddress("bb.gem.m5.clust.clustv_adc", bb_gem_m5_clust_clustv_adc, &b_bb_gem_m5_clust_clustv_adc);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.clust.clustv_pos", &Ndata_bb_gem_m5_clust_clustv_pos, &b_Ndata_bb_gem_m5_clust_clustv_pos);
   fChain->SetBranchAddress("bb.gem.m5.clust.clustv_pos", bb_gem_m5_clust_clustv_pos, &b_bb_gem_m5_clust_clustv_pos);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.clust.clustv_strips", &Ndata_bb_gem_m5_clust_clustv_strips, &b_Ndata_bb_gem_m5_clust_clustv_strips);
   fChain->SetBranchAddress("bb.gem.m5.clust.clustv_strips", bb_gem_m5_clust_clustv_strips, &b_bb_gem_m5_clust_clustv_strips);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.clust.clustv_time", &Ndata_bb_gem_m5_clust_clustv_time, &b_Ndata_bb_gem_m5_clust_clustv_time);
   fChain->SetBranchAddress("bb.gem.m5.clust.clustv_time", bb_gem_m5_clust_clustv_time, &b_bb_gem_m5_clust_clustv_time);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.clust.isnegativeU", &Ndata_bb_gem_m5_clust_isnegativeU, &b_Ndata_bb_gem_m5_clust_isnegativeU);
   fChain->SetBranchAddress("bb.gem.m5.clust.isnegativeU", bb_gem_m5_clust_isnegativeU, &b_bb_gem_m5_clust_isnegativeU);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.clust.isnegativeV", &Ndata_bb_gem_m5_clust_isnegativeV, &b_Ndata_bb_gem_m5_clust_isnegativeV);
   fChain->SetBranchAddress("bb.gem.m5.clust.isnegativeV", bb_gem_m5_clust_isnegativeV, &b_bb_gem_m5_clust_isnegativeV);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.clust.isnegontrackU", &Ndata_bb_gem_m5_clust_isnegontrackU, &b_Ndata_bb_gem_m5_clust_isnegontrackU);
   fChain->SetBranchAddress("bb.gem.m5.clust.isnegontrackU", bb_gem_m5_clust_isnegontrackU, &b_bb_gem_m5_clust_isnegontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.clust.isnegontrackV", &Ndata_bb_gem_m5_clust_isnegontrackV, &b_Ndata_bb_gem_m5_clust_isnegontrackV);
   fChain->SetBranchAddress("bb.gem.m5.clust.isnegontrackV", bb_gem_m5_clust_isnegontrackV, &b_bb_gem_m5_clust_isnegontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.ADCavg", &Ndata_bb_gem_m5_strip_ADCavg, &b_Ndata_bb_gem_m5_strip_ADCavg);
   fChain->SetBranchAddress("bb.gem.m5.strip.ADCavg", bb_gem_m5_strip_ADCavg, &b_bb_gem_m5_strip_ADCavg);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.ADCmax", &Ndata_bb_gem_m5_strip_ADCmax, &b_Ndata_bb_gem_m5_strip_ADCmax);
   fChain->SetBranchAddress("bb.gem.m5.strip.ADCmax", bb_gem_m5_strip_ADCmax, &b_bb_gem_m5_strip_ADCmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.ADCsamples", &Ndata_bb_gem_m5_strip_ADCsamples, &b_Ndata_bb_gem_m5_strip_ADCsamples);
   fChain->SetBranchAddress("bb.gem.m5.strip.ADCsamples", bb_gem_m5_strip_ADCsamples, &b_bb_gem_m5_strip_ADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.ADCsum", &Ndata_bb_gem_m5_strip_ADCsum, &b_Ndata_bb_gem_m5_strip_ADCsum);
   fChain->SetBranchAddress("bb.gem.m5.strip.ADCsum", bb_gem_m5_strip_ADCsum, &b_bb_gem_m5_strip_ADCsum);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.BUILD_ALL_SAMPLES", &Ndata_bb_gem_m5_strip_BUILD_ALL_SAMPLES, &b_Ndata_bb_gem_m5_strip_BUILD_ALL_SAMPLES);
   fChain->SetBranchAddress("bb.gem.m5.strip.BUILD_ALL_SAMPLES", bb_gem_m5_strip_BUILD_ALL_SAMPLES, &b_bb_gem_m5_strip_BUILD_ALL_SAMPLES);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.CM_GOOD", &Ndata_bb_gem_m5_strip_CM_GOOD, &b_Ndata_bb_gem_m5_strip_CM_GOOD);
   fChain->SetBranchAddress("bb.gem.m5.strip.CM_GOOD", bb_gem_m5_strip_CM_GOOD, &b_bb_gem_m5_strip_CM_GOOD);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.CorrCoeff", &Ndata_bb_gem_m5_strip_CorrCoeff, &b_Ndata_bb_gem_m5_strip_CorrCoeff);
   fChain->SetBranchAddress("bb.gem.m5.strip.CorrCoeff", bb_gem_m5_strip_CorrCoeff, &b_bb_gem_m5_strip_CorrCoeff);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.DeconvADCmax", &Ndata_bb_gem_m5_strip_DeconvADCmax, &b_Ndata_bb_gem_m5_strip_DeconvADCmax);
   fChain->SetBranchAddress("bb.gem.m5.strip.DeconvADCmax", bb_gem_m5_strip_DeconvADCmax, &b_bb_gem_m5_strip_DeconvADCmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.DeconvADCmaxCombo", &Ndata_bb_gem_m5_strip_DeconvADCmaxCombo, &b_Ndata_bb_gem_m5_strip_DeconvADCmaxCombo);
   fChain->SetBranchAddress("bb.gem.m5.strip.DeconvADCmaxCombo", bb_gem_m5_strip_DeconvADCmaxCombo, &b_bb_gem_m5_strip_DeconvADCmaxCombo);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.DeconvADCsamples", &Ndata_bb_gem_m5_strip_DeconvADCsamples, &b_Ndata_bb_gem_m5_strip_DeconvADCsamples);
   fChain->SetBranchAddress("bb.gem.m5.strip.DeconvADCsamples", bb_gem_m5_strip_DeconvADCsamples, &b_bb_gem_m5_strip_DeconvADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.DeconvADCsum", &Ndata_bb_gem_m5_strip_DeconvADCsum, &b_Ndata_bb_gem_m5_strip_DeconvADCsum);
   fChain->SetBranchAddress("bb.gem.m5.strip.DeconvADCsum", bb_gem_m5_strip_DeconvADCsum, &b_bb_gem_m5_strip_DeconvADCsum);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.ENABLE_CM", &Ndata_bb_gem_m5_strip_ENABLE_CM, &b_Ndata_bb_gem_m5_strip_ENABLE_CM);
   fChain->SetBranchAddress("bb.gem.m5.strip.ENABLE_CM", bb_gem_m5_strip_ENABLE_CM, &b_bb_gem_m5_strip_ENABLE_CM);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.IsU", &Ndata_bb_gem_m5_strip_IsU, &b_Ndata_bb_gem_m5_strip_IsU);
   fChain->SetBranchAddress("bb.gem.m5.strip.IsU", bb_gem_m5_strip_IsU, &b_bb_gem_m5_strip_IsU);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.IsV", &Ndata_bb_gem_m5_strip_IsV, &b_Ndata_bb_gem_m5_strip_IsV);
   fChain->SetBranchAddress("bb.gem.m5.strip.IsV", bb_gem_m5_strip_IsV, &b_bb_gem_m5_strip_IsV);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.TSchi2", &Ndata_bb_gem_m5_strip_TSchi2, &b_Ndata_bb_gem_m5_strip_TSchi2);
   fChain->SetBranchAddress("bb.gem.m5.strip.TSchi2", bb_gem_m5_strip_TSchi2, &b_bb_gem_m5_strip_TSchi2);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.TSprob", &Ndata_bb_gem_m5_strip_TSprob, &b_Ndata_bb_gem_m5_strip_TSprob);
   fChain->SetBranchAddress("bb.gem.m5.strip.TSprob", bb_gem_m5_strip_TSprob, &b_bb_gem_m5_strip_TSprob);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.Tcorr", &Ndata_bb_gem_m5_strip_Tcorr, &b_Ndata_bb_gem_m5_strip_Tcorr);
   fChain->SetBranchAddress("bb.gem.m5.strip.Tcorr", bb_gem_m5_strip_Tcorr, &b_bb_gem_m5_strip_Tcorr);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.Tdiff", &Ndata_bb_gem_m5_strip_Tdiff, &b_Ndata_bb_gem_m5_strip_Tdiff);
   fChain->SetBranchAddress("bb.gem.m5.strip.Tdiff", bb_gem_m5_strip_Tdiff, &b_bb_gem_m5_strip_Tdiff);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.Tfit", &Ndata_bb_gem_m5_strip_Tfit, &b_Ndata_bb_gem_m5_strip_Tfit);
   fChain->SetBranchAddress("bb.gem.m5.strip.Tfit", bb_gem_m5_strip_Tfit, &b_bb_gem_m5_strip_Tfit);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.Tmean", &Ndata_bb_gem_m5_strip_Tmean, &b_Ndata_bb_gem_m5_strip_Tmean);
   fChain->SetBranchAddress("bb.gem.m5.strip.Tmean", bb_gem_m5_strip_Tmean, &b_bb_gem_m5_strip_Tmean);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.TmeanDeconv", &Ndata_bb_gem_m5_strip_TmeanDeconv, &b_Ndata_bb_gem_m5_strip_TmeanDeconv);
   fChain->SetBranchAddress("bb.gem.m5.strip.TmeanDeconv", bb_gem_m5_strip_TmeanDeconv, &b_bb_gem_m5_strip_TmeanDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.Tsigma", &Ndata_bb_gem_m5_strip_Tsigma, &b_Ndata_bb_gem_m5_strip_Tsigma);
   fChain->SetBranchAddress("bb.gem.m5.strip.Tsigma", bb_gem_m5_strip_Tsigma, &b_bb_gem_m5_strip_Tsigma);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.adc_id", &Ndata_bb_gem_m5_strip_adc_id, &b_Ndata_bb_gem_m5_strip_adc_id);
   fChain->SetBranchAddress("bb.gem.m5.strip.adc_id", bb_gem_m5_strip_adc_id, &b_bb_gem_m5_strip_adc_id);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.crate", &Ndata_bb_gem_m5_strip_crate, &b_Ndata_bb_gem_m5_strip_crate);
   fChain->SetBranchAddress("bb.gem.m5.strip.crate", bb_gem_m5_strip_crate, &b_bb_gem_m5_strip_crate);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.isampmax", &Ndata_bb_gem_m5_strip_isampmax, &b_Ndata_bb_gem_m5_strip_isampmax);
   fChain->SetBranchAddress("bb.gem.m5.strip.isampmax", bb_gem_m5_strip_isampmax, &b_bb_gem_m5_strip_isampmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.isampmaxDeconv", &Ndata_bb_gem_m5_strip_isampmaxDeconv, &b_Ndata_bb_gem_m5_strip_isampmaxDeconv);
   fChain->SetBranchAddress("bb.gem.m5.strip.isampmaxDeconv", bb_gem_m5_strip_isampmaxDeconv, &b_bb_gem_m5_strip_isampmaxDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.isampmaxDeconvCombo", &Ndata_bb_gem_m5_strip_isampmaxDeconvCombo, &b_Ndata_bb_gem_m5_strip_isampmaxDeconvCombo);
   fChain->SetBranchAddress("bb.gem.m5.strip.isampmaxDeconvCombo", bb_gem_m5_strip_isampmaxDeconvCombo, &b_bb_gem_m5_strip_isampmaxDeconvCombo);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.isnegative", &Ndata_bb_gem_m5_strip_isnegative, &b_Ndata_bb_gem_m5_strip_isnegative);
   fChain->SetBranchAddress("bb.gem.m5.strip.isnegative", bb_gem_m5_strip_isnegative, &b_bb_gem_m5_strip_isnegative);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.isnegativeU", &Ndata_bb_gem_m5_strip_isnegativeU, &b_Ndata_bb_gem_m5_strip_isnegativeU);
   fChain->SetBranchAddress("bb.gem.m5.strip.isnegativeU", bb_gem_m5_strip_isnegativeU, &b_bb_gem_m5_strip_isnegativeU);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.isnegativeV", &Ndata_bb_gem_m5_strip_isnegativeV, &b_Ndata_bb_gem_m5_strip_isnegativeV);
   fChain->SetBranchAddress("bb.gem.m5.strip.isnegativeV", bb_gem_m5_strip_isnegativeV, &b_bb_gem_m5_strip_isnegativeV);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.isnegontrack", &Ndata_bb_gem_m5_strip_isnegontrack, &b_Ndata_bb_gem_m5_strip_isnegontrack);
   fChain->SetBranchAddress("bb.gem.m5.strip.isnegontrack", bb_gem_m5_strip_isnegontrack, &b_bb_gem_m5_strip_isnegontrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.isnegontrackU", &Ndata_bb_gem_m5_strip_isnegontrackU, &b_Ndata_bb_gem_m5_strip_isnegontrackU);
   fChain->SetBranchAddress("bb.gem.m5.strip.isnegontrackU", bb_gem_m5_strip_isnegontrackU, &b_bb_gem_m5_strip_isnegontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.isnegontrackV", &Ndata_bb_gem_m5_strip_isnegontrackV, &b_Ndata_bb_gem_m5_strip_isnegontrackV);
   fChain->SetBranchAddress("bb.gem.m5.strip.isnegontrackV", bb_gem_m5_strip_isnegontrackV, &b_bb_gem_m5_strip_isnegontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.istrip", &Ndata_bb_gem_m5_strip_istrip, &b_Ndata_bb_gem_m5_strip_istrip);
   fChain->SetBranchAddress("bb.gem.m5.strip.istrip", bb_gem_m5_strip_istrip, &b_bb_gem_m5_strip_istrip);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.itrack", &Ndata_bb_gem_m5_strip_itrack, &b_Ndata_bb_gem_m5_strip_itrack);
   fChain->SetBranchAddress("bb.gem.m5.strip.itrack", bb_gem_m5_strip_itrack, &b_bb_gem_m5_strip_itrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.mpd", &Ndata_bb_gem_m5_strip_mpd, &b_Ndata_bb_gem_m5_strip_mpd);
   fChain->SetBranchAddress("bb.gem.m5.strip.mpd", bb_gem_m5_strip_mpd, &b_bb_gem_m5_strip_mpd);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.ontrack", &Ndata_bb_gem_m5_strip_ontrack, &b_Ndata_bb_gem_m5_strip_ontrack);
   fChain->SetBranchAddress("bb.gem.m5.strip.ontrack", bb_gem_m5_strip_ontrack, &b_bb_gem_m5_strip_ontrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.ontrackU", &Ndata_bb_gem_m5_strip_ontrackU, &b_Ndata_bb_gem_m5_strip_ontrackU);
   fChain->SetBranchAddress("bb.gem.m5.strip.ontrackU", bb_gem_m5_strip_ontrackU, &b_bb_gem_m5_strip_ontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.ontrackV", &Ndata_bb_gem_m5_strip_ontrackV, &b_Ndata_bb_gem_m5_strip_ontrackV);
   fChain->SetBranchAddress("bb.gem.m5.strip.ontrackV", bb_gem_m5_strip_ontrackV, &b_bb_gem_m5_strip_ontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m5.strip.rawADCsamples", &Ndata_bb_gem_m5_strip_rawADCsamples, &b_Ndata_bb_gem_m5_strip_rawADCsamples);
   fChain->SetBranchAddress("bb.gem.m5.strip.rawADCsamples", bb_gem_m5_strip_rawADCsamples, &b_bb_gem_m5_strip_rawADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.clust.clustu_adc", &Ndata_bb_gem_m6_clust_clustu_adc, &b_Ndata_bb_gem_m6_clust_clustu_adc);
   fChain->SetBranchAddress("bb.gem.m6.clust.clustu_adc", bb_gem_m6_clust_clustu_adc, &b_bb_gem_m6_clust_clustu_adc);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.clust.clustu_pos", &Ndata_bb_gem_m6_clust_clustu_pos, &b_Ndata_bb_gem_m6_clust_clustu_pos);
   fChain->SetBranchAddress("bb.gem.m6.clust.clustu_pos", bb_gem_m6_clust_clustu_pos, &b_bb_gem_m6_clust_clustu_pos);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.clust.clustu_strips", &Ndata_bb_gem_m6_clust_clustu_strips, &b_Ndata_bb_gem_m6_clust_clustu_strips);
   fChain->SetBranchAddress("bb.gem.m6.clust.clustu_strips", bb_gem_m6_clust_clustu_strips, &b_bb_gem_m6_clust_clustu_strips);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.clust.clustu_time", &Ndata_bb_gem_m6_clust_clustu_time, &b_Ndata_bb_gem_m6_clust_clustu_time);
   fChain->SetBranchAddress("bb.gem.m6.clust.clustu_time", bb_gem_m6_clust_clustu_time, &b_bb_gem_m6_clust_clustu_time);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.clust.clustv_adc", &Ndata_bb_gem_m6_clust_clustv_adc, &b_Ndata_bb_gem_m6_clust_clustv_adc);
   fChain->SetBranchAddress("bb.gem.m6.clust.clustv_adc", bb_gem_m6_clust_clustv_adc, &b_bb_gem_m6_clust_clustv_adc);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.clust.clustv_pos", &Ndata_bb_gem_m6_clust_clustv_pos, &b_Ndata_bb_gem_m6_clust_clustv_pos);
   fChain->SetBranchAddress("bb.gem.m6.clust.clustv_pos", bb_gem_m6_clust_clustv_pos, &b_bb_gem_m6_clust_clustv_pos);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.clust.clustv_strips", &Ndata_bb_gem_m6_clust_clustv_strips, &b_Ndata_bb_gem_m6_clust_clustv_strips);
   fChain->SetBranchAddress("bb.gem.m6.clust.clustv_strips", bb_gem_m6_clust_clustv_strips, &b_bb_gem_m6_clust_clustv_strips);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.clust.clustv_time", &Ndata_bb_gem_m6_clust_clustv_time, &b_Ndata_bb_gem_m6_clust_clustv_time);
   fChain->SetBranchAddress("bb.gem.m6.clust.clustv_time", bb_gem_m6_clust_clustv_time, &b_bb_gem_m6_clust_clustv_time);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.clust.isnegativeU", &Ndata_bb_gem_m6_clust_isnegativeU, &b_Ndata_bb_gem_m6_clust_isnegativeU);
   fChain->SetBranchAddress("bb.gem.m6.clust.isnegativeU", bb_gem_m6_clust_isnegativeU, &b_bb_gem_m6_clust_isnegativeU);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.clust.isnegativeV", &Ndata_bb_gem_m6_clust_isnegativeV, &b_Ndata_bb_gem_m6_clust_isnegativeV);
   fChain->SetBranchAddress("bb.gem.m6.clust.isnegativeV", bb_gem_m6_clust_isnegativeV, &b_bb_gem_m6_clust_isnegativeV);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.clust.isnegontrackU", &Ndata_bb_gem_m6_clust_isnegontrackU, &b_Ndata_bb_gem_m6_clust_isnegontrackU);
   fChain->SetBranchAddress("bb.gem.m6.clust.isnegontrackU", bb_gem_m6_clust_isnegontrackU, &b_bb_gem_m6_clust_isnegontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.clust.isnegontrackV", &Ndata_bb_gem_m6_clust_isnegontrackV, &b_Ndata_bb_gem_m6_clust_isnegontrackV);
   fChain->SetBranchAddress("bb.gem.m6.clust.isnegontrackV", bb_gem_m6_clust_isnegontrackV, &b_bb_gem_m6_clust_isnegontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.ADCavg", &Ndata_bb_gem_m6_strip_ADCavg, &b_Ndata_bb_gem_m6_strip_ADCavg);
   fChain->SetBranchAddress("bb.gem.m6.strip.ADCavg", bb_gem_m6_strip_ADCavg, &b_bb_gem_m6_strip_ADCavg);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.ADCmax", &Ndata_bb_gem_m6_strip_ADCmax, &b_Ndata_bb_gem_m6_strip_ADCmax);
   fChain->SetBranchAddress("bb.gem.m6.strip.ADCmax", bb_gem_m6_strip_ADCmax, &b_bb_gem_m6_strip_ADCmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.ADCsamples", &Ndata_bb_gem_m6_strip_ADCsamples, &b_Ndata_bb_gem_m6_strip_ADCsamples);
   fChain->SetBranchAddress("bb.gem.m6.strip.ADCsamples", bb_gem_m6_strip_ADCsamples, &b_bb_gem_m6_strip_ADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.ADCsum", &Ndata_bb_gem_m6_strip_ADCsum, &b_Ndata_bb_gem_m6_strip_ADCsum);
   fChain->SetBranchAddress("bb.gem.m6.strip.ADCsum", bb_gem_m6_strip_ADCsum, &b_bb_gem_m6_strip_ADCsum);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.BUILD_ALL_SAMPLES", &Ndata_bb_gem_m6_strip_BUILD_ALL_SAMPLES, &b_Ndata_bb_gem_m6_strip_BUILD_ALL_SAMPLES);
   fChain->SetBranchAddress("bb.gem.m6.strip.BUILD_ALL_SAMPLES", bb_gem_m6_strip_BUILD_ALL_SAMPLES, &b_bb_gem_m6_strip_BUILD_ALL_SAMPLES);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.CM_GOOD", &Ndata_bb_gem_m6_strip_CM_GOOD, &b_Ndata_bb_gem_m6_strip_CM_GOOD);
   fChain->SetBranchAddress("bb.gem.m6.strip.CM_GOOD", bb_gem_m6_strip_CM_GOOD, &b_bb_gem_m6_strip_CM_GOOD);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.CorrCoeff", &Ndata_bb_gem_m6_strip_CorrCoeff, &b_Ndata_bb_gem_m6_strip_CorrCoeff);
   fChain->SetBranchAddress("bb.gem.m6.strip.CorrCoeff", bb_gem_m6_strip_CorrCoeff, &b_bb_gem_m6_strip_CorrCoeff);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.DeconvADCmax", &Ndata_bb_gem_m6_strip_DeconvADCmax, &b_Ndata_bb_gem_m6_strip_DeconvADCmax);
   fChain->SetBranchAddress("bb.gem.m6.strip.DeconvADCmax", bb_gem_m6_strip_DeconvADCmax, &b_bb_gem_m6_strip_DeconvADCmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.DeconvADCmaxCombo", &Ndata_bb_gem_m6_strip_DeconvADCmaxCombo, &b_Ndata_bb_gem_m6_strip_DeconvADCmaxCombo);
   fChain->SetBranchAddress("bb.gem.m6.strip.DeconvADCmaxCombo", bb_gem_m6_strip_DeconvADCmaxCombo, &b_bb_gem_m6_strip_DeconvADCmaxCombo);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.DeconvADCsamples", &Ndata_bb_gem_m6_strip_DeconvADCsamples, &b_Ndata_bb_gem_m6_strip_DeconvADCsamples);
   fChain->SetBranchAddress("bb.gem.m6.strip.DeconvADCsamples", bb_gem_m6_strip_DeconvADCsamples, &b_bb_gem_m6_strip_DeconvADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.DeconvADCsum", &Ndata_bb_gem_m6_strip_DeconvADCsum, &b_Ndata_bb_gem_m6_strip_DeconvADCsum);
   fChain->SetBranchAddress("bb.gem.m6.strip.DeconvADCsum", bb_gem_m6_strip_DeconvADCsum, &b_bb_gem_m6_strip_DeconvADCsum);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.ENABLE_CM", &Ndata_bb_gem_m6_strip_ENABLE_CM, &b_Ndata_bb_gem_m6_strip_ENABLE_CM);
   fChain->SetBranchAddress("bb.gem.m6.strip.ENABLE_CM", bb_gem_m6_strip_ENABLE_CM, &b_bb_gem_m6_strip_ENABLE_CM);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.IsU", &Ndata_bb_gem_m6_strip_IsU, &b_Ndata_bb_gem_m6_strip_IsU);
   fChain->SetBranchAddress("bb.gem.m6.strip.IsU", bb_gem_m6_strip_IsU, &b_bb_gem_m6_strip_IsU);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.IsV", &Ndata_bb_gem_m6_strip_IsV, &b_Ndata_bb_gem_m6_strip_IsV);
   fChain->SetBranchAddress("bb.gem.m6.strip.IsV", bb_gem_m6_strip_IsV, &b_bb_gem_m6_strip_IsV);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.TSchi2", &Ndata_bb_gem_m6_strip_TSchi2, &b_Ndata_bb_gem_m6_strip_TSchi2);
   fChain->SetBranchAddress("bb.gem.m6.strip.TSchi2", bb_gem_m6_strip_TSchi2, &b_bb_gem_m6_strip_TSchi2);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.TSprob", &Ndata_bb_gem_m6_strip_TSprob, &b_Ndata_bb_gem_m6_strip_TSprob);
   fChain->SetBranchAddress("bb.gem.m6.strip.TSprob", bb_gem_m6_strip_TSprob, &b_bb_gem_m6_strip_TSprob);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.Tcorr", &Ndata_bb_gem_m6_strip_Tcorr, &b_Ndata_bb_gem_m6_strip_Tcorr);
   fChain->SetBranchAddress("bb.gem.m6.strip.Tcorr", bb_gem_m6_strip_Tcorr, &b_bb_gem_m6_strip_Tcorr);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.Tdiff", &Ndata_bb_gem_m6_strip_Tdiff, &b_Ndata_bb_gem_m6_strip_Tdiff);
   fChain->SetBranchAddress("bb.gem.m6.strip.Tdiff", bb_gem_m6_strip_Tdiff, &b_bb_gem_m6_strip_Tdiff);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.Tfit", &Ndata_bb_gem_m6_strip_Tfit, &b_Ndata_bb_gem_m6_strip_Tfit);
   fChain->SetBranchAddress("bb.gem.m6.strip.Tfit", bb_gem_m6_strip_Tfit, &b_bb_gem_m6_strip_Tfit);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.Tmean", &Ndata_bb_gem_m6_strip_Tmean, &b_Ndata_bb_gem_m6_strip_Tmean);
   fChain->SetBranchAddress("bb.gem.m6.strip.Tmean", bb_gem_m6_strip_Tmean, &b_bb_gem_m6_strip_Tmean);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.TmeanDeconv", &Ndata_bb_gem_m6_strip_TmeanDeconv, &b_Ndata_bb_gem_m6_strip_TmeanDeconv);
   fChain->SetBranchAddress("bb.gem.m6.strip.TmeanDeconv", bb_gem_m6_strip_TmeanDeconv, &b_bb_gem_m6_strip_TmeanDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.Tsigma", &Ndata_bb_gem_m6_strip_Tsigma, &b_Ndata_bb_gem_m6_strip_Tsigma);
   fChain->SetBranchAddress("bb.gem.m6.strip.Tsigma", bb_gem_m6_strip_Tsigma, &b_bb_gem_m6_strip_Tsigma);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.adc_id", &Ndata_bb_gem_m6_strip_adc_id, &b_Ndata_bb_gem_m6_strip_adc_id);
   fChain->SetBranchAddress("bb.gem.m6.strip.adc_id", bb_gem_m6_strip_adc_id, &b_bb_gem_m6_strip_adc_id);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.crate", &Ndata_bb_gem_m6_strip_crate, &b_Ndata_bb_gem_m6_strip_crate);
   fChain->SetBranchAddress("bb.gem.m6.strip.crate", bb_gem_m6_strip_crate, &b_bb_gem_m6_strip_crate);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.isampmax", &Ndata_bb_gem_m6_strip_isampmax, &b_Ndata_bb_gem_m6_strip_isampmax);
   fChain->SetBranchAddress("bb.gem.m6.strip.isampmax", bb_gem_m6_strip_isampmax, &b_bb_gem_m6_strip_isampmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.isampmaxDeconv", &Ndata_bb_gem_m6_strip_isampmaxDeconv, &b_Ndata_bb_gem_m6_strip_isampmaxDeconv);
   fChain->SetBranchAddress("bb.gem.m6.strip.isampmaxDeconv", bb_gem_m6_strip_isampmaxDeconv, &b_bb_gem_m6_strip_isampmaxDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.isampmaxDeconvCombo", &Ndata_bb_gem_m6_strip_isampmaxDeconvCombo, &b_Ndata_bb_gem_m6_strip_isampmaxDeconvCombo);
   fChain->SetBranchAddress("bb.gem.m6.strip.isampmaxDeconvCombo", bb_gem_m6_strip_isampmaxDeconvCombo, &b_bb_gem_m6_strip_isampmaxDeconvCombo);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.isnegative", &Ndata_bb_gem_m6_strip_isnegative, &b_Ndata_bb_gem_m6_strip_isnegative);
   fChain->SetBranchAddress("bb.gem.m6.strip.isnegative", bb_gem_m6_strip_isnegative, &b_bb_gem_m6_strip_isnegative);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.isnegativeU", &Ndata_bb_gem_m6_strip_isnegativeU, &b_Ndata_bb_gem_m6_strip_isnegativeU);
   fChain->SetBranchAddress("bb.gem.m6.strip.isnegativeU", bb_gem_m6_strip_isnegativeU, &b_bb_gem_m6_strip_isnegativeU);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.isnegativeV", &Ndata_bb_gem_m6_strip_isnegativeV, &b_Ndata_bb_gem_m6_strip_isnegativeV);
   fChain->SetBranchAddress("bb.gem.m6.strip.isnegativeV", bb_gem_m6_strip_isnegativeV, &b_bb_gem_m6_strip_isnegativeV);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.isnegontrack", &Ndata_bb_gem_m6_strip_isnegontrack, &b_Ndata_bb_gem_m6_strip_isnegontrack);
   fChain->SetBranchAddress("bb.gem.m6.strip.isnegontrack", bb_gem_m6_strip_isnegontrack, &b_bb_gem_m6_strip_isnegontrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.isnegontrackU", &Ndata_bb_gem_m6_strip_isnegontrackU, &b_Ndata_bb_gem_m6_strip_isnegontrackU);
   fChain->SetBranchAddress("bb.gem.m6.strip.isnegontrackU", bb_gem_m6_strip_isnegontrackU, &b_bb_gem_m6_strip_isnegontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.isnegontrackV", &Ndata_bb_gem_m6_strip_isnegontrackV, &b_Ndata_bb_gem_m6_strip_isnegontrackV);
   fChain->SetBranchAddress("bb.gem.m6.strip.isnegontrackV", bb_gem_m6_strip_isnegontrackV, &b_bb_gem_m6_strip_isnegontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.istrip", &Ndata_bb_gem_m6_strip_istrip, &b_Ndata_bb_gem_m6_strip_istrip);
   fChain->SetBranchAddress("bb.gem.m6.strip.istrip", bb_gem_m6_strip_istrip, &b_bb_gem_m6_strip_istrip);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.itrack", &Ndata_bb_gem_m6_strip_itrack, &b_Ndata_bb_gem_m6_strip_itrack);
   fChain->SetBranchAddress("bb.gem.m6.strip.itrack", bb_gem_m6_strip_itrack, &b_bb_gem_m6_strip_itrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.mpd", &Ndata_bb_gem_m6_strip_mpd, &b_Ndata_bb_gem_m6_strip_mpd);
   fChain->SetBranchAddress("bb.gem.m6.strip.mpd", bb_gem_m6_strip_mpd, &b_bb_gem_m6_strip_mpd);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.ontrack", &Ndata_bb_gem_m6_strip_ontrack, &b_Ndata_bb_gem_m6_strip_ontrack);
   fChain->SetBranchAddress("bb.gem.m6.strip.ontrack", bb_gem_m6_strip_ontrack, &b_bb_gem_m6_strip_ontrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.ontrackU", &Ndata_bb_gem_m6_strip_ontrackU, &b_Ndata_bb_gem_m6_strip_ontrackU);
   fChain->SetBranchAddress("bb.gem.m6.strip.ontrackU", bb_gem_m6_strip_ontrackU, &b_bb_gem_m6_strip_ontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.ontrackV", &Ndata_bb_gem_m6_strip_ontrackV, &b_Ndata_bb_gem_m6_strip_ontrackV);
   fChain->SetBranchAddress("bb.gem.m6.strip.ontrackV", bb_gem_m6_strip_ontrackV, &b_bb_gem_m6_strip_ontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m6.strip.rawADCsamples", &Ndata_bb_gem_m6_strip_rawADCsamples, &b_Ndata_bb_gem_m6_strip_rawADCsamples);
   fChain->SetBranchAddress("bb.gem.m6.strip.rawADCsamples", bb_gem_m6_strip_rawADCsamples, &b_bb_gem_m6_strip_rawADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.clust.clustu_adc", &Ndata_bb_gem_m7_clust_clustu_adc, &b_Ndata_bb_gem_m7_clust_clustu_adc);
   fChain->SetBranchAddress("bb.gem.m7.clust.clustu_adc", bb_gem_m7_clust_clustu_adc, &b_bb_gem_m7_clust_clustu_adc);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.clust.clustu_pos", &Ndata_bb_gem_m7_clust_clustu_pos, &b_Ndata_bb_gem_m7_clust_clustu_pos);
   fChain->SetBranchAddress("bb.gem.m7.clust.clustu_pos", bb_gem_m7_clust_clustu_pos, &b_bb_gem_m7_clust_clustu_pos);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.clust.clustu_strips", &Ndata_bb_gem_m7_clust_clustu_strips, &b_Ndata_bb_gem_m7_clust_clustu_strips);
   fChain->SetBranchAddress("bb.gem.m7.clust.clustu_strips", bb_gem_m7_clust_clustu_strips, &b_bb_gem_m7_clust_clustu_strips);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.clust.clustu_time", &Ndata_bb_gem_m7_clust_clustu_time, &b_Ndata_bb_gem_m7_clust_clustu_time);
   fChain->SetBranchAddress("bb.gem.m7.clust.clustu_time", bb_gem_m7_clust_clustu_time, &b_bb_gem_m7_clust_clustu_time);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.clust.clustv_adc", &Ndata_bb_gem_m7_clust_clustv_adc, &b_Ndata_bb_gem_m7_clust_clustv_adc);
   fChain->SetBranchAddress("bb.gem.m7.clust.clustv_adc", bb_gem_m7_clust_clustv_adc, &b_bb_gem_m7_clust_clustv_adc);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.clust.clustv_pos", &Ndata_bb_gem_m7_clust_clustv_pos, &b_Ndata_bb_gem_m7_clust_clustv_pos);
   fChain->SetBranchAddress("bb.gem.m7.clust.clustv_pos", bb_gem_m7_clust_clustv_pos, &b_bb_gem_m7_clust_clustv_pos);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.clust.clustv_strips", &Ndata_bb_gem_m7_clust_clustv_strips, &b_Ndata_bb_gem_m7_clust_clustv_strips);
   fChain->SetBranchAddress("bb.gem.m7.clust.clustv_strips", bb_gem_m7_clust_clustv_strips, &b_bb_gem_m7_clust_clustv_strips);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.clust.clustv_time", &Ndata_bb_gem_m7_clust_clustv_time, &b_Ndata_bb_gem_m7_clust_clustv_time);
   fChain->SetBranchAddress("bb.gem.m7.clust.clustv_time", bb_gem_m7_clust_clustv_time, &b_bb_gem_m7_clust_clustv_time);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.clust.isnegativeU", &Ndata_bb_gem_m7_clust_isnegativeU, &b_Ndata_bb_gem_m7_clust_isnegativeU);
   fChain->SetBranchAddress("bb.gem.m7.clust.isnegativeU", bb_gem_m7_clust_isnegativeU, &b_bb_gem_m7_clust_isnegativeU);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.clust.isnegativeV", &Ndata_bb_gem_m7_clust_isnegativeV, &b_Ndata_bb_gem_m7_clust_isnegativeV);
   fChain->SetBranchAddress("bb.gem.m7.clust.isnegativeV", bb_gem_m7_clust_isnegativeV, &b_bb_gem_m7_clust_isnegativeV);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.clust.isnegontrackU", &Ndata_bb_gem_m7_clust_isnegontrackU, &b_Ndata_bb_gem_m7_clust_isnegontrackU);
   fChain->SetBranchAddress("bb.gem.m7.clust.isnegontrackU", bb_gem_m7_clust_isnegontrackU, &b_bb_gem_m7_clust_isnegontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.clust.isnegontrackV", &Ndata_bb_gem_m7_clust_isnegontrackV, &b_Ndata_bb_gem_m7_clust_isnegontrackV);
   fChain->SetBranchAddress("bb.gem.m7.clust.isnegontrackV", bb_gem_m7_clust_isnegontrackV, &b_bb_gem_m7_clust_isnegontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.ADCavg", &Ndata_bb_gem_m7_strip_ADCavg, &b_Ndata_bb_gem_m7_strip_ADCavg);
   fChain->SetBranchAddress("bb.gem.m7.strip.ADCavg", bb_gem_m7_strip_ADCavg, &b_bb_gem_m7_strip_ADCavg);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.ADCmax", &Ndata_bb_gem_m7_strip_ADCmax, &b_Ndata_bb_gem_m7_strip_ADCmax);
   fChain->SetBranchAddress("bb.gem.m7.strip.ADCmax", bb_gem_m7_strip_ADCmax, &b_bb_gem_m7_strip_ADCmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.ADCsamples", &Ndata_bb_gem_m7_strip_ADCsamples, &b_Ndata_bb_gem_m7_strip_ADCsamples);
   fChain->SetBranchAddress("bb.gem.m7.strip.ADCsamples", bb_gem_m7_strip_ADCsamples, &b_bb_gem_m7_strip_ADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.ADCsum", &Ndata_bb_gem_m7_strip_ADCsum, &b_Ndata_bb_gem_m7_strip_ADCsum);
   fChain->SetBranchAddress("bb.gem.m7.strip.ADCsum", bb_gem_m7_strip_ADCsum, &b_bb_gem_m7_strip_ADCsum);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.BUILD_ALL_SAMPLES", &Ndata_bb_gem_m7_strip_BUILD_ALL_SAMPLES, &b_Ndata_bb_gem_m7_strip_BUILD_ALL_SAMPLES);
   fChain->SetBranchAddress("bb.gem.m7.strip.BUILD_ALL_SAMPLES", bb_gem_m7_strip_BUILD_ALL_SAMPLES, &b_bb_gem_m7_strip_BUILD_ALL_SAMPLES);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.CM_GOOD", &Ndata_bb_gem_m7_strip_CM_GOOD, &b_Ndata_bb_gem_m7_strip_CM_GOOD);
   fChain->SetBranchAddress("bb.gem.m7.strip.CM_GOOD", bb_gem_m7_strip_CM_GOOD, &b_bb_gem_m7_strip_CM_GOOD);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.CorrCoeff", &Ndata_bb_gem_m7_strip_CorrCoeff, &b_Ndata_bb_gem_m7_strip_CorrCoeff);
   fChain->SetBranchAddress("bb.gem.m7.strip.CorrCoeff", bb_gem_m7_strip_CorrCoeff, &b_bb_gem_m7_strip_CorrCoeff);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.DeconvADCmax", &Ndata_bb_gem_m7_strip_DeconvADCmax, &b_Ndata_bb_gem_m7_strip_DeconvADCmax);
   fChain->SetBranchAddress("bb.gem.m7.strip.DeconvADCmax", bb_gem_m7_strip_DeconvADCmax, &b_bb_gem_m7_strip_DeconvADCmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.DeconvADCmaxCombo", &Ndata_bb_gem_m7_strip_DeconvADCmaxCombo, &b_Ndata_bb_gem_m7_strip_DeconvADCmaxCombo);
   fChain->SetBranchAddress("bb.gem.m7.strip.DeconvADCmaxCombo", bb_gem_m7_strip_DeconvADCmaxCombo, &b_bb_gem_m7_strip_DeconvADCmaxCombo);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.DeconvADCsamples", &Ndata_bb_gem_m7_strip_DeconvADCsamples, &b_Ndata_bb_gem_m7_strip_DeconvADCsamples);
   fChain->SetBranchAddress("bb.gem.m7.strip.DeconvADCsamples", bb_gem_m7_strip_DeconvADCsamples, &b_bb_gem_m7_strip_DeconvADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.DeconvADCsum", &Ndata_bb_gem_m7_strip_DeconvADCsum, &b_Ndata_bb_gem_m7_strip_DeconvADCsum);
   fChain->SetBranchAddress("bb.gem.m7.strip.DeconvADCsum", bb_gem_m7_strip_DeconvADCsum, &b_bb_gem_m7_strip_DeconvADCsum);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.ENABLE_CM", &Ndata_bb_gem_m7_strip_ENABLE_CM, &b_Ndata_bb_gem_m7_strip_ENABLE_CM);
   fChain->SetBranchAddress("bb.gem.m7.strip.ENABLE_CM", bb_gem_m7_strip_ENABLE_CM, &b_bb_gem_m7_strip_ENABLE_CM);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.IsU", &Ndata_bb_gem_m7_strip_IsU, &b_Ndata_bb_gem_m7_strip_IsU);
   fChain->SetBranchAddress("bb.gem.m7.strip.IsU", bb_gem_m7_strip_IsU, &b_bb_gem_m7_strip_IsU);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.IsV", &Ndata_bb_gem_m7_strip_IsV, &b_Ndata_bb_gem_m7_strip_IsV);
   fChain->SetBranchAddress("bb.gem.m7.strip.IsV", bb_gem_m7_strip_IsV, &b_bb_gem_m7_strip_IsV);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.TSchi2", &Ndata_bb_gem_m7_strip_TSchi2, &b_Ndata_bb_gem_m7_strip_TSchi2);
   fChain->SetBranchAddress("bb.gem.m7.strip.TSchi2", bb_gem_m7_strip_TSchi2, &b_bb_gem_m7_strip_TSchi2);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.TSprob", &Ndata_bb_gem_m7_strip_TSprob, &b_Ndata_bb_gem_m7_strip_TSprob);
   fChain->SetBranchAddress("bb.gem.m7.strip.TSprob", bb_gem_m7_strip_TSprob, &b_bb_gem_m7_strip_TSprob);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.Tcorr", &Ndata_bb_gem_m7_strip_Tcorr, &b_Ndata_bb_gem_m7_strip_Tcorr);
   fChain->SetBranchAddress("bb.gem.m7.strip.Tcorr", bb_gem_m7_strip_Tcorr, &b_bb_gem_m7_strip_Tcorr);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.Tdiff", &Ndata_bb_gem_m7_strip_Tdiff, &b_Ndata_bb_gem_m7_strip_Tdiff);
   fChain->SetBranchAddress("bb.gem.m7.strip.Tdiff", bb_gem_m7_strip_Tdiff, &b_bb_gem_m7_strip_Tdiff);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.Tfit", &Ndata_bb_gem_m7_strip_Tfit, &b_Ndata_bb_gem_m7_strip_Tfit);
   fChain->SetBranchAddress("bb.gem.m7.strip.Tfit", bb_gem_m7_strip_Tfit, &b_bb_gem_m7_strip_Tfit);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.Tmean", &Ndata_bb_gem_m7_strip_Tmean, &b_Ndata_bb_gem_m7_strip_Tmean);
   fChain->SetBranchAddress("bb.gem.m7.strip.Tmean", bb_gem_m7_strip_Tmean, &b_bb_gem_m7_strip_Tmean);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.TmeanDeconv", &Ndata_bb_gem_m7_strip_TmeanDeconv, &b_Ndata_bb_gem_m7_strip_TmeanDeconv);
   fChain->SetBranchAddress("bb.gem.m7.strip.TmeanDeconv", bb_gem_m7_strip_TmeanDeconv, &b_bb_gem_m7_strip_TmeanDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.Tsigma", &Ndata_bb_gem_m7_strip_Tsigma, &b_Ndata_bb_gem_m7_strip_Tsigma);
   fChain->SetBranchAddress("bb.gem.m7.strip.Tsigma", bb_gem_m7_strip_Tsigma, &b_bb_gem_m7_strip_Tsigma);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.adc_id", &Ndata_bb_gem_m7_strip_adc_id, &b_Ndata_bb_gem_m7_strip_adc_id);
   fChain->SetBranchAddress("bb.gem.m7.strip.adc_id", bb_gem_m7_strip_adc_id, &b_bb_gem_m7_strip_adc_id);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.crate", &Ndata_bb_gem_m7_strip_crate, &b_Ndata_bb_gem_m7_strip_crate);
   fChain->SetBranchAddress("bb.gem.m7.strip.crate", bb_gem_m7_strip_crate, &b_bb_gem_m7_strip_crate);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.isampmax", &Ndata_bb_gem_m7_strip_isampmax, &b_Ndata_bb_gem_m7_strip_isampmax);
   fChain->SetBranchAddress("bb.gem.m7.strip.isampmax", bb_gem_m7_strip_isampmax, &b_bb_gem_m7_strip_isampmax);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.isampmaxDeconv", &Ndata_bb_gem_m7_strip_isampmaxDeconv, &b_Ndata_bb_gem_m7_strip_isampmaxDeconv);
   fChain->SetBranchAddress("bb.gem.m7.strip.isampmaxDeconv", bb_gem_m7_strip_isampmaxDeconv, &b_bb_gem_m7_strip_isampmaxDeconv);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.isampmaxDeconvCombo", &Ndata_bb_gem_m7_strip_isampmaxDeconvCombo, &b_Ndata_bb_gem_m7_strip_isampmaxDeconvCombo);
   fChain->SetBranchAddress("bb.gem.m7.strip.isampmaxDeconvCombo", bb_gem_m7_strip_isampmaxDeconvCombo, &b_bb_gem_m7_strip_isampmaxDeconvCombo);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.isnegative", &Ndata_bb_gem_m7_strip_isnegative, &b_Ndata_bb_gem_m7_strip_isnegative);
   fChain->SetBranchAddress("bb.gem.m7.strip.isnegative", bb_gem_m7_strip_isnegative, &b_bb_gem_m7_strip_isnegative);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.isnegativeU", &Ndata_bb_gem_m7_strip_isnegativeU, &b_Ndata_bb_gem_m7_strip_isnegativeU);
   fChain->SetBranchAddress("bb.gem.m7.strip.isnegativeU", bb_gem_m7_strip_isnegativeU, &b_bb_gem_m7_strip_isnegativeU);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.isnegativeV", &Ndata_bb_gem_m7_strip_isnegativeV, &b_Ndata_bb_gem_m7_strip_isnegativeV);
   fChain->SetBranchAddress("bb.gem.m7.strip.isnegativeV", bb_gem_m7_strip_isnegativeV, &b_bb_gem_m7_strip_isnegativeV);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.isnegontrack", &Ndata_bb_gem_m7_strip_isnegontrack, &b_Ndata_bb_gem_m7_strip_isnegontrack);
   fChain->SetBranchAddress("bb.gem.m7.strip.isnegontrack", bb_gem_m7_strip_isnegontrack, &b_bb_gem_m7_strip_isnegontrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.isnegontrackU", &Ndata_bb_gem_m7_strip_isnegontrackU, &b_Ndata_bb_gem_m7_strip_isnegontrackU);
   fChain->SetBranchAddress("bb.gem.m7.strip.isnegontrackU", bb_gem_m7_strip_isnegontrackU, &b_bb_gem_m7_strip_isnegontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.isnegontrackV", &Ndata_bb_gem_m7_strip_isnegontrackV, &b_Ndata_bb_gem_m7_strip_isnegontrackV);
   fChain->SetBranchAddress("bb.gem.m7.strip.isnegontrackV", bb_gem_m7_strip_isnegontrackV, &b_bb_gem_m7_strip_isnegontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.istrip", &Ndata_bb_gem_m7_strip_istrip, &b_Ndata_bb_gem_m7_strip_istrip);
   fChain->SetBranchAddress("bb.gem.m7.strip.istrip", bb_gem_m7_strip_istrip, &b_bb_gem_m7_strip_istrip);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.itrack", &Ndata_bb_gem_m7_strip_itrack, &b_Ndata_bb_gem_m7_strip_itrack);
   fChain->SetBranchAddress("bb.gem.m7.strip.itrack", bb_gem_m7_strip_itrack, &b_bb_gem_m7_strip_itrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.mpd", &Ndata_bb_gem_m7_strip_mpd, &b_Ndata_bb_gem_m7_strip_mpd);
   fChain->SetBranchAddress("bb.gem.m7.strip.mpd", bb_gem_m7_strip_mpd, &b_bb_gem_m7_strip_mpd);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.ontrack", &Ndata_bb_gem_m7_strip_ontrack, &b_Ndata_bb_gem_m7_strip_ontrack);
   fChain->SetBranchAddress("bb.gem.m7.strip.ontrack", bb_gem_m7_strip_ontrack, &b_bb_gem_m7_strip_ontrack);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.ontrackU", &Ndata_bb_gem_m7_strip_ontrackU, &b_Ndata_bb_gem_m7_strip_ontrackU);
   fChain->SetBranchAddress("bb.gem.m7.strip.ontrackU", bb_gem_m7_strip_ontrackU, &b_bb_gem_m7_strip_ontrackU);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.ontrackV", &Ndata_bb_gem_m7_strip_ontrackV, &b_Ndata_bb_gem_m7_strip_ontrackV);
   fChain->SetBranchAddress("bb.gem.m7.strip.ontrackV", bb_gem_m7_strip_ontrackV, &b_bb_gem_m7_strip_ontrackV);
   fChain->SetBranchAddress("Ndata.bb.gem.m7.strip.rawADCsamples", &Ndata_bb_gem_m7_strip_rawADCsamples, &b_Ndata_bb_gem_m7_strip_rawADCsamples);
   fChain->SetBranchAddress("bb.gem.m7.strip.rawADCsamples", bb_gem_m7_strip_rawADCsamples, &b_bb_gem_m7_strip_rawADCsamples);
   fChain->SetBranchAddress("Ndata.bb.gem.n2Dhit_layer", &Ndata_bb_gem_n2Dhit_layer, &b_Ndata_bb_gem_n2Dhit_layer);
   fChain->SetBranchAddress("bb.gem.n2Dhit_layer", bb_gem_n2Dhit_layer, &b_bb_gem_n2Dhit_layer);
   fChain->SetBranchAddress("Ndata.bb.gem.nclustu_layer", &Ndata_bb_gem_nclustu_layer, &b_Ndata_bb_gem_nclustu_layer);
   fChain->SetBranchAddress("bb.gem.nclustu_layer", bb_gem_nclustu_layer, &b_bb_gem_nclustu_layer);
   fChain->SetBranchAddress("Ndata.bb.gem.nclustv_layer", &Ndata_bb_gem_nclustv_layer, &b_Ndata_bb_gem_nclustv_layer);
   fChain->SetBranchAddress("bb.gem.nclustv_layer", bb_gem_nclustv_layer, &b_bb_gem_nclustv_layer);
   fChain->SetBranchAddress("Ndata.bb.gem.nstripsu_layer", &Ndata_bb_gem_nstripsu_layer, &b_Ndata_bb_gem_nstripsu_layer);
   fChain->SetBranchAddress("bb.gem.nstripsu_layer", bb_gem_nstripsu_layer, &b_bb_gem_nstripsu_layer);
   fChain->SetBranchAddress("Ndata.bb.gem.nstripsv_layer", &Ndata_bb_gem_nstripsv_layer, &b_Ndata_bb_gem_nstripsv_layer);
   fChain->SetBranchAddress("bb.gem.nstripsv_layer", bb_gem_nstripsv_layer, &b_bb_gem_nstripsv_layer);
   fChain->SetBranchAddress("Ndata.bb.gem.track.chi2ndf", &Ndata_bb_gem_track_chi2ndf, &b_Ndata_bb_gem_track_chi2ndf);
   fChain->SetBranchAddress("bb.gem.track.chi2ndf", bb_gem_track_chi2ndf, &b_bb_gem_track_chi2ndf);
   fChain->SetBranchAddress("Ndata.bb.gem.track.chi2ndf_hitquality", &Ndata_bb_gem_track_chi2ndf_hitquality, &b_Ndata_bb_gem_track_chi2ndf_hitquality);
   fChain->SetBranchAddress("bb.gem.track.chi2ndf_hitquality", bb_gem_track_chi2ndf_hitquality, &b_bb_gem_track_chi2ndf_hitquality);
   fChain->SetBranchAddress("Ndata.bb.gem.track.ngoodhits", &Ndata_bb_gem_track_ngoodhits, &b_Ndata_bb_gem_track_ngoodhits);
   fChain->SetBranchAddress("bb.gem.track.ngoodhits", bb_gem_track_ngoodhits, &b_bb_gem_track_ngoodhits);
   fChain->SetBranchAddress("Ndata.bb.gem.track.nhits", &Ndata_bb_gem_track_nhits, &b_Ndata_bb_gem_track_nhits);
   fChain->SetBranchAddress("bb.gem.track.nhits", bb_gem_track_nhits, &b_bb_gem_track_nhits);
   fChain->SetBranchAddress("Ndata.bb.gem.track.t0", &Ndata_bb_gem_track_t0, &b_Ndata_bb_gem_track_t0);
   fChain->SetBranchAddress("bb.gem.track.t0", bb_gem_track_t0, &b_bb_gem_track_t0);
   fChain->SetBranchAddress("Ndata.bb.gem.track.x", &Ndata_bb_gem_track_x, &b_Ndata_bb_gem_track_x);
   fChain->SetBranchAddress("bb.gem.track.x", bb_gem_track_x, &b_bb_gem_track_x);
   fChain->SetBranchAddress("Ndata.bb.gem.track.xp", &Ndata_bb_gem_track_xp, &b_Ndata_bb_gem_track_xp);
   fChain->SetBranchAddress("bb.gem.track.xp", bb_gem_track_xp, &b_bb_gem_track_xp);
   fChain->SetBranchAddress("Ndata.bb.gem.track.y", &Ndata_bb_gem_track_y, &b_Ndata_bb_gem_track_y);
   fChain->SetBranchAddress("bb.gem.track.y", bb_gem_track_y, &b_bb_gem_track_y);
   fChain->SetBranchAddress("Ndata.bb.gem.track.yp", &Ndata_bb_gem_track_yp, &b_Ndata_bb_gem_track_yp);
   fChain->SetBranchAddress("bb.gem.track.yp", bb_gem_track_yp, &b_bb_gem_track_yp);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.allclus.adc", &Ndata_bb_grinch_tdc_allclus_adc, &b_Ndata_bb_grinch_tdc_allclus_adc);
   fChain->SetBranchAddress("bb.grinch_tdc.allclus.adc", bb_grinch_tdc_allclus_adc, &b_bb_grinch_tdc_allclus_adc);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.allclus.dx", &Ndata_bb_grinch_tdc_allclus_dx, &b_Ndata_bb_grinch_tdc_allclus_dx);
   fChain->SetBranchAddress("bb.grinch_tdc.allclus.dx", bb_grinch_tdc_allclus_dx, &b_bb_grinch_tdc_allclus_dx);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.allclus.dy", &Ndata_bb_grinch_tdc_allclus_dy, &b_Ndata_bb_grinch_tdc_allclus_dy);
   fChain->SetBranchAddress("bb.grinch_tdc.allclus.dy", bb_grinch_tdc_allclus_dy, &b_bb_grinch_tdc_allclus_dy);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.allclus.mirrorindex", &Ndata_bb_grinch_tdc_allclus_mirrorindex, &b_Ndata_bb_grinch_tdc_allclus_mirrorindex);
   fChain->SetBranchAddress("bb.grinch_tdc.allclus.mirrorindex", bb_grinch_tdc_allclus_mirrorindex, &b_bb_grinch_tdc_allclus_mirrorindex);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.allclus.size", &Ndata_bb_grinch_tdc_allclus_size, &b_Ndata_bb_grinch_tdc_allclus_size);
   fChain->SetBranchAddress("bb.grinch_tdc.allclus.size", bb_grinch_tdc_allclus_size, &b_bb_grinch_tdc_allclus_size);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.allclus.t_mean", &Ndata_bb_grinch_tdc_allclus_t_mean, &b_Ndata_bb_grinch_tdc_allclus_t_mean);
   fChain->SetBranchAddress("bb.grinch_tdc.allclus.t_mean", bb_grinch_tdc_allclus_t_mean, &b_bb_grinch_tdc_allclus_t_mean);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.allclus.t_rms", &Ndata_bb_grinch_tdc_allclus_t_rms, &b_Ndata_bb_grinch_tdc_allclus_t_rms);
   fChain->SetBranchAddress("bb.grinch_tdc.allclus.t_rms", bb_grinch_tdc_allclus_t_rms, &b_bb_grinch_tdc_allclus_t_rms);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.allclus.tot_mean", &Ndata_bb_grinch_tdc_allclus_tot_mean, &b_Ndata_bb_grinch_tdc_allclus_tot_mean);
   fChain->SetBranchAddress("bb.grinch_tdc.allclus.tot_mean", bb_grinch_tdc_allclus_tot_mean, &b_bb_grinch_tdc_allclus_tot_mean);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.allclus.trackindex", &Ndata_bb_grinch_tdc_allclus_trackindex, &b_Ndata_bb_grinch_tdc_allclus_trackindex);
   fChain->SetBranchAddress("bb.grinch_tdc.allclus.trackindex", bb_grinch_tdc_allclus_trackindex, &b_bb_grinch_tdc_allclus_trackindex);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.allclus.x_mean", &Ndata_bb_grinch_tdc_allclus_x_mean, &b_Ndata_bb_grinch_tdc_allclus_x_mean);
   fChain->SetBranchAddress("bb.grinch_tdc.allclus.x_mean", bb_grinch_tdc_allclus_x_mean, &b_bb_grinch_tdc_allclus_x_mean);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.allclus.y_mean", &Ndata_bb_grinch_tdc_allclus_y_mean, &b_Ndata_bb_grinch_tdc_allclus_y_mean);
   fChain->SetBranchAddress("bb.grinch_tdc.allclus.y_mean", bb_grinch_tdc_allclus_y_mean, &b_bb_grinch_tdc_allclus_y_mean);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.hit.amp", &Ndata_bb_grinch_tdc_hit_amp, &b_Ndata_bb_grinch_tdc_hit_amp);
   fChain->SetBranchAddress("bb.grinch_tdc.hit.amp", bb_grinch_tdc_hit_amp, &b_bb_grinch_tdc_hit_amp);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.hit.clustindex", &Ndata_bb_grinch_tdc_hit_clustindex, &b_Ndata_bb_grinch_tdc_hit_clustindex);
   fChain->SetBranchAddress("bb.grinch_tdc.hit.clustindex", bb_grinch_tdc_hit_clustindex, &b_bb_grinch_tdc_hit_clustindex);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.hit.col", &Ndata_bb_grinch_tdc_hit_col, &b_Ndata_bb_grinch_tdc_hit_col);
   fChain->SetBranchAddress("bb.grinch_tdc.hit.col", bb_grinch_tdc_hit_col, &b_bb_grinch_tdc_hit_col);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.hit.pmtnum", &Ndata_bb_grinch_tdc_hit_pmtnum, &b_Ndata_bb_grinch_tdc_hit_pmtnum);
   fChain->SetBranchAddress("bb.grinch_tdc.hit.pmtnum", bb_grinch_tdc_hit_pmtnum, &b_bb_grinch_tdc_hit_pmtnum);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.hit.row", &Ndata_bb_grinch_tdc_hit_row, &b_Ndata_bb_grinch_tdc_hit_row);
   fChain->SetBranchAddress("bb.grinch_tdc.hit.row", bb_grinch_tdc_hit_row, &b_bb_grinch_tdc_hit_row);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.hit.time", &Ndata_bb_grinch_tdc_hit_time, &b_Ndata_bb_grinch_tdc_hit_time);
   fChain->SetBranchAddress("bb.grinch_tdc.hit.time", bb_grinch_tdc_hit_time, &b_bb_grinch_tdc_hit_time);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.hit.trackindex", &Ndata_bb_grinch_tdc_hit_trackindex, &b_Ndata_bb_grinch_tdc_hit_trackindex);
   fChain->SetBranchAddress("bb.grinch_tdc.hit.trackindex", bb_grinch_tdc_hit_trackindex, &b_bb_grinch_tdc_hit_trackindex);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.hit.xhit", &Ndata_bb_grinch_tdc_hit_xhit, &b_Ndata_bb_grinch_tdc_hit_xhit);
   fChain->SetBranchAddress("bb.grinch_tdc.hit.xhit", bb_grinch_tdc_hit_xhit, &b_bb_grinch_tdc_hit_xhit);
   fChain->SetBranchAddress("Ndata.bb.grinch_tdc.hit.yhit", &Ndata_bb_grinch_tdc_hit_yhit, &b_Ndata_bb_grinch_tdc_hit_yhit);
   fChain->SetBranchAddress("bb.grinch_tdc.hit.yhit", bb_grinch_tdc_hit_yhit, &b_bb_grinch_tdc_hit_yhit);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.Ref.tdc", &Ndata_bb_hodotdc_Ref_tdc, &b_Ndata_bb_hodotdc_Ref_tdc);
   fChain->SetBranchAddress("bb.hodotdc.Ref.tdc", bb_hodotdc_Ref_tdc, &b_bb_hodotdc_Ref_tdc);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.Ref.tdc_mult", &Ndata_bb_hodotdc_Ref_tdc_mult, &b_Ndata_bb_hodotdc_Ref_tdc_mult);
   fChain->SetBranchAddress("bb.hodotdc.Ref.tdc_mult", bb_hodotdc_Ref_tdc_mult, &b_bb_hodotdc_Ref_tdc_mult);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.Ref.tdc_te", &Ndata_bb_hodotdc_Ref_tdc_te, &b_Ndata_bb_hodotdc_Ref_tdc_te);
   fChain->SetBranchAddress("bb.hodotdc.Ref.tdc_te", bb_hodotdc_Ref_tdc_te, &b_bb_hodotdc_Ref_tdc_te);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.Ref.tdc_tot", &Ndata_bb_hodotdc_Ref_tdc_tot, &b_Ndata_bb_hodotdc_Ref_tdc_tot);
   fChain->SetBranchAddress("bb.hodotdc.Ref.tdc_tot", bb_hodotdc_Ref_tdc_tot, &b_bb_hodotdc_Ref_tdc_tot);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.Ref.tdcelemID", &Ndata_bb_hodotdc_Ref_tdcelemID, &b_Ndata_bb_hodotdc_Ref_tdcelemID);
   fChain->SetBranchAddress("bb.hodotdc.Ref.tdcelemID", bb_hodotdc_Ref_tdcelemID, &b_bb_hodotdc_Ref_tdcelemID);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.bar.tdc.id", &Ndata_bb_hodotdc_clus_bar_tdc_id, &b_Ndata_bb_hodotdc_clus_bar_tdc_id);
   fChain->SetBranchAddress("bb.hodotdc.clus.bar.tdc.id", bb_hodotdc_clus_bar_tdc_id, &b_bb_hodotdc_clus_bar_tdc_id);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.bar.tdc.itrack", &Ndata_bb_hodotdc_clus_bar_tdc_itrack, &b_Ndata_bb_hodotdc_clus_bar_tdc_itrack);
   fChain->SetBranchAddress("bb.hodotdc.clus.bar.tdc.itrack", bb_hodotdc_clus_bar_tdc_itrack, &b_bb_hodotdc_clus_bar_tdc_itrack);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.bar.tdc.meantime", &Ndata_bb_hodotdc_clus_bar_tdc_meantime, &b_Ndata_bb_hodotdc_clus_bar_tdc_meantime);
   fChain->SetBranchAddress("bb.hodotdc.clus.bar.tdc.meantime", bb_hodotdc_clus_bar_tdc_meantime, &b_bb_hodotdc_clus_bar_tdc_meantime);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.bar.tdc.meantot", &Ndata_bb_hodotdc_clus_bar_tdc_meantot, &b_Ndata_bb_hodotdc_clus_bar_tdc_meantot);
   fChain->SetBranchAddress("bb.hodotdc.clus.bar.tdc.meantot", bb_hodotdc_clus_bar_tdc_meantot, &b_bb_hodotdc_clus_bar_tdc_meantot);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.bar.tdc.timediff", &Ndata_bb_hodotdc_clus_bar_tdc_timediff, &b_Ndata_bb_hodotdc_clus_bar_tdc_timediff);
   fChain->SetBranchAddress("bb.hodotdc.clus.bar.tdc.timediff", bb_hodotdc_clus_bar_tdc_timediff, &b_bb_hodotdc_clus_bar_tdc_timediff);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.bar.tdc.timehitpos", &Ndata_bb_hodotdc_clus_bar_tdc_timehitpos, &b_Ndata_bb_hodotdc_clus_bar_tdc_timehitpos);
   fChain->SetBranchAddress("bb.hodotdc.clus.bar.tdc.timehitpos", bb_hodotdc_clus_bar_tdc_timehitpos, &b_bb_hodotdc_clus_bar_tdc_timehitpos);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.bar.tdc.tleft", &Ndata_bb_hodotdc_clus_bar_tdc_tleft, &b_Ndata_bb_hodotdc_clus_bar_tdc_tleft);
   fChain->SetBranchAddress("bb.hodotdc.clus.bar.tdc.tleft", bb_hodotdc_clus_bar_tdc_tleft, &b_bb_hodotdc_clus_bar_tdc_tleft);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.bar.tdc.totleft", &Ndata_bb_hodotdc_clus_bar_tdc_totleft, &b_Ndata_bb_hodotdc_clus_bar_tdc_totleft);
   fChain->SetBranchAddress("bb.hodotdc.clus.bar.tdc.totleft", bb_hodotdc_clus_bar_tdc_totleft, &b_bb_hodotdc_clus_bar_tdc_totleft);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.bar.tdc.totright", &Ndata_bb_hodotdc_clus_bar_tdc_totright, &b_Ndata_bb_hodotdc_clus_bar_tdc_totright);
   fChain->SetBranchAddress("bb.hodotdc.clus.bar.tdc.totright", bb_hodotdc_clus_bar_tdc_totright, &b_bb_hodotdc_clus_bar_tdc_totright);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.bar.tdc.tright", &Ndata_bb_hodotdc_clus_bar_tdc_tright, &b_Ndata_bb_hodotdc_clus_bar_tdc_tright);
   fChain->SetBranchAddress("bb.hodotdc.clus.bar.tdc.tright", bb_hodotdc_clus_bar_tdc_tright, &b_bb_hodotdc_clus_bar_tdc_tright);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.bar.tdc.vpos", &Ndata_bb_hodotdc_clus_bar_tdc_vpos, &b_Ndata_bb_hodotdc_clus_bar_tdc_vpos);
   fChain->SetBranchAddress("bb.hodotdc.clus.bar.tdc.vpos", bb_hodotdc_clus_bar_tdc_vpos, &b_bb_hodotdc_clus_bar_tdc_vpos);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.id", &Ndata_bb_hodotdc_clus_id, &b_Ndata_bb_hodotdc_clus_id);
   fChain->SetBranchAddress("bb.hodotdc.clus.id", bb_hodotdc_clus_id, &b_bb_hodotdc_clus_id);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.size", &Ndata_bb_hodotdc_clus_size, &b_Ndata_bb_hodotdc_clus_size);
   fChain->SetBranchAddress("bb.hodotdc.clus.size", bb_hodotdc_clus_size, &b_bb_hodotdc_clus_size);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.tdiff", &Ndata_bb_hodotdc_clus_tdiff, &b_Ndata_bb_hodotdc_clus_tdiff);
   fChain->SetBranchAddress("bb.hodotdc.clus.tdiff", bb_hodotdc_clus_tdiff, &b_bb_hodotdc_clus_tdiff);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.tleft", &Ndata_bb_hodotdc_clus_tleft, &b_Ndata_bb_hodotdc_clus_tleft);
   fChain->SetBranchAddress("bb.hodotdc.clus.tleft", bb_hodotdc_clus_tleft, &b_bb_hodotdc_clus_tleft);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.tmean", &Ndata_bb_hodotdc_clus_tmean, &b_Ndata_bb_hodotdc_clus_tmean);
   fChain->SetBranchAddress("bb.hodotdc.clus.tmean", bb_hodotdc_clus_tmean, &b_bb_hodotdc_clus_tmean);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.totleft", &Ndata_bb_hodotdc_clus_totleft, &b_Ndata_bb_hodotdc_clus_totleft);
   fChain->SetBranchAddress("bb.hodotdc.clus.totleft", bb_hodotdc_clus_totleft, &b_bb_hodotdc_clus_totleft);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.totmean", &Ndata_bb_hodotdc_clus_totmean, &b_Ndata_bb_hodotdc_clus_totmean);
   fChain->SetBranchAddress("bb.hodotdc.clus.totmean", bb_hodotdc_clus_totmean, &b_bb_hodotdc_clus_totmean);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.totright", &Ndata_bb_hodotdc_clus_totright, &b_Ndata_bb_hodotdc_clus_totright);
   fChain->SetBranchAddress("bb.hodotdc.clus.totright", bb_hodotdc_clus_totright, &b_bb_hodotdc_clus_totright);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.trackindex", &Ndata_bb_hodotdc_clus_trackindex, &b_Ndata_bb_hodotdc_clus_trackindex);
   fChain->SetBranchAddress("bb.hodotdc.clus.trackindex", bb_hodotdc_clus_trackindex, &b_bb_hodotdc_clus_trackindex);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.tright", &Ndata_bb_hodotdc_clus_tright, &b_Ndata_bb_hodotdc_clus_tright);
   fChain->SetBranchAddress("bb.hodotdc.clus.tright", bb_hodotdc_clus_tright, &b_bb_hodotdc_clus_tright);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.xmean", &Ndata_bb_hodotdc_clus_xmean, &b_Ndata_bb_hodotdc_clus_xmean);
   fChain->SetBranchAddress("bb.hodotdc.clus.xmean", bb_hodotdc_clus_xmean, &b_bb_hodotdc_clus_xmean);
   fChain->SetBranchAddress("Ndata.bb.hodotdc.clus.ymean", &Ndata_bb_hodotdc_clus_ymean, &b_Ndata_bb_hodotdc_clus_ymean);
   fChain->SetBranchAddress("bb.hodotdc.clus.ymean", bb_hodotdc_clus_ymean, &b_bb_hodotdc_clus_ymean);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.again", &Ndata_bb_ps_clus_again, &b_Ndata_bb_ps_clus_again);
   fChain->SetBranchAddress("bb.ps.clus.again", bb_ps_clus_again, &b_bb_ps_clus_again);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.atime", &Ndata_bb_ps_clus_atime, &b_Ndata_bb_ps_clus_atime);
   fChain->SetBranchAddress("bb.ps.clus.atime", bb_ps_clus_atime, &b_bb_ps_clus_atime);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.col", &Ndata_bb_ps_clus_col, &b_Ndata_bb_ps_clus_col);
   fChain->SetBranchAddress("bb.ps.clus.col", bb_ps_clus_col, &b_bb_ps_clus_col);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.e", &Ndata_bb_ps_clus_e, &b_Ndata_bb_ps_clus_e);
   fChain->SetBranchAddress("bb.ps.clus.e", bb_ps_clus_e, &b_bb_ps_clus_e);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.eblk", &Ndata_bb_ps_clus_eblk, &b_Ndata_bb_ps_clus_eblk);
   fChain->SetBranchAddress("bb.ps.clus.eblk", bb_ps_clus_eblk, &b_bb_ps_clus_eblk);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.id", &Ndata_bb_ps_clus_id, &b_Ndata_bb_ps_clus_id);
   fChain->SetBranchAddress("bb.ps.clus.id", bb_ps_clus_id, &b_bb_ps_clus_id);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.nblk", &Ndata_bb_ps_clus_nblk, &b_Ndata_bb_ps_clus_nblk);
   fChain->SetBranchAddress("bb.ps.clus.nblk", bb_ps_clus_nblk, &b_bb_ps_clus_nblk);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.row", &Ndata_bb_ps_clus_row, &b_Ndata_bb_ps_clus_row);
   fChain->SetBranchAddress("bb.ps.clus.row", bb_ps_clus_row, &b_bb_ps_clus_row);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.tdctime", &Ndata_bb_ps_clus_tdctime, &b_Ndata_bb_ps_clus_tdctime);
   fChain->SetBranchAddress("bb.ps.clus.tdctime", bb_ps_clus_tdctime, &b_bb_ps_clus_tdctime);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.x", &Ndata_bb_ps_clus_x, &b_Ndata_bb_ps_clus_x);
   fChain->SetBranchAddress("bb.ps.clus.x", bb_ps_clus_x, &b_bb_ps_clus_x);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.y", &Ndata_bb_ps_clus_y, &b_Ndata_bb_ps_clus_y);
   fChain->SetBranchAddress("bb.ps.clus.y", bb_ps_clus_y, &b_bb_ps_clus_y);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.again", &Ndata_bb_ps_clus_blk_again, &b_Ndata_bb_ps_clus_blk_again);
   fChain->SetBranchAddress("bb.ps.clus_blk.again", bb_ps_clus_blk_again, &b_bb_ps_clus_blk_again);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.atime", &Ndata_bb_ps_clus_blk_atime, &b_Ndata_bb_ps_clus_blk_atime);
   fChain->SetBranchAddress("bb.ps.clus_blk.atime", bb_ps_clus_blk_atime, &b_bb_ps_clus_blk_atime);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.col", &Ndata_bb_ps_clus_blk_col, &b_Ndata_bb_ps_clus_blk_col);
   fChain->SetBranchAddress("bb.ps.clus_blk.col", bb_ps_clus_blk_col, &b_bb_ps_clus_blk_col);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.e", &Ndata_bb_ps_clus_blk_e, &b_Ndata_bb_ps_clus_blk_e);
   fChain->SetBranchAddress("bb.ps.clus_blk.e", bb_ps_clus_blk_e, &b_bb_ps_clus_blk_e);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.id", &Ndata_bb_ps_clus_blk_id, &b_Ndata_bb_ps_clus_blk_id);
   fChain->SetBranchAddress("bb.ps.clus_blk.id", bb_ps_clus_blk_id, &b_bb_ps_clus_blk_id);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.row", &Ndata_bb_ps_clus_blk_row, &b_Ndata_bb_ps_clus_blk_row);
   fChain->SetBranchAddress("bb.ps.clus_blk.row", bb_ps_clus_blk_row, &b_bb_ps_clus_blk_row);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.tdctime", &Ndata_bb_ps_clus_blk_tdctime, &b_Ndata_bb_ps_clus_blk_tdctime);
   fChain->SetBranchAddress("bb.ps.clus_blk.tdctime", bb_ps_clus_blk_tdctime, &b_bb_ps_clus_blk_tdctime);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.x", &Ndata_bb_ps_clus_blk_x, &b_Ndata_bb_ps_clus_blk_x);
   fChain->SetBranchAddress("bb.ps.clus_blk.x", bb_ps_clus_blk_x, &b_bb_ps_clus_blk_x);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.y", &Ndata_bb_ps_clus_blk_y, &b_Ndata_bb_ps_clus_blk_y);
   fChain->SetBranchAddress("bb.ps.clus_blk.y", bb_ps_clus_blk_y, &b_bb_ps_clus_blk_y);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.again", &Ndata_bb_sh_clus_again, &b_Ndata_bb_sh_clus_again);
   fChain->SetBranchAddress("bb.sh.clus.again", bb_sh_clus_again, &b_bb_sh_clus_again);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.atime", &Ndata_bb_sh_clus_atime, &b_Ndata_bb_sh_clus_atime);
   fChain->SetBranchAddress("bb.sh.clus.atime", bb_sh_clus_atime, &b_bb_sh_clus_atime);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.col", &Ndata_bb_sh_clus_col, &b_Ndata_bb_sh_clus_col);
   fChain->SetBranchAddress("bb.sh.clus.col", bb_sh_clus_col, &b_bb_sh_clus_col);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.e", &Ndata_bb_sh_clus_e, &b_Ndata_bb_sh_clus_e);
   fChain->SetBranchAddress("bb.sh.clus.e", bb_sh_clus_e, &b_bb_sh_clus_e);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.eblk", &Ndata_bb_sh_clus_eblk, &b_Ndata_bb_sh_clus_eblk);
   fChain->SetBranchAddress("bb.sh.clus.eblk", bb_sh_clus_eblk, &b_bb_sh_clus_eblk);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.id", &Ndata_bb_sh_clus_id, &b_Ndata_bb_sh_clus_id);
   fChain->SetBranchAddress("bb.sh.clus.id", bb_sh_clus_id, &b_bb_sh_clus_id);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.nblk", &Ndata_bb_sh_clus_nblk, &b_Ndata_bb_sh_clus_nblk);
   fChain->SetBranchAddress("bb.sh.clus.nblk", bb_sh_clus_nblk, &b_bb_sh_clus_nblk);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.row", &Ndata_bb_sh_clus_row, &b_Ndata_bb_sh_clus_row);
   fChain->SetBranchAddress("bb.sh.clus.row", bb_sh_clus_row, &b_bb_sh_clus_row);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.tdctime", &Ndata_bb_sh_clus_tdctime, &b_Ndata_bb_sh_clus_tdctime);
   fChain->SetBranchAddress("bb.sh.clus.tdctime", bb_sh_clus_tdctime, &b_bb_sh_clus_tdctime);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.x", &Ndata_bb_sh_clus_x, &b_Ndata_bb_sh_clus_x);
   fChain->SetBranchAddress("bb.sh.clus.x", bb_sh_clus_x, &b_bb_sh_clus_x);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.y", &Ndata_bb_sh_clus_y, &b_Ndata_bb_sh_clus_y);
   fChain->SetBranchAddress("bb.sh.clus.y", bb_sh_clus_y, &b_bb_sh_clus_y);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.again", &Ndata_bb_sh_clus_blk_again, &b_Ndata_bb_sh_clus_blk_again);
   fChain->SetBranchAddress("bb.sh.clus_blk.again", bb_sh_clus_blk_again, &b_bb_sh_clus_blk_again);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.atime", &Ndata_bb_sh_clus_blk_atime, &b_Ndata_bb_sh_clus_blk_atime);
   fChain->SetBranchAddress("bb.sh.clus_blk.atime", bb_sh_clus_blk_atime, &b_bb_sh_clus_blk_atime);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.col", &Ndata_bb_sh_clus_blk_col, &b_Ndata_bb_sh_clus_blk_col);
   fChain->SetBranchAddress("bb.sh.clus_blk.col", bb_sh_clus_blk_col, &b_bb_sh_clus_blk_col);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.e", &Ndata_bb_sh_clus_blk_e, &b_Ndata_bb_sh_clus_blk_e);
   fChain->SetBranchAddress("bb.sh.clus_blk.e", bb_sh_clus_blk_e, &b_bb_sh_clus_blk_e);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.id", &Ndata_bb_sh_clus_blk_id, &b_Ndata_bb_sh_clus_blk_id);
   fChain->SetBranchAddress("bb.sh.clus_blk.id", bb_sh_clus_blk_id, &b_bb_sh_clus_blk_id);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.row", &Ndata_bb_sh_clus_blk_row, &b_Ndata_bb_sh_clus_blk_row);
   fChain->SetBranchAddress("bb.sh.clus_blk.row", bb_sh_clus_blk_row, &b_bb_sh_clus_blk_row);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.tdctime", &Ndata_bb_sh_clus_blk_tdctime, &b_Ndata_bb_sh_clus_blk_tdctime);
   fChain->SetBranchAddress("bb.sh.clus_blk.tdctime", bb_sh_clus_blk_tdctime, &b_bb_sh_clus_blk_tdctime);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.x", &Ndata_bb_sh_clus_blk_x, &b_Ndata_bb_sh_clus_blk_x);
   fChain->SetBranchAddress("bb.sh.clus_blk.x", bb_sh_clus_blk_x, &b_bb_sh_clus_blk_x);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.y", &Ndata_bb_sh_clus_blk_y, &b_Ndata_bb_sh_clus_blk_y);
   fChain->SetBranchAddress("bb.sh.clus_blk.y", bb_sh_clus_blk_y, &b_bb_sh_clus_blk_y);
   fChain->SetBranchAddress("Ndata.bb.tdctrig.tdc", &Ndata_bb_tdctrig_tdc, &b_Ndata_bb_tdctrig_tdc);
   fChain->SetBranchAddress("bb.tdctrig.tdc", bb_tdctrig_tdc, &b_bb_tdctrig_tdc);
   fChain->SetBranchAddress("Ndata.bb.tdctrig.tdcelemID", &Ndata_bb_tdctrig_tdcelemID, &b_Ndata_bb_tdctrig_tdcelemID);
   fChain->SetBranchAddress("bb.tdctrig.tdcelemID", bb_tdctrig_tdcelemID, &b_bb_tdctrig_tdcelemID);
   fChain->SetBranchAddress("Ndata.bb.tr.beta", &Ndata_bb_tr_beta, &b_Ndata_bb_tr_beta);
   fChain->SetBranchAddress("bb.tr.beta", bb_tr_beta, &b_bb_tr_beta);
   fChain->SetBranchAddress("Ndata.bb.tr.chi2", &Ndata_bb_tr_chi2, &b_Ndata_bb_tr_chi2);
   fChain->SetBranchAddress("bb.tr.chi2", bb_tr_chi2, &b_bb_tr_chi2);
   fChain->SetBranchAddress("Ndata.bb.tr.d_ph", &Ndata_bb_tr_d_ph, &b_Ndata_bb_tr_d_ph);
   fChain->SetBranchAddress("bb.tr.d_ph", bb_tr_d_ph, &b_bb_tr_d_ph);
   fChain->SetBranchAddress("Ndata.bb.tr.d_th", &Ndata_bb_tr_d_th, &b_Ndata_bb_tr_d_th);
   fChain->SetBranchAddress("bb.tr.d_th", bb_tr_d_th, &b_bb_tr_d_th);
   fChain->SetBranchAddress("Ndata.bb.tr.d_x", &Ndata_bb_tr_d_x, &b_Ndata_bb_tr_d_x);
   fChain->SetBranchAddress("bb.tr.d_x", bb_tr_d_x, &b_bb_tr_d_x);
   fChain->SetBranchAddress("Ndata.bb.tr.d_y", &Ndata_bb_tr_d_y, &b_Ndata_bb_tr_d_y);
   fChain->SetBranchAddress("bb.tr.d_y", bb_tr_d_y, &b_bb_tr_d_y);
   fChain->SetBranchAddress("Ndata.bb.tr.dbeta", &Ndata_bb_tr_dbeta, &b_Ndata_bb_tr_dbeta);
   fChain->SetBranchAddress("bb.tr.dbeta", bb_tr_dbeta, &b_bb_tr_dbeta);
   fChain->SetBranchAddress("Ndata.bb.tr.dtime", &Ndata_bb_tr_dtime, &b_Ndata_bb_tr_dtime);
   fChain->SetBranchAddress("bb.tr.dtime", bb_tr_dtime, &b_bb_tr_dtime);
   fChain->SetBranchAddress("Ndata.bb.tr.flag", &Ndata_bb_tr_flag, &b_Ndata_bb_tr_flag);
   fChain->SetBranchAddress("bb.tr.flag", bb_tr_flag, &b_bb_tr_flag);
   fChain->SetBranchAddress("Ndata.bb.tr.ndof", &Ndata_bb_tr_ndof, &b_Ndata_bb_tr_ndof);
   fChain->SetBranchAddress("bb.tr.ndof", bb_tr_ndof, &b_bb_tr_ndof);
   fChain->SetBranchAddress("Ndata.bb.tr.p", &Ndata_bb_tr_p, &b_Ndata_bb_tr_p);
   fChain->SetBranchAddress("bb.tr.p", bb_tr_p, &b_bb_tr_p);
   fChain->SetBranchAddress("Ndata.bb.tr.pathl", &Ndata_bb_tr_pathl, &b_Ndata_bb_tr_pathl);
   fChain->SetBranchAddress("bb.tr.pathl", bb_tr_pathl, &b_bb_tr_pathl);
   fChain->SetBranchAddress("Ndata.bb.tr.ph", &Ndata_bb_tr_ph, &b_Ndata_bb_tr_ph);
   fChain->SetBranchAddress("bb.tr.ph", bb_tr_ph, &b_bb_tr_ph);
   fChain->SetBranchAddress("Ndata.bb.tr.px", &Ndata_bb_tr_px, &b_Ndata_bb_tr_px);
   fChain->SetBranchAddress("bb.tr.px", bb_tr_px, &b_bb_tr_px);
   fChain->SetBranchAddress("Ndata.bb.tr.py", &Ndata_bb_tr_py, &b_Ndata_bb_tr_py);
   fChain->SetBranchAddress("bb.tr.py", bb_tr_py, &b_bb_tr_py);
   fChain->SetBranchAddress("Ndata.bb.tr.pz", &Ndata_bb_tr_pz, &b_Ndata_bb_tr_pz);
   fChain->SetBranchAddress("bb.tr.pz", bb_tr_pz, &b_bb_tr_pz);
   fChain->SetBranchAddress("Ndata.bb.tr.r_ph", &Ndata_bb_tr_r_ph, &b_Ndata_bb_tr_r_ph);
   fChain->SetBranchAddress("bb.tr.r_ph", bb_tr_r_ph, &b_bb_tr_r_ph);
   fChain->SetBranchAddress("Ndata.bb.tr.r_th", &Ndata_bb_tr_r_th, &b_Ndata_bb_tr_r_th);
   fChain->SetBranchAddress("bb.tr.r_th", bb_tr_r_th, &b_bb_tr_r_th);
   fChain->SetBranchAddress("Ndata.bb.tr.r_x", &Ndata_bb_tr_r_x, &b_Ndata_bb_tr_r_x);
   fChain->SetBranchAddress("bb.tr.r_x", bb_tr_r_x, &b_bb_tr_r_x);
   fChain->SetBranchAddress("Ndata.bb.tr.r_y", &Ndata_bb_tr_r_y, &b_Ndata_bb_tr_r_y);
   fChain->SetBranchAddress("bb.tr.r_y", bb_tr_r_y, &b_bb_tr_r_y);
   fChain->SetBranchAddress("Ndata.bb.tr.tg_dp", &Ndata_bb_tr_tg_dp, &b_Ndata_bb_tr_tg_dp);
   fChain->SetBranchAddress("bb.tr.tg_dp", bb_tr_tg_dp, &b_bb_tr_tg_dp);
   fChain->SetBranchAddress("Ndata.bb.tr.tg_ph", &Ndata_bb_tr_tg_ph, &b_Ndata_bb_tr_tg_ph);
   fChain->SetBranchAddress("bb.tr.tg_ph", bb_tr_tg_ph, &b_bb_tr_tg_ph);
   fChain->SetBranchAddress("Ndata.bb.tr.tg_th", &Ndata_bb_tr_tg_th, &b_Ndata_bb_tr_tg_th);
   fChain->SetBranchAddress("bb.tr.tg_th", bb_tr_tg_th, &b_bb_tr_tg_th);
   fChain->SetBranchAddress("Ndata.bb.tr.tg_x", &Ndata_bb_tr_tg_x, &b_Ndata_bb_tr_tg_x);
   fChain->SetBranchAddress("bb.tr.tg_x", bb_tr_tg_x, &b_bb_tr_tg_x);
   fChain->SetBranchAddress("Ndata.bb.tr.tg_y", &Ndata_bb_tr_tg_y, &b_Ndata_bb_tr_tg_y);
   fChain->SetBranchAddress("bb.tr.tg_y", bb_tr_tg_y, &b_bb_tr_tg_y);
   fChain->SetBranchAddress("Ndata.bb.tr.th", &Ndata_bb_tr_th, &b_Ndata_bb_tr_th);
   fChain->SetBranchAddress("bb.tr.th", bb_tr_th, &b_bb_tr_th);
   fChain->SetBranchAddress("Ndata.bb.tr.time", &Ndata_bb_tr_time, &b_Ndata_bb_tr_time);
   fChain->SetBranchAddress("bb.tr.time", bb_tr_time, &b_bb_tr_time);
   fChain->SetBranchAddress("Ndata.bb.tr.vx", &Ndata_bb_tr_vx, &b_Ndata_bb_tr_vx);
   fChain->SetBranchAddress("bb.tr.vx", bb_tr_vx, &b_bb_tr_vx);
   fChain->SetBranchAddress("Ndata.bb.tr.vy", &Ndata_bb_tr_vy, &b_Ndata_bb_tr_vy);
   fChain->SetBranchAddress("bb.tr.vy", bb_tr_vy, &b_bb_tr_vy);
   fChain->SetBranchAddress("Ndata.bb.tr.vz", &Ndata_bb_tr_vz, &b_Ndata_bb_tr_vz);
   fChain->SetBranchAddress("bb.tr.vz", bb_tr_vz, &b_bb_tr_vz);
   fChain->SetBranchAddress("Ndata.bb.tr.x", &Ndata_bb_tr_x, &b_Ndata_bb_tr_x);
   fChain->SetBranchAddress("bb.tr.x", bb_tr_x, &b_bb_tr_x);
   fChain->SetBranchAddress("Ndata.bb.tr.y", &Ndata_bb_tr_y, &b_Ndata_bb_tr_y);
   fChain->SetBranchAddress("bb.tr.y", bb_tr_y, &b_bb_tr_y);
   fChain->SetBranchAddress("Ndata.bb.x_bcp", &Ndata_bb_x_bcp, &b_Ndata_bb_x_bcp);
   fChain->SetBranchAddress("bb.x_bcp", bb_x_bcp, &b_bb_x_bcp);
   fChain->SetBranchAddress("Ndata.bb.x_fcp", &Ndata_bb_x_fcp, &b_Ndata_bb_x_fcp);
   fChain->SetBranchAddress("bb.x_fcp", bb_x_fcp, &b_bb_x_fcp);
   fChain->SetBranchAddress("Ndata.bb.y_bcp", &Ndata_bb_y_bcp, &b_Ndata_bb_y_bcp);
   fChain->SetBranchAddress("bb.y_bcp", bb_y_bcp, &b_bb_y_bcp);
   fChain->SetBranchAddress("Ndata.bb.y_fcp", &Ndata_bb_y_fcp, &b_Ndata_bb_y_fcp);
   fChain->SetBranchAddress("bb.y_fcp", bb_y_fcp, &b_bb_y_fcp);
   fChain->SetBranchAddress("Ndata.bb.z_bcp", &Ndata_bb_z_bcp, &b_Ndata_bb_z_bcp);
   fChain->SetBranchAddress("bb.z_bcp", bb_z_bcp, &b_bb_z_bcp);
   fChain->SetBranchAddress("Ndata.bb.z_fcp", &Ndata_bb_z_fcp, &b_Ndata_bb_z_fcp);
   fChain->SetBranchAddress("bb.z_fcp", bb_z_fcp, &b_bb_z_fcp);
   fChain->SetBranchAddress("Ndata.sbs.hcal.Ref.tdc", &Ndata_sbs_hcal_Ref_tdc, &b_Ndata_sbs_hcal_Ref_tdc);
   fChain->SetBranchAddress("sbs.hcal.Ref.tdc", sbs_hcal_Ref_tdc, &b_sbs_hcal_Ref_tdc);
   fChain->SetBranchAddress("Ndata.sbs.hcal.Ref.tdc_mult", &Ndata_sbs_hcal_Ref_tdc_mult, &b_Ndata_sbs_hcal_Ref_tdc_mult);
   fChain->SetBranchAddress("sbs.hcal.Ref.tdc_mult", sbs_hcal_Ref_tdc_mult, &b_sbs_hcal_Ref_tdc_mult);
   fChain->SetBranchAddress("Ndata.sbs.hcal.Ref.tdcelemID", &Ndata_sbs_hcal_Ref_tdcelemID, &b_Ndata_sbs_hcal_Ref_tdcelemID);
   fChain->SetBranchAddress("sbs.hcal.Ref.tdcelemID", sbs_hcal_Ref_tdcelemID, &b_sbs_hcal_Ref_tdcelemID);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.again", &Ndata_sbs_hcal_clus_again, &b_Ndata_sbs_hcal_clus_again);
   fChain->SetBranchAddress("sbs.hcal.clus.again", sbs_hcal_clus_again, &b_sbs_hcal_clus_again);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.atime", &Ndata_sbs_hcal_clus_atime, &b_Ndata_sbs_hcal_clus_atime);
   fChain->SetBranchAddress("sbs.hcal.clus.atime", sbs_hcal_clus_atime, &b_sbs_hcal_clus_atime);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.col", &Ndata_sbs_hcal_clus_col, &b_Ndata_sbs_hcal_clus_col);
   fChain->SetBranchAddress("sbs.hcal.clus.col", sbs_hcal_clus_col, &b_sbs_hcal_clus_col);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.e", &Ndata_sbs_hcal_clus_e, &b_Ndata_sbs_hcal_clus_e);
   fChain->SetBranchAddress("sbs.hcal.clus.e", sbs_hcal_clus_e, &b_sbs_hcal_clus_e);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.eblk", &Ndata_sbs_hcal_clus_eblk, &b_Ndata_sbs_hcal_clus_eblk);
   fChain->SetBranchAddress("sbs.hcal.clus.eblk", sbs_hcal_clus_eblk, &b_sbs_hcal_clus_eblk);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.id", &Ndata_sbs_hcal_clus_id, &b_Ndata_sbs_hcal_clus_id);
   fChain->SetBranchAddress("sbs.hcal.clus.id", sbs_hcal_clus_id, &b_sbs_hcal_clus_id);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.nblk", &Ndata_sbs_hcal_clus_nblk, &b_Ndata_sbs_hcal_clus_nblk);
   fChain->SetBranchAddress("sbs.hcal.clus.nblk", sbs_hcal_clus_nblk, &b_sbs_hcal_clus_nblk);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.row", &Ndata_sbs_hcal_clus_row, &b_Ndata_sbs_hcal_clus_row);
   fChain->SetBranchAddress("sbs.hcal.clus.row", sbs_hcal_clus_row, &b_sbs_hcal_clus_row);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.tdctime", &Ndata_sbs_hcal_clus_tdctime, &b_Ndata_sbs_hcal_clus_tdctime);
   fChain->SetBranchAddress("sbs.hcal.clus.tdctime", sbs_hcal_clus_tdctime, &b_sbs_hcal_clus_tdctime);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.x", &Ndata_sbs_hcal_clus_x, &b_Ndata_sbs_hcal_clus_x);
   fChain->SetBranchAddress("sbs.hcal.clus.x", sbs_hcal_clus_x, &b_sbs_hcal_clus_x);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.y", &Ndata_sbs_hcal_clus_y, &b_Ndata_sbs_hcal_clus_y);
   fChain->SetBranchAddress("sbs.hcal.clus.y", sbs_hcal_clus_y, &b_sbs_hcal_clus_y);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.again", &Ndata_sbs_hcal_clus_blk_again, &b_Ndata_sbs_hcal_clus_blk_again);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.again", sbs_hcal_clus_blk_again, &b_sbs_hcal_clus_blk_again);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.atime", &Ndata_sbs_hcal_clus_blk_atime, &b_Ndata_sbs_hcal_clus_blk_atime);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.atime", sbs_hcal_clus_blk_atime, &b_sbs_hcal_clus_blk_atime);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.col", &Ndata_sbs_hcal_clus_blk_col, &b_Ndata_sbs_hcal_clus_blk_col);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.col", sbs_hcal_clus_blk_col, &b_sbs_hcal_clus_blk_col);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.e", &Ndata_sbs_hcal_clus_blk_e, &b_Ndata_sbs_hcal_clus_blk_e);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.e", sbs_hcal_clus_blk_e, &b_sbs_hcal_clus_blk_e);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.id", &Ndata_sbs_hcal_clus_blk_id, &b_Ndata_sbs_hcal_clus_blk_id);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.id", sbs_hcal_clus_blk_id, &b_sbs_hcal_clus_blk_id);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.row", &Ndata_sbs_hcal_clus_blk_row, &b_Ndata_sbs_hcal_clus_blk_row);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.row", sbs_hcal_clus_blk_row, &b_sbs_hcal_clus_blk_row);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.tdctime", &Ndata_sbs_hcal_clus_blk_tdctime, &b_Ndata_sbs_hcal_clus_blk_tdctime);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.tdctime", sbs_hcal_clus_blk_tdctime, &b_sbs_hcal_clus_blk_tdctime);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.x", &Ndata_sbs_hcal_clus_blk_x, &b_Ndata_sbs_hcal_clus_blk_x);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.x", sbs_hcal_clus_blk_x, &b_sbs_hcal_clus_blk_x);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.y", &Ndata_sbs_hcal_clus_blk_y, &b_Ndata_sbs_hcal_clus_blk_y);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.y", sbs_hcal_clus_blk_y, &b_sbs_hcal_clus_blk_y);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.atime", &Ndata_sbs_hcal_goodblock_atime, &b_Ndata_sbs_hcal_goodblock_atime);
   fChain->SetBranchAddress("sbs.hcal.goodblock.atime", sbs_hcal_goodblock_atime, &b_sbs_hcal_goodblock_atime);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.cid", &Ndata_sbs_hcal_goodblock_cid, &b_Ndata_sbs_hcal_goodblock_cid);
   fChain->SetBranchAddress("sbs.hcal.goodblock.cid", sbs_hcal_goodblock_cid, &b_sbs_hcal_goodblock_cid);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.col", &Ndata_sbs_hcal_goodblock_col, &b_Ndata_sbs_hcal_goodblock_col);
   fChain->SetBranchAddress("sbs.hcal.goodblock.col", sbs_hcal_goodblock_col, &b_sbs_hcal_goodblock_col);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.e", &Ndata_sbs_hcal_goodblock_e, &b_Ndata_sbs_hcal_goodblock_e);
   fChain->SetBranchAddress("sbs.hcal.goodblock.e", sbs_hcal_goodblock_e, &b_sbs_hcal_goodblock_e);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.id", &Ndata_sbs_hcal_goodblock_id, &b_Ndata_sbs_hcal_goodblock_id);
   fChain->SetBranchAddress("sbs.hcal.goodblock.id", sbs_hcal_goodblock_id, &b_sbs_hcal_goodblock_id);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.row", &Ndata_sbs_hcal_goodblock_row, &b_Ndata_sbs_hcal_goodblock_row);
   fChain->SetBranchAddress("sbs.hcal.goodblock.row", sbs_hcal_goodblock_row, &b_sbs_hcal_goodblock_row);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.tdctime", &Ndata_sbs_hcal_goodblock_tdctime, &b_Ndata_sbs_hcal_goodblock_tdctime);
   fChain->SetBranchAddress("sbs.hcal.goodblock.tdctime", sbs_hcal_goodblock_tdctime, &b_sbs_hcal_goodblock_tdctime);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.x", &Ndata_sbs_hcal_goodblock_x, &b_Ndata_sbs_hcal_goodblock_x);
   fChain->SetBranchAddress("sbs.hcal.goodblock.x", sbs_hcal_goodblock_x, &b_sbs_hcal_goodblock_x);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.y", &Ndata_sbs_hcal_goodblock_y, &b_Ndata_sbs_hcal_goodblock_y);
   fChain->SetBranchAddress("sbs.hcal.goodblock.y", sbs_hcal_goodblock_y, &b_sbs_hcal_goodblock_y);
   fChain->SetBranchAddress("Ndata.sbs.tr.beta", &Ndata_sbs_tr_beta, &b_Ndata_sbs_tr_beta);
   fChain->SetBranchAddress("sbs.tr.beta", &sbs_tr_beta, &b_sbs_tr_beta);
   fChain->SetBranchAddress("Ndata.sbs.tr.chi2", &Ndata_sbs_tr_chi2, &b_Ndata_sbs_tr_chi2);
   fChain->SetBranchAddress("sbs.tr.chi2", &sbs_tr_chi2, &b_sbs_tr_chi2);
   fChain->SetBranchAddress("Ndata.sbs.tr.d_ph", &Ndata_sbs_tr_d_ph, &b_Ndata_sbs_tr_d_ph);
   fChain->SetBranchAddress("sbs.tr.d_ph", &sbs_tr_d_ph, &b_sbs_tr_d_ph);
   fChain->SetBranchAddress("Ndata.sbs.tr.d_th", &Ndata_sbs_tr_d_th, &b_Ndata_sbs_tr_d_th);
   fChain->SetBranchAddress("sbs.tr.d_th", &sbs_tr_d_th, &b_sbs_tr_d_th);
   fChain->SetBranchAddress("Ndata.sbs.tr.d_x", &Ndata_sbs_tr_d_x, &b_Ndata_sbs_tr_d_x);
   fChain->SetBranchAddress("sbs.tr.d_x", &sbs_tr_d_x, &b_sbs_tr_d_x);
   fChain->SetBranchAddress("Ndata.sbs.tr.d_y", &Ndata_sbs_tr_d_y, &b_Ndata_sbs_tr_d_y);
   fChain->SetBranchAddress("sbs.tr.d_y", &sbs_tr_d_y, &b_sbs_tr_d_y);
   fChain->SetBranchAddress("Ndata.sbs.tr.dbeta", &Ndata_sbs_tr_dbeta, &b_Ndata_sbs_tr_dbeta);
   fChain->SetBranchAddress("sbs.tr.dbeta", &sbs_tr_dbeta, &b_sbs_tr_dbeta);
   fChain->SetBranchAddress("Ndata.sbs.tr.dtime", &Ndata_sbs_tr_dtime, &b_Ndata_sbs_tr_dtime);
   fChain->SetBranchAddress("sbs.tr.dtime", &sbs_tr_dtime, &b_sbs_tr_dtime);
   fChain->SetBranchAddress("Ndata.sbs.tr.flag", &Ndata_sbs_tr_flag, &b_Ndata_sbs_tr_flag);
   fChain->SetBranchAddress("sbs.tr.flag", &sbs_tr_flag, &b_sbs_tr_flag);
   fChain->SetBranchAddress("Ndata.sbs.tr.ndof", &Ndata_sbs_tr_ndof, &b_Ndata_sbs_tr_ndof);
   fChain->SetBranchAddress("sbs.tr.ndof", &sbs_tr_ndof, &b_sbs_tr_ndof);
   fChain->SetBranchAddress("Ndata.sbs.tr.p", &Ndata_sbs_tr_p, &b_Ndata_sbs_tr_p);
   fChain->SetBranchAddress("sbs.tr.p", &sbs_tr_p, &b_sbs_tr_p);
   fChain->SetBranchAddress("Ndata.sbs.tr.pathl", &Ndata_sbs_tr_pathl, &b_Ndata_sbs_tr_pathl);
   fChain->SetBranchAddress("sbs.tr.pathl", &sbs_tr_pathl, &b_sbs_tr_pathl);
   fChain->SetBranchAddress("Ndata.sbs.tr.ph", &Ndata_sbs_tr_ph, &b_Ndata_sbs_tr_ph);
   fChain->SetBranchAddress("sbs.tr.ph", &sbs_tr_ph, &b_sbs_tr_ph);
   fChain->SetBranchAddress("Ndata.sbs.tr.px", &Ndata_sbs_tr_px, &b_Ndata_sbs_tr_px);
   fChain->SetBranchAddress("sbs.tr.px", &sbs_tr_px, &b_sbs_tr_px);
   fChain->SetBranchAddress("Ndata.sbs.tr.py", &Ndata_sbs_tr_py, &b_Ndata_sbs_tr_py);
   fChain->SetBranchAddress("sbs.tr.py", &sbs_tr_py, &b_sbs_tr_py);
   fChain->SetBranchAddress("Ndata.sbs.tr.pz", &Ndata_sbs_tr_pz, &b_Ndata_sbs_tr_pz);
   fChain->SetBranchAddress("sbs.tr.pz", &sbs_tr_pz, &b_sbs_tr_pz);
   fChain->SetBranchAddress("Ndata.sbs.tr.r_ph", &Ndata_sbs_tr_r_ph, &b_Ndata_sbs_tr_r_ph);
   fChain->SetBranchAddress("sbs.tr.r_ph", &sbs_tr_r_ph, &b_sbs_tr_r_ph);
   fChain->SetBranchAddress("Ndata.sbs.tr.r_th", &Ndata_sbs_tr_r_th, &b_Ndata_sbs_tr_r_th);
   fChain->SetBranchAddress("sbs.tr.r_th", &sbs_tr_r_th, &b_sbs_tr_r_th);
   fChain->SetBranchAddress("Ndata.sbs.tr.r_x", &Ndata_sbs_tr_r_x, &b_Ndata_sbs_tr_r_x);
   fChain->SetBranchAddress("sbs.tr.r_x", &sbs_tr_r_x, &b_sbs_tr_r_x);
   fChain->SetBranchAddress("Ndata.sbs.tr.r_y", &Ndata_sbs_tr_r_y, &b_Ndata_sbs_tr_r_y);
   fChain->SetBranchAddress("sbs.tr.r_y", &sbs_tr_r_y, &b_sbs_tr_r_y);
   fChain->SetBranchAddress("Ndata.sbs.tr.tg_dp", &Ndata_sbs_tr_tg_dp, &b_Ndata_sbs_tr_tg_dp);
   fChain->SetBranchAddress("sbs.tr.tg_dp", &sbs_tr_tg_dp, &b_sbs_tr_tg_dp);
   fChain->SetBranchAddress("Ndata.sbs.tr.tg_ph", &Ndata_sbs_tr_tg_ph, &b_Ndata_sbs_tr_tg_ph);
   fChain->SetBranchAddress("sbs.tr.tg_ph", &sbs_tr_tg_ph, &b_sbs_tr_tg_ph);
   fChain->SetBranchAddress("Ndata.sbs.tr.tg_th", &Ndata_sbs_tr_tg_th, &b_Ndata_sbs_tr_tg_th);
   fChain->SetBranchAddress("sbs.tr.tg_th", &sbs_tr_tg_th, &b_sbs_tr_tg_th);
   fChain->SetBranchAddress("Ndata.sbs.tr.tg_y", &Ndata_sbs_tr_tg_y, &b_Ndata_sbs_tr_tg_y);
   fChain->SetBranchAddress("sbs.tr.tg_y", &sbs_tr_tg_y, &b_sbs_tr_tg_y);
   fChain->SetBranchAddress("Ndata.sbs.tr.th", &Ndata_sbs_tr_th, &b_Ndata_sbs_tr_th);
   fChain->SetBranchAddress("sbs.tr.th", &sbs_tr_th, &b_sbs_tr_th);
   fChain->SetBranchAddress("Ndata.sbs.tr.time", &Ndata_sbs_tr_time, &b_Ndata_sbs_tr_time);
   fChain->SetBranchAddress("sbs.tr.time", &sbs_tr_time, &b_sbs_tr_time);
   fChain->SetBranchAddress("Ndata.sbs.tr.vx", &Ndata_sbs_tr_vx, &b_Ndata_sbs_tr_vx);
   fChain->SetBranchAddress("sbs.tr.vx", &sbs_tr_vx, &b_sbs_tr_vx);
   fChain->SetBranchAddress("Ndata.sbs.tr.vy", &Ndata_sbs_tr_vy, &b_Ndata_sbs_tr_vy);
   fChain->SetBranchAddress("sbs.tr.vy", &sbs_tr_vy, &b_sbs_tr_vy);
   fChain->SetBranchAddress("Ndata.sbs.tr.vz", &Ndata_sbs_tr_vz, &b_Ndata_sbs_tr_vz);
   fChain->SetBranchAddress("sbs.tr.vz", &sbs_tr_vz, &b_sbs_tr_vz);
   fChain->SetBranchAddress("Ndata.sbs.tr.x", &Ndata_sbs_tr_x, &b_Ndata_sbs_tr_x);
   fChain->SetBranchAddress("sbs.tr.x", &sbs_tr_x, &b_sbs_tr_x);
   fChain->SetBranchAddress("Ndata.sbs.tr.y", &Ndata_sbs_tr_y, &b_Ndata_sbs_tr_y);
   fChain->SetBranchAddress("sbs.tr.y", &sbs_tr_y, &b_sbs_tr_y);
   fChain->SetBranchAddress("Ndata.sbs.trig.a_amp_p", &Ndata_sbs_trig_a_amp_p, &b_Ndata_sbs_trig_a_amp_p);
   fChain->SetBranchAddress("sbs.trig.a_amp_p", sbs_trig_a_amp_p, &b_sbs_trig_a_amp_p);
   fChain->SetBranchAddress("Ndata.sbs.trig.a_p", &Ndata_sbs_trig_a_p, &b_Ndata_sbs_trig_a_p);
   fChain->SetBranchAddress("sbs.trig.a_p", sbs_trig_a_p, &b_sbs_trig_a_p);
   fChain->SetBranchAddress("Ndata.sbs.trig.a_time", &Ndata_sbs_trig_a_time, &b_Ndata_sbs_trig_a_time);
   fChain->SetBranchAddress("sbs.trig.a_time", sbs_trig_a_time, &b_sbs_trig_a_time);
   fChain->SetBranchAddress("Ndata.sbs.trig.adcelemID", &Ndata_sbs_trig_adcelemID, &b_Ndata_sbs_trig_adcelemID);
   fChain->SetBranchAddress("sbs.trig.adcelemID", sbs_trig_adcelemID, &b_sbs_trig_adcelemID);
   fChain->SetBranchAddress("Ndata.sbs.x_bcp", &Ndata_sbs_x_bcp, &b_Ndata_sbs_x_bcp);
   fChain->SetBranchAddress("sbs.x_bcp", sbs_x_bcp, &b_sbs_x_bcp);
   fChain->SetBranchAddress("Ndata.sbs.x_fcp", &Ndata_sbs_x_fcp, &b_Ndata_sbs_x_fcp);
   fChain->SetBranchAddress("sbs.x_fcp", sbs_x_fcp, &b_sbs_x_fcp);
   fChain->SetBranchAddress("Ndata.sbs.y_bcp", &Ndata_sbs_y_bcp, &b_Ndata_sbs_y_bcp);
   fChain->SetBranchAddress("sbs.y_bcp", sbs_y_bcp, &b_sbs_y_bcp);
   fChain->SetBranchAddress("Ndata.sbs.y_fcp", &Ndata_sbs_y_fcp, &b_Ndata_sbs_y_fcp);
   fChain->SetBranchAddress("sbs.y_fcp", sbs_y_fcp, &b_sbs_y_fcp);
   fChain->SetBranchAddress("Ndata.sbs.z_bcp", &Ndata_sbs_z_bcp, &b_Ndata_sbs_z_bcp);
   fChain->SetBranchAddress("sbs.z_bcp", sbs_z_bcp, &b_sbs_z_bcp);
   fChain->SetBranchAddress("Ndata.sbs.z_fcp", &Ndata_sbs_z_fcp, &b_Ndata_sbs_z_fcp);
   fChain->SetBranchAddress("sbs.z_fcp", sbs_z_fcp, &b_sbs_z_fcp);
   fChain->SetBranchAddress("BB.gold.beta", &BB_gold_beta, &b_BB_gold_beta);
   fChain->SetBranchAddress("BB.gold.dp", &BB_gold_dp, &b_BB_gold_dp);
   fChain->SetBranchAddress("BB.gold.index", &BB_gold_index, &b_BB_gold_index);
   fChain->SetBranchAddress("BB.gold.ok", &BB_gold_ok, &b_BB_gold_ok);
   fChain->SetBranchAddress("BB.gold.p", &BB_gold_p, &b_BB_gold_p);
   fChain->SetBranchAddress("BB.gold.ph", &BB_gold_ph, &b_BB_gold_ph);
   fChain->SetBranchAddress("BB.gold.px", &BB_gold_px, &b_BB_gold_px);
   fChain->SetBranchAddress("BB.gold.py", &BB_gold_py, &b_BB_gold_py);
   fChain->SetBranchAddress("BB.gold.pz", &BB_gold_pz, &b_BB_gold_pz);
   fChain->SetBranchAddress("BB.gold.th", &BB_gold_th, &b_BB_gold_th);
   fChain->SetBranchAddress("BB.gold.x", &BB_gold_x, &b_BB_gold_x);
   fChain->SetBranchAddress("BB.gold.y", &BB_gold_y, &b_BB_gold_y);
   fChain->SetBranchAddress("Lrb.BPMA.rawcur.1", &Lrb_BPMA_rawcur_1, &b_Lrb_BPMA_rawcur_1);
   fChain->SetBranchAddress("Lrb.BPMA.rawcur.2", &Lrb_BPMA_rawcur_2, &b_Lrb_BPMA_rawcur_2);
   fChain->SetBranchAddress("Lrb.BPMA.rawcur.3", &Lrb_BPMA_rawcur_3, &b_Lrb_BPMA_rawcur_3);
   fChain->SetBranchAddress("Lrb.BPMA.rawcur.4", &Lrb_BPMA_rawcur_4, &b_Lrb_BPMA_rawcur_4);
   fChain->SetBranchAddress("Lrb.BPMA.rotpos1", &Lrb_BPMA_rotpos1, &b_Lrb_BPMA_rotpos1);
   fChain->SetBranchAddress("Lrb.BPMA.rotpos2", &Lrb_BPMA_rotpos2, &b_Lrb_BPMA_rotpos2);
   fChain->SetBranchAddress("Lrb.BPMA.x", &Lrb_BPMA_x, &b_Lrb_BPMA_x);
   fChain->SetBranchAddress("Lrb.BPMA.y", &Lrb_BPMA_y, &b_Lrb_BPMA_y);
   fChain->SetBranchAddress("Lrb.BPMA.z", &Lrb_BPMA_z, &b_Lrb_BPMA_z);
   fChain->SetBranchAddress("Lrb.BPMB.rawcur.1", &Lrb_BPMB_rawcur_1, &b_Lrb_BPMB_rawcur_1);
   fChain->SetBranchAddress("Lrb.BPMB.rawcur.2", &Lrb_BPMB_rawcur_2, &b_Lrb_BPMB_rawcur_2);
   fChain->SetBranchAddress("Lrb.BPMB.rawcur.3", &Lrb_BPMB_rawcur_3, &b_Lrb_BPMB_rawcur_3);
   fChain->SetBranchAddress("Lrb.BPMB.rawcur.4", &Lrb_BPMB_rawcur_4, &b_Lrb_BPMB_rawcur_4);
   fChain->SetBranchAddress("Lrb.BPMB.rotpos1", &Lrb_BPMB_rotpos1, &b_Lrb_BPMB_rotpos1);
   fChain->SetBranchAddress("Lrb.BPMB.rotpos2", &Lrb_BPMB_rotpos2, &b_Lrb_BPMB_rotpos2);
   fChain->SetBranchAddress("Lrb.BPMB.x", &Lrb_BPMB_x, &b_Lrb_BPMB_x);
   fChain->SetBranchAddress("Lrb.BPMB.y", &Lrb_BPMB_y, &b_Lrb_BPMB_y);
   fChain->SetBranchAddress("Lrb.BPMB.z", &Lrb_BPMB_z, &b_Lrb_BPMB_z);
   fChain->SetBranchAddress("Lrb.Raster.bpma.x", &Lrb_Raster_bpma_x, &b_Lrb_Raster_bpma_x);
   fChain->SetBranchAddress("Lrb.Raster.bpma.y", &Lrb_Raster_bpma_y, &b_Lrb_Raster_bpma_y);
   fChain->SetBranchAddress("Lrb.Raster.bpma.z", &Lrb_Raster_bpma_z, &b_Lrb_Raster_bpma_z);
   fChain->SetBranchAddress("Lrb.Raster.bpmb.x", &Lrb_Raster_bpmb_x, &b_Lrb_Raster_bpmb_x);
   fChain->SetBranchAddress("Lrb.Raster.bpmb.y", &Lrb_Raster_bpmb_y, &b_Lrb_Raster_bpmb_y);
   fChain->SetBranchAddress("Lrb.Raster.bpmb.z", &Lrb_Raster_bpmb_z, &b_Lrb_Raster_bpmb_z);
   fChain->SetBranchAddress("Lrb.Raster.rawcur.x", &Lrb_Raster_rawcur_x, &b_Lrb_Raster_rawcur_x);
   fChain->SetBranchAddress("Lrb.Raster.rawcur.y", &Lrb_Raster_rawcur_y, &b_Lrb_Raster_rawcur_y);
   fChain->SetBranchAddress("Lrb.Raster.rawslope.x", &Lrb_Raster_rawslope_x, &b_Lrb_Raster_rawslope_x);
   fChain->SetBranchAddress("Lrb.Raster.rawslope.y", &Lrb_Raster_rawslope_y, &b_Lrb_Raster_rawslope_y);
   fChain->SetBranchAddress("Lrb.Raster.target.dir.x", &Lrb_Raster_target_dir_x, &b_Lrb_Raster_target_dir_x);
   fChain->SetBranchAddress("Lrb.Raster.target.dir.y", &Lrb_Raster_target_dir_y, &b_Lrb_Raster_target_dir_y);
   fChain->SetBranchAddress("Lrb.Raster.target.dir.z", &Lrb_Raster_target_dir_z, &b_Lrb_Raster_target_dir_z);
   fChain->SetBranchAddress("Lrb.Raster.target.x", &Lrb_Raster_target_x, &b_Lrb_Raster_target_x);
   fChain->SetBranchAddress("Lrb.Raster.target.y", &Lrb_Raster_target_y, &b_Lrb_Raster_target_y);
   fChain->SetBranchAddress("Lrb.Raster.target.z", &Lrb_Raster_target_z, &b_Lrb_Raster_target_z);
   fChain->SetBranchAddress("Lrb.Raster2.bpma.x", &Lrb_Raster2_bpma_x, &b_Lrb_Raster2_bpma_x);
   fChain->SetBranchAddress("Lrb.Raster2.bpma.y", &Lrb_Raster2_bpma_y, &b_Lrb_Raster2_bpma_y);
   fChain->SetBranchAddress("Lrb.Raster2.bpma.z", &Lrb_Raster2_bpma_z, &b_Lrb_Raster2_bpma_z);
   fChain->SetBranchAddress("Lrb.Raster2.bpmb.x", &Lrb_Raster2_bpmb_x, &b_Lrb_Raster2_bpmb_x);
   fChain->SetBranchAddress("Lrb.Raster2.bpmb.y", &Lrb_Raster2_bpmb_y, &b_Lrb_Raster2_bpmb_y);
   fChain->SetBranchAddress("Lrb.Raster2.bpmb.z", &Lrb_Raster2_bpmb_z, &b_Lrb_Raster2_bpmb_z);
   fChain->SetBranchAddress("Lrb.Raster2.rawcur.x", &Lrb_Raster2_rawcur_x, &b_Lrb_Raster2_rawcur_x);
   fChain->SetBranchAddress("Lrb.Raster2.rawcur.y", &Lrb_Raster2_rawcur_y, &b_Lrb_Raster2_rawcur_y);
   fChain->SetBranchAddress("Lrb.Raster2.rawslope.x", &Lrb_Raster2_rawslope_x, &b_Lrb_Raster2_rawslope_x);
   fChain->SetBranchAddress("Lrb.Raster2.rawslope.y", &Lrb_Raster2_rawslope_y, &b_Lrb_Raster2_rawslope_y);
   fChain->SetBranchAddress("Lrb.Raster2.target.dir.x", &Lrb_Raster2_target_dir_x, &b_Lrb_Raster2_target_dir_x);
   fChain->SetBranchAddress("Lrb.Raster2.target.dir.y", &Lrb_Raster2_target_dir_y, &b_Lrb_Raster2_target_dir_y);
   fChain->SetBranchAddress("Lrb.Raster2.target.dir.z", &Lrb_Raster2_target_dir_z, &b_Lrb_Raster2_target_dir_z);
   fChain->SetBranchAddress("Lrb.Raster2.target.x", &Lrb_Raster2_target_x, &b_Lrb_Raster2_target_x);
   fChain->SetBranchAddress("Lrb.Raster2.target.y", &Lrb_Raster2_target_y, &b_Lrb_Raster2_target_y);
   fChain->SetBranchAddress("Lrb.Raster2.target.z", &Lrb_Raster2_target_z, &b_Lrb_Raster2_target_z);
   fChain->SetBranchAddress("SBSrb.BPMA.rawcur.1", &SBSrb_BPMA_rawcur_1, &b_SBSrb_BPMA_rawcur_1);
   fChain->SetBranchAddress("SBSrb.BPMA.rawcur.2", &SBSrb_BPMA_rawcur_2, &b_SBSrb_BPMA_rawcur_2);
   fChain->SetBranchAddress("SBSrb.BPMA.rawcur.3", &SBSrb_BPMA_rawcur_3, &b_SBSrb_BPMA_rawcur_3);
   fChain->SetBranchAddress("SBSrb.BPMA.rawcur.4", &SBSrb_BPMA_rawcur_4, &b_SBSrb_BPMA_rawcur_4);
   fChain->SetBranchAddress("SBSrb.BPMA.rotpos1", &SBSrb_BPMA_rotpos1, &b_SBSrb_BPMA_rotpos1);
   fChain->SetBranchAddress("SBSrb.BPMA.rotpos2", &SBSrb_BPMA_rotpos2, &b_SBSrb_BPMA_rotpos2);
   fChain->SetBranchAddress("SBSrb.BPMA.x", &SBSrb_BPMA_x, &b_SBSrb_BPMA_x);
   fChain->SetBranchAddress("SBSrb.BPMA.y", &SBSrb_BPMA_y, &b_SBSrb_BPMA_y);
   fChain->SetBranchAddress("SBSrb.BPMA.z", &SBSrb_BPMA_z, &b_SBSrb_BPMA_z);
   fChain->SetBranchAddress("SBSrb.BPMB.rawcur.1", &SBSrb_BPMB_rawcur_1, &b_SBSrb_BPMB_rawcur_1);
   fChain->SetBranchAddress("SBSrb.BPMB.rawcur.2", &SBSrb_BPMB_rawcur_2, &b_SBSrb_BPMB_rawcur_2);
   fChain->SetBranchAddress("SBSrb.BPMB.rawcur.3", &SBSrb_BPMB_rawcur_3, &b_SBSrb_BPMB_rawcur_3);
   fChain->SetBranchAddress("SBSrb.BPMB.rawcur.4", &SBSrb_BPMB_rawcur_4, &b_SBSrb_BPMB_rawcur_4);
   fChain->SetBranchAddress("SBSrb.BPMB.rotpos1", &SBSrb_BPMB_rotpos1, &b_SBSrb_BPMB_rotpos1);
   fChain->SetBranchAddress("SBSrb.BPMB.rotpos2", &SBSrb_BPMB_rotpos2, &b_SBSrb_BPMB_rotpos2);
   fChain->SetBranchAddress("SBSrb.BPMB.x", &SBSrb_BPMB_x, &b_SBSrb_BPMB_x);
   fChain->SetBranchAddress("SBSrb.BPMB.y", &SBSrb_BPMB_y, &b_SBSrb_BPMB_y);
   fChain->SetBranchAddress("SBSrb.BPMB.z", &SBSrb_BPMB_z, &b_SBSrb_BPMB_z);
   fChain->SetBranchAddress("SBSrb.Raster.bpma.x", &SBSrb_Raster_bpma_x, &b_SBSrb_Raster_bpma_x);
   fChain->SetBranchAddress("SBSrb.Raster.bpma.y", &SBSrb_Raster_bpma_y, &b_SBSrb_Raster_bpma_y);
   fChain->SetBranchAddress("SBSrb.Raster.bpma.z", &SBSrb_Raster_bpma_z, &b_SBSrb_Raster_bpma_z);
   fChain->SetBranchAddress("SBSrb.Raster.bpmb.x", &SBSrb_Raster_bpmb_x, &b_SBSrb_Raster_bpmb_x);
   fChain->SetBranchAddress("SBSrb.Raster.bpmb.y", &SBSrb_Raster_bpmb_y, &b_SBSrb_Raster_bpmb_y);
   fChain->SetBranchAddress("SBSrb.Raster.bpmb.z", &SBSrb_Raster_bpmb_z, &b_SBSrb_Raster_bpmb_z);
   fChain->SetBranchAddress("SBSrb.Raster.rawcur.x", &SBSrb_Raster_rawcur_x, &b_SBSrb_Raster_rawcur_x);
   fChain->SetBranchAddress("SBSrb.Raster.rawcur.y", &SBSrb_Raster_rawcur_y, &b_SBSrb_Raster_rawcur_y);
   fChain->SetBranchAddress("SBSrb.Raster.rawslope.x", &SBSrb_Raster_rawslope_x, &b_SBSrb_Raster_rawslope_x);
   fChain->SetBranchAddress("SBSrb.Raster.rawslope.y", &SBSrb_Raster_rawslope_y, &b_SBSrb_Raster_rawslope_y);
   fChain->SetBranchAddress("SBSrb.Raster.target.dir.x", &SBSrb_Raster_target_dir_x, &b_SBSrb_Raster_target_dir_x);
   fChain->SetBranchAddress("SBSrb.Raster.target.dir.y", &SBSrb_Raster_target_dir_y, &b_SBSrb_Raster_target_dir_y);
   fChain->SetBranchAddress("SBSrb.Raster.target.dir.z", &SBSrb_Raster_target_dir_z, &b_SBSrb_Raster_target_dir_z);
   fChain->SetBranchAddress("SBSrb.Raster.target.x", &SBSrb_Raster_target_x, &b_SBSrb_Raster_target_x);
   fChain->SetBranchAddress("SBSrb.Raster.target.y", &SBSrb_Raster_target_y, &b_SBSrb_Raster_target_y);
   fChain->SetBranchAddress("SBSrb.Raster.target.z", &SBSrb_Raster_target_z, &b_SBSrb_Raster_target_z);
   fChain->SetBranchAddress("SBSrb.Raster2.bpma.x", &SBSrb_Raster2_bpma_x, &b_SBSrb_Raster2_bpma_x);
   fChain->SetBranchAddress("SBSrb.Raster2.bpma.y", &SBSrb_Raster2_bpma_y, &b_SBSrb_Raster2_bpma_y);
   fChain->SetBranchAddress("SBSrb.Raster2.bpma.z", &SBSrb_Raster2_bpma_z, &b_SBSrb_Raster2_bpma_z);
   fChain->SetBranchAddress("SBSrb.Raster2.bpmb.x", &SBSrb_Raster2_bpmb_x, &b_SBSrb_Raster2_bpmb_x);
   fChain->SetBranchAddress("SBSrb.Raster2.bpmb.y", &SBSrb_Raster2_bpmb_y, &b_SBSrb_Raster2_bpmb_y);
   fChain->SetBranchAddress("SBSrb.Raster2.bpmb.z", &SBSrb_Raster2_bpmb_z, &b_SBSrb_Raster2_bpmb_z);
   fChain->SetBranchAddress("SBSrb.Raster2.rawcur.x", &SBSrb_Raster2_rawcur_x, &b_SBSrb_Raster2_rawcur_x);
   fChain->SetBranchAddress("SBSrb.Raster2.rawcur.y", &SBSrb_Raster2_rawcur_y, &b_SBSrb_Raster2_rawcur_y);
   fChain->SetBranchAddress("SBSrb.Raster2.rawslope.x", &SBSrb_Raster2_rawslope_x, &b_SBSrb_Raster2_rawslope_x);
   fChain->SetBranchAddress("SBSrb.Raster2.rawslope.y", &SBSrb_Raster2_rawslope_y, &b_SBSrb_Raster2_rawslope_y);
   fChain->SetBranchAddress("SBSrb.Raster2.target.dir.x", &SBSrb_Raster2_target_dir_x, &b_SBSrb_Raster2_target_dir_x);
   fChain->SetBranchAddress("SBSrb.Raster2.target.dir.y", &SBSrb_Raster2_target_dir_y, &b_SBSrb_Raster2_target_dir_y);
   fChain->SetBranchAddress("SBSrb.Raster2.target.dir.z", &SBSrb_Raster2_target_dir_z, &b_SBSrb_Raster2_target_dir_z);
   fChain->SetBranchAddress("SBSrb.Raster2.target.x", &SBSrb_Raster2_target_x, &b_SBSrb_Raster2_target_x);
   fChain->SetBranchAddress("SBSrb.Raster2.target.y", &SBSrb_Raster2_target_y, &b_SBSrb_Raster2_target_y);
   fChain->SetBranchAddress("SBSrb.Raster2.target.z", &SBSrb_Raster2_target_z, &b_SBSrb_Raster2_target_z);
   fChain->SetBranchAddress("bb.gem.hit.ngoodhits", &bb_gem_hit_ngoodhits, &b_bb_gem_hit_ngoodhits);
   fChain->SetBranchAddress("bb.gem.m0.clust.nclustu", &bb_gem_m0_clust_nclustu, &b_bb_gem_m0_clust_nclustu);
   fChain->SetBranchAddress("bb.gem.m0.clust.nclustu_neg", &bb_gem_m0_clust_nclustu_neg, &b_bb_gem_m0_clust_nclustu_neg);
   fChain->SetBranchAddress("bb.gem.m0.clust.nclustu_tot", &bb_gem_m0_clust_nclustu_tot, &b_bb_gem_m0_clust_nclustu_tot);
   fChain->SetBranchAddress("bb.gem.m0.clust.nclustv", &bb_gem_m0_clust_nclustv, &b_bb_gem_m0_clust_nclustv);
   fChain->SetBranchAddress("bb.gem.m0.clust.nclustv_neg", &bb_gem_m0_clust_nclustv_neg, &b_bb_gem_m0_clust_nclustv_neg);
   fChain->SetBranchAddress("bb.gem.m0.clust.nclustv_tot", &bb_gem_m0_clust_nclustv_tot, &b_bb_gem_m0_clust_nclustv_tot);
   fChain->SetBranchAddress("bb.gem.m0.strip.nstrips_keep", &bb_gem_m0_strip_nstrips_keep, &b_bb_gem_m0_strip_nstrips_keep);
   fChain->SetBranchAddress("bb.gem.m0.strip.nstrips_keepU", &bb_gem_m0_strip_nstrips_keepU, &b_bb_gem_m0_strip_nstrips_keepU);
   fChain->SetBranchAddress("bb.gem.m0.strip.nstrips_keepV", &bb_gem_m0_strip_nstrips_keepV, &b_bb_gem_m0_strip_nstrips_keepV);
   fChain->SetBranchAddress("bb.gem.m0.strip.nstrips_keep_lmax", &bb_gem_m0_strip_nstrips_keep_lmax, &b_bb_gem_m0_strip_nstrips_keep_lmax);
   fChain->SetBranchAddress("bb.gem.m0.strip.nstrips_keep_lmaxU", &bb_gem_m0_strip_nstrips_keep_lmaxU, &b_bb_gem_m0_strip_nstrips_keep_lmaxU);
   fChain->SetBranchAddress("bb.gem.m0.strip.nstrips_keep_lmaxV", &bb_gem_m0_strip_nstrips_keep_lmaxV, &b_bb_gem_m0_strip_nstrips_keep_lmaxV);
   fChain->SetBranchAddress("bb.gem.m0.strip.nstripsfired", &bb_gem_m0_strip_nstripsfired, &b_bb_gem_m0_strip_nstripsfired);
   fChain->SetBranchAddress("bb.gem.m0.strip.nstripsfired_neg", &bb_gem_m0_strip_nstripsfired_neg, &b_bb_gem_m0_strip_nstripsfired_neg);
   fChain->SetBranchAddress("bb.gem.m1.clust.nclustu", &bb_gem_m1_clust_nclustu, &b_bb_gem_m1_clust_nclustu);
   fChain->SetBranchAddress("bb.gem.m1.clust.nclustu_neg", &bb_gem_m1_clust_nclustu_neg, &b_bb_gem_m1_clust_nclustu_neg);
   fChain->SetBranchAddress("bb.gem.m1.clust.nclustu_tot", &bb_gem_m1_clust_nclustu_tot, &b_bb_gem_m1_clust_nclustu_tot);
   fChain->SetBranchAddress("bb.gem.m1.clust.nclustv", &bb_gem_m1_clust_nclustv, &b_bb_gem_m1_clust_nclustv);
   fChain->SetBranchAddress("bb.gem.m1.clust.nclustv_neg", &bb_gem_m1_clust_nclustv_neg, &b_bb_gem_m1_clust_nclustv_neg);
   fChain->SetBranchAddress("bb.gem.m1.clust.nclustv_tot", &bb_gem_m1_clust_nclustv_tot, &b_bb_gem_m1_clust_nclustv_tot);
   fChain->SetBranchAddress("bb.gem.m1.strip.nstrips_keep", &bb_gem_m1_strip_nstrips_keep, &b_bb_gem_m1_strip_nstrips_keep);
   fChain->SetBranchAddress("bb.gem.m1.strip.nstrips_keepU", &bb_gem_m1_strip_nstrips_keepU, &b_bb_gem_m1_strip_nstrips_keepU);
   fChain->SetBranchAddress("bb.gem.m1.strip.nstrips_keepV", &bb_gem_m1_strip_nstrips_keepV, &b_bb_gem_m1_strip_nstrips_keepV);
   fChain->SetBranchAddress("bb.gem.m1.strip.nstrips_keep_lmax", &bb_gem_m1_strip_nstrips_keep_lmax, &b_bb_gem_m1_strip_nstrips_keep_lmax);
   fChain->SetBranchAddress("bb.gem.m1.strip.nstrips_keep_lmaxU", &bb_gem_m1_strip_nstrips_keep_lmaxU, &b_bb_gem_m1_strip_nstrips_keep_lmaxU);
   fChain->SetBranchAddress("bb.gem.m1.strip.nstrips_keep_lmaxV", &bb_gem_m1_strip_nstrips_keep_lmaxV, &b_bb_gem_m1_strip_nstrips_keep_lmaxV);
   fChain->SetBranchAddress("bb.gem.m1.strip.nstripsfired", &bb_gem_m1_strip_nstripsfired, &b_bb_gem_m1_strip_nstripsfired);
   fChain->SetBranchAddress("bb.gem.m1.strip.nstripsfired_neg", &bb_gem_m1_strip_nstripsfired_neg, &b_bb_gem_m1_strip_nstripsfired_neg);
   fChain->SetBranchAddress("bb.gem.m2.clust.nclustu", &bb_gem_m2_clust_nclustu, &b_bb_gem_m2_clust_nclustu);
   fChain->SetBranchAddress("bb.gem.m2.clust.nclustu_neg", &bb_gem_m2_clust_nclustu_neg, &b_bb_gem_m2_clust_nclustu_neg);
   fChain->SetBranchAddress("bb.gem.m2.clust.nclustu_tot", &bb_gem_m2_clust_nclustu_tot, &b_bb_gem_m2_clust_nclustu_tot);
   fChain->SetBranchAddress("bb.gem.m2.clust.nclustv", &bb_gem_m2_clust_nclustv, &b_bb_gem_m2_clust_nclustv);
   fChain->SetBranchAddress("bb.gem.m2.clust.nclustv_neg", &bb_gem_m2_clust_nclustv_neg, &b_bb_gem_m2_clust_nclustv_neg);
   fChain->SetBranchAddress("bb.gem.m2.clust.nclustv_tot", &bb_gem_m2_clust_nclustv_tot, &b_bb_gem_m2_clust_nclustv_tot);
   fChain->SetBranchAddress("bb.gem.m2.strip.nstrips_keep", &bb_gem_m2_strip_nstrips_keep, &b_bb_gem_m2_strip_nstrips_keep);
   fChain->SetBranchAddress("bb.gem.m2.strip.nstrips_keepU", &bb_gem_m2_strip_nstrips_keepU, &b_bb_gem_m2_strip_nstrips_keepU);
   fChain->SetBranchAddress("bb.gem.m2.strip.nstrips_keepV", &bb_gem_m2_strip_nstrips_keepV, &b_bb_gem_m2_strip_nstrips_keepV);
   fChain->SetBranchAddress("bb.gem.m2.strip.nstrips_keep_lmax", &bb_gem_m2_strip_nstrips_keep_lmax, &b_bb_gem_m2_strip_nstrips_keep_lmax);
   fChain->SetBranchAddress("bb.gem.m2.strip.nstrips_keep_lmaxU", &bb_gem_m2_strip_nstrips_keep_lmaxU, &b_bb_gem_m2_strip_nstrips_keep_lmaxU);
   fChain->SetBranchAddress("bb.gem.m2.strip.nstrips_keep_lmaxV", &bb_gem_m2_strip_nstrips_keep_lmaxV, &b_bb_gem_m2_strip_nstrips_keep_lmaxV);
   fChain->SetBranchAddress("bb.gem.m2.strip.nstripsfired", &bb_gem_m2_strip_nstripsfired, &b_bb_gem_m2_strip_nstripsfired);
   fChain->SetBranchAddress("bb.gem.m2.strip.nstripsfired_neg", &bb_gem_m2_strip_nstripsfired_neg, &b_bb_gem_m2_strip_nstripsfired_neg);
   fChain->SetBranchAddress("bb.gem.m3.clust.nclustu", &bb_gem_m3_clust_nclustu, &b_bb_gem_m3_clust_nclustu);
   fChain->SetBranchAddress("bb.gem.m3.clust.nclustu_neg", &bb_gem_m3_clust_nclustu_neg, &b_bb_gem_m3_clust_nclustu_neg);
   fChain->SetBranchAddress("bb.gem.m3.clust.nclustu_tot", &bb_gem_m3_clust_nclustu_tot, &b_bb_gem_m3_clust_nclustu_tot);
   fChain->SetBranchAddress("bb.gem.m3.clust.nclustv", &bb_gem_m3_clust_nclustv, &b_bb_gem_m3_clust_nclustv);
   fChain->SetBranchAddress("bb.gem.m3.clust.nclustv_neg", &bb_gem_m3_clust_nclustv_neg, &b_bb_gem_m3_clust_nclustv_neg);
   fChain->SetBranchAddress("bb.gem.m3.clust.nclustv_tot", &bb_gem_m3_clust_nclustv_tot, &b_bb_gem_m3_clust_nclustv_tot);
   fChain->SetBranchAddress("bb.gem.m3.strip.nstrips_keep", &bb_gem_m3_strip_nstrips_keep, &b_bb_gem_m3_strip_nstrips_keep);
   fChain->SetBranchAddress("bb.gem.m3.strip.nstrips_keepU", &bb_gem_m3_strip_nstrips_keepU, &b_bb_gem_m3_strip_nstrips_keepU);
   fChain->SetBranchAddress("bb.gem.m3.strip.nstrips_keepV", &bb_gem_m3_strip_nstrips_keepV, &b_bb_gem_m3_strip_nstrips_keepV);
   fChain->SetBranchAddress("bb.gem.m3.strip.nstrips_keep_lmax", &bb_gem_m3_strip_nstrips_keep_lmax, &b_bb_gem_m3_strip_nstrips_keep_lmax);
   fChain->SetBranchAddress("bb.gem.m3.strip.nstrips_keep_lmaxU", &bb_gem_m3_strip_nstrips_keep_lmaxU, &b_bb_gem_m3_strip_nstrips_keep_lmaxU);
   fChain->SetBranchAddress("bb.gem.m3.strip.nstrips_keep_lmaxV", &bb_gem_m3_strip_nstrips_keep_lmaxV, &b_bb_gem_m3_strip_nstrips_keep_lmaxV);
   fChain->SetBranchAddress("bb.gem.m3.strip.nstripsfired", &bb_gem_m3_strip_nstripsfired, &b_bb_gem_m3_strip_nstripsfired);
   fChain->SetBranchAddress("bb.gem.m3.strip.nstripsfired_neg", &bb_gem_m3_strip_nstripsfired_neg, &b_bb_gem_m3_strip_nstripsfired_neg);
   fChain->SetBranchAddress("bb.gem.m4.clust.nclustu", &bb_gem_m4_clust_nclustu, &b_bb_gem_m4_clust_nclustu);
   fChain->SetBranchAddress("bb.gem.m4.clust.nclustu_neg", &bb_gem_m4_clust_nclustu_neg, &b_bb_gem_m4_clust_nclustu_neg);
   fChain->SetBranchAddress("bb.gem.m4.clust.nclustu_tot", &bb_gem_m4_clust_nclustu_tot, &b_bb_gem_m4_clust_nclustu_tot);
   fChain->SetBranchAddress("bb.gem.m4.clust.nclustv", &bb_gem_m4_clust_nclustv, &b_bb_gem_m4_clust_nclustv);
   fChain->SetBranchAddress("bb.gem.m4.clust.nclustv_neg", &bb_gem_m4_clust_nclustv_neg, &b_bb_gem_m4_clust_nclustv_neg);
   fChain->SetBranchAddress("bb.gem.m4.clust.nclustv_tot", &bb_gem_m4_clust_nclustv_tot, &b_bb_gem_m4_clust_nclustv_tot);
   fChain->SetBranchAddress("bb.gem.m4.strip.nstrips_keep", &bb_gem_m4_strip_nstrips_keep, &b_bb_gem_m4_strip_nstrips_keep);
   fChain->SetBranchAddress("bb.gem.m4.strip.nstrips_keepU", &bb_gem_m4_strip_nstrips_keepU, &b_bb_gem_m4_strip_nstrips_keepU);
   fChain->SetBranchAddress("bb.gem.m4.strip.nstrips_keepV", &bb_gem_m4_strip_nstrips_keepV, &b_bb_gem_m4_strip_nstrips_keepV);
   fChain->SetBranchAddress("bb.gem.m4.strip.nstrips_keep_lmax", &bb_gem_m4_strip_nstrips_keep_lmax, &b_bb_gem_m4_strip_nstrips_keep_lmax);
   fChain->SetBranchAddress("bb.gem.m4.strip.nstrips_keep_lmaxU", &bb_gem_m4_strip_nstrips_keep_lmaxU, &b_bb_gem_m4_strip_nstrips_keep_lmaxU);
   fChain->SetBranchAddress("bb.gem.m4.strip.nstrips_keep_lmaxV", &bb_gem_m4_strip_nstrips_keep_lmaxV, &b_bb_gem_m4_strip_nstrips_keep_lmaxV);
   fChain->SetBranchAddress("bb.gem.m4.strip.nstripsfired", &bb_gem_m4_strip_nstripsfired, &b_bb_gem_m4_strip_nstripsfired);
   fChain->SetBranchAddress("bb.gem.m4.strip.nstripsfired_neg", &bb_gem_m4_strip_nstripsfired_neg, &b_bb_gem_m4_strip_nstripsfired_neg);
   fChain->SetBranchAddress("bb.gem.m5.clust.nclustu", &bb_gem_m5_clust_nclustu, &b_bb_gem_m5_clust_nclustu);
   fChain->SetBranchAddress("bb.gem.m5.clust.nclustu_neg", &bb_gem_m5_clust_nclustu_neg, &b_bb_gem_m5_clust_nclustu_neg);
   fChain->SetBranchAddress("bb.gem.m5.clust.nclustu_tot", &bb_gem_m5_clust_nclustu_tot, &b_bb_gem_m5_clust_nclustu_tot);
   fChain->SetBranchAddress("bb.gem.m5.clust.nclustv", &bb_gem_m5_clust_nclustv, &b_bb_gem_m5_clust_nclustv);
   fChain->SetBranchAddress("bb.gem.m5.clust.nclustv_neg", &bb_gem_m5_clust_nclustv_neg, &b_bb_gem_m5_clust_nclustv_neg);
   fChain->SetBranchAddress("bb.gem.m5.clust.nclustv_tot", &bb_gem_m5_clust_nclustv_tot, &b_bb_gem_m5_clust_nclustv_tot);
   fChain->SetBranchAddress("bb.gem.m5.strip.nstrips_keep", &bb_gem_m5_strip_nstrips_keep, &b_bb_gem_m5_strip_nstrips_keep);
   fChain->SetBranchAddress("bb.gem.m5.strip.nstrips_keepU", &bb_gem_m5_strip_nstrips_keepU, &b_bb_gem_m5_strip_nstrips_keepU);
   fChain->SetBranchAddress("bb.gem.m5.strip.nstrips_keepV", &bb_gem_m5_strip_nstrips_keepV, &b_bb_gem_m5_strip_nstrips_keepV);
   fChain->SetBranchAddress("bb.gem.m5.strip.nstrips_keep_lmax", &bb_gem_m5_strip_nstrips_keep_lmax, &b_bb_gem_m5_strip_nstrips_keep_lmax);
   fChain->SetBranchAddress("bb.gem.m5.strip.nstrips_keep_lmaxU", &bb_gem_m5_strip_nstrips_keep_lmaxU, &b_bb_gem_m5_strip_nstrips_keep_lmaxU);
   fChain->SetBranchAddress("bb.gem.m5.strip.nstrips_keep_lmaxV", &bb_gem_m5_strip_nstrips_keep_lmaxV, &b_bb_gem_m5_strip_nstrips_keep_lmaxV);
   fChain->SetBranchAddress("bb.gem.m5.strip.nstripsfired", &bb_gem_m5_strip_nstripsfired, &b_bb_gem_m5_strip_nstripsfired);
   fChain->SetBranchAddress("bb.gem.m5.strip.nstripsfired_neg", &bb_gem_m5_strip_nstripsfired_neg, &b_bb_gem_m5_strip_nstripsfired_neg);
   fChain->SetBranchAddress("bb.gem.m6.clust.nclustu", &bb_gem_m6_clust_nclustu, &b_bb_gem_m6_clust_nclustu);
   fChain->SetBranchAddress("bb.gem.m6.clust.nclustu_neg", &bb_gem_m6_clust_nclustu_neg, &b_bb_gem_m6_clust_nclustu_neg);
   fChain->SetBranchAddress("bb.gem.m6.clust.nclustu_tot", &bb_gem_m6_clust_nclustu_tot, &b_bb_gem_m6_clust_nclustu_tot);
   fChain->SetBranchAddress("bb.gem.m6.clust.nclustv", &bb_gem_m6_clust_nclustv, &b_bb_gem_m6_clust_nclustv);
   fChain->SetBranchAddress("bb.gem.m6.clust.nclustv_neg", &bb_gem_m6_clust_nclustv_neg, &b_bb_gem_m6_clust_nclustv_neg);
   fChain->SetBranchAddress("bb.gem.m6.clust.nclustv_tot", &bb_gem_m6_clust_nclustv_tot, &b_bb_gem_m6_clust_nclustv_tot);
   fChain->SetBranchAddress("bb.gem.m6.strip.nstrips_keep", &bb_gem_m6_strip_nstrips_keep, &b_bb_gem_m6_strip_nstrips_keep);
   fChain->SetBranchAddress("bb.gem.m6.strip.nstrips_keepU", &bb_gem_m6_strip_nstrips_keepU, &b_bb_gem_m6_strip_nstrips_keepU);
   fChain->SetBranchAddress("bb.gem.m6.strip.nstrips_keepV", &bb_gem_m6_strip_nstrips_keepV, &b_bb_gem_m6_strip_nstrips_keepV);
   fChain->SetBranchAddress("bb.gem.m6.strip.nstrips_keep_lmax", &bb_gem_m6_strip_nstrips_keep_lmax, &b_bb_gem_m6_strip_nstrips_keep_lmax);
   fChain->SetBranchAddress("bb.gem.m6.strip.nstrips_keep_lmaxU", &bb_gem_m6_strip_nstrips_keep_lmaxU, &b_bb_gem_m6_strip_nstrips_keep_lmaxU);
   fChain->SetBranchAddress("bb.gem.m6.strip.nstrips_keep_lmaxV", &bb_gem_m6_strip_nstrips_keep_lmaxV, &b_bb_gem_m6_strip_nstrips_keep_lmaxV);
   fChain->SetBranchAddress("bb.gem.m6.strip.nstripsfired", &bb_gem_m6_strip_nstripsfired, &b_bb_gem_m6_strip_nstripsfired);
   fChain->SetBranchAddress("bb.gem.m6.strip.nstripsfired_neg", &bb_gem_m6_strip_nstripsfired_neg, &b_bb_gem_m6_strip_nstripsfired_neg);
   fChain->SetBranchAddress("bb.gem.m7.clust.nclustu", &bb_gem_m7_clust_nclustu, &b_bb_gem_m7_clust_nclustu);
   fChain->SetBranchAddress("bb.gem.m7.clust.nclustu_neg", &bb_gem_m7_clust_nclustu_neg, &b_bb_gem_m7_clust_nclustu_neg);
   fChain->SetBranchAddress("bb.gem.m7.clust.nclustu_tot", &bb_gem_m7_clust_nclustu_tot, &b_bb_gem_m7_clust_nclustu_tot);
   fChain->SetBranchAddress("bb.gem.m7.clust.nclustv", &bb_gem_m7_clust_nclustv, &b_bb_gem_m7_clust_nclustv);
   fChain->SetBranchAddress("bb.gem.m7.clust.nclustv_neg", &bb_gem_m7_clust_nclustv_neg, &b_bb_gem_m7_clust_nclustv_neg);
   fChain->SetBranchAddress("bb.gem.m7.clust.nclustv_tot", &bb_gem_m7_clust_nclustv_tot, &b_bb_gem_m7_clust_nclustv_tot);
   fChain->SetBranchAddress("bb.gem.m7.strip.nstrips_keep", &bb_gem_m7_strip_nstrips_keep, &b_bb_gem_m7_strip_nstrips_keep);
   fChain->SetBranchAddress("bb.gem.m7.strip.nstrips_keepU", &bb_gem_m7_strip_nstrips_keepU, &b_bb_gem_m7_strip_nstrips_keepU);
   fChain->SetBranchAddress("bb.gem.m7.strip.nstrips_keepV", &bb_gem_m7_strip_nstrips_keepV, &b_bb_gem_m7_strip_nstrips_keepV);
   fChain->SetBranchAddress("bb.gem.m7.strip.nstrips_keep_lmax", &bb_gem_m7_strip_nstrips_keep_lmax, &b_bb_gem_m7_strip_nstrips_keep_lmax);
   fChain->SetBranchAddress("bb.gem.m7.strip.nstrips_keep_lmaxU", &bb_gem_m7_strip_nstrips_keep_lmaxU, &b_bb_gem_m7_strip_nstrips_keep_lmaxU);
   fChain->SetBranchAddress("bb.gem.m7.strip.nstrips_keep_lmaxV", &bb_gem_m7_strip_nstrips_keep_lmaxV, &b_bb_gem_m7_strip_nstrips_keep_lmaxV);
   fChain->SetBranchAddress("bb.gem.m7.strip.nstripsfired", &bb_gem_m7_strip_nstripsfired, &b_bb_gem_m7_strip_nstripsfired);
   fChain->SetBranchAddress("bb.gem.m7.strip.nstripsfired_neg", &bb_gem_m7_strip_nstripsfired_neg, &b_bb_gem_m7_strip_nstripsfired_neg);
   fChain->SetBranchAddress("bb.gem.nlayershit", &bb_gem_nlayershit, &b_bb_gem_nlayershit);
   fChain->SetBranchAddress("bb.gem.nlayershitu", &bb_gem_nlayershitu, &b_bb_gem_nlayershitu);
   fChain->SetBranchAddress("bb.gem.nlayershituv", &bb_gem_nlayershituv, &b_bb_gem_nlayershituv);
   fChain->SetBranchAddress("bb.gem.nlayershitv", &bb_gem_nlayershitv, &b_bb_gem_nlayershitv);
   fChain->SetBranchAddress("bb.gem.track.besttrack", &bb_gem_track_besttrack, &b_bb_gem_track_besttrack);
   fChain->SetBranchAddress("bb.gem.track.ntrack", &bb_gem_track_ntrack, &b_bb_gem_track_ntrack);
   fChain->SetBranchAddress("bb.gem.trigtime", &bb_gem_trigtime, &b_bb_gem_trigtime);
   fChain->SetBranchAddress("bb.grinch_tdc.bestcluster", &bb_grinch_tdc_bestcluster, &b_bb_grinch_tdc_bestcluster);
   fChain->SetBranchAddress("bb.grinch_tdc.clus.adc", &bb_grinch_tdc_clus_adc, &b_bb_grinch_tdc_clus_adc);
   fChain->SetBranchAddress("bb.grinch_tdc.clus.dx", &bb_grinch_tdc_clus_dx, &b_bb_grinch_tdc_clus_dx);
   fChain->SetBranchAddress("bb.grinch_tdc.clus.dy", &bb_grinch_tdc_clus_dy, &b_bb_grinch_tdc_clus_dy);
   fChain->SetBranchAddress("bb.grinch_tdc.clus.mirrorindex", &bb_grinch_tdc_clus_mirrorindex, &b_bb_grinch_tdc_clus_mirrorindex);
   fChain->SetBranchAddress("bb.grinch_tdc.clus.size", &bb_grinch_tdc_clus_size, &b_bb_grinch_tdc_clus_size);
   fChain->SetBranchAddress("bb.grinch_tdc.clus.t_mean", &bb_grinch_tdc_clus_t_mean, &b_bb_grinch_tdc_clus_t_mean);
   fChain->SetBranchAddress("bb.grinch_tdc.clus.t_rms", &bb_grinch_tdc_clus_t_rms, &b_bb_grinch_tdc_clus_t_rms);
   fChain->SetBranchAddress("bb.grinch_tdc.clus.tot_mean", &bb_grinch_tdc_clus_tot_mean, &b_bb_grinch_tdc_clus_tot_mean);
   fChain->SetBranchAddress("bb.grinch_tdc.clus.trackindex", &bb_grinch_tdc_clus_trackindex, &b_bb_grinch_tdc_clus_trackindex);
   fChain->SetBranchAddress("bb.grinch_tdc.clus.x_mean", &bb_grinch_tdc_clus_x_mean, &b_bb_grinch_tdc_clus_x_mean);
   fChain->SetBranchAddress("bb.grinch_tdc.clus.y_mean", &bb_grinch_tdc_clus_y_mean, &b_bb_grinch_tdc_clus_y_mean);
   fChain->SetBranchAddress("bb.grinch_tdc.nclus", &bb_grinch_tdc_nclus, &b_bb_grinch_tdc_nclus);
   fChain->SetBranchAddress("bb.grinch_tdc.ngoodhits", &bb_grinch_tdc_ngoodhits, &b_bb_grinch_tdc_ngoodhits);
   fChain->SetBranchAddress("bb.grinch_tdc.ntrackmatch", &bb_grinch_tdc_ntrackmatch, &b_bb_grinch_tdc_ntrackmatch);
   fChain->SetBranchAddress("bb.hodotdc.nclus", &bb_hodotdc_nclus, &b_bb_hodotdc_nclus);
   fChain->SetBranchAddress("bb.ps.againblk", &bb_ps_againblk, &b_bb_ps_againblk);
   fChain->SetBranchAddress("bb.ps.atimeblk", &bb_ps_atimeblk, &b_bb_ps_atimeblk);
   fChain->SetBranchAddress("bb.ps.colblk", &bb_ps_colblk, &b_bb_ps_colblk);
   fChain->SetBranchAddress("bb.ps.e", &bb_ps_e, &b_bb_ps_e);
   fChain->SetBranchAddress("bb.ps.eblk", &bb_ps_eblk, &b_bb_ps_eblk);
   fChain->SetBranchAddress("bb.ps.idblk", &bb_ps_idblk, &b_bb_ps_idblk);
   fChain->SetBranchAddress("bb.ps.index", &bb_ps_index, &b_bb_ps_index);
   fChain->SetBranchAddress("bb.ps.nblk", &bb_ps_nblk, &b_bb_ps_nblk);
   fChain->SetBranchAddress("bb.ps.nclus", &bb_ps_nclus, &b_bb_ps_nclus);
   fChain->SetBranchAddress("bb.ps.ngoodADChits", &bb_ps_ngoodADChits, &b_bb_ps_ngoodADChits);
   fChain->SetBranchAddress("bb.ps.rowblk", &bb_ps_rowblk, &b_bb_ps_rowblk);
   fChain->SetBranchAddress("bb.ps.x", &bb_ps_x, &b_bb_ps_x);
   fChain->SetBranchAddress("bb.ps.y", &bb_ps_y, &b_bb_ps_y);
   fChain->SetBranchAddress("bb.sh.againblk", &bb_sh_againblk, &b_bb_sh_againblk);
   fChain->SetBranchAddress("bb.sh.atimeblk", &bb_sh_atimeblk, &b_bb_sh_atimeblk);
   fChain->SetBranchAddress("bb.sh.colblk", &bb_sh_colblk, &b_bb_sh_colblk);
   fChain->SetBranchAddress("bb.sh.e", &bb_sh_e, &b_bb_sh_e);
   fChain->SetBranchAddress("bb.sh.eblk", &bb_sh_eblk, &b_bb_sh_eblk);
   fChain->SetBranchAddress("bb.sh.idblk", &bb_sh_idblk, &b_bb_sh_idblk);
   fChain->SetBranchAddress("bb.sh.index", &bb_sh_index, &b_bb_sh_index);
   fChain->SetBranchAddress("bb.sh.nblk", &bb_sh_nblk, &b_bb_sh_nblk);
   fChain->SetBranchAddress("bb.sh.nclus", &bb_sh_nclus, &b_bb_sh_nclus);
   fChain->SetBranchAddress("bb.sh.ngoodADChits", &bb_sh_ngoodADChits, &b_bb_sh_ngoodADChits);
   fChain->SetBranchAddress("bb.sh.rowblk", &bb_sh_rowblk, &b_bb_sh_rowblk);
   fChain->SetBranchAddress("bb.sh.x", &bb_sh_x, &b_bb_sh_x);
   fChain->SetBranchAddress("bb.sh.y", &bb_sh_y, &b_bb_sh_y);
   fChain->SetBranchAddress("bb.tr.n", &bb_tr_n, &b_bb_tr_n);
   fChain->SetBranchAddress("bb.ts.over_threshold", &bb_ts_over_threshold, &b_bb_ts_over_threshold);
   fChain->SetBranchAddress("e.kine.Q2", &e_kine_Q2, &b_e_kine_Q2);
   fChain->SetBranchAddress("e.kine.W2", &e_kine_W2, &b_e_kine_W2);
   fChain->SetBranchAddress("e.kine.angle", &e_kine_angle, &b_e_kine_angle);
   fChain->SetBranchAddress("e.kine.epsilon", &e_kine_epsilon, &b_e_kine_epsilon);
   fChain->SetBranchAddress("e.kine.nu", &e_kine_nu, &b_e_kine_nu);
   fChain->SetBranchAddress("e.kine.omega", &e_kine_omega, &b_e_kine_omega);
   fChain->SetBranchAddress("e.kine.ph_q", &e_kine_ph_q, &b_e_kine_ph_q);
   fChain->SetBranchAddress("e.kine.q3m", &e_kine_q3m, &b_e_kine_q3m);
   fChain->SetBranchAddress("e.kine.q_x", &e_kine_q_x, &b_e_kine_q_x);
   fChain->SetBranchAddress("e.kine.q_y", &e_kine_q_y, &b_e_kine_q_y);
   fChain->SetBranchAddress("e.kine.q_z", &e_kine_q_z, &b_e_kine_q_z);
   fChain->SetBranchAddress("e.kine.th_q", &e_kine_th_q, &b_e_kine_th_q);
   fChain->SetBranchAddress("e.kine.x_bj", &e_kine_x_bj, &b_e_kine_x_bj);
   fChain->SetBranchAddress("g.datatype", &g_datatype, &b_g_datatype);
   fChain->SetBranchAddress("g.evlen", &g_evlen, &b_g_evlen);
   fChain->SetBranchAddress("g.evnum", &g_evnum, &b_g_evnum);
   fChain->SetBranchAddress("g.evtime", &g_evtime, &b_g_evtime);
   fChain->SetBranchAddress("g.evtyp", &g_evtyp, &b_g_evtyp);
   fChain->SetBranchAddress("g.runnum", &g_runnum, &b_g_runnum);
   fChain->SetBranchAddress("g.runtime", &g_runtime, &b_g_runtime);
   fChain->SetBranchAddress("g.runtype", &g_runtype, &b_g_runtype);
   fChain->SetBranchAddress("g.trigbits", &g_trigbits, &b_g_trigbits);
   fChain->SetBranchAddress("sbs.HCALdir_x", &sbs_HCALdir_x, &b_sbs_HCALdir_x);
   fChain->SetBranchAddress("sbs.HCALdir_y", &sbs_HCALdir_y, &b_sbs_HCALdir_y);
   fChain->SetBranchAddress("sbs.HCALdir_z", &sbs_HCALdir_z, &b_sbs_HCALdir_z);
   fChain->SetBranchAddress("sbs.HCALph_n", &sbs_HCALph_n, &b_sbs_HCALph_n);
   fChain->SetBranchAddress("sbs.HCALth_n", &sbs_HCALth_n, &b_sbs_HCALth_n);
   fChain->SetBranchAddress("sbs.hcal.againblk", &sbs_hcal_againblk, &b_sbs_hcal_againblk);
   fChain->SetBranchAddress("sbs.hcal.atimeblk", &sbs_hcal_atimeblk, &b_sbs_hcal_atimeblk);
   fChain->SetBranchAddress("sbs.hcal.colblk", &sbs_hcal_colblk, &b_sbs_hcal_colblk);
   fChain->SetBranchAddress("sbs.hcal.e", &sbs_hcal_e, &b_sbs_hcal_e);
   fChain->SetBranchAddress("sbs.hcal.eblk", &sbs_hcal_eblk, &b_sbs_hcal_eblk);
   fChain->SetBranchAddress("sbs.hcal.idblk", &sbs_hcal_idblk, &b_sbs_hcal_idblk);
   fChain->SetBranchAddress("sbs.hcal.index", &sbs_hcal_index, &b_sbs_hcal_index);
   fChain->SetBranchAddress("sbs.hcal.nblk", &sbs_hcal_nblk, &b_sbs_hcal_nblk);
   fChain->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus, &b_sbs_hcal_nclus);
   fChain->SetBranchAddress("sbs.hcal.ngoodADChits", &sbs_hcal_ngoodADChits, &b_sbs_hcal_ngoodADChits);
   fChain->SetBranchAddress("sbs.hcal.rowblk", &sbs_hcal_rowblk, &b_sbs_hcal_rowblk);
   fChain->SetBranchAddress("sbs.hcal.tdctimeblk", &sbs_hcal_tdctimeblk, &b_sbs_hcal_tdctimeblk);
   fChain->SetBranchAddress("sbs.hcal.x", &sbs_hcal_x, &b_sbs_hcal_x);
   fChain->SetBranchAddress("sbs.hcal.y", &sbs_hcal_y, &b_sbs_hcal_y);
   fChain->SetBranchAddress("sbs.tr.n", &sbs_tr_n, &b_sbs_tr_n);
   fChain->SetBranchAddress("singletrack", &singletrack, &b_singletrack);
   fChain->SetBranchAddress("anytrack", &anytrack, &b_anytrack);
   fChain->SetBranchAddress("HALLA_p", &HALLA_p, &b_HALLA_p);
   fChain->SetBranchAddress("fEvtHdr.fEvtTime", &fEvtHdr_fEvtTime, &b_Event_Branch_fEvtHdr_fEvtTime);
   fChain->SetBranchAddress("fEvtHdr.fEvtNum", &fEvtHdr_fEvtNum, &b_Event_Branch_fEvtHdr_fEvtNum);
   fChain->SetBranchAddress("fEvtHdr.fEvtType", &fEvtHdr_fEvtType, &b_Event_Branch_fEvtHdr_fEvtType);
   fChain->SetBranchAddress("fEvtHdr.fEvtLen", &fEvtHdr_fEvtLen, &b_Event_Branch_fEvtHdr_fEvtLen);
   fChain->SetBranchAddress("fEvtHdr.fHelicity", &fEvtHdr_fHelicity, &b_Event_Branch_fEvtHdr_fHelicity);
   fChain->SetBranchAddress("fEvtHdr.fTrigBits", &fEvtHdr_fTrigBits, &b_Event_Branch_fEvtHdr_fTrigBits);
   fChain->SetBranchAddress("fEvtHdr.fRun", &fEvtHdr_fRun, &b_Event_Branch_fEvtHdr_fRun);
   Notify();
}

Bool_t gmn_tree_allstripoutputs::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gmn_tree_allstripoutputs::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gmn_tree_allstripoutputs::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gmn_tree_allstripoutputs_cxx
