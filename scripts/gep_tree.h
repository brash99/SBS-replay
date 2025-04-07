//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr  2 10:09:18 2025 by ROOT version 6.30/04
// from TChain T/
//////////////////////////////////////////////////////////

#ifndef gep_tree_h
#define gep_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "THaEvent.h"
#include "THaEvent.h"

class gep_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Ndata_MC_bbgemhit_ptridx;
   Double_t        MC_bbgemhit_ptridx[1];   //[Ndata.MC.bbgemhit_ptridx]
   Int_t           Ndata_MC_bbgemhit_sdtridx;
   Double_t        MC_bbgemhit_sdtridx[1];   //[Ndata.MC.bbgemhit_sdtridx]
   Int_t           Ndata_MC_bbgemtrack_ptridx;
   Double_t        MC_bbgemtrack_ptridx[1];   //[Ndata.MC.bbgemtrack_ptridx]
   Int_t           Ndata_MC_bbgemtrack_sdtridx;
   Double_t        MC_bbgemtrack_sdtridx[1];   //[Ndata.MC.bbgemtrack_sdtridx]
   Int_t           Ndata_MC_bbhodohit_ptridx;
   Double_t        MC_bbhodohit_ptridx[1];   //[Ndata.MC.bbhodohit_ptridx]
   Int_t           Ndata_MC_bbhodohit_sdtridx;
   Double_t        MC_bbhodohit_sdtridx[1];   //[Ndata.MC.bbhodohit_sdtridx]
   Int_t           Ndata_MC_bbpshit_ptridx;
   Double_t        MC_bbpshit_ptridx[1];   //[Ndata.MC.bbpshit_ptridx]
   Int_t           Ndata_MC_bbpshit_sdtridx;
   Double_t        MC_bbpshit_sdtridx[1];   //[Ndata.MC.bbpshit_sdtridx]
   Int_t           Ndata_MC_bbshhit_ptridx;
   Double_t        MC_bbshhit_ptridx[1];   //[Ndata.MC.bbshhit_ptridx]
   Int_t           Ndata_MC_bbshhit_sdtridx;
   Double_t        MC_bbshhit_sdtridx[1];   //[Ndata.MC.bbshhit_sdtridx]
   Int_t           Ndata_MC_hcalhit_ptridx;
   Double_t        MC_hcalhit_ptridx[1];   //[Ndata.MC.hcalhit_ptridx]
   Int_t           Ndata_MC_hcalhit_sdtridx;
   Double_t        MC_hcalhit_sdtridx[1];   //[Ndata.MC.hcalhit_sdtridx]
   Int_t           Ndata_MC_ptrack_etot;
   Double_t        MC_ptrack_etot[1];   //[Ndata.MC.ptrack_etot]
   Int_t           Ndata_MC_ptrack_momx;
   Double_t        MC_ptrack_momx[1];   //[Ndata.MC.ptrack_momx]
   Int_t           Ndata_MC_ptrack_momy;
   Double_t        MC_ptrack_momy[1];   //[Ndata.MC.ptrack_momy]
   Int_t           Ndata_MC_ptrack_momz;
   Double_t        MC_ptrack_momz[1];   //[Ndata.MC.ptrack_momz]
   Int_t           Ndata_MC_ptrack_pid;
   Double_t        MC_ptrack_pid[1];   //[Ndata.MC.ptrack_pid]
   Int_t           Ndata_MC_ptrack_polx;
   Double_t        MC_ptrack_polx[1];   //[Ndata.MC.ptrack_polx]
   Int_t           Ndata_MC_ptrack_poly;
   Double_t        MC_ptrack_poly[1];   //[Ndata.MC.ptrack_poly]
   Int_t           Ndata_MC_ptrack_polz;
   Double_t        MC_ptrack_polz[1];   //[Ndata.MC.ptrack_polz]
   Int_t           Ndata_MC_ptrack_posx;
   Double_t        MC_ptrack_posx[1];   //[Ndata.MC.ptrack_posx]
   Int_t           Ndata_MC_ptrack_posy;
   Double_t        MC_ptrack_posy[1];   //[Ndata.MC.ptrack_posy]
   Int_t           Ndata_MC_ptrack_posz;
   Double_t        MC_ptrack_posz[1];   //[Ndata.MC.ptrack_posz]
   Int_t           Ndata_MC_ptrack_t;
   Double_t        MC_ptrack_t[1];   //[Ndata.MC.ptrack_t]
   Int_t           Ndata_MC_ptrack_tid;
   Double_t        MC_ptrack_tid[1];   //[Ndata.MC.ptrack_tid]
   Int_t           Ndata_MC_sdtrack_etot;
   Double_t        MC_sdtrack_etot[1];   //[Ndata.MC.sdtrack_etot]
   Int_t           Ndata_MC_sdtrack_mid;
   Double_t        MC_sdtrack_mid[1];   //[Ndata.MC.sdtrack_mid]
   Int_t           Ndata_MC_sdtrack_momx;
   Double_t        MC_sdtrack_momx[1];   //[Ndata.MC.sdtrack_momx]
   Int_t           Ndata_MC_sdtrack_momy;
   Double_t        MC_sdtrack_momy[1];   //[Ndata.MC.sdtrack_momy]
   Int_t           Ndata_MC_sdtrack_momz;
   Double_t        MC_sdtrack_momz[1];   //[Ndata.MC.sdtrack_momz]
   Int_t           Ndata_MC_sdtrack_pid;
   Double_t        MC_sdtrack_pid[1];   //[Ndata.MC.sdtrack_pid]
   Int_t           Ndata_MC_sdtrack_polx;
   Double_t        MC_sdtrack_polx[1];   //[Ndata.MC.sdtrack_polx]
   Int_t           Ndata_MC_sdtrack_poly;
   Double_t        MC_sdtrack_poly[1];   //[Ndata.MC.sdtrack_poly]
   Int_t           Ndata_MC_sdtrack_polz;
   Double_t        MC_sdtrack_polz[1];   //[Ndata.MC.sdtrack_polz]
   Int_t           Ndata_MC_sdtrack_posx;
   Double_t        MC_sdtrack_posx[1];   //[Ndata.MC.sdtrack_posx]
   Int_t           Ndata_MC_sdtrack_posy;
   Double_t        MC_sdtrack_posy[1];   //[Ndata.MC.sdtrack_posy]
   Int_t           Ndata_MC_sdtrack_posz;
   Double_t        MC_sdtrack_posz[1];   //[Ndata.MC.sdtrack_posz]
   Int_t           Ndata_MC_sdtrack_t;
   Double_t        MC_sdtrack_t[1];   //[Ndata.MC.sdtrack_t]
   Int_t           Ndata_MC_sdtrack_tid;
   Double_t        MC_sdtrack_tid[1];   //[Ndata.MC.sdtrack_tid]
   Int_t           Ndata_MC_sdtrack_vEkin;
   Double_t        MC_sdtrack_vEkin[1];   //[Ndata.MC.sdtrack_vEkin]
   Int_t           Ndata_MC_sdtrack_vnx;
   Double_t        MC_sdtrack_vnx[1];   //[Ndata.MC.sdtrack_vnx]
   Int_t           Ndata_MC_sdtrack_vny;
   Double_t        MC_sdtrack_vny[1];   //[Ndata.MC.sdtrack_vny]
   Int_t           Ndata_MC_sdtrack_vnz;
   Double_t        MC_sdtrack_vnz[1];   //[Ndata.MC.sdtrack_vnz]
   Int_t           Ndata_MC_sdtrack_vx;
   Double_t        MC_sdtrack_vx[1];   //[Ndata.MC.sdtrack_vx]
   Int_t           Ndata_MC_sdtrack_vy;
   Double_t        MC_sdtrack_vy[1];   //[Ndata.MC.sdtrack_vy]
   Int_t           Ndata_MC_sdtrack_vz;
   Double_t        MC_sdtrack_vz[1];   //[Ndata.MC.sdtrack_vz]
   Int_t           Ndata_earm_ecal_a_amp_p;
   Double_t        earm_ecal_a_amp_p[26];   //[Ndata.earm.ecal.a_amp_p]
   Int_t           Ndata_earm_ecal_a_p;
   Double_t        earm_ecal_a_p[26];   //[Ndata.earm.ecal.a_p]
   Int_t           Ndata_earm_ecal_a_time;
   Double_t        earm_ecal_a_time[26];   //[Ndata.earm.ecal.a_time]
   Int_t           Ndata_earm_ecal_adccol;
   Double_t        earm_ecal_adccol[26];   //[Ndata.earm.ecal.adccol]
   Int_t           Ndata_earm_ecal_adcelemID;
   Double_t        earm_ecal_adcelemID[26];   //[Ndata.earm.ecal.adcelemID]
   Int_t           Ndata_earm_ecal_adcrow;
   Double_t        earm_ecal_adcrow[26];   //[Ndata.earm.ecal.adcrow]
   Int_t           Ndata_earm_ecal_adcxpos;
   Double_t        earm_ecal_adcxpos[26];   //[Ndata.earm.ecal.adcxpos]
   Int_t           Ndata_earm_ecal_adcypos;
   Double_t        earm_ecal_adcypos[26];   //[Ndata.earm.ecal.adcypos]
   Int_t           Ndata_earm_ecal_clus_adctime;
   Double_t        earm_ecal_clus_adctime[7];   //[Ndata.earm.ecal.clus.adctime]
   Int_t           Ndata_earm_ecal_clus_again;
   Double_t        earm_ecal_clus_again[7];   //[Ndata.earm.ecal.clus.again]
   Int_t           Ndata_earm_ecal_clus_atimeblk;
   Double_t        earm_ecal_clus_atimeblk[7];   //[Ndata.earm.ecal.clus.atimeblk]
   Int_t           Ndata_earm_ecal_clus_col;
   Double_t        earm_ecal_clus_col[7];   //[Ndata.earm.ecal.clus.col]
   Int_t           Ndata_earm_ecal_clus_col_goodtdc;
   Double_t        earm_ecal_clus_col_goodtdc[7];   //[Ndata.earm.ecal.clus.col_goodtdc]
   Int_t           Ndata_earm_ecal_clus_e;
   Double_t        earm_ecal_clus_e[7];   //[Ndata.earm.ecal.clus.e]
   Int_t           Ndata_earm_ecal_clus_e_goodtdc;
   Double_t        earm_ecal_clus_e_goodtdc[7];   //[Ndata.earm.ecal.clus.e_goodtdc]
   Int_t           Ndata_earm_ecal_clus_eblk;
   Double_t        earm_ecal_clus_eblk[7];   //[Ndata.earm.ecal.clus.eblk]
   Int_t           Ndata_earm_ecal_clus_eblk_goodtdc;
   Double_t        earm_ecal_clus_eblk_goodtdc[7];   //[Ndata.earm.ecal.clus.eblk_goodtdc]
   Int_t           Ndata_earm_ecal_clus_id;
   Double_t        earm_ecal_clus_id[7];   //[Ndata.earm.ecal.clus.id]
   Int_t           Ndata_earm_ecal_clus_id_goodtdc;
   Double_t        earm_ecal_clus_id_goodtdc[7];   //[Ndata.earm.ecal.clus.id_goodtdc]
   Int_t           Ndata_earm_ecal_clus_nblk;
   Double_t        earm_ecal_clus_nblk[7];   //[Ndata.earm.ecal.clus.nblk]
   Int_t           Ndata_earm_ecal_clus_nblk_goodtdc;
   Double_t        earm_ecal_clus_nblk_goodtdc[7];   //[Ndata.earm.ecal.clus.nblk_goodtdc]
   Int_t           Ndata_earm_ecal_clus_row;
   Double_t        earm_ecal_clus_row[7];   //[Ndata.earm.ecal.clus.row]
   Int_t           Ndata_earm_ecal_clus_row_goodtdc;
   Double_t        earm_ecal_clus_row_goodtdc[7];   //[Ndata.earm.ecal.clus.row_goodtdc]
   Int_t           Ndata_earm_ecal_clus_tdctime;
   Double_t        earm_ecal_clus_tdctime[7];   //[Ndata.earm.ecal.clus.tdctime]
   Int_t           Ndata_earm_ecal_clus_tdctime_tw;
   Double_t        earm_ecal_clus_tdctime_tw[7];   //[Ndata.earm.ecal.clus.tdctime_tw]
   Int_t           Ndata_earm_ecal_clus_tdctimeblk;
   Double_t        earm_ecal_clus_tdctimeblk[7];   //[Ndata.earm.ecal.clus.tdctimeblk]
   Int_t           Ndata_earm_ecal_clus_tdctimeblk_tw;
   Double_t        earm_ecal_clus_tdctimeblk_tw[7];   //[Ndata.earm.ecal.clus.tdctimeblk_tw]
   Int_t           Ndata_earm_ecal_clus_x;
   Double_t        earm_ecal_clus_x[7];   //[Ndata.earm.ecal.clus.x]
   Int_t           Ndata_earm_ecal_clus_y;
   Double_t        earm_ecal_clus_y[7];   //[Ndata.earm.ecal.clus.y]
   Int_t           Ndata_earm_ecal_clus_blk_again;
   Double_t        earm_ecal_clus_blk_again[21];   //[Ndata.earm.ecal.clus_blk.again]
   Int_t           Ndata_earm_ecal_clus_blk_atime;
   Double_t        earm_ecal_clus_blk_atime[21];   //[Ndata.earm.ecal.clus_blk.atime]
   Int_t           Ndata_earm_ecal_clus_blk_col;
   Double_t        earm_ecal_clus_blk_col[21];   //[Ndata.earm.ecal.clus_blk.col]
   Int_t           Ndata_earm_ecal_clus_blk_e;
   Double_t        earm_ecal_clus_blk_e[21];   //[Ndata.earm.ecal.clus_blk.e]
   Int_t           Ndata_earm_ecal_clus_blk_id;
   Double_t        earm_ecal_clus_blk_id[21];   //[Ndata.earm.ecal.clus_blk.id]
   Int_t           Ndata_earm_ecal_clus_blk_row;
   Double_t        earm_ecal_clus_blk_row[21];   //[Ndata.earm.ecal.clus_blk.row]
   Int_t           Ndata_earm_ecal_clus_blk_tdctime;
   Double_t        earm_ecal_clus_blk_tdctime[21];   //[Ndata.earm.ecal.clus_blk.tdctime]
   Int_t           Ndata_earm_ecal_clus_blk_tdctime_tw;
   Double_t        earm_ecal_clus_blk_tdctime_tw[21];   //[Ndata.earm.ecal.clus_blk.tdctime_tw]
   Int_t           Ndata_earm_ecal_clus_blk_x;
   Double_t        earm_ecal_clus_blk_x[21];   //[Ndata.earm.ecal.clus_blk.x]
   Int_t           Ndata_earm_ecal_clus_blk_y;
   Double_t        earm_ecal_clus_blk_y[21];   //[Ndata.earm.ecal.clus_blk.y]
   Int_t           Ndata_earm_ecal_ped;
   Double_t        earm_ecal_ped[26];   //[Ndata.earm.ecal.ped]
   Int_t           Ndata_earm_tr_beta;
   Double_t        earm_tr_beta[1];   //[Ndata.earm.tr.beta]
   Int_t           Ndata_earm_tr_chi2;
   Double_t        earm_tr_chi2[1];   //[Ndata.earm.tr.chi2]
   Int_t           Ndata_earm_tr_d_ph;
   Double_t        earm_tr_d_ph[1];   //[Ndata.earm.tr.d_ph]
   Int_t           Ndata_earm_tr_d_th;
   Double_t        earm_tr_d_th[1];   //[Ndata.earm.tr.d_th]
   Int_t           Ndata_earm_tr_d_x;
   Double_t        earm_tr_d_x[1];   //[Ndata.earm.tr.d_x]
   Int_t           Ndata_earm_tr_d_y;
   Double_t        earm_tr_d_y[1];   //[Ndata.earm.tr.d_y]
   Int_t           Ndata_earm_tr_dbeta;
   Double_t        earm_tr_dbeta[1];   //[Ndata.earm.tr.dbeta]
   Int_t           Ndata_earm_tr_dtime;
   Double_t        earm_tr_dtime[1];   //[Ndata.earm.tr.dtime]
   Int_t           Ndata_earm_tr_flag;
   Double_t        earm_tr_flag[1];   //[Ndata.earm.tr.flag]
   Int_t           Ndata_earm_tr_ndof;
   Double_t        earm_tr_ndof[1];   //[Ndata.earm.tr.ndof]
   Int_t           Ndata_earm_tr_p;
   Double_t        earm_tr_p[1];   //[Ndata.earm.tr.p]
   Int_t           Ndata_earm_tr_pathl;
   Double_t        earm_tr_pathl[1];   //[Ndata.earm.tr.pathl]
   Int_t           Ndata_earm_tr_ph;
   Double_t        earm_tr_ph[1];   //[Ndata.earm.tr.ph]
   Int_t           Ndata_earm_tr_px;
   Double_t        earm_tr_px[1];   //[Ndata.earm.tr.px]
   Int_t           Ndata_earm_tr_py;
   Double_t        earm_tr_py[1];   //[Ndata.earm.tr.py]
   Int_t           Ndata_earm_tr_pz;
   Double_t        earm_tr_pz[1];   //[Ndata.earm.tr.pz]
   Int_t           Ndata_earm_tr_r_ph;
   Double_t        earm_tr_r_ph[1];   //[Ndata.earm.tr.r_ph]
   Int_t           Ndata_earm_tr_r_th;
   Double_t        earm_tr_r_th[1];   //[Ndata.earm.tr.r_th]
   Int_t           Ndata_earm_tr_r_x;
   Double_t        earm_tr_r_x[1];   //[Ndata.earm.tr.r_x]
   Int_t           Ndata_earm_tr_r_y;
   Double_t        earm_tr_r_y[1];   //[Ndata.earm.tr.r_y]
   Int_t           Ndata_earm_tr_tg_dp;
   Double_t        earm_tr_tg_dp[1];   //[Ndata.earm.tr.tg_dp]
   Int_t           Ndata_earm_tr_tg_ph;
   Double_t        earm_tr_tg_ph[1];   //[Ndata.earm.tr.tg_ph]
   Int_t           Ndata_earm_tr_tg_th;
   Double_t        earm_tr_tg_th[1];   //[Ndata.earm.tr.tg_th]
   Int_t           Ndata_earm_tr_tg_y;
   Double_t        earm_tr_tg_y[1];   //[Ndata.earm.tr.tg_y]
   Int_t           Ndata_earm_tr_th;
   Double_t        earm_tr_th[1];   //[Ndata.earm.tr.th]
   Int_t           Ndata_earm_tr_time;
   Double_t        earm_tr_time[1];   //[Ndata.earm.tr.time]
   Int_t           Ndata_earm_tr_vx;
   Double_t        earm_tr_vx[1];   //[Ndata.earm.tr.vx]
   Int_t           Ndata_earm_tr_vy;
   Double_t        earm_tr_vy[1];   //[Ndata.earm.tr.vy]
   Int_t           Ndata_earm_tr_vz;
   Double_t        earm_tr_vz[1];   //[Ndata.earm.tr.vz]
   Int_t           Ndata_earm_tr_x;
   Double_t        earm_tr_x[1];   //[Ndata.earm.tr.x]
   Int_t           Ndata_earm_tr_y;
   Double_t        earm_tr_y[1];   //[Ndata.earm.tr.y]
   Int_t           Ndata_sbs_gemFPP_hit_ADCU;
   Double_t        sbs_gemFPP_hit_ADCU[11];   //[Ndata.sbs.gemFPP.hit.ADCU]
   Int_t           Ndata_sbs_gemFPP_hit_ADCU_deconv;
   Double_t        sbs_gemFPP_hit_ADCU_deconv[11];   //[Ndata.sbs.gemFPP.hit.ADCU_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_ADCV;
   Double_t        sbs_gemFPP_hit_ADCV[11];   //[Ndata.sbs.gemFPP.hit.ADCV]
   Int_t           Ndata_sbs_gemFPP_hit_ADCV_deconv;
   Double_t        sbs_gemFPP_hit_ADCV_deconv[11];   //[Ndata.sbs.gemFPP.hit.ADCV_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_ADCasym;
   Double_t        sbs_gemFPP_hit_ADCasym[11];   //[Ndata.sbs.gemFPP.hit.ADCasym]
   Int_t           Ndata_sbs_gemFPP_hit_ADCasym_deconv;
   Double_t        sbs_gemFPP_hit_ADCasym_deconv[11];   //[Ndata.sbs.gemFPP.hit.ADCasym_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_ADCavg;
   Double_t        sbs_gemFPP_hit_ADCavg[11];   //[Ndata.sbs.gemFPP.hit.ADCavg]
   Int_t           Ndata_sbs_gemFPP_hit_ADCavg_deconv;
   Double_t        sbs_gemFPP_hit_ADCavg_deconv[11];   //[Ndata.sbs.gemFPP.hit.ADCavg_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac0_Umax;
   Double_t        sbs_gemFPP_hit_ADCfrac0_Umax[11];   //[Ndata.sbs.gemFPP.hit.ADCfrac0_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac0_Vmax;
   Double_t        sbs_gemFPP_hit_ADCfrac0_Vmax[11];   //[Ndata.sbs.gemFPP.hit.ADCfrac0_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac1_Umax;
   Double_t        sbs_gemFPP_hit_ADCfrac1_Umax[11];   //[Ndata.sbs.gemFPP.hit.ADCfrac1_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac1_Vmax;
   Double_t        sbs_gemFPP_hit_ADCfrac1_Vmax[11];   //[Ndata.sbs.gemFPP.hit.ADCfrac1_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac2_Umax;
   Double_t        sbs_gemFPP_hit_ADCfrac2_Umax[11];   //[Ndata.sbs.gemFPP.hit.ADCfrac2_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac2_Vmax;
   Double_t        sbs_gemFPP_hit_ADCfrac2_Vmax[11];   //[Ndata.sbs.gemFPP.hit.ADCfrac2_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac3_Umax;
   Double_t        sbs_gemFPP_hit_ADCfrac3_Umax[11];   //[Ndata.sbs.gemFPP.hit.ADCfrac3_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac3_Vmax;
   Double_t        sbs_gemFPP_hit_ADCfrac3_Vmax[11];   //[Ndata.sbs.gemFPP.hit.ADCfrac3_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac4_Umax;
   Double_t        sbs_gemFPP_hit_ADCfrac4_Umax[11];   //[Ndata.sbs.gemFPP.hit.ADCfrac4_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac4_Vmax;
   Double_t        sbs_gemFPP_hit_ADCfrac4_Vmax[11];   //[Ndata.sbs.gemFPP.hit.ADCfrac4_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac5_Umax;
   Double_t        sbs_gemFPP_hit_ADCfrac5_Umax[11];   //[Ndata.sbs.gemFPP.hit.ADCfrac5_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCfrac5_Vmax;
   Double_t        sbs_gemFPP_hit_ADCfrac5_Vmax[11];   //[Ndata.sbs.gemFPP.hit.ADCfrac5_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxcomboUclust_deconv;
   Double_t        sbs_gemFPP_hit_ADCmaxcomboUclust_deconv[11];   //[Ndata.sbs.gemFPP.hit.ADCmaxcomboUclust_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxcomboVclust_deconv;
   Double_t        sbs_gemFPP_hit_ADCmaxcomboVclust_deconv[11];   //[Ndata.sbs.gemFPP.hit.ADCmaxcomboVclust_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxsampU;
   Double_t        sbs_gemFPP_hit_ADCmaxsampU[11];   //[Ndata.sbs.gemFPP.hit.ADCmaxsampU]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxsampUclust;
   Double_t        sbs_gemFPP_hit_ADCmaxsampUclust[11];   //[Ndata.sbs.gemFPP.hit.ADCmaxsampUclust]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxsampUclust_deconv;
   Double_t        sbs_gemFPP_hit_ADCmaxsampUclust_deconv[11];   //[Ndata.sbs.gemFPP.hit.ADCmaxsampUclust_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxsampV;
   Double_t        sbs_gemFPP_hit_ADCmaxsampV[11];   //[Ndata.sbs.gemFPP.hit.ADCmaxsampV]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxsampVclust;
   Double_t        sbs_gemFPP_hit_ADCmaxsampVclust[11];   //[Ndata.sbs.gemFPP.hit.ADCmaxsampVclust]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxsampVclust_deconv;
   Double_t        sbs_gemFPP_hit_ADCmaxsampVclust_deconv[11];   //[Ndata.sbs.gemFPP.hit.ADCmaxsampVclust_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxstripU;
   Double_t        sbs_gemFPP_hit_ADCmaxstripU[11];   //[Ndata.sbs.gemFPP.hit.ADCmaxstripU]
   Int_t           Ndata_sbs_gemFPP_hit_ADCmaxstripV;
   Double_t        sbs_gemFPP_hit_ADCmaxstripV[11];   //[Ndata.sbs.gemFPP.hit.ADCmaxstripV]
   Int_t           Ndata_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_U;
   Double_t        sbs_gemFPP_hit_BUILD_ALL_SAMPLES_U[11];   //[Ndata.sbs.gemFPP.hit.BUILD_ALL_SAMPLES_U]
   Int_t           Ndata_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_V;
   Double_t        sbs_gemFPP_hit_BUILD_ALL_SAMPLES_V[11];   //[Ndata.sbs.gemFPP.hit.BUILD_ALL_SAMPLES_V]
   Int_t           Ndata_sbs_gemFPP_hit_CM_GOOD_U;
   Double_t        sbs_gemFPP_hit_CM_GOOD_U[11];   //[Ndata.sbs.gemFPP.hit.CM_GOOD_U]
   Int_t           Ndata_sbs_gemFPP_hit_CM_GOOD_V;
   Double_t        sbs_gemFPP_hit_CM_GOOD_V[11];   //[Ndata.sbs.gemFPP.hit.CM_GOOD_V]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC0_Umax;
   Double_t        sbs_gemFPP_hit_DeconvADC0_Umax[11];   //[Ndata.sbs.gemFPP.hit.DeconvADC0_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC0_Vmax;
   Double_t        sbs_gemFPP_hit_DeconvADC0_Vmax[11];   //[Ndata.sbs.gemFPP.hit.DeconvADC0_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC1_Umax;
   Double_t        sbs_gemFPP_hit_DeconvADC1_Umax[11];   //[Ndata.sbs.gemFPP.hit.DeconvADC1_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC1_Vmax;
   Double_t        sbs_gemFPP_hit_DeconvADC1_Vmax[11];   //[Ndata.sbs.gemFPP.hit.DeconvADC1_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC2_Umax;
   Double_t        sbs_gemFPP_hit_DeconvADC2_Umax[11];   //[Ndata.sbs.gemFPP.hit.DeconvADC2_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC2_Vmax;
   Double_t        sbs_gemFPP_hit_DeconvADC2_Vmax[11];   //[Ndata.sbs.gemFPP.hit.DeconvADC2_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC3_Umax;
   Double_t        sbs_gemFPP_hit_DeconvADC3_Umax[11];   //[Ndata.sbs.gemFPP.hit.DeconvADC3_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC3_Vmax;
   Double_t        sbs_gemFPP_hit_DeconvADC3_Vmax[11];   //[Ndata.sbs.gemFPP.hit.DeconvADC3_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC4_Umax;
   Double_t        sbs_gemFPP_hit_DeconvADC4_Umax[11];   //[Ndata.sbs.gemFPP.hit.DeconvADC4_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC4_Vmax;
   Double_t        sbs_gemFPP_hit_DeconvADC4_Vmax[11];   //[Ndata.sbs.gemFPP.hit.DeconvADC4_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC5_Umax;
   Double_t        sbs_gemFPP_hit_DeconvADC5_Umax[11];   //[Ndata.sbs.gemFPP.hit.DeconvADC5_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADC5_Vmax;
   Double_t        sbs_gemFPP_hit_DeconvADC5_Vmax[11];   //[Ndata.sbs.gemFPP.hit.DeconvADC5_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADCmaxcomboU;
   Double_t        sbs_gemFPP_hit_DeconvADCmaxcomboU[11];   //[Ndata.sbs.gemFPP.hit.DeconvADCmaxcomboU]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADCmaxcomboV;
   Double_t        sbs_gemFPP_hit_DeconvADCmaxcomboV[11];   //[Ndata.sbs.gemFPP.hit.DeconvADCmaxcomboV]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADCmaxsampU;
   Double_t        sbs_gemFPP_hit_DeconvADCmaxsampU[11];   //[Ndata.sbs.gemFPP.hit.DeconvADCmaxsampU]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADCmaxsampV;
   Double_t        sbs_gemFPP_hit_DeconvADCmaxsampV[11];   //[Ndata.sbs.gemFPP.hit.DeconvADCmaxsampV]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADCmaxstripU;
   Double_t        sbs_gemFPP_hit_DeconvADCmaxstripU[11];   //[Ndata.sbs.gemFPP.hit.DeconvADCmaxstripU]
   Int_t           Ndata_sbs_gemFPP_hit_DeconvADCmaxstripV;
   Double_t        sbs_gemFPP_hit_DeconvADCmaxstripV[11];   //[Ndata.sbs.gemFPP.hit.DeconvADCmaxstripV]
   Int_t           Ndata_sbs_gemFPP_hit_ENABLE_CM_U;
   Double_t        sbs_gemFPP_hit_ENABLE_CM_U[11];   //[Ndata.sbs.gemFPP.hit.ENABLE_CM_U]
   Int_t           Ndata_sbs_gemFPP_hit_ENABLE_CM_V;
   Double_t        sbs_gemFPP_hit_ENABLE_CM_V[11];   //[Ndata.sbs.gemFPP.hit.ENABLE_CM_V]
   Int_t           Ndata_sbs_gemFPP_hit_TSchi2_Umax;
   Double_t        sbs_gemFPP_hit_TSchi2_Umax[11];   //[Ndata.sbs.gemFPP.hit.TSchi2_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_TSchi2_Vmax;
   Double_t        sbs_gemFPP_hit_TSchi2_Vmax[11];   //[Ndata.sbs.gemFPP.hit.TSchi2_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_TSprob_Umax;
   Double_t        sbs_gemFPP_hit_TSprob_Umax[11];   //[Ndata.sbs.gemFPP.hit.TSprob_Umax]
   Int_t           Ndata_sbs_gemFPP_hit_TSprob_Vmax;
   Double_t        sbs_gemFPP_hit_TSprob_Vmax[11];   //[Ndata.sbs.gemFPP.hit.TSprob_Vmax]
   Int_t           Ndata_sbs_gemFPP_hit_Tavg;
   Double_t        sbs_gemFPP_hit_Tavg[11];   //[Ndata.sbs.gemFPP.hit.Tavg]
   Int_t           Ndata_sbs_gemFPP_hit_Tavg_corr;
   Double_t        sbs_gemFPP_hit_Tavg_corr[11];   //[Ndata.sbs.gemFPP.hit.Tavg_corr]
   Int_t           Ndata_sbs_gemFPP_hit_Tavg_deconv;
   Double_t        sbs_gemFPP_hit_Tavg_deconv[11];   //[Ndata.sbs.gemFPP.hit.Tavg_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_Tavg_fit;
   Double_t        sbs_gemFPP_hit_Tavg_fit[11];   //[Ndata.sbs.gemFPP.hit.Tavg_fit]
   Int_t           Ndata_sbs_gemFPP_hit_Ugain;
   Double_t        sbs_gemFPP_hit_Ugain[11];   //[Ndata.sbs.gemFPP.hit.Ugain]
   Int_t           Ndata_sbs_gemFPP_hit_Utime;
   Double_t        sbs_gemFPP_hit_Utime[11];   //[Ndata.sbs.gemFPP.hit.Utime]
   Int_t           Ndata_sbs_gemFPP_hit_UtimeDeconv;
   Double_t        sbs_gemFPP_hit_UtimeDeconv[11];   //[Ndata.sbs.gemFPP.hit.UtimeDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_UtimeFit;
   Double_t        sbs_gemFPP_hit_UtimeFit[11];   //[Ndata.sbs.gemFPP.hit.UtimeFit]
   Int_t           Ndata_sbs_gemFPP_hit_UtimeMaxStrip;
   Double_t        sbs_gemFPP_hit_UtimeMaxStrip[11];   //[Ndata.sbs.gemFPP.hit.UtimeMaxStrip]
   Int_t           Ndata_sbs_gemFPP_hit_UtimeMaxStripDeconv;
   Double_t        sbs_gemFPP_hit_UtimeMaxStripDeconv[11];   //[Ndata.sbs.gemFPP.hit.UtimeMaxStripDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_UtimeMaxStripFit;
   Double_t        sbs_gemFPP_hit_UtimeMaxStripFit[11];   //[Ndata.sbs.gemFPP.hit.UtimeMaxStripFit]
   Int_t           Ndata_sbs_gemFPP_hit_Vgain;
   Double_t        sbs_gemFPP_hit_Vgain[11];   //[Ndata.sbs.gemFPP.hit.Vgain]
   Int_t           Ndata_sbs_gemFPP_hit_Vtime;
   Double_t        sbs_gemFPP_hit_Vtime[11];   //[Ndata.sbs.gemFPP.hit.Vtime]
   Int_t           Ndata_sbs_gemFPP_hit_VtimeDeconv;
   Double_t        sbs_gemFPP_hit_VtimeDeconv[11];   //[Ndata.sbs.gemFPP.hit.VtimeDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_VtimeFit;
   Double_t        sbs_gemFPP_hit_VtimeFit[11];   //[Ndata.sbs.gemFPP.hit.VtimeFit]
   Int_t           Ndata_sbs_gemFPP_hit_VtimeMaxStrip;
   Double_t        sbs_gemFPP_hit_VtimeMaxStrip[11];   //[Ndata.sbs.gemFPP.hit.VtimeMaxStrip]
   Int_t           Ndata_sbs_gemFPP_hit_VtimeMaxStripDeconv;
   Double_t        sbs_gemFPP_hit_VtimeMaxStripDeconv[11];   //[Ndata.sbs.gemFPP.hit.VtimeMaxStripDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_VtimeMaxStripFit;
   Double_t        sbs_gemFPP_hit_VtimeMaxStripFit[11];   //[Ndata.sbs.gemFPP.hit.VtimeMaxStripFit]
   Int_t           Ndata_sbs_gemFPP_hit_adc_id_U;
   Double_t        sbs_gemFPP_hit_adc_id_U[11];   //[Ndata.sbs.gemFPP.hit.adc_id_U]
   Int_t           Ndata_sbs_gemFPP_hit_adc_id_V;
   Double_t        sbs_gemFPP_hit_adc_id_V[11];   //[Ndata.sbs.gemFPP.hit.adc_id_V]
   Int_t           Ndata_sbs_gemFPP_hit_ccor_clust;
   Double_t        sbs_gemFPP_hit_ccor_clust[11];   //[Ndata.sbs.gemFPP.hit.ccor_clust]
   Int_t           Ndata_sbs_gemFPP_hit_ccor_clust_deconv;
   Double_t        sbs_gemFPP_hit_ccor_clust_deconv[11];   //[Ndata.sbs.gemFPP.hit.ccor_clust_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_ccor_strip;
   Double_t        sbs_gemFPP_hit_ccor_strip[11];   //[Ndata.sbs.gemFPP.hit.ccor_strip]
   Int_t           Ndata_sbs_gemFPP_hit_ccor_strip_deconv;
   Double_t        sbs_gemFPP_hit_ccor_strip_deconv[11];   //[Ndata.sbs.gemFPP.hit.ccor_strip_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_crate_U;
   Double_t        sbs_gemFPP_hit_crate_U[11];   //[Ndata.sbs.gemFPP.hit.crate_U]
   Int_t           Ndata_sbs_gemFPP_hit_crate_V;
   Double_t        sbs_gemFPP_hit_crate_V[11];   //[Ndata.sbs.gemFPP.hit.crate_V]
   Int_t           Ndata_sbs_gemFPP_hit_deltat;
   Double_t        sbs_gemFPP_hit_deltat[11];   //[Ndata.sbs.gemFPP.hit.deltat]
   Int_t           Ndata_sbs_gemFPP_hit_deltat_deconv;
   Double_t        sbs_gemFPP_hit_deltat_deconv[11];   //[Ndata.sbs.gemFPP.hit.deltat_deconv]
   Int_t           Ndata_sbs_gemFPP_hit_deltat_fit;
   Double_t        sbs_gemFPP_hit_deltat_fit[11];   //[Ndata.sbs.gemFPP.hit.deltat_fit]
   Int_t           Ndata_sbs_gemFPP_hit_eresidu;
   Double_t        sbs_gemFPP_hit_eresidu[11];   //[Ndata.sbs.gemFPP.hit.eresidu]
   Int_t           Ndata_sbs_gemFPP_hit_eresidv;
   Double_t        sbs_gemFPP_hit_eresidv[11];   //[Ndata.sbs.gemFPP.hit.eresidv]
   Int_t           Ndata_sbs_gemFPP_hit_icombomaxUclustDeconv;
   Double_t        sbs_gemFPP_hit_icombomaxUclustDeconv[11];   //[Ndata.sbs.gemFPP.hit.icombomaxUclustDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_icombomaxUstripDeconv;
   Double_t        sbs_gemFPP_hit_icombomaxUstripDeconv[11];   //[Ndata.sbs.gemFPP.hit.icombomaxUstripDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_icombomaxVclustDeconv;
   Double_t        sbs_gemFPP_hit_icombomaxVclustDeconv[11];   //[Ndata.sbs.gemFPP.hit.icombomaxVclustDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_icombomaxVstripDeconv;
   Double_t        sbs_gemFPP_hit_icombomaxVstripDeconv[11];   //[Ndata.sbs.gemFPP.hit.icombomaxVstripDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_isampmaxUclust;
   Double_t        sbs_gemFPP_hit_isampmaxUclust[11];   //[Ndata.sbs.gemFPP.hit.isampmaxUclust]
   Int_t           Ndata_sbs_gemFPP_hit_isampmaxUclustDeconv;
   Double_t        sbs_gemFPP_hit_isampmaxUclustDeconv[11];   //[Ndata.sbs.gemFPP.hit.isampmaxUclustDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_isampmaxUstrip;
   Double_t        sbs_gemFPP_hit_isampmaxUstrip[11];   //[Ndata.sbs.gemFPP.hit.isampmaxUstrip]
   Int_t           Ndata_sbs_gemFPP_hit_isampmaxUstripDeconv;
   Double_t        sbs_gemFPP_hit_isampmaxUstripDeconv[11];   //[Ndata.sbs.gemFPP.hit.isampmaxUstripDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_isampmaxVclust;
   Double_t        sbs_gemFPP_hit_isampmaxVclust[11];   //[Ndata.sbs.gemFPP.hit.isampmaxVclust]
   Int_t           Ndata_sbs_gemFPP_hit_isampmaxVclustDeconv;
   Double_t        sbs_gemFPP_hit_isampmaxVclustDeconv[11];   //[Ndata.sbs.gemFPP.hit.isampmaxVclustDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_isampmaxVstrip;
   Double_t        sbs_gemFPP_hit_isampmaxVstrip[11];   //[Ndata.sbs.gemFPP.hit.isampmaxVstrip]
   Int_t           Ndata_sbs_gemFPP_hit_isampmaxVstripDeconv;
   Double_t        sbs_gemFPP_hit_isampmaxVstripDeconv[11];   //[Ndata.sbs.gemFPP.hit.isampmaxVstripDeconv]
   Int_t           Ndata_sbs_gemFPP_hit_layer;
   Double_t        sbs_gemFPP_hit_layer[11];   //[Ndata.sbs.gemFPP.hit.layer]
   Int_t           Ndata_sbs_gemFPP_hit_module;
   Double_t        sbs_gemFPP_hit_module[11];   //[Ndata.sbs.gemFPP.hit.module]
   Int_t           Ndata_sbs_gemFPP_hit_mpd_U;
   Double_t        sbs_gemFPP_hit_mpd_U[11];   //[Ndata.sbs.gemFPP.hit.mpd_U]
   Int_t           Ndata_sbs_gemFPP_hit_mpd_V;
   Double_t        sbs_gemFPP_hit_mpd_V[11];   //[Ndata.sbs.gemFPP.hit.mpd_V]
   Int_t           Ndata_sbs_gemFPP_hit_nstripu;
   Double_t        sbs_gemFPP_hit_nstripu[11];   //[Ndata.sbs.gemFPP.hit.nstripu]
   Int_t           Ndata_sbs_gemFPP_hit_nstripv;
   Double_t        sbs_gemFPP_hit_nstripv[11];   //[Ndata.sbs.gemFPP.hit.nstripv]
   Int_t           Ndata_sbs_gemFPP_hit_residu;
   Double_t        sbs_gemFPP_hit_residu[11];   //[Ndata.sbs.gemFPP.hit.residu]
   Int_t           Ndata_sbs_gemFPP_hit_residv;
   Double_t        sbs_gemFPP_hit_residv[11];   //[Ndata.sbs.gemFPP.hit.residv]
   Int_t           Ndata_sbs_gemFPP_hit_trackindex;
   Double_t        sbs_gemFPP_hit_trackindex[11];   //[Ndata.sbs.gemFPP.hit.trackindex]
   Int_t           Ndata_sbs_gemFPP_hit_u;
   Double_t        sbs_gemFPP_hit_u[11];   //[Ndata.sbs.gemFPP.hit.u]
   Int_t           Ndata_sbs_gemFPP_hit_umoment;
   Double_t        sbs_gemFPP_hit_umoment[11];   //[Ndata.sbs.gemFPP.hit.umoment]
   Int_t           Ndata_sbs_gemFPP_hit_usigma;
   Double_t        sbs_gemFPP_hit_usigma[11];   //[Ndata.sbs.gemFPP.hit.usigma]
   Int_t           Ndata_sbs_gemFPP_hit_ustriphi;
   Double_t        sbs_gemFPP_hit_ustriphi[11];   //[Ndata.sbs.gemFPP.hit.ustriphi]
   Int_t           Ndata_sbs_gemFPP_hit_ustriplo;
   Double_t        sbs_gemFPP_hit_ustriplo[11];   //[Ndata.sbs.gemFPP.hit.ustriplo]
   Int_t           Ndata_sbs_gemFPP_hit_ustripmax;
   Double_t        sbs_gemFPP_hit_ustripmax[11];   //[Ndata.sbs.gemFPP.hit.ustripmax]
   Int_t           Ndata_sbs_gemFPP_hit_v;
   Double_t        sbs_gemFPP_hit_v[11];   //[Ndata.sbs.gemFPP.hit.v]
   Int_t           Ndata_sbs_gemFPP_hit_vmoment;
   Double_t        sbs_gemFPP_hit_vmoment[11];   //[Ndata.sbs.gemFPP.hit.vmoment]
   Int_t           Ndata_sbs_gemFPP_hit_vsigma;
   Double_t        sbs_gemFPP_hit_vsigma[11];   //[Ndata.sbs.gemFPP.hit.vsigma]
   Int_t           Ndata_sbs_gemFPP_hit_vstriphi;
   Double_t        sbs_gemFPP_hit_vstriphi[11];   //[Ndata.sbs.gemFPP.hit.vstriphi]
   Int_t           Ndata_sbs_gemFPP_hit_vstriplo;
   Double_t        sbs_gemFPP_hit_vstriplo[11];   //[Ndata.sbs.gemFPP.hit.vstriplo]
   Int_t           Ndata_sbs_gemFPP_hit_vstripmax;
   Double_t        sbs_gemFPP_hit_vstripmax[11];   //[Ndata.sbs.gemFPP.hit.vstripmax]
   Int_t           Ndata_sbs_gemFPP_hit_xglobal;
   Double_t        sbs_gemFPP_hit_xglobal[11];   //[Ndata.sbs.gemFPP.hit.xglobal]
   Int_t           Ndata_sbs_gemFPP_hit_xlocal;
   Double_t        sbs_gemFPP_hit_xlocal[11];   //[Ndata.sbs.gemFPP.hit.xlocal]
   Int_t           Ndata_sbs_gemFPP_hit_yglobal;
   Double_t        sbs_gemFPP_hit_yglobal[11];   //[Ndata.sbs.gemFPP.hit.yglobal]
   Int_t           Ndata_sbs_gemFPP_hit_ylocal;
   Double_t        sbs_gemFPP_hit_ylocal[11];   //[Ndata.sbs.gemFPP.hit.ylocal]
   Int_t           Ndata_sbs_gemFPP_hit_zglobal;
   Double_t        sbs_gemFPP_hit_zglobal[11];   //[Ndata.sbs.gemFPP.hit.zglobal]
   Int_t           Ndata_sbs_gemFPP_n2Dhit_layer;
   Double_t        sbs_gemFPP_n2Dhit_layer[8];   //[Ndata.sbs.gemFPP.n2Dhit_layer]
   Int_t           Ndata_sbs_gemFPP_nclustu_layer;
   Double_t        sbs_gemFPP_nclustu_layer[8];   //[Ndata.sbs.gemFPP.nclustu_layer]
   Int_t           Ndata_sbs_gemFPP_nclustv_layer;
   Double_t        sbs_gemFPP_nclustv_layer[8];   //[Ndata.sbs.gemFPP.nclustv_layer]
   Int_t           Ndata_sbs_gemFPP_nstripsu_layer;
   Double_t        sbs_gemFPP_nstripsu_layer[8];   //[Ndata.sbs.gemFPP.nstripsu_layer]
   Int_t           Ndata_sbs_gemFPP_nstripsv_layer;
   Double_t        sbs_gemFPP_nstripsv_layer[8];   //[Ndata.sbs.gemFPP.nstripsv_layer]
   Int_t           Ndata_sbs_gemFPP_track_chi2ndf;
   Double_t        sbs_gemFPP_track_chi2ndf[2];   //[Ndata.sbs.gemFPP.track.chi2ndf]
   Int_t           Ndata_sbs_gemFPP_track_chi2ndf_hitquality;
   Double_t        sbs_gemFPP_track_chi2ndf_hitquality[2];   //[Ndata.sbs.gemFPP.track.chi2ndf_hitquality]
   Int_t           Ndata_sbs_gemFPP_track_ngoodhits;
   Double_t        sbs_gemFPP_track_ngoodhits[2];   //[Ndata.sbs.gemFPP.track.ngoodhits]
   Int_t           Ndata_sbs_gemFPP_track_nhits;
   Double_t        sbs_gemFPP_track_nhits[2];   //[Ndata.sbs.gemFPP.track.nhits]
   Int_t           Ndata_sbs_gemFPP_track_phi;
   Double_t        sbs_gemFPP_track_phi[2];   //[Ndata.sbs.gemFPP.track.phi]
   Int_t           Ndata_sbs_gemFPP_track_sclose;
   Double_t        sbs_gemFPP_track_sclose[2];   //[Ndata.sbs.gemFPP.track.sclose]
   Int_t           Ndata_sbs_gemFPP_track_t0;
   Double_t        sbs_gemFPP_track_t0[2];   //[Ndata.sbs.gemFPP.track.t0]
   Int_t           Ndata_sbs_gemFPP_track_theta;
   Double_t        sbs_gemFPP_track_theta[2];   //[Ndata.sbs.gemFPP.track.theta]
   Int_t           Ndata_sbs_gemFPP_track_x;
   Double_t        sbs_gemFPP_track_x[2];   //[Ndata.sbs.gemFPP.track.x]
   Int_t           Ndata_sbs_gemFPP_track_xp;
   Double_t        sbs_gemFPP_track_xp[2];   //[Ndata.sbs.gemFPP.track.xp]
   Int_t           Ndata_sbs_gemFPP_track_y;
   Double_t        sbs_gemFPP_track_y[2];   //[Ndata.sbs.gemFPP.track.y]
   Int_t           Ndata_sbs_gemFPP_track_yp;
   Double_t        sbs_gemFPP_track_yp[2];   //[Ndata.sbs.gemFPP.track.yp]
   Int_t           Ndata_sbs_gemFPP_track_zclose;
   Double_t        sbs_gemFPP_track_zclose[2];   //[Ndata.sbs.gemFPP.track.zclose]
   Int_t           Ndata_sbs_gemFT_hit_ADCU;
   Double_t        sbs_gemFT_hit_ADCU[8];   //[Ndata.sbs.gemFT.hit.ADCU]
   Int_t           Ndata_sbs_gemFT_hit_ADCU_deconv;
   Double_t        sbs_gemFT_hit_ADCU_deconv[8];   //[Ndata.sbs.gemFT.hit.ADCU_deconv]
   Int_t           Ndata_sbs_gemFT_hit_ADCV;
   Double_t        sbs_gemFT_hit_ADCV[8];   //[Ndata.sbs.gemFT.hit.ADCV]
   Int_t           Ndata_sbs_gemFT_hit_ADCV_deconv;
   Double_t        sbs_gemFT_hit_ADCV_deconv[8];   //[Ndata.sbs.gemFT.hit.ADCV_deconv]
   Int_t           Ndata_sbs_gemFT_hit_ADCasym;
   Double_t        sbs_gemFT_hit_ADCasym[8];   //[Ndata.sbs.gemFT.hit.ADCasym]
   Int_t           Ndata_sbs_gemFT_hit_ADCasym_deconv;
   Double_t        sbs_gemFT_hit_ADCasym_deconv[8];   //[Ndata.sbs.gemFT.hit.ADCasym_deconv]
   Int_t           Ndata_sbs_gemFT_hit_ADCavg;
   Double_t        sbs_gemFT_hit_ADCavg[8];   //[Ndata.sbs.gemFT.hit.ADCavg]
   Int_t           Ndata_sbs_gemFT_hit_ADCavg_deconv;
   Double_t        sbs_gemFT_hit_ADCavg_deconv[8];   //[Ndata.sbs.gemFT.hit.ADCavg_deconv]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac0_Umax;
   Double_t        sbs_gemFT_hit_ADCfrac0_Umax[8];   //[Ndata.sbs.gemFT.hit.ADCfrac0_Umax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac0_Vmax;
   Double_t        sbs_gemFT_hit_ADCfrac0_Vmax[8];   //[Ndata.sbs.gemFT.hit.ADCfrac0_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac1_Umax;
   Double_t        sbs_gemFT_hit_ADCfrac1_Umax[8];   //[Ndata.sbs.gemFT.hit.ADCfrac1_Umax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac1_Vmax;
   Double_t        sbs_gemFT_hit_ADCfrac1_Vmax[8];   //[Ndata.sbs.gemFT.hit.ADCfrac1_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac2_Umax;
   Double_t        sbs_gemFT_hit_ADCfrac2_Umax[8];   //[Ndata.sbs.gemFT.hit.ADCfrac2_Umax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac2_Vmax;
   Double_t        sbs_gemFT_hit_ADCfrac2_Vmax[8];   //[Ndata.sbs.gemFT.hit.ADCfrac2_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac3_Umax;
   Double_t        sbs_gemFT_hit_ADCfrac3_Umax[8];   //[Ndata.sbs.gemFT.hit.ADCfrac3_Umax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac3_Vmax;
   Double_t        sbs_gemFT_hit_ADCfrac3_Vmax[8];   //[Ndata.sbs.gemFT.hit.ADCfrac3_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac4_Umax;
   Double_t        sbs_gemFT_hit_ADCfrac4_Umax[8];   //[Ndata.sbs.gemFT.hit.ADCfrac4_Umax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac4_Vmax;
   Double_t        sbs_gemFT_hit_ADCfrac4_Vmax[8];   //[Ndata.sbs.gemFT.hit.ADCfrac4_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac5_Umax;
   Double_t        sbs_gemFT_hit_ADCfrac5_Umax[8];   //[Ndata.sbs.gemFT.hit.ADCfrac5_Umax]
   Int_t           Ndata_sbs_gemFT_hit_ADCfrac5_Vmax;
   Double_t        sbs_gemFT_hit_ADCfrac5_Vmax[8];   //[Ndata.sbs.gemFT.hit.ADCfrac5_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxcomboUclust_deconv;
   Double_t        sbs_gemFT_hit_ADCmaxcomboUclust_deconv[8];   //[Ndata.sbs.gemFT.hit.ADCmaxcomboUclust_deconv]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxcomboVclust_deconv;
   Double_t        sbs_gemFT_hit_ADCmaxcomboVclust_deconv[8];   //[Ndata.sbs.gemFT.hit.ADCmaxcomboVclust_deconv]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxsampU;
   Double_t        sbs_gemFT_hit_ADCmaxsampU[8];   //[Ndata.sbs.gemFT.hit.ADCmaxsampU]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxsampUclust;
   Double_t        sbs_gemFT_hit_ADCmaxsampUclust[8];   //[Ndata.sbs.gemFT.hit.ADCmaxsampUclust]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxsampUclust_deconv;
   Double_t        sbs_gemFT_hit_ADCmaxsampUclust_deconv[8];   //[Ndata.sbs.gemFT.hit.ADCmaxsampUclust_deconv]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxsampV;
   Double_t        sbs_gemFT_hit_ADCmaxsampV[8];   //[Ndata.sbs.gemFT.hit.ADCmaxsampV]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxsampVclust;
   Double_t        sbs_gemFT_hit_ADCmaxsampVclust[8];   //[Ndata.sbs.gemFT.hit.ADCmaxsampVclust]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxsampVclust_deconv;
   Double_t        sbs_gemFT_hit_ADCmaxsampVclust_deconv[8];   //[Ndata.sbs.gemFT.hit.ADCmaxsampVclust_deconv]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxstripU;
   Double_t        sbs_gemFT_hit_ADCmaxstripU[8];   //[Ndata.sbs.gemFT.hit.ADCmaxstripU]
   Int_t           Ndata_sbs_gemFT_hit_ADCmaxstripV;
   Double_t        sbs_gemFT_hit_ADCmaxstripV[8];   //[Ndata.sbs.gemFT.hit.ADCmaxstripV]
   Int_t           Ndata_sbs_gemFT_hit_BUILD_ALL_SAMPLES_U;
   Double_t        sbs_gemFT_hit_BUILD_ALL_SAMPLES_U[8];   //[Ndata.sbs.gemFT.hit.BUILD_ALL_SAMPLES_U]
   Int_t           Ndata_sbs_gemFT_hit_BUILD_ALL_SAMPLES_V;
   Double_t        sbs_gemFT_hit_BUILD_ALL_SAMPLES_V[8];   //[Ndata.sbs.gemFT.hit.BUILD_ALL_SAMPLES_V]
   Int_t           Ndata_sbs_gemFT_hit_CM_GOOD_U;
   Double_t        sbs_gemFT_hit_CM_GOOD_U[8];   //[Ndata.sbs.gemFT.hit.CM_GOOD_U]
   Int_t           Ndata_sbs_gemFT_hit_CM_GOOD_V;
   Double_t        sbs_gemFT_hit_CM_GOOD_V[8];   //[Ndata.sbs.gemFT.hit.CM_GOOD_V]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC0_Umax;
   Double_t        sbs_gemFT_hit_DeconvADC0_Umax[8];   //[Ndata.sbs.gemFT.hit.DeconvADC0_Umax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC0_Vmax;
   Double_t        sbs_gemFT_hit_DeconvADC0_Vmax[8];   //[Ndata.sbs.gemFT.hit.DeconvADC0_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC1_Umax;
   Double_t        sbs_gemFT_hit_DeconvADC1_Umax[8];   //[Ndata.sbs.gemFT.hit.DeconvADC1_Umax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC1_Vmax;
   Double_t        sbs_gemFT_hit_DeconvADC1_Vmax[8];   //[Ndata.sbs.gemFT.hit.DeconvADC1_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC2_Umax;
   Double_t        sbs_gemFT_hit_DeconvADC2_Umax[8];   //[Ndata.sbs.gemFT.hit.DeconvADC2_Umax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC2_Vmax;
   Double_t        sbs_gemFT_hit_DeconvADC2_Vmax[8];   //[Ndata.sbs.gemFT.hit.DeconvADC2_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC3_Umax;
   Double_t        sbs_gemFT_hit_DeconvADC3_Umax[8];   //[Ndata.sbs.gemFT.hit.DeconvADC3_Umax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC3_Vmax;
   Double_t        sbs_gemFT_hit_DeconvADC3_Vmax[8];   //[Ndata.sbs.gemFT.hit.DeconvADC3_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC4_Umax;
   Double_t        sbs_gemFT_hit_DeconvADC4_Umax[8];   //[Ndata.sbs.gemFT.hit.DeconvADC4_Umax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC4_Vmax;
   Double_t        sbs_gemFT_hit_DeconvADC4_Vmax[8];   //[Ndata.sbs.gemFT.hit.DeconvADC4_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC5_Umax;
   Double_t        sbs_gemFT_hit_DeconvADC5_Umax[8];   //[Ndata.sbs.gemFT.hit.DeconvADC5_Umax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADC5_Vmax;
   Double_t        sbs_gemFT_hit_DeconvADC5_Vmax[8];   //[Ndata.sbs.gemFT.hit.DeconvADC5_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADCmaxcomboU;
   Double_t        sbs_gemFT_hit_DeconvADCmaxcomboU[8];   //[Ndata.sbs.gemFT.hit.DeconvADCmaxcomboU]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADCmaxcomboV;
   Double_t        sbs_gemFT_hit_DeconvADCmaxcomboV[8];   //[Ndata.sbs.gemFT.hit.DeconvADCmaxcomboV]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADCmaxsampU;
   Double_t        sbs_gemFT_hit_DeconvADCmaxsampU[8];   //[Ndata.sbs.gemFT.hit.DeconvADCmaxsampU]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADCmaxsampV;
   Double_t        sbs_gemFT_hit_DeconvADCmaxsampV[8];   //[Ndata.sbs.gemFT.hit.DeconvADCmaxsampV]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADCmaxstripU;
   Double_t        sbs_gemFT_hit_DeconvADCmaxstripU[8];   //[Ndata.sbs.gemFT.hit.DeconvADCmaxstripU]
   Int_t           Ndata_sbs_gemFT_hit_DeconvADCmaxstripV;
   Double_t        sbs_gemFT_hit_DeconvADCmaxstripV[8];   //[Ndata.sbs.gemFT.hit.DeconvADCmaxstripV]
   Int_t           Ndata_sbs_gemFT_hit_ENABLE_CM_U;
   Double_t        sbs_gemFT_hit_ENABLE_CM_U[8];   //[Ndata.sbs.gemFT.hit.ENABLE_CM_U]
   Int_t           Ndata_sbs_gemFT_hit_ENABLE_CM_V;
   Double_t        sbs_gemFT_hit_ENABLE_CM_V[8];   //[Ndata.sbs.gemFT.hit.ENABLE_CM_V]
   Int_t           Ndata_sbs_gemFT_hit_TSchi2_Umax;
   Double_t        sbs_gemFT_hit_TSchi2_Umax[8];   //[Ndata.sbs.gemFT.hit.TSchi2_Umax]
   Int_t           Ndata_sbs_gemFT_hit_TSchi2_Vmax;
   Double_t        sbs_gemFT_hit_TSchi2_Vmax[8];   //[Ndata.sbs.gemFT.hit.TSchi2_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_TSprob_Umax;
   Double_t        sbs_gemFT_hit_TSprob_Umax[8];   //[Ndata.sbs.gemFT.hit.TSprob_Umax]
   Int_t           Ndata_sbs_gemFT_hit_TSprob_Vmax;
   Double_t        sbs_gemFT_hit_TSprob_Vmax[8];   //[Ndata.sbs.gemFT.hit.TSprob_Vmax]
   Int_t           Ndata_sbs_gemFT_hit_Tavg;
   Double_t        sbs_gemFT_hit_Tavg[8];   //[Ndata.sbs.gemFT.hit.Tavg]
   Int_t           Ndata_sbs_gemFT_hit_Tavg_corr;
   Double_t        sbs_gemFT_hit_Tavg_corr[8];   //[Ndata.sbs.gemFT.hit.Tavg_corr]
   Int_t           Ndata_sbs_gemFT_hit_Tavg_deconv;
   Double_t        sbs_gemFT_hit_Tavg_deconv[8];   //[Ndata.sbs.gemFT.hit.Tavg_deconv]
   Int_t           Ndata_sbs_gemFT_hit_Tavg_fit;
   Double_t        sbs_gemFT_hit_Tavg_fit[8];   //[Ndata.sbs.gemFT.hit.Tavg_fit]
   Int_t           Ndata_sbs_gemFT_hit_Ugain;
   Double_t        sbs_gemFT_hit_Ugain[8];   //[Ndata.sbs.gemFT.hit.Ugain]
   Int_t           Ndata_sbs_gemFT_hit_Utime;
   Double_t        sbs_gemFT_hit_Utime[8];   //[Ndata.sbs.gemFT.hit.Utime]
   Int_t           Ndata_sbs_gemFT_hit_UtimeDeconv;
   Double_t        sbs_gemFT_hit_UtimeDeconv[8];   //[Ndata.sbs.gemFT.hit.UtimeDeconv]
   Int_t           Ndata_sbs_gemFT_hit_UtimeFit;
   Double_t        sbs_gemFT_hit_UtimeFit[8];   //[Ndata.sbs.gemFT.hit.UtimeFit]
   Int_t           Ndata_sbs_gemFT_hit_UtimeMaxStrip;
   Double_t        sbs_gemFT_hit_UtimeMaxStrip[8];   //[Ndata.sbs.gemFT.hit.UtimeMaxStrip]
   Int_t           Ndata_sbs_gemFT_hit_UtimeMaxStripDeconv;
   Double_t        sbs_gemFT_hit_UtimeMaxStripDeconv[8];   //[Ndata.sbs.gemFT.hit.UtimeMaxStripDeconv]
   Int_t           Ndata_sbs_gemFT_hit_UtimeMaxStripFit;
   Double_t        sbs_gemFT_hit_UtimeMaxStripFit[8];   //[Ndata.sbs.gemFT.hit.UtimeMaxStripFit]
   Int_t           Ndata_sbs_gemFT_hit_Vgain;
   Double_t        sbs_gemFT_hit_Vgain[8];   //[Ndata.sbs.gemFT.hit.Vgain]
   Int_t           Ndata_sbs_gemFT_hit_Vtime;
   Double_t        sbs_gemFT_hit_Vtime[8];   //[Ndata.sbs.gemFT.hit.Vtime]
   Int_t           Ndata_sbs_gemFT_hit_VtimeDeconv;
   Double_t        sbs_gemFT_hit_VtimeDeconv[8];   //[Ndata.sbs.gemFT.hit.VtimeDeconv]
   Int_t           Ndata_sbs_gemFT_hit_VtimeFit;
   Double_t        sbs_gemFT_hit_VtimeFit[8];   //[Ndata.sbs.gemFT.hit.VtimeFit]
   Int_t           Ndata_sbs_gemFT_hit_VtimeMaxStrip;
   Double_t        sbs_gemFT_hit_VtimeMaxStrip[8];   //[Ndata.sbs.gemFT.hit.VtimeMaxStrip]
   Int_t           Ndata_sbs_gemFT_hit_VtimeMaxStripDeconv;
   Double_t        sbs_gemFT_hit_VtimeMaxStripDeconv[8];   //[Ndata.sbs.gemFT.hit.VtimeMaxStripDeconv]
   Int_t           Ndata_sbs_gemFT_hit_VtimeMaxStripFit;
   Double_t        sbs_gemFT_hit_VtimeMaxStripFit[8];   //[Ndata.sbs.gemFT.hit.VtimeMaxStripFit]
   Int_t           Ndata_sbs_gemFT_hit_adc_id_U;
   Double_t        sbs_gemFT_hit_adc_id_U[8];   //[Ndata.sbs.gemFT.hit.adc_id_U]
   Int_t           Ndata_sbs_gemFT_hit_adc_id_V;
   Double_t        sbs_gemFT_hit_adc_id_V[8];   //[Ndata.sbs.gemFT.hit.adc_id_V]
   Int_t           Ndata_sbs_gemFT_hit_ccor_clust;
   Double_t        sbs_gemFT_hit_ccor_clust[8];   //[Ndata.sbs.gemFT.hit.ccor_clust]
   Int_t           Ndata_sbs_gemFT_hit_ccor_clust_deconv;
   Double_t        sbs_gemFT_hit_ccor_clust_deconv[8];   //[Ndata.sbs.gemFT.hit.ccor_clust_deconv]
   Int_t           Ndata_sbs_gemFT_hit_ccor_strip;
   Double_t        sbs_gemFT_hit_ccor_strip[8];   //[Ndata.sbs.gemFT.hit.ccor_strip]
   Int_t           Ndata_sbs_gemFT_hit_ccor_strip_deconv;
   Double_t        sbs_gemFT_hit_ccor_strip_deconv[8];   //[Ndata.sbs.gemFT.hit.ccor_strip_deconv]
   Int_t           Ndata_sbs_gemFT_hit_crate_U;
   Double_t        sbs_gemFT_hit_crate_U[8];   //[Ndata.sbs.gemFT.hit.crate_U]
   Int_t           Ndata_sbs_gemFT_hit_crate_V;
   Double_t        sbs_gemFT_hit_crate_V[8];   //[Ndata.sbs.gemFT.hit.crate_V]
   Int_t           Ndata_sbs_gemFT_hit_deltat;
   Double_t        sbs_gemFT_hit_deltat[8];   //[Ndata.sbs.gemFT.hit.deltat]
   Int_t           Ndata_sbs_gemFT_hit_deltat_deconv;
   Double_t        sbs_gemFT_hit_deltat_deconv[8];   //[Ndata.sbs.gemFT.hit.deltat_deconv]
   Int_t           Ndata_sbs_gemFT_hit_deltat_fit;
   Double_t        sbs_gemFT_hit_deltat_fit[8];   //[Ndata.sbs.gemFT.hit.deltat_fit]
   Int_t           Ndata_sbs_gemFT_hit_eresidu;
   Double_t        sbs_gemFT_hit_eresidu[8];   //[Ndata.sbs.gemFT.hit.eresidu]
   Int_t           Ndata_sbs_gemFT_hit_eresidv;
   Double_t        sbs_gemFT_hit_eresidv[8];   //[Ndata.sbs.gemFT.hit.eresidv]
   Int_t           Ndata_sbs_gemFT_hit_icombomaxUclustDeconv;
   Double_t        sbs_gemFT_hit_icombomaxUclustDeconv[8];   //[Ndata.sbs.gemFT.hit.icombomaxUclustDeconv]
   Int_t           Ndata_sbs_gemFT_hit_icombomaxUstripDeconv;
   Double_t        sbs_gemFT_hit_icombomaxUstripDeconv[8];   //[Ndata.sbs.gemFT.hit.icombomaxUstripDeconv]
   Int_t           Ndata_sbs_gemFT_hit_icombomaxVclustDeconv;
   Double_t        sbs_gemFT_hit_icombomaxVclustDeconv[8];   //[Ndata.sbs.gemFT.hit.icombomaxVclustDeconv]
   Int_t           Ndata_sbs_gemFT_hit_icombomaxVstripDeconv;
   Double_t        sbs_gemFT_hit_icombomaxVstripDeconv[8];   //[Ndata.sbs.gemFT.hit.icombomaxVstripDeconv]
   Int_t           Ndata_sbs_gemFT_hit_isampmaxUclust;
   Double_t        sbs_gemFT_hit_isampmaxUclust[8];   //[Ndata.sbs.gemFT.hit.isampmaxUclust]
   Int_t           Ndata_sbs_gemFT_hit_isampmaxUclustDeconv;
   Double_t        sbs_gemFT_hit_isampmaxUclustDeconv[8];   //[Ndata.sbs.gemFT.hit.isampmaxUclustDeconv]
   Int_t           Ndata_sbs_gemFT_hit_isampmaxUstrip;
   Double_t        sbs_gemFT_hit_isampmaxUstrip[8];   //[Ndata.sbs.gemFT.hit.isampmaxUstrip]
   Int_t           Ndata_sbs_gemFT_hit_isampmaxUstripDeconv;
   Double_t        sbs_gemFT_hit_isampmaxUstripDeconv[8];   //[Ndata.sbs.gemFT.hit.isampmaxUstripDeconv]
   Int_t           Ndata_sbs_gemFT_hit_isampmaxVclust;
   Double_t        sbs_gemFT_hit_isampmaxVclust[8];   //[Ndata.sbs.gemFT.hit.isampmaxVclust]
   Int_t           Ndata_sbs_gemFT_hit_isampmaxVclustDeconv;
   Double_t        sbs_gemFT_hit_isampmaxVclustDeconv[8];   //[Ndata.sbs.gemFT.hit.isampmaxVclustDeconv]
   Int_t           Ndata_sbs_gemFT_hit_isampmaxVstrip;
   Double_t        sbs_gemFT_hit_isampmaxVstrip[8];   //[Ndata.sbs.gemFT.hit.isampmaxVstrip]
   Int_t           Ndata_sbs_gemFT_hit_isampmaxVstripDeconv;
   Double_t        sbs_gemFT_hit_isampmaxVstripDeconv[8];   //[Ndata.sbs.gemFT.hit.isampmaxVstripDeconv]
   Int_t           Ndata_sbs_gemFT_hit_layer;
   Double_t        sbs_gemFT_hit_layer[8];   //[Ndata.sbs.gemFT.hit.layer]
   Int_t           Ndata_sbs_gemFT_hit_module;
   Double_t        sbs_gemFT_hit_module[8];   //[Ndata.sbs.gemFT.hit.module]
   Int_t           Ndata_sbs_gemFT_hit_mpd_U;
   Double_t        sbs_gemFT_hit_mpd_U[8];   //[Ndata.sbs.gemFT.hit.mpd_U]
   Int_t           Ndata_sbs_gemFT_hit_mpd_V;
   Double_t        sbs_gemFT_hit_mpd_V[8];   //[Ndata.sbs.gemFT.hit.mpd_V]
   Int_t           Ndata_sbs_gemFT_hit_nstripu;
   Double_t        sbs_gemFT_hit_nstripu[8];   //[Ndata.sbs.gemFT.hit.nstripu]
   Int_t           Ndata_sbs_gemFT_hit_nstripv;
   Double_t        sbs_gemFT_hit_nstripv[8];   //[Ndata.sbs.gemFT.hit.nstripv]
   Int_t           Ndata_sbs_gemFT_hit_residu;
   Double_t        sbs_gemFT_hit_residu[8];   //[Ndata.sbs.gemFT.hit.residu]
   Int_t           Ndata_sbs_gemFT_hit_residv;
   Double_t        sbs_gemFT_hit_residv[8];   //[Ndata.sbs.gemFT.hit.residv]
   Int_t           Ndata_sbs_gemFT_hit_trackindex;
   Double_t        sbs_gemFT_hit_trackindex[8];   //[Ndata.sbs.gemFT.hit.trackindex]
   Int_t           Ndata_sbs_gemFT_hit_u;
   Double_t        sbs_gemFT_hit_u[8];   //[Ndata.sbs.gemFT.hit.u]
   Int_t           Ndata_sbs_gemFT_hit_umoment;
   Double_t        sbs_gemFT_hit_umoment[8];   //[Ndata.sbs.gemFT.hit.umoment]
   Int_t           Ndata_sbs_gemFT_hit_usigma;
   Double_t        sbs_gemFT_hit_usigma[8];   //[Ndata.sbs.gemFT.hit.usigma]
   Int_t           Ndata_sbs_gemFT_hit_ustriphi;
   Double_t        sbs_gemFT_hit_ustriphi[8];   //[Ndata.sbs.gemFT.hit.ustriphi]
   Int_t           Ndata_sbs_gemFT_hit_ustriplo;
   Double_t        sbs_gemFT_hit_ustriplo[8];   //[Ndata.sbs.gemFT.hit.ustriplo]
   Int_t           Ndata_sbs_gemFT_hit_ustripmax;
   Double_t        sbs_gemFT_hit_ustripmax[8];   //[Ndata.sbs.gemFT.hit.ustripmax]
   Int_t           Ndata_sbs_gemFT_hit_v;
   Double_t        sbs_gemFT_hit_v[8];   //[Ndata.sbs.gemFT.hit.v]
   Int_t           Ndata_sbs_gemFT_hit_vmoment;
   Double_t        sbs_gemFT_hit_vmoment[8];   //[Ndata.sbs.gemFT.hit.vmoment]
   Int_t           Ndata_sbs_gemFT_hit_vsigma;
   Double_t        sbs_gemFT_hit_vsigma[8];   //[Ndata.sbs.gemFT.hit.vsigma]
   Int_t           Ndata_sbs_gemFT_hit_vstriphi;
   Double_t        sbs_gemFT_hit_vstriphi[8];   //[Ndata.sbs.gemFT.hit.vstriphi]
   Int_t           Ndata_sbs_gemFT_hit_vstriplo;
   Double_t        sbs_gemFT_hit_vstriplo[8];   //[Ndata.sbs.gemFT.hit.vstriplo]
   Int_t           Ndata_sbs_gemFT_hit_vstripmax;
   Double_t        sbs_gemFT_hit_vstripmax[8];   //[Ndata.sbs.gemFT.hit.vstripmax]
   Int_t           Ndata_sbs_gemFT_hit_xglobal;
   Double_t        sbs_gemFT_hit_xglobal[8];   //[Ndata.sbs.gemFT.hit.xglobal]
   Int_t           Ndata_sbs_gemFT_hit_xlocal;
   Double_t        sbs_gemFT_hit_xlocal[8];   //[Ndata.sbs.gemFT.hit.xlocal]
   Int_t           Ndata_sbs_gemFT_hit_yglobal;
   Double_t        sbs_gemFT_hit_yglobal[8];   //[Ndata.sbs.gemFT.hit.yglobal]
   Int_t           Ndata_sbs_gemFT_hit_ylocal;
   Double_t        sbs_gemFT_hit_ylocal[8];   //[Ndata.sbs.gemFT.hit.ylocal]
   Int_t           Ndata_sbs_gemFT_hit_zglobal;
   Double_t        sbs_gemFT_hit_zglobal[8];   //[Ndata.sbs.gemFT.hit.zglobal]
   Int_t           Ndata_sbs_gemFT_n2Dhit_layer;
   Double_t        sbs_gemFT_n2Dhit_layer[8];   //[Ndata.sbs.gemFT.n2Dhit_layer]
   Int_t           Ndata_sbs_gemFT_nclustu_layer;
   Double_t        sbs_gemFT_nclustu_layer[8];   //[Ndata.sbs.gemFT.nclustu_layer]
   Int_t           Ndata_sbs_gemFT_nclustv_layer;
   Double_t        sbs_gemFT_nclustv_layer[8];   //[Ndata.sbs.gemFT.nclustv_layer]
   Int_t           Ndata_sbs_gemFT_nstripsu_layer;
   Double_t        sbs_gemFT_nstripsu_layer[8];   //[Ndata.sbs.gemFT.nstripsu_layer]
   Int_t           Ndata_sbs_gemFT_nstripsv_layer;
   Double_t        sbs_gemFT_nstripsv_layer[8];   //[Ndata.sbs.gemFT.nstripsv_layer]
   Int_t           Ndata_sbs_gemFT_track_chi2ndf;
   Double_t        sbs_gemFT_track_chi2ndf[1];   //[Ndata.sbs.gemFT.track.chi2ndf]
   Int_t           Ndata_sbs_gemFT_track_chi2ndf_hitquality;
   Double_t        sbs_gemFT_track_chi2ndf_hitquality[1];   //[Ndata.sbs.gemFT.track.chi2ndf_hitquality]
   Int_t           Ndata_sbs_gemFT_track_ngoodhits;
   Double_t        sbs_gemFT_track_ngoodhits[1];   //[Ndata.sbs.gemFT.track.ngoodhits]
   Int_t           Ndata_sbs_gemFT_track_nhits;
   Double_t        sbs_gemFT_track_nhits[1];   //[Ndata.sbs.gemFT.track.nhits]
   Int_t           Ndata_sbs_gemFT_track_t0;
   Double_t        sbs_gemFT_track_t0[1];   //[Ndata.sbs.gemFT.track.t0]
   Int_t           Ndata_sbs_gemFT_track_x;
   Double_t        sbs_gemFT_track_x[1];   //[Ndata.sbs.gemFT.track.x]
   Int_t           Ndata_sbs_gemFT_track_xp;
   Double_t        sbs_gemFT_track_xp[1];   //[Ndata.sbs.gemFT.track.xp]
   Int_t           Ndata_sbs_gemFT_track_y;
   Double_t        sbs_gemFT_track_y[1];   //[Ndata.sbs.gemFT.track.y]
   Int_t           Ndata_sbs_gemFT_track_yp;
   Double_t        sbs_gemFT_track_yp[1];   //[Ndata.sbs.gemFT.track.yp]
   Int_t           Ndata_sbs_hcal_clus_adctime;
   Double_t        sbs_hcal_clus_adctime[6];   //[Ndata.sbs.hcal.clus.adctime]
   Int_t           Ndata_sbs_hcal_clus_again;
   Double_t        sbs_hcal_clus_again[6];   //[Ndata.sbs.hcal.clus.again]
   Int_t           Ndata_sbs_hcal_clus_atimeblk;
   Double_t        sbs_hcal_clus_atimeblk[6];   //[Ndata.sbs.hcal.clus.atimeblk]
   Int_t           Ndata_sbs_hcal_clus_col;
   Double_t        sbs_hcal_clus_col[6];   //[Ndata.sbs.hcal.clus.col]
   Int_t           Ndata_sbs_hcal_clus_col_goodtdc;
   Double_t        sbs_hcal_clus_col_goodtdc[6];   //[Ndata.sbs.hcal.clus.col_goodtdc]
   Int_t           Ndata_sbs_hcal_clus_e;
   Double_t        sbs_hcal_clus_e[6];   //[Ndata.sbs.hcal.clus.e]
   Int_t           Ndata_sbs_hcal_clus_e_goodtdc;
   Double_t        sbs_hcal_clus_e_goodtdc[6];   //[Ndata.sbs.hcal.clus.e_goodtdc]
   Int_t           Ndata_sbs_hcal_clus_eblk;
   Double_t        sbs_hcal_clus_eblk[6];   //[Ndata.sbs.hcal.clus.eblk]
   Int_t           Ndata_sbs_hcal_clus_eblk_goodtdc;
   Double_t        sbs_hcal_clus_eblk_goodtdc[6];   //[Ndata.sbs.hcal.clus.eblk_goodtdc]
   Int_t           Ndata_sbs_hcal_clus_id;
   Double_t        sbs_hcal_clus_id[6];   //[Ndata.sbs.hcal.clus.id]
   Int_t           Ndata_sbs_hcal_clus_id_goodtdc;
   Double_t        sbs_hcal_clus_id_goodtdc[6];   //[Ndata.sbs.hcal.clus.id_goodtdc]
   Int_t           Ndata_sbs_hcal_clus_nblk;
   Double_t        sbs_hcal_clus_nblk[6];   //[Ndata.sbs.hcal.clus.nblk]
   Int_t           Ndata_sbs_hcal_clus_nblk_goodtdc;
   Double_t        sbs_hcal_clus_nblk_goodtdc[6];   //[Ndata.sbs.hcal.clus.nblk_goodtdc]
   Int_t           Ndata_sbs_hcal_clus_row;
   Double_t        sbs_hcal_clus_row[6];   //[Ndata.sbs.hcal.clus.row]
   Int_t           Ndata_sbs_hcal_clus_row_goodtdc;
   Double_t        sbs_hcal_clus_row_goodtdc[6];   //[Ndata.sbs.hcal.clus.row_goodtdc]
   Int_t           Ndata_sbs_hcal_clus_tdctime;
   Double_t        sbs_hcal_clus_tdctime[6];   //[Ndata.sbs.hcal.clus.tdctime]
   Int_t           Ndata_sbs_hcal_clus_tdctime_tw;
   Double_t        sbs_hcal_clus_tdctime_tw[6];   //[Ndata.sbs.hcal.clus.tdctime_tw]
   Int_t           Ndata_sbs_hcal_clus_tdctimeblk;
   Double_t        sbs_hcal_clus_tdctimeblk[6];   //[Ndata.sbs.hcal.clus.tdctimeblk]
   Int_t           Ndata_sbs_hcal_clus_tdctimeblk_tw;
   Double_t        sbs_hcal_clus_tdctimeblk_tw[6];   //[Ndata.sbs.hcal.clus.tdctimeblk_tw]
   Int_t           Ndata_sbs_hcal_clus_x;
   Double_t        sbs_hcal_clus_x[6];   //[Ndata.sbs.hcal.clus.x]
   Int_t           Ndata_sbs_hcal_clus_y;
   Double_t        sbs_hcal_clus_y[6];   //[Ndata.sbs.hcal.clus.y]
   Int_t           Ndata_sbs_hcal_clus_blk_again;
   Double_t        sbs_hcal_clus_blk_again[15];   //[Ndata.sbs.hcal.clus_blk.again]
   Int_t           Ndata_sbs_hcal_clus_blk_atime;
   Double_t        sbs_hcal_clus_blk_atime[15];   //[Ndata.sbs.hcal.clus_blk.atime]
   Int_t           Ndata_sbs_hcal_clus_blk_col;
   Double_t        sbs_hcal_clus_blk_col[15];   //[Ndata.sbs.hcal.clus_blk.col]
   Int_t           Ndata_sbs_hcal_clus_blk_e;
   Double_t        sbs_hcal_clus_blk_e[15];   //[Ndata.sbs.hcal.clus_blk.e]
   Int_t           Ndata_sbs_hcal_clus_blk_id;
   Double_t        sbs_hcal_clus_blk_id[15];   //[Ndata.sbs.hcal.clus_blk.id]
   Int_t           Ndata_sbs_hcal_clus_blk_row;
   Double_t        sbs_hcal_clus_blk_row[15];   //[Ndata.sbs.hcal.clus_blk.row]
   Int_t           Ndata_sbs_hcal_clus_blk_tdctime;
   Double_t        sbs_hcal_clus_blk_tdctime[15];   //[Ndata.sbs.hcal.clus_blk.tdctime]
   Int_t           Ndata_sbs_hcal_clus_blk_tdctime_tw;
   Double_t        sbs_hcal_clus_blk_tdctime_tw[15];   //[Ndata.sbs.hcal.clus_blk.tdctime_tw]
   Int_t           Ndata_sbs_hcal_clus_blk_x;
   Double_t        sbs_hcal_clus_blk_x[15];   //[Ndata.sbs.hcal.clus_blk.x]
   Int_t           Ndata_sbs_hcal_clus_blk_y;
   Double_t        sbs_hcal_clus_blk_y[15];   //[Ndata.sbs.hcal.clus_blk.y]
   Int_t           Ndata_sbs_hcal_goodblock_atime;
   Double_t        sbs_hcal_goodblock_atime[29];   //[Ndata.sbs.hcal.goodblock.atime]
   Int_t           Ndata_sbs_hcal_goodblock_cid;
   Double_t        sbs_hcal_goodblock_cid[29];   //[Ndata.sbs.hcal.goodblock.cid]
   Int_t           Ndata_sbs_hcal_goodblock_col;
   Double_t        sbs_hcal_goodblock_col[29];   //[Ndata.sbs.hcal.goodblock.col]
   Int_t           Ndata_sbs_hcal_goodblock_e;
   Double_t        sbs_hcal_goodblock_e[29];   //[Ndata.sbs.hcal.goodblock.e]
   Int_t           Ndata_sbs_hcal_goodblock_id;
   Double_t        sbs_hcal_goodblock_id[29];   //[Ndata.sbs.hcal.goodblock.id]
   Int_t           Ndata_sbs_hcal_goodblock_row;
   Double_t        sbs_hcal_goodblock_row[29];   //[Ndata.sbs.hcal.goodblock.row]
   Int_t           Ndata_sbs_hcal_goodblock_tdctime;
   Double_t        sbs_hcal_goodblock_tdctime[29];   //[Ndata.sbs.hcal.goodblock.tdctime]
   Int_t           Ndata_sbs_hcal_goodblock_tdctime_tw;
   Double_t        sbs_hcal_goodblock_tdctime_tw[29];   //[Ndata.sbs.hcal.goodblock.tdctime_tw]
   Int_t           Ndata_sbs_hcal_goodblock_x;
   Double_t        sbs_hcal_goodblock_x[29];   //[Ndata.sbs.hcal.goodblock.x]
   Int_t           Ndata_sbs_hcal_goodblock_y;
   Double_t        sbs_hcal_goodblock_y[29];   //[Ndata.sbs.hcal.goodblock.y]
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
   Int_t           Ndata_sbs_x_bcp;
   Double_t        sbs_x_bcp[35];   //[Ndata.sbs.x_bcp]
   Int_t           Ndata_sbs_x_fcp;
   Double_t        sbs_x_fcp[35];   //[Ndata.sbs.x_fcp]
   Int_t           Ndata_sbs_y_bcp;
   Double_t        sbs_y_bcp[35];   //[Ndata.sbs.y_bcp]
   Int_t           Ndata_sbs_y_fcp;
   Double_t        sbs_y_fcp[35];   //[Ndata.sbs.y_fcp]
   Int_t           Ndata_sbs_z_bcp;
   Double_t        sbs_z_bcp[35];   //[Ndata.sbs.z_bcp]
   Int_t           Ndata_sbs_z_fcp;
   Double_t        sbs_z_fcp[35];   //[Ndata.sbs.z_fcp]
   Double_t        FTROI_ECAL_energy;
   Double_t        FTROI_Eprime;
   Double_t        FTROI_ephi;
   Double_t        FTROI_etheta;
   Double_t        FTROI_good;
   Double_t        FTROI_pp;
   Double_t        FTROI_pphi;
   Double_t        FTROI_ptheta;
   Double_t        FTROI_xECAL_global;
   Double_t        FTROI_yECAL_global;
   Double_t        FTROI_zECAL_global;
   Double_t        MC_mc_ep;
   Double_t        MC_mc_epx;
   Double_t        MC_mc_epy;
   Double_t        MC_mc_epz;
   Double_t        MC_mc_fnucl;
   Double_t        MC_mc_np;
   Double_t        MC_mc_npx;
   Double_t        MC_mc_npy;
   Double_t        MC_mc_npz;
   Double_t        MC_mc_nucl;
   Double_t        MC_mc_omega;
   Double_t        MC_mc_sigma;
   Double_t        MC_mc_vx;
   Double_t        MC_mc_vy;
   Double_t        MC_mc_vz;
   Double_t        earm_ecal_againblk;
   Double_t        earm_ecal_atimeblk;
   Double_t        earm_ecal_colblk;
   Double_t        earm_ecal_e;
   Double_t        earm_ecal_eblk;
   Double_t        earm_ecal_idblk;
   Double_t        earm_ecal_index;
   Double_t        earm_ecal_nblk;
   Double_t        earm_ecal_nclus;
   Double_t        earm_ecal_ngoodADChits;
   Double_t        earm_ecal_rowblk;
   Double_t        earm_ecal_x;
   Double_t        earm_ecal_y;
   Double_t        earm_tr_n;
   Double_t        sbs_HCALdir_x;
   Double_t        sbs_HCALdir_y;
   Double_t        sbs_HCALdir_z;
   Double_t        sbs_HCALph_n;
   Double_t        sbs_HCALth_n;
   Double_t        sbs_gemFPP_hit_ngoodhits;
   Double_t        sbs_gemFPP_nlayershit;
   Double_t        sbs_gemFPP_nlayershitu;
   Double_t        sbs_gemFPP_nlayershituv;
   Double_t        sbs_gemFPP_nlayershitv;
   Double_t        sbs_gemFPP_track_besttrack;
   Double_t        sbs_gemFPP_track_ntrack;
   Double_t        sbs_gemFT_hit_ngoodhits;
   Double_t        sbs_gemFT_nlayershit;
   Double_t        sbs_gemFT_nlayershitu;
   Double_t        sbs_gemFT_nlayershituv;
   Double_t        sbs_gemFT_nlayershitv;
   Double_t        sbs_gemFT_track_besttrack;
   Double_t        sbs_gemFT_track_ntrack;
   Double_t        sbs_hcal_adctime;
   Double_t        sbs_hcal_againblk;
   Double_t        sbs_hcal_atimeblk;
   Double_t        sbs_hcal_colblk;
   Double_t        sbs_hcal_colblk_goodtdc;
   Double_t        sbs_hcal_e;
   Double_t        sbs_hcal_e_goodtdc;
   Double_t        sbs_hcal_eblk;
   Double_t        sbs_hcal_eblk_goodtdc;
   Double_t        sbs_hcal_idblk;
   Double_t        sbs_hcal_idblk_goodtdc;
   Double_t        sbs_hcal_index;
   Double_t        sbs_hcal_nblk;
   Double_t        sbs_hcal_nblk_goodtdc;
   Double_t        sbs_hcal_nclus;
   Double_t        sbs_hcal_ngoodADChits;
   Double_t        sbs_hcal_rowblk;
   Double_t        sbs_hcal_rowblk_goodtdc;
   Double_t        sbs_hcal_tdctime;
   Double_t        sbs_hcal_tdctime_tw;
   Double_t        sbs_hcal_tdctimeblk;
   Double_t        sbs_hcal_tdctimeblk_tw;
   Double_t        sbs_hcal_x;
   Double_t        sbs_hcal_y;
   Double_t        sbs_tr_n;
   Double_t        singletrack;
   Double_t        anytrack;
   Double_t        singletrack;
   Double_t        anytrack;
 //THaEvent        *Event_Branch;
   ULong64_t       fEvtHdr_fEvtTime;
   UInt_t          fEvtHdr_fEvtNum;
   UInt_t          fEvtHdr_fEvtType;
   UInt_t          fEvtHdr_fEvtLen;
   Int_t           fEvtHdr_fHelicity;
   UInt_t          fEvtHdr_fTrigBits;
   UInt_t          fEvtHdr_fRun;

   // List of branches
   TBranch        *b_Ndata_MC_bbgemhit_ptridx;   //!
   TBranch        *b_MC_bbgemhit_ptridx;   //!
   TBranch        *b_Ndata_MC_bbgemhit_sdtridx;   //!
   TBranch        *b_MC_bbgemhit_sdtridx;   //!
   TBranch        *b_Ndata_MC_bbgemtrack_ptridx;   //!
   TBranch        *b_MC_bbgemtrack_ptridx;   //!
   TBranch        *b_Ndata_MC_bbgemtrack_sdtridx;   //!
   TBranch        *b_MC_bbgemtrack_sdtridx;   //!
   TBranch        *b_Ndata_MC_bbhodohit_ptridx;   //!
   TBranch        *b_MC_bbhodohit_ptridx;   //!
   TBranch        *b_Ndata_MC_bbhodohit_sdtridx;   //!
   TBranch        *b_MC_bbhodohit_sdtridx;   //!
   TBranch        *b_Ndata_MC_bbpshit_ptridx;   //!
   TBranch        *b_MC_bbpshit_ptridx;   //!
   TBranch        *b_Ndata_MC_bbpshit_sdtridx;   //!
   TBranch        *b_MC_bbpshit_sdtridx;   //!
   TBranch        *b_Ndata_MC_bbshhit_ptridx;   //!
   TBranch        *b_MC_bbshhit_ptridx;   //!
   TBranch        *b_Ndata_MC_bbshhit_sdtridx;   //!
   TBranch        *b_MC_bbshhit_sdtridx;   //!
   TBranch        *b_Ndata_MC_hcalhit_ptridx;   //!
   TBranch        *b_MC_hcalhit_ptridx;   //!
   TBranch        *b_Ndata_MC_hcalhit_sdtridx;   //!
   TBranch        *b_MC_hcalhit_sdtridx;   //!
   TBranch        *b_Ndata_MC_ptrack_etot;   //!
   TBranch        *b_MC_ptrack_etot;   //!
   TBranch        *b_Ndata_MC_ptrack_momx;   //!
   TBranch        *b_MC_ptrack_momx;   //!
   TBranch        *b_Ndata_MC_ptrack_momy;   //!
   TBranch        *b_MC_ptrack_momy;   //!
   TBranch        *b_Ndata_MC_ptrack_momz;   //!
   TBranch        *b_MC_ptrack_momz;   //!
   TBranch        *b_Ndata_MC_ptrack_pid;   //!
   TBranch        *b_MC_ptrack_pid;   //!
   TBranch        *b_Ndata_MC_ptrack_polx;   //!
   TBranch        *b_MC_ptrack_polx;   //!
   TBranch        *b_Ndata_MC_ptrack_poly;   //!
   TBranch        *b_MC_ptrack_poly;   //!
   TBranch        *b_Ndata_MC_ptrack_polz;   //!
   TBranch        *b_MC_ptrack_polz;   //!
   TBranch        *b_Ndata_MC_ptrack_posx;   //!
   TBranch        *b_MC_ptrack_posx;   //!
   TBranch        *b_Ndata_MC_ptrack_posy;   //!
   TBranch        *b_MC_ptrack_posy;   //!
   TBranch        *b_Ndata_MC_ptrack_posz;   //!
   TBranch        *b_MC_ptrack_posz;   //!
   TBranch        *b_Ndata_MC_ptrack_t;   //!
   TBranch        *b_MC_ptrack_t;   //!
   TBranch        *b_Ndata_MC_ptrack_tid;   //!
   TBranch        *b_MC_ptrack_tid;   //!
   TBranch        *b_Ndata_MC_sdtrack_etot;   //!
   TBranch        *b_MC_sdtrack_etot;   //!
   TBranch        *b_Ndata_MC_sdtrack_mid;   //!
   TBranch        *b_MC_sdtrack_mid;   //!
   TBranch        *b_Ndata_MC_sdtrack_momx;   //!
   TBranch        *b_MC_sdtrack_momx;   //!
   TBranch        *b_Ndata_MC_sdtrack_momy;   //!
   TBranch        *b_MC_sdtrack_momy;   //!
   TBranch        *b_Ndata_MC_sdtrack_momz;   //!
   TBranch        *b_MC_sdtrack_momz;   //!
   TBranch        *b_Ndata_MC_sdtrack_pid;   //!
   TBranch        *b_MC_sdtrack_pid;   //!
   TBranch        *b_Ndata_MC_sdtrack_polx;   //!
   TBranch        *b_MC_sdtrack_polx;   //!
   TBranch        *b_Ndata_MC_sdtrack_poly;   //!
   TBranch        *b_MC_sdtrack_poly;   //!
   TBranch        *b_Ndata_MC_sdtrack_polz;   //!
   TBranch        *b_MC_sdtrack_polz;   //!
   TBranch        *b_Ndata_MC_sdtrack_posx;   //!
   TBranch        *b_MC_sdtrack_posx;   //!
   TBranch        *b_Ndata_MC_sdtrack_posy;   //!
   TBranch        *b_MC_sdtrack_posy;   //!
   TBranch        *b_Ndata_MC_sdtrack_posz;   //!
   TBranch        *b_MC_sdtrack_posz;   //!
   TBranch        *b_Ndata_MC_sdtrack_t;   //!
   TBranch        *b_MC_sdtrack_t;   //!
   TBranch        *b_Ndata_MC_sdtrack_tid;   //!
   TBranch        *b_MC_sdtrack_tid;   //!
   TBranch        *b_Ndata_MC_sdtrack_vEkin;   //!
   TBranch        *b_MC_sdtrack_vEkin;   //!
   TBranch        *b_Ndata_MC_sdtrack_vnx;   //!
   TBranch        *b_MC_sdtrack_vnx;   //!
   TBranch        *b_Ndata_MC_sdtrack_vny;   //!
   TBranch        *b_MC_sdtrack_vny;   //!
   TBranch        *b_Ndata_MC_sdtrack_vnz;   //!
   TBranch        *b_MC_sdtrack_vnz;   //!
   TBranch        *b_Ndata_MC_sdtrack_vx;   //!
   TBranch        *b_MC_sdtrack_vx;   //!
   TBranch        *b_Ndata_MC_sdtrack_vy;   //!
   TBranch        *b_MC_sdtrack_vy;   //!
   TBranch        *b_Ndata_MC_sdtrack_vz;   //!
   TBranch        *b_MC_sdtrack_vz;   //!
   TBranch        *b_Ndata_earm_ecal_a_amp_p;   //!
   TBranch        *b_earm_ecal_a_amp_p;   //!
   TBranch        *b_Ndata_earm_ecal_a_p;   //!
   TBranch        *b_earm_ecal_a_p;   //!
   TBranch        *b_Ndata_earm_ecal_a_time;   //!
   TBranch        *b_earm_ecal_a_time;   //!
   TBranch        *b_Ndata_earm_ecal_adccol;   //!
   TBranch        *b_earm_ecal_adccol;   //!
   TBranch        *b_Ndata_earm_ecal_adcelemID;   //!
   TBranch        *b_earm_ecal_adcelemID;   //!
   TBranch        *b_Ndata_earm_ecal_adcrow;   //!
   TBranch        *b_earm_ecal_adcrow;   //!
   TBranch        *b_Ndata_earm_ecal_adcxpos;   //!
   TBranch        *b_earm_ecal_adcxpos;   //!
   TBranch        *b_Ndata_earm_ecal_adcypos;   //!
   TBranch        *b_earm_ecal_adcypos;   //!
   TBranch        *b_Ndata_earm_ecal_clus_adctime;   //!
   TBranch        *b_earm_ecal_clus_adctime;   //!
   TBranch        *b_Ndata_earm_ecal_clus_again;   //!
   TBranch        *b_earm_ecal_clus_again;   //!
   TBranch        *b_Ndata_earm_ecal_clus_atimeblk;   //!
   TBranch        *b_earm_ecal_clus_atimeblk;   //!
   TBranch        *b_Ndata_earm_ecal_clus_col;   //!
   TBranch        *b_earm_ecal_clus_col;   //!
   TBranch        *b_Ndata_earm_ecal_clus_col_goodtdc;   //!
   TBranch        *b_earm_ecal_clus_col_goodtdc;   //!
   TBranch        *b_Ndata_earm_ecal_clus_e;   //!
   TBranch        *b_earm_ecal_clus_e;   //!
   TBranch        *b_Ndata_earm_ecal_clus_e_goodtdc;   //!
   TBranch        *b_earm_ecal_clus_e_goodtdc;   //!
   TBranch        *b_Ndata_earm_ecal_clus_eblk;   //!
   TBranch        *b_earm_ecal_clus_eblk;   //!
   TBranch        *b_Ndata_earm_ecal_clus_eblk_goodtdc;   //!
   TBranch        *b_earm_ecal_clus_eblk_goodtdc;   //!
   TBranch        *b_Ndata_earm_ecal_clus_id;   //!
   TBranch        *b_earm_ecal_clus_id;   //!
   TBranch        *b_Ndata_earm_ecal_clus_id_goodtdc;   //!
   TBranch        *b_earm_ecal_clus_id_goodtdc;   //!
   TBranch        *b_Ndata_earm_ecal_clus_nblk;   //!
   TBranch        *b_earm_ecal_clus_nblk;   //!
   TBranch        *b_Ndata_earm_ecal_clus_nblk_goodtdc;   //!
   TBranch        *b_earm_ecal_clus_nblk_goodtdc;   //!
   TBranch        *b_Ndata_earm_ecal_clus_row;   //!
   TBranch        *b_earm_ecal_clus_row;   //!
   TBranch        *b_Ndata_earm_ecal_clus_row_goodtdc;   //!
   TBranch        *b_earm_ecal_clus_row_goodtdc;   //!
   TBranch        *b_Ndata_earm_ecal_clus_tdctime;   //!
   TBranch        *b_earm_ecal_clus_tdctime;   //!
   TBranch        *b_Ndata_earm_ecal_clus_tdctime_tw;   //!
   TBranch        *b_earm_ecal_clus_tdctime_tw;   //!
   TBranch        *b_Ndata_earm_ecal_clus_tdctimeblk;   //!
   TBranch        *b_earm_ecal_clus_tdctimeblk;   //!
   TBranch        *b_Ndata_earm_ecal_clus_tdctimeblk_tw;   //!
   TBranch        *b_earm_ecal_clus_tdctimeblk_tw;   //!
   TBranch        *b_Ndata_earm_ecal_clus_x;   //!
   TBranch        *b_earm_ecal_clus_x;   //!
   TBranch        *b_Ndata_earm_ecal_clus_y;   //!
   TBranch        *b_earm_ecal_clus_y;   //!
   TBranch        *b_Ndata_earm_ecal_clus_blk_again;   //!
   TBranch        *b_earm_ecal_clus_blk_again;   //!
   TBranch        *b_Ndata_earm_ecal_clus_blk_atime;   //!
   TBranch        *b_earm_ecal_clus_blk_atime;   //!
   TBranch        *b_Ndata_earm_ecal_clus_blk_col;   //!
   TBranch        *b_earm_ecal_clus_blk_col;   //!
   TBranch        *b_Ndata_earm_ecal_clus_blk_e;   //!
   TBranch        *b_earm_ecal_clus_blk_e;   //!
   TBranch        *b_Ndata_earm_ecal_clus_blk_id;   //!
   TBranch        *b_earm_ecal_clus_blk_id;   //!
   TBranch        *b_Ndata_earm_ecal_clus_blk_row;   //!
   TBranch        *b_earm_ecal_clus_blk_row;   //!
   TBranch        *b_Ndata_earm_ecal_clus_blk_tdctime;   //!
   TBranch        *b_earm_ecal_clus_blk_tdctime;   //!
   TBranch        *b_Ndata_earm_ecal_clus_blk_tdctime_tw;   //!
   TBranch        *b_earm_ecal_clus_blk_tdctime_tw;   //!
   TBranch        *b_Ndata_earm_ecal_clus_blk_x;   //!
   TBranch        *b_earm_ecal_clus_blk_x;   //!
   TBranch        *b_Ndata_earm_ecal_clus_blk_y;   //!
   TBranch        *b_earm_ecal_clus_blk_y;   //!
   TBranch        *b_Ndata_earm_ecal_ped;   //!
   TBranch        *b_earm_ecal_ped;   //!
   TBranch        *b_Ndata_earm_tr_beta;   //!
   TBranch        *b_earm_tr_beta;   //!
   TBranch        *b_Ndata_earm_tr_chi2;   //!
   TBranch        *b_earm_tr_chi2;   //!
   TBranch        *b_Ndata_earm_tr_d_ph;   //!
   TBranch        *b_earm_tr_d_ph;   //!
   TBranch        *b_Ndata_earm_tr_d_th;   //!
   TBranch        *b_earm_tr_d_th;   //!
   TBranch        *b_Ndata_earm_tr_d_x;   //!
   TBranch        *b_earm_tr_d_x;   //!
   TBranch        *b_Ndata_earm_tr_d_y;   //!
   TBranch        *b_earm_tr_d_y;   //!
   TBranch        *b_Ndata_earm_tr_dbeta;   //!
   TBranch        *b_earm_tr_dbeta;   //!
   TBranch        *b_Ndata_earm_tr_dtime;   //!
   TBranch        *b_earm_tr_dtime;   //!
   TBranch        *b_Ndata_earm_tr_flag;   //!
   TBranch        *b_earm_tr_flag;   //!
   TBranch        *b_Ndata_earm_tr_ndof;   //!
   TBranch        *b_earm_tr_ndof;   //!
   TBranch        *b_Ndata_earm_tr_p;   //!
   TBranch        *b_earm_tr_p;   //!
   TBranch        *b_Ndata_earm_tr_pathl;   //!
   TBranch        *b_earm_tr_pathl;   //!
   TBranch        *b_Ndata_earm_tr_ph;   //!
   TBranch        *b_earm_tr_ph;   //!
   TBranch        *b_Ndata_earm_tr_px;   //!
   TBranch        *b_earm_tr_px;   //!
   TBranch        *b_Ndata_earm_tr_py;   //!
   TBranch        *b_earm_tr_py;   //!
   TBranch        *b_Ndata_earm_tr_pz;   //!
   TBranch        *b_earm_tr_pz;   //!
   TBranch        *b_Ndata_earm_tr_r_ph;   //!
   TBranch        *b_earm_tr_r_ph;   //!
   TBranch        *b_Ndata_earm_tr_r_th;   //!
   TBranch        *b_earm_tr_r_th;   //!
   TBranch        *b_Ndata_earm_tr_r_x;   //!
   TBranch        *b_earm_tr_r_x;   //!
   TBranch        *b_Ndata_earm_tr_r_y;   //!
   TBranch        *b_earm_tr_r_y;   //!
   TBranch        *b_Ndata_earm_tr_tg_dp;   //!
   TBranch        *b_earm_tr_tg_dp;   //!
   TBranch        *b_Ndata_earm_tr_tg_ph;   //!
   TBranch        *b_earm_tr_tg_ph;   //!
   TBranch        *b_Ndata_earm_tr_tg_th;   //!
   TBranch        *b_earm_tr_tg_th;   //!
   TBranch        *b_Ndata_earm_tr_tg_y;   //!
   TBranch        *b_earm_tr_tg_y;   //!
   TBranch        *b_Ndata_earm_tr_th;   //!
   TBranch        *b_earm_tr_th;   //!
   TBranch        *b_Ndata_earm_tr_time;   //!
   TBranch        *b_earm_tr_time;   //!
   TBranch        *b_Ndata_earm_tr_vx;   //!
   TBranch        *b_earm_tr_vx;   //!
   TBranch        *b_Ndata_earm_tr_vy;   //!
   TBranch        *b_earm_tr_vy;   //!
   TBranch        *b_Ndata_earm_tr_vz;   //!
   TBranch        *b_earm_tr_vz;   //!
   TBranch        *b_Ndata_earm_tr_x;   //!
   TBranch        *b_earm_tr_x;   //!
   TBranch        *b_Ndata_earm_tr_y;   //!
   TBranch        *b_earm_tr_y;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCU;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCU;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCU_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCU_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCV;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCV;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCV_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCV_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCasym;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCasym;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCasym_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCasym_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCavg;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCavg;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCavg_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCavg_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac0_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac0_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac0_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac0_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac1_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac1_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac1_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac1_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac2_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac2_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac2_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac2_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac3_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac3_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac3_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac3_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac4_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac4_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac4_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac4_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac5_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac5_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCfrac5_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCfrac5_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxcomboUclust_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxcomboUclust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxcomboVclust_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxcomboVclust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxsampU;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxsampU;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxsampUclust;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxsampUclust;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxsampUclust_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxsampUclust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxsampV;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxsampV;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxsampVclust;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxsampVclust;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxsampVclust_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxsampVclust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxstripU;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxstripU;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ADCmaxstripV;   //!
   TBranch        *b_sbs_gemFPP_hit_ADCmaxstripV;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_U;   //!
   TBranch        *b_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_U;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_V;   //!
   TBranch        *b_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_V;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_CM_GOOD_U;   //!
   TBranch        *b_sbs_gemFPP_hit_CM_GOOD_U;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_CM_GOOD_V;   //!
   TBranch        *b_sbs_gemFPP_hit_CM_GOOD_V;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC0_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC0_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC0_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC0_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC1_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC1_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC1_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC1_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC2_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC2_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC2_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC2_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC3_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC3_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC3_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC3_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC4_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC4_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC4_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC4_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC5_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC5_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADC5_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADC5_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADCmaxcomboU;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADCmaxcomboU;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADCmaxcomboV;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADCmaxcomboV;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADCmaxsampU;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADCmaxsampU;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADCmaxsampV;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADCmaxsampV;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADCmaxstripU;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADCmaxstripU;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_DeconvADCmaxstripV;   //!
   TBranch        *b_sbs_gemFPP_hit_DeconvADCmaxstripV;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ENABLE_CM_U;   //!
   TBranch        *b_sbs_gemFPP_hit_ENABLE_CM_U;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ENABLE_CM_V;   //!
   TBranch        *b_sbs_gemFPP_hit_ENABLE_CM_V;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_TSchi2_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_TSchi2_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_TSchi2_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_TSchi2_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_TSprob_Umax;   //!
   TBranch        *b_sbs_gemFPP_hit_TSprob_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_TSprob_Vmax;   //!
   TBranch        *b_sbs_gemFPP_hit_TSprob_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_Tavg;   //!
   TBranch        *b_sbs_gemFPP_hit_Tavg;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_Tavg_corr;   //!
   TBranch        *b_sbs_gemFPP_hit_Tavg_corr;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_Tavg_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_Tavg_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_Tavg_fit;   //!
   TBranch        *b_sbs_gemFPP_hit_Tavg_fit;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_Ugain;   //!
   TBranch        *b_sbs_gemFPP_hit_Ugain;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_Utime;   //!
   TBranch        *b_sbs_gemFPP_hit_Utime;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_UtimeDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_UtimeDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_UtimeFit;   //!
   TBranch        *b_sbs_gemFPP_hit_UtimeFit;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_UtimeMaxStrip;   //!
   TBranch        *b_sbs_gemFPP_hit_UtimeMaxStrip;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_UtimeMaxStripDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_UtimeMaxStripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_UtimeMaxStripFit;   //!
   TBranch        *b_sbs_gemFPP_hit_UtimeMaxStripFit;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_Vgain;   //!
   TBranch        *b_sbs_gemFPP_hit_Vgain;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_Vtime;   //!
   TBranch        *b_sbs_gemFPP_hit_Vtime;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_VtimeDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_VtimeDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_VtimeFit;   //!
   TBranch        *b_sbs_gemFPP_hit_VtimeFit;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_VtimeMaxStrip;   //!
   TBranch        *b_sbs_gemFPP_hit_VtimeMaxStrip;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_VtimeMaxStripDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_VtimeMaxStripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_VtimeMaxStripFit;   //!
   TBranch        *b_sbs_gemFPP_hit_VtimeMaxStripFit;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_adc_id_U;   //!
   TBranch        *b_sbs_gemFPP_hit_adc_id_U;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_adc_id_V;   //!
   TBranch        *b_sbs_gemFPP_hit_adc_id_V;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ccor_clust;   //!
   TBranch        *b_sbs_gemFPP_hit_ccor_clust;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ccor_clust_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ccor_clust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ccor_strip;   //!
   TBranch        *b_sbs_gemFPP_hit_ccor_strip;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ccor_strip_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_ccor_strip_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_crate_U;   //!
   TBranch        *b_sbs_gemFPP_hit_crate_U;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_crate_V;   //!
   TBranch        *b_sbs_gemFPP_hit_crate_V;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_deltat;   //!
   TBranch        *b_sbs_gemFPP_hit_deltat;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_deltat_deconv;   //!
   TBranch        *b_sbs_gemFPP_hit_deltat_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_deltat_fit;   //!
   TBranch        *b_sbs_gemFPP_hit_deltat_fit;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_eresidu;   //!
   TBranch        *b_sbs_gemFPP_hit_eresidu;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_eresidv;   //!
   TBranch        *b_sbs_gemFPP_hit_eresidv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_icombomaxUclustDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_icombomaxUclustDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_icombomaxUstripDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_icombomaxUstripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_icombomaxVclustDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_icombomaxVclustDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_icombomaxVstripDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_icombomaxVstripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_isampmaxUclust;   //!
   TBranch        *b_sbs_gemFPP_hit_isampmaxUclust;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_isampmaxUclustDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_isampmaxUclustDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_isampmaxUstrip;   //!
   TBranch        *b_sbs_gemFPP_hit_isampmaxUstrip;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_isampmaxUstripDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_isampmaxUstripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_isampmaxVclust;   //!
   TBranch        *b_sbs_gemFPP_hit_isampmaxVclust;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_isampmaxVclustDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_isampmaxVclustDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_isampmaxVstrip;   //!
   TBranch        *b_sbs_gemFPP_hit_isampmaxVstrip;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_isampmaxVstripDeconv;   //!
   TBranch        *b_sbs_gemFPP_hit_isampmaxVstripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_layer;   //!
   TBranch        *b_sbs_gemFPP_hit_layer;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_module;   //!
   TBranch        *b_sbs_gemFPP_hit_module;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_mpd_U;   //!
   TBranch        *b_sbs_gemFPP_hit_mpd_U;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_mpd_V;   //!
   TBranch        *b_sbs_gemFPP_hit_mpd_V;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_nstripu;   //!
   TBranch        *b_sbs_gemFPP_hit_nstripu;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_nstripv;   //!
   TBranch        *b_sbs_gemFPP_hit_nstripv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_residu;   //!
   TBranch        *b_sbs_gemFPP_hit_residu;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_residv;   //!
   TBranch        *b_sbs_gemFPP_hit_residv;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_trackindex;   //!
   TBranch        *b_sbs_gemFPP_hit_trackindex;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_u;   //!
   TBranch        *b_sbs_gemFPP_hit_u;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_umoment;   //!
   TBranch        *b_sbs_gemFPP_hit_umoment;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_usigma;   //!
   TBranch        *b_sbs_gemFPP_hit_usigma;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ustriphi;   //!
   TBranch        *b_sbs_gemFPP_hit_ustriphi;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ustriplo;   //!
   TBranch        *b_sbs_gemFPP_hit_ustriplo;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ustripmax;   //!
   TBranch        *b_sbs_gemFPP_hit_ustripmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_v;   //!
   TBranch        *b_sbs_gemFPP_hit_v;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_vmoment;   //!
   TBranch        *b_sbs_gemFPP_hit_vmoment;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_vsigma;   //!
   TBranch        *b_sbs_gemFPP_hit_vsigma;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_vstriphi;   //!
   TBranch        *b_sbs_gemFPP_hit_vstriphi;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_vstriplo;   //!
   TBranch        *b_sbs_gemFPP_hit_vstriplo;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_vstripmax;   //!
   TBranch        *b_sbs_gemFPP_hit_vstripmax;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_xglobal;   //!
   TBranch        *b_sbs_gemFPP_hit_xglobal;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_xlocal;   //!
   TBranch        *b_sbs_gemFPP_hit_xlocal;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_yglobal;   //!
   TBranch        *b_sbs_gemFPP_hit_yglobal;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_ylocal;   //!
   TBranch        *b_sbs_gemFPP_hit_ylocal;   //!
   TBranch        *b_Ndata_sbs_gemFPP_hit_zglobal;   //!
   TBranch        *b_sbs_gemFPP_hit_zglobal;   //!
   TBranch        *b_Ndata_sbs_gemFPP_n2Dhit_layer;   //!
   TBranch        *b_sbs_gemFPP_n2Dhit_layer;   //!
   TBranch        *b_Ndata_sbs_gemFPP_nclustu_layer;   //!
   TBranch        *b_sbs_gemFPP_nclustu_layer;   //!
   TBranch        *b_Ndata_sbs_gemFPP_nclustv_layer;   //!
   TBranch        *b_sbs_gemFPP_nclustv_layer;   //!
   TBranch        *b_Ndata_sbs_gemFPP_nstripsu_layer;   //!
   TBranch        *b_sbs_gemFPP_nstripsu_layer;   //!
   TBranch        *b_Ndata_sbs_gemFPP_nstripsv_layer;   //!
   TBranch        *b_sbs_gemFPP_nstripsv_layer;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_chi2ndf;   //!
   TBranch        *b_sbs_gemFPP_track_chi2ndf;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_chi2ndf_hitquality;   //!
   TBranch        *b_sbs_gemFPP_track_chi2ndf_hitquality;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_ngoodhits;   //!
   TBranch        *b_sbs_gemFPP_track_ngoodhits;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_nhits;   //!
   TBranch        *b_sbs_gemFPP_track_nhits;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_phi;   //!
   TBranch        *b_sbs_gemFPP_track_phi;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_sclose;   //!
   TBranch        *b_sbs_gemFPP_track_sclose;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_t0;   //!
   TBranch        *b_sbs_gemFPP_track_t0;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_theta;   //!
   TBranch        *b_sbs_gemFPP_track_theta;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_x;   //!
   TBranch        *b_sbs_gemFPP_track_x;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_xp;   //!
   TBranch        *b_sbs_gemFPP_track_xp;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_y;   //!
   TBranch        *b_sbs_gemFPP_track_y;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_yp;   //!
   TBranch        *b_sbs_gemFPP_track_yp;   //!
   TBranch        *b_Ndata_sbs_gemFPP_track_zclose;   //!
   TBranch        *b_sbs_gemFPP_track_zclose;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCU;   //!
   TBranch        *b_sbs_gemFT_hit_ADCU;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCU_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ADCU_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCV;   //!
   TBranch        *b_sbs_gemFT_hit_ADCV;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCV_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ADCV_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCasym;   //!
   TBranch        *b_sbs_gemFT_hit_ADCasym;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCasym_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ADCasym_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCavg;   //!
   TBranch        *b_sbs_gemFT_hit_ADCavg;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCavg_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ADCavg_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac0_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac0_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac0_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac0_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac1_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac1_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac1_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac1_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac2_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac2_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac2_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac2_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac3_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac3_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac3_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac3_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac4_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac4_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac4_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac4_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac5_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac5_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCfrac5_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_ADCfrac5_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxcomboUclust_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxcomboUclust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxcomboVclust_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxcomboVclust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxsampU;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxsampU;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxsampUclust;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxsampUclust;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxsampUclust_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxsampUclust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxsampV;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxsampV;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxsampVclust;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxsampVclust;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxsampVclust_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxsampVclust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxstripU;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxstripU;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ADCmaxstripV;   //!
   TBranch        *b_sbs_gemFT_hit_ADCmaxstripV;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_BUILD_ALL_SAMPLES_U;   //!
   TBranch        *b_sbs_gemFT_hit_BUILD_ALL_SAMPLES_U;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_BUILD_ALL_SAMPLES_V;   //!
   TBranch        *b_sbs_gemFT_hit_BUILD_ALL_SAMPLES_V;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_CM_GOOD_U;   //!
   TBranch        *b_sbs_gemFT_hit_CM_GOOD_U;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_CM_GOOD_V;   //!
   TBranch        *b_sbs_gemFT_hit_CM_GOOD_V;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC0_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC0_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC0_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC0_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC1_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC1_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC1_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC1_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC2_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC2_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC2_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC2_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC3_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC3_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC3_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC3_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC4_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC4_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC4_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC4_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC5_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC5_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADC5_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADC5_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADCmaxcomboU;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADCmaxcomboU;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADCmaxcomboV;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADCmaxcomboV;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADCmaxsampU;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADCmaxsampU;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADCmaxsampV;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADCmaxsampV;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADCmaxstripU;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADCmaxstripU;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_DeconvADCmaxstripV;   //!
   TBranch        *b_sbs_gemFT_hit_DeconvADCmaxstripV;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ENABLE_CM_U;   //!
   TBranch        *b_sbs_gemFT_hit_ENABLE_CM_U;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ENABLE_CM_V;   //!
   TBranch        *b_sbs_gemFT_hit_ENABLE_CM_V;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_TSchi2_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_TSchi2_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_TSchi2_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_TSchi2_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_TSprob_Umax;   //!
   TBranch        *b_sbs_gemFT_hit_TSprob_Umax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_TSprob_Vmax;   //!
   TBranch        *b_sbs_gemFT_hit_TSprob_Vmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_Tavg;   //!
   TBranch        *b_sbs_gemFT_hit_Tavg;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_Tavg_corr;   //!
   TBranch        *b_sbs_gemFT_hit_Tavg_corr;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_Tavg_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_Tavg_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_Tavg_fit;   //!
   TBranch        *b_sbs_gemFT_hit_Tavg_fit;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_Ugain;   //!
   TBranch        *b_sbs_gemFT_hit_Ugain;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_Utime;   //!
   TBranch        *b_sbs_gemFT_hit_Utime;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_UtimeDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_UtimeDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_UtimeFit;   //!
   TBranch        *b_sbs_gemFT_hit_UtimeFit;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_UtimeMaxStrip;   //!
   TBranch        *b_sbs_gemFT_hit_UtimeMaxStrip;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_UtimeMaxStripDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_UtimeMaxStripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_UtimeMaxStripFit;   //!
   TBranch        *b_sbs_gemFT_hit_UtimeMaxStripFit;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_Vgain;   //!
   TBranch        *b_sbs_gemFT_hit_Vgain;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_Vtime;   //!
   TBranch        *b_sbs_gemFT_hit_Vtime;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_VtimeDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_VtimeDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_VtimeFit;   //!
   TBranch        *b_sbs_gemFT_hit_VtimeFit;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_VtimeMaxStrip;   //!
   TBranch        *b_sbs_gemFT_hit_VtimeMaxStrip;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_VtimeMaxStripDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_VtimeMaxStripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_VtimeMaxStripFit;   //!
   TBranch        *b_sbs_gemFT_hit_VtimeMaxStripFit;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_adc_id_U;   //!
   TBranch        *b_sbs_gemFT_hit_adc_id_U;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_adc_id_V;   //!
   TBranch        *b_sbs_gemFT_hit_adc_id_V;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ccor_clust;   //!
   TBranch        *b_sbs_gemFT_hit_ccor_clust;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ccor_clust_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ccor_clust_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ccor_strip;   //!
   TBranch        *b_sbs_gemFT_hit_ccor_strip;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ccor_strip_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_ccor_strip_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_crate_U;   //!
   TBranch        *b_sbs_gemFT_hit_crate_U;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_crate_V;   //!
   TBranch        *b_sbs_gemFT_hit_crate_V;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_deltat;   //!
   TBranch        *b_sbs_gemFT_hit_deltat;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_deltat_deconv;   //!
   TBranch        *b_sbs_gemFT_hit_deltat_deconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_deltat_fit;   //!
   TBranch        *b_sbs_gemFT_hit_deltat_fit;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_eresidu;   //!
   TBranch        *b_sbs_gemFT_hit_eresidu;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_eresidv;   //!
   TBranch        *b_sbs_gemFT_hit_eresidv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_icombomaxUclustDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_icombomaxUclustDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_icombomaxUstripDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_icombomaxUstripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_icombomaxVclustDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_icombomaxVclustDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_icombomaxVstripDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_icombomaxVstripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_isampmaxUclust;   //!
   TBranch        *b_sbs_gemFT_hit_isampmaxUclust;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_isampmaxUclustDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_isampmaxUclustDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_isampmaxUstrip;   //!
   TBranch        *b_sbs_gemFT_hit_isampmaxUstrip;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_isampmaxUstripDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_isampmaxUstripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_isampmaxVclust;   //!
   TBranch        *b_sbs_gemFT_hit_isampmaxVclust;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_isampmaxVclustDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_isampmaxVclustDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_isampmaxVstrip;   //!
   TBranch        *b_sbs_gemFT_hit_isampmaxVstrip;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_isampmaxVstripDeconv;   //!
   TBranch        *b_sbs_gemFT_hit_isampmaxVstripDeconv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_layer;   //!
   TBranch        *b_sbs_gemFT_hit_layer;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_module;   //!
   TBranch        *b_sbs_gemFT_hit_module;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_mpd_U;   //!
   TBranch        *b_sbs_gemFT_hit_mpd_U;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_mpd_V;   //!
   TBranch        *b_sbs_gemFT_hit_mpd_V;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_nstripu;   //!
   TBranch        *b_sbs_gemFT_hit_nstripu;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_nstripv;   //!
   TBranch        *b_sbs_gemFT_hit_nstripv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_residu;   //!
   TBranch        *b_sbs_gemFT_hit_residu;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_residv;   //!
   TBranch        *b_sbs_gemFT_hit_residv;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_trackindex;   //!
   TBranch        *b_sbs_gemFT_hit_trackindex;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_u;   //!
   TBranch        *b_sbs_gemFT_hit_u;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_umoment;   //!
   TBranch        *b_sbs_gemFT_hit_umoment;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_usigma;   //!
   TBranch        *b_sbs_gemFT_hit_usigma;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ustriphi;   //!
   TBranch        *b_sbs_gemFT_hit_ustriphi;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ustriplo;   //!
   TBranch        *b_sbs_gemFT_hit_ustriplo;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ustripmax;   //!
   TBranch        *b_sbs_gemFT_hit_ustripmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_v;   //!
   TBranch        *b_sbs_gemFT_hit_v;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_vmoment;   //!
   TBranch        *b_sbs_gemFT_hit_vmoment;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_vsigma;   //!
   TBranch        *b_sbs_gemFT_hit_vsigma;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_vstriphi;   //!
   TBranch        *b_sbs_gemFT_hit_vstriphi;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_vstriplo;   //!
   TBranch        *b_sbs_gemFT_hit_vstriplo;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_vstripmax;   //!
   TBranch        *b_sbs_gemFT_hit_vstripmax;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_xglobal;   //!
   TBranch        *b_sbs_gemFT_hit_xglobal;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_xlocal;   //!
   TBranch        *b_sbs_gemFT_hit_xlocal;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_yglobal;   //!
   TBranch        *b_sbs_gemFT_hit_yglobal;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_ylocal;   //!
   TBranch        *b_sbs_gemFT_hit_ylocal;   //!
   TBranch        *b_Ndata_sbs_gemFT_hit_zglobal;   //!
   TBranch        *b_sbs_gemFT_hit_zglobal;   //!
   TBranch        *b_Ndata_sbs_gemFT_n2Dhit_layer;   //!
   TBranch        *b_sbs_gemFT_n2Dhit_layer;   //!
   TBranch        *b_Ndata_sbs_gemFT_nclustu_layer;   //!
   TBranch        *b_sbs_gemFT_nclustu_layer;   //!
   TBranch        *b_Ndata_sbs_gemFT_nclustv_layer;   //!
   TBranch        *b_sbs_gemFT_nclustv_layer;   //!
   TBranch        *b_Ndata_sbs_gemFT_nstripsu_layer;   //!
   TBranch        *b_sbs_gemFT_nstripsu_layer;   //!
   TBranch        *b_Ndata_sbs_gemFT_nstripsv_layer;   //!
   TBranch        *b_sbs_gemFT_nstripsv_layer;   //!
   TBranch        *b_Ndata_sbs_gemFT_track_chi2ndf;   //!
   TBranch        *b_sbs_gemFT_track_chi2ndf;   //!
   TBranch        *b_Ndata_sbs_gemFT_track_chi2ndf_hitquality;   //!
   TBranch        *b_sbs_gemFT_track_chi2ndf_hitquality;   //!
   TBranch        *b_Ndata_sbs_gemFT_track_ngoodhits;   //!
   TBranch        *b_sbs_gemFT_track_ngoodhits;   //!
   TBranch        *b_Ndata_sbs_gemFT_track_nhits;   //!
   TBranch        *b_sbs_gemFT_track_nhits;   //!
   TBranch        *b_Ndata_sbs_gemFT_track_t0;   //!
   TBranch        *b_sbs_gemFT_track_t0;   //!
   TBranch        *b_Ndata_sbs_gemFT_track_x;   //!
   TBranch        *b_sbs_gemFT_track_x;   //!
   TBranch        *b_Ndata_sbs_gemFT_track_xp;   //!
   TBranch        *b_sbs_gemFT_track_xp;   //!
   TBranch        *b_Ndata_sbs_gemFT_track_y;   //!
   TBranch        *b_sbs_gemFT_track_y;   //!
   TBranch        *b_Ndata_sbs_gemFT_track_yp;   //!
   TBranch        *b_sbs_gemFT_track_yp;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_adctime;   //!
   TBranch        *b_sbs_hcal_clus_adctime;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_again;   //!
   TBranch        *b_sbs_hcal_clus_again;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_atimeblk;   //!
   TBranch        *b_sbs_hcal_clus_atimeblk;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_col;   //!
   TBranch        *b_sbs_hcal_clus_col;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_col_goodtdc;   //!
   TBranch        *b_sbs_hcal_clus_col_goodtdc;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_e;   //!
   TBranch        *b_sbs_hcal_clus_e;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_e_goodtdc;   //!
   TBranch        *b_sbs_hcal_clus_e_goodtdc;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_eblk;   //!
   TBranch        *b_sbs_hcal_clus_eblk;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_eblk_goodtdc;   //!
   TBranch        *b_sbs_hcal_clus_eblk_goodtdc;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_id;   //!
   TBranch        *b_sbs_hcal_clus_id;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_id_goodtdc;   //!
   TBranch        *b_sbs_hcal_clus_id_goodtdc;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_nblk;   //!
   TBranch        *b_sbs_hcal_clus_nblk;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_nblk_goodtdc;   //!
   TBranch        *b_sbs_hcal_clus_nblk_goodtdc;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_row;   //!
   TBranch        *b_sbs_hcal_clus_row;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_row_goodtdc;   //!
   TBranch        *b_sbs_hcal_clus_row_goodtdc;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_tdctime;   //!
   TBranch        *b_sbs_hcal_clus_tdctime;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_tdctime_tw;   //!
   TBranch        *b_sbs_hcal_clus_tdctime_tw;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_tdctimeblk;   //!
   TBranch        *b_sbs_hcal_clus_tdctimeblk;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_tdctimeblk_tw;   //!
   TBranch        *b_sbs_hcal_clus_tdctimeblk_tw;   //!
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
   TBranch        *b_Ndata_sbs_hcal_clus_blk_tdctime_tw;   //!
   TBranch        *b_sbs_hcal_clus_blk_tdctime_tw;   //!
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
   TBranch        *b_Ndata_sbs_hcal_goodblock_tdctime_tw;   //!
   TBranch        *b_sbs_hcal_goodblock_tdctime_tw;   //!
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
   TBranch        *b_FTROI_ECAL_energy;   //!
   TBranch        *b_FTROI_Eprime;   //!
   TBranch        *b_FTROI_ephi;   //!
   TBranch        *b_FTROI_etheta;   //!
   TBranch        *b_FTROI_good;   //!
   TBranch        *b_FTROI_pp;   //!
   TBranch        *b_FTROI_pphi;   //!
   TBranch        *b_FTROI_ptheta;   //!
   TBranch        *b_FTROI_xECAL_global;   //!
   TBranch        *b_FTROI_yECAL_global;   //!
   TBranch        *b_FTROI_zECAL_global;   //!
   TBranch        *b_MC_mc_ep;   //!
   TBranch        *b_MC_mc_epx;   //!
   TBranch        *b_MC_mc_epy;   //!
   TBranch        *b_MC_mc_epz;   //!
   TBranch        *b_MC_mc_fnucl;   //!
   TBranch        *b_MC_mc_np;   //!
   TBranch        *b_MC_mc_npx;   //!
   TBranch        *b_MC_mc_npy;   //!
   TBranch        *b_MC_mc_npz;   //!
   TBranch        *b_MC_mc_nucl;   //!
   TBranch        *b_MC_mc_omega;   //!
   TBranch        *b_MC_mc_sigma;   //!
   TBranch        *b_MC_mc_vx;   //!
   TBranch        *b_MC_mc_vy;   //!
   TBranch        *b_MC_mc_vz;   //!
   TBranch        *b_earm_ecal_againblk;   //!
   TBranch        *b_earm_ecal_atimeblk;   //!
   TBranch        *b_earm_ecal_colblk;   //!
   TBranch        *b_earm_ecal_e;   //!
   TBranch        *b_earm_ecal_eblk;   //!
   TBranch        *b_earm_ecal_idblk;   //!
   TBranch        *b_earm_ecal_index;   //!
   TBranch        *b_earm_ecal_nblk;   //!
   TBranch        *b_earm_ecal_nclus;   //!
   TBranch        *b_earm_ecal_ngoodADChits;   //!
   TBranch        *b_earm_ecal_rowblk;   //!
   TBranch        *b_earm_ecal_x;   //!
   TBranch        *b_earm_ecal_y;   //!
   TBranch        *b_earm_tr_n;   //!
   TBranch        *b_sbs_HCALdir_x;   //!
   TBranch        *b_sbs_HCALdir_y;   //!
   TBranch        *b_sbs_HCALdir_z;   //!
   TBranch        *b_sbs_HCALph_n;   //!
   TBranch        *b_sbs_HCALth_n;   //!
   TBranch        *b_sbs_gemFPP_hit_ngoodhits;   //!
   TBranch        *b_sbs_gemFPP_nlayershit;   //!
   TBranch        *b_sbs_gemFPP_nlayershitu;   //!
   TBranch        *b_sbs_gemFPP_nlayershituv;   //!
   TBranch        *b_sbs_gemFPP_nlayershitv;   //!
   TBranch        *b_sbs_gemFPP_track_besttrack;   //!
   TBranch        *b_sbs_gemFPP_track_ntrack;   //!
   TBranch        *b_sbs_gemFT_hit_ngoodhits;   //!
   TBranch        *b_sbs_gemFT_nlayershit;   //!
   TBranch        *b_sbs_gemFT_nlayershitu;   //!
   TBranch        *b_sbs_gemFT_nlayershituv;   //!
   TBranch        *b_sbs_gemFT_nlayershitv;   //!
   TBranch        *b_sbs_gemFT_track_besttrack;   //!
   TBranch        *b_sbs_gemFT_track_ntrack;   //!
   TBranch        *b_sbs_hcal_adctime;   //!
   TBranch        *b_sbs_hcal_againblk;   //!
   TBranch        *b_sbs_hcal_atimeblk;   //!
   TBranch        *b_sbs_hcal_colblk;   //!
   TBranch        *b_sbs_hcal_colblk_goodtdc;   //!
   TBranch        *b_sbs_hcal_e;   //!
   TBranch        *b_sbs_hcal_e_goodtdc;   //!
   TBranch        *b_sbs_hcal_eblk;   //!
   TBranch        *b_sbs_hcal_eblk_goodtdc;   //!
   TBranch        *b_sbs_hcal_idblk;   //!
   TBranch        *b_sbs_hcal_idblk_goodtdc;   //!
   TBranch        *b_sbs_hcal_index;   //!
   TBranch        *b_sbs_hcal_nblk;   //!
   TBranch        *b_sbs_hcal_nblk_goodtdc;   //!
   TBranch        *b_sbs_hcal_nclus;   //!
   TBranch        *b_sbs_hcal_ngoodADChits;   //!
   TBranch        *b_sbs_hcal_rowblk;   //!
   TBranch        *b_sbs_hcal_rowblk_goodtdc;   //!
   TBranch        *b_sbs_hcal_tdctime;   //!
   TBranch        *b_sbs_hcal_tdctime_tw;   //!
   TBranch        *b_sbs_hcal_tdctimeblk;   //!
   TBranch        *b_sbs_hcal_tdctimeblk_tw;   //!
   TBranch        *b_sbs_hcal_x;   //!
   TBranch        *b_sbs_hcal_y;   //!
   TBranch        *b_sbs_tr_n;   //!
   TBranch        *b_singletrack;   //!
   TBranch        *b_anytrack;   //!
   TBranch        *b_singletrack;   //!
   TBranch        *b_anytrack;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtTime;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtNum;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtType;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtLen;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fHelicity;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fTrigBits;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fRun;   //!

   gep_tree(TTree *tree=0);
   virtual ~gep_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef gep_tree_cxx
gep_tree::gep_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("T",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("T","");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job1.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job10.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job100.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job11.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job12.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job13.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job14.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job15.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job16.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job17.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job18.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job19.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job2.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job20.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job21.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job22.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job23.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job24.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job25.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job26.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job27.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job28.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job29.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job3.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job30.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job31.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job32.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job33.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job34.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job35.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job36.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job37.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job38.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job39.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job4.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job40.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job41.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job42.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job43.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job44.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job45.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job46.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job47.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job48.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job49.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job5.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job50.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job51.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job52.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job53.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job54.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job55.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job56.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job57.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job58.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job59.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job6.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job60.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job61.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job62.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job63.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job64.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job65.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job66.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job67.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job68.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job69.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job7.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job70.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job71.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job72.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job73.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job74.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job75.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job76.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job77.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job78.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job79.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job8.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job80.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job81.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job82.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job83.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job84.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job85.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job86.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job87.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job88.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job89.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job9.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job90.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job91.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job92.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job93.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job94.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job95.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job96.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job97.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job98.root/T");
      chain->Add("$OUT_DIR/replayed_GEP1_elastic_job99.root/T");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

gep_tree::~gep_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gep_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gep_tree::LoadTree(Long64_t entry)
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

void gep_tree::Init(TTree *tree)
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

   fChain->SetBranchAddress("Ndata.MC.bbgemhit_ptridx", &Ndata_MC_bbgemhit_ptridx, &b_Ndata_MC_bbgemhit_ptridx);
   fChain->SetBranchAddress("MC.bbgemhit_ptridx", &MC_bbgemhit_ptridx, &b_MC_bbgemhit_ptridx);
   fChain->SetBranchAddress("Ndata.MC.bbgemhit_sdtridx", &Ndata_MC_bbgemhit_sdtridx, &b_Ndata_MC_bbgemhit_sdtridx);
   fChain->SetBranchAddress("MC.bbgemhit_sdtridx", &MC_bbgemhit_sdtridx, &b_MC_bbgemhit_sdtridx);
   fChain->SetBranchAddress("Ndata.MC.bbgemtrack_ptridx", &Ndata_MC_bbgemtrack_ptridx, &b_Ndata_MC_bbgemtrack_ptridx);
   fChain->SetBranchAddress("MC.bbgemtrack_ptridx", &MC_bbgemtrack_ptridx, &b_MC_bbgemtrack_ptridx);
   fChain->SetBranchAddress("Ndata.MC.bbgemtrack_sdtridx", &Ndata_MC_bbgemtrack_sdtridx, &b_Ndata_MC_bbgemtrack_sdtridx);
   fChain->SetBranchAddress("MC.bbgemtrack_sdtridx", &MC_bbgemtrack_sdtridx, &b_MC_bbgemtrack_sdtridx);
   fChain->SetBranchAddress("Ndata.MC.bbhodohit_ptridx", &Ndata_MC_bbhodohit_ptridx, &b_Ndata_MC_bbhodohit_ptridx);
   fChain->SetBranchAddress("MC.bbhodohit_ptridx", &MC_bbhodohit_ptridx, &b_MC_bbhodohit_ptridx);
   fChain->SetBranchAddress("Ndata.MC.bbhodohit_sdtridx", &Ndata_MC_bbhodohit_sdtridx, &b_Ndata_MC_bbhodohit_sdtridx);
   fChain->SetBranchAddress("MC.bbhodohit_sdtridx", &MC_bbhodohit_sdtridx, &b_MC_bbhodohit_sdtridx);
   fChain->SetBranchAddress("Ndata.MC.bbpshit_ptridx", &Ndata_MC_bbpshit_ptridx, &b_Ndata_MC_bbpshit_ptridx);
   fChain->SetBranchAddress("MC.bbpshit_ptridx", &MC_bbpshit_ptridx, &b_MC_bbpshit_ptridx);
   fChain->SetBranchAddress("Ndata.MC.bbpshit_sdtridx", &Ndata_MC_bbpshit_sdtridx, &b_Ndata_MC_bbpshit_sdtridx);
   fChain->SetBranchAddress("MC.bbpshit_sdtridx", &MC_bbpshit_sdtridx, &b_MC_bbpshit_sdtridx);
   fChain->SetBranchAddress("Ndata.MC.bbshhit_ptridx", &Ndata_MC_bbshhit_ptridx, &b_Ndata_MC_bbshhit_ptridx);
   fChain->SetBranchAddress("MC.bbshhit_ptridx", &MC_bbshhit_ptridx, &b_MC_bbshhit_ptridx);
   fChain->SetBranchAddress("Ndata.MC.bbshhit_sdtridx", &Ndata_MC_bbshhit_sdtridx, &b_Ndata_MC_bbshhit_sdtridx);
   fChain->SetBranchAddress("MC.bbshhit_sdtridx", &MC_bbshhit_sdtridx, &b_MC_bbshhit_sdtridx);
   fChain->SetBranchAddress("Ndata.MC.hcalhit_ptridx", &Ndata_MC_hcalhit_ptridx, &b_Ndata_MC_hcalhit_ptridx);
   fChain->SetBranchAddress("MC.hcalhit_ptridx", &MC_hcalhit_ptridx, &b_MC_hcalhit_ptridx);
   fChain->SetBranchAddress("Ndata.MC.hcalhit_sdtridx", &Ndata_MC_hcalhit_sdtridx, &b_Ndata_MC_hcalhit_sdtridx);
   fChain->SetBranchAddress("MC.hcalhit_sdtridx", &MC_hcalhit_sdtridx, &b_MC_hcalhit_sdtridx);
   fChain->SetBranchAddress("Ndata.MC.ptrack_etot", &Ndata_MC_ptrack_etot, &b_Ndata_MC_ptrack_etot);
   fChain->SetBranchAddress("MC.ptrack_etot", &MC_ptrack_etot, &b_MC_ptrack_etot);
   fChain->SetBranchAddress("Ndata.MC.ptrack_momx", &Ndata_MC_ptrack_momx, &b_Ndata_MC_ptrack_momx);
   fChain->SetBranchAddress("MC.ptrack_momx", &MC_ptrack_momx, &b_MC_ptrack_momx);
   fChain->SetBranchAddress("Ndata.MC.ptrack_momy", &Ndata_MC_ptrack_momy, &b_Ndata_MC_ptrack_momy);
   fChain->SetBranchAddress("MC.ptrack_momy", &MC_ptrack_momy, &b_MC_ptrack_momy);
   fChain->SetBranchAddress("Ndata.MC.ptrack_momz", &Ndata_MC_ptrack_momz, &b_Ndata_MC_ptrack_momz);
   fChain->SetBranchAddress("MC.ptrack_momz", &MC_ptrack_momz, &b_MC_ptrack_momz);
   fChain->SetBranchAddress("Ndata.MC.ptrack_pid", &Ndata_MC_ptrack_pid, &b_Ndata_MC_ptrack_pid);
   fChain->SetBranchAddress("MC.ptrack_pid", &MC_ptrack_pid, &b_MC_ptrack_pid);
   fChain->SetBranchAddress("Ndata.MC.ptrack_polx", &Ndata_MC_ptrack_polx, &b_Ndata_MC_ptrack_polx);
   fChain->SetBranchAddress("MC.ptrack_polx", &MC_ptrack_polx, &b_MC_ptrack_polx);
   fChain->SetBranchAddress("Ndata.MC.ptrack_poly", &Ndata_MC_ptrack_poly, &b_Ndata_MC_ptrack_poly);
   fChain->SetBranchAddress("MC.ptrack_poly", &MC_ptrack_poly, &b_MC_ptrack_poly);
   fChain->SetBranchAddress("Ndata.MC.ptrack_polz", &Ndata_MC_ptrack_polz, &b_Ndata_MC_ptrack_polz);
   fChain->SetBranchAddress("MC.ptrack_polz", &MC_ptrack_polz, &b_MC_ptrack_polz);
   fChain->SetBranchAddress("Ndata.MC.ptrack_posx", &Ndata_MC_ptrack_posx, &b_Ndata_MC_ptrack_posx);
   fChain->SetBranchAddress("MC.ptrack_posx", &MC_ptrack_posx, &b_MC_ptrack_posx);
   fChain->SetBranchAddress("Ndata.MC.ptrack_posy", &Ndata_MC_ptrack_posy, &b_Ndata_MC_ptrack_posy);
   fChain->SetBranchAddress("MC.ptrack_posy", &MC_ptrack_posy, &b_MC_ptrack_posy);
   fChain->SetBranchAddress("Ndata.MC.ptrack_posz", &Ndata_MC_ptrack_posz, &b_Ndata_MC_ptrack_posz);
   fChain->SetBranchAddress("MC.ptrack_posz", &MC_ptrack_posz, &b_MC_ptrack_posz);
   fChain->SetBranchAddress("Ndata.MC.ptrack_t", &Ndata_MC_ptrack_t, &b_Ndata_MC_ptrack_t);
   fChain->SetBranchAddress("MC.ptrack_t", &MC_ptrack_t, &b_MC_ptrack_t);
   fChain->SetBranchAddress("Ndata.MC.ptrack_tid", &Ndata_MC_ptrack_tid, &b_Ndata_MC_ptrack_tid);
   fChain->SetBranchAddress("MC.ptrack_tid", &MC_ptrack_tid, &b_MC_ptrack_tid);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_etot", &Ndata_MC_sdtrack_etot, &b_Ndata_MC_sdtrack_etot);
   fChain->SetBranchAddress("MC.sdtrack_etot", &MC_sdtrack_etot, &b_MC_sdtrack_etot);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_mid", &Ndata_MC_sdtrack_mid, &b_Ndata_MC_sdtrack_mid);
   fChain->SetBranchAddress("MC.sdtrack_mid", &MC_sdtrack_mid, &b_MC_sdtrack_mid);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_momx", &Ndata_MC_sdtrack_momx, &b_Ndata_MC_sdtrack_momx);
   fChain->SetBranchAddress("MC.sdtrack_momx", &MC_sdtrack_momx, &b_MC_sdtrack_momx);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_momy", &Ndata_MC_sdtrack_momy, &b_Ndata_MC_sdtrack_momy);
   fChain->SetBranchAddress("MC.sdtrack_momy", &MC_sdtrack_momy, &b_MC_sdtrack_momy);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_momz", &Ndata_MC_sdtrack_momz, &b_Ndata_MC_sdtrack_momz);
   fChain->SetBranchAddress("MC.sdtrack_momz", &MC_sdtrack_momz, &b_MC_sdtrack_momz);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_pid", &Ndata_MC_sdtrack_pid, &b_Ndata_MC_sdtrack_pid);
   fChain->SetBranchAddress("MC.sdtrack_pid", &MC_sdtrack_pid, &b_MC_sdtrack_pid);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_polx", &Ndata_MC_sdtrack_polx, &b_Ndata_MC_sdtrack_polx);
   fChain->SetBranchAddress("MC.sdtrack_polx", &MC_sdtrack_polx, &b_MC_sdtrack_polx);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_poly", &Ndata_MC_sdtrack_poly, &b_Ndata_MC_sdtrack_poly);
   fChain->SetBranchAddress("MC.sdtrack_poly", &MC_sdtrack_poly, &b_MC_sdtrack_poly);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_polz", &Ndata_MC_sdtrack_polz, &b_Ndata_MC_sdtrack_polz);
   fChain->SetBranchAddress("MC.sdtrack_polz", &MC_sdtrack_polz, &b_MC_sdtrack_polz);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_posx", &Ndata_MC_sdtrack_posx, &b_Ndata_MC_sdtrack_posx);
   fChain->SetBranchAddress("MC.sdtrack_posx", &MC_sdtrack_posx, &b_MC_sdtrack_posx);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_posy", &Ndata_MC_sdtrack_posy, &b_Ndata_MC_sdtrack_posy);
   fChain->SetBranchAddress("MC.sdtrack_posy", &MC_sdtrack_posy, &b_MC_sdtrack_posy);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_posz", &Ndata_MC_sdtrack_posz, &b_Ndata_MC_sdtrack_posz);
   fChain->SetBranchAddress("MC.sdtrack_posz", &MC_sdtrack_posz, &b_MC_sdtrack_posz);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_t", &Ndata_MC_sdtrack_t, &b_Ndata_MC_sdtrack_t);
   fChain->SetBranchAddress("MC.sdtrack_t", &MC_sdtrack_t, &b_MC_sdtrack_t);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_tid", &Ndata_MC_sdtrack_tid, &b_Ndata_MC_sdtrack_tid);
   fChain->SetBranchAddress("MC.sdtrack_tid", &MC_sdtrack_tid, &b_MC_sdtrack_tid);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_vEkin", &Ndata_MC_sdtrack_vEkin, &b_Ndata_MC_sdtrack_vEkin);
   fChain->SetBranchAddress("MC.sdtrack_vEkin", &MC_sdtrack_vEkin, &b_MC_sdtrack_vEkin);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_vnx", &Ndata_MC_sdtrack_vnx, &b_Ndata_MC_sdtrack_vnx);
   fChain->SetBranchAddress("MC.sdtrack_vnx", &MC_sdtrack_vnx, &b_MC_sdtrack_vnx);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_vny", &Ndata_MC_sdtrack_vny, &b_Ndata_MC_sdtrack_vny);
   fChain->SetBranchAddress("MC.sdtrack_vny", &MC_sdtrack_vny, &b_MC_sdtrack_vny);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_vnz", &Ndata_MC_sdtrack_vnz, &b_Ndata_MC_sdtrack_vnz);
   fChain->SetBranchAddress("MC.sdtrack_vnz", &MC_sdtrack_vnz, &b_MC_sdtrack_vnz);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_vx", &Ndata_MC_sdtrack_vx, &b_Ndata_MC_sdtrack_vx);
   fChain->SetBranchAddress("MC.sdtrack_vx", &MC_sdtrack_vx, &b_MC_sdtrack_vx);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_vy", &Ndata_MC_sdtrack_vy, &b_Ndata_MC_sdtrack_vy);
   fChain->SetBranchAddress("MC.sdtrack_vy", &MC_sdtrack_vy, &b_MC_sdtrack_vy);
   fChain->SetBranchAddress("Ndata.MC.sdtrack_vz", &Ndata_MC_sdtrack_vz, &b_Ndata_MC_sdtrack_vz);
   fChain->SetBranchAddress("MC.sdtrack_vz", &MC_sdtrack_vz, &b_MC_sdtrack_vz);
   fChain->SetBranchAddress("Ndata.earm.ecal.a_amp_p", &Ndata_earm_ecal_a_amp_p, &b_Ndata_earm_ecal_a_amp_p);
   fChain->SetBranchAddress("earm.ecal.a_amp_p", earm_ecal_a_amp_p, &b_earm_ecal_a_amp_p);
   fChain->SetBranchAddress("Ndata.earm.ecal.a_p", &Ndata_earm_ecal_a_p, &b_Ndata_earm_ecal_a_p);
   fChain->SetBranchAddress("earm.ecal.a_p", earm_ecal_a_p, &b_earm_ecal_a_p);
   fChain->SetBranchAddress("Ndata.earm.ecal.a_time", &Ndata_earm_ecal_a_time, &b_Ndata_earm_ecal_a_time);
   fChain->SetBranchAddress("earm.ecal.a_time", earm_ecal_a_time, &b_earm_ecal_a_time);
   fChain->SetBranchAddress("Ndata.earm.ecal.adccol", &Ndata_earm_ecal_adccol, &b_Ndata_earm_ecal_adccol);
   fChain->SetBranchAddress("earm.ecal.adccol", earm_ecal_adccol, &b_earm_ecal_adccol);
   fChain->SetBranchAddress("Ndata.earm.ecal.adcelemID", &Ndata_earm_ecal_adcelemID, &b_Ndata_earm_ecal_adcelemID);
   fChain->SetBranchAddress("earm.ecal.adcelemID", earm_ecal_adcelemID, &b_earm_ecal_adcelemID);
   fChain->SetBranchAddress("Ndata.earm.ecal.adcrow", &Ndata_earm_ecal_adcrow, &b_Ndata_earm_ecal_adcrow);
   fChain->SetBranchAddress("earm.ecal.adcrow", earm_ecal_adcrow, &b_earm_ecal_adcrow);
   fChain->SetBranchAddress("Ndata.earm.ecal.adcxpos", &Ndata_earm_ecal_adcxpos, &b_Ndata_earm_ecal_adcxpos);
   fChain->SetBranchAddress("earm.ecal.adcxpos", earm_ecal_adcxpos, &b_earm_ecal_adcxpos);
   fChain->SetBranchAddress("Ndata.earm.ecal.adcypos", &Ndata_earm_ecal_adcypos, &b_Ndata_earm_ecal_adcypos);
   fChain->SetBranchAddress("earm.ecal.adcypos", earm_ecal_adcypos, &b_earm_ecal_adcypos);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.adctime", &Ndata_earm_ecal_clus_adctime, &b_Ndata_earm_ecal_clus_adctime);
   fChain->SetBranchAddress("earm.ecal.clus.adctime", earm_ecal_clus_adctime, &b_earm_ecal_clus_adctime);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.again", &Ndata_earm_ecal_clus_again, &b_Ndata_earm_ecal_clus_again);
   fChain->SetBranchAddress("earm.ecal.clus.again", earm_ecal_clus_again, &b_earm_ecal_clus_again);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.atimeblk", &Ndata_earm_ecal_clus_atimeblk, &b_Ndata_earm_ecal_clus_atimeblk);
   fChain->SetBranchAddress("earm.ecal.clus.atimeblk", earm_ecal_clus_atimeblk, &b_earm_ecal_clus_atimeblk);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.col", &Ndata_earm_ecal_clus_col, &b_Ndata_earm_ecal_clus_col);
   fChain->SetBranchAddress("earm.ecal.clus.col", earm_ecal_clus_col, &b_earm_ecal_clus_col);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.col_goodtdc", &Ndata_earm_ecal_clus_col_goodtdc, &b_Ndata_earm_ecal_clus_col_goodtdc);
   fChain->SetBranchAddress("earm.ecal.clus.col_goodtdc", earm_ecal_clus_col_goodtdc, &b_earm_ecal_clus_col_goodtdc);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.e", &Ndata_earm_ecal_clus_e, &b_Ndata_earm_ecal_clus_e);
   fChain->SetBranchAddress("earm.ecal.clus.e", earm_ecal_clus_e, &b_earm_ecal_clus_e);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.e_goodtdc", &Ndata_earm_ecal_clus_e_goodtdc, &b_Ndata_earm_ecal_clus_e_goodtdc);
   fChain->SetBranchAddress("earm.ecal.clus.e_goodtdc", earm_ecal_clus_e_goodtdc, &b_earm_ecal_clus_e_goodtdc);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.eblk", &Ndata_earm_ecal_clus_eblk, &b_Ndata_earm_ecal_clus_eblk);
   fChain->SetBranchAddress("earm.ecal.clus.eblk", earm_ecal_clus_eblk, &b_earm_ecal_clus_eblk);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.eblk_goodtdc", &Ndata_earm_ecal_clus_eblk_goodtdc, &b_Ndata_earm_ecal_clus_eblk_goodtdc);
   fChain->SetBranchAddress("earm.ecal.clus.eblk_goodtdc", earm_ecal_clus_eblk_goodtdc, &b_earm_ecal_clus_eblk_goodtdc);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.id", &Ndata_earm_ecal_clus_id, &b_Ndata_earm_ecal_clus_id);
   fChain->SetBranchAddress("earm.ecal.clus.id", earm_ecal_clus_id, &b_earm_ecal_clus_id);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.id_goodtdc", &Ndata_earm_ecal_clus_id_goodtdc, &b_Ndata_earm_ecal_clus_id_goodtdc);
   fChain->SetBranchAddress("earm.ecal.clus.id_goodtdc", earm_ecal_clus_id_goodtdc, &b_earm_ecal_clus_id_goodtdc);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.nblk", &Ndata_earm_ecal_clus_nblk, &b_Ndata_earm_ecal_clus_nblk);
   fChain->SetBranchAddress("earm.ecal.clus.nblk", earm_ecal_clus_nblk, &b_earm_ecal_clus_nblk);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.nblk_goodtdc", &Ndata_earm_ecal_clus_nblk_goodtdc, &b_Ndata_earm_ecal_clus_nblk_goodtdc);
   fChain->SetBranchAddress("earm.ecal.clus.nblk_goodtdc", earm_ecal_clus_nblk_goodtdc, &b_earm_ecal_clus_nblk_goodtdc);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.row", &Ndata_earm_ecal_clus_row, &b_Ndata_earm_ecal_clus_row);
   fChain->SetBranchAddress("earm.ecal.clus.row", earm_ecal_clus_row, &b_earm_ecal_clus_row);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.row_goodtdc", &Ndata_earm_ecal_clus_row_goodtdc, &b_Ndata_earm_ecal_clus_row_goodtdc);
   fChain->SetBranchAddress("earm.ecal.clus.row_goodtdc", earm_ecal_clus_row_goodtdc, &b_earm_ecal_clus_row_goodtdc);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.tdctime", &Ndata_earm_ecal_clus_tdctime, &b_Ndata_earm_ecal_clus_tdctime);
   fChain->SetBranchAddress("earm.ecal.clus.tdctime", earm_ecal_clus_tdctime, &b_earm_ecal_clus_tdctime);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.tdctime_tw", &Ndata_earm_ecal_clus_tdctime_tw, &b_Ndata_earm_ecal_clus_tdctime_tw);
   fChain->SetBranchAddress("earm.ecal.clus.tdctime_tw", earm_ecal_clus_tdctime_tw, &b_earm_ecal_clus_tdctime_tw);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.tdctimeblk", &Ndata_earm_ecal_clus_tdctimeblk, &b_Ndata_earm_ecal_clus_tdctimeblk);
   fChain->SetBranchAddress("earm.ecal.clus.tdctimeblk", earm_ecal_clus_tdctimeblk, &b_earm_ecal_clus_tdctimeblk);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.tdctimeblk_tw", &Ndata_earm_ecal_clus_tdctimeblk_tw, &b_Ndata_earm_ecal_clus_tdctimeblk_tw);
   fChain->SetBranchAddress("earm.ecal.clus.tdctimeblk_tw", earm_ecal_clus_tdctimeblk_tw, &b_earm_ecal_clus_tdctimeblk_tw);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.x", &Ndata_earm_ecal_clus_x, &b_Ndata_earm_ecal_clus_x);
   fChain->SetBranchAddress("earm.ecal.clus.x", earm_ecal_clus_x, &b_earm_ecal_clus_x);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus.y", &Ndata_earm_ecal_clus_y, &b_Ndata_earm_ecal_clus_y);
   fChain->SetBranchAddress("earm.ecal.clus.y", earm_ecal_clus_y, &b_earm_ecal_clus_y);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus_blk.again", &Ndata_earm_ecal_clus_blk_again, &b_Ndata_earm_ecal_clus_blk_again);
   fChain->SetBranchAddress("earm.ecal.clus_blk.again", earm_ecal_clus_blk_again, &b_earm_ecal_clus_blk_again);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus_blk.atime", &Ndata_earm_ecal_clus_blk_atime, &b_Ndata_earm_ecal_clus_blk_atime);
   fChain->SetBranchAddress("earm.ecal.clus_blk.atime", earm_ecal_clus_blk_atime, &b_earm_ecal_clus_blk_atime);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus_blk.col", &Ndata_earm_ecal_clus_blk_col, &b_Ndata_earm_ecal_clus_blk_col);
   fChain->SetBranchAddress("earm.ecal.clus_blk.col", earm_ecal_clus_blk_col, &b_earm_ecal_clus_blk_col);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus_blk.e", &Ndata_earm_ecal_clus_blk_e, &b_Ndata_earm_ecal_clus_blk_e);
   fChain->SetBranchAddress("earm.ecal.clus_blk.e", earm_ecal_clus_blk_e, &b_earm_ecal_clus_blk_e);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus_blk.id", &Ndata_earm_ecal_clus_blk_id, &b_Ndata_earm_ecal_clus_blk_id);
   fChain->SetBranchAddress("earm.ecal.clus_blk.id", earm_ecal_clus_blk_id, &b_earm_ecal_clus_blk_id);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus_blk.row", &Ndata_earm_ecal_clus_blk_row, &b_Ndata_earm_ecal_clus_blk_row);
   fChain->SetBranchAddress("earm.ecal.clus_blk.row", earm_ecal_clus_blk_row, &b_earm_ecal_clus_blk_row);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus_blk.tdctime", &Ndata_earm_ecal_clus_blk_tdctime, &b_Ndata_earm_ecal_clus_blk_tdctime);
   fChain->SetBranchAddress("earm.ecal.clus_blk.tdctime", earm_ecal_clus_blk_tdctime, &b_earm_ecal_clus_blk_tdctime);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus_blk.tdctime_tw", &Ndata_earm_ecal_clus_blk_tdctime_tw, &b_Ndata_earm_ecal_clus_blk_tdctime_tw);
   fChain->SetBranchAddress("earm.ecal.clus_blk.tdctime_tw", earm_ecal_clus_blk_tdctime_tw, &b_earm_ecal_clus_blk_tdctime_tw);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus_blk.x", &Ndata_earm_ecal_clus_blk_x, &b_Ndata_earm_ecal_clus_blk_x);
   fChain->SetBranchAddress("earm.ecal.clus_blk.x", earm_ecal_clus_blk_x, &b_earm_ecal_clus_blk_x);
   fChain->SetBranchAddress("Ndata.earm.ecal.clus_blk.y", &Ndata_earm_ecal_clus_blk_y, &b_Ndata_earm_ecal_clus_blk_y);
   fChain->SetBranchAddress("earm.ecal.clus_blk.y", earm_ecal_clus_blk_y, &b_earm_ecal_clus_blk_y);
   fChain->SetBranchAddress("Ndata.earm.ecal.ped", &Ndata_earm_ecal_ped, &b_Ndata_earm_ecal_ped);
   fChain->SetBranchAddress("earm.ecal.ped", earm_ecal_ped, &b_earm_ecal_ped);
   fChain->SetBranchAddress("Ndata.earm.tr.beta", &Ndata_earm_tr_beta, &b_Ndata_earm_tr_beta);
   fChain->SetBranchAddress("earm.tr.beta", earm_tr_beta, &b_earm_tr_beta);
   fChain->SetBranchAddress("Ndata.earm.tr.chi2", &Ndata_earm_tr_chi2, &b_Ndata_earm_tr_chi2);
   fChain->SetBranchAddress("earm.tr.chi2", earm_tr_chi2, &b_earm_tr_chi2);
   fChain->SetBranchAddress("Ndata.earm.tr.d_ph", &Ndata_earm_tr_d_ph, &b_Ndata_earm_tr_d_ph);
   fChain->SetBranchAddress("earm.tr.d_ph", earm_tr_d_ph, &b_earm_tr_d_ph);
   fChain->SetBranchAddress("Ndata.earm.tr.d_th", &Ndata_earm_tr_d_th, &b_Ndata_earm_tr_d_th);
   fChain->SetBranchAddress("earm.tr.d_th", earm_tr_d_th, &b_earm_tr_d_th);
   fChain->SetBranchAddress("Ndata.earm.tr.d_x", &Ndata_earm_tr_d_x, &b_Ndata_earm_tr_d_x);
   fChain->SetBranchAddress("earm.tr.d_x", earm_tr_d_x, &b_earm_tr_d_x);
   fChain->SetBranchAddress("Ndata.earm.tr.d_y", &Ndata_earm_tr_d_y, &b_Ndata_earm_tr_d_y);
   fChain->SetBranchAddress("earm.tr.d_y", earm_tr_d_y, &b_earm_tr_d_y);
   fChain->SetBranchAddress("Ndata.earm.tr.dbeta", &Ndata_earm_tr_dbeta, &b_Ndata_earm_tr_dbeta);
   fChain->SetBranchAddress("earm.tr.dbeta", earm_tr_dbeta, &b_earm_tr_dbeta);
   fChain->SetBranchAddress("Ndata.earm.tr.dtime", &Ndata_earm_tr_dtime, &b_Ndata_earm_tr_dtime);
   fChain->SetBranchAddress("earm.tr.dtime", earm_tr_dtime, &b_earm_tr_dtime);
   fChain->SetBranchAddress("Ndata.earm.tr.flag", &Ndata_earm_tr_flag, &b_Ndata_earm_tr_flag);
   fChain->SetBranchAddress("earm.tr.flag", earm_tr_flag, &b_earm_tr_flag);
   fChain->SetBranchAddress("Ndata.earm.tr.ndof", &Ndata_earm_tr_ndof, &b_Ndata_earm_tr_ndof);
   fChain->SetBranchAddress("earm.tr.ndof", earm_tr_ndof, &b_earm_tr_ndof);
   fChain->SetBranchAddress("Ndata.earm.tr.p", &Ndata_earm_tr_p, &b_Ndata_earm_tr_p);
   fChain->SetBranchAddress("earm.tr.p", earm_tr_p, &b_earm_tr_p);
   fChain->SetBranchAddress("Ndata.earm.tr.pathl", &Ndata_earm_tr_pathl, &b_Ndata_earm_tr_pathl);
   fChain->SetBranchAddress("earm.tr.pathl", earm_tr_pathl, &b_earm_tr_pathl);
   fChain->SetBranchAddress("Ndata.earm.tr.ph", &Ndata_earm_tr_ph, &b_Ndata_earm_tr_ph);
   fChain->SetBranchAddress("earm.tr.ph", earm_tr_ph, &b_earm_tr_ph);
   fChain->SetBranchAddress("Ndata.earm.tr.px", &Ndata_earm_tr_px, &b_Ndata_earm_tr_px);
   fChain->SetBranchAddress("earm.tr.px", earm_tr_px, &b_earm_tr_px);
   fChain->SetBranchAddress("Ndata.earm.tr.py", &Ndata_earm_tr_py, &b_Ndata_earm_tr_py);
   fChain->SetBranchAddress("earm.tr.py", earm_tr_py, &b_earm_tr_py);
   fChain->SetBranchAddress("Ndata.earm.tr.pz", &Ndata_earm_tr_pz, &b_Ndata_earm_tr_pz);
   fChain->SetBranchAddress("earm.tr.pz", earm_tr_pz, &b_earm_tr_pz);
   fChain->SetBranchAddress("Ndata.earm.tr.r_ph", &Ndata_earm_tr_r_ph, &b_Ndata_earm_tr_r_ph);
   fChain->SetBranchAddress("earm.tr.r_ph", earm_tr_r_ph, &b_earm_tr_r_ph);
   fChain->SetBranchAddress("Ndata.earm.tr.r_th", &Ndata_earm_tr_r_th, &b_Ndata_earm_tr_r_th);
   fChain->SetBranchAddress("earm.tr.r_th", earm_tr_r_th, &b_earm_tr_r_th);
   fChain->SetBranchAddress("Ndata.earm.tr.r_x", &Ndata_earm_tr_r_x, &b_Ndata_earm_tr_r_x);
   fChain->SetBranchAddress("earm.tr.r_x", earm_tr_r_x, &b_earm_tr_r_x);
   fChain->SetBranchAddress("Ndata.earm.tr.r_y", &Ndata_earm_tr_r_y, &b_Ndata_earm_tr_r_y);
   fChain->SetBranchAddress("earm.tr.r_y", earm_tr_r_y, &b_earm_tr_r_y);
   fChain->SetBranchAddress("Ndata.earm.tr.tg_dp", &Ndata_earm_tr_tg_dp, &b_Ndata_earm_tr_tg_dp);
   fChain->SetBranchAddress("earm.tr.tg_dp", earm_tr_tg_dp, &b_earm_tr_tg_dp);
   fChain->SetBranchAddress("Ndata.earm.tr.tg_ph", &Ndata_earm_tr_tg_ph, &b_Ndata_earm_tr_tg_ph);
   fChain->SetBranchAddress("earm.tr.tg_ph", earm_tr_tg_ph, &b_earm_tr_tg_ph);
   fChain->SetBranchAddress("Ndata.earm.tr.tg_th", &Ndata_earm_tr_tg_th, &b_Ndata_earm_tr_tg_th);
   fChain->SetBranchAddress("earm.tr.tg_th", earm_tr_tg_th, &b_earm_tr_tg_th);
   fChain->SetBranchAddress("Ndata.earm.tr.tg_y", &Ndata_earm_tr_tg_y, &b_Ndata_earm_tr_tg_y);
   fChain->SetBranchAddress("earm.tr.tg_y", earm_tr_tg_y, &b_earm_tr_tg_y);
   fChain->SetBranchAddress("Ndata.earm.tr.th", &Ndata_earm_tr_th, &b_Ndata_earm_tr_th);
   fChain->SetBranchAddress("earm.tr.th", earm_tr_th, &b_earm_tr_th);
   fChain->SetBranchAddress("Ndata.earm.tr.time", &Ndata_earm_tr_time, &b_Ndata_earm_tr_time);
   fChain->SetBranchAddress("earm.tr.time", earm_tr_time, &b_earm_tr_time);
   fChain->SetBranchAddress("Ndata.earm.tr.vx", &Ndata_earm_tr_vx, &b_Ndata_earm_tr_vx);
   fChain->SetBranchAddress("earm.tr.vx", earm_tr_vx, &b_earm_tr_vx);
   fChain->SetBranchAddress("Ndata.earm.tr.vy", &Ndata_earm_tr_vy, &b_Ndata_earm_tr_vy);
   fChain->SetBranchAddress("earm.tr.vy", earm_tr_vy, &b_earm_tr_vy);
   fChain->SetBranchAddress("Ndata.earm.tr.vz", &Ndata_earm_tr_vz, &b_Ndata_earm_tr_vz);
   fChain->SetBranchAddress("earm.tr.vz", earm_tr_vz, &b_earm_tr_vz);
   fChain->SetBranchAddress("Ndata.earm.tr.x", &Ndata_earm_tr_x, &b_Ndata_earm_tr_x);
   fChain->SetBranchAddress("earm.tr.x", earm_tr_x, &b_earm_tr_x);
   fChain->SetBranchAddress("Ndata.earm.tr.y", &Ndata_earm_tr_y, &b_Ndata_earm_tr_y);
   fChain->SetBranchAddress("earm.tr.y", earm_tr_y, &b_earm_tr_y);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCU", &Ndata_sbs_gemFPP_hit_ADCU, &b_Ndata_sbs_gemFPP_hit_ADCU);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCU", sbs_gemFPP_hit_ADCU, &b_sbs_gemFPP_hit_ADCU);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCU_deconv", &Ndata_sbs_gemFPP_hit_ADCU_deconv, &b_Ndata_sbs_gemFPP_hit_ADCU_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCU_deconv", sbs_gemFPP_hit_ADCU_deconv, &b_sbs_gemFPP_hit_ADCU_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCV", &Ndata_sbs_gemFPP_hit_ADCV, &b_Ndata_sbs_gemFPP_hit_ADCV);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCV", sbs_gemFPP_hit_ADCV, &b_sbs_gemFPP_hit_ADCV);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCV_deconv", &Ndata_sbs_gemFPP_hit_ADCV_deconv, &b_Ndata_sbs_gemFPP_hit_ADCV_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCV_deconv", sbs_gemFPP_hit_ADCV_deconv, &b_sbs_gemFPP_hit_ADCV_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCasym", &Ndata_sbs_gemFPP_hit_ADCasym, &b_Ndata_sbs_gemFPP_hit_ADCasym);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCasym", sbs_gemFPP_hit_ADCasym, &b_sbs_gemFPP_hit_ADCasym);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCasym_deconv", &Ndata_sbs_gemFPP_hit_ADCasym_deconv, &b_Ndata_sbs_gemFPP_hit_ADCasym_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCasym_deconv", sbs_gemFPP_hit_ADCasym_deconv, &b_sbs_gemFPP_hit_ADCasym_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCavg", &Ndata_sbs_gemFPP_hit_ADCavg, &b_Ndata_sbs_gemFPP_hit_ADCavg);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCavg", sbs_gemFPP_hit_ADCavg, &b_sbs_gemFPP_hit_ADCavg);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCavg_deconv", &Ndata_sbs_gemFPP_hit_ADCavg_deconv, &b_Ndata_sbs_gemFPP_hit_ADCavg_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCavg_deconv", sbs_gemFPP_hit_ADCavg_deconv, &b_sbs_gemFPP_hit_ADCavg_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac0_Umax", &Ndata_sbs_gemFPP_hit_ADCfrac0_Umax, &b_Ndata_sbs_gemFPP_hit_ADCfrac0_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac0_Umax", sbs_gemFPP_hit_ADCfrac0_Umax, &b_sbs_gemFPP_hit_ADCfrac0_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac0_Vmax", &Ndata_sbs_gemFPP_hit_ADCfrac0_Vmax, &b_Ndata_sbs_gemFPP_hit_ADCfrac0_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac0_Vmax", sbs_gemFPP_hit_ADCfrac0_Vmax, &b_sbs_gemFPP_hit_ADCfrac0_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac1_Umax", &Ndata_sbs_gemFPP_hit_ADCfrac1_Umax, &b_Ndata_sbs_gemFPP_hit_ADCfrac1_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac1_Umax", sbs_gemFPP_hit_ADCfrac1_Umax, &b_sbs_gemFPP_hit_ADCfrac1_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac1_Vmax", &Ndata_sbs_gemFPP_hit_ADCfrac1_Vmax, &b_Ndata_sbs_gemFPP_hit_ADCfrac1_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac1_Vmax", sbs_gemFPP_hit_ADCfrac1_Vmax, &b_sbs_gemFPP_hit_ADCfrac1_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac2_Umax", &Ndata_sbs_gemFPP_hit_ADCfrac2_Umax, &b_Ndata_sbs_gemFPP_hit_ADCfrac2_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac2_Umax", sbs_gemFPP_hit_ADCfrac2_Umax, &b_sbs_gemFPP_hit_ADCfrac2_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac2_Vmax", &Ndata_sbs_gemFPP_hit_ADCfrac2_Vmax, &b_Ndata_sbs_gemFPP_hit_ADCfrac2_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac2_Vmax", sbs_gemFPP_hit_ADCfrac2_Vmax, &b_sbs_gemFPP_hit_ADCfrac2_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac3_Umax", &Ndata_sbs_gemFPP_hit_ADCfrac3_Umax, &b_Ndata_sbs_gemFPP_hit_ADCfrac3_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac3_Umax", sbs_gemFPP_hit_ADCfrac3_Umax, &b_sbs_gemFPP_hit_ADCfrac3_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac3_Vmax", &Ndata_sbs_gemFPP_hit_ADCfrac3_Vmax, &b_Ndata_sbs_gemFPP_hit_ADCfrac3_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac3_Vmax", sbs_gemFPP_hit_ADCfrac3_Vmax, &b_sbs_gemFPP_hit_ADCfrac3_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac4_Umax", &Ndata_sbs_gemFPP_hit_ADCfrac4_Umax, &b_Ndata_sbs_gemFPP_hit_ADCfrac4_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac4_Umax", sbs_gemFPP_hit_ADCfrac4_Umax, &b_sbs_gemFPP_hit_ADCfrac4_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac4_Vmax", &Ndata_sbs_gemFPP_hit_ADCfrac4_Vmax, &b_Ndata_sbs_gemFPP_hit_ADCfrac4_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac4_Vmax", sbs_gemFPP_hit_ADCfrac4_Vmax, &b_sbs_gemFPP_hit_ADCfrac4_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac5_Umax", &Ndata_sbs_gemFPP_hit_ADCfrac5_Umax, &b_Ndata_sbs_gemFPP_hit_ADCfrac5_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac5_Umax", sbs_gemFPP_hit_ADCfrac5_Umax, &b_sbs_gemFPP_hit_ADCfrac5_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCfrac5_Vmax", &Ndata_sbs_gemFPP_hit_ADCfrac5_Vmax, &b_Ndata_sbs_gemFPP_hit_ADCfrac5_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCfrac5_Vmax", sbs_gemFPP_hit_ADCfrac5_Vmax, &b_sbs_gemFPP_hit_ADCfrac5_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxcomboUclust_deconv", &Ndata_sbs_gemFPP_hit_ADCmaxcomboUclust_deconv, &b_Ndata_sbs_gemFPP_hit_ADCmaxcomboUclust_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxcomboUclust_deconv", sbs_gemFPP_hit_ADCmaxcomboUclust_deconv, &b_sbs_gemFPP_hit_ADCmaxcomboUclust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxcomboVclust_deconv", &Ndata_sbs_gemFPP_hit_ADCmaxcomboVclust_deconv, &b_Ndata_sbs_gemFPP_hit_ADCmaxcomboVclust_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxcomboVclust_deconv", sbs_gemFPP_hit_ADCmaxcomboVclust_deconv, &b_sbs_gemFPP_hit_ADCmaxcomboVclust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxsampU", &Ndata_sbs_gemFPP_hit_ADCmaxsampU, &b_Ndata_sbs_gemFPP_hit_ADCmaxsampU);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxsampU", sbs_gemFPP_hit_ADCmaxsampU, &b_sbs_gemFPP_hit_ADCmaxsampU);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxsampUclust", &Ndata_sbs_gemFPP_hit_ADCmaxsampUclust, &b_Ndata_sbs_gemFPP_hit_ADCmaxsampUclust);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxsampUclust", sbs_gemFPP_hit_ADCmaxsampUclust, &b_sbs_gemFPP_hit_ADCmaxsampUclust);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxsampUclust_deconv", &Ndata_sbs_gemFPP_hit_ADCmaxsampUclust_deconv, &b_Ndata_sbs_gemFPP_hit_ADCmaxsampUclust_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxsampUclust_deconv", sbs_gemFPP_hit_ADCmaxsampUclust_deconv, &b_sbs_gemFPP_hit_ADCmaxsampUclust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxsampV", &Ndata_sbs_gemFPP_hit_ADCmaxsampV, &b_Ndata_sbs_gemFPP_hit_ADCmaxsampV);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxsampV", sbs_gemFPP_hit_ADCmaxsampV, &b_sbs_gemFPP_hit_ADCmaxsampV);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxsampVclust", &Ndata_sbs_gemFPP_hit_ADCmaxsampVclust, &b_Ndata_sbs_gemFPP_hit_ADCmaxsampVclust);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxsampVclust", sbs_gemFPP_hit_ADCmaxsampVclust, &b_sbs_gemFPP_hit_ADCmaxsampVclust);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxsampVclust_deconv", &Ndata_sbs_gemFPP_hit_ADCmaxsampVclust_deconv, &b_Ndata_sbs_gemFPP_hit_ADCmaxsampVclust_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxsampVclust_deconv", sbs_gemFPP_hit_ADCmaxsampVclust_deconv, &b_sbs_gemFPP_hit_ADCmaxsampVclust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxstripU", &Ndata_sbs_gemFPP_hit_ADCmaxstripU, &b_Ndata_sbs_gemFPP_hit_ADCmaxstripU);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxstripU", sbs_gemFPP_hit_ADCmaxstripU, &b_sbs_gemFPP_hit_ADCmaxstripU);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ADCmaxstripV", &Ndata_sbs_gemFPP_hit_ADCmaxstripV, &b_Ndata_sbs_gemFPP_hit_ADCmaxstripV);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ADCmaxstripV", sbs_gemFPP_hit_ADCmaxstripV, &b_sbs_gemFPP_hit_ADCmaxstripV);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.BUILD_ALL_SAMPLES_U", &Ndata_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_U, &b_Ndata_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_U);
   fChain->SetBranchAddress("sbs.gemFPP.hit.BUILD_ALL_SAMPLES_U", sbs_gemFPP_hit_BUILD_ALL_SAMPLES_U, &b_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.BUILD_ALL_SAMPLES_V", &Ndata_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_V, &b_Ndata_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_V);
   fChain->SetBranchAddress("sbs.gemFPP.hit.BUILD_ALL_SAMPLES_V", sbs_gemFPP_hit_BUILD_ALL_SAMPLES_V, &b_sbs_gemFPP_hit_BUILD_ALL_SAMPLES_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.CM_GOOD_U", &Ndata_sbs_gemFPP_hit_CM_GOOD_U, &b_Ndata_sbs_gemFPP_hit_CM_GOOD_U);
   fChain->SetBranchAddress("sbs.gemFPP.hit.CM_GOOD_U", sbs_gemFPP_hit_CM_GOOD_U, &b_sbs_gemFPP_hit_CM_GOOD_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.CM_GOOD_V", &Ndata_sbs_gemFPP_hit_CM_GOOD_V, &b_Ndata_sbs_gemFPP_hit_CM_GOOD_V);
   fChain->SetBranchAddress("sbs.gemFPP.hit.CM_GOOD_V", sbs_gemFPP_hit_CM_GOOD_V, &b_sbs_gemFPP_hit_CM_GOOD_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC0_Umax", &Ndata_sbs_gemFPP_hit_DeconvADC0_Umax, &b_Ndata_sbs_gemFPP_hit_DeconvADC0_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC0_Umax", sbs_gemFPP_hit_DeconvADC0_Umax, &b_sbs_gemFPP_hit_DeconvADC0_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC0_Vmax", &Ndata_sbs_gemFPP_hit_DeconvADC0_Vmax, &b_Ndata_sbs_gemFPP_hit_DeconvADC0_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC0_Vmax", sbs_gemFPP_hit_DeconvADC0_Vmax, &b_sbs_gemFPP_hit_DeconvADC0_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC1_Umax", &Ndata_sbs_gemFPP_hit_DeconvADC1_Umax, &b_Ndata_sbs_gemFPP_hit_DeconvADC1_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC1_Umax", sbs_gemFPP_hit_DeconvADC1_Umax, &b_sbs_gemFPP_hit_DeconvADC1_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC1_Vmax", &Ndata_sbs_gemFPP_hit_DeconvADC1_Vmax, &b_Ndata_sbs_gemFPP_hit_DeconvADC1_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC1_Vmax", sbs_gemFPP_hit_DeconvADC1_Vmax, &b_sbs_gemFPP_hit_DeconvADC1_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC2_Umax", &Ndata_sbs_gemFPP_hit_DeconvADC2_Umax, &b_Ndata_sbs_gemFPP_hit_DeconvADC2_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC2_Umax", sbs_gemFPP_hit_DeconvADC2_Umax, &b_sbs_gemFPP_hit_DeconvADC2_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC2_Vmax", &Ndata_sbs_gemFPP_hit_DeconvADC2_Vmax, &b_Ndata_sbs_gemFPP_hit_DeconvADC2_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC2_Vmax", sbs_gemFPP_hit_DeconvADC2_Vmax, &b_sbs_gemFPP_hit_DeconvADC2_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC3_Umax", &Ndata_sbs_gemFPP_hit_DeconvADC3_Umax, &b_Ndata_sbs_gemFPP_hit_DeconvADC3_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC3_Umax", sbs_gemFPP_hit_DeconvADC3_Umax, &b_sbs_gemFPP_hit_DeconvADC3_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC3_Vmax", &Ndata_sbs_gemFPP_hit_DeconvADC3_Vmax, &b_Ndata_sbs_gemFPP_hit_DeconvADC3_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC3_Vmax", sbs_gemFPP_hit_DeconvADC3_Vmax, &b_sbs_gemFPP_hit_DeconvADC3_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC4_Umax", &Ndata_sbs_gemFPP_hit_DeconvADC4_Umax, &b_Ndata_sbs_gemFPP_hit_DeconvADC4_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC4_Umax", sbs_gemFPP_hit_DeconvADC4_Umax, &b_sbs_gemFPP_hit_DeconvADC4_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC4_Vmax", &Ndata_sbs_gemFPP_hit_DeconvADC4_Vmax, &b_Ndata_sbs_gemFPP_hit_DeconvADC4_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC4_Vmax", sbs_gemFPP_hit_DeconvADC4_Vmax, &b_sbs_gemFPP_hit_DeconvADC4_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC5_Umax", &Ndata_sbs_gemFPP_hit_DeconvADC5_Umax, &b_Ndata_sbs_gemFPP_hit_DeconvADC5_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC5_Umax", sbs_gemFPP_hit_DeconvADC5_Umax, &b_sbs_gemFPP_hit_DeconvADC5_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADC5_Vmax", &Ndata_sbs_gemFPP_hit_DeconvADC5_Vmax, &b_Ndata_sbs_gemFPP_hit_DeconvADC5_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADC5_Vmax", sbs_gemFPP_hit_DeconvADC5_Vmax, &b_sbs_gemFPP_hit_DeconvADC5_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADCmaxcomboU", &Ndata_sbs_gemFPP_hit_DeconvADCmaxcomboU, &b_Ndata_sbs_gemFPP_hit_DeconvADCmaxcomboU);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADCmaxcomboU", sbs_gemFPP_hit_DeconvADCmaxcomboU, &b_sbs_gemFPP_hit_DeconvADCmaxcomboU);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADCmaxcomboV", &Ndata_sbs_gemFPP_hit_DeconvADCmaxcomboV, &b_Ndata_sbs_gemFPP_hit_DeconvADCmaxcomboV);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADCmaxcomboV", sbs_gemFPP_hit_DeconvADCmaxcomboV, &b_sbs_gemFPP_hit_DeconvADCmaxcomboV);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADCmaxsampU", &Ndata_sbs_gemFPP_hit_DeconvADCmaxsampU, &b_Ndata_sbs_gemFPP_hit_DeconvADCmaxsampU);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADCmaxsampU", sbs_gemFPP_hit_DeconvADCmaxsampU, &b_sbs_gemFPP_hit_DeconvADCmaxsampU);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADCmaxsampV", &Ndata_sbs_gemFPP_hit_DeconvADCmaxsampV, &b_Ndata_sbs_gemFPP_hit_DeconvADCmaxsampV);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADCmaxsampV", sbs_gemFPP_hit_DeconvADCmaxsampV, &b_sbs_gemFPP_hit_DeconvADCmaxsampV);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADCmaxstripU", &Ndata_sbs_gemFPP_hit_DeconvADCmaxstripU, &b_Ndata_sbs_gemFPP_hit_DeconvADCmaxstripU);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADCmaxstripU", sbs_gemFPP_hit_DeconvADCmaxstripU, &b_sbs_gemFPP_hit_DeconvADCmaxstripU);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.DeconvADCmaxstripV", &Ndata_sbs_gemFPP_hit_DeconvADCmaxstripV, &b_Ndata_sbs_gemFPP_hit_DeconvADCmaxstripV);
   fChain->SetBranchAddress("sbs.gemFPP.hit.DeconvADCmaxstripV", sbs_gemFPP_hit_DeconvADCmaxstripV, &b_sbs_gemFPP_hit_DeconvADCmaxstripV);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ENABLE_CM_U", &Ndata_sbs_gemFPP_hit_ENABLE_CM_U, &b_Ndata_sbs_gemFPP_hit_ENABLE_CM_U);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ENABLE_CM_U", sbs_gemFPP_hit_ENABLE_CM_U, &b_sbs_gemFPP_hit_ENABLE_CM_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ENABLE_CM_V", &Ndata_sbs_gemFPP_hit_ENABLE_CM_V, &b_Ndata_sbs_gemFPP_hit_ENABLE_CM_V);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ENABLE_CM_V", sbs_gemFPP_hit_ENABLE_CM_V, &b_sbs_gemFPP_hit_ENABLE_CM_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.TSchi2_Umax", &Ndata_sbs_gemFPP_hit_TSchi2_Umax, &b_Ndata_sbs_gemFPP_hit_TSchi2_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.TSchi2_Umax", sbs_gemFPP_hit_TSchi2_Umax, &b_sbs_gemFPP_hit_TSchi2_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.TSchi2_Vmax", &Ndata_sbs_gemFPP_hit_TSchi2_Vmax, &b_Ndata_sbs_gemFPP_hit_TSchi2_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.TSchi2_Vmax", sbs_gemFPP_hit_TSchi2_Vmax, &b_sbs_gemFPP_hit_TSchi2_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.TSprob_Umax", &Ndata_sbs_gemFPP_hit_TSprob_Umax, &b_Ndata_sbs_gemFPP_hit_TSprob_Umax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.TSprob_Umax", sbs_gemFPP_hit_TSprob_Umax, &b_sbs_gemFPP_hit_TSprob_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.TSprob_Vmax", &Ndata_sbs_gemFPP_hit_TSprob_Vmax, &b_Ndata_sbs_gemFPP_hit_TSprob_Vmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.TSprob_Vmax", sbs_gemFPP_hit_TSprob_Vmax, &b_sbs_gemFPP_hit_TSprob_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.Tavg", &Ndata_sbs_gemFPP_hit_Tavg, &b_Ndata_sbs_gemFPP_hit_Tavg);
   fChain->SetBranchAddress("sbs.gemFPP.hit.Tavg", sbs_gemFPP_hit_Tavg, &b_sbs_gemFPP_hit_Tavg);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.Tavg_corr", &Ndata_sbs_gemFPP_hit_Tavg_corr, &b_Ndata_sbs_gemFPP_hit_Tavg_corr);
   fChain->SetBranchAddress("sbs.gemFPP.hit.Tavg_corr", sbs_gemFPP_hit_Tavg_corr, &b_sbs_gemFPP_hit_Tavg_corr);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.Tavg_deconv", &Ndata_sbs_gemFPP_hit_Tavg_deconv, &b_Ndata_sbs_gemFPP_hit_Tavg_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.Tavg_deconv", sbs_gemFPP_hit_Tavg_deconv, &b_sbs_gemFPP_hit_Tavg_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.Tavg_fit", &Ndata_sbs_gemFPP_hit_Tavg_fit, &b_Ndata_sbs_gemFPP_hit_Tavg_fit);
   fChain->SetBranchAddress("sbs.gemFPP.hit.Tavg_fit", sbs_gemFPP_hit_Tavg_fit, &b_sbs_gemFPP_hit_Tavg_fit);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.Ugain", &Ndata_sbs_gemFPP_hit_Ugain, &b_Ndata_sbs_gemFPP_hit_Ugain);
   fChain->SetBranchAddress("sbs.gemFPP.hit.Ugain", sbs_gemFPP_hit_Ugain, &b_sbs_gemFPP_hit_Ugain);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.Utime", &Ndata_sbs_gemFPP_hit_Utime, &b_Ndata_sbs_gemFPP_hit_Utime);
   fChain->SetBranchAddress("sbs.gemFPP.hit.Utime", sbs_gemFPP_hit_Utime, &b_sbs_gemFPP_hit_Utime);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.UtimeDeconv", &Ndata_sbs_gemFPP_hit_UtimeDeconv, &b_Ndata_sbs_gemFPP_hit_UtimeDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.UtimeDeconv", sbs_gemFPP_hit_UtimeDeconv, &b_sbs_gemFPP_hit_UtimeDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.UtimeFit", &Ndata_sbs_gemFPP_hit_UtimeFit, &b_Ndata_sbs_gemFPP_hit_UtimeFit);
   fChain->SetBranchAddress("sbs.gemFPP.hit.UtimeFit", sbs_gemFPP_hit_UtimeFit, &b_sbs_gemFPP_hit_UtimeFit);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.UtimeMaxStrip", &Ndata_sbs_gemFPP_hit_UtimeMaxStrip, &b_Ndata_sbs_gemFPP_hit_UtimeMaxStrip);
   fChain->SetBranchAddress("sbs.gemFPP.hit.UtimeMaxStrip", sbs_gemFPP_hit_UtimeMaxStrip, &b_sbs_gemFPP_hit_UtimeMaxStrip);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.UtimeMaxStripDeconv", &Ndata_sbs_gemFPP_hit_UtimeMaxStripDeconv, &b_Ndata_sbs_gemFPP_hit_UtimeMaxStripDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.UtimeMaxStripDeconv", sbs_gemFPP_hit_UtimeMaxStripDeconv, &b_sbs_gemFPP_hit_UtimeMaxStripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.UtimeMaxStripFit", &Ndata_sbs_gemFPP_hit_UtimeMaxStripFit, &b_Ndata_sbs_gemFPP_hit_UtimeMaxStripFit);
   fChain->SetBranchAddress("sbs.gemFPP.hit.UtimeMaxStripFit", sbs_gemFPP_hit_UtimeMaxStripFit, &b_sbs_gemFPP_hit_UtimeMaxStripFit);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.Vgain", &Ndata_sbs_gemFPP_hit_Vgain, &b_Ndata_sbs_gemFPP_hit_Vgain);
   fChain->SetBranchAddress("sbs.gemFPP.hit.Vgain", sbs_gemFPP_hit_Vgain, &b_sbs_gemFPP_hit_Vgain);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.Vtime", &Ndata_sbs_gemFPP_hit_Vtime, &b_Ndata_sbs_gemFPP_hit_Vtime);
   fChain->SetBranchAddress("sbs.gemFPP.hit.Vtime", sbs_gemFPP_hit_Vtime, &b_sbs_gemFPP_hit_Vtime);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.VtimeDeconv", &Ndata_sbs_gemFPP_hit_VtimeDeconv, &b_Ndata_sbs_gemFPP_hit_VtimeDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.VtimeDeconv", sbs_gemFPP_hit_VtimeDeconv, &b_sbs_gemFPP_hit_VtimeDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.VtimeFit", &Ndata_sbs_gemFPP_hit_VtimeFit, &b_Ndata_sbs_gemFPP_hit_VtimeFit);
   fChain->SetBranchAddress("sbs.gemFPP.hit.VtimeFit", sbs_gemFPP_hit_VtimeFit, &b_sbs_gemFPP_hit_VtimeFit);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.VtimeMaxStrip", &Ndata_sbs_gemFPP_hit_VtimeMaxStrip, &b_Ndata_sbs_gemFPP_hit_VtimeMaxStrip);
   fChain->SetBranchAddress("sbs.gemFPP.hit.VtimeMaxStrip", sbs_gemFPP_hit_VtimeMaxStrip, &b_sbs_gemFPP_hit_VtimeMaxStrip);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.VtimeMaxStripDeconv", &Ndata_sbs_gemFPP_hit_VtimeMaxStripDeconv, &b_Ndata_sbs_gemFPP_hit_VtimeMaxStripDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.VtimeMaxStripDeconv", sbs_gemFPP_hit_VtimeMaxStripDeconv, &b_sbs_gemFPP_hit_VtimeMaxStripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.VtimeMaxStripFit", &Ndata_sbs_gemFPP_hit_VtimeMaxStripFit, &b_Ndata_sbs_gemFPP_hit_VtimeMaxStripFit);
   fChain->SetBranchAddress("sbs.gemFPP.hit.VtimeMaxStripFit", sbs_gemFPP_hit_VtimeMaxStripFit, &b_sbs_gemFPP_hit_VtimeMaxStripFit);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.adc_id_U", &Ndata_sbs_gemFPP_hit_adc_id_U, &b_Ndata_sbs_gemFPP_hit_adc_id_U);
   fChain->SetBranchAddress("sbs.gemFPP.hit.adc_id_U", sbs_gemFPP_hit_adc_id_U, &b_sbs_gemFPP_hit_adc_id_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.adc_id_V", &Ndata_sbs_gemFPP_hit_adc_id_V, &b_Ndata_sbs_gemFPP_hit_adc_id_V);
   fChain->SetBranchAddress("sbs.gemFPP.hit.adc_id_V", sbs_gemFPP_hit_adc_id_V, &b_sbs_gemFPP_hit_adc_id_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ccor_clust", &Ndata_sbs_gemFPP_hit_ccor_clust, &b_Ndata_sbs_gemFPP_hit_ccor_clust);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ccor_clust", sbs_gemFPP_hit_ccor_clust, &b_sbs_gemFPP_hit_ccor_clust);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ccor_clust_deconv", &Ndata_sbs_gemFPP_hit_ccor_clust_deconv, &b_Ndata_sbs_gemFPP_hit_ccor_clust_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ccor_clust_deconv", sbs_gemFPP_hit_ccor_clust_deconv, &b_sbs_gemFPP_hit_ccor_clust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ccor_strip", &Ndata_sbs_gemFPP_hit_ccor_strip, &b_Ndata_sbs_gemFPP_hit_ccor_strip);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ccor_strip", sbs_gemFPP_hit_ccor_strip, &b_sbs_gemFPP_hit_ccor_strip);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ccor_strip_deconv", &Ndata_sbs_gemFPP_hit_ccor_strip_deconv, &b_Ndata_sbs_gemFPP_hit_ccor_strip_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ccor_strip_deconv", sbs_gemFPP_hit_ccor_strip_deconv, &b_sbs_gemFPP_hit_ccor_strip_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.crate_U", &Ndata_sbs_gemFPP_hit_crate_U, &b_Ndata_sbs_gemFPP_hit_crate_U);
   fChain->SetBranchAddress("sbs.gemFPP.hit.crate_U", sbs_gemFPP_hit_crate_U, &b_sbs_gemFPP_hit_crate_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.crate_V", &Ndata_sbs_gemFPP_hit_crate_V, &b_Ndata_sbs_gemFPP_hit_crate_V);
   fChain->SetBranchAddress("sbs.gemFPP.hit.crate_V", sbs_gemFPP_hit_crate_V, &b_sbs_gemFPP_hit_crate_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.deltat", &Ndata_sbs_gemFPP_hit_deltat, &b_Ndata_sbs_gemFPP_hit_deltat);
   fChain->SetBranchAddress("sbs.gemFPP.hit.deltat", sbs_gemFPP_hit_deltat, &b_sbs_gemFPP_hit_deltat);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.deltat_deconv", &Ndata_sbs_gemFPP_hit_deltat_deconv, &b_Ndata_sbs_gemFPP_hit_deltat_deconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.deltat_deconv", sbs_gemFPP_hit_deltat_deconv, &b_sbs_gemFPP_hit_deltat_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.deltat_fit", &Ndata_sbs_gemFPP_hit_deltat_fit, &b_Ndata_sbs_gemFPP_hit_deltat_fit);
   fChain->SetBranchAddress("sbs.gemFPP.hit.deltat_fit", sbs_gemFPP_hit_deltat_fit, &b_sbs_gemFPP_hit_deltat_fit);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.eresidu", &Ndata_sbs_gemFPP_hit_eresidu, &b_Ndata_sbs_gemFPP_hit_eresidu);
   fChain->SetBranchAddress("sbs.gemFPP.hit.eresidu", sbs_gemFPP_hit_eresidu, &b_sbs_gemFPP_hit_eresidu);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.eresidv", &Ndata_sbs_gemFPP_hit_eresidv, &b_Ndata_sbs_gemFPP_hit_eresidv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.eresidv", sbs_gemFPP_hit_eresidv, &b_sbs_gemFPP_hit_eresidv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.icombomaxUclustDeconv", &Ndata_sbs_gemFPP_hit_icombomaxUclustDeconv, &b_Ndata_sbs_gemFPP_hit_icombomaxUclustDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.icombomaxUclustDeconv", sbs_gemFPP_hit_icombomaxUclustDeconv, &b_sbs_gemFPP_hit_icombomaxUclustDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.icombomaxUstripDeconv", &Ndata_sbs_gemFPP_hit_icombomaxUstripDeconv, &b_Ndata_sbs_gemFPP_hit_icombomaxUstripDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.icombomaxUstripDeconv", sbs_gemFPP_hit_icombomaxUstripDeconv, &b_sbs_gemFPP_hit_icombomaxUstripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.icombomaxVclustDeconv", &Ndata_sbs_gemFPP_hit_icombomaxVclustDeconv, &b_Ndata_sbs_gemFPP_hit_icombomaxVclustDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.icombomaxVclustDeconv", sbs_gemFPP_hit_icombomaxVclustDeconv, &b_sbs_gemFPP_hit_icombomaxVclustDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.icombomaxVstripDeconv", &Ndata_sbs_gemFPP_hit_icombomaxVstripDeconv, &b_Ndata_sbs_gemFPP_hit_icombomaxVstripDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.icombomaxVstripDeconv", sbs_gemFPP_hit_icombomaxVstripDeconv, &b_sbs_gemFPP_hit_icombomaxVstripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.isampmaxUclust", &Ndata_sbs_gemFPP_hit_isampmaxUclust, &b_Ndata_sbs_gemFPP_hit_isampmaxUclust);
   fChain->SetBranchAddress("sbs.gemFPP.hit.isampmaxUclust", sbs_gemFPP_hit_isampmaxUclust, &b_sbs_gemFPP_hit_isampmaxUclust);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.isampmaxUclustDeconv", &Ndata_sbs_gemFPP_hit_isampmaxUclustDeconv, &b_Ndata_sbs_gemFPP_hit_isampmaxUclustDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.isampmaxUclustDeconv", sbs_gemFPP_hit_isampmaxUclustDeconv, &b_sbs_gemFPP_hit_isampmaxUclustDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.isampmaxUstrip", &Ndata_sbs_gemFPP_hit_isampmaxUstrip, &b_Ndata_sbs_gemFPP_hit_isampmaxUstrip);
   fChain->SetBranchAddress("sbs.gemFPP.hit.isampmaxUstrip", sbs_gemFPP_hit_isampmaxUstrip, &b_sbs_gemFPP_hit_isampmaxUstrip);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.isampmaxUstripDeconv", &Ndata_sbs_gemFPP_hit_isampmaxUstripDeconv, &b_Ndata_sbs_gemFPP_hit_isampmaxUstripDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.isampmaxUstripDeconv", sbs_gemFPP_hit_isampmaxUstripDeconv, &b_sbs_gemFPP_hit_isampmaxUstripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.isampmaxVclust", &Ndata_sbs_gemFPP_hit_isampmaxVclust, &b_Ndata_sbs_gemFPP_hit_isampmaxVclust);
   fChain->SetBranchAddress("sbs.gemFPP.hit.isampmaxVclust", sbs_gemFPP_hit_isampmaxVclust, &b_sbs_gemFPP_hit_isampmaxVclust);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.isampmaxVclustDeconv", &Ndata_sbs_gemFPP_hit_isampmaxVclustDeconv, &b_Ndata_sbs_gemFPP_hit_isampmaxVclustDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.isampmaxVclustDeconv", sbs_gemFPP_hit_isampmaxVclustDeconv, &b_sbs_gemFPP_hit_isampmaxVclustDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.isampmaxVstrip", &Ndata_sbs_gemFPP_hit_isampmaxVstrip, &b_Ndata_sbs_gemFPP_hit_isampmaxVstrip);
   fChain->SetBranchAddress("sbs.gemFPP.hit.isampmaxVstrip", sbs_gemFPP_hit_isampmaxVstrip, &b_sbs_gemFPP_hit_isampmaxVstrip);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.isampmaxVstripDeconv", &Ndata_sbs_gemFPP_hit_isampmaxVstripDeconv, &b_Ndata_sbs_gemFPP_hit_isampmaxVstripDeconv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.isampmaxVstripDeconv", sbs_gemFPP_hit_isampmaxVstripDeconv, &b_sbs_gemFPP_hit_isampmaxVstripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.layer", &Ndata_sbs_gemFPP_hit_layer, &b_Ndata_sbs_gemFPP_hit_layer);
   fChain->SetBranchAddress("sbs.gemFPP.hit.layer", sbs_gemFPP_hit_layer, &b_sbs_gemFPP_hit_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.module", &Ndata_sbs_gemFPP_hit_module, &b_Ndata_sbs_gemFPP_hit_module);
   fChain->SetBranchAddress("sbs.gemFPP.hit.module", sbs_gemFPP_hit_module, &b_sbs_gemFPP_hit_module);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.mpd_U", &Ndata_sbs_gemFPP_hit_mpd_U, &b_Ndata_sbs_gemFPP_hit_mpd_U);
   fChain->SetBranchAddress("sbs.gemFPP.hit.mpd_U", sbs_gemFPP_hit_mpd_U, &b_sbs_gemFPP_hit_mpd_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.mpd_V", &Ndata_sbs_gemFPP_hit_mpd_V, &b_Ndata_sbs_gemFPP_hit_mpd_V);
   fChain->SetBranchAddress("sbs.gemFPP.hit.mpd_V", sbs_gemFPP_hit_mpd_V, &b_sbs_gemFPP_hit_mpd_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.nstripu", &Ndata_sbs_gemFPP_hit_nstripu, &b_Ndata_sbs_gemFPP_hit_nstripu);
   fChain->SetBranchAddress("sbs.gemFPP.hit.nstripu", sbs_gemFPP_hit_nstripu, &b_sbs_gemFPP_hit_nstripu);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.nstripv", &Ndata_sbs_gemFPP_hit_nstripv, &b_Ndata_sbs_gemFPP_hit_nstripv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.nstripv", sbs_gemFPP_hit_nstripv, &b_sbs_gemFPP_hit_nstripv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.residu", &Ndata_sbs_gemFPP_hit_residu, &b_Ndata_sbs_gemFPP_hit_residu);
   fChain->SetBranchAddress("sbs.gemFPP.hit.residu", sbs_gemFPP_hit_residu, &b_sbs_gemFPP_hit_residu);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.residv", &Ndata_sbs_gemFPP_hit_residv, &b_Ndata_sbs_gemFPP_hit_residv);
   fChain->SetBranchAddress("sbs.gemFPP.hit.residv", sbs_gemFPP_hit_residv, &b_sbs_gemFPP_hit_residv);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.trackindex", &Ndata_sbs_gemFPP_hit_trackindex, &b_Ndata_sbs_gemFPP_hit_trackindex);
   fChain->SetBranchAddress("sbs.gemFPP.hit.trackindex", sbs_gemFPP_hit_trackindex, &b_sbs_gemFPP_hit_trackindex);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.u", &Ndata_sbs_gemFPP_hit_u, &b_Ndata_sbs_gemFPP_hit_u);
   fChain->SetBranchAddress("sbs.gemFPP.hit.u", sbs_gemFPP_hit_u, &b_sbs_gemFPP_hit_u);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.umoment", &Ndata_sbs_gemFPP_hit_umoment, &b_Ndata_sbs_gemFPP_hit_umoment);
   fChain->SetBranchAddress("sbs.gemFPP.hit.umoment", sbs_gemFPP_hit_umoment, &b_sbs_gemFPP_hit_umoment);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.usigma", &Ndata_sbs_gemFPP_hit_usigma, &b_Ndata_sbs_gemFPP_hit_usigma);
   fChain->SetBranchAddress("sbs.gemFPP.hit.usigma", sbs_gemFPP_hit_usigma, &b_sbs_gemFPP_hit_usigma);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ustriphi", &Ndata_sbs_gemFPP_hit_ustriphi, &b_Ndata_sbs_gemFPP_hit_ustriphi);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ustriphi", sbs_gemFPP_hit_ustriphi, &b_sbs_gemFPP_hit_ustriphi);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ustriplo", &Ndata_sbs_gemFPP_hit_ustriplo, &b_Ndata_sbs_gemFPP_hit_ustriplo);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ustriplo", sbs_gemFPP_hit_ustriplo, &b_sbs_gemFPP_hit_ustriplo);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ustripmax", &Ndata_sbs_gemFPP_hit_ustripmax, &b_Ndata_sbs_gemFPP_hit_ustripmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ustripmax", sbs_gemFPP_hit_ustripmax, &b_sbs_gemFPP_hit_ustripmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.v", &Ndata_sbs_gemFPP_hit_v, &b_Ndata_sbs_gemFPP_hit_v);
   fChain->SetBranchAddress("sbs.gemFPP.hit.v", sbs_gemFPP_hit_v, &b_sbs_gemFPP_hit_v);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.vmoment", &Ndata_sbs_gemFPP_hit_vmoment, &b_Ndata_sbs_gemFPP_hit_vmoment);
   fChain->SetBranchAddress("sbs.gemFPP.hit.vmoment", sbs_gemFPP_hit_vmoment, &b_sbs_gemFPP_hit_vmoment);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.vsigma", &Ndata_sbs_gemFPP_hit_vsigma, &b_Ndata_sbs_gemFPP_hit_vsigma);
   fChain->SetBranchAddress("sbs.gemFPP.hit.vsigma", sbs_gemFPP_hit_vsigma, &b_sbs_gemFPP_hit_vsigma);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.vstriphi", &Ndata_sbs_gemFPP_hit_vstriphi, &b_Ndata_sbs_gemFPP_hit_vstriphi);
   fChain->SetBranchAddress("sbs.gemFPP.hit.vstriphi", sbs_gemFPP_hit_vstriphi, &b_sbs_gemFPP_hit_vstriphi);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.vstriplo", &Ndata_sbs_gemFPP_hit_vstriplo, &b_Ndata_sbs_gemFPP_hit_vstriplo);
   fChain->SetBranchAddress("sbs.gemFPP.hit.vstriplo", sbs_gemFPP_hit_vstriplo, &b_sbs_gemFPP_hit_vstriplo);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.vstripmax", &Ndata_sbs_gemFPP_hit_vstripmax, &b_Ndata_sbs_gemFPP_hit_vstripmax);
   fChain->SetBranchAddress("sbs.gemFPP.hit.vstripmax", sbs_gemFPP_hit_vstripmax, &b_sbs_gemFPP_hit_vstripmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.xglobal", &Ndata_sbs_gemFPP_hit_xglobal, &b_Ndata_sbs_gemFPP_hit_xglobal);
   fChain->SetBranchAddress("sbs.gemFPP.hit.xglobal", sbs_gemFPP_hit_xglobal, &b_sbs_gemFPP_hit_xglobal);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.xlocal", &Ndata_sbs_gemFPP_hit_xlocal, &b_Ndata_sbs_gemFPP_hit_xlocal);
   fChain->SetBranchAddress("sbs.gemFPP.hit.xlocal", sbs_gemFPP_hit_xlocal, &b_sbs_gemFPP_hit_xlocal);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.yglobal", &Ndata_sbs_gemFPP_hit_yglobal, &b_Ndata_sbs_gemFPP_hit_yglobal);
   fChain->SetBranchAddress("sbs.gemFPP.hit.yglobal", sbs_gemFPP_hit_yglobal, &b_sbs_gemFPP_hit_yglobal);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.ylocal", &Ndata_sbs_gemFPP_hit_ylocal, &b_Ndata_sbs_gemFPP_hit_ylocal);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ylocal", sbs_gemFPP_hit_ylocal, &b_sbs_gemFPP_hit_ylocal);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.hit.zglobal", &Ndata_sbs_gemFPP_hit_zglobal, &b_Ndata_sbs_gemFPP_hit_zglobal);
   fChain->SetBranchAddress("sbs.gemFPP.hit.zglobal", sbs_gemFPP_hit_zglobal, &b_sbs_gemFPP_hit_zglobal);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.n2Dhit_layer", &Ndata_sbs_gemFPP_n2Dhit_layer, &b_Ndata_sbs_gemFPP_n2Dhit_layer);
   fChain->SetBranchAddress("sbs.gemFPP.n2Dhit_layer", sbs_gemFPP_n2Dhit_layer, &b_sbs_gemFPP_n2Dhit_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.nclustu_layer", &Ndata_sbs_gemFPP_nclustu_layer, &b_Ndata_sbs_gemFPP_nclustu_layer);
   fChain->SetBranchAddress("sbs.gemFPP.nclustu_layer", sbs_gemFPP_nclustu_layer, &b_sbs_gemFPP_nclustu_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.nclustv_layer", &Ndata_sbs_gemFPP_nclustv_layer, &b_Ndata_sbs_gemFPP_nclustv_layer);
   fChain->SetBranchAddress("sbs.gemFPP.nclustv_layer", sbs_gemFPP_nclustv_layer, &b_sbs_gemFPP_nclustv_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.nstripsu_layer", &Ndata_sbs_gemFPP_nstripsu_layer, &b_Ndata_sbs_gemFPP_nstripsu_layer);
   fChain->SetBranchAddress("sbs.gemFPP.nstripsu_layer", sbs_gemFPP_nstripsu_layer, &b_sbs_gemFPP_nstripsu_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.nstripsv_layer", &Ndata_sbs_gemFPP_nstripsv_layer, &b_Ndata_sbs_gemFPP_nstripsv_layer);
   fChain->SetBranchAddress("sbs.gemFPP.nstripsv_layer", sbs_gemFPP_nstripsv_layer, &b_sbs_gemFPP_nstripsv_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.chi2ndf", &Ndata_sbs_gemFPP_track_chi2ndf, &b_Ndata_sbs_gemFPP_track_chi2ndf);
   fChain->SetBranchAddress("sbs.gemFPP.track.chi2ndf", sbs_gemFPP_track_chi2ndf, &b_sbs_gemFPP_track_chi2ndf);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.chi2ndf_hitquality", &Ndata_sbs_gemFPP_track_chi2ndf_hitquality, &b_Ndata_sbs_gemFPP_track_chi2ndf_hitquality);
   fChain->SetBranchAddress("sbs.gemFPP.track.chi2ndf_hitquality", sbs_gemFPP_track_chi2ndf_hitquality, &b_sbs_gemFPP_track_chi2ndf_hitquality);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.ngoodhits", &Ndata_sbs_gemFPP_track_ngoodhits, &b_Ndata_sbs_gemFPP_track_ngoodhits);
   fChain->SetBranchAddress("sbs.gemFPP.track.ngoodhits", sbs_gemFPP_track_ngoodhits, &b_sbs_gemFPP_track_ngoodhits);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.nhits", &Ndata_sbs_gemFPP_track_nhits, &b_Ndata_sbs_gemFPP_track_nhits);
   fChain->SetBranchAddress("sbs.gemFPP.track.nhits", sbs_gemFPP_track_nhits, &b_sbs_gemFPP_track_nhits);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.phi", &Ndata_sbs_gemFPP_track_phi, &b_Ndata_sbs_gemFPP_track_phi);
   fChain->SetBranchAddress("sbs.gemFPP.track.phi", sbs_gemFPP_track_phi, &b_sbs_gemFPP_track_phi);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.sclose", &Ndata_sbs_gemFPP_track_sclose, &b_Ndata_sbs_gemFPP_track_sclose);
   fChain->SetBranchAddress("sbs.gemFPP.track.sclose", sbs_gemFPP_track_sclose, &b_sbs_gemFPP_track_sclose);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.t0", &Ndata_sbs_gemFPP_track_t0, &b_Ndata_sbs_gemFPP_track_t0);
   fChain->SetBranchAddress("sbs.gemFPP.track.t0", sbs_gemFPP_track_t0, &b_sbs_gemFPP_track_t0);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.theta", &Ndata_sbs_gemFPP_track_theta, &b_Ndata_sbs_gemFPP_track_theta);
   fChain->SetBranchAddress("sbs.gemFPP.track.theta", sbs_gemFPP_track_theta, &b_sbs_gemFPP_track_theta);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.x", &Ndata_sbs_gemFPP_track_x, &b_Ndata_sbs_gemFPP_track_x);
   fChain->SetBranchAddress("sbs.gemFPP.track.x", sbs_gemFPP_track_x, &b_sbs_gemFPP_track_x);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.xp", &Ndata_sbs_gemFPP_track_xp, &b_Ndata_sbs_gemFPP_track_xp);
   fChain->SetBranchAddress("sbs.gemFPP.track.xp", sbs_gemFPP_track_xp, &b_sbs_gemFPP_track_xp);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.y", &Ndata_sbs_gemFPP_track_y, &b_Ndata_sbs_gemFPP_track_y);
   fChain->SetBranchAddress("sbs.gemFPP.track.y", sbs_gemFPP_track_y, &b_sbs_gemFPP_track_y);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.yp", &Ndata_sbs_gemFPP_track_yp, &b_Ndata_sbs_gemFPP_track_yp);
   fChain->SetBranchAddress("sbs.gemFPP.track.yp", sbs_gemFPP_track_yp, &b_sbs_gemFPP_track_yp);
   fChain->SetBranchAddress("Ndata.sbs.gemFPP.track.zclose", &Ndata_sbs_gemFPP_track_zclose, &b_Ndata_sbs_gemFPP_track_zclose);
   fChain->SetBranchAddress("sbs.gemFPP.track.zclose", sbs_gemFPP_track_zclose, &b_sbs_gemFPP_track_zclose);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCU", &Ndata_sbs_gemFT_hit_ADCU, &b_Ndata_sbs_gemFT_hit_ADCU);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCU", sbs_gemFT_hit_ADCU, &b_sbs_gemFT_hit_ADCU);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCU_deconv", &Ndata_sbs_gemFT_hit_ADCU_deconv, &b_Ndata_sbs_gemFT_hit_ADCU_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCU_deconv", sbs_gemFT_hit_ADCU_deconv, &b_sbs_gemFT_hit_ADCU_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCV", &Ndata_sbs_gemFT_hit_ADCV, &b_Ndata_sbs_gemFT_hit_ADCV);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCV", sbs_gemFT_hit_ADCV, &b_sbs_gemFT_hit_ADCV);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCV_deconv", &Ndata_sbs_gemFT_hit_ADCV_deconv, &b_Ndata_sbs_gemFT_hit_ADCV_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCV_deconv", sbs_gemFT_hit_ADCV_deconv, &b_sbs_gemFT_hit_ADCV_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCasym", &Ndata_sbs_gemFT_hit_ADCasym, &b_Ndata_sbs_gemFT_hit_ADCasym);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCasym", sbs_gemFT_hit_ADCasym, &b_sbs_gemFT_hit_ADCasym);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCasym_deconv", &Ndata_sbs_gemFT_hit_ADCasym_deconv, &b_Ndata_sbs_gemFT_hit_ADCasym_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCasym_deconv", sbs_gemFT_hit_ADCasym_deconv, &b_sbs_gemFT_hit_ADCasym_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCavg", &Ndata_sbs_gemFT_hit_ADCavg, &b_Ndata_sbs_gemFT_hit_ADCavg);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCavg", sbs_gemFT_hit_ADCavg, &b_sbs_gemFT_hit_ADCavg);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCavg_deconv", &Ndata_sbs_gemFT_hit_ADCavg_deconv, &b_Ndata_sbs_gemFT_hit_ADCavg_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCavg_deconv", sbs_gemFT_hit_ADCavg_deconv, &b_sbs_gemFT_hit_ADCavg_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac0_Umax", &Ndata_sbs_gemFT_hit_ADCfrac0_Umax, &b_Ndata_sbs_gemFT_hit_ADCfrac0_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac0_Umax", sbs_gemFT_hit_ADCfrac0_Umax, &b_sbs_gemFT_hit_ADCfrac0_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac0_Vmax", &Ndata_sbs_gemFT_hit_ADCfrac0_Vmax, &b_Ndata_sbs_gemFT_hit_ADCfrac0_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac0_Vmax", sbs_gemFT_hit_ADCfrac0_Vmax, &b_sbs_gemFT_hit_ADCfrac0_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac1_Umax", &Ndata_sbs_gemFT_hit_ADCfrac1_Umax, &b_Ndata_sbs_gemFT_hit_ADCfrac1_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac1_Umax", sbs_gemFT_hit_ADCfrac1_Umax, &b_sbs_gemFT_hit_ADCfrac1_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac1_Vmax", &Ndata_sbs_gemFT_hit_ADCfrac1_Vmax, &b_Ndata_sbs_gemFT_hit_ADCfrac1_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac1_Vmax", sbs_gemFT_hit_ADCfrac1_Vmax, &b_sbs_gemFT_hit_ADCfrac1_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac2_Umax", &Ndata_sbs_gemFT_hit_ADCfrac2_Umax, &b_Ndata_sbs_gemFT_hit_ADCfrac2_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac2_Umax", sbs_gemFT_hit_ADCfrac2_Umax, &b_sbs_gemFT_hit_ADCfrac2_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac2_Vmax", &Ndata_sbs_gemFT_hit_ADCfrac2_Vmax, &b_Ndata_sbs_gemFT_hit_ADCfrac2_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac2_Vmax", sbs_gemFT_hit_ADCfrac2_Vmax, &b_sbs_gemFT_hit_ADCfrac2_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac3_Umax", &Ndata_sbs_gemFT_hit_ADCfrac3_Umax, &b_Ndata_sbs_gemFT_hit_ADCfrac3_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac3_Umax", sbs_gemFT_hit_ADCfrac3_Umax, &b_sbs_gemFT_hit_ADCfrac3_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac3_Vmax", &Ndata_sbs_gemFT_hit_ADCfrac3_Vmax, &b_Ndata_sbs_gemFT_hit_ADCfrac3_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac3_Vmax", sbs_gemFT_hit_ADCfrac3_Vmax, &b_sbs_gemFT_hit_ADCfrac3_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac4_Umax", &Ndata_sbs_gemFT_hit_ADCfrac4_Umax, &b_Ndata_sbs_gemFT_hit_ADCfrac4_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac4_Umax", sbs_gemFT_hit_ADCfrac4_Umax, &b_sbs_gemFT_hit_ADCfrac4_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac4_Vmax", &Ndata_sbs_gemFT_hit_ADCfrac4_Vmax, &b_Ndata_sbs_gemFT_hit_ADCfrac4_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac4_Vmax", sbs_gemFT_hit_ADCfrac4_Vmax, &b_sbs_gemFT_hit_ADCfrac4_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac5_Umax", &Ndata_sbs_gemFT_hit_ADCfrac5_Umax, &b_Ndata_sbs_gemFT_hit_ADCfrac5_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac5_Umax", sbs_gemFT_hit_ADCfrac5_Umax, &b_sbs_gemFT_hit_ADCfrac5_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCfrac5_Vmax", &Ndata_sbs_gemFT_hit_ADCfrac5_Vmax, &b_Ndata_sbs_gemFT_hit_ADCfrac5_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCfrac5_Vmax", sbs_gemFT_hit_ADCfrac5_Vmax, &b_sbs_gemFT_hit_ADCfrac5_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxcomboUclust_deconv", &Ndata_sbs_gemFT_hit_ADCmaxcomboUclust_deconv, &b_Ndata_sbs_gemFT_hit_ADCmaxcomboUclust_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxcomboUclust_deconv", sbs_gemFT_hit_ADCmaxcomboUclust_deconv, &b_sbs_gemFT_hit_ADCmaxcomboUclust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxcomboVclust_deconv", &Ndata_sbs_gemFT_hit_ADCmaxcomboVclust_deconv, &b_Ndata_sbs_gemFT_hit_ADCmaxcomboVclust_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxcomboVclust_deconv", sbs_gemFT_hit_ADCmaxcomboVclust_deconv, &b_sbs_gemFT_hit_ADCmaxcomboVclust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxsampU", &Ndata_sbs_gemFT_hit_ADCmaxsampU, &b_Ndata_sbs_gemFT_hit_ADCmaxsampU);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxsampU", sbs_gemFT_hit_ADCmaxsampU, &b_sbs_gemFT_hit_ADCmaxsampU);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxsampUclust", &Ndata_sbs_gemFT_hit_ADCmaxsampUclust, &b_Ndata_sbs_gemFT_hit_ADCmaxsampUclust);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxsampUclust", sbs_gemFT_hit_ADCmaxsampUclust, &b_sbs_gemFT_hit_ADCmaxsampUclust);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxsampUclust_deconv", &Ndata_sbs_gemFT_hit_ADCmaxsampUclust_deconv, &b_Ndata_sbs_gemFT_hit_ADCmaxsampUclust_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxsampUclust_deconv", sbs_gemFT_hit_ADCmaxsampUclust_deconv, &b_sbs_gemFT_hit_ADCmaxsampUclust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxsampV", &Ndata_sbs_gemFT_hit_ADCmaxsampV, &b_Ndata_sbs_gemFT_hit_ADCmaxsampV);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxsampV", sbs_gemFT_hit_ADCmaxsampV, &b_sbs_gemFT_hit_ADCmaxsampV);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxsampVclust", &Ndata_sbs_gemFT_hit_ADCmaxsampVclust, &b_Ndata_sbs_gemFT_hit_ADCmaxsampVclust);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxsampVclust", sbs_gemFT_hit_ADCmaxsampVclust, &b_sbs_gemFT_hit_ADCmaxsampVclust);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxsampVclust_deconv", &Ndata_sbs_gemFT_hit_ADCmaxsampVclust_deconv, &b_Ndata_sbs_gemFT_hit_ADCmaxsampVclust_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxsampVclust_deconv", sbs_gemFT_hit_ADCmaxsampVclust_deconv, &b_sbs_gemFT_hit_ADCmaxsampVclust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxstripU", &Ndata_sbs_gemFT_hit_ADCmaxstripU, &b_Ndata_sbs_gemFT_hit_ADCmaxstripU);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxstripU", sbs_gemFT_hit_ADCmaxstripU, &b_sbs_gemFT_hit_ADCmaxstripU);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ADCmaxstripV", &Ndata_sbs_gemFT_hit_ADCmaxstripV, &b_Ndata_sbs_gemFT_hit_ADCmaxstripV);
   fChain->SetBranchAddress("sbs.gemFT.hit.ADCmaxstripV", sbs_gemFT_hit_ADCmaxstripV, &b_sbs_gemFT_hit_ADCmaxstripV);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.BUILD_ALL_SAMPLES_U", &Ndata_sbs_gemFT_hit_BUILD_ALL_SAMPLES_U, &b_Ndata_sbs_gemFT_hit_BUILD_ALL_SAMPLES_U);
   fChain->SetBranchAddress("sbs.gemFT.hit.BUILD_ALL_SAMPLES_U", sbs_gemFT_hit_BUILD_ALL_SAMPLES_U, &b_sbs_gemFT_hit_BUILD_ALL_SAMPLES_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.BUILD_ALL_SAMPLES_V", &Ndata_sbs_gemFT_hit_BUILD_ALL_SAMPLES_V, &b_Ndata_sbs_gemFT_hit_BUILD_ALL_SAMPLES_V);
   fChain->SetBranchAddress("sbs.gemFT.hit.BUILD_ALL_SAMPLES_V", sbs_gemFT_hit_BUILD_ALL_SAMPLES_V, &b_sbs_gemFT_hit_BUILD_ALL_SAMPLES_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.CM_GOOD_U", &Ndata_sbs_gemFT_hit_CM_GOOD_U, &b_Ndata_sbs_gemFT_hit_CM_GOOD_U);
   fChain->SetBranchAddress("sbs.gemFT.hit.CM_GOOD_U", sbs_gemFT_hit_CM_GOOD_U, &b_sbs_gemFT_hit_CM_GOOD_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.CM_GOOD_V", &Ndata_sbs_gemFT_hit_CM_GOOD_V, &b_Ndata_sbs_gemFT_hit_CM_GOOD_V);
   fChain->SetBranchAddress("sbs.gemFT.hit.CM_GOOD_V", sbs_gemFT_hit_CM_GOOD_V, &b_sbs_gemFT_hit_CM_GOOD_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC0_Umax", &Ndata_sbs_gemFT_hit_DeconvADC0_Umax, &b_Ndata_sbs_gemFT_hit_DeconvADC0_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC0_Umax", sbs_gemFT_hit_DeconvADC0_Umax, &b_sbs_gemFT_hit_DeconvADC0_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC0_Vmax", &Ndata_sbs_gemFT_hit_DeconvADC0_Vmax, &b_Ndata_sbs_gemFT_hit_DeconvADC0_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC0_Vmax", sbs_gemFT_hit_DeconvADC0_Vmax, &b_sbs_gemFT_hit_DeconvADC0_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC1_Umax", &Ndata_sbs_gemFT_hit_DeconvADC1_Umax, &b_Ndata_sbs_gemFT_hit_DeconvADC1_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC1_Umax", sbs_gemFT_hit_DeconvADC1_Umax, &b_sbs_gemFT_hit_DeconvADC1_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC1_Vmax", &Ndata_sbs_gemFT_hit_DeconvADC1_Vmax, &b_Ndata_sbs_gemFT_hit_DeconvADC1_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC1_Vmax", sbs_gemFT_hit_DeconvADC1_Vmax, &b_sbs_gemFT_hit_DeconvADC1_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC2_Umax", &Ndata_sbs_gemFT_hit_DeconvADC2_Umax, &b_Ndata_sbs_gemFT_hit_DeconvADC2_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC2_Umax", sbs_gemFT_hit_DeconvADC2_Umax, &b_sbs_gemFT_hit_DeconvADC2_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC2_Vmax", &Ndata_sbs_gemFT_hit_DeconvADC2_Vmax, &b_Ndata_sbs_gemFT_hit_DeconvADC2_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC2_Vmax", sbs_gemFT_hit_DeconvADC2_Vmax, &b_sbs_gemFT_hit_DeconvADC2_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC3_Umax", &Ndata_sbs_gemFT_hit_DeconvADC3_Umax, &b_Ndata_sbs_gemFT_hit_DeconvADC3_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC3_Umax", sbs_gemFT_hit_DeconvADC3_Umax, &b_sbs_gemFT_hit_DeconvADC3_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC3_Vmax", &Ndata_sbs_gemFT_hit_DeconvADC3_Vmax, &b_Ndata_sbs_gemFT_hit_DeconvADC3_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC3_Vmax", sbs_gemFT_hit_DeconvADC3_Vmax, &b_sbs_gemFT_hit_DeconvADC3_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC4_Umax", &Ndata_sbs_gemFT_hit_DeconvADC4_Umax, &b_Ndata_sbs_gemFT_hit_DeconvADC4_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC4_Umax", sbs_gemFT_hit_DeconvADC4_Umax, &b_sbs_gemFT_hit_DeconvADC4_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC4_Vmax", &Ndata_sbs_gemFT_hit_DeconvADC4_Vmax, &b_Ndata_sbs_gemFT_hit_DeconvADC4_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC4_Vmax", sbs_gemFT_hit_DeconvADC4_Vmax, &b_sbs_gemFT_hit_DeconvADC4_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC5_Umax", &Ndata_sbs_gemFT_hit_DeconvADC5_Umax, &b_Ndata_sbs_gemFT_hit_DeconvADC5_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC5_Umax", sbs_gemFT_hit_DeconvADC5_Umax, &b_sbs_gemFT_hit_DeconvADC5_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADC5_Vmax", &Ndata_sbs_gemFT_hit_DeconvADC5_Vmax, &b_Ndata_sbs_gemFT_hit_DeconvADC5_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADC5_Vmax", sbs_gemFT_hit_DeconvADC5_Vmax, &b_sbs_gemFT_hit_DeconvADC5_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADCmaxcomboU", &Ndata_sbs_gemFT_hit_DeconvADCmaxcomboU, &b_Ndata_sbs_gemFT_hit_DeconvADCmaxcomboU);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADCmaxcomboU", sbs_gemFT_hit_DeconvADCmaxcomboU, &b_sbs_gemFT_hit_DeconvADCmaxcomboU);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADCmaxcomboV", &Ndata_sbs_gemFT_hit_DeconvADCmaxcomboV, &b_Ndata_sbs_gemFT_hit_DeconvADCmaxcomboV);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADCmaxcomboV", sbs_gemFT_hit_DeconvADCmaxcomboV, &b_sbs_gemFT_hit_DeconvADCmaxcomboV);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADCmaxsampU", &Ndata_sbs_gemFT_hit_DeconvADCmaxsampU, &b_Ndata_sbs_gemFT_hit_DeconvADCmaxsampU);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADCmaxsampU", sbs_gemFT_hit_DeconvADCmaxsampU, &b_sbs_gemFT_hit_DeconvADCmaxsampU);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADCmaxsampV", &Ndata_sbs_gemFT_hit_DeconvADCmaxsampV, &b_Ndata_sbs_gemFT_hit_DeconvADCmaxsampV);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADCmaxsampV", sbs_gemFT_hit_DeconvADCmaxsampV, &b_sbs_gemFT_hit_DeconvADCmaxsampV);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADCmaxstripU", &Ndata_sbs_gemFT_hit_DeconvADCmaxstripU, &b_Ndata_sbs_gemFT_hit_DeconvADCmaxstripU);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADCmaxstripU", sbs_gemFT_hit_DeconvADCmaxstripU, &b_sbs_gemFT_hit_DeconvADCmaxstripU);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.DeconvADCmaxstripV", &Ndata_sbs_gemFT_hit_DeconvADCmaxstripV, &b_Ndata_sbs_gemFT_hit_DeconvADCmaxstripV);
   fChain->SetBranchAddress("sbs.gemFT.hit.DeconvADCmaxstripV", sbs_gemFT_hit_DeconvADCmaxstripV, &b_sbs_gemFT_hit_DeconvADCmaxstripV);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ENABLE_CM_U", &Ndata_sbs_gemFT_hit_ENABLE_CM_U, &b_Ndata_sbs_gemFT_hit_ENABLE_CM_U);
   fChain->SetBranchAddress("sbs.gemFT.hit.ENABLE_CM_U", sbs_gemFT_hit_ENABLE_CM_U, &b_sbs_gemFT_hit_ENABLE_CM_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ENABLE_CM_V", &Ndata_sbs_gemFT_hit_ENABLE_CM_V, &b_Ndata_sbs_gemFT_hit_ENABLE_CM_V);
   fChain->SetBranchAddress("sbs.gemFT.hit.ENABLE_CM_V", sbs_gemFT_hit_ENABLE_CM_V, &b_sbs_gemFT_hit_ENABLE_CM_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.TSchi2_Umax", &Ndata_sbs_gemFT_hit_TSchi2_Umax, &b_Ndata_sbs_gemFT_hit_TSchi2_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.TSchi2_Umax", sbs_gemFT_hit_TSchi2_Umax, &b_sbs_gemFT_hit_TSchi2_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.TSchi2_Vmax", &Ndata_sbs_gemFT_hit_TSchi2_Vmax, &b_Ndata_sbs_gemFT_hit_TSchi2_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.TSchi2_Vmax", sbs_gemFT_hit_TSchi2_Vmax, &b_sbs_gemFT_hit_TSchi2_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.TSprob_Umax", &Ndata_sbs_gemFT_hit_TSprob_Umax, &b_Ndata_sbs_gemFT_hit_TSprob_Umax);
   fChain->SetBranchAddress("sbs.gemFT.hit.TSprob_Umax", sbs_gemFT_hit_TSprob_Umax, &b_sbs_gemFT_hit_TSprob_Umax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.TSprob_Vmax", &Ndata_sbs_gemFT_hit_TSprob_Vmax, &b_Ndata_sbs_gemFT_hit_TSprob_Vmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.TSprob_Vmax", sbs_gemFT_hit_TSprob_Vmax, &b_sbs_gemFT_hit_TSprob_Vmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.Tavg", &Ndata_sbs_gemFT_hit_Tavg, &b_Ndata_sbs_gemFT_hit_Tavg);
   fChain->SetBranchAddress("sbs.gemFT.hit.Tavg", sbs_gemFT_hit_Tavg, &b_sbs_gemFT_hit_Tavg);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.Tavg_corr", &Ndata_sbs_gemFT_hit_Tavg_corr, &b_Ndata_sbs_gemFT_hit_Tavg_corr);
   fChain->SetBranchAddress("sbs.gemFT.hit.Tavg_corr", sbs_gemFT_hit_Tavg_corr, &b_sbs_gemFT_hit_Tavg_corr);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.Tavg_deconv", &Ndata_sbs_gemFT_hit_Tavg_deconv, &b_Ndata_sbs_gemFT_hit_Tavg_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.Tavg_deconv", sbs_gemFT_hit_Tavg_deconv, &b_sbs_gemFT_hit_Tavg_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.Tavg_fit", &Ndata_sbs_gemFT_hit_Tavg_fit, &b_Ndata_sbs_gemFT_hit_Tavg_fit);
   fChain->SetBranchAddress("sbs.gemFT.hit.Tavg_fit", sbs_gemFT_hit_Tavg_fit, &b_sbs_gemFT_hit_Tavg_fit);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.Ugain", &Ndata_sbs_gemFT_hit_Ugain, &b_Ndata_sbs_gemFT_hit_Ugain);
   fChain->SetBranchAddress("sbs.gemFT.hit.Ugain", sbs_gemFT_hit_Ugain, &b_sbs_gemFT_hit_Ugain);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.Utime", &Ndata_sbs_gemFT_hit_Utime, &b_Ndata_sbs_gemFT_hit_Utime);
   fChain->SetBranchAddress("sbs.gemFT.hit.Utime", sbs_gemFT_hit_Utime, &b_sbs_gemFT_hit_Utime);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.UtimeDeconv", &Ndata_sbs_gemFT_hit_UtimeDeconv, &b_Ndata_sbs_gemFT_hit_UtimeDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.UtimeDeconv", sbs_gemFT_hit_UtimeDeconv, &b_sbs_gemFT_hit_UtimeDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.UtimeFit", &Ndata_sbs_gemFT_hit_UtimeFit, &b_Ndata_sbs_gemFT_hit_UtimeFit);
   fChain->SetBranchAddress("sbs.gemFT.hit.UtimeFit", sbs_gemFT_hit_UtimeFit, &b_sbs_gemFT_hit_UtimeFit);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.UtimeMaxStrip", &Ndata_sbs_gemFT_hit_UtimeMaxStrip, &b_Ndata_sbs_gemFT_hit_UtimeMaxStrip);
   fChain->SetBranchAddress("sbs.gemFT.hit.UtimeMaxStrip", sbs_gemFT_hit_UtimeMaxStrip, &b_sbs_gemFT_hit_UtimeMaxStrip);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.UtimeMaxStripDeconv", &Ndata_sbs_gemFT_hit_UtimeMaxStripDeconv, &b_Ndata_sbs_gemFT_hit_UtimeMaxStripDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.UtimeMaxStripDeconv", sbs_gemFT_hit_UtimeMaxStripDeconv, &b_sbs_gemFT_hit_UtimeMaxStripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.UtimeMaxStripFit", &Ndata_sbs_gemFT_hit_UtimeMaxStripFit, &b_Ndata_sbs_gemFT_hit_UtimeMaxStripFit);
   fChain->SetBranchAddress("sbs.gemFT.hit.UtimeMaxStripFit", sbs_gemFT_hit_UtimeMaxStripFit, &b_sbs_gemFT_hit_UtimeMaxStripFit);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.Vgain", &Ndata_sbs_gemFT_hit_Vgain, &b_Ndata_sbs_gemFT_hit_Vgain);
   fChain->SetBranchAddress("sbs.gemFT.hit.Vgain", sbs_gemFT_hit_Vgain, &b_sbs_gemFT_hit_Vgain);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.Vtime", &Ndata_sbs_gemFT_hit_Vtime, &b_Ndata_sbs_gemFT_hit_Vtime);
   fChain->SetBranchAddress("sbs.gemFT.hit.Vtime", sbs_gemFT_hit_Vtime, &b_sbs_gemFT_hit_Vtime);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.VtimeDeconv", &Ndata_sbs_gemFT_hit_VtimeDeconv, &b_Ndata_sbs_gemFT_hit_VtimeDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.VtimeDeconv", sbs_gemFT_hit_VtimeDeconv, &b_sbs_gemFT_hit_VtimeDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.VtimeFit", &Ndata_sbs_gemFT_hit_VtimeFit, &b_Ndata_sbs_gemFT_hit_VtimeFit);
   fChain->SetBranchAddress("sbs.gemFT.hit.VtimeFit", sbs_gemFT_hit_VtimeFit, &b_sbs_gemFT_hit_VtimeFit);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.VtimeMaxStrip", &Ndata_sbs_gemFT_hit_VtimeMaxStrip, &b_Ndata_sbs_gemFT_hit_VtimeMaxStrip);
   fChain->SetBranchAddress("sbs.gemFT.hit.VtimeMaxStrip", sbs_gemFT_hit_VtimeMaxStrip, &b_sbs_gemFT_hit_VtimeMaxStrip);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.VtimeMaxStripDeconv", &Ndata_sbs_gemFT_hit_VtimeMaxStripDeconv, &b_Ndata_sbs_gemFT_hit_VtimeMaxStripDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.VtimeMaxStripDeconv", sbs_gemFT_hit_VtimeMaxStripDeconv, &b_sbs_gemFT_hit_VtimeMaxStripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.VtimeMaxStripFit", &Ndata_sbs_gemFT_hit_VtimeMaxStripFit, &b_Ndata_sbs_gemFT_hit_VtimeMaxStripFit);
   fChain->SetBranchAddress("sbs.gemFT.hit.VtimeMaxStripFit", sbs_gemFT_hit_VtimeMaxStripFit, &b_sbs_gemFT_hit_VtimeMaxStripFit);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.adc_id_U", &Ndata_sbs_gemFT_hit_adc_id_U, &b_Ndata_sbs_gemFT_hit_adc_id_U);
   fChain->SetBranchAddress("sbs.gemFT.hit.adc_id_U", sbs_gemFT_hit_adc_id_U, &b_sbs_gemFT_hit_adc_id_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.adc_id_V", &Ndata_sbs_gemFT_hit_adc_id_V, &b_Ndata_sbs_gemFT_hit_adc_id_V);
   fChain->SetBranchAddress("sbs.gemFT.hit.adc_id_V", sbs_gemFT_hit_adc_id_V, &b_sbs_gemFT_hit_adc_id_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ccor_clust", &Ndata_sbs_gemFT_hit_ccor_clust, &b_Ndata_sbs_gemFT_hit_ccor_clust);
   fChain->SetBranchAddress("sbs.gemFT.hit.ccor_clust", sbs_gemFT_hit_ccor_clust, &b_sbs_gemFT_hit_ccor_clust);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ccor_clust_deconv", &Ndata_sbs_gemFT_hit_ccor_clust_deconv, &b_Ndata_sbs_gemFT_hit_ccor_clust_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ccor_clust_deconv", sbs_gemFT_hit_ccor_clust_deconv, &b_sbs_gemFT_hit_ccor_clust_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ccor_strip", &Ndata_sbs_gemFT_hit_ccor_strip, &b_Ndata_sbs_gemFT_hit_ccor_strip);
   fChain->SetBranchAddress("sbs.gemFT.hit.ccor_strip", sbs_gemFT_hit_ccor_strip, &b_sbs_gemFT_hit_ccor_strip);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ccor_strip_deconv", &Ndata_sbs_gemFT_hit_ccor_strip_deconv, &b_Ndata_sbs_gemFT_hit_ccor_strip_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.ccor_strip_deconv", sbs_gemFT_hit_ccor_strip_deconv, &b_sbs_gemFT_hit_ccor_strip_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.crate_U", &Ndata_sbs_gemFT_hit_crate_U, &b_Ndata_sbs_gemFT_hit_crate_U);
   fChain->SetBranchAddress("sbs.gemFT.hit.crate_U", sbs_gemFT_hit_crate_U, &b_sbs_gemFT_hit_crate_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.crate_V", &Ndata_sbs_gemFT_hit_crate_V, &b_Ndata_sbs_gemFT_hit_crate_V);
   fChain->SetBranchAddress("sbs.gemFT.hit.crate_V", sbs_gemFT_hit_crate_V, &b_sbs_gemFT_hit_crate_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.deltat", &Ndata_sbs_gemFT_hit_deltat, &b_Ndata_sbs_gemFT_hit_deltat);
   fChain->SetBranchAddress("sbs.gemFT.hit.deltat", sbs_gemFT_hit_deltat, &b_sbs_gemFT_hit_deltat);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.deltat_deconv", &Ndata_sbs_gemFT_hit_deltat_deconv, &b_Ndata_sbs_gemFT_hit_deltat_deconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.deltat_deconv", sbs_gemFT_hit_deltat_deconv, &b_sbs_gemFT_hit_deltat_deconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.deltat_fit", &Ndata_sbs_gemFT_hit_deltat_fit, &b_Ndata_sbs_gemFT_hit_deltat_fit);
   fChain->SetBranchAddress("sbs.gemFT.hit.deltat_fit", sbs_gemFT_hit_deltat_fit, &b_sbs_gemFT_hit_deltat_fit);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.eresidu", &Ndata_sbs_gemFT_hit_eresidu, &b_Ndata_sbs_gemFT_hit_eresidu);
   fChain->SetBranchAddress("sbs.gemFT.hit.eresidu", sbs_gemFT_hit_eresidu, &b_sbs_gemFT_hit_eresidu);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.eresidv", &Ndata_sbs_gemFT_hit_eresidv, &b_Ndata_sbs_gemFT_hit_eresidv);
   fChain->SetBranchAddress("sbs.gemFT.hit.eresidv", sbs_gemFT_hit_eresidv, &b_sbs_gemFT_hit_eresidv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.icombomaxUclustDeconv", &Ndata_sbs_gemFT_hit_icombomaxUclustDeconv, &b_Ndata_sbs_gemFT_hit_icombomaxUclustDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.icombomaxUclustDeconv", sbs_gemFT_hit_icombomaxUclustDeconv, &b_sbs_gemFT_hit_icombomaxUclustDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.icombomaxUstripDeconv", &Ndata_sbs_gemFT_hit_icombomaxUstripDeconv, &b_Ndata_sbs_gemFT_hit_icombomaxUstripDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.icombomaxUstripDeconv", sbs_gemFT_hit_icombomaxUstripDeconv, &b_sbs_gemFT_hit_icombomaxUstripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.icombomaxVclustDeconv", &Ndata_sbs_gemFT_hit_icombomaxVclustDeconv, &b_Ndata_sbs_gemFT_hit_icombomaxVclustDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.icombomaxVclustDeconv", sbs_gemFT_hit_icombomaxVclustDeconv, &b_sbs_gemFT_hit_icombomaxVclustDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.icombomaxVstripDeconv", &Ndata_sbs_gemFT_hit_icombomaxVstripDeconv, &b_Ndata_sbs_gemFT_hit_icombomaxVstripDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.icombomaxVstripDeconv", sbs_gemFT_hit_icombomaxVstripDeconv, &b_sbs_gemFT_hit_icombomaxVstripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.isampmaxUclust", &Ndata_sbs_gemFT_hit_isampmaxUclust, &b_Ndata_sbs_gemFT_hit_isampmaxUclust);
   fChain->SetBranchAddress("sbs.gemFT.hit.isampmaxUclust", sbs_gemFT_hit_isampmaxUclust, &b_sbs_gemFT_hit_isampmaxUclust);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.isampmaxUclustDeconv", &Ndata_sbs_gemFT_hit_isampmaxUclustDeconv, &b_Ndata_sbs_gemFT_hit_isampmaxUclustDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.isampmaxUclustDeconv", sbs_gemFT_hit_isampmaxUclustDeconv, &b_sbs_gemFT_hit_isampmaxUclustDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.isampmaxUstrip", &Ndata_sbs_gemFT_hit_isampmaxUstrip, &b_Ndata_sbs_gemFT_hit_isampmaxUstrip);
   fChain->SetBranchAddress("sbs.gemFT.hit.isampmaxUstrip", sbs_gemFT_hit_isampmaxUstrip, &b_sbs_gemFT_hit_isampmaxUstrip);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.isampmaxUstripDeconv", &Ndata_sbs_gemFT_hit_isampmaxUstripDeconv, &b_Ndata_sbs_gemFT_hit_isampmaxUstripDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.isampmaxUstripDeconv", sbs_gemFT_hit_isampmaxUstripDeconv, &b_sbs_gemFT_hit_isampmaxUstripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.isampmaxVclust", &Ndata_sbs_gemFT_hit_isampmaxVclust, &b_Ndata_sbs_gemFT_hit_isampmaxVclust);
   fChain->SetBranchAddress("sbs.gemFT.hit.isampmaxVclust", sbs_gemFT_hit_isampmaxVclust, &b_sbs_gemFT_hit_isampmaxVclust);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.isampmaxVclustDeconv", &Ndata_sbs_gemFT_hit_isampmaxVclustDeconv, &b_Ndata_sbs_gemFT_hit_isampmaxVclustDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.isampmaxVclustDeconv", sbs_gemFT_hit_isampmaxVclustDeconv, &b_sbs_gemFT_hit_isampmaxVclustDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.isampmaxVstrip", &Ndata_sbs_gemFT_hit_isampmaxVstrip, &b_Ndata_sbs_gemFT_hit_isampmaxVstrip);
   fChain->SetBranchAddress("sbs.gemFT.hit.isampmaxVstrip", sbs_gemFT_hit_isampmaxVstrip, &b_sbs_gemFT_hit_isampmaxVstrip);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.isampmaxVstripDeconv", &Ndata_sbs_gemFT_hit_isampmaxVstripDeconv, &b_Ndata_sbs_gemFT_hit_isampmaxVstripDeconv);
   fChain->SetBranchAddress("sbs.gemFT.hit.isampmaxVstripDeconv", sbs_gemFT_hit_isampmaxVstripDeconv, &b_sbs_gemFT_hit_isampmaxVstripDeconv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.layer", &Ndata_sbs_gemFT_hit_layer, &b_Ndata_sbs_gemFT_hit_layer);
   fChain->SetBranchAddress("sbs.gemFT.hit.layer", sbs_gemFT_hit_layer, &b_sbs_gemFT_hit_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.module", &Ndata_sbs_gemFT_hit_module, &b_Ndata_sbs_gemFT_hit_module);
   fChain->SetBranchAddress("sbs.gemFT.hit.module", sbs_gemFT_hit_module, &b_sbs_gemFT_hit_module);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.mpd_U", &Ndata_sbs_gemFT_hit_mpd_U, &b_Ndata_sbs_gemFT_hit_mpd_U);
   fChain->SetBranchAddress("sbs.gemFT.hit.mpd_U", sbs_gemFT_hit_mpd_U, &b_sbs_gemFT_hit_mpd_U);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.mpd_V", &Ndata_sbs_gemFT_hit_mpd_V, &b_Ndata_sbs_gemFT_hit_mpd_V);
   fChain->SetBranchAddress("sbs.gemFT.hit.mpd_V", sbs_gemFT_hit_mpd_V, &b_sbs_gemFT_hit_mpd_V);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.nstripu", &Ndata_sbs_gemFT_hit_nstripu, &b_Ndata_sbs_gemFT_hit_nstripu);
   fChain->SetBranchAddress("sbs.gemFT.hit.nstripu", sbs_gemFT_hit_nstripu, &b_sbs_gemFT_hit_nstripu);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.nstripv", &Ndata_sbs_gemFT_hit_nstripv, &b_Ndata_sbs_gemFT_hit_nstripv);
   fChain->SetBranchAddress("sbs.gemFT.hit.nstripv", sbs_gemFT_hit_nstripv, &b_sbs_gemFT_hit_nstripv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.residu", &Ndata_sbs_gemFT_hit_residu, &b_Ndata_sbs_gemFT_hit_residu);
   fChain->SetBranchAddress("sbs.gemFT.hit.residu", sbs_gemFT_hit_residu, &b_sbs_gemFT_hit_residu);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.residv", &Ndata_sbs_gemFT_hit_residv, &b_Ndata_sbs_gemFT_hit_residv);
   fChain->SetBranchAddress("sbs.gemFT.hit.residv", sbs_gemFT_hit_residv, &b_sbs_gemFT_hit_residv);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.trackindex", &Ndata_sbs_gemFT_hit_trackindex, &b_Ndata_sbs_gemFT_hit_trackindex);
   fChain->SetBranchAddress("sbs.gemFT.hit.trackindex", sbs_gemFT_hit_trackindex, &b_sbs_gemFT_hit_trackindex);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.u", &Ndata_sbs_gemFT_hit_u, &b_Ndata_sbs_gemFT_hit_u);
   fChain->SetBranchAddress("sbs.gemFT.hit.u", sbs_gemFT_hit_u, &b_sbs_gemFT_hit_u);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.umoment", &Ndata_sbs_gemFT_hit_umoment, &b_Ndata_sbs_gemFT_hit_umoment);
   fChain->SetBranchAddress("sbs.gemFT.hit.umoment", sbs_gemFT_hit_umoment, &b_sbs_gemFT_hit_umoment);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.usigma", &Ndata_sbs_gemFT_hit_usigma, &b_Ndata_sbs_gemFT_hit_usigma);
   fChain->SetBranchAddress("sbs.gemFT.hit.usigma", sbs_gemFT_hit_usigma, &b_sbs_gemFT_hit_usigma);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ustriphi", &Ndata_sbs_gemFT_hit_ustriphi, &b_Ndata_sbs_gemFT_hit_ustriphi);
   fChain->SetBranchAddress("sbs.gemFT.hit.ustriphi", sbs_gemFT_hit_ustriphi, &b_sbs_gemFT_hit_ustriphi);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ustriplo", &Ndata_sbs_gemFT_hit_ustriplo, &b_Ndata_sbs_gemFT_hit_ustriplo);
   fChain->SetBranchAddress("sbs.gemFT.hit.ustriplo", sbs_gemFT_hit_ustriplo, &b_sbs_gemFT_hit_ustriplo);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ustripmax", &Ndata_sbs_gemFT_hit_ustripmax, &b_Ndata_sbs_gemFT_hit_ustripmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.ustripmax", sbs_gemFT_hit_ustripmax, &b_sbs_gemFT_hit_ustripmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.v", &Ndata_sbs_gemFT_hit_v, &b_Ndata_sbs_gemFT_hit_v);
   fChain->SetBranchAddress("sbs.gemFT.hit.v", sbs_gemFT_hit_v, &b_sbs_gemFT_hit_v);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.vmoment", &Ndata_sbs_gemFT_hit_vmoment, &b_Ndata_sbs_gemFT_hit_vmoment);
   fChain->SetBranchAddress("sbs.gemFT.hit.vmoment", sbs_gemFT_hit_vmoment, &b_sbs_gemFT_hit_vmoment);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.vsigma", &Ndata_sbs_gemFT_hit_vsigma, &b_Ndata_sbs_gemFT_hit_vsigma);
   fChain->SetBranchAddress("sbs.gemFT.hit.vsigma", sbs_gemFT_hit_vsigma, &b_sbs_gemFT_hit_vsigma);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.vstriphi", &Ndata_sbs_gemFT_hit_vstriphi, &b_Ndata_sbs_gemFT_hit_vstriphi);
   fChain->SetBranchAddress("sbs.gemFT.hit.vstriphi", sbs_gemFT_hit_vstriphi, &b_sbs_gemFT_hit_vstriphi);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.vstriplo", &Ndata_sbs_gemFT_hit_vstriplo, &b_Ndata_sbs_gemFT_hit_vstriplo);
   fChain->SetBranchAddress("sbs.gemFT.hit.vstriplo", sbs_gemFT_hit_vstriplo, &b_sbs_gemFT_hit_vstriplo);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.vstripmax", &Ndata_sbs_gemFT_hit_vstripmax, &b_Ndata_sbs_gemFT_hit_vstripmax);
   fChain->SetBranchAddress("sbs.gemFT.hit.vstripmax", sbs_gemFT_hit_vstripmax, &b_sbs_gemFT_hit_vstripmax);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.xglobal", &Ndata_sbs_gemFT_hit_xglobal, &b_Ndata_sbs_gemFT_hit_xglobal);
   fChain->SetBranchAddress("sbs.gemFT.hit.xglobal", sbs_gemFT_hit_xglobal, &b_sbs_gemFT_hit_xglobal);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.xlocal", &Ndata_sbs_gemFT_hit_xlocal, &b_Ndata_sbs_gemFT_hit_xlocal);
   fChain->SetBranchAddress("sbs.gemFT.hit.xlocal", sbs_gemFT_hit_xlocal, &b_sbs_gemFT_hit_xlocal);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.yglobal", &Ndata_sbs_gemFT_hit_yglobal, &b_Ndata_sbs_gemFT_hit_yglobal);
   fChain->SetBranchAddress("sbs.gemFT.hit.yglobal", sbs_gemFT_hit_yglobal, &b_sbs_gemFT_hit_yglobal);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.ylocal", &Ndata_sbs_gemFT_hit_ylocal, &b_Ndata_sbs_gemFT_hit_ylocal);
   fChain->SetBranchAddress("sbs.gemFT.hit.ylocal", sbs_gemFT_hit_ylocal, &b_sbs_gemFT_hit_ylocal);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.hit.zglobal", &Ndata_sbs_gemFT_hit_zglobal, &b_Ndata_sbs_gemFT_hit_zglobal);
   fChain->SetBranchAddress("sbs.gemFT.hit.zglobal", sbs_gemFT_hit_zglobal, &b_sbs_gemFT_hit_zglobal);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.n2Dhit_layer", &Ndata_sbs_gemFT_n2Dhit_layer, &b_Ndata_sbs_gemFT_n2Dhit_layer);
   fChain->SetBranchAddress("sbs.gemFT.n2Dhit_layer", sbs_gemFT_n2Dhit_layer, &b_sbs_gemFT_n2Dhit_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.nclustu_layer", &Ndata_sbs_gemFT_nclustu_layer, &b_Ndata_sbs_gemFT_nclustu_layer);
   fChain->SetBranchAddress("sbs.gemFT.nclustu_layer", sbs_gemFT_nclustu_layer, &b_sbs_gemFT_nclustu_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.nclustv_layer", &Ndata_sbs_gemFT_nclustv_layer, &b_Ndata_sbs_gemFT_nclustv_layer);
   fChain->SetBranchAddress("sbs.gemFT.nclustv_layer", sbs_gemFT_nclustv_layer, &b_sbs_gemFT_nclustv_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.nstripsu_layer", &Ndata_sbs_gemFT_nstripsu_layer, &b_Ndata_sbs_gemFT_nstripsu_layer);
   fChain->SetBranchAddress("sbs.gemFT.nstripsu_layer", sbs_gemFT_nstripsu_layer, &b_sbs_gemFT_nstripsu_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.nstripsv_layer", &Ndata_sbs_gemFT_nstripsv_layer, &b_Ndata_sbs_gemFT_nstripsv_layer);
   fChain->SetBranchAddress("sbs.gemFT.nstripsv_layer", sbs_gemFT_nstripsv_layer, &b_sbs_gemFT_nstripsv_layer);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.track.chi2ndf", &Ndata_sbs_gemFT_track_chi2ndf, &b_Ndata_sbs_gemFT_track_chi2ndf);
   fChain->SetBranchAddress("sbs.gemFT.track.chi2ndf", sbs_gemFT_track_chi2ndf, &b_sbs_gemFT_track_chi2ndf);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.track.chi2ndf_hitquality", &Ndata_sbs_gemFT_track_chi2ndf_hitquality, &b_Ndata_sbs_gemFT_track_chi2ndf_hitquality);
   fChain->SetBranchAddress("sbs.gemFT.track.chi2ndf_hitquality", sbs_gemFT_track_chi2ndf_hitquality, &b_sbs_gemFT_track_chi2ndf_hitquality);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.track.ngoodhits", &Ndata_sbs_gemFT_track_ngoodhits, &b_Ndata_sbs_gemFT_track_ngoodhits);
   fChain->SetBranchAddress("sbs.gemFT.track.ngoodhits", sbs_gemFT_track_ngoodhits, &b_sbs_gemFT_track_ngoodhits);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.track.nhits", &Ndata_sbs_gemFT_track_nhits, &b_Ndata_sbs_gemFT_track_nhits);
   fChain->SetBranchAddress("sbs.gemFT.track.nhits", sbs_gemFT_track_nhits, &b_sbs_gemFT_track_nhits);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.track.t0", &Ndata_sbs_gemFT_track_t0, &b_Ndata_sbs_gemFT_track_t0);
   fChain->SetBranchAddress("sbs.gemFT.track.t0", sbs_gemFT_track_t0, &b_sbs_gemFT_track_t0);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.track.x", &Ndata_sbs_gemFT_track_x, &b_Ndata_sbs_gemFT_track_x);
   fChain->SetBranchAddress("sbs.gemFT.track.x", sbs_gemFT_track_x, &b_sbs_gemFT_track_x);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.track.xp", &Ndata_sbs_gemFT_track_xp, &b_Ndata_sbs_gemFT_track_xp);
   fChain->SetBranchAddress("sbs.gemFT.track.xp", sbs_gemFT_track_xp, &b_sbs_gemFT_track_xp);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.track.y", &Ndata_sbs_gemFT_track_y, &b_Ndata_sbs_gemFT_track_y);
   fChain->SetBranchAddress("sbs.gemFT.track.y", sbs_gemFT_track_y, &b_sbs_gemFT_track_y);
   fChain->SetBranchAddress("Ndata.sbs.gemFT.track.yp", &Ndata_sbs_gemFT_track_yp, &b_Ndata_sbs_gemFT_track_yp);
   fChain->SetBranchAddress("sbs.gemFT.track.yp", sbs_gemFT_track_yp, &b_sbs_gemFT_track_yp);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.adctime", &Ndata_sbs_hcal_clus_adctime, &b_Ndata_sbs_hcal_clus_adctime);
   fChain->SetBranchAddress("sbs.hcal.clus.adctime", sbs_hcal_clus_adctime, &b_sbs_hcal_clus_adctime);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.again", &Ndata_sbs_hcal_clus_again, &b_Ndata_sbs_hcal_clus_again);
   fChain->SetBranchAddress("sbs.hcal.clus.again", sbs_hcal_clus_again, &b_sbs_hcal_clus_again);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.atimeblk", &Ndata_sbs_hcal_clus_atimeblk, &b_Ndata_sbs_hcal_clus_atimeblk);
   fChain->SetBranchAddress("sbs.hcal.clus.atimeblk", sbs_hcal_clus_atimeblk, &b_sbs_hcal_clus_atimeblk);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.col", &Ndata_sbs_hcal_clus_col, &b_Ndata_sbs_hcal_clus_col);
   fChain->SetBranchAddress("sbs.hcal.clus.col", sbs_hcal_clus_col, &b_sbs_hcal_clus_col);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.col_goodtdc", &Ndata_sbs_hcal_clus_col_goodtdc, &b_Ndata_sbs_hcal_clus_col_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.clus.col_goodtdc", sbs_hcal_clus_col_goodtdc, &b_sbs_hcal_clus_col_goodtdc);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.e", &Ndata_sbs_hcal_clus_e, &b_Ndata_sbs_hcal_clus_e);
   fChain->SetBranchAddress("sbs.hcal.clus.e", sbs_hcal_clus_e, &b_sbs_hcal_clus_e);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.e_goodtdc", &Ndata_sbs_hcal_clus_e_goodtdc, &b_Ndata_sbs_hcal_clus_e_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.clus.e_goodtdc", sbs_hcal_clus_e_goodtdc, &b_sbs_hcal_clus_e_goodtdc);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.eblk", &Ndata_sbs_hcal_clus_eblk, &b_Ndata_sbs_hcal_clus_eblk);
   fChain->SetBranchAddress("sbs.hcal.clus.eblk", sbs_hcal_clus_eblk, &b_sbs_hcal_clus_eblk);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.eblk_goodtdc", &Ndata_sbs_hcal_clus_eblk_goodtdc, &b_Ndata_sbs_hcal_clus_eblk_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.clus.eblk_goodtdc", sbs_hcal_clus_eblk_goodtdc, &b_sbs_hcal_clus_eblk_goodtdc);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.id", &Ndata_sbs_hcal_clus_id, &b_Ndata_sbs_hcal_clus_id);
   fChain->SetBranchAddress("sbs.hcal.clus.id", sbs_hcal_clus_id, &b_sbs_hcal_clus_id);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.id_goodtdc", &Ndata_sbs_hcal_clus_id_goodtdc, &b_Ndata_sbs_hcal_clus_id_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.clus.id_goodtdc", sbs_hcal_clus_id_goodtdc, &b_sbs_hcal_clus_id_goodtdc);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.nblk", &Ndata_sbs_hcal_clus_nblk, &b_Ndata_sbs_hcal_clus_nblk);
   fChain->SetBranchAddress("sbs.hcal.clus.nblk", sbs_hcal_clus_nblk, &b_sbs_hcal_clus_nblk);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.nblk_goodtdc", &Ndata_sbs_hcal_clus_nblk_goodtdc, &b_Ndata_sbs_hcal_clus_nblk_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.clus.nblk_goodtdc", sbs_hcal_clus_nblk_goodtdc, &b_sbs_hcal_clus_nblk_goodtdc);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.row", &Ndata_sbs_hcal_clus_row, &b_Ndata_sbs_hcal_clus_row);
   fChain->SetBranchAddress("sbs.hcal.clus.row", sbs_hcal_clus_row, &b_sbs_hcal_clus_row);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.row_goodtdc", &Ndata_sbs_hcal_clus_row_goodtdc, &b_Ndata_sbs_hcal_clus_row_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.clus.row_goodtdc", sbs_hcal_clus_row_goodtdc, &b_sbs_hcal_clus_row_goodtdc);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.tdctime", &Ndata_sbs_hcal_clus_tdctime, &b_Ndata_sbs_hcal_clus_tdctime);
   fChain->SetBranchAddress("sbs.hcal.clus.tdctime", sbs_hcal_clus_tdctime, &b_sbs_hcal_clus_tdctime);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.tdctime_tw", &Ndata_sbs_hcal_clus_tdctime_tw, &b_Ndata_sbs_hcal_clus_tdctime_tw);
   fChain->SetBranchAddress("sbs.hcal.clus.tdctime_tw", sbs_hcal_clus_tdctime_tw, &b_sbs_hcal_clus_tdctime_tw);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.tdctimeblk", &Ndata_sbs_hcal_clus_tdctimeblk, &b_Ndata_sbs_hcal_clus_tdctimeblk);
   fChain->SetBranchAddress("sbs.hcal.clus.tdctimeblk", sbs_hcal_clus_tdctimeblk, &b_sbs_hcal_clus_tdctimeblk);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.tdctimeblk_tw", &Ndata_sbs_hcal_clus_tdctimeblk_tw, &b_Ndata_sbs_hcal_clus_tdctimeblk_tw);
   fChain->SetBranchAddress("sbs.hcal.clus.tdctimeblk_tw", sbs_hcal_clus_tdctimeblk_tw, &b_sbs_hcal_clus_tdctimeblk_tw);
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
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.tdctime_tw", &Ndata_sbs_hcal_clus_blk_tdctime_tw, &b_Ndata_sbs_hcal_clus_blk_tdctime_tw);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.tdctime_tw", sbs_hcal_clus_blk_tdctime_tw, &b_sbs_hcal_clus_blk_tdctime_tw);
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
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.tdctime_tw", &Ndata_sbs_hcal_goodblock_tdctime_tw, &b_Ndata_sbs_hcal_goodblock_tdctime_tw);
   fChain->SetBranchAddress("sbs.hcal.goodblock.tdctime_tw", sbs_hcal_goodblock_tdctime_tw, &b_sbs_hcal_goodblock_tdctime_tw);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.x", &Ndata_sbs_hcal_goodblock_x, &b_Ndata_sbs_hcal_goodblock_x);
   fChain->SetBranchAddress("sbs.hcal.goodblock.x", sbs_hcal_goodblock_x, &b_sbs_hcal_goodblock_x);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.y", &Ndata_sbs_hcal_goodblock_y, &b_Ndata_sbs_hcal_goodblock_y);
   fChain->SetBranchAddress("sbs.hcal.goodblock.y", sbs_hcal_goodblock_y, &b_sbs_hcal_goodblock_y);
   fChain->SetBranchAddress("Ndata.sbs.tr.beta", &Ndata_sbs_tr_beta, &b_Ndata_sbs_tr_beta);
   fChain->SetBranchAddress("sbs.tr.beta", sbs_tr_beta, &b_sbs_tr_beta);
   fChain->SetBranchAddress("Ndata.sbs.tr.chi2", &Ndata_sbs_tr_chi2, &b_Ndata_sbs_tr_chi2);
   fChain->SetBranchAddress("sbs.tr.chi2", sbs_tr_chi2, &b_sbs_tr_chi2);
   fChain->SetBranchAddress("Ndata.sbs.tr.d_ph", &Ndata_sbs_tr_d_ph, &b_Ndata_sbs_tr_d_ph);
   fChain->SetBranchAddress("sbs.tr.d_ph", sbs_tr_d_ph, &b_sbs_tr_d_ph);
   fChain->SetBranchAddress("Ndata.sbs.tr.d_th", &Ndata_sbs_tr_d_th, &b_Ndata_sbs_tr_d_th);
   fChain->SetBranchAddress("sbs.tr.d_th", sbs_tr_d_th, &b_sbs_tr_d_th);
   fChain->SetBranchAddress("Ndata.sbs.tr.d_x", &Ndata_sbs_tr_d_x, &b_Ndata_sbs_tr_d_x);
   fChain->SetBranchAddress("sbs.tr.d_x", sbs_tr_d_x, &b_sbs_tr_d_x);
   fChain->SetBranchAddress("Ndata.sbs.tr.d_y", &Ndata_sbs_tr_d_y, &b_Ndata_sbs_tr_d_y);
   fChain->SetBranchAddress("sbs.tr.d_y", sbs_tr_d_y, &b_sbs_tr_d_y);
   fChain->SetBranchAddress("Ndata.sbs.tr.dbeta", &Ndata_sbs_tr_dbeta, &b_Ndata_sbs_tr_dbeta);
   fChain->SetBranchAddress("sbs.tr.dbeta", sbs_tr_dbeta, &b_sbs_tr_dbeta);
   fChain->SetBranchAddress("Ndata.sbs.tr.dtime", &Ndata_sbs_tr_dtime, &b_Ndata_sbs_tr_dtime);
   fChain->SetBranchAddress("sbs.tr.dtime", sbs_tr_dtime, &b_sbs_tr_dtime);
   fChain->SetBranchAddress("Ndata.sbs.tr.flag", &Ndata_sbs_tr_flag, &b_Ndata_sbs_tr_flag);
   fChain->SetBranchAddress("sbs.tr.flag", sbs_tr_flag, &b_sbs_tr_flag);
   fChain->SetBranchAddress("Ndata.sbs.tr.ndof", &Ndata_sbs_tr_ndof, &b_Ndata_sbs_tr_ndof);
   fChain->SetBranchAddress("sbs.tr.ndof", sbs_tr_ndof, &b_sbs_tr_ndof);
   fChain->SetBranchAddress("Ndata.sbs.tr.p", &Ndata_sbs_tr_p, &b_Ndata_sbs_tr_p);
   fChain->SetBranchAddress("sbs.tr.p", sbs_tr_p, &b_sbs_tr_p);
   fChain->SetBranchAddress("Ndata.sbs.tr.pathl", &Ndata_sbs_tr_pathl, &b_Ndata_sbs_tr_pathl);
   fChain->SetBranchAddress("sbs.tr.pathl", sbs_tr_pathl, &b_sbs_tr_pathl);
   fChain->SetBranchAddress("Ndata.sbs.tr.ph", &Ndata_sbs_tr_ph, &b_Ndata_sbs_tr_ph);
   fChain->SetBranchAddress("sbs.tr.ph", sbs_tr_ph, &b_sbs_tr_ph);
   fChain->SetBranchAddress("Ndata.sbs.tr.px", &Ndata_sbs_tr_px, &b_Ndata_sbs_tr_px);
   fChain->SetBranchAddress("sbs.tr.px", sbs_tr_px, &b_sbs_tr_px);
   fChain->SetBranchAddress("Ndata.sbs.tr.py", &Ndata_sbs_tr_py, &b_Ndata_sbs_tr_py);
   fChain->SetBranchAddress("sbs.tr.py", sbs_tr_py, &b_sbs_tr_py);
   fChain->SetBranchAddress("Ndata.sbs.tr.pz", &Ndata_sbs_tr_pz, &b_Ndata_sbs_tr_pz);
   fChain->SetBranchAddress("sbs.tr.pz", sbs_tr_pz, &b_sbs_tr_pz);
   fChain->SetBranchAddress("Ndata.sbs.tr.r_ph", &Ndata_sbs_tr_r_ph, &b_Ndata_sbs_tr_r_ph);
   fChain->SetBranchAddress("sbs.tr.r_ph", sbs_tr_r_ph, &b_sbs_tr_r_ph);
   fChain->SetBranchAddress("Ndata.sbs.tr.r_th", &Ndata_sbs_tr_r_th, &b_Ndata_sbs_tr_r_th);
   fChain->SetBranchAddress("sbs.tr.r_th", sbs_tr_r_th, &b_sbs_tr_r_th);
   fChain->SetBranchAddress("Ndata.sbs.tr.r_x", &Ndata_sbs_tr_r_x, &b_Ndata_sbs_tr_r_x);
   fChain->SetBranchAddress("sbs.tr.r_x", sbs_tr_r_x, &b_sbs_tr_r_x);
   fChain->SetBranchAddress("Ndata.sbs.tr.r_y", &Ndata_sbs_tr_r_y, &b_Ndata_sbs_tr_r_y);
   fChain->SetBranchAddress("sbs.tr.r_y", sbs_tr_r_y, &b_sbs_tr_r_y);
   fChain->SetBranchAddress("Ndata.sbs.tr.tg_dp", &Ndata_sbs_tr_tg_dp, &b_Ndata_sbs_tr_tg_dp);
   fChain->SetBranchAddress("sbs.tr.tg_dp", sbs_tr_tg_dp, &b_sbs_tr_tg_dp);
   fChain->SetBranchAddress("Ndata.sbs.tr.tg_ph", &Ndata_sbs_tr_tg_ph, &b_Ndata_sbs_tr_tg_ph);
   fChain->SetBranchAddress("sbs.tr.tg_ph", sbs_tr_tg_ph, &b_sbs_tr_tg_ph);
   fChain->SetBranchAddress("Ndata.sbs.tr.tg_th", &Ndata_sbs_tr_tg_th, &b_Ndata_sbs_tr_tg_th);
   fChain->SetBranchAddress("sbs.tr.tg_th", sbs_tr_tg_th, &b_sbs_tr_tg_th);
   fChain->SetBranchAddress("Ndata.sbs.tr.tg_y", &Ndata_sbs_tr_tg_y, &b_Ndata_sbs_tr_tg_y);
   fChain->SetBranchAddress("sbs.tr.tg_y", sbs_tr_tg_y, &b_sbs_tr_tg_y);
   fChain->SetBranchAddress("Ndata.sbs.tr.th", &Ndata_sbs_tr_th, &b_Ndata_sbs_tr_th);
   fChain->SetBranchAddress("sbs.tr.th", sbs_tr_th, &b_sbs_tr_th);
   fChain->SetBranchAddress("Ndata.sbs.tr.time", &Ndata_sbs_tr_time, &b_Ndata_sbs_tr_time);
   fChain->SetBranchAddress("sbs.tr.time", sbs_tr_time, &b_sbs_tr_time);
   fChain->SetBranchAddress("Ndata.sbs.tr.vx", &Ndata_sbs_tr_vx, &b_Ndata_sbs_tr_vx);
   fChain->SetBranchAddress("sbs.tr.vx", sbs_tr_vx, &b_sbs_tr_vx);
   fChain->SetBranchAddress("Ndata.sbs.tr.vy", &Ndata_sbs_tr_vy, &b_Ndata_sbs_tr_vy);
   fChain->SetBranchAddress("sbs.tr.vy", sbs_tr_vy, &b_sbs_tr_vy);
   fChain->SetBranchAddress("Ndata.sbs.tr.vz", &Ndata_sbs_tr_vz, &b_Ndata_sbs_tr_vz);
   fChain->SetBranchAddress("sbs.tr.vz", sbs_tr_vz, &b_sbs_tr_vz);
   fChain->SetBranchAddress("Ndata.sbs.tr.x", &Ndata_sbs_tr_x, &b_Ndata_sbs_tr_x);
   fChain->SetBranchAddress("sbs.tr.x", sbs_tr_x, &b_sbs_tr_x);
   fChain->SetBranchAddress("Ndata.sbs.tr.y", &Ndata_sbs_tr_y, &b_Ndata_sbs_tr_y);
   fChain->SetBranchAddress("sbs.tr.y", sbs_tr_y, &b_sbs_tr_y);
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
   fChain->SetBranchAddress("FTROI.ECAL_energy", &FTROI_ECAL_energy, &b_FTROI_ECAL_energy);
   fChain->SetBranchAddress("FTROI.Eprime", &FTROI_Eprime, &b_FTROI_Eprime);
   fChain->SetBranchAddress("FTROI.ephi", &FTROI_ephi, &b_FTROI_ephi);
   fChain->SetBranchAddress("FTROI.etheta", &FTROI_etheta, &b_FTROI_etheta);
   fChain->SetBranchAddress("FTROI.good", &FTROI_good, &b_FTROI_good);
   fChain->SetBranchAddress("FTROI.pp", &FTROI_pp, &b_FTROI_pp);
   fChain->SetBranchAddress("FTROI.pphi", &FTROI_pphi, &b_FTROI_pphi);
   fChain->SetBranchAddress("FTROI.ptheta", &FTROI_ptheta, &b_FTROI_ptheta);
   fChain->SetBranchAddress("FTROI.xECAL_global", &FTROI_xECAL_global, &b_FTROI_xECAL_global);
   fChain->SetBranchAddress("FTROI.yECAL_global", &FTROI_yECAL_global, &b_FTROI_yECAL_global);
   fChain->SetBranchAddress("FTROI.zECAL_global", &FTROI_zECAL_global, &b_FTROI_zECAL_global);
   fChain->SetBranchAddress("MC.mc_ep", &MC_mc_ep, &b_MC_mc_ep);
   fChain->SetBranchAddress("MC.mc_epx", &MC_mc_epx, &b_MC_mc_epx);
   fChain->SetBranchAddress("MC.mc_epy", &MC_mc_epy, &b_MC_mc_epy);
   fChain->SetBranchAddress("MC.mc_epz", &MC_mc_epz, &b_MC_mc_epz);
   fChain->SetBranchAddress("MC.mc_fnucl", &MC_mc_fnucl, &b_MC_mc_fnucl);
   fChain->SetBranchAddress("MC.mc_np", &MC_mc_np, &b_MC_mc_np);
   fChain->SetBranchAddress("MC.mc_npx", &MC_mc_npx, &b_MC_mc_npx);
   fChain->SetBranchAddress("MC.mc_npy", &MC_mc_npy, &b_MC_mc_npy);
   fChain->SetBranchAddress("MC.mc_npz", &MC_mc_npz, &b_MC_mc_npz);
   fChain->SetBranchAddress("MC.mc_nucl", &MC_mc_nucl, &b_MC_mc_nucl);
   fChain->SetBranchAddress("MC.mc_omega", &MC_mc_omega, &b_MC_mc_omega);
   fChain->SetBranchAddress("MC.mc_sigma", &MC_mc_sigma, &b_MC_mc_sigma);
   fChain->SetBranchAddress("MC.mc_vx", &MC_mc_vx, &b_MC_mc_vx);
   fChain->SetBranchAddress("MC.mc_vy", &MC_mc_vy, &b_MC_mc_vy);
   fChain->SetBranchAddress("MC.mc_vz", &MC_mc_vz, &b_MC_mc_vz);
   fChain->SetBranchAddress("earm.ecal.againblk", &earm_ecal_againblk, &b_earm_ecal_againblk);
   fChain->SetBranchAddress("earm.ecal.atimeblk", &earm_ecal_atimeblk, &b_earm_ecal_atimeblk);
   fChain->SetBranchAddress("earm.ecal.colblk", &earm_ecal_colblk, &b_earm_ecal_colblk);
   fChain->SetBranchAddress("earm.ecal.e", &earm_ecal_e, &b_earm_ecal_e);
   fChain->SetBranchAddress("earm.ecal.eblk", &earm_ecal_eblk, &b_earm_ecal_eblk);
   fChain->SetBranchAddress("earm.ecal.idblk", &earm_ecal_idblk, &b_earm_ecal_idblk);
   fChain->SetBranchAddress("earm.ecal.index", &earm_ecal_index, &b_earm_ecal_index);
   fChain->SetBranchAddress("earm.ecal.nblk", &earm_ecal_nblk, &b_earm_ecal_nblk);
   fChain->SetBranchAddress("earm.ecal.nclus", &earm_ecal_nclus, &b_earm_ecal_nclus);
   fChain->SetBranchAddress("earm.ecal.ngoodADChits", &earm_ecal_ngoodADChits, &b_earm_ecal_ngoodADChits);
   fChain->SetBranchAddress("earm.ecal.rowblk", &earm_ecal_rowblk, &b_earm_ecal_rowblk);
   fChain->SetBranchAddress("earm.ecal.x", &earm_ecal_x, &b_earm_ecal_x);
   fChain->SetBranchAddress("earm.ecal.y", &earm_ecal_y, &b_earm_ecal_y);
   fChain->SetBranchAddress("earm.tr.n", &earm_tr_n, &b_earm_tr_n);
   fChain->SetBranchAddress("sbs.HCALdir_x", &sbs_HCALdir_x, &b_sbs_HCALdir_x);
   fChain->SetBranchAddress("sbs.HCALdir_y", &sbs_HCALdir_y, &b_sbs_HCALdir_y);
   fChain->SetBranchAddress("sbs.HCALdir_z", &sbs_HCALdir_z, &b_sbs_HCALdir_z);
   fChain->SetBranchAddress("sbs.HCALph_n", &sbs_HCALph_n, &b_sbs_HCALph_n);
   fChain->SetBranchAddress("sbs.HCALth_n", &sbs_HCALth_n, &b_sbs_HCALth_n);
   fChain->SetBranchAddress("sbs.gemFPP.hit.ngoodhits", &sbs_gemFPP_hit_ngoodhits, &b_sbs_gemFPP_hit_ngoodhits);
   fChain->SetBranchAddress("sbs.gemFPP.nlayershit", &sbs_gemFPP_nlayershit, &b_sbs_gemFPP_nlayershit);
   fChain->SetBranchAddress("sbs.gemFPP.nlayershitu", &sbs_gemFPP_nlayershitu, &b_sbs_gemFPP_nlayershitu);
   fChain->SetBranchAddress("sbs.gemFPP.nlayershituv", &sbs_gemFPP_nlayershituv, &b_sbs_gemFPP_nlayershituv);
   fChain->SetBranchAddress("sbs.gemFPP.nlayershitv", &sbs_gemFPP_nlayershitv, &b_sbs_gemFPP_nlayershitv);
   fChain->SetBranchAddress("sbs.gemFPP.track.besttrack", &sbs_gemFPP_track_besttrack, &b_sbs_gemFPP_track_besttrack);
   fChain->SetBranchAddress("sbs.gemFPP.track.ntrack", &sbs_gemFPP_track_ntrack, &b_sbs_gemFPP_track_ntrack);
   fChain->SetBranchAddress("sbs.gemFT.hit.ngoodhits", &sbs_gemFT_hit_ngoodhits, &b_sbs_gemFT_hit_ngoodhits);
   fChain->SetBranchAddress("sbs.gemFT.nlayershit", &sbs_gemFT_nlayershit, &b_sbs_gemFT_nlayershit);
   fChain->SetBranchAddress("sbs.gemFT.nlayershitu", &sbs_gemFT_nlayershitu, &b_sbs_gemFT_nlayershitu);
   fChain->SetBranchAddress("sbs.gemFT.nlayershituv", &sbs_gemFT_nlayershituv, &b_sbs_gemFT_nlayershituv);
   fChain->SetBranchAddress("sbs.gemFT.nlayershitv", &sbs_gemFT_nlayershitv, &b_sbs_gemFT_nlayershitv);
   fChain->SetBranchAddress("sbs.gemFT.track.besttrack", &sbs_gemFT_track_besttrack, &b_sbs_gemFT_track_besttrack);
   fChain->SetBranchAddress("sbs.gemFT.track.ntrack", &sbs_gemFT_track_ntrack, &b_sbs_gemFT_track_ntrack);
   fChain->SetBranchAddress("sbs.hcal.adctime", &sbs_hcal_adctime, &b_sbs_hcal_adctime);
   fChain->SetBranchAddress("sbs.hcal.againblk", &sbs_hcal_againblk, &b_sbs_hcal_againblk);
   fChain->SetBranchAddress("sbs.hcal.atimeblk", &sbs_hcal_atimeblk, &b_sbs_hcal_atimeblk);
   fChain->SetBranchAddress("sbs.hcal.colblk", &sbs_hcal_colblk, &b_sbs_hcal_colblk);
   fChain->SetBranchAddress("sbs.hcal.colblk_goodtdc", &sbs_hcal_colblk_goodtdc, &b_sbs_hcal_colblk_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.e", &sbs_hcal_e, &b_sbs_hcal_e);
   fChain->SetBranchAddress("sbs.hcal.e_goodtdc", &sbs_hcal_e_goodtdc, &b_sbs_hcal_e_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.eblk", &sbs_hcal_eblk, &b_sbs_hcal_eblk);
   fChain->SetBranchAddress("sbs.hcal.eblk_goodtdc", &sbs_hcal_eblk_goodtdc, &b_sbs_hcal_eblk_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.idblk", &sbs_hcal_idblk, &b_sbs_hcal_idblk);
   fChain->SetBranchAddress("sbs.hcal.idblk_goodtdc", &sbs_hcal_idblk_goodtdc, &b_sbs_hcal_idblk_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.index", &sbs_hcal_index, &b_sbs_hcal_index);
   fChain->SetBranchAddress("sbs.hcal.nblk", &sbs_hcal_nblk, &b_sbs_hcal_nblk);
   fChain->SetBranchAddress("sbs.hcal.nblk_goodtdc", &sbs_hcal_nblk_goodtdc, &b_sbs_hcal_nblk_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus, &b_sbs_hcal_nclus);
   fChain->SetBranchAddress("sbs.hcal.ngoodADChits", &sbs_hcal_ngoodADChits, &b_sbs_hcal_ngoodADChits);
   fChain->SetBranchAddress("sbs.hcal.rowblk", &sbs_hcal_rowblk, &b_sbs_hcal_rowblk);
   fChain->SetBranchAddress("sbs.hcal.rowblk_goodtdc", &sbs_hcal_rowblk_goodtdc, &b_sbs_hcal_rowblk_goodtdc);
   fChain->SetBranchAddress("sbs.hcal.tdctime", &sbs_hcal_tdctime, &b_sbs_hcal_tdctime);
   fChain->SetBranchAddress("sbs.hcal.tdctime_tw", &sbs_hcal_tdctime_tw, &b_sbs_hcal_tdctime_tw);
   fChain->SetBranchAddress("sbs.hcal.tdctimeblk", &sbs_hcal_tdctimeblk, &b_sbs_hcal_tdctimeblk);
   fChain->SetBranchAddress("sbs.hcal.tdctimeblk_tw", &sbs_hcal_tdctimeblk_tw, &b_sbs_hcal_tdctimeblk_tw);
   fChain->SetBranchAddress("sbs.hcal.x", &sbs_hcal_x, &b_sbs_hcal_x);
   fChain->SetBranchAddress("sbs.hcal.y", &sbs_hcal_y, &b_sbs_hcal_y);
   fChain->SetBranchAddress("sbs.tr.n", &sbs_tr_n, &b_sbs_tr_n);
   fChain->SetBranchAddress("singletrack", &singletrack, &b_singletrack);
   fChain->SetBranchAddress("anytrack", &anytrack, &b_anytrack);
//    fChain->SetBranchAddress("singletrack", &singletrack, &b_singletrack);
//    fChain->SetBranchAddress("anytrack", &anytrack, &b_anytrack);
   fChain->SetBranchAddress("fEvtHdr.fEvtTime", &fEvtHdr_fEvtTime, &b_Event_Branch_fEvtHdr_fEvtTime);
   fChain->SetBranchAddress("fEvtHdr.fEvtNum", &fEvtHdr_fEvtNum, &b_Event_Branch_fEvtHdr_fEvtNum);
   fChain->SetBranchAddress("fEvtHdr.fEvtType", &fEvtHdr_fEvtType, &b_Event_Branch_fEvtHdr_fEvtType);
   fChain->SetBranchAddress("fEvtHdr.fEvtLen", &fEvtHdr_fEvtLen, &b_Event_Branch_fEvtHdr_fEvtLen);
   fChain->SetBranchAddress("fEvtHdr.fHelicity", &fEvtHdr_fHelicity, &b_Event_Branch_fEvtHdr_fHelicity);
   fChain->SetBranchAddress("fEvtHdr.fTrigBits", &fEvtHdr_fTrigBits, &b_Event_Branch_fEvtHdr_fTrigBits);
   fChain->SetBranchAddress("fEvtHdr.fRun", &fEvtHdr_fRun, &b_Event_Branch_fEvtHdr_fRun);
   Notify();
}

Bool_t gep_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gep_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gep_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gep_tree_cxx
