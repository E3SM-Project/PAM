
#ifndef __VARS_H__
#define __VARS_H__

#include "const.h"


//////////////////////////////////////////////////////////////////////////////////
// These arrays use non-1 lower bounds in the Fortran code
// They must be indexed differently in the C++ code
//////////////////////////////////////////////////////////////////////////////////
extern umgReal4d u    ; // Index as u     (k , offy_u    +j , offx_u    +i , icrm )
extern umgReal4d v    ; // Index as v     (k , offy_v    +j , offx_v    +i , icrm )
extern umgReal4d w    ; // Index as w     (k , offy_w    +j , offx_w    +i , icrm )
extern umgReal4d t    ; // Index as t     (k , offy_t    +j , offx_t    +i , icrm )
extern umgReal4d p    ; // Index as p     (k , offy_p    +j , offx_p    +i , icrm )
extern umgReal4d tke2 ; // Index as tke2  (k , offy_tke2 +j , offx_tke2 +i , icrm )
extern umgReal4d tk2  ; // Index as tk2   (k , offy_tk2  +j , offx_tk2  +i , icrm )
extern umgReal3d sstxy; // Index as sstxy (    offy_sstxy+j , offx_sstxy+i , icrm )
extern umgReal2d fcory; // Index as fcory (    offy_fcory+j                , icrm )


extern "C" void wrap_arrays(real *u_p, real *v_p, real *w_p, real *t_p, real *p_p, real *tabs_p, real *qv_p, real *qcl_p, real *qpl_p,
                            real *qci_p, real *qpi_p, real *tke2_p, real *tk2_p, real *dudt_p, real *dvdt_p, real *dwdt_p, real *misc_p, real *fluxbu_p,
                            real *fluxbv_p, real *fluxbt_p, real *fluxbq_p, real *fluxtu_p, real *fluxtv_p, real *fluxtt_p, real *fluxtq_p,
                            real *fzero_p, real *precsfc_p, real *precssfc_p, real *t0_p, real *q0_p, real *qv0_p, real *tabs0_p, real *tv0_p,
                            real *u0_p, real *v0_p, real *tg0_p, real *qg0_p, real *ug0_p, real *vg0_p, real *p0_p, real *tke0_p, real *t01_p, real *q01_p,
                            real *qp0_p, real *qn0_p, real *prespot_p, real *rho_p, real *rhow_p, real *bet_p, real *gamaz_p, real *wsub_p,
                            real *qtend_p, real *ttend_p, real *utend_p, real *vtend_p, real *sstxy_p, real *fcory_p, real *fcorzy_p, real *latitude_p,
                            real *longitude_p, real *prec_xy_p, real *pw_xy_p, real *cw_xy_p, real *iw_xy_p, real *cld_xy_p, real *u200_xy_p,
                            real *usfc_xy_p, real *v200_xy_p, real *vsfc_xy_p, real *w500_xy_p, real *w_max_p, real *u_max_p, real *twsb_p,
                            real *precflux_p, real *uwle_p, real *uwsb_p, real *vwle_p, real *vwsb_p, real *tkelediss_p, real *tdiff_p, real *tlat_p,
                            real *tlatqi_p, real *qifall_p, real *qpfall_p, real *total_water_evap_p, real *total_water_prec_p, real *CF3D_p,
                            real *u850_xy_p, real *v850_xy_p, real *psfc_xy_p, real *swvp_xy_p, real *cloudtopheight_p, real *echotopheight_p,
                            real *cloudtoptemp_p, real *fcorz_p, real *fcor_p, real *longitude0_p, real *latitude0_p, real *z0_p, real *uhl_p,
                            real *vhl_p, real *taux0_p, real *tauy0_p, real *z_p, real *pres_p, real *zi_p, real *presi_p, real *adz_p, real *adzw_p,
                            real *dt3_p, real *dz_p);


extern int  nstep                    ;
extern int  ncycle                   ;
extern int  icycle                   ;
extern int  na, nb, nc               ;
extern real at, bt, ct               ;
extern real dtn                      ;
extern real dtfactor                 ;
extern int  rank                     ;
extern int  ranknn                   ;
extern int  rankss                   ;
extern int  rankee                   ;
extern int  rankww                   ;
extern int  rankne                   ;
extern int  ranknw                   ;
extern int  rankse                   ;
extern int  ranksw                   ;
extern bool dompi                    ;
extern bool masterproc               ;
extern bool dostatis                 ;
extern bool dostatisrad              ;
extern int  nstatis                  ;
extern bool compute_reffc            ;
extern bool compute_reffi            ;
extern bool notopened2D              ;
extern bool notopened3D              ;
extern bool notopenedmom             ;
extern real dx                       ;
extern real dy                       ;
extern bool doconstdz                ;
extern int  nstop                    ;
extern int  nelapse                  ;
extern real dt                       ;
extern real day0                     ;
extern int  nrad                     ;
extern int  nrestart                 ;
extern int  nstat                    ;
extern int  nstatfrq                 ;
extern bool restart_sep              ;
extern int  nrestart_skip            ;
extern bool output_sep               ;
extern bool doisccp                  ;
extern bool domodis                  ;
extern bool domisr                   ;
extern bool dosimfilesout            ;
extern bool doSAMconditionals        ;
extern bool dosatupdnconditionals    ;
extern bool doscamiopdata            ;
extern bool dozero_out_day0          ;
extern int  nsave3Dstart             ;
extern int  nsave3Dend               ;
extern bool save3Dbin                ;
extern bool save3Dsep                ;
extern real qnsave3D                 ;
extern bool dogzip3D                 ;
extern bool rad3Dout                 ;
extern int  nsave2D                  ;
extern int  nsave2Dstart             ;
extern int  nsave2Dend               ;
extern bool save2Dbin                ;
extern bool save2Dsep                ;
extern bool save2Davg                ;
extern bool dogzip2D                 ;
extern int  nstatmom                 ;
extern int  nstatmomstart            ;
extern int  nstatmomend              ;
extern bool savemomsep               ;
extern bool savemombin               ;
extern int  nmovie                   ;
extern int  nmoviestart              ;
extern int  nmovieend                ;
extern bool isInitialized_scamiopdata;
extern bool wgls_holds_omega         ;

extern real epsv            ;
extern bool dosubsidence    ;
extern real ug              ;
extern real vg              ;
extern bool les             ;
extern bool sfc_flx_fxd     ;
extern bool sfc_tau_fxd     ;
extern bool dodamping       ;
extern bool docloud         ;
extern bool docam_sfc_fluxes;
extern bool doprecip        ;
extern bool dosgs           ;
extern bool docoriolis      ;
extern bool dosurface       ;
extern bool dowallx         ;
extern bool dowally         ;
extern bool docolumn        ;
extern bool dotracers       ;
extern bool dosmoke         ;

extern umgReal4d tabs            ;
extern umgReal4d qv              ;
extern umgReal4d qcl             ;
extern umgReal4d qpl             ;
extern umgReal4d qci             ;
extern umgReal4d qpi             ;
extern umgReal5d dudt            ;
extern umgReal5d dvdt            ;
extern umgReal5d dwdt            ;
extern umgReal4d misc            ;
extern umgReal3d fluxbu          ;
extern umgReal3d fluxbv          ;
extern umgReal3d fluxbt          ;
extern umgReal3d fluxbq          ;
extern umgReal3d fluxtu          ;
extern umgReal3d fluxtv          ;
extern umgReal3d fluxtt          ;
extern umgReal3d fluxtq          ;
extern umgReal3d fzero           ;
extern umgReal3d precsfc         ;
extern umgReal3d precssfc        ;
extern umgReal2d t0              ;
extern umgReal2d q0              ;
extern umgReal2d qv0             ;
extern umgReal2d tabs0           ;
extern umgReal2d tv0             ;
extern umgReal2d u0              ;
extern umgReal2d v0              ;
extern umgReal2d tg0             ;
extern umgReal2d qg0             ;
extern umgReal2d ug0             ;
extern umgReal2d vg0             ;
extern umgReal2d p0              ;
extern umgReal2d tke0            ;
extern umgReal2d t01             ;
extern umgReal2d q01             ;
extern umgReal2d qp0             ;
extern umgReal2d qn0             ;
extern umgReal2d prespot         ;
extern umgReal2d rho             ;
extern umgReal2d rhow            ;
extern umgReal2d bet             ;
extern umgReal2d gamaz           ;
extern umgReal2d wsub            ;
extern umgReal2d qtend           ;
extern umgReal2d ttend           ;
extern umgReal2d utend           ;
extern umgReal2d vtend           ;
extern umgReal2d fcorzy          ;
extern umgReal3d latitude        ;
extern umgReal3d longitude       ;
extern umgReal3d prec_xy         ;
extern umgReal3d pw_xy           ;
extern umgReal3d cw_xy           ;
extern umgReal3d iw_xy           ;
extern umgReal3d cld_xy          ;
extern umgReal3d u200_xy         ;
extern umgReal3d usfc_xy         ;
extern umgReal3d v200_xy         ;
extern umgReal3d vsfc_xy         ;
extern umgReal3d w500_xy         ;
extern umgReal1d w_max           ;
extern umgReal1d u_max           ;
extern umgReal2d twsb            ;
extern umgReal2d precflux        ;
extern umgReal2d uwle            ;
extern umgReal2d uwsb            ;
extern umgReal2d vwle            ;
extern umgReal2d vwsb            ;
extern umgReal2d tkelediss       ;
extern umgReal2d tdiff           ;
extern umgReal2d tlat            ;
extern umgReal2d tlatqi          ;
extern umgReal2d qifall          ;
extern umgReal2d qpfall          ;
extern umgReal1d total_water_evap;
extern umgReal1d total_water_prec;
extern umgReal4d CF3D            ;
extern umgReal3d u850_xy         ;
extern umgReal3d v850_xy         ;
extern umgReal3d psfc_xy         ;
extern umgReal3d swvp_xy         ;
extern umgReal3d cloudtopheight  ;
extern umgReal3d echotopheight   ;
extern umgReal3d cloudtoptemp    ;

extern umgReal1d fcorz           ;
extern umgReal1d fcor            ;
extern umgReal1d longitude0      ;
extern umgReal1d latitude0       ;
extern umgReal1d z0              ;
extern umgReal1d uhl             ;
extern umgReal1d vhl             ;
extern umgReal1d taux0           ;
extern umgReal1d tauy0           ;

extern umgReal2d z               ;
extern umgReal2d pres            ;
extern umgReal2d zi              ;
extern umgReal2d presi           ;
extern umgReal2d adz             ;
extern umgReal2d adzw            ;
extern umgReal1d dt3             ;
extern umgReal1d dz              ;

#endif

