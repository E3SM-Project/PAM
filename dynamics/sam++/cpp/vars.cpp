
#include "vars.h"



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
                            real *dt3_p, real *dz_p) {
  
  u                = umgReal4d(u_p               ,ncrms,dimx_u,dimy_u,nzm  ); 
  v                = umgReal4d(v_p               ,ncrms,dimx_v,dimy_v,nzm  ); 
  w                = umgReal4d(w_p               ,ncrms,dimx_w,dimy_w,nz   ); 
  t                = umgReal4d(t_p               ,ncrms,dimx_s,dimy_s,nzm  ); 
  p                = umgReal4d(p_p               ,ncrms,nxp1,dimy_p, nzm   ); 
  tabs             = umgReal4d(tabs_p            ,ncrms,nx, ny, nzm        ); 
  qv               = umgReal4d(qv_p              ,ncrms,nx, ny, nzm        ); 
  qcl              = umgReal4d(qcl_p             ,ncrms,nx, ny, nzm        ); 
  qpl              = umgReal4d(qpl_p             ,ncrms,nx, ny, nzm        ); 
  qci              = umgReal4d(qci_p             ,ncrms,nx, ny, nzm        ); 
  qpi              = umgReal4d(qpi_p             ,ncrms,nx, ny, nzm        ); 
  tke2             = umgReal4d(tke2_p            ,ncrms,dimx_s,dimy_s, nzm ); 
  tk2              = umgReal4d(tk2_p             ,ncrms,nxp2,dimy_tk2, nzm ); 
  dudt             = umgReal5d(dudt_p            ,ncrms,nxp1, ny, nzm, 3   ); 
  dvdt             = umgReal5d(dvdt_p            ,ncrms,nx, nyp1, nzm, 3   ); 
  dwdt             = umgReal5d(dwdt_p            ,ncrms,nx, ny  , nz,  3   ); 
  misc             = umgReal4d(misc_p            ,ncrms,nx, ny, nz         ); 
  fluxbu           = umgReal3d(fluxbu_p          ,ncrms,nx,ny              ); 
  fluxbv           = umgReal3d(fluxbv_p          ,ncrms,nx,ny              ); 
  fluxbt           = umgReal3d(fluxbt_p          ,ncrms,nx,ny              ); 
  fluxbq           = umgReal3d(fluxbq_p          ,ncrms,nx,ny              ); 
  fluxtu           = umgReal3d(fluxtu_p          ,ncrms,nx,ny              ); 
  fluxtv           = umgReal3d(fluxtv_p          ,ncrms,nx,ny              ); 
  fluxtt           = umgReal3d(fluxtt_p          ,ncrms,nx,ny              ); 
  fluxtq           = umgReal3d(fluxtq_p          ,ncrms,nx,ny              ); 
  fzero            = umgReal3d(fzero_p           ,ncrms,nx,ny              ); 
  precsfc          = umgReal3d(precsfc_p         ,ncrms,nx,ny              ); 
  precssfc         = umgReal3d(precssfc_p        ,ncrms,nx,ny              ); 
  t0               = umgReal2d(t0_p              ,ncrms,nzm                ); 
  q0               = umgReal2d(q0_p              ,ncrms,nzm                ); 
  qv0              = umgReal2d(qv0_p             ,ncrms,nzm                ); 
  tabs0            = umgReal2d(tabs0_p           ,ncrms,nzm                ); 
  tv0              = umgReal2d(tv0_p             ,ncrms,nzm                ); 
  u0               = umgReal2d(u0_p              ,ncrms,nzm                ); 
  v0               = umgReal2d(v0_p              ,ncrms,nzm                ); 
  tg0              = umgReal2d(tg0_p             ,ncrms,nzm                ); 
  qg0              = umgReal2d(qg0_p             ,ncrms,nzm                ); 
  ug0              = umgReal2d(ug0_p             ,ncrms,nzm                ); 
  vg0              = umgReal2d(vg0_p             ,ncrms,nzm                ); 
  p0               = umgReal2d(p0_p              ,ncrms,nzm                ); 
  tke0             = umgReal2d(tke0_p            ,ncrms,nzm                ); 
  t01              = umgReal2d(t01_p             ,ncrms,nzm                ); 
  q01              = umgReal2d(q01_p             ,ncrms,nzm                ); 
  qp0              = umgReal2d(qp0_p             ,ncrms,nzm                ); 
  qn0              = umgReal2d(qn0_p             ,ncrms,nzm                ); 
  prespot          = umgReal2d(prespot_p         ,ncrms,nzm                ); 
  rho              = umgReal2d(rho_p             ,ncrms,nzm                ); 
  rhow             = umgReal2d(rhow_p            ,ncrms,nz                 ); 
  bet              = umgReal2d(bet_p             ,ncrms,nzm                ); 
  gamaz            = umgReal2d(gamaz_p           ,ncrms,nzm                ); 
  wsub             = umgReal2d(wsub_p            ,ncrms,nz                 ); 
  qtend            = umgReal2d(qtend_p           ,ncrms,nzm                ); 
  ttend            = umgReal2d(ttend_p           ,ncrms,nzm                ); 
  utend            = umgReal2d(utend_p           ,ncrms,nzm                ); 
  vtend            = umgReal2d(vtend_p           ,ncrms,nzm                ); 
  sstxy            = umgReal3d(sstxy_p           ,ncrms,nxp1,dimy_sstxy    ); 
  fcory            = umgReal2d(fcory_p           ,ncrms,nyp1               ); 
  fcorzy           = umgReal2d(fcorzy_p          ,ncrms,ny                 ); 
  latitude         = umgReal3d(latitude_p        ,ncrms,nx,ny              ); 
  longitude        = umgReal3d(longitude_p       ,ncrms,nx,ny              ); 
  prec_xy          = umgReal3d(prec_xy_p         ,ncrms,nx,ny              ); 
  pw_xy            = umgReal3d(pw_xy_p           ,ncrms,nx,ny              ); 
  cw_xy            = umgReal3d(cw_xy_p           ,ncrms,nx,ny              ); 
  iw_xy            = umgReal3d(iw_xy_p           ,ncrms,nx,ny              ); 
  cld_xy           = umgReal3d(cld_xy_p          ,ncrms,nx,ny              ); 
  u200_xy          = umgReal3d(u200_xy_p         ,ncrms,nx,ny              ); 
  usfc_xy          = umgReal3d(usfc_xy_p         ,ncrms,nx,ny              ); 
  v200_xy          = umgReal3d(v200_xy_p         ,ncrms,nx,ny              ); 
  vsfc_xy          = umgReal3d(vsfc_xy_p         ,ncrms,nx,ny              ); 
  w500_xy          = umgReal3d(w500_xy_p         ,ncrms,nx,ny              ); 
  w_max            = umgReal1d(w_max_p           ,ncrms,nz                 ); 
  u_max            = umgReal1d(u_max_p           ,ncrms,nz                 ); 
  twsb             = umgReal2d(twsb_p            ,ncrms,nz                 ); 
  precflux         = umgReal2d(precflux_p        ,ncrms,nz                 ); 
  uwle             = umgReal2d(uwle_p            ,ncrms,nz                 ); 
  uwsb             = umgReal2d(uwsb_p            ,ncrms,nz                 ); 
  vwle             = umgReal2d(vwle_p            ,ncrms,nz                 ); 
  vwsb             = umgReal2d(vwsb_p            ,ncrms,nz                 ); 
  tkelediss        = umgReal2d(tkelediss_p       ,ncrms,nz                 ); 
  tdiff            = umgReal2d(tdiff_p           ,ncrms,nz                 ); 
  tlat             = umgReal2d(tlat_p            ,ncrms,nz                 ); 
  tlatqi           = umgReal2d(tlatqi_p          ,ncrms,nz                 ); 
  qifall           = umgReal2d(qifall_p          ,ncrms,nx,ny,nzm          ); 
  qpfall           = umgReal2d(qpfall_p          ,ncrms,nx,ny              ); 
  total_water_evap = umgReal1d(total_water_evap_p,ncrms,nx,ny              ); 
  total_water_prec = umgReal1d(total_water_prec_p,ncrms,nx,ny              ); 
  CF3D             = umgReal4d(CF3D_p            ,ncrms,nx,ny              ); 
  u850_xy          = umgReal3d(u850_xy_p         ,ncrms,nx,ny              ); 
  v850_xy          = umgReal3d(v850_xy_p         ,ncrms,nx,ny              ); 
  psfc_xy          = umgReal3d(psfc_xy_p         ,ncrms,nx,ny              ); 
  swvp_xy          = umgReal3d(swvp_xy_p         ,ncrms                    ); 
  cloudtopheight   = umgReal3d(cloudtopheight_p  ,ncrms                    ); 
  echotopheight    = umgReal3d(echotopheight_p   ,ncrms                    ); 
  cloudtoptemp     = umgReal3d(cloudtoptemp_p    ,ncrms                    ); 
  fcorz            = umgReal1d(fcorz_p           ,ncrms                    ); 
  fcor             = umgReal1d(fcor_p            ,ncrms                    ); 
  longitude0       = umgReal1d(longitude0_p      ,ncrms                    ); 
  latitude0        = umgReal1d(latitude0_p       ,ncrms                    ); 
  z0               = umgReal1d(z0_p              ,ncrms                    ); 
  uhl              = umgReal1d(uhl_p             ,ncrms                    ); 
  vhl              = umgReal1d(vhl_p             ,ncrms                    ); 
  taux0            = umgReal1d(taux0_p           ,ncrms                    ); 
  tauy0            = umgReal1d(tauy0_p           ,ncrms                    ); 
  z                = umgReal2d(z_p               ,ncrms,nz                 ); 
  pres             = umgReal2d(pres_p            ,ncrms,nzm                ); 
  zi               = umgReal2d(zi_p              ,ncrms,nz                 ); 
  presi            = umgReal2d(presi_p           ,ncrms,nz                 ); 
  adz              = umgReal2d(adz_p             ,ncrms,nzm                ); 
  adzw             = umgReal2d(adzw_p            ,ncrms,nz                 ); 
  dt3              = umgReal1d(dt3_p             ,3                        ); 
  dz               = umgReal1d(dz_p              ,ncrms                    ); 
}



umgReal4d u               ;
umgReal4d v               ;
umgReal4d w               ;
umgReal4d t               ;
umgReal4d p               ;
umgReal4d tabs            ;
umgReal4d qv              ;
umgReal4d qcl             ;
umgReal4d qpl             ;
umgReal4d qci             ;
umgReal4d qpi             ;
umgReal4d tke2            ;
umgReal4d tk2             ;
umgReal5d dudt            ;
umgReal5d dvdt            ;
umgReal5d dwdt            ;
umgReal4d misc            ;
umgReal3d fluxbu          ;
umgReal3d fluxbv          ;
umgReal3d fluxbt          ;
umgReal3d fluxbq          ;
umgReal3d fluxtu          ;
umgReal3d fluxtv          ;
umgReal3d fluxtt          ;
umgReal3d fluxtq          ;
umgReal3d fzero           ;
umgReal3d precsfc         ;
umgReal3d precssfc        ;
umgReal2d t0              ;
umgReal2d q0              ;
umgReal2d qv0             ;
umgReal2d tabs0           ;
umgReal2d tv0             ;
umgReal2d u0              ;
umgReal2d v0              ;
umgReal2d tg0             ;
umgReal2d qg0             ;
umgReal2d ug0             ;
umgReal2d vg0             ;
umgReal2d p0              ;
umgReal2d tke0            ;
umgReal2d t01             ;
umgReal2d q01             ;
umgReal2d qp0             ;
umgReal2d qn0             ;
umgReal2d prespot         ;
umgReal2d rho             ;
umgReal2d rhow            ;
umgReal2d bet             ;
umgReal2d gamaz           ;
umgReal2d wsub            ;
umgReal2d qtend           ;
umgReal2d ttend           ;
umgReal2d utend           ;
umgReal2d vtend           ;
umgReal3d sstxy           ;
umgReal2d fcory           ;
umgReal2d fcorzy          ;
umgReal3d latitude        ;
umgReal3d longitude       ;
umgReal3d prec_xy         ;
umgReal3d pw_xy           ;
umgReal3d cw_xy           ;
umgReal3d iw_xy           ;
umgReal3d cld_xy          ;
umgReal3d u200_xy         ;
umgReal3d usfc_xy         ;
umgReal3d v200_xy         ;
umgReal3d vsfc_xy         ;
umgReal3d w500_xy         ;
umgReal1d w_max           ;
umgReal1d u_max           ;
umgReal2d twsb            ;
umgReal2d precflux        ;
umgReal2d uwle            ;
umgReal2d uwsb            ;
umgReal2d vwle            ;
umgReal2d vwsb            ;
umgReal2d tkelediss       ;
umgReal2d tdiff           ;
umgReal2d tlat            ;
umgReal2d tlatqi          ;
umgReal2d qifall          ;
umgReal2d qpfall          ;
umgReal1d total_water_evap;
umgReal1d total_water_prec;
umgReal4d CF3D            ;
umgReal3d u850_xy         ;
umgReal3d v850_xy         ;
umgReal3d psfc_xy         ;
umgReal3d swvp_xy         ;
umgReal3d cloudtopheight  ;
umgReal3d echotopheight   ;
umgReal3d cloudtoptemp    ;

umgReal1d fcorz           ;
umgReal1d fcor            ;
umgReal1d longitude0      ;
umgReal1d latitude0       ;
umgReal1d z0              ;
umgReal1d uhl             ;
umgReal1d vhl             ;
umgReal1d taux0           ;
umgReal1d tauy0           ;

umgReal2d z               ;
umgReal2d pres            ;
umgReal2d zi              ;
umgReal2d presi           ;
umgReal2d adz             ;
umgReal2d adzw            ;
umgReal1d dt3             ;
umgReal1d dz              ;







