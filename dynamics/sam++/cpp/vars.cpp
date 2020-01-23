
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
                            real *dt3_p, real *dz_p, real *sgs_field_p, real *sgs_field_diag_p, real *grdf_x_p, real *grdf_y_p, real *grdf_z_p,
                            real *tkesbbuoy_p, real *tkesbshear_p, real *tkesbdiss_p) {
  
  u                = umgReal4d( "u               " , u_p                     , nzm , dimy_u     , dimx_u , ncrms ); 
  v                = umgReal4d( "v               " , v_p                     , nzm , dimy_v     , dimx_v , ncrms ); 
  w                = umgReal4d( "w               " , w_p                     , nz  , dimy_w     , dimx_w , ncrms ); 
  t                = umgReal4d( "t               " , t_p                     , nzm , dimy_s     , dimx_s , ncrms ); 
  p                = umgReal4d( "p               " , p_p                     , nzm , dimy_p     , nxp1   , ncrms ); 
  tabs             = umgReal4d( "tabs            " , tabs_p                  , nzm , ny         , nx     , ncrms ); 
  qv               = umgReal4d( "qv              " , qv_p                    , nzm , ny         , nx     , ncrms ); 
  qcl              = umgReal4d( "qcl             " , qcl_p                   , nzm , ny         , nx     , ncrms ); 
  qpl              = umgReal4d( "qpl             " , qpl_p                   , nzm , ny         , nx     , ncrms ); 
  qci              = umgReal4d( "qci             " , qci_p                   , nzm , ny         , nx     , ncrms ); 
  qpi              = umgReal4d( "qpi             " , qpi_p                   , nzm , ny         , nx     , ncrms ); 
  tke2             = umgReal4d( "tke2            " , tke2_p                  , nzm , dimy_s     , dimx_s , ncrms ); 
  tk2              = umgReal4d( "tk2             " , tk2_p                   , nzm , dimy_tk2   , nxp2   , ncrms ); 
  dudt             = umgReal5d( "dudt            " , dudt_p              , 3 , nzm , ny         , nxp1   , ncrms ); 
  dvdt             = umgReal5d( "dvdt            " , dvdt_p              , 3 , nzm , nyp1       , nx     , ncrms ); 
  dwdt             = umgReal5d( "dwdt            " , dwdt_p              , 3 , nz  , ny         , nx     , ncrms ); 
  misc             = umgReal4d( "misc            " , misc_p                  , nz  , ny         , nx     , ncrms ); 
  fluxbu           = umgReal3d( "fluxbu          " , fluxbu_p                      , ny         , nx     , ncrms ); 
  fluxbv           = umgReal3d( "fluxbv          " , fluxbv_p                      , ny         , nx     , ncrms ); 
  fluxbt           = umgReal3d( "fluxbt          " , fluxbt_p                      , ny         , nx     , ncrms ); 
  fluxbq           = umgReal3d( "fluxbq          " , fluxbq_p                      , ny         , nx     , ncrms ); 
  fluxtu           = umgReal3d( "fluxtu          " , fluxtu_p                      , ny         , nx     , ncrms ); 
  fluxtv           = umgReal3d( "fluxtv          " , fluxtv_p                      , ny         , nx     , ncrms ); 
  fluxtt           = umgReal3d( "fluxtt          " , fluxtt_p                      , ny         , nx     , ncrms ); 
  fluxtq           = umgReal3d( "fluxtq          " , fluxtq_p                      , ny         , nx     , ncrms ); 
  fzero            = umgReal3d( "fzero           " , fzero_p                       , ny         , nx     , ncrms ); 
  precsfc          = umgReal3d( "precsfc         " , precsfc_p                     , ny         , nx     , ncrms ); 
  precssfc         = umgReal3d( "precssfc        " , precssfc_p                    , ny         , nx     , ncrms ); 
  t0               = umgReal2d( "t0              " , t0_p                                       , nzm    , ncrms ); 
  q0               = umgReal2d( "q0              " , q0_p                                       , nzm    , ncrms ); 
  qv0              = umgReal2d( "qv0             " , qv0_p                                      , nzm    , ncrms ); 
  tabs0            = umgReal2d( "tabs0           " , tabs0_p                                    , nzm    , ncrms ); 
  tv0              = umgReal2d( "tv0             " , tv0_p                                      , nzm    , ncrms ); 
  u0               = umgReal2d( "u0              " , u0_p                                       , nzm    , ncrms ); 
  v0               = umgReal2d( "v0              " , v0_p                                       , nzm    , ncrms ); 
  tg0              = umgReal2d( "tg0             " , tg0_p                                      , nzm    , ncrms ); 
  qg0              = umgReal2d( "qg0             " , qg0_p                                      , nzm    , ncrms ); 
  ug0              = umgReal2d( "ug0             " , ug0_p                                      , nzm    , ncrms ); 
  vg0              = umgReal2d( "vg0             " , vg0_p                                      , nzm    , ncrms ); 
  p0               = umgReal2d( "p0              " , p0_p                                       , nzm    , ncrms ); 
  tke0             = umgReal2d( "tke0            " , tke0_p                                     , nzm    , ncrms ); 
  t01              = umgReal2d( "t01             " , t01_p                                      , nzm    , ncrms ); 
  q01              = umgReal2d( "q01             " , q01_p                                      , nzm    , ncrms ); 
  qp0              = umgReal2d( "qp0             " , qp0_p                                      , nzm    , ncrms ); 
  qn0              = umgReal2d( "qn0             " , qn0_p                                      , nzm    , ncrms ); 
  prespot          = umgReal2d( "prespot         " , prespot_p                                  , nzm    , ncrms ); 
  rho              = umgReal2d( "rho             " , rho_p                                      , nzm    , ncrms ); 
  rhow             = umgReal2d( "rhow            " , rhow_p                                     , nz     , ncrms ); 
  bet              = umgReal2d( "bet             " , bet_p                                      , nzm    , ncrms ); 
  gamaz            = umgReal2d( "gamaz           " , gamaz_p                                    , nzm    , ncrms ); 
  wsub             = umgReal2d( "wsub            " , wsub_p                                     , nz     , ncrms ); 
  qtend            = umgReal2d( "qtend           " , qtend_p                                    , nzm    , ncrms ); 
  ttend            = umgReal2d( "ttend           " , ttend_p                                    , nzm    , ncrms ); 
  utend            = umgReal2d( "utend           " , utend_p                                    , nzm    , ncrms ); 
  vtend            = umgReal2d( "vtend           " , vtend_p                                    , nzm    , ncrms ); 
  sstxy            = umgReal3d( "sstxy           " , sstxy_p                       , dimy_sstxy , nxp1   , ncrms ); 
  fcory            = umgReal2d( "fcory           " , fcory_p                                    , nyp1   , ncrms ); 
  fcorzy           = umgReal2d( "fcorzy          " , fcorzy_p                                   , ny     , ncrms ); 
  latitude         = umgReal3d( "latitude        " , latitude_p                    , ny         , nx     , ncrms ); 
  longitude        = umgReal3d( "longitude       " , longitude_p                   , ny         , nx     , ncrms ); 
  prec_xy          = umgReal3d( "prec_xy         " , prec_xy_p                     , ny         , nx     , ncrms ); 
  pw_xy            = umgReal3d( "pw_xy           " , pw_xy_p                       , ny         , nx     , ncrms ); 
  cw_xy            = umgReal3d( "cw_xy           " , cw_xy_p                       , ny         , nx     , ncrms ); 
  iw_xy            = umgReal3d( "iw_xy           " , iw_xy_p                       , ny         , nx     , ncrms ); 
  cld_xy           = umgReal3d( "cld_xy          " , cld_xy_p                      , ny         , nx     , ncrms ); 
  u200_xy          = umgReal3d( "u200_xy         " , u200_xy_p                     , ny         , nx     , ncrms ); 
  usfc_xy          = umgReal3d( "usfc_xy         " , usfc_xy_p                     , ny         , nx     , ncrms ); 
  v200_xy          = umgReal3d( "v200_xy         " , v200_xy_p                     , ny         , nx     , ncrms ); 
  vsfc_xy          = umgReal3d( "vsfc_xy         " , vsfc_xy_p                     , ny         , nx     , ncrms ); 
  w500_xy          = umgReal3d( "w500_xy         " , w500_xy_p                     , ny         , nx     , ncrms ); 
  w_max            = umgReal1d( "w_max           " , w_max_p                                    , nz     , ncrms ); 
  u_max            = umgReal1d( "u_max           " , u_max_p                                    , nz     , ncrms ); 
  twsb             = umgReal2d( "twsb            " , twsb_p                                     , nz     , ncrms ); 
  precflux         = umgReal2d( "precflux        " , precflux_p                                 , nz     , ncrms ); 
  uwle             = umgReal2d( "uwle            " , uwle_p                                     , nz     , ncrms ); 
  uwsb             = umgReal2d( "uwsb            " , uwsb_p                                     , nz     , ncrms ); 
  vwle             = umgReal2d( "vwle            " , vwle_p                                     , nz     , ncrms ); 
  vwsb             = umgReal2d( "vwsb            " , vwsb_p                                     , nz     , ncrms ); 
  tkelediss        = umgReal2d( "tkelediss       " , tkelediss_p                                , nz     , ncrms ); 
  tdiff            = umgReal2d( "tdiff           " , tdiff_p                                    , nz     , ncrms ); 
  tlat             = umgReal2d( "tlat            " , tlat_p                                     , nz     , ncrms ); 
  tlatqi           = umgReal2d( "tlatqi          " , tlatqi_p                                   , nz     , ncrms ); 
  qifall           = umgReal2d( "qifall          " , qifall_p                , nzm , ny         , nx     , ncrms ); 
  qpfall           = umgReal2d( "qpfall          " , qpfall_p                      , ny         , nx     , ncrms ); 
  total_water_evap = umgReal1d( "total_water_evap" , total_water_evap_p            , ny         , nx     , ncrms ); 
  total_water_prec = umgReal1d( "total_water_prec" , total_water_prec_p            , ny         , nx     , ncrms ); 
  CF3D             = umgReal4d( "CF3D            " , CF3D_p                        , ny         , nx     , ncrms ); 
  u850_xy          = umgReal3d( "u850_xy         " , u850_xy_p                     , ny         , nx     , ncrms ); 
  v850_xy          = umgReal3d( "v850_xy         " , v850_xy_p                     , ny         , nx     , ncrms ); 
  psfc_xy          = umgReal3d( "psfc_xy         " , psfc_xy_p                     , ny         , nx     , ncrms ); 
  swvp_xy          = umgReal3d( "swvp_xy         " , swvp_xy_p                                           , ncrms ); 
  cloudtopheight   = umgReal3d( "cloudtopheight  " , cloudtopheight_p                                    , ncrms ); 
  echotopheight    = umgReal3d( "echotopheight   " , echotopheight_p                                     , ncrms ); 
  cloudtoptemp     = umgReal3d( "cloudtoptemp    " , cloudtoptemp_p                                      , ncrms ); 
  fcorz            = umgReal1d( "fcorz           " , fcorz_p                                             , ncrms ); 
  fcor             = umgReal1d( "fcor            " , fcor_p                                              , ncrms ); 
  longitude0       = umgReal1d( "longitude0      " , longitude0_p                                        , ncrms ); 
  latitude0        = umgReal1d( "latitude0       " , latitude0_p                                         , ncrms ); 
  z0               = umgReal1d( "z0              " , z0_p                                                , ncrms ); 
  uhl              = umgReal1d( "uhl             " , uhl_p                                               , ncrms ); 
  vhl              = umgReal1d( "vhl             " , vhl_p                                               , ncrms ); 
  taux0            = umgReal1d( "taux0           " , taux0_p                                             , ncrms ); 
  tauy0            = umgReal1d( "tauy0           " , tauy0_p                                             , ncrms ); 
  z                = umgReal2d( "z               " , z_p                                        , nz     , ncrms ); 
  pres             = umgReal2d( "pres            " , pres_p                                     , nzm    , ncrms ); 
  zi               = umgReal2d( "zi              " , zi_p                                       , nz     , ncrms ); 
  presi            = umgReal2d( "presi           " , presi_p                                    , nz     , ncrms ); 
  adz              = umgReal2d( "adz             " , adz_p                                      , nzm    , ncrms ); 
  adzw             = umgReal2d( "adzw            " , adzw_p                                     , nz     , ncrms ); 
  dz               = umgReal1d( "dz              " , dz_p                                                , ncrms ); 
  dt3              = umgReal1d( "dt3             " , dt3_p               , 3                                     ); 
  sgs_field        = umgReal5d( "sgs_field       " , sgs_field_p      , nsgs_fields      , nzm , dimy_s , dimx_s , ncrms );
  sgs_field_diag   = umgReal5d( "sgs_field_diag  " , sgs_field_diag_p , nsgs_fields_diag , nzm , dimy_d , dimx_d , ncrms );
  grdf_x           = umgReal2d( "grdf_x          " , grdf_x_p                            , nzm                   , ncrms );
  grdf_y           = umgReal2d( "grdf_y          " , grdf_y_p                            , nzm                   , ncrms );
  grdf_z           = umgReal2d( "grdf_z          " , grdf_z_p                            , nzm                   , ncrms );
  tkesbbuoy        = umgReal2d( "tkesbbuoy       " , tkesbbuoy_p                         , nz                    , ncrms );
  tkesbshear       = umgReal2d( "tkesbshear      " , tkesbshear_p                        , nz                    , ncrms );
  tkesbdiss        = umgReal2d( "tkesbdiss       " , tkesbdiss_p                         , nz                    , ncrms );
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

umgReal5d sgs_field       ;
umgReal5d sgs_field_diag  ;
umgReal2d grdf_x          ;
umgReal2d grdf_y          ;
umgReal2d grdf_z          ;
umgReal2d tkesbbuoy       ;
umgReal2d tkesbshear      ;
umgReal2d tkesbdiss       ;







