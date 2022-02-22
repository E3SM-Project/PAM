
#pragma once

#include "ArrayIR.h"

namespace pam {

  void call_shoc_from_pam( int ncol, int nlev, int nlevi, double dtime, int nadv, int num_qtracers,
                           ArrayIR<double,1,IRMemDevice> shoc_host_dx    ,
                           ArrayIR<double,1,IRMemDevice> shoc_host_dy    ,
                           ArrayIR<double,2,IRMemDevice> shoc_thv        ,
                           ArrayIR<double,2,IRMemDevice> shoc_zt_grid    ,
                           ArrayIR<double,2,IRMemDevice> shoc_zi_grid    ,
                           ArrayIR<double,2,IRMemDevice> shoc_pres       ,
                           ArrayIR<double,2,IRMemDevice> shoc_presi      ,
                           ArrayIR<double,2,IRMemDevice> shoc_pdel       ,
                           ArrayIR<double,1,IRMemDevice> shoc_wthl_sfc   ,
                           ArrayIR<double,1,IRMemDevice> shoc_wqw_sfc    ,
                           ArrayIR<double,1,IRMemDevice> shoc_uw_sfc     ,
                           ArrayIR<double,1,IRMemDevice> shoc_vw_sfc     ,
                           ArrayIR<double,2,IRMemDevice> shoc_wtracer_sfc,
                           ArrayIR<double,2,IRMemDevice> shoc_w_field    ,
                           ArrayIR<double,2,IRMemDevice> shoc_inv_exner  ,
                           ArrayIR<double,1,IRMemDevice> shoc_phis       ,
                           ArrayIR<double,2,IRMemDevice> shoc_host_dse   ,
                           ArrayIR<double,2,IRMemDevice> shoc_tke        ,
                           ArrayIR<double,2,IRMemDevice> shoc_thetal     ,
                           ArrayIR<double,2,IRMemDevice> shoc_qw         ,
                           ArrayIR<double,2,IRMemDevice> shoc_u_wind     ,
                           ArrayIR<double,2,IRMemDevice> shoc_v_wind     ,
                           ArrayIR<double,3,IRMemDevice> shoc_qtracers   ,
                           ArrayIR<double,2,IRMemDevice> shoc_wthv_sec   ,
                           ArrayIR<double,2,IRMemDevice> shoc_tkh        ,
                           ArrayIR<double,2,IRMemDevice> shoc_tk         ,
                           ArrayIR<double,2,IRMemDevice> shoc_ql         ,
                           ArrayIR<double,2,IRMemDevice> shoc_cldfrac    ,
                           ArrayIR<double,1,IRMemDevice> shoc_pblh       ,
                           ArrayIR<double,2,IRMemDevice> shoc_mix        ,
                           ArrayIR<double,2,IRMemDevice> shoc_isotropy   ,
                           ArrayIR<double,2,IRMemDevice> shoc_w_sec      ,
                           ArrayIR<double,2,IRMemDevice> shoc_thl_sec    ,
                           ArrayIR<double,2,IRMemDevice> shoc_qw_sec     ,
                           ArrayIR<double,2,IRMemDevice> shoc_qwthl_sec  ,
                           ArrayIR<double,2,IRMemDevice> shoc_wthl_sec   ,
                           ArrayIR<double,2,IRMemDevice> shoc_wqw_sec    ,
                           ArrayIR<double,2,IRMemDevice> shoc_wtke_sec   ,
                           ArrayIR<double,2,IRMemDevice> shoc_uw_sec     ,
                           ArrayIR<double,2,IRMemDevice> shoc_vw_sec     ,
                           ArrayIR<double,2,IRMemDevice> shoc_w3         ,
                           ArrayIR<double,2,IRMemDevice> shoc_wqls_sec   ,
                           ArrayIR<double,2,IRMemDevice> shoc_brunt      ,
                           ArrayIR<double,2,IRMemDevice> shoc_ql2        );

}


