
#pragma once

#include "ArrayIR.h"

namespace pam {

  void call_p3_main_from_pam(double dt , int it , int its , int ite , int kts , int kte , bool do_predict_nc ,
                             bool do_prescribed_CCN , double elapsed_s ,
                             ArrayIR<double,2,IRMemDevice> input_qc                 ,
                             ArrayIR<double,2,IRMemDevice> input_nc                 ,
                             ArrayIR<double,2,IRMemDevice> input_qr                 ,
                             ArrayIR<double,2,IRMemDevice> input_nr                 ,
                             ArrayIR<double,2,IRMemDevice> input_theta              ,
                             ArrayIR<double,2,IRMemDevice> input_qv                 ,
                             ArrayIR<double,2,IRMemDevice> input_qi                 ,
                             ArrayIR<double,2,IRMemDevice> input_qm                 ,
                             ArrayIR<double,2,IRMemDevice> input_ni                 ,
                             ArrayIR<double,2,IRMemDevice> input_bm                 ,
                             ArrayIR<double,2,IRMemDevice> input_pressure           ,
                             ArrayIR<double,2,IRMemDevice> input_dz                 ,
                             ArrayIR<double,2,IRMemDevice> input_nc_nuceat_tend     ,
                             ArrayIR<double,2,IRMemDevice> input_nccn_prescribed    ,
                             ArrayIR<double,2,IRMemDevice> input_ni_activated       ,
                             ArrayIR<double,2,IRMemDevice> input_inv_qc_relvar      ,
                             ArrayIR<double,1,IRMemDevice> input_precip_liq_surf    ,
                             ArrayIR<double,1,IRMemDevice> input_precip_ice_surf    ,
                             ArrayIR<double,2,IRMemDevice> input_diag_eff_radius_qc ,
                             ArrayIR<double,2,IRMemDevice> input_diag_eff_radius_qi ,
                             ArrayIR<double,2,IRMemDevice> input_bulk_qi            ,
                             ArrayIR<double,2,IRMemDevice> input_dpres              ,
                             ArrayIR<double,2,IRMemDevice> input_inv_exner          ,
                             ArrayIR<double,2,IRMemDevice> input_qv2qi_depos_tend   ,
                             ArrayIR<double,2,IRMemDevice> input_precip_total_tend  ,
                             ArrayIR<double,2,IRMemDevice> input_nevapr             ,
                             ArrayIR<double,2,IRMemDevice> input_qr_evap_tend       ,
                             ArrayIR<double,2,IRMemDevice> input_precip_liq_flux    ,
                             ArrayIR<double,2,IRMemDevice> input_precip_ice_flux    ,
                             ArrayIR<double,2,IRMemDevice> input_cld_frac_r         ,
                             ArrayIR<double,2,IRMemDevice> input_cld_frac_l         ,
                             ArrayIR<double,2,IRMemDevice> input_cld_frac_i         ,
                             ArrayIR<double,3,IRMemDevice> input_p3_tend_out        ,
                             ArrayIR<double,2,IRMemDevice> input_mu_c               ,
                             ArrayIR<double,2,IRMemDevice> input_lamc               ,
                             ArrayIR<double,2,IRMemDevice> input_liq_ice_exchange   ,
                             ArrayIR<double,2,IRMemDevice> input_vap_liq_exchange   ,
                             ArrayIR<double,2,IRMemDevice> input_vap_ice_exchange   ,
                             ArrayIR<double,2,IRMemDevice> input_qv_prev            ,
                             ArrayIR<double,2,IRMemDevice> input_t_prev             ,
                             ArrayIR<double,2,IRMemDevice> input_col_location       );

}


