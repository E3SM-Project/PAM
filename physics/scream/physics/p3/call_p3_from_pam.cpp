
#include "call_p3_from_pam.h"
#include "p3_functions.hpp"
#include "ArrayIR_kokkos.h"

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
                             ArrayIR<double,2,IRMemDevice> input_col_location       ) {
    // Create some convenient types using ekat and shoc Functions
    typedef double                                          Scalar             ;
    typedef ekat::DefaultDevice                             Device             ;
    typedef typename scream::p3::Functions<double,Device>   P3                 ;
    typedef typename P3::Spack                              Spack              ;
    typedef typename P3::P3PrognosticState                  P3PrognosticState  ;
    typedef typename P3::P3DiagnosticInputs                 P3DiagnosticInputs ;
    typedef typename P3::P3DiagnosticOutputs                P3DiagnosticOutputs;
    typedef typename P3::P3Infrastructure                   P3Infrastructure   ;
    typedef typename P3::P3HistoryOnly                      P3HistoryOnly      ;
    typedef typename P3::P3LookupTables                     P3LookupTables     ;
    typedef typename ekat::WorkspaceManager<Spack,Device>   WorkspaceManager   ;

    int nlev = kte - kts + 1;
    int ncol = ite - its + 1;

    static_assert( SCREAM_SMALL_PACK_SIZE == 1 ,
                   "ERROR: PAM's p3 integration isn't setup to deal with SCREAM_SMALL_PACK_SIZE > 1" );

    // Transform the PAM ArrayIR metadata and data pointers into Kokkos View objects of the appropriate dimensions
    auto qc                  = arrayIR_to_kokkos_view( input_qc                 );
    auto nc                  = arrayIR_to_kokkos_view( input_nc                 );
    auto qr                  = arrayIR_to_kokkos_view( input_qr                 );
    auto nr                  = arrayIR_to_kokkos_view( input_nr                 );
    auto theta               = arrayIR_to_kokkos_view( input_theta              );
    auto qv                  = arrayIR_to_kokkos_view( input_qv                 );
    auto qi                  = arrayIR_to_kokkos_view( input_qi                 );
    auto qm                  = arrayIR_to_kokkos_view( input_qm                 );
    auto ni                  = arrayIR_to_kokkos_view( input_ni                 );
    auto bm                  = arrayIR_to_kokkos_view( input_bm                 );
    auto pressure            = arrayIR_to_kokkos_view( input_pressure           );
    auto dz                  = arrayIR_to_kokkos_view( input_dz                 );
    auto nc_nuceat_tend      = arrayIR_to_kokkos_view( input_nc_nuceat_tend     );
    auto nccn_prescribed     = arrayIR_to_kokkos_view( input_nccn_prescribed    );
    auto ni_activated        = arrayIR_to_kokkos_view( input_ni_activated       );
    auto inv_qc_relvar       = arrayIR_to_kokkos_view( input_inv_qc_relvar      );
    auto precip_liq_surf     = arrayIR_to_kokkos_view( input_precip_liq_surf    );
    auto precip_ice_surf     = arrayIR_to_kokkos_view( input_precip_ice_surf    );
    auto diag_eff_radius_qc  = arrayIR_to_kokkos_view( input_diag_eff_radius_qc );
    auto diag_eff_radius_qi  = arrayIR_to_kokkos_view( input_diag_eff_radius_qi );
    auto bulk_qi             = arrayIR_to_kokkos_view( input_bulk_qi            );
    auto dpres               = arrayIR_to_kokkos_view( input_dpres              );
    auto inv_exner           = arrayIR_to_kokkos_view( input_inv_exner          );
    auto qv2qi_depos_tend    = arrayIR_to_kokkos_view( input_qv2qi_depos_tend   );
    auto precip_total_tend   = arrayIR_to_kokkos_view( input_precip_total_tend  );
    auto nevapr              = arrayIR_to_kokkos_view( input_nevapr             );
    auto qr_evap_tend        = arrayIR_to_kokkos_view( input_qr_evap_tend       );
    auto precip_liq_flux     = arrayIR_to_kokkos_view( input_precip_liq_flux    );
    auto precip_ice_flux     = arrayIR_to_kokkos_view( input_precip_ice_flux    );
    auto cld_frac_r          = arrayIR_to_kokkos_view( input_cld_frac_r         );
    auto cld_frac_l          = arrayIR_to_kokkos_view( input_cld_frac_l         );
    auto cld_frac_i          = arrayIR_to_kokkos_view( input_cld_frac_i         );
    auto p3_tend_out         = arrayIR_to_kokkos_view( input_p3_tend_out        );
    auto mu_c                = arrayIR_to_kokkos_view( input_mu_c               );
    auto lamc                = arrayIR_to_kokkos_view( input_lamc               );
    auto liq_ice_exchange    = arrayIR_to_kokkos_view( input_liq_ice_exchange   );
    auto vap_liq_exchange    = arrayIR_to_kokkos_view( input_vap_liq_exchange   );
    auto vap_ice_exchange    = arrayIR_to_kokkos_view( input_vap_ice_exchange   );
    auto qv_prev             = arrayIR_to_kokkos_view( input_qv_prev            );
    auto t_prev              = arrayIR_to_kokkos_view( input_t_prev             );
    auto col_location        = arrayIR_to_kokkos_view( input_col_location       );

    typedef typename P3::view_1d_table       view_1d_table     ;
    typedef typename P3::view_2d_table       view_2d_table     ;
    typedef typename P3::view_ice_table      view_ice_table    ;
    typedef typename P3::view_collect_table  view_collect_table;
    typedef typename P3::view_dnu_table      view_dnu_table    ;
    view_1d_table      mu_r_table_vals;
    view_2d_table      vn_table_vals, vm_table_vals, revap_table_vals;
    view_ice_table     ice_table_vals;
    view_collect_table collect_table_vals;
    view_dnu_table     dnu_table_vals;
    P3::init_kokkos_ice_lookup_tables(ice_table_vals, collect_table_vals);
    P3::init_kokkos_tables(vn_table_vals, vm_table_vals, revap_table_vals, mu_r_table_vals, dnu_table_vals);

    P3LookupTables lookup_tables{mu_r_table_vals, vn_table_vals, vm_table_vals, revap_table_vals,
                                 ice_table_vals, collect_table_vals, dnu_table_vals};

    const int nlev_pack = ekat::npack<Spack>(nlev);
    const auto policy = ekat::ExeSpaceUtils<ekat::KokkosTypes<Device>::ExeSpace>::get_default_team_policy(ncol, nlev_pack);
    ekat::WorkspaceManager<Spack,Device> workspace_mgr(nlev_pack, 52, policy);
  }

}


