
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

    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3PrognosticState_qc                  ("qc                ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3PrognosticState_nc                  ("nc                ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3PrognosticState_qr                  ("qr                ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3PrognosticState_nr                  ("nr                ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3PrognosticState_qi                  ("qi                ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3PrognosticState_qm                  ("qm                ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3PrognosticState_ni                  ("ni                ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3PrognosticState_bm                  ("bm                ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3PrognosticState_qv                  ("qv                ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3PrognosticState_th                  ("th                ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3DiagnosticInputs_nc_nuceat_tend     ("nc_nuceat_tend    ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3DiagnosticInputs_nccn               ("nccn              ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3DiagnosticInputs_ni_activated       ("ni_activated      ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3DiagnosticInputs_inv_qc_relvar      ("inv_qc_relvar     ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3DiagnosticInputs_cld_frac_i         ("cld_frac_i        ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3DiagnosticInputs_cld_frac_l         ("cld_frac_l        ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3DiagnosticInputs_cld_frac_r         ("cld_frac_r        ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3DiagnosticInputs_pres               ("pres              ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3DiagnosticInputs_dz                 ("dz                ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3DiagnosticInputs_dpres              ("dpres             ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3DiagnosticInputs_inv_exner          ("inv_exner         ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3DiagnosticInputs_qv_prev            ("qv_prev           ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3DiagnosticInputs_t_prev             ("t_prev            ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P2DiagnosticOutputs_qv2qi_depos_tend  ("qv2qi_depos_tend  ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_1d<Scalar> P3DiagnosticOutputs_precip_liq_surf   ("precip_liq_surf   ",ncol       );
    typename ekat::KokkosTypes<Device>::view_1d<Scalar> P3DiagnosticOutputs_precip_ice_surf   ("precip_ice_surf   ",ncol       );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3DiagnosticOutputs_diag_eff_radius_qc("diag_eff_radius_qc",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3DiagnosticOutputs_diag_eff_radius_qi("diag_eff_radius_qi",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3DiagnosticOutputs_rho_qi            ("rho_qi            ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3DiagnosticOutputs_precip_liq_flux   ("precip_liq_flux   ",ncol,nlev+1);
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3DiagnosticOutputs_precip_ice_flux   ("precip_ice_flux   ",ncol,nlev+1);
    double                                              P3Infrastructure_dt                                                     ;
    int                                                 P3Infrastructure_it                                                     ;
    int                                                 P3Infrastructure_its                                                    ;
    int                                                 P3Infrastructure_ite                                                    ;
    int                                                 P3Infrastructure_kts                                                    ;
    int                                                 P3Infrastructure_kte                                                    ;
    bool                                                P3Infrastructure_predictNc                                              ;
    bool                                                P3Infrastructure_prescribedCCN                                          ;
    typename ekat::KokkosTypes<Device>::view_2d<Scalar> P3Infrastructure_col_location         ("col_location      ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3HistoryOnly_liq_ice_exchange        ("liq_ice_exchange  ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3HistoryOnly_vap_liq_exchange        ("vap_liq_exchange  ",ncol,nlev  );
    typename ekat::KokkosTypes<Device>::view_2d<Spack > P3HistoryOnly_vap_ice_exchange        ("vap_ice_exchange  ",ncol,nlev  );

    // Initialize the inputs
    P3Infrastructure_dt            = dt               ;
    P3Infrastructure_it            = it               ;
    P3Infrastructure_its           = its              ;
    P3Infrastructure_ite           = ite              ;
    P3Infrastructure_kts           = kts              ;
    P3Infrastructure_kte           = kte              ;
    P3Infrastructure_predictNc     = do_predict_nc    ;
    P3Infrastructure_prescribedCCN = do_prescribed_CCN;
    Kokkos::parallel_for( Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0},{ncol,nlev}) , KOKKOS_LAMBDA (int i, int k) {
      P3PrognosticState_qc             (i,k) = qc             (k,i);
      P3PrognosticState_nc             (i,k) = nc             (k,i);
      P3PrognosticState_qr             (i,k) = qr             (k,i);
      P3PrognosticState_nr             (i,k) = nr             (k,i);
      P3PrognosticState_qi             (i,k) = qi             (k,i);
      P3PrognosticState_qm             (i,k) = qm             (k,i);
      P3PrognosticState_ni             (i,k) = ni             (k,i);
      P3PrognosticState_bm             (i,k) = bm             (k,i);
      P3PrognosticState_qv             (i,k) = qv             (k,i);
      P3PrognosticState_th             (i,k) = theta          (k,i);
      P3DiagnosticInputs_nc_nuceat_tend(i,k) = nc_nuceat_tend (k,i);
      P3DiagnosticInputs_nccn          (i,k) = nccn_prescribed(k,i);
      P3DiagnosticInputs_ni_activated  (i,k) = ni_activated   (k,i);
      P3DiagnosticInputs_inv_qc_relvar (i,k) = inv_qc_relvar  (k,i);
      P3DiagnosticInputs_cld_frac_i    (i,k) = cld_frac_i     (k,i);
      P3DiagnosticInputs_cld_frac_l    (i,k) = cld_frac_l     (k,i);
      P3DiagnosticInputs_cld_frac_r    (i,k) = cld_frac_r     (k,i);
      P3DiagnosticInputs_pres          (i,k) = pressure       (k,i);
      P3DiagnosticInputs_dz            (i,k) = dz             (k,i);
      P3DiagnosticInputs_dpres         (i,k) = dpres          (k,i);
      P3DiagnosticInputs_inv_exner     (i,k) = inv_exner      (k,i);
      P3DiagnosticInputs_qv_prev       (i,k) = qv_prev        (k,i);
      P3DiagnosticInputs_t_prev        (i,k) = t_prev         (k,i);
      P3Infrastructure_col_location    (i,k) = col_location   (k,i);
    });

    Kokkos::fence();

    P3PrognosticState   prog_state;
    P3DiagnosticInputs  diag_inputs;
    P3DiagnosticOutputs diag_outputs;
    P3Infrastructure    infrastructure;
    P3HistoryOnly       history_only;

    prog_state.qc                   = P3PrognosticState_qc                  ;
    prog_state.nc                   = P3PrognosticState_nc                  ;
    prog_state.qr                   = P3PrognosticState_qr                  ;
    prog_state.nr                   = P3PrognosticState_nr                  ;
    prog_state.qi                   = P3PrognosticState_qi                  ;
    prog_state.qm                   = P3PrognosticState_qm                  ;
    prog_state.ni                   = P3PrognosticState_ni                  ;
    prog_state.bm                   = P3PrognosticState_bm                  ;
    prog_state.qv                   = P3PrognosticState_qv                  ;
    prog_state.th                   = P3PrognosticState_th                  ;
    diag_inputs.nc_nuceat_tend      = P3DiagnosticInputs_nc_nuceat_tend     ;
    diag_inputs.nccn                = P3DiagnosticInputs_nccn               ;
    diag_inputs.ni_activated        = P3DiagnosticInputs_ni_activated       ;
    diag_inputs.inv_qc_relvar       = P3DiagnosticInputs_inv_qc_relvar      ;
    diag_inputs.cld_frac_i          = P3DiagnosticInputs_cld_frac_i         ;
    diag_inputs.cld_frac_l          = P3DiagnosticInputs_cld_frac_l         ;
    diag_inputs.cld_frac_r          = P3DiagnosticInputs_cld_frac_r         ;
    diag_inputs.pres                = P3DiagnosticInputs_pres               ;
    diag_inputs.dz                  = P3DiagnosticInputs_dz                 ;
    diag_inputs.dpres               = P3DiagnosticInputs_dpres              ;
    diag_inputs.inv_exner           = P3DiagnosticInputs_inv_exner          ;
    diag_inputs.qv_prev             = P3DiagnosticInputs_qv_prev            ;
    diag_inputs.t_prev              = P3DiagnosticInputs_t_prev             ;
    diag_outputs.qv2qi_depos_tend   = P2DiagnosticOutputs_qv2qi_depos_tend  ;
    diag_outputs.precip_liq_surf    = P3DiagnosticOutputs_precip_liq_surf   ;
    diag_outputs.precip_ice_surf    = P3DiagnosticOutputs_precip_ice_surf   ;
    diag_outputs.diag_eff_radius_qc = P3DiagnosticOutputs_diag_eff_radius_qc;
    diag_outputs.diag_eff_radius_qi = P3DiagnosticOutputs_diag_eff_radius_qi;
    diag_outputs.rho_qi             = P3DiagnosticOutputs_rho_qi            ;
    diag_outputs.precip_liq_flux    = P3DiagnosticOutputs_precip_liq_flux   ;
    diag_outputs.precip_ice_flux    = P3DiagnosticOutputs_precip_ice_flux   ;
    infrastructure.dt               = P3Infrastructure_dt                   ;
    infrastructure.it               = P3Infrastructure_it                   ;
    infrastructure.its              = P3Infrastructure_its                  ;
    infrastructure.ite              = P3Infrastructure_ite                  ;
    infrastructure.kts              = P3Infrastructure_kts                  ;
    infrastructure.kte              = P3Infrastructure_kte                  ;
    infrastructure.predictNc        = P3Infrastructure_predictNc            ;
    infrastructure.prescribedCCN    = P3Infrastructure_prescribedCCN        ;
    infrastructure.col_location     = P3Infrastructure_col_location         ;
    history_only.liq_ice_exchange   = P3HistoryOnly_liq_ice_exchange        ;
    history_only.vap_liq_exchange   = P3HistoryOnly_vap_liq_exchange        ;
    history_only.vap_ice_exchange   = P3HistoryOnly_vap_ice_exchange        ;

    auto elapsed_microsec = P3::p3_main(prog_state, diag_inputs, diag_outputs, infrastructure,
                                        history_only, lookup_tables, workspace_mgr, ncol, nlev);



  }

}


