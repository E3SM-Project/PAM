
#include "EKAT_utils.h"
#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"
#include "shoc_f90.hpp"

namespace pam {

  void shoc_main_cxx(int                         const & shcol,        // in
                     int                         const & nlev,         // in
                     int                         const & nlevi,        // in
                     double                      const & dtime,        // in
                     int                         const & nadv,         // in
                     array_ir::ArrayIR<double,1> const & host_dx,      // in
                     array_ir::ArrayIR<double,1> const & host_dy,      // in
                     array_ir::ArrayIR<double,2> const & thv,          // in
                     array_ir::ArrayIR<double,2> const & zt_grid,      // in
                     array_ir::ArrayIR<double,2> const & zi_grid,      // in
                     array_ir::ArrayIR<double,2> const & pres,         // in
                     array_ir::ArrayIR<double,2> const & presi,        // in
                     array_ir::ArrayIR<double,2> const & pdel,         // in
                     array_ir::ArrayIR<double,1> const & wthl_sfc,     // in
                     array_ir::ArrayIR<double,1> const & wqw_sfc,      // in
                     array_ir::ArrayIR<double,1> const & uw_sfc,       // in
                     array_ir::ArrayIR<double,1> const & vw_sfc,       // in
                     array_ir::ArrayIR<double,2> const & wtracer_sfc,  // in
                     int                         const & num_qtracers, // in
                     array_ir::ArrayIR<double,2> const & w_field,      // in
                     array_ir::ArrayIR<double,2> const & inv_exner,    // in
                     array_ir::ArrayIR<double,1> const & phis,         // in
                     array_ir::ArrayIR<double,2> const & host_dse,     // inout
                     array_ir::ArrayIR<double,2> const & tke,          // inout
                     array_ir::ArrayIR<double,2> const & thetal,       // inout
                     array_ir::ArrayIR<double,2> const & qw,           // inout
                     array_ir::ArrayIR<double,3> const & shoc_hwind,   // inout - (col,2,lev) - u is index zero - v is index one
                     array_ir::ArrayIR<double,3> const & qtracers,     // inout - (col,qtr,lev)
                     array_ir::ArrayIR<double,2> const & wthv_sec,     // inout
                     array_ir::ArrayIR<double,2> const & tk,           // inout
                     array_ir::ArrayIR<double,2> const & shoc_ql,      // inout
                     array_ir::ArrayIR<double,2> const & shoc_cldfrac, // inout
                     array_ir::ArrayIR<double,1> const & pblh,         //   out
                     array_ir::ArrayIR<double,2> const & shoc_mix,     //   out
                     array_ir::ArrayIR<double,2> const & isotropy,     //   out
                     array_ir::ArrayIR<double,2> const & w_sec,        //   out
                     array_ir::ArrayIR<double,2> const & thl_sec,      //   out
                     array_ir::ArrayIR<double,2> const & qw_sec,       //   out
                     array_ir::ArrayIR<double,2> const & qwthl_sec,    //   out
                     array_ir::ArrayIR<double,2> const & wthl_sec,     //   out
                     array_ir::ArrayIR<double,2> const & wqw_sec,      //   out
                     array_ir::ArrayIR<double,2> const & wtke_sec,     //   out
                     array_ir::ArrayIR<double,2> const & uw_sec,       //   out
                     array_ir::ArrayIR<double,2> const & vw_sec,       //   out
                     array_ir::ArrayIR<double,2> const & w3,           //   out
                     array_ir::ArrayIR<double,2> const & wqls_sec,     //   out
                     array_ir::ArrayIR<double,2> const & brunt,        //   out
                     array_ir::ArrayIR<double,2> const & shoc_ql2,     //   out
                     array_ir::ArrayIR<double,2> const & shoc_tkh ) {  //   out
    using ScreamCXX::ArrayIR_to_View_of_Packs;
    using ScreamCXX::ArrayIR_to_View;
    using namespace scream;
    using namespace scream::shoc;
    using SHOC       = shoc::Functions<Real, DefaultDevice>;
    using Spack      = typename SHOC::Spack;

    if (! Kokkos::is_initialized()) { Kokkos::initialize(); }

    const int ncol   = shcol;
    const int npack  = ekat::npack<Spack>(nlev);
    const int nipack = ekat::npack<Spack>(nlevi);

    //--------------------------------------------------------------------------
    // Input Variables

    auto host_dx_1d     = ArrayIR_to_View         (host_dx    );
    auto host_dy_1d     = ArrayIR_to_View         (host_dy    );
    auto wthl_sfc_1d    = ArrayIR_to_View         (wthl_sfc   );
    auto wqw_sfc_1d     = ArrayIR_to_View         (wqw_sfc    );
    auto uw_sfc_1d      = ArrayIR_to_View         (uw_sfc     );
    auto vw_sfc_1d      = ArrayIR_to_View         (vw_sfc     );
    auto phis_1d        = ArrayIR_to_View         (phis       );
    auto zt_grid_2d     = ArrayIR_to_View_of_Packs(zt_grid    );
    auto zi_grid_2d     = ArrayIR_to_View_of_Packs(zi_grid    );
    auto pres_2d        = ArrayIR_to_View_of_Packs(pres       );
    auto presi_2d       = ArrayIR_to_View_of_Packs(presi      );
    auto pdel_2d        = ArrayIR_to_View_of_Packs(pdel       );
    auto thv_2d         = ArrayIR_to_View_of_Packs(thv        );
    auto w_field_2d     = ArrayIR_to_View_of_Packs(w_field    );
    auto wtracer_sfc_2d = ArrayIR_to_View_of_Packs(wtracer_sfc);
    auto inv_exner_2d   = ArrayIR_to_View_of_Packs(inv_exner  );

    SHOC::SHOCInput shoc_input;
    shoc_input.dx          = host_dx_1d;
    shoc_input.dy          = host_dy_1d;
    shoc_input.zt_grid     = zt_grid_2d;
    shoc_input.zi_grid     = zi_grid_2d;
    shoc_input.pres        = pres_2d;
    shoc_input.presi       = presi_2d;
    shoc_input.pdel        = pdel_2d;
    shoc_input.thv         = thv_2d;
    shoc_input.w_field     = w_field_2d;     // wm_zt;
    shoc_input.wthl_sfc    = wthl_sfc_1d;    // wpthlp_sfc;
    shoc_input.wqw_sfc     = wqw_sfc_1d;     // wprtp_sfc;
    shoc_input.uw_sfc      = uw_sfc_1d;      // upwp_sfc;
    shoc_input.vw_sfc      = vw_sfc_1d;      // vpwp_sfc;
    shoc_input.wtracer_sfc = wtracer_sfc_2d;
    shoc_input.inv_exner   = inv_exner_2d;
    shoc_input.phis        = phis_1d;

    //--------------------------------------------------------------------------
    // Input/Output Variables

    auto host_dse_2d     = ArrayIR_to_View_of_Packs(host_dse    );
    auto tke_2d          = ArrayIR_to_View_of_Packs(tke         );
    auto thetal_2d       = ArrayIR_to_View_of_Packs(thetal      );
    auto qw_2d           = ArrayIR_to_View_of_Packs(qw          );
    auto shoc_ql_2d      = ArrayIR_to_View_of_Packs(shoc_ql     );
    auto wthv_sec_2d     = ArrayIR_to_View_of_Packs(wthv_sec    );
    auto tk_2d           = ArrayIR_to_View_of_Packs(tk          );
    auto shoc_cldfrac_2d = ArrayIR_to_View_of_Packs(shoc_cldfrac);
    auto shoc_hwind_3d   = ArrayIR_to_View_of_Packs(shoc_hwind  );
    auto qtracers_3d     = ArrayIR_to_View_of_Packs(qtracers    );

    SHOC::SHOCInputOutput shoc_input_output;
    shoc_input_output.host_dse     = host_dse_2d;
    shoc_input_output.tke          = tke_2d;
    shoc_input_output.thetal       = thetal_2d;
    shoc_input_output.qw           = qw_2d;
    shoc_input_output.horiz_wind   = shoc_hwind_3d;
    shoc_input_output.wthv_sec     = wthv_sec_2d;
    shoc_input_output.qtracers     = qtracers_3d;
    shoc_input_output.tk           = tk_2d;
    shoc_input_output.shoc_cldfrac = shoc_cldfrac_2d;
    shoc_input_output.shoc_ql      = shoc_ql_2d;

    //--------------------------------------------------------------------------
    // Output Variables

    auto pblh_1d     = ArrayIR_to_View         (pblh    );
    auto shoc_ql2_2d = ArrayIR_to_View_of_Packs(shoc_ql2);
    auto shoc_tkh_2d = ArrayIR_to_View_of_Packs(shoc_tkh);

    SHOC::SHOCOutput shoc_output;
    shoc_output.pblh     = pblh_1d;
    shoc_output.shoc_ql2 = shoc_ql2_2d;
    shoc_output.tkh      = shoc_tkh_2d;


    //--------------------------------------------------------------------------
    // Diagnostic Output

    auto shoc_mix_2d  = ArrayIR_to_View_of_Packs(shoc_mix );
    auto w_sec_2d     = ArrayIR_to_View_of_Packs(w_sec    );
    auto thl_sec_2d   = ArrayIR_to_View_of_Packs(thl_sec  );
    auto qw_sec_2d    = ArrayIR_to_View_of_Packs(qw_sec   );
    auto qwthl_sec_2d = ArrayIR_to_View_of_Packs(qwthl_sec);
    auto wthl_sec_2d  = ArrayIR_to_View_of_Packs(wthl_sec );
    auto wqw_sec_2d   = ArrayIR_to_View_of_Packs(wqw_sec  );
    auto wtke_sec_2d  = ArrayIR_to_View_of_Packs(wtke_sec );
    auto uw_sec_2d    = ArrayIR_to_View_of_Packs(uw_sec   );
    auto vw_sec_2d    = ArrayIR_to_View_of_Packs(vw_sec   );
    auto w3_2d        = ArrayIR_to_View_of_Packs(w3       );
    auto wqls_sec_2d  = ArrayIR_to_View_of_Packs(wqls_sec );
    auto brunt_2d     = ArrayIR_to_View_of_Packs(brunt    );
    auto isotropy_2d  = ArrayIR_to_View_of_Packs(isotropy );

    SHOC::SHOCHistoryOutput shoc_history_output;
    shoc_history_output.shoc_mix  = shoc_mix_2d;
    shoc_history_output.isotropy  = isotropy_2d;
    shoc_history_output.w_sec     = w_sec_2d;
    shoc_history_output.thl_sec   = thl_sec_2d;
    shoc_history_output.qw_sec    = qw_sec_2d;
    shoc_history_output.qwthl_sec = qwthl_sec_2d;
    shoc_history_output.wthl_sec  = wthl_sec_2d;
    shoc_history_output.wqw_sec   = wqw_sec_2d;
    shoc_history_output.wtke_sec  = wtke_sec_2d;
    shoc_history_output.uw_sec    = uw_sec_2d;
    shoc_history_output.vw_sec    = vw_sec_2d;
    shoc_history_output.w3        = w3_2d;
    shoc_history_output.wqls_sec  = wqls_sec_2d;
    shoc_history_output.brunt     = brunt_2d;

    //--------------------------------------------------------------------------

    // hardcode runtime options to match common scream settings for now
    SHOC::SHOCRuntime shoc_runtime_options {
        0.001, // lambda_low
        0.04, // lambda_high
        2.65, // lambda_slope
        0.02, // lambda_thresh
        1.0, // thl2tune
        1.0, // qw2tune
        1.0, // qwthl2tune
        1.0, // w2tune
        0.5, // length_fac
        7.0, // c_diag_3rd_mom
        0.1, // Ckh
        0.1 // Ckm
    };

    const int nwind = ekat::npack<Spack>(2)*Spack::n;
    const int ntrac = ekat::npack<Spack>(num_qtracers+3)*Spack::n;
    const auto policy = ekat::ExeSpaceUtils<SHOC::KT::ExeSpace>::get_default_team_policy(ncol, npack);
    ekat::WorkspaceManager<Spack, SHOC::KT::Device> workspace_mgr(nipack, 14+(nwind+ntrac), policy);

    const auto elapsed_microsec = SHOC::shoc_main(shcol, nlev, nlevi, nlev, nadv, num_qtracers, dtime, workspace_mgr,
                                                  shoc_runtime_options, shoc_input, shoc_input_output, shoc_output, shoc_history_output);
  }

}


