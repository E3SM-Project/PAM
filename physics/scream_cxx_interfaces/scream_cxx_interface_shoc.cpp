
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
                     array_ir::ArrayIR<double,2> const & shoc_ql2 ) {  //   out
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

    SHOC::SHOCInput shoc_input{host_dx_1d, host_dy_1d, zt_grid_2d, zi_grid_2d,
                               pres_2d, presi_2d, pdel_2d, thv_2d,
                               w_field_2d, wthl_sfc_1d, wqw_sfc_1d, uw_sfc_1d,
                               vw_sfc_1d, wtracer_sfc_2d, inv_exner_2d, phis_1d};

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

    SHOC::SHOCInputOutput shoc_input_output{host_dse_2d, tke_2d, thetal_2d, qw_2d,
                                           shoc_hwind_3d, wthv_sec_2d, qtracers_3d,
                                           tk_2d, shoc_cldfrac_2d, shoc_ql_2d};

    auto pblh_1d     = ArrayIR_to_View         (pblh    );
    auto shoc_ql2_2d = ArrayIR_to_View_of_Packs(shoc_ql2);

    SHOC::SHOCOutput shoc_output{pblh_1d, shoc_ql2_2d};

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

    SHOC::SHOCHistoryOutput shoc_history_output{shoc_mix_2d, w_sec_2d, thl_sec_2d, qw_sec_2d,
                                                qwthl_sec_2d, wthl_sec_2d, wqw_sec_2d, wtke_sec_2d,
                                                uw_sec_2d, vw_sec_2d, w3_2d, wqls_sec_2d, brunt_2d, isotropy_2d};

    const int nwind = ekat::npack<Spack>(2)*Spack::n;
    const int ntrac = ekat::npack<Spack>(num_qtracers+3)*Spack::n;
    const auto policy = ekat::ExeSpaceUtils<SHOC::KT::ExeSpace>::get_default_team_policy(ncol, npack);
    ekat::WorkspaceManager<Spack, SHOC::KT::Device> workspace_mgr(nipack, 13+(nwind+ntrac), policy);

    const auto elapsed_microsec = SHOC::shoc_main(shcol, nlev, nlevi, nlev, nadv, num_qtracers, dtime, workspace_mgr,
                                                  shoc_input, shoc_input_output, shoc_output, shoc_history_output);
  }

}


