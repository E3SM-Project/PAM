
#pragma once

#include "ArrayIR.h"

namespace pam {

  void p3_init_lookup_tables();

  void p3_main_cxx(array_ir::ArrayIR<double,2> const & qc,                 // inout
                   array_ir::ArrayIR<double,2> const & nc,                 // inout
                   array_ir::ArrayIR<double,2> const & qr,                 // inout
                   array_ir::ArrayIR<double,2> const & nr,                 // inout
                   array_ir::ArrayIR<double,2> const & th_atm,             // inout
                   array_ir::ArrayIR<double,2> const & qv,                 // inout
                   double                      const & dt,                 // inout
                   array_ir::ArrayIR<double,2> const & qi,                 // inout
                   array_ir::ArrayIR<double,2> const & qm,                 // inout
                   array_ir::ArrayIR<double,2> const & ni,                 // inout
                   array_ir::ArrayIR<double,2> const & bm,                 // inout
                   array_ir::ArrayIR<double,2> const & pres,               // in
                   array_ir::ArrayIR<double,2> const & dz,                 // in
                   array_ir::ArrayIR<double,2> const & nc_nuceat_tend,     // in
                   array_ir::ArrayIR<double,2> const & nccn_prescribed,    // in
                   array_ir::ArrayIR<double,2> const & ni_activated,       // in
                   array_ir::ArrayIR<double,2> const & inv_qc_relvar,      // in
                   int                         const & it,                 // in
                   array_ir::ArrayIR<double,1> const & precip_liq_surf,    //   out
                   array_ir::ArrayIR<double,1> const & precip_ice_surf,    //   out
                   int                         const & its,                // in
                   int                         const & ite,                // in
                   int                         const & kts,                // in
                   int                         const & kte,                // in
                   array_ir::ArrayIR<double,2> const & diag_eff_radius_qc, //   out
                   array_ir::ArrayIR<double,2> const & diag_eff_radius_qi, //   out
                   array_ir::ArrayIR<double,2> const & rho_qi,             //   out
                   bool                        const & do_predict_nc,      // in
                   bool                        const & do_prescribed_CCN,  // in
                   array_ir::ArrayIR<double,2> const & dpres,              // in
                   array_ir::ArrayIR<double,2> const & inv_exner,          // in
                   array_ir::ArrayIR<double,2> const & qv2qi_depos_tend,   //   out
                   array_ir::ArrayIR<double,2> const & precip_liq_flux,    //   out
                   array_ir::ArrayIR<double,2> const & precip_ice_flux,    //   out
                   array_ir::ArrayIR<double,2> const & cld_frac_r,         // in
                   array_ir::ArrayIR<double,2> const & cld_frac_l,         // in
                   array_ir::ArrayIR<double,2> const & cld_frac_i,         // in
                   array_ir::ArrayIR<double,2> const & liq_ice_exchange,   //   out
                   array_ir::ArrayIR<double,2> const & vap_liq_exchange,   //   out
                   array_ir::ArrayIR<double,2> const & vap_ice_exchange,   //   out
                   array_ir::ArrayIR<double,2> const & qv_prev,            // in
                   array_ir::ArrayIR<double,2> const & t_prev,             // in
                   array_ir::ArrayIR<double,2> const & col_location,       // in
                   double                              *elapsed_s );       //   out

}


