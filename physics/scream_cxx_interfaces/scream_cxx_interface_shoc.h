
#pragma once 

#include "ArrayIR.h"

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
                     array_ir::ArrayIR<double,2> const & shoc_tkh );   //   out

}


