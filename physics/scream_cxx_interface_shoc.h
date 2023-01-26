
#pragma once 

#include "ArrayIR.h"

void shoc_main_cxx(int                         & shcol,        // in
                   int                         & nlev,         // in
                   int                         & nlevi,        // in
                   double                      & dtime,        // in
                   int                         & nadv,         // in
                   array_ir::ArrayIR<double,1> & host_dx,      // in
                   array_ir::ArrayIR<double,1> & host_dy,      // in
                   array_ir::ArrayIR<double,2> & thv,          // in
                   array_ir::ArrayIR<double,2> & zt_grid,      // in
                   array_ir::ArrayIR<double,2> & zi_grid,      // in
                   array_ir::ArrayIR<double,2> & pres,         // in
                   array_ir::ArrayIR<double,2> & presi,        // in
                   array_ir::ArrayIR<double,2> & pdel,         // in
                   array_ir::ArrayIR<double,1> & wthl_sfc,     // in
                   array_ir::ArrayIR<double,1> & wqw_sfc,      // in
                   array_ir::ArrayIR<double,1> & uw_sfc,       // in
                   array_ir::ArrayIR<double,1> & vw_sfc,       // in
                   array_ir::ArrayIR<double,2> & wtracer_sfc,  // in
                   int                         & num_qtracers, // in
                   array_ir::ArrayIR<double,2> & w_field,      // in
                   array_ir::ArrayIR<double,2> & inv_exner,    // in
                   array_ir::ArrayIR<double,1> & phis,         // in
                   array_ir::ArrayIR<double,2> & host_dse,     // inout
                   array_ir::ArrayIR<double,2> & tke,          // inout
                   array_ir::ArrayIR<double,2> & thetal,       // inout
                   array_ir::ArrayIR<double,2> & qw,           // inout
                   array_ir::ArrayIR<double,3> & shoc_hwind,   // inout - (col,2,lev) - u is index zero - v is index one
                   array_ir::ArrayIR<double,3> & qtracers,     // inout - (col,qtr,lev)
                   array_ir::ArrayIR<double,2> & wthv_sec,     // inout
                   array_ir::ArrayIR<double,2> & tk,           // inout
                   array_ir::ArrayIR<double,2> & shoc_ql,      // inout
                   array_ir::ArrayIR<double,2> & shoc_cldfrac, // inout
                   array_ir::ArrayIR<double,1> & pblh,         //   out
                   array_ir::ArrayIR<double,2> & shoc_mix,     //   out
                   array_ir::ArrayIR<double,2> & isotropy,     //   out
                   array_ir::ArrayIR<double,2> & w_sec,        //   out
                   array_ir::ArrayIR<double,2> & thl_sec,      //   out
                   array_ir::ArrayIR<double,2> & qw_sec,       //   out
                   array_ir::ArrayIR<double,2> & qwthl_sec,    //   out
                   array_ir::ArrayIR<double,2> & wthl_sec,     //   out
                   array_ir::ArrayIR<double,2> & wqw_sec,      //   out
                   array_ir::ArrayIR<double,2> & wtke_sec,     //   out
                   array_ir::ArrayIR<double,2> & uw_sec,       //   out
                   array_ir::ArrayIR<double,2> & vw_sec,       //   out
                   array_ir::ArrayIR<double,2> & w3,           //   out
                   array_ir::ArrayIR<double,2> & wqls_sec,     //   out
                   array_ir::ArrayIR<double,2> & brunt,        //   out
                   array_ir::ArrayIR<double,2> & shoc_ql2 );   //   out

