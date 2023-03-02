
#include "EKAT_utils.h"
#include "p3_functions.hpp"
#include "p3_functions_f90.hpp"
#include "p3_f90.hpp"

namespace pam {

  extern scream::p3::Functions<scream::Real,scream::DefaultDevice>::view_1d_table      mu_r_table_vals;
  extern scream::p3::Functions<scream::Real,scream::DefaultDevice>::view_2d_table      vn_table_vals;
  extern scream::p3::Functions<scream::Real,scream::DefaultDevice>::view_2d_table      vm_table_vals;
  extern scream::p3::Functions<scream::Real,scream::DefaultDevice>::view_2d_table      revap_table_vals;
  extern scream::p3::Functions<scream::Real,scream::DefaultDevice>::view_ice_table     ice_table_vals;
  extern scream::p3::Functions<scream::Real,scream::DefaultDevice>::view_collect_table collect_table_vals;
  extern scream::p3::Functions<scream::Real,scream::DefaultDevice>::view_dnu_table     dnu_table_vals;

  void deallocate_scream_cxx_globals() {
    mu_r_table_vals    = scream::p3::Functions<scream::Real,scream::DefaultDevice>::view_1d_table     ();
    vn_table_vals      = scream::p3::Functions<scream::Real,scream::DefaultDevice>::view_2d_table     ();
    vm_table_vals      = scream::p3::Functions<scream::Real,scream::DefaultDevice>::view_2d_table     ();
    revap_table_vals   = scream::p3::Functions<scream::Real,scream::DefaultDevice>::view_2d_table     ();
    ice_table_vals     = scream::p3::Functions<scream::Real,scream::DefaultDevice>::view_ice_table    ();
    collect_table_vals = scream::p3::Functions<scream::Real,scream::DefaultDevice>::view_collect_table();
    dnu_table_vals     = scream::p3::Functions<scream::Real,scream::DefaultDevice>::view_dnu_table    ();
  }

  void call_kokkos_finalize() {
    Kokkos::finalize();
  }

}


