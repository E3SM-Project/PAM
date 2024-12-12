
#include "EKAT_utils.h"
#include "p3_functions.hpp"

namespace pam {

  extern scream::p3::Functions<scream::Real,scream::DefaultDevice>::P3LookupTables     lookup_tables;

  void deallocate_scream_cxx_globals() {
    lookup_tables = decltype(lookup_tables)();
  }

  void call_kokkos_finalize() {
    Kokkos::finalize();
  }

}


