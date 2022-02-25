
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

namespace pam {

  void kokkos_initialize() {
    Kokkos::initialize();
  }


  void kokkos_finalize() {
    Kokkos::finalize();
  }

}


