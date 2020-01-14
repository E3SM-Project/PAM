
#include "Kokkos_Core.hpp"

extern "C" void cpp_finalize() {
  Kokkos::finalize();
}


