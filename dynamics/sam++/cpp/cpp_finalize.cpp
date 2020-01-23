
#include "YAKL.h"

extern "C" void cpp_finalize() {
  yakl::finalize();
}


