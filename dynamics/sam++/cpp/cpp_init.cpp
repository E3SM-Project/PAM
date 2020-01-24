
#include "YAKL.h"
#include "vars.h"

extern "C" void cpp_init() {
  yakl::init( ncrms*crm_nx*crm_ny*crm_nz*1500 );
  // yakl::init();
}


