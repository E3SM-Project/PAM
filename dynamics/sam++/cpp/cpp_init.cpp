
#include "YAKL.h"

extern "C" void cpp_init() {
  yakl::init( NCRMS*CRM_NX*CRM_NY*CRM_NZ*1500 );
}


