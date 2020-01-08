
#include <string>
#include <iostream>
#include <cmath>
#include "abcoefs.h"



// Computes a relative difference between v1 and v2, where v1 is the normalizer.
// Compares the relative difference against the tolerance
//    If rel diff <= tol, returns 0
//    Otherwise returns -1
// If v1 is zero, then the absolute difference is compared against tolerance to avoid
// division by zero
inline template <class T> int compareScalar( std::string const &label , T v1 , T v2 , T tol ) {
  double diff = std::abs( (double) v1 - (double) v2);
  if (v1 != 0) { diff /= (double) v1; }
  bool pass = diff <= tol;
  std::cout << std::scientific;
  std::cout << label << ": " << (pass ? "PASS" : "FAIL") << std::endl;
  if (!pass) {
    std::cout << "Rel Diff:  " << diff << std::endl;
    std::cout << "v1:        " << v1   << std::endl;
    std::cout << "v2:        " << v2   << std::endl;
    std::cout << "Tolerance: " << tol  << std::endl;
  }
  std::cout << std::endl;
  return pass ? 0 : -1;
}



inline void initDomainStep3NoCycle( Domain &dom ) {
  dom.dx      = 500.;
  dom.dy      = 500.;
  dom.ncycle  = 1;
  dom.icycle  = 1;
  dom.dt      = 1.;
  dom.dtn     = 5.;
  dom.dtfactor = 1.;
  dom.na      = 1;
  dom.nb      = 2;
  dom.nc      = 3;
  dom.dt3(0)  = 1.;
  dom.dt3(1)  = 1.;
  dom.dt3(2)  = 1.;
  dom.nstep   = 3;
  dom.nstop   = 100;
  abcoefs( dom );  // Initializes at, bt, and ct

  for (int icrm = 0; icrm < ncrms; icrm++) {
    
  }
}




