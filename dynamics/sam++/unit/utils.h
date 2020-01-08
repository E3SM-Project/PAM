
#include <string>
#include <iostream>
#include <cmath>

// Computes a relative difference between v1 and v2, where v1 is the normalizer.
// Compares the relative difference against the tolerance
//    If rel diff <= tol, returns 0
//    Otherwise returns -1
// If v1 is zero, then the absolute difference is compared against tolerance to avoid
// division by zero
template <class T> int compareScalar( std::string const &label , T v1 , T v2 , T tol ) {
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
