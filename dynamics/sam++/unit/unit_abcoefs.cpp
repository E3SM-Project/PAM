
#include "abcoefs.h"

int main() {
  real1d dt3 = real1d("dt3",3);
  real   at, bt, ct;
  int    na, nb, nc;

  dt3(0) = 4.2;
  dt3(1) = 8.4;
  dt3(2) = 5.1;

  // Test 1
  na = 1 ; nb = 2 ; nc = 3;
  abcoefs(dt3, na, nb, nc, at, bt, ct);
  std::cout << "abcoefs Fortran Test 1\n";
  std::cout << "at: " << at << "\n";
  std::cout << "bt: " << bt << "\n";
  std::cout << "ct: " << ct << "\n\n";

  // Test 2
  na = 3 ; nb = 1 ; nc = 2;
  abcoefs(dt3, na, nb, nc, at, bt, ct);
  std::cout << "abcoefs Fortran Test 2\n";
  std::cout << "at: " << at << "\n";
  std::cout << "bt: " << bt << "\n";
  std::cout << "ct: " << ct << "\n\n";

  // Test 3
  na = 2 ; nb = 3 ; nc = 1;
  abcoefs(dt3, na, nb, nc, at, bt, ct);
  std::cout << "abcoefs Fortran Test 3\n";
  std::cout << "at: " << at << "\n";
  std::cout << "bt: " << bt << "\n";
  std::cout << "ct: " << ct << "\n\n";
}

