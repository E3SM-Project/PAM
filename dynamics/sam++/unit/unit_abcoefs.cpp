
#include "abcoefs.h"

extern "C" void abcoefs_f90(real *dt3, int nstep, int na, int nb, int nc, real &at, real &bt, real &ct);

int main() {
  real1d dt3 = real1d("dt3",3);
  real   at, bt, ct;
  real   at_f90, bt_f90, ct_f90;
  int    na, nb, nc, nstep;

  nstep = 1;

  std::cout << "Got here 0" << std::endl;

  dt3(0) = 4.1;
  dt3(1) = 4.1;
  dt3(2) = 4.1;

  // Test 1
  na = 1 ; nb = 2 ; nc = 3;
  std::cout << "Got here 1" << std::endl;
  abcoefs    (dt3       , nstep, na, nb, nc, at    , bt    , ct    );
  std::cout << "Got here 2" << std::endl;
  abcoefs_f90(dt3.data(), nstep, na, nb, nc, at_f90, bt_f90, ct_f90);
  std::cout << "abcoefs Test 1\n";
  std::cout << "at: " << fabs(at-at_f90) << " " << at << " " << at_f90 << "\n";
  std::cout << "bt: " << fabs(bt-bt_f90) << " " << bt << " " << bt_f90 << "\n";
  std::cout << "ct: " << fabs(ct-ct_f90) << " " << ct << " " << ct_f90 << "\n\n";

  // Test 2
  na = 3 ; nb = 1 ; nc = 2;
  abcoefs    (dt3       , nstep, na, nb, nc, at    , bt    , ct    );
  abcoefs_f90(dt3.data(), nstep, na, nb, nc, at_f90, bt_f90, ct_f90);
  std::cout << "abcoefs Test 2\n";
  std::cout << "at: " << fabs(at-at_f90) << " " << at << " " << at_f90 << "\n";
  std::cout << "bt: " << fabs(bt-bt_f90) << " " << bt << " " << bt_f90 << "\n";
  std::cout << "ct: " << fabs(ct-ct_f90) << " " << ct << " " << ct_f90 << "\n\n";

  // Test 3
  na = 2 ; nb = 3 ; nc = 1;
  abcoefs    (dt3       , nstep, na, nb, nc, at    , bt    , ct    );
  abcoefs_f90(dt3.data(), nstep, na, nb, nc, at_f90, bt_f90, ct_f90);
  std::cout << "abcoefs Test 3\n";
  std::cout << "at: " << fabs(at-at_f90) << " " << at << " " << at_f90 << "\n";
  std::cout << "bt: " << fabs(bt-bt_f90) << " " << bt << " " << bt_f90 << "\n";
  std::cout << "ct: " << fabs(ct-ct_f90) << " " << ct << " " << ct_f90 << "\n\n";

  nstep = 2;

  // Test 4
  na = 1 ; nb = 2 ; nc = 3;
  abcoefs    (dt3       , nstep, na, nb, nc, at    , bt    , ct    );
  abcoefs_f90(dt3.data(), nstep, na, nb, nc, at_f90, bt_f90, ct_f90);
  std::cout << "abcoefs Test 4\n";
  std::cout << "at: " << fabs(at-at_f90) << " " << at << " " << at_f90 << "\n";
  std::cout << "bt: " << fabs(bt-bt_f90) << " " << bt << " " << bt_f90 << "\n";
  std::cout << "ct: " << fabs(ct-ct_f90) << " " << ct << " " << ct_f90 << "\n\n";

  // Test 5
  na = 3 ; nb = 1 ; nc = 2;
  abcoefs    (dt3       , nstep, na, nb, nc, at    , bt    , ct    );
  abcoefs_f90(dt3.data(), nstep, na, nb, nc, at_f90, bt_f90, ct_f90);
  std::cout << "abcoefs Test 5\n";
  std::cout << "at: " << fabs(at-at_f90) << " " << at << " " << at_f90 << "\n";
  std::cout << "bt: " << fabs(bt-bt_f90) << " " << bt << " " << bt_f90 << "\n";
  std::cout << "ct: " << fabs(ct-ct_f90) << " " << ct << " " << ct_f90 << "\n\n";

  // Test 6
  na = 2 ; nb = 3 ; nc = 1;
  abcoefs    (dt3       , nstep, na, nb, nc, at    , bt    , ct    );
  abcoefs_f90(dt3.data(), nstep, na, nb, nc, at_f90, bt_f90, ct_f90);
  std::cout << "abcoefs Test 6\n";
  std::cout << "at: " << fabs(at-at_f90) << " " << at << " " << at_f90 << "\n";
  std::cout << "bt: " << fabs(bt-bt_f90) << " " << bt << " " << bt_f90 << "\n";
  std::cout << "ct: " << fabs(ct-ct_f90) << " " << ct << " " << ct_f90 << "\n\n";

  nstep = 3;

  // Test 7
  na = 1 ; nb = 2 ; nc = 3;
  abcoefs    (dt3       , nstep, na, nb, nc, at    , bt    , ct    );
  abcoefs_f90(dt3.data(), nstep, na, nb, nc, at_f90, bt_f90, ct_f90);
  std::cout << "abcoefs Test 7\n";
  std::cout << "at: " << fabs(at-at_f90) << " " << at << " " << at_f90 << "\n";
  std::cout << "bt: " << fabs(bt-bt_f90) << " " << bt << " " << bt_f90 << "\n";
  std::cout << "ct: " << fabs(ct-ct_f90) << " " << ct << " " << ct_f90 << "\n\n";

  // Test 8
  na = 3 ; nb = 1 ; nc = 2;
  abcoefs    (dt3       , nstep, na, nb, nc, at    , bt    , ct    );
  abcoefs_f90(dt3.data(), nstep, na, nb, nc, at_f90, bt_f90, ct_f90);
  std::cout << "abcoefs Test 8\n";
  std::cout << "at: " << fabs(at-at_f90) << " " << at << " " << at_f90 << "\n";
  std::cout << "bt: " << fabs(bt-bt_f90) << " " << bt << " " << bt_f90 << "\n";
  std::cout << "ct: " << fabs(ct-ct_f90) << " " << ct << " " << ct_f90 << "\n\n";

  // Test 9
  na = 2 ; nb = 3 ; nc = 1;
  abcoefs    (dt3       , nstep, na, nb, nc, at    , bt    , ct    );
  abcoefs_f90(dt3.data(), nstep, na, nb, nc, at_f90, bt_f90, ct_f90);
  std::cout << "abcoefs Test 9\n";
  std::cout << "at: " << fabs(at-at_f90) << " " << at << " " << at_f90 << "\n";
  std::cout << "bt: " << fabs(bt-bt_f90) << " " << bt << " " << bt_f90 << "\n";
  std::cout << "ct: " << fabs(ct-ct_f90) << " " << ct << " " << ct_f90 << "\n\n";
}

