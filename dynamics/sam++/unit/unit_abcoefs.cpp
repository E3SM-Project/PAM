
#include "abcoefs.h"
#include "utils.h"

extern "C" void abcoefs_f90(real *dt3, int nstep, int na, int nb, int nc, real &at, real &bt, real &ct);


int main() {
  real1d dt3 = real1d("dt3",3);
  real   at, bt, ct;
  real   at_f90, bt_f90, ct_f90;
  int    na, nb, nc, nstep;
  real   tol = 1.e-14;
  int    ret = 0;

  nstep = 1;

  dt3(0) = 4.1;
  dt3(1) = 4.2;
  dt3(2) = 4.3;

  std::cout << std::scientific;

  // Test 1
  na = 1 ; nb = 2 ; nc = 3;
  abcoefs    (dt3       , nstep, na, nb, nc, at    , bt    , ct    );
  abcoefs_f90(dt3.data(), nstep, na, nb, nc, at_f90, bt_f90, ct_f90);
  ret += compareScalar( "Test 1, at" , at_f90 , at , tol );
  ret += compareScalar( "Test 1, bt" , bt_f90 , bt , tol );
  ret += compareScalar( "Test 1, ct" , ct_f90 , ct , tol );

  // Test 2
  na = 3 ; nb = 1 ; nc = 2;
  abcoefs    (dt3       , nstep, na, nb, nc, at    , bt    , ct    );
  abcoefs_f90(dt3.data(), nstep, na, nb, nc, at_f90, bt_f90, ct_f90);
  ret += compareScalar( "Test 2, at" , at_f90 , at , tol );
  ret += compareScalar( "Test 2, bt" , bt_f90 , bt , tol );
  ret += compareScalar( "Test 2, ct" , ct_f90 , ct , tol );

  // Test 3
  na = 2 ; nb = 3 ; nc = 1;
  abcoefs    (dt3       , nstep, na, nb, nc, at    , bt    , ct    );
  abcoefs_f90(dt3.data(), nstep, na, nb, nc, at_f90, bt_f90, ct_f90);
  ret += compareScalar( "Test 3, at" , at_f90 , at , tol );
  ret += compareScalar( "Test 3, bt" , bt_f90 , bt , tol );
  ret += compareScalar( "Test 3, ct" , ct_f90 , ct , tol );

  nstep = 2;

  // Test 4
  na = 1 ; nb = 2 ; nc = 3;
  abcoefs    (dt3       , nstep, na, nb, nc, at    , bt    , ct    );
  abcoefs_f90(dt3.data(), nstep, na, nb, nc, at_f90, bt_f90, ct_f90);
  ret += compareScalar( "Test 4, at" , at_f90 , at , tol );
  ret += compareScalar( "Test 4, bt" , bt_f90 , bt , tol );
  ret += compareScalar( "Test 4, ct" , ct_f90 , ct , tol );

  // Test 5
  na = 3 ; nb = 1 ; nc = 2;
  abcoefs    (dt3       , nstep, na, nb, nc, at    , bt    , ct    );
  abcoefs_f90(dt3.data(), nstep, na, nb, nc, at_f90, bt_f90, ct_f90);
  ret += compareScalar( "Test 5, at" , at_f90 , at , tol );
  ret += compareScalar( "Test 5, bt" , bt_f90 , bt , tol );
  ret += compareScalar( "Test 5, ct" , ct_f90 , ct , tol );

  // Test 6
  na = 2 ; nb = 3 ; nc = 1;
  abcoefs    (dt3       , nstep, na, nb, nc, at    , bt    , ct    );
  abcoefs_f90(dt3.data(), nstep, na, nb, nc, at_f90, bt_f90, ct_f90);
  ret += compareScalar( "Test 6, at" , at_f90 , at , tol );
  ret += compareScalar( "Test 6, bt" , bt_f90 , bt , tol );
  ret += compareScalar( "Test 6, ct" , ct_f90 , ct , tol );

  nstep = 3;

  // Test 7
  na = 1 ; nb = 2 ; nc = 3;
  abcoefs    (dt3       , nstep, na, nb, nc, at    , bt    , ct    );
  abcoefs_f90(dt3.data(), nstep, na, nb, nc, at_f90, bt_f90, ct_f90);
  ret += compareScalar( "Test 7, at" , at_f90 , at , tol );
  ret += compareScalar( "Test 7, bt" , bt_f90 , bt , tol );
  ret += compareScalar( "Test 7, ct" , ct_f90 , ct , tol );

  // Test 8
  na = 3 ; nb = 1 ; nc = 2;
  abcoefs    (dt3       , nstep, na, nb, nc, at    , bt    , ct    );
  abcoefs_f90(dt3.data(), nstep, na, nb, nc, at_f90, bt_f90, ct_f90);
  ret += compareScalar( "Test 8, at" , at_f90 , at , tol );
  ret += compareScalar( "Test 8, bt" , bt_f90 , bt , tol );
  ret += compareScalar( "Test 8, ct" , ct_f90 , ct , tol );

  // Test 9
  na = 2 ; nb = 3 ; nc = 1;
  abcoefs    (dt3       , nstep, na, nb, nc, at    , bt    , ct    );
  abcoefs_f90(dt3.data(), nstep, na, nb, nc, at_f90, bt_f90, ct_f90);
  ret += compareScalar( "Test 9, at" , at_f90 , at , tol );
  ret += compareScalar( "Test 9, bt" , bt_f90 , bt , tol );
  ret += compareScalar( "Test 9, ct" , ct_f90 , ct , tol );

  return ret;
}

