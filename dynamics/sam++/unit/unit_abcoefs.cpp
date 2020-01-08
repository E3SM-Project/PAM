
#include "abcoefs.h"
#include "utils.h"

extern "C" void abcoefs_f90(real *dt3, int nstep, int na, int nb, int nc, real &at, real &bt, real &ct);


int main() {
  int  constexpr ncrms = 1;
  real constexpr tol = 1.e-14;
  Domain dom(ncrms);
  real   at_f90, bt_f90, ct_f90;
  int    ret = 0;

  dom.nstep = 1;

  dom.dt3(0) = 4.1;
  dom.dt3(1) = 4.2;
  dom.dt3(2) = 4.3;

  std::cout << std::scientific;

  // Test 1
  dom.na = 1 ; dom.nb = 2 ; dom.nc = 3;
  abcoefs    (dom);
  abcoefs_f90(dom.dt3.data(), dom.nstep, dom.na, dom.nb, dom.nc, at_f90, bt_f90, ct_f90);
  ret += compareScalar( "Test 1, at" , at_f90 , dom.at , tol );
  ret += compareScalar( "Test 1, bt" , bt_f90 , dom.bt , tol );
  ret += compareScalar( "Test 1, ct" , ct_f90 , dom.ct , tol );

  // Test 2
  dom.na = 3 ; dom.nb = 1 ; dom.nc = 2;
  abcoefs    (dom);
  abcoefs_f90(dom.dt3.data(), dom.nstep, dom.na, dom.nb, dom.nc, at_f90, bt_f90, ct_f90);
  ret += compareScalar( "Test 2, at" , at_f90 , dom.at , tol );
  ret += compareScalar( "Test 2, bt" , bt_f90 , dom.bt , tol );
  ret += compareScalar( "Test 2, ct" , ct_f90 , dom.ct , tol );

  // Test 3
  dom.na = 2 ; dom.nb = 3 ; dom.nc = 1;
  abcoefs    (dom);
  abcoefs_f90(dom.dt3.data(), dom.nstep, dom.na, dom.nb, dom.nc, at_f90, bt_f90, ct_f90);
  ret += compareScalar( "Test 3, at" , at_f90 , dom.at , tol );
  ret += compareScalar( "Test 3, bt" , bt_f90 , dom.bt , tol );
  ret += compareScalar( "Test 3, ct" , ct_f90 , dom.ct , tol );

  dom.nstep = 2;

  // Test 4
  dom.na = 1 ; dom.nb = 2 ; dom.nc = 3;
  abcoefs    (dom);
  abcoefs_f90(dom.dt3.data(), dom.nstep, dom.na, dom.nb, dom.nc, at_f90, bt_f90, ct_f90);
  ret += compareScalar( "Test 4, at" , at_f90 , dom.at , tol );
  ret += compareScalar( "Test 4, bt" , bt_f90 , dom.bt , tol );
  ret += compareScalar( "Test 4, ct" , ct_f90 , dom.ct , tol );

  // Test 5
  dom.na = 3 ; dom.nb = 1 ; dom.nc = 2;
  abcoefs    (dom);
  abcoefs_f90(dom.dt3.data(), dom.nstep, dom.na, dom.nb, dom.nc, at_f90, bt_f90, ct_f90);
  ret += compareScalar( "Test 5, at" , at_f90 , dom.at , tol );
  ret += compareScalar( "Test 5, bt" , bt_f90 , dom.bt , tol );
  ret += compareScalar( "Test 5, ct" , ct_f90 , dom.ct , tol );

  // Test 6
  dom.na = 2 ; dom.nb = 3 ; dom.nc = 1;
  abcoefs    (dom);
  abcoefs_f90(dom.dt3.data(), dom.nstep, dom.na, dom.nb, dom.nc, at_f90, bt_f90, ct_f90);
  ret += compareScalar( "Test 6, at" , at_f90 , dom.at , tol );
  ret += compareScalar( "Test 6, bt" , bt_f90 , dom.bt , tol );
  ret += compareScalar( "Test 6, ct" , ct_f90 , dom.ct , tol );

  dom.nstep = 3;

  // Test 7
  dom.na = 1 ; dom.nb = 2 ; dom.nc = 3;
  abcoefs    (dom);
  abcoefs_f90(dom.dt3.data(), dom.nstep, dom.na, dom.nb, dom.nc, at_f90, bt_f90, ct_f90);
  ret += compareScalar( "Test 7, at" , at_f90 , dom.at , tol );
  ret += compareScalar( "Test 7, bt" , bt_f90 , dom.bt , tol );
  ret += compareScalar( "Test 7, ct" , ct_f90 , dom.ct , tol );

  // Test 8
  dom.na = 3 ; dom.nb = 1 ; dom.nc = 2;
  abcoefs    (dom);
  abcoefs_f90(dom.dt3.data(), dom.nstep, dom.na, dom.nb, dom.nc, at_f90, bt_f90, ct_f90);
  ret += compareScalar( "Test 8, at" , at_f90 , dom.at , tol );
  ret += compareScalar( "Test 8, bt" , bt_f90 , dom.bt , tol );
  ret += compareScalar( "Test 8, ct" , ct_f90 , dom.ct , tol );

  // Test 9
  dom.na = 2 ; dom.nb = 3 ; dom.nc = 1;
  abcoefs    (dom);
  abcoefs_f90(dom.dt3.data(), dom.nstep, dom.na, dom.nb, dom.nc, at_f90, bt_f90, ct_f90);
  ret += compareScalar( "Test 9, at" , at_f90 , dom.at , tol );
  ret += compareScalar( "Test 9, bt" , bt_f90 , dom.bt , tol );
  ret += compareScalar( "Test 9, ct" , ct_f90 , dom.ct , tol );

  return ret;
}

