
#include "TransformMatrices.h"
#include "WenoLimiter.h"


real func(real xloc) {
  return cos(2*M_PI*xloc-M_PI/10);
}


template <int ord>
void test_convergence() {
  int constexpr hs = (ord-1)/2;

  std::cout << "Order: " << ord << "\n";

  // Setup transformation matrices and WENO
  SArray<real,1,ord> gllPts_ord;
  SArray<real,1,ord> gllWts_ord;
  SArray<real,2,ord,ord> s2c;
  SArray<real,2,ord,ord> c2g;
  SArray<real,3,ord,ord,ord> wenoRecon;
  SArray<real,1,hs+2> idl;
  real sigma;

  TransformMatrices::weno_sten_to_coefs(wenoRecon);
  weno::wenoSetIdealSigma<ord>(idl,sigma);

  TransformMatrices::sten_to_coefs(s2c);
  TransformMatrices::coefs_to_gll (c2g);

  TransformMatrices::get_gll_points (gllPts_ord);
  TransformMatrices::get_gll_weights(gllWts_ord);

  /******************************************************************
   nx=5
  *******************************************************************/
  int nx1 = 5;
  // Initialize a stencil of a cosine wave centered at zero
  SArray<real,1,ord> stencil;
  real dx = 1./nx1;
  for (int i = -hs; i <= hs; i++) {
    stencil(hs+i) = 0;
    for (int ii=0; ii < ord; ii++) {
      real xloc = (i+0.5)*dx + gllPts_ord(ii)*dx;
      stencil(hs+i) += func(xloc) * gllWts_ord(ii);
    }
  }

  // Reconstruct GLL points using the stencil information
  SArray<real,1,ord> gll = c2g * s2c * stencil;

  // Integrate the error of the approximation using reconstruction
  real error1_nolim = 0;
  for (int ii=0; ii < ord; ii++) {
    real xloc = (0.5)*dx + gllPts_ord(ii)*dx;
    error1_nolim += abs(gll(ii) - func(xloc)) * gllWts_ord(ii);
  }
  error1_nolim /= nx1;

  // WENO reconstruct GLL points using the stencil information
  SArray<real,1,ord> wenoCoefs;
  weno::compute_weno_coefs<ord>( wenoRecon , stencil , wenoCoefs , idl , sigma );
  gll = c2g * wenoCoefs;

  // Integrate the error of the approximation using reconstruction
  real error1_weno = 0;
  for (int ii=0; ii < ord; ii++) {
    real xloc = (0.5)*dx + gllPts_ord(ii)*dx;
    error1_weno += abs(gll(ii) - func(xloc)) * gllWts_ord(ii);
  }
  error1_weno /= nx1;

  /******************************************************************
   nx=6
  *******************************************************************/
  int nx2 = 6;
  // Initialize a stencil of a cosine wave centered at zero
  dx = 1./nx2;
  for (int i = -hs; i <= hs; i++) {
    stencil(hs+i) = 0;
    for (int ii=0; ii < ord; ii++) {
      real xloc = (i+0.5)*dx + gllPts_ord(ii)*dx;
      stencil(hs+i) += func(xloc) * gllWts_ord(ii);
    }
  }

  // Reconstruct GLL points using the stencil information
  gll = c2g * s2c * stencil;

  // Integrate the error of the approximation using reconstruction
  real error2_nolim = 0;
  for (int ii=0; ii < ord; ii++) {
    real xloc = (0.5)*dx + gllPts_ord(ii)*dx;
    error2_nolim += abs(gll(ii) - func(xloc)) * gllWts_ord(ii);
  }
  error2_nolim /= nx2;

  // WENO reconstruct GLL points using the stencil information
  weno::compute_weno_coefs<ord>( wenoRecon , stencil , wenoCoefs , idl , sigma );
  gll = c2g * wenoCoefs;

  // Integrate the error of the approximation using reconstruction
  real error2_weno = 0;
  for (int ii=0; ii < ord; ii++) {
    real xloc = (0.5)*dx + gllPts_ord(ii)*dx;
    error2_weno += abs(gll(ii) - func(xloc)) * gllWts_ord(ii);
  }
  error2_weno /= nx2;

  real conv = log(error1_nolim/error2_nolim) / log((double) nx2/ (double) nx1);
  std::cout << "Regular: err1 err2 conv: " << error1_nolim << " " << error2_nolim << " " << conv << "\n";
  if (conv < ord-0.1) {
    std::cerr << "ERROR: Wrong convergence for regular interpolation\n";
    exit(-1);
  }

  conv = log(error1_weno/error2_weno) / log((double) nx2/ (double) nx1);
  std::cout << "WENO: err1 err2 conv: "  << error1_weno << " " << error2_weno << " " << conv << "\n";
  if (conv < ord-0.1) {
    std::cerr << "ERROR: Wrong convergence for WENO interpolation\n";
    exit(-1);
  }



  // WENO reconstruct GLL points using discontinuous data
  for (int ii=0; ii < hs; ii++) { stencil(ii) = 0; }
  stencil(hs) = 0.85;
  for (int ii=hs+1; ii < ord; ii++) { stencil(ii) = 1; }

  weno::compute_weno_coefs<ord>( wenoRecon , stencil , wenoCoefs , idl , sigma );
  gll = c2g * wenoCoefs;
  real over_weno = yakl::intrinsics::maxval(gll);

  std::cout << "GLL values (disc WENO): ";
  for (int ii=0; ii < ord; ii++) {
    std::cout << gll(ii) << "  ";
  }
  std::cout << "\n";

  gll = c2g * s2c * stencil;
  real over_nolim = yakl::intrinsics::maxval(gll);

  std::cout << "GLL values (disc nolim): ";
  for (int ii=0; ii < ord; ii++) {
    std::cout << gll(ii) << "  ";
  }
  std::cout << "\n\n";

  if (over_weno >= over_nolim) {
    std::cerr << "ERROR: WENO does not reduce the overshoot magnitude\n";
    exit(-1);
  }
}

int main() {
  test_convergence<3>();
  test_convergence<5>();
  test_convergence<7>();
  test_convergence<9>();
}

