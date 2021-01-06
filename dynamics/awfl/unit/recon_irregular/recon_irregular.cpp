
#include "TransformMatrices_variable.h"
#include "TransformMatrices.h"


real func(real xloc) {
  return cos(2*M_PI*xloc-M_PI/10);
}


template <int ord>
void test_convergence() {
  int constexpr hs = (ord-1)/2;

  std::cout << "Order: " << ord << "\n";

  SArray<real,1,ord> gllPts_ord;
  SArray<real,1,ord> gllWts_ord;
  SArray<double,2,ord,ord> s2c_tmp;
  SArray<real,2,ord,ord> s2c;
  SArray<real,2,ord,ord> c2g;

  SArray<double,1,ord> dx;
  dx(0) = 1;
  for (int ii=1; ii < ord; ii++) { dx(ii) = dx(ii-1)*1.5; }
  real norm = dx(hs);
  for (int ii=0; ii < ord; ii++) { dx(ii) /= norm; }

  SArray<double,1,ord+1> locs;
  locs(0) = 0;
  for (int ii=1; ii < ord+1; ii++) { locs(ii) = locs(ii-1) + dx(ii-1); }
  real offset = ( locs(hs) + locs(hs+1) ) / 2;
  for (int ii=0; ii < ord+1; ii++) { locs(ii) -= offset; }

  std::cout << "dx: ";
  for (int ii=0; ii < ord; ii++) {
    std::cout << dx(ii) << "  ";
  }
  std::cout << "\n";

  std::cout << "locs: ";
  for (int ii=0; ii < ord+1; ii++) {
    std::cout << locs(ii) << "  ";
  }
  std::cout << "\n";


  TransformMatrices_variable::sten_to_coefs_variable<ord>(locs,s2c_tmp);
  for (int jj=0; jj < ord; jj++) {
    for (int ii=0; ii < ord; ii++) {
      s2c(jj,ii) = s2c_tmp(jj,ii);
    }
  }
  TransformMatrices::coefs_to_gll(c2g);

  TransformMatrices::get_gll_points (gllPts_ord);
  TransformMatrices::get_gll_weights(gllWts_ord);

  /******************************************************************
   nx=20
   ******************************************************************/
  int nx1 = 10;
  // Initialize a stencil of a cosine wave centered at zero
  SArray<real,1,ord> stencil;
  real dxnorm = 1./nx1;
  for (int i = -hs; i <= hs; i++) {
    stencil(hs+i) = 0;
    for (int ii=0; ii < ord; ii++) {
      real xloc = (locs(hs+i)+dx(hs+i)/2)*dxnorm + dx(hs+i)*gllPts_ord(ii)*dxnorm;
      stencil(hs+i) += func(xloc) * gllWts_ord(ii);
    }
  }

  // Reconstruct GLL points using the stencil information
  SArray<real,1,ord> gll = c2g * s2c * stencil;

  // Integrate the error of the approximation using reconstruction
  real error1_nolim = 0;
  for (int ii=0; ii < ord; ii++) {
    real xloc = gllPts_ord(ii)*dxnorm;
    error1_nolim += abs(gll(ii) - func(xloc)) * gllWts_ord(ii);
  }
  error1_nolim /= nx1;

  /******************************************************************
   nx=30
   ******************************************************************/
  int nx2 = 30;
  // Initialize a stencil of a cosine wave centered at zero
  dxnorm = 1./nx2;
  for (int i = -hs; i <= hs; i++) {
    stencil(hs+i) = 0;
    for (int ii=0; ii < ord; ii++) {
      real xloc = (locs(hs+i)+dx(hs+i)/2)*dxnorm + dx(hs+i)*gllPts_ord(ii)*dxnorm;
      stencil(hs+i) += func(xloc) * gllWts_ord(ii);
    }
  }

  // Reconstruct GLL points using the stencil information
  gll = c2g * s2c * stencil;

  // Integrate the error of the approximation using reconstruction
  real error2_nolim = 0;
  for (int ii=0; ii < ord; ii++) {
    real xloc = gllPts_ord(ii)*dxnorm;
    error2_nolim += abs(gll(ii) - func(xloc)) * gllWts_ord(ii);
  }
  error2_nolim /= nx2;

  real conv = log(error1_nolim/error2_nolim) / log((double) nx2/ (double) nx1);
  std::cout << "Regular: err1 err2 conv: " << error1_nolim << " " << error2_nolim << " " << conv << "\n";
  if (conv < ord-1) {
    std::cerr << "ERROR: Wrong convergence for regular interpolation\n";
    exit(-1);
  }

  std::cout << "\n";

}

int main() {
  test_convergence<3>();
  test_convergence<5>();
  test_convergence<7>();
  test_convergence<9>();
}

