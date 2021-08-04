
#include "TransformMatrices.h"
#include "WenoLimiter.h"
#include <string>


int constexpr nx = 1024*1024*10;


template <int ord>
void matvec_standard() {
  int constexpr hs = (ord-1)/2;

  SArray<real,2,ord,ord> mat;
  memset( mat , 1._fp );

  real1d arr("arr",nx+2*hs);
  memset( arr , 1._fp );

  real2d rslt("rslt",ord,nx);
  memset( arr , 1._fp );

  parallel_for( (std::string("standard")+std::to_string(ord)).c_str() , SimpleBounds<1>(nx) , YAKL_LAMBDA (int i) {
    SArray<real,1,ord> stencil;
    for (int ii=0; ii < ord; ii++) {
      stencil(ii) = arr(i+ii);
    }
    SArray<real,1,ord> result;
    for (int ii=0; ii < ord; ii++) {
      real tmp = 0;
      for (int s=0; s < ord; s++) {
        tmp += mat(ii,s) * stencil(s);
      }
      rslt(ii,i) = tmp;
    }
  });
}


template <int ord>
void matvec_nonstandard() {
  int constexpr hs = (ord-1)/2;

  SArray<real,2,ord,ord> mat;
  memset( mat , 1._fp );

  real1d arr("arr",nx+2*hs);
  memset( arr , 1._fp );

  real2d rslt("rslt",ord,nx);
  memset( arr , 1._fp );

  parallel_for( (std::string("non-standard")+std::to_string(ord)).c_str() , SimpleBounds<1>(nx) , YAKL_LAMBDA (int i) {
    SArray<real,1,ord> stencil;
    for (int ii=0; ii < ord; ii++) {
      stencil(ii) = arr(i+ii);
    }
    SArray<real,1,ord> result;
    for (int ii=0; ii < ord; ii++) {
      real tmp = 0;
      for (int s=0; s < ord; s++) {
        tmp += mat(s,ii) * stencil(s);
      }
      rslt(ii,i) = tmp;
    }
  });
}

int main() {
  yakl::init();
  {
    matvec_standard<3>();
    matvec_standard<5>();
    matvec_standard<7>();
    matvec_standard<9>();

    matvec_nonstandard<3>();
    matvec_nonstandard<5>();
    matvec_nonstandard<7>();
    matvec_nonstandard<9>();
  }
  yakl::finalize();
}

