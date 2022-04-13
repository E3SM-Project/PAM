
#pragma once

#include "YAKL.h"
#include "yaml-cpp/yaml.h"
#include <thread>

using yakl::SArray;
using yakl::memDevice;
using yakl::memHost;
using yakl::styleC;
using yakl::Array;

#ifndef PAM_ORD
  #define PAM_ORD 5
#endif

int constexpr pam_ord = PAM_ORD;

typedef double real;

int constexpr max_fields = 50;

YAKL_INLINE real constexpr operator"" _fp( long double x ) {
  return static_cast<real>(x);
}

typedef yakl::Array<real,1,yakl::memDevice,yakl::styleC> real1d;
typedef yakl::Array<real,2,yakl::memDevice,yakl::styleC> real2d;
typedef yakl::Array<real,3,yakl::memDevice,yakl::styleC> real3d;
typedef yakl::Array<real,4,yakl::memDevice,yakl::styleC> real4d;
typedef yakl::Array<real,5,yakl::memDevice,yakl::styleC> real5d;
typedef yakl::Array<real,6,yakl::memDevice,yakl::styleC> real6d;
typedef yakl::Array<real,7,yakl::memDevice,yakl::styleC> real7d;
typedef yakl::Array<real,8,yakl::memDevice,yakl::styleC> real8d;

typedef yakl::Array<int,1,yakl::memDevice,yakl::styleC> int1d;
typedef yakl::Array<int,2,yakl::memDevice,yakl::styleC> int2d;
typedef yakl::Array<int,3,yakl::memDevice,yakl::styleC> int3d;
typedef yakl::Array<int,4,yakl::memDevice,yakl::styleC> int4d;
typedef yakl::Array<int,5,yakl::memDevice,yakl::styleC> int5d;
typedef yakl::Array<int,6,yakl::memDevice,yakl::styleC> int6d;
typedef yakl::Array<int,7,yakl::memDevice,yakl::styleC> int7d;
typedef yakl::Array<int,8,yakl::memDevice,yakl::styleC> int8d;

typedef yakl::Array<bool,1,yakl::memDevice,yakl::styleC> bool1d;
typedef yakl::Array<bool,2,yakl::memDevice,yakl::styleC> bool2d;
typedef yakl::Array<bool,3,yakl::memDevice,yakl::styleC> bool3d;
typedef yakl::Array<bool,4,yakl::memDevice,yakl::styleC> bool4d;
typedef yakl::Array<bool,5,yakl::memDevice,yakl::styleC> bool5d;
typedef yakl::Array<bool,6,yakl::memDevice,yakl::styleC> bool6d;
typedef yakl::Array<bool,7,yakl::memDevice,yakl::styleC> bool7d;
typedef yakl::Array<bool,8,yakl::memDevice,yakl::styleC> bool8d;

typedef yakl::Array<real,1,yakl::memHost,yakl::styleC> realHost1d;
typedef yakl::Array<real,2,yakl::memHost,yakl::styleC> realHost2d;
typedef yakl::Array<real,3,yakl::memHost,yakl::styleC> realHost3d;
typedef yakl::Array<real,4,yakl::memHost,yakl::styleC> realHost4d;
typedef yakl::Array<real,5,yakl::memHost,yakl::styleC> realHost5d;
typedef yakl::Array<real,6,yakl::memHost,yakl::styleC> realHost6d;
typedef yakl::Array<real,7,yakl::memHost,yakl::styleC> realHost7d;
typedef yakl::Array<real,8,yakl::memHost,yakl::styleC> realHost8d;

typedef yakl::Array<int,1,yakl::memHost,yakl::styleC> intHost1d;
typedef yakl::Array<int,2,yakl::memHost,yakl::styleC> intHost2d;
typedef yakl::Array<int,3,yakl::memHost,yakl::styleC> intHost3d;
typedef yakl::Array<int,4,yakl::memHost,yakl::styleC> intHost4d;
typedef yakl::Array<int,5,yakl::memHost,yakl::styleC> intHost5d;
typedef yakl::Array<int,6,yakl::memHost,yakl::styleC> intHost6d;
typedef yakl::Array<int,7,yakl::memHost,yakl::styleC> intHost7d;
typedef yakl::Array<int,8,yakl::memHost,yakl::styleC> intHost8d;

typedef yakl::Array<bool,1,yakl::memHost,yakl::styleC> boolHost1d;
typedef yakl::Array<bool,2,yakl::memHost,yakl::styleC> boolHost2d;
typedef yakl::Array<bool,3,yakl::memHost,yakl::styleC> boolHost3d;
typedef yakl::Array<bool,4,yakl::memHost,yakl::styleC> boolHost4d;
typedef yakl::Array<bool,5,yakl::memHost,yakl::styleC> boolHost5d;
typedef yakl::Array<bool,6,yakl::memHost,yakl::styleC> boolHost6d;
typedef yakl::Array<bool,7,yakl::memHost,yakl::styleC> boolHost7d;
typedef yakl::Array<bool,8,yakl::memHost,yakl::styleC> boolHost8d;

typedef yakl::Array<real,1,yakl::memDevice,yakl::styleFortran> F_real1d;
typedef yakl::Array<real,2,yakl::memDevice,yakl::styleFortran> F_real2d;
typedef yakl::Array<real,3,yakl::memDevice,yakl::styleFortran> F_real3d;
typedef yakl::Array<real,4,yakl::memDevice,yakl::styleFortran> F_real4d;
typedef yakl::Array<real,5,yakl::memDevice,yakl::styleFortran> F_real5d;
typedef yakl::Array<real,6,yakl::memDevice,yakl::styleFortran> F_real6d;
typedef yakl::Array<real,7,yakl::memDevice,yakl::styleFortran> F_real7d;
typedef yakl::Array<real,8,yakl::memDevice,yakl::styleFortran> F_real8d;

typedef yakl::Array<int,1,yakl::memDevice,yakl::styleFortran> F_int1d;
typedef yakl::Array<int,2,yakl::memDevice,yakl::styleFortran> F_int2d;
typedef yakl::Array<int,3,yakl::memDevice,yakl::styleFortran> F_int3d;
typedef yakl::Array<int,4,yakl::memDevice,yakl::styleFortran> F_int4d;
typedef yakl::Array<int,5,yakl::memDevice,yakl::styleFortran> F_int5d;
typedef yakl::Array<int,6,yakl::memDevice,yakl::styleFortran> F_int6d;
typedef yakl::Array<int,7,yakl::memDevice,yakl::styleFortran> F_int7d;
typedef yakl::Array<int,8,yakl::memDevice,yakl::styleFortran> F_int8d;

typedef yakl::Array<bool,1,yakl::memDevice,yakl::styleFortran> F_bool1d;
typedef yakl::Array<bool,2,yakl::memDevice,yakl::styleFortran> F_bool2d;
typedef yakl::Array<bool,3,yakl::memDevice,yakl::styleFortran> F_bool3d;
typedef yakl::Array<bool,4,yakl::memDevice,yakl::styleFortran> F_bool4d;
typedef yakl::Array<bool,5,yakl::memDevice,yakl::styleFortran> F_bool5d;
typedef yakl::Array<bool,6,yakl::memDevice,yakl::styleFortran> F_bool6d;
typedef yakl::Array<bool,7,yakl::memDevice,yakl::styleFortran> F_bool7d;
typedef yakl::Array<bool,8,yakl::memDevice,yakl::styleFortran> F_bool8d;

typedef yakl::Array<real,1,yakl::memHost,yakl::styleFortran> F_realHost1d;
typedef yakl::Array<real,2,yakl::memHost,yakl::styleFortran> F_realHost2d;
typedef yakl::Array<real,3,yakl::memHost,yakl::styleFortran> F_realHost3d;
typedef yakl::Array<real,4,yakl::memHost,yakl::styleFortran> F_realHost4d;
typedef yakl::Array<real,5,yakl::memHost,yakl::styleFortran> F_realHost5d;
typedef yakl::Array<real,6,yakl::memHost,yakl::styleFortran> F_realHost6d;
typedef yakl::Array<real,7,yakl::memHost,yakl::styleFortran> F_realHost7d;
typedef yakl::Array<real,8,yakl::memHost,yakl::styleFortran> F_realHost8d;

typedef yakl::Array<int,1,yakl::memHost,yakl::styleFortran> F_intHost1d;
typedef yakl::Array<int,2,yakl::memHost,yakl::styleFortran> F_intHost2d;
typedef yakl::Array<int,3,yakl::memHost,yakl::styleFortran> F_intHost3d;
typedef yakl::Array<int,4,yakl::memHost,yakl::styleFortran> F_intHost4d;
typedef yakl::Array<int,5,yakl::memHost,yakl::styleFortran> F_intHost5d;
typedef yakl::Array<int,6,yakl::memHost,yakl::styleFortran> F_intHost6d;
typedef yakl::Array<int,7,yakl::memHost,yakl::styleFortran> F_intHost7d;
typedef yakl::Array<int,8,yakl::memHost,yakl::styleFortran> F_intHost8d;

typedef yakl::Array<bool,1,yakl::memHost,yakl::styleFortran> F_boolHost1d;
typedef yakl::Array<bool,2,yakl::memHost,yakl::styleFortran> F_boolHost2d;
typedef yakl::Array<bool,3,yakl::memHost,yakl::styleFortran> F_boolHost3d;
typedef yakl::Array<bool,4,yakl::memHost,yakl::styleFortran> F_boolHost4d;
typedef yakl::Array<bool,5,yakl::memHost,yakl::styleFortran> F_boolHost5d;
typedef yakl::Array<bool,6,yakl::memHost,yakl::styleFortran> F_boolHost6d;
typedef yakl::Array<bool,7,yakl::memHost,yakl::styleFortran> F_boolHost7d;
typedef yakl::Array<bool,8,yakl::memHost,yakl::styleFortran> F_boolHost8d;



typedef yakl::Array<real const,1,yakl::memDevice,yakl::styleC> realConst1d;
typedef yakl::Array<real const,2,yakl::memDevice,yakl::styleC> realConst2d;
typedef yakl::Array<real const,3,yakl::memDevice,yakl::styleC> realConst3d;
typedef yakl::Array<real const,4,yakl::memDevice,yakl::styleC> realConst4d;
typedef yakl::Array<real const,5,yakl::memDevice,yakl::styleC> realConst5d;
typedef yakl::Array<real const,6,yakl::memDevice,yakl::styleC> realConst6d;
typedef yakl::Array<real const,7,yakl::memDevice,yakl::styleC> realConst7d;
typedef yakl::Array<real const,8,yakl::memDevice,yakl::styleC> realConst8d;

typedef yakl::Array<int const,1,yakl::memDevice,yakl::styleC> intConst1d;
typedef yakl::Array<int const,2,yakl::memDevice,yakl::styleC> intConst2d;
typedef yakl::Array<int const,3,yakl::memDevice,yakl::styleC> intConst3d;
typedef yakl::Array<int const,4,yakl::memDevice,yakl::styleC> intConst4d;
typedef yakl::Array<int const,5,yakl::memDevice,yakl::styleC> intConst5d;
typedef yakl::Array<int const,6,yakl::memDevice,yakl::styleC> intConst6d;
typedef yakl::Array<int const,7,yakl::memDevice,yakl::styleC> intConst7d;
typedef yakl::Array<int const,8,yakl::memDevice,yakl::styleC> intConst8d;

typedef yakl::Array<bool const,1,yakl::memDevice,yakl::styleC> boolConst1d;
typedef yakl::Array<bool const,2,yakl::memDevice,yakl::styleC> boolConst2d;
typedef yakl::Array<bool const,3,yakl::memDevice,yakl::styleC> boolConst3d;
typedef yakl::Array<bool const,4,yakl::memDevice,yakl::styleC> boolConst4d;
typedef yakl::Array<bool const,5,yakl::memDevice,yakl::styleC> boolConst5d;
typedef yakl::Array<bool const,6,yakl::memDevice,yakl::styleC> boolConst6d;
typedef yakl::Array<bool const,7,yakl::memDevice,yakl::styleC> boolConst7d;
typedef yakl::Array<bool const,8,yakl::memDevice,yakl::styleC> boolConst8d;

typedef yakl::Array<real const,1,yakl::memHost,yakl::styleC> realConstHost1d;
typedef yakl::Array<real const,2,yakl::memHost,yakl::styleC> realConstHost2d;
typedef yakl::Array<real const,3,yakl::memHost,yakl::styleC> realConstHost3d;
typedef yakl::Array<real const,4,yakl::memHost,yakl::styleC> realConstHost4d;
typedef yakl::Array<real const,5,yakl::memHost,yakl::styleC> realConstHost5d;
typedef yakl::Array<real const,6,yakl::memHost,yakl::styleC> realConstHost6d;
typedef yakl::Array<real const,7,yakl::memHost,yakl::styleC> realConstHost7d;
typedef yakl::Array<real const,8,yakl::memHost,yakl::styleC> realConstHost8d;

typedef yakl::Array<int const,1,yakl::memHost,yakl::styleC> intConstHost1d;
typedef yakl::Array<int const,2,yakl::memHost,yakl::styleC> intConstHost2d;
typedef yakl::Array<int const,3,yakl::memHost,yakl::styleC> intConstHost3d;
typedef yakl::Array<int const,4,yakl::memHost,yakl::styleC> intConstHost4d;
typedef yakl::Array<int const,5,yakl::memHost,yakl::styleC> intConstHost5d;
typedef yakl::Array<int const,6,yakl::memHost,yakl::styleC> intConstHost6d;
typedef yakl::Array<int const,7,yakl::memHost,yakl::styleC> intConstHost7d;
typedef yakl::Array<int const,8,yakl::memHost,yakl::styleC> intConstHost8d;

typedef yakl::Array<bool const,1,yakl::memHost,yakl::styleC> boolConstHost1d;
typedef yakl::Array<bool const,2,yakl::memHost,yakl::styleC> boolConstHost2d;
typedef yakl::Array<bool const,3,yakl::memHost,yakl::styleC> boolConstHost3d;
typedef yakl::Array<bool const,4,yakl::memHost,yakl::styleC> boolConstHost4d;
typedef yakl::Array<bool const,5,yakl::memHost,yakl::styleC> boolConstHost5d;
typedef yakl::Array<bool const,6,yakl::memHost,yakl::styleC> boolConstHost6d;
typedef yakl::Array<bool const,7,yakl::memHost,yakl::styleC> boolConstHost7d;
typedef yakl::Array<bool const,8,yakl::memHost,yakl::styleC> boolConstHost8d;

typedef yakl::Array<real const,1,yakl::memDevice,yakl::styleFortran> F_realConst1d;
typedef yakl::Array<real const,2,yakl::memDevice,yakl::styleFortran> F_realConst2d;
typedef yakl::Array<real const,3,yakl::memDevice,yakl::styleFortran> F_realConst3d;
typedef yakl::Array<real const,4,yakl::memDevice,yakl::styleFortran> F_realConst4d;
typedef yakl::Array<real const,5,yakl::memDevice,yakl::styleFortran> F_realConst5d;
typedef yakl::Array<real const,6,yakl::memDevice,yakl::styleFortran> F_realConst6d;
typedef yakl::Array<real const,7,yakl::memDevice,yakl::styleFortran> F_realConst7d;
typedef yakl::Array<real const,8,yakl::memDevice,yakl::styleFortran> F_realConst8d;

typedef yakl::Array<int const,1,yakl::memDevice,yakl::styleFortran> F_intConst1d;
typedef yakl::Array<int const,2,yakl::memDevice,yakl::styleFortran> F_intConst2d;
typedef yakl::Array<int const,3,yakl::memDevice,yakl::styleFortran> F_intConst3d;
typedef yakl::Array<int const,4,yakl::memDevice,yakl::styleFortran> F_intConst4d;
typedef yakl::Array<int const,5,yakl::memDevice,yakl::styleFortran> F_intConst5d;
typedef yakl::Array<int const,6,yakl::memDevice,yakl::styleFortran> F_intConst6d;
typedef yakl::Array<int const,7,yakl::memDevice,yakl::styleFortran> F_intConst7d;
typedef yakl::Array<int const,8,yakl::memDevice,yakl::styleFortran> F_intConst8d;

typedef yakl::Array<bool const,1,yakl::memDevice,yakl::styleFortran> F_boolConst1d;
typedef yakl::Array<bool const,2,yakl::memDevice,yakl::styleFortran> F_boolConst2d;
typedef yakl::Array<bool const,3,yakl::memDevice,yakl::styleFortran> F_boolConst3d;
typedef yakl::Array<bool const,4,yakl::memDevice,yakl::styleFortran> F_boolConst4d;
typedef yakl::Array<bool const,5,yakl::memDevice,yakl::styleFortran> F_boolConst5d;
typedef yakl::Array<bool const,6,yakl::memDevice,yakl::styleFortran> F_boolConst6d;
typedef yakl::Array<bool const,7,yakl::memDevice,yakl::styleFortran> F_boolConst7d;
typedef yakl::Array<bool const,8,yakl::memDevice,yakl::styleFortran> F_boolConst8d;

typedef yakl::Array<real const,1,yakl::memHost,yakl::styleFortran> F_realConstHost1d;
typedef yakl::Array<real const,2,yakl::memHost,yakl::styleFortran> F_realConstHost2d;
typedef yakl::Array<real const,3,yakl::memHost,yakl::styleFortran> F_realConstHost3d;
typedef yakl::Array<real const,4,yakl::memHost,yakl::styleFortran> F_realConstHost4d;
typedef yakl::Array<real const,5,yakl::memHost,yakl::styleFortran> F_realConstHost5d;
typedef yakl::Array<real const,6,yakl::memHost,yakl::styleFortran> F_realConstHost6d;
typedef yakl::Array<real const,7,yakl::memHost,yakl::styleFortran> F_realConstHost7d;
typedef yakl::Array<real const,8,yakl::memHost,yakl::styleFortran> F_realConstHost8d;

typedef yakl::Array<int const,1,yakl::memHost,yakl::styleFortran> F_intConstHost1d;
typedef yakl::Array<int const,2,yakl::memHost,yakl::styleFortran> F_intConstHost2d;
typedef yakl::Array<int const,3,yakl::memHost,yakl::styleFortran> F_intConstHost3d;
typedef yakl::Array<int const,4,yakl::memHost,yakl::styleFortran> F_intConstHost4d;
typedef yakl::Array<int const,5,yakl::memHost,yakl::styleFortran> F_intConstHost5d;
typedef yakl::Array<int const,6,yakl::memHost,yakl::styleFortran> F_intConstHost6d;
typedef yakl::Array<int const,7,yakl::memHost,yakl::styleFortran> F_intConstHost7d;
typedef yakl::Array<int const,8,yakl::memHost,yakl::styleFortran> F_intConstHost8d;

typedef yakl::Array<bool const,1,yakl::memHost,yakl::styleFortran> F_boolConstHost1d;
typedef yakl::Array<bool const,2,yakl::memHost,yakl::styleFortran> F_boolConstHost2d;
typedef yakl::Array<bool const,3,yakl::memHost,yakl::styleFortran> F_boolConstHost3d;
typedef yakl::Array<bool const,4,yakl::memHost,yakl::styleFortran> F_boolConstHost4d;
typedef yakl::Array<bool const,5,yakl::memHost,yakl::styleFortran> F_boolConstHost5d;
typedef yakl::Array<bool const,6,yakl::memHost,yakl::styleFortran> F_boolConstHost6d;
typedef yakl::Array<bool const,7,yakl::memHost,yakl::styleFortran> F_boolConstHost7d;
typedef yakl::Array<bool const,8,yakl::memHost,yakl::styleFortran> F_boolConstHost8d;


inline void endrun(std::string err="") {
  std::cerr << err << std::endl;
  throw err;
}


template <class T, int N, int MEM, int STYLE>
inline void validate_array_nan( yakl::Array<T,N,MEM,STYLE> const &arr) {
  using yakl::c::parallel_for;
  yakl::ScalarLiveOut<bool> nan_found(false);
  parallel_for( "validate nan" , arr.totElems() , YAKL_LAMBDA (int i) {
    if ( std::isnan(arr.myData[i]) ) { nan_found = true; }
  });
  if (nan_found.hostRead()) { std::cout << "WARNING: NaN encountered in array: " << std::endl; endrun(); }
}


template <class T, int N, int MEM, int STYLE>
inline void validate_array_inf( yakl::Array<T,N,MEM,STYLE> const &arr) {
  using yakl::c::parallel_for;
  yakl::ScalarLiveOut<bool> inf_found(false);
  parallel_for( "validate inf" , arr.totElems() , YAKL_LAMBDA (int i) {
    if ( std::isinf(arr.myData[i]) ) { inf_found = true; }
  });
  if (inf_found.hostRead()) { std::cout << "WARNING: inf encountered in array" << std::endl; endrun(); }
}


template <class T, int N, int MEM, int STYLE>
inline void validate_array_inf_nan( yakl::Array<T,N,MEM,STYLE> const &arr) {
  using yakl::c::parallel_for;
  yakl::ScalarLiveOut<bool> nan_found(false);
  yakl::ScalarLiveOut<bool> inf_found(false);
  parallel_for( "validate inf nan" , arr.totElems() , YAKL_LAMBDA (int i) {
    if ( std::isnan(arr.myData[i]) ) { nan_found = true; }
    if ( std::isinf(arr.myData[i]) ) { inf_found = true; }
  });
  if (nan_found.hostRead()) { std::cout << "WARNING: NaN encountered in array" << std::endl; endrun(); }
  if (inf_found.hostRead()) { std::cout << "WARNING: inf encountered in array" << std::endl; endrun(); }
}


template <class T, int N, int MEM, int STYLE>
inline void validate_array_positive( yakl::Array<T,N,MEM,STYLE> const &arr) {
  using yakl::c::parallel_for;
  yakl::ScalarLiveOut<bool> neg_found(false);
  parallel_for( "validate positive" , arr.totElems() , YAKL_LAMBDA (int i) {
    if ( arr.myData[i] < 0. ) { neg_found = true; }
  });
  if (neg_found.hostRead()) { std::cout << "WARNING: negative value encountered in array" << yakl::intrinsics::minval(arr) << std::endl; }
}


