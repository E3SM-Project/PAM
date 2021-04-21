
#pragma once

#include "YAKL.h"
#include "yaml-cpp/yaml.h"

using yakl::SArray;
using yakl::c::parallel_for;
using yakl::c::SimpleBounds;
using yakl::c::Bounds;
using yakl::memDevice;
using yakl::memHost;
using yakl::styleC;
using yakl::index_t;

typedef double real;

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


inline void endrun(std::string err="") {
  std::cerr << err << std::endl;
  throw err;
}


template <class T, int N, int MEM, int STYLE>
inline void validate_array_nan( yakl::Array<T,N,MEM,STYLE> const &arr) {
  auto arrHost = arr.createHostCopy();
  for (index_t i=0; i < arrHost.totElems(); i++) {
    if ( std::isnan(arrHost.myData[i]) ) {
      std::cout << "WARNING: NaN encountered at global index: " << i << std::endl;
    }
  }
}


template <class T, int N, int MEM, int STYLE>
inline void validate_array_inf( yakl::Array<T,N,MEM,STYLE> const &arr) {
  auto arrHost = arr.createHostCopy();
  for (index_t i=0; i < arrHost.totElems(); i++) {
    if ( std::isinf(arrHost.myData[i]) ) {
      std::cout << "WARNING: inf encountered at global index: " << i << std::endl;
    }
  }
}


template <class T, int N, int MEM, int STYLE>
inline void validate_array_inf_nan( yakl::Array<T,N,MEM,STYLE> const &arr) {
  auto arrHost = arr.createHostCopy();
  for (index_t i=0; i < arrHost.totElems(); i++) {
    if ( std::isinf(arrHost.myData[i]) ) {
      std::cout << "WARNING: inf encountered at global index: " << i << std::endl;
    }
    if ( std::isnan(arrHost.myData[i]) ) {
      std::cout << "WARNING: NaN encountered at global index: " << i << std::endl;
    }
  }
}


template <class T, int N, int MEM, int STYLE>
inline void validate_array_positive( yakl::Array<T,N,MEM,STYLE> const &arr) {
  auto arrHost = arr.createHostCopy();
  for (index_t i=0; i < arrHost.totElems(); i++) {
    if (arrHost.myData[i] < 0. ) {
      std::cout << "WARNING: negative value encountered: " << arrHost.myData[i] << " at global index: " << i << std::endl;
    }
  }
}


