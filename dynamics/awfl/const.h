
#pragma once

#include "YAKL.h"
#include "yaml-cpp/yaml.h"
#include "YAKL_netcdf.h"

using yakl::c::parallel_for;
using yakl::c::Bounds;
using yakl::fence;
using yakl::min;
using yakl::max;
using yakl::SArray;
using yakl::memDevice;
using yakl::memHost;
using yakl::memset;

#ifndef ORD
  #define ORD 5
#endif

#ifndef NGLL
  #define NGLL 3
#endif

typedef float real;

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

int constexpr ord  = ORD;
int constexpr ngll = NGLL;

static_assert(ngll <= ord , "ERROR: ngll must be <= ord");

template <class T> void endrun(T err) {
  std::cerr << err << std::endl;
  throw err;
}

