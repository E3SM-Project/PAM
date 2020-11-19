
#pragma once

#include "YAKL.h"
#include "yaml-cpp/yaml.h"
#include "YAKL_netcdf.h"
#include "YAKL_tridiagonal.h"

using yakl::c::parallel_for;
using yakl::c::SimpleBounds;
using yakl::c::Bounds;
using yakl::fence;
using yakl::min;
using yakl::max;
using yakl::SArray;
using yakl::memDevice;
using yakl::memHost;
using yakl::memset;
using yakl::styleC;

#ifndef ORD
  #define ORD 5
#endif

#ifndef NGLL
  #define NGLL 3
#endif

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

typedef yakl::Array<bool,1,yakl::memHost,yakl::styleC> boolHost1d;
typedef yakl::Array<bool,2,yakl::memHost,yakl::styleC> boolHost2d;
typedef yakl::Array<bool,3,yakl::memHost,yakl::styleC> boolHost3d;
typedef yakl::Array<bool,4,yakl::memHost,yakl::styleC> boolHost4d;
typedef yakl::Array<bool,5,yakl::memHost,yakl::styleC> boolHost5d;
typedef yakl::Array<bool,6,yakl::memHost,yakl::styleC> boolHost6d;
typedef yakl::Array<bool,7,yakl::memHost,yakl::styleC> boolHost7d;
typedef yakl::Array<bool,8,yakl::memHost,yakl::styleC> boolHost8d;

int constexpr ord  = ORD;
int constexpr ngll = NGLL;

static_assert(ngll <= ord , "ERROR: ngll must be <= ord");

template <class T> void endrun(T err) {
  std::cerr << err << std::endl;
  throw err;
}

