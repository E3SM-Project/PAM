
#pragma once

#include "pam_const.h"

namespace sam1mom {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;
  using yakl::fortran::Bounds;
  using yakl::intrinsics::shape;
  using yakl::intrinsics::size;
  using yakl::min;
  using yakl::max;
  using yakl::abs;

  typedef yakl::Array<real,1,yakl::memDevice,yakl::styleFortran> real1d;
  typedef yakl::Array<real,2,yakl::memDevice,yakl::styleFortran> real2d;
  typedef yakl::Array<real,3,yakl::memDevice,yakl::styleFortran> real3d;
  typedef yakl::Array<real,4,yakl::memDevice,yakl::styleFortran> real4d;
  typedef yakl::Array<real,5,yakl::memDevice,yakl::styleFortran> real5d;
  typedef yakl::Array<real,6,yakl::memDevice,yakl::styleFortran> real6d;
  typedef yakl::Array<real,7,yakl::memDevice,yakl::styleFortran> real7d;
  typedef yakl::Array<real,8,yakl::memDevice,yakl::styleFortran> real8d;

  typedef yakl::Array<int,1,yakl::memDevice,yakl::styleFortran> int1d;
  typedef yakl::Array<int,2,yakl::memDevice,yakl::styleFortran> int2d;
  typedef yakl::Array<int,3,yakl::memDevice,yakl::styleFortran> int3d;
  typedef yakl::Array<int,4,yakl::memDevice,yakl::styleFortran> int4d;
  typedef yakl::Array<int,5,yakl::memDevice,yakl::styleFortran> int5d;
  typedef yakl::Array<int,6,yakl::memDevice,yakl::styleFortran> int6d;
  typedef yakl::Array<int,7,yakl::memDevice,yakl::styleFortran> int7d;
  typedef yakl::Array<int,8,yakl::memDevice,yakl::styleFortran> int8d;
}


