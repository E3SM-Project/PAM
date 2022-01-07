
#pragma once

#include "pam_const.h"

#include "YAKL.h"
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
using yakl::COLON;
using yakl::ScalarLiveOut;

#ifndef PAM_A_ORD
  #define PAM_A_ORD 5
#endif

#ifndef PAM_A_TORD
  #define PAM_A_TORD 3
#endif


int constexpr ord  = PAM_A_ORD;
int constexpr ngll = PAM_A_TORD;

static_assert(ngll <= ord , "ERROR: ngll must be <= ord");
