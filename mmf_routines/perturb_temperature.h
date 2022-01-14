
#pragma once

#include "pam_const.h"
#include "pam_coupler.h"


// It is best if id is globally unique for each batch of CRMs
inline void perturb_temperature( PamCoupler &coupler , int id ) {
  int constexpr num_levels = 5;
  real constexpr magnitude = 1.;

  size_t seed = static_cast<size_t>(std::clock() + id*nz*nx*ny*nens);

  // ny*nx*nens can all be globbed together for this routine
  auto temp = coupler.dm.get_lev_col<real>("temp");
  int nz   = temp.dimensions[0];
  int ncol = temp.dimensions[1];

  parallel_for( "perturb temperature" , SimpleBounds<2>(num_levels,ncol) , YAKL_LAMBDA (int k, int icol) {
    yakl::Random prng(seed+k*ncol+i);  // seed + k*ncol + i  is a globally unique identifier
    real rand = prng.genFP<real>()*2._fp - 1._fp;  // Random number in [-1,1]
    real scaling = ( num_levels - static_cast<real>(k) ) / num_levels;  // Less effect at higher levels
    temp(k,j,i,iens) += rand * magnitude * scaling;
  });
}


