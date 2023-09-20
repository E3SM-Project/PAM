
#pragma once

#include "pam_const.h"
#include "pam_coupler.h"

namespace modules {

  // It is best if id is globally unique for each batch of CRMs
  inline void perturb_temperature( pam::PamCoupler &coupler , intConst1d id , real magnitude = 0.1 ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    using yakl::atomicAdd;

    int nz   = coupler.get_nz  ();
    int ny   = coupler.get_ny  ();
    int nx   = coupler.get_nx  ();
    int nens = coupler.get_nens();

    if (id.size() != nens) endrun("ERROR: size of id array must be the same as nens");

    int  num_levels = nz/4;

    // ny*nx*nens can all be globbed together for this routine
    auto &dm = coupler.get_data_manager_device_readwrite();
    auto temp = dm.get<real,4>("temp");

    real2d temp_hmean1("temp_hmean1",nz,nens);
    real2d temp_hmean2("temp_hmean1",nz,nens);
    parallel_for( YAKL_AUTO_LABEL(),
                  SimpleBounds<2>(nz,nens),
                  YAKL_LAMBDA (int k, int iens) {
      temp_hmean1(k,iens) = 0;
      temp_hmean2(k,iens) = 0;
    });

    real r_nx_ny  = 1._fp/(nx*ny);  // precompute reciprocal to avoid costly divisions
    parallel_for( YAKL_AUTO_LABEL(), SimpleBounds<4>(nz,ny,nx,nens), 
                  YAKL_LAMBDA (int k, int j, int i, int iens) {
      atomicAdd( temp_hmean1(k,iens), temp(k,j,i,iens) * r_nx_ny );
    });

    parallel_for( YAKL_AUTO_LABEL(),
                  SimpleBounds<4>(num_levels,ny,nx,nens),
                  YAKL_LAMBDA (int k, int j, int i, int iens) {
      int seed = id(iens)*num_levels*ny*nx + k*ny*nx + j*nx + i;
      yakl::Random prng(seed*1e3);
      real rand = prng.genFP<real>()*2._fp - 1._fp;  // Random number in [-1,1]
      rand = std::min( rand, 1.0 );
      rand = std::max( rand,-1.0 );
      real scaling = ( num_levels - static_cast<real>(k) ) / num_levels;  // Less effect at higher levels
      temp(k,j,i,iens) += rand * magnitude * scaling;
      atomicAdd( temp_hmean2(k,iens), temp(k,j,i,iens) * r_nx_ny );
    });

    // energy conservation scaling
    parallel_for( YAKL_AUTO_LABEL(),
                  SimpleBounds<4>(num_levels,ny,nx,nens),
                  YAKL_LAMBDA (int k, int j, int i, int iens) {
      temp(k,j,i,iens) = temp(k,j,i,iens) * temp_hmean1(k,iens) / temp_hmean2(k,iens);
    });

  }

}


