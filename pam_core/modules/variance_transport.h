
#pragma once

#include "pam_coupler.h"

// Enforce that this CRM's variance w/r to the column mean be forced toward the input variance
// Since the growth / decay factor will scale the same regardless, I use standard deviation to reduce computations

namespace modules {

  class VarianceTransport {
    real4d temp_factor ;
    real4d rho_v_factor;
    int    nz, ny, nx, nens;


    void init( pam::PamCoupler const &coupler ) {
      using yakl::c::parallel_for;
      using yakl::c::SimpleBounds;

      this->nz             = coupler.get_nz  ();
      this->ny             = coupler.get_ny  ();
      this->nx             = coupler.get_nx  ();
      this->nens           = coupler.get_nens();
      this->temp_factor    = real4d("temp_vt_factor" ,nz,ny,nx,nens);
      this->rho_v_factor   = real4d("rho_v_vt_factor",nz,ny,nx,nens);

      YAKL_SCOPE( temp_factor  , this->temp_factor  );
      YAKL_SCOPE( rho_v_factor , this->rho_v_factor );

      auto dm = coupler.get_data_manager_readonly();

      // Compute current variance
      auto temp  = dm.get<real const,4>("temp       ");
      auto rho_v = dm.get<real const,4>("water_vapor");

      real2d temp_var ("temp_var" ,nz,nens);
      real2d rho_v_var("rho_v_var",nz,nens);

      compute_variance( temp , rho_v , temp_var , rho_v_var );

      // Compute factor by which variance must increase
      auto input_temp_var  = dm.get<real const,2>("intput_temp_variance" );
      auto input_rho_v_var = dm.get<real const,2>("intput_rho_v_variance");

      real constexpr eps = 1.e-20;

      parallel_for( Bounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
        temp_factor (k,iens) = input_temp_var (k,iens) / ( temp_var (k,iens) + eps );
        rho_v_factor(k,iens) = input_rho_v_var(k,iens) / ( rho_v_var(k,iens) + eps );
      });

    }


    void apply_forcing( pam::PamCoupler &coupler , real dt ) {
      using yakl::c::parallel_for;
      using yakl::c::SimpleBounds;

      YAKL_SCOPE( temp_factor    , this->temp_factor    );
      YAKL_SCOPE( rho_v_factor   , this->rho_v_factor   );

      auto &dm = coupler.get_data_manager_readwrite();

      auto temp  = dm.get<real,3>("temp"       );
      auto rho_v = dm.get<real,3>("water_vapor");

      real fact = std::exp( dt / coupler.get_dt_gcm() );

      real pert_scale_min = 0.9;
      real pert_scale_max = 1.1;

      parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        temp (k,j,i,iens) *= std::max( pert_scale_min , std::min( pert_scale_max , temp_factor (k,iens)*fact ) );
        rho_v(k,j,i,iens) *= std::max( pert_scale_min , std::min( pert_scale_max , rho_v_factor(k,iens)*fact ) );
      });
    }


    void finalize( pam::PamCoupler &coupler ) {
      auto &dm = coupler.get_data_manager_readwrite();

      auto temp  = dm.get<real const,4>("temp       ");
      auto rho_v = dm.get<real const,4>("water_vapor");

      dm.register_and_allocate<real>( "temperature_variance" , "" , {nz,nens} );
      dm.register_and_allocate<real>( "water_vapor_variance" , "" , {nz,nens} );

      auto temp_var  = dm.get<real,2>("temperature_variance");
      auto rho_v_var = dm.get<real,2>("water_vapor_variance");

      compute_variance( temp , rho_v , temp_var , rho_v_var );
    }


    void compute_variance( realConst4d   temp     , realConst4d   rho_v     ,
                           real2d const &temp_var , real2d const &rho_v_var ) {
      using yakl::c::parallel_for;
      using yakl::c::SimpleBounds;

      // Initialize mean and variance to zero in the same kernel to reduce kernel count
      real2d temp_mean ("temp_mean" ,nz,nens);
      real2d rho_v_mean("rho_v_mean",nz,nens);
      real2d temp_var  ("temp_var"  ,nz,nens);
      real2d rho_v_var ("rho_v_var" ,nz,nens);

      parallel_for( Bounds<4>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
        temp_mean (k,iens) = 0;
        rho_v_mean(k,iens) = 0;
        temp_var  (k,iens) = 0;
        rho_v_var (k,iens) = 0;
      });

      real factor = 1._fp / (nx*ny);

      // Calculate the column means
      parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        yakl::atomicAdd( temp_mean (k,iens) , temp (k,j,i,iens)*factor );
        yakl::atomicAdd( rho_v_mean(k,iens) , rho_v(k,j,i,iens)*factor );
      });

      // Compute and sum the column variances
      parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        yakl::atomicAdd( temp_var (k,iens) , std::abs( temp (k,j,i,iens) - temp_mean (k,iens) )*factor );
        yakl::atomicAdd( rho_v_var(k,iens) , std::abs( rho_v(k,j,i,iens) - rho_v_mean(k,iens) )*factor );
      });
    }

  };

}


