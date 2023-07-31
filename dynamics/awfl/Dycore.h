
#pragma once

#include "awfl_const.h"
#include "TransformMatrices.h"
#include "TransformMatrices_variable.h"
#include "WenoLimiter.h"
#include "idealized_profiles.h"
#include "MultipleFields.h"
#include "pam_coupler.h"


class Dycore {
  public:

  // Order of accuracy (numerical convergence for smooth flows) for the dynamical core
  #ifndef MW_ORD
    int  static constexpr ord = 5;
  #else
    int  static constexpr ord = MW_ORD;
  #endif
  // This is one value higher than normal to facilitate removing the Riemann solver phase
  int  static constexpr hs  = (ord+1)/2; // Number of halo cells ("hs" == "halo size")
  int  static constexpr num_state = 5;

  // IDs for the variables in the state vector
  int  static constexpr idR = 0;  // Density
  int  static constexpr idU = 1;  // u-momentum
  int  static constexpr idV = 2;  // v-momentum
  int  static constexpr idW = 3;  // w-momentum
  int  static constexpr idT = 4;  // Density * potential temperature

  SArray<real,3,hs,hs,hs> weno_recon_lower;



  // Compute the maximum stable time step using very conservative assumptions about max wind speed
  real compute_time_step( pam::PamCoupler const &coupler , real cfl = 0.5 ) const {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    auto nens     = coupler.get_nens();
    auto nx       = coupler.get_nx();
    auto ny       = coupler.get_ny();
    auto nz       = coupler.get_nz();
    auto dx       = coupler.get_dx();
    auto dy       = coupler.get_dy();
    auto R_d      = coupler.get_option<real>("R_d");
    auto R_v      = coupler.get_option<real>("R_v");
    auto gamma_d  = coupler.get_option<real>("gamma_d");
    auto &dm      = coupler.get_data_manager_device_readonly();
    auto dm_rho_d = dm.get<real const,4>("density_dry");
    auto dm_uvel  = dm.get<real const,4>("uvel"       );
    auto dm_vvel  = dm.get<real const,4>("vvel"       );
    auto dm_wvel  = dm.get<real const,4>("wvel"       );
    auto dm_temp  = dm.get<real const,4>("temp"       );
    auto dm_rho_v = dm.get<real const,4>("water_vapor");
    auto dz       = dm.get<real const,2>("vertical_cell_dz");
    real4d dt4d("dt4d",nz,ny,nx,nens);
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      real rho_d = dm_rho_d(k,j,i,iens);
      real u     = dm_uvel (k,j,i,iens);
      real v     = dm_vvel (k,j,i,iens);
      real w     = dm_wvel (k,j,i,iens);
      real temp  = dm_temp (k,j,i,iens);
      real rho_v = dm_rho_v(k,j,i,iens);
      real rho   = rho_d + rho_v;
      real p = ( rho_d*R_d + rho_v*R_v ) * temp;
      real cs = sqrt(gamma_d*p/rho);
      real dtx = cfl * dx         / (std::abs(u)+cs);
      real dty = cfl * dy         / (std::abs(v)+cs);
      real dtz = cfl * dz(k,iens) / (std::abs(w)+cs);
      dt4d(k,j,i,iens) = std::min( std::min( dtx , dty ) , dtz );
    });
    return yakl::intrinsics::minval( dt4d );
  }



  // Perform a single time step using SSPRK3 time stepping
  void timeStep(pam::PamCoupler &coupler) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    using yakl::intrinsics::maxval;
    using yakl::intrinsics::abs;

    auto num_tracers     = coupler.get_num_tracers();
    auto nens            = coupler.get_nens();
    auto nx              = coupler.get_nx();
    auto ny              = coupler.get_ny();
    auto nz              = coupler.get_nz();
    auto tracer_positive = coupler.get_data_manager_device_readonly().get<bool const,1>("tracer_positive");

    real dt_phys = coupler.get_option<real>("crm_dt");

    // Create arrays to hold state and tracers with halos on the left and right of the domain
    // Cells [0:hs-1] are the left halos, and cells [nx+hs:nx+2*hs-1] are the right halos
    real5d state  ("state"  ,num_state  ,nz+2*hs,ny+2*hs,nx+2*hs,nens);
    real5d tracers("tracers",num_tracers,nz+2*hs,ny+2*hs,nx+2*hs,nens);

    // Populate the state and tracers arrays using data from the coupler, convert to the dycore's desired state
    convert_coupler_to_dynamics( coupler , state , tracers );

    // Get the max stable time step for the dynamics. dt_phys might be > dt_dyn, meaning we would need to sub-cycle
    real dt_dyn = compute_time_step( coupler );

    // Get the number of sub-cycles we need, and set the dynamics time step accordingly
    int ncycles = (int) std::ceil( dt_phys / dt_dyn );
    dt_dyn = dt_phys / ncycles;

    for (int icycle = 0; icycle < ncycles; icycle++) {
      // SSPRK3 requires temporary arrays to hold intermediate state and tracers arrays
      real5d state_tmp   ("state_tmp"   ,num_state  ,nz+2*hs,ny+2*hs,nx+2*hs,nens);
      real5d state_tend  ("state_tend"  ,num_state  ,nz     ,ny     ,nx     ,nens);
      real5d tracers_tmp ("tracers_tmp" ,num_tracers,nz+2*hs,ny+2*hs,nx+2*hs,nens);
      real5d tracers_tend("tracers_tend",num_tracers,nz     ,ny     ,nx     ,nens);
      //////////////
      // Stage 1
      //////////////
      compute_tendencies( coupler , state     , state_tend , tracers     , tracers_tend , dt_dyn );
      // Apply tendencies
      parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        for (int l = 0; l < num_state  ; l++) {
          state_tmp  (l,hs+k,hs+j,hs+i,iens) = state  (l,hs+k,hs+j,hs+i,iens) + dt_dyn * state_tend  (l,k,j,i,iens);
        }
        for (int l = 0; l < num_tracers; l++) {
          tracers_tmp(l,hs+k,hs+j,hs+i,iens) = tracers(l,hs+k,hs+j,hs+i,iens) + dt_dyn * tracers_tend(l,k,j,i,iens);
          // For machine precision negative values after FCT-enforced positivity application
          if (tracer_positive(l)) {
            tracers_tmp(l,hs+k,hs+j,hs+i,iens) = std::max( 0._fp , tracers_tmp(l,hs+k,hs+j,hs+i,iens) );
          }
        }
      });
      //////////////
      // Stage 2
      //////////////
      compute_tendencies( coupler , state_tmp , state_tend , tracers_tmp , tracers_tend , (1._fp/4._fp) * dt_dyn );
      // Apply tendencies
      parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        for (int l = 0; l < num_state  ; l++) {
          state_tmp  (l,hs+k,hs+j,hs+i,iens) = (3._fp/4._fp) * state      (l,hs+k,hs+j,hs+i,iens) + 
                                               (1._fp/4._fp) * state_tmp  (l,hs+k,hs+j,hs+i,iens) +
                                               (1._fp/4._fp) * dt_dyn * state_tend  (l,k,j,i,iens);
        }
        for (int l = 0; l < num_tracers; l++) {
          tracers_tmp(l,hs+k,hs+j,hs+i,iens) = (3._fp/4._fp) * tracers    (l,hs+k,hs+j,hs+i,iens) + 
                                               (1._fp/4._fp) * tracers_tmp(l,hs+k,hs+j,hs+i,iens) +
                                               (1._fp/4._fp) * dt_dyn * tracers_tend(l,k,j,i,iens);
          // For machine precision negative values after FCT-enforced positivity application
          if (tracer_positive(l)) {
            tracers_tmp(l,hs+k,hs+j,hs+i,iens) = std::max( 0._fp , tracers_tmp(l,hs+k,hs+j,hs+i,iens) );
          }
        }
      });
      //////////////
      // Stage 3
      //////////////
      compute_tendencies( coupler , state_tmp , state_tend , tracers_tmp , tracers_tend , (2._fp/3._fp) * dt_dyn );
      // Apply tendencies
      parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        for (int l = 0; l < num_state  ; l++) {
          state      (l,hs+k,hs+j,hs+i,iens) = (1._fp/3._fp) * state      (l,hs+k,hs+j,hs+i,iens) +
                                               (2._fp/3._fp) * state_tmp  (l,hs+k,hs+j,hs+i,iens) +
                                               (2._fp/3._fp) * dt_dyn * state_tend  (l,k,j,i,iens);
        }
        for (int l = 0; l < num_tracers; l++) {
          tracers    (l,hs+k,hs+j,hs+i,iens) = (1._fp/3._fp) * tracers    (l,hs+k,hs+j,hs+i,iens) +
                                               (2._fp/3._fp) * tracers_tmp(l,hs+k,hs+j,hs+i,iens) +
                                               (2._fp/3._fp) * dt_dyn * tracers_tend(l,k,j,i,iens);
          // For machine precision negative values after FCT-enforced positivity application
          if (tracer_positive(l)) {
            tracers    (l,hs+k,hs+j,hs+i,iens) = std::max( 0._fp , tracers    (l,hs+k,hs+j,hs+i,iens) );
          }
        }
      });
    }

    // Convert the dycore's state back to the coupler's state
    convert_dynamics_to_coupler( coupler , state , tracers );
  }



  // Compute the tendencies for state and tracers for one semi-discretized step inside the RK integrator
  // Tendencies are the time rate of change for a quantity
  // Coupler is non-const because we are writing to the flux variables
  void compute_tendencies( pam::PamCoupler &coupler   ,
                           real5d const &state        ,
                           real5d const &state_tend   ,
                           real5d const &tracers      ,
                           real5d const &tracers_tend ,
                           real dt                    ) const {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    using std::min;
    using std::max;

    auto nens         = coupler.get_nens();
    auto nx           = coupler.get_nx();
    auto ny           = coupler.get_ny();
    auto nz           = coupler.get_nz();
    auto dx           = coupler.get_dx();
    auto dy           = coupler.get_dy();
    auto sim2d        = ny == 1;
    auto C0           = coupler.get_option<real>("C0"     );
    auto gamma_d      = coupler.get_option<real>("gamma_d");
    auto grav         = coupler.get_option<real>("grav"   );
    auto num_tracers  = coupler.get_num_tracers();
    auto grav_balance = coupler.get_option<bool>("balance_hydrostasis_with_gravity");

    // The store a single values flux at cell edges
    auto &dm                   = coupler.get_data_manager_device_readwrite();
    auto tracer_positive       = dm.get<bool const,1>("tracer_positive");
    auto dz                    = dm.get<real const,2>("vertical_cell_dz");
    auto vert_weno_recon_lower = dm.get<real const,5>("vert_weno_recon_lower");
    auto vert_sten_to_coefs    = dm.get<real const,4>("vert_sten_to_coefs");
    auto hy_dens_cells         = dm.get<real const,2>("hy_dens_cells");
    auto hy_pressure_cells     = dm.get<real const,2>("hy_pressure_cells");
    auto grav_var              = dm.get<real const,2>("variable_gravity");

    YAKL_SCOPE( weno_recon_lower , this->weno_recon_lower );

    // Use TransformMatrices class to create matrices & GLL points to convert degrees of freedom as needed
    SArray<real,2,ord,ord> sten_to_coefs;
    SArray<real,2,ord,2  > coefs_to_gll;
    TransformMatrices::coefs_to_gll_lower(coefs_to_gll );
    TransformMatrices::sten_to_coefs     (sten_to_coefs);
    SArray<real,1,hs+1> idl;
    real                sigma;
    weno::wenoSetIdealSigma<ord>(idl,sigma);

    real4d pressure("pressure",nz+2*hs,ny+2*hs,nx+2*hs,nens);

    // Compute pressure perturbation, density perturbation, and divide density from all other quantities before interpolation
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      if (grav_balance) {
        pressure (hs+k,hs+j,hs+i,iens) = C0*std::pow(state(idT,hs+k,hs+j,hs+i,iens),gamma_d);
      } else {
        pressure (hs+k,hs+j,hs+i,iens) = C0*std::pow(state(idT,hs+k,hs+j,hs+i,iens),gamma_d) - hy_pressure_cells(k,iens);
      }
      state(idU,hs+k,hs+j,hs+i,iens) /= state(idR,hs+k,hs+j,hs+i,iens);
      state(idV,hs+k,hs+j,hs+i,iens) /= state(idR,hs+k,hs+j,hs+i,iens);
      state(idW,hs+k,hs+j,hs+i,iens) /= state(idR,hs+k,hs+j,hs+i,iens);
      state(idT,hs+k,hs+j,hs+i,iens) /= state(idR,hs+k,hs+j,hs+i,iens);
      for (int tr=0; tr < num_tracers; tr++) { tracers(tr,hs+k,hs+j,hs+i,iens) /= state(idR,hs+k,hs+j,hs+i,iens); }
    });

    halo_exchange( coupler , state , tracers , pressure );

    real5d state_flux_x  ("state_flux_x"  ,num_state  ,nz,ny,nx+1,nens);
    real5d state_flux_y  ("state_flux_y"  ,num_state  ,nz,ny+1,nx,nens);
    real5d state_flux_z  ("state_flux_z"  ,num_state  ,nz+1,ny,nx,nens);
    real5d tracers_flux_x("tracers_flux_x",num_tracers,nz,ny,nx+1,nens);
    real5d tracers_flux_y("tracers_flux_y",num_tracers,nz,ny+1,nx,nens);
    real5d tracers_flux_z("tracers_flux_z",num_tracers,nz+1,ny,nx,nens);

    // Compute samples of state and tracers at cell edges using cell-centered reconstructions at high-order with WENO
    // At the end of this, we will have two samples per cell edge in each dimension, one from each adjacent cell.
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz+1,ny+1,nx+1,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      real constexpr cs = 350;
      ////////////////////////////////////////////////////////
      // X-direction
      ////////////////////////////////////////////////////////
      if (j < ny && k < nz) {
        SArray<real,1,ord> stencil;
        // ACOUSTIC
        real ru, pp;
        {
          // rho*u (left estimate)
          int i_upw = 0;
          for (int s=0; s < ord; s++) { stencil(s) = state(idR,hs+k,hs+j,i+i_upw+s,iens)*state(idU,hs+k,hs+j,i+i_upw+s,iens); }
          real ru_L = reconstruct(stencil,coefs_to_gll,sten_to_coefs,weno_recon_lower,idl,sigma,1-i_upw);
          // rho*u (right estimate)
          i_upw = 1;
          for (int s=0; s < ord; s++) { stencil(s) = state(idR,hs+k,hs+j,i+i_upw+s,iens)*state(idU,hs+k,hs+j,i+i_upw+s,iens); }
          real ru_R = reconstruct(stencil,coefs_to_gll,sten_to_coefs,weno_recon_lower,idl,sigma,1-i_upw);
          // pressure perturbation (left estimate)
          i_upw = 0;
          for (int s=0; s < ord; s++) { stencil(s) = pressure(hs+k,hs+j,i+i_upw+s,iens); }
          real pp_L = reconstruct(stencil,coefs_to_gll,sten_to_coefs,weno_recon_lower,idl,sigma,1-i_upw);
          // pressure perturbation (right estimate)
          i_upw = 1;
          for (int s=0; s < ord; s++) { stencil(s) = pressure(hs+k,hs+j,i+i_upw+s,iens); }
          real pp_R = reconstruct(stencil,coefs_to_gll,sten_to_coefs,weno_recon_lower,idl,sigma,1-i_upw);
          // Characteristics & upwind values
          real w1 = 0.5_fp * (pp_R-cs*ru_R);
          real w2 = 0.5_fp * (pp_L+cs*ru_L);
          pp = w1+w2;
          ru = (w2-w1)/cs;
          state_flux_x(idR,k,j,i,iens) = ru;
        }
        // ADVECTIVE
        int i_upw = ru > 0 ? 0 : 1;
        // u-velocity
        for (int s=0; s < ord; s++) { stencil(s) = state(idU,hs+k,hs+j,i+i_upw+s,iens); }
        state_flux_x(idU,k,j,i,iens) = ru * reconstruct(stencil,coefs_to_gll,sten_to_coefs,weno_recon_lower,idl,sigma,1-i_upw) + pp;
        // v-velocity
        for (int s=0; s < ord; s++) { stencil(s) = state(idV,hs+k,hs+j,i+i_upw+s,iens); }
        state_flux_x(idV,k,j,i,iens) = ru * reconstruct(stencil,coefs_to_gll,sten_to_coefs,weno_recon_lower,idl,sigma,1-i_upw);
        // w-velocity
        for (int s=0; s < ord; s++) { stencil(s) = state(idW,hs+k,hs+j,i+i_upw+s,iens); }
        state_flux_x(idW,k,j,i,iens) = ru * reconstruct(stencil,coefs_to_gll,sten_to_coefs,weno_recon_lower,idl,sigma,1-i_upw);
        // theta
        for (int s=0; s < ord; s++) { stencil(s) = state(idT,hs+k,hs+j,i+i_upw+s,iens); }
        state_flux_x(idT,k,j,i,iens) = ru * reconstruct(stencil,coefs_to_gll,sten_to_coefs,weno_recon_lower,idl,sigma,1-i_upw);
        // tracers
        for (int tr=0; tr < num_tracers; tr++) {
          for (int s=0; s < ord; s++) { stencil(s) = tracers(tr,hs+k,hs+j,i+i_upw+s,iens); }
          tracers_flux_x(tr,k,j,i,iens) = ru * reconstruct(stencil,coefs_to_gll,sten_to_coefs,weno_recon_lower,idl,sigma,1-i_upw);
        }
      }

      ////////////////////////////////////////////////////////
      // Y-direction
      ////////////////////////////////////////////////////////
      if (i < nx && k < nz) {
        if (! sim2d) {
          SArray<real,1,ord> stencil;
          // ACOUSTIC
          real rv, pp;
          {
            // rho*v (left estimate)
            int j_upw = 0;
            for (int s=0; s < ord; s++) { stencil(s) = state(idR,hs+k,j+j_upw+s,hs+i,iens)*state(idV,hs+k,j+j_upw+s,hs+i,iens); }
            real rv_L = reconstruct(stencil,coefs_to_gll,sten_to_coefs,weno_recon_lower,idl,sigma,1-j_upw);
            // rho*v (right estimate)
            j_upw = 1;
            for (int s=0; s < ord; s++) { stencil(s) = state(idR,hs+k,j+j_upw+s,hs+i,iens)*state(idV,hs+k,j+j_upw+s,hs+i,iens); }
            real rv_R = reconstruct(stencil,coefs_to_gll,sten_to_coefs,weno_recon_lower,idl,sigma,1-j_upw);
            // pressure perturbation (left estimate)
            j_upw = 0;
            for (int s=0; s < ord; s++) { stencil(s) = pressure(hs+k,j+j_upw+s,hs+i,iens); }
            real pp_L = reconstruct(stencil,coefs_to_gll,sten_to_coefs,weno_recon_lower,idl,sigma,1-j_upw);
            // pressure perturbation (right estimate)
            j_upw = 1;
            for (int s=0; s < ord; s++) { stencil(s) = pressure(hs+k,j+j_upw+s,hs+i,iens); }
            real pp_R = reconstruct(stencil,coefs_to_gll,sten_to_coefs,weno_recon_lower,idl,sigma,1-j_upw);
            // Characteristics & upwind values
            real w1 = 0.5_fp * (pp_R-cs*rv_R);
            real w2 = 0.5_fp * (pp_L+cs*rv_L);
            pp = w1+w2;
            rv = (w2-w1)/cs;
            state_flux_y(idR,k,j,i,iens) = rv;
          }
          // ADVECTIVE
          int j_upw = rv > 0 ? 0 : 1;
          // u-velocity
          for (int s=0; s < ord; s++) { stencil(s) = state(idU,hs+k,j+j_upw+s,hs+i,iens); }
          state_flux_y(idU,k,j,i,iens) = rv * reconstruct(stencil,coefs_to_gll,sten_to_coefs,weno_recon_lower,idl,sigma,1-j_upw);
          // v-velocity
          for (int s=0; s < ord; s++) { stencil(s) = state(idV,hs+k,j+j_upw+s,hs+i,iens); }
          state_flux_y(idV,k,j,i,iens) = rv * reconstruct(stencil,coefs_to_gll,sten_to_coefs,weno_recon_lower,idl,sigma,1-j_upw) + pp;
          // w-velocity
          for (int s=0; s < ord; s++) { stencil(s) = state(idW,hs+k,j+j_upw+s,hs+i,iens); }
          state_flux_y(idW,k,j,i,iens) = rv * reconstruct(stencil,coefs_to_gll,sten_to_coefs,weno_recon_lower,idl,sigma,1-j_upw);
          // theta
          for (int s=0; s < ord; s++) { stencil(s) = state(idT,hs+k,j+j_upw+s,hs+i,iens); }
          state_flux_y(idT,k,j,i,iens) = rv * reconstruct(stencil,coefs_to_gll,sten_to_coefs,weno_recon_lower,idl,sigma,1-j_upw);
          // tracers
          for (int tr=0; tr < num_tracers; tr++) {
            for (int s=0; s < ord; s++) { stencil(s) = tracers(tr,hs+k,j+j_upw+s,hs+i,iens); }
            tracers_flux_y(tr,k,j,i,iens) = rv * reconstruct(stencil,coefs_to_gll,sten_to_coefs,weno_recon_lower,idl,sigma,1-j_upw);
          }
        } else {
          state_flux_y(idR,k,j,i,iens) = 0;
          state_flux_y(idU,k,j,i,iens) = 0;
          state_flux_y(idV,k,j,i,iens) = 0;
          state_flux_y(idW,k,j,i,iens) = 0;
          state_flux_y(idT,k,j,i,iens) = 0;
          for (int tr=0; tr < num_tracers; tr++) { tracers_flux_y(tr,k,j,i,iens) = 0; }
        }
      }

      ////////////////////////////////////////////////////////
      // Z-direction
      ////////////////////////////////////////////////////////
      if (i < nx && j < ny) {
        SArray<real,1,ord> stencil;
        SArray<real,2,ord,ord>  s2c_loc[2];
        SArray<real,3,hs,hs,hs> wrl_loc[2];
        for (int i1=0; i1 < ord; i1++) {
          for (int i2=0; i2 < ord; i2++) {
            s2c_loc[0](i1,i2) = vert_sten_to_coefs(k  ,i1,i2,iens);
            s2c_loc[1](i1,i2) = vert_sten_to_coefs(k+1,i1,i2,iens);
          }
        }
        for (int i1=0; i1 < hs; i1++) {
          for (int i2=0; i2 < hs; i2++) {
            for (int i3=0; i3 < hs; i3++) {
              wrl_loc[0](i1,i2,i3) = vert_weno_recon_lower(k  ,i1,i2,i3,iens);
              wrl_loc[1](i1,i2,i3) = vert_weno_recon_lower(k+1,i1,i2,i3,iens);
            }
          }
        }
        // ACOUSTIC
        real rw, pp;
        {
          // rho*w (left estimate)
          int k_upw = 0;
          for (int s=0; s < ord; s++) { stencil(s) = state(idR,k+k_upw+s,hs+j,hs+i,iens)*state(idW,k+k_upw+s,hs+j,hs+i,iens); }
          real rw_L = reconstruct(stencil,coefs_to_gll,s2c_loc[k_upw],wrl_loc[k_upw],idl,sigma,1-k_upw);
          if ((k == 0 || k == nz)) rw_L = 0;
          // rho*w (right estimate)
          k_upw = 1;
          for (int s=0; s < ord; s++) { stencil(s) = state(idR,k+k_upw+s,hs+j,hs+i,iens)*state(idW,k+k_upw+s,hs+j,hs+i,iens); }
          real rw_R = reconstruct(stencil,coefs_to_gll,s2c_loc[k_upw],wrl_loc[k_upw],idl,sigma,1-k_upw);
          if ((k == 0 || k == nz)) rw_R = 0;
          // pressure perturbation (left estimate)
          k_upw = 0;
          for (int s=0; s < ord; s++) { stencil(s) = pressure(k+k_upw+s,hs+j,hs+i,iens); }
          real pp_L = reconstruct(stencil,coefs_to_gll,s2c_loc[k_upw],wrl_loc[k_upw],idl,sigma,1-k_upw);
          // pressure perturbation (right estimate)
          k_upw = 1;
          for (int s=0; s < ord; s++) { stencil(s) = pressure(k+k_upw+s,hs+j,hs+i,iens); }
          real pp_R = reconstruct(stencil,coefs_to_gll,s2c_loc[k_upw],wrl_loc[k_upw],idl,sigma,1-k_upw);
          // Characteristics & upwind values
          real w1 = 0.5_fp * (pp_R-cs*rw_R);
          real w2 = 0.5_fp * (pp_L+cs*rw_L);
          pp = w1+w2;
          rw = (w2-w1)/cs;
          if ((k == 0 || k == nz)) rw = 0;
          state_flux_z(idR,k,j,i,iens) = rw;
        }
        // ADVECTIVE
        int k_upw = rw > 0 ? 0 : 1;
        // u-velocity
        for (int s=0; s < ord; s++) { stencil(s) = state(idU,k+k_upw+s,hs+j,hs+i,iens); }
        state_flux_z(idU,k,j,i,iens) = rw * reconstruct(stencil,coefs_to_gll,s2c_loc[k_upw],wrl_loc[k_upw],idl,sigma,1-k_upw);
        // v-velocity
        for (int s=0; s < ord; s++) { stencil(s) = state(idV,k+k_upw+s,hs+j,hs+i,iens); }
        state_flux_z(idV,k,j,i,iens) = rw * reconstruct(stencil,coefs_to_gll,s2c_loc[k_upw],wrl_loc[k_upw],idl,sigma,1-k_upw);
        // w-velocity
        for (int s=0; s < ord; s++) { stencil(s) = state(idW,k+k_upw+s,hs+j,hs+i,iens); }
        state_flux_z(idW,k,j,i,iens) = rw * reconstruct(stencil,coefs_to_gll,s2c_loc[k_upw],wrl_loc[k_upw],idl,sigma,1-k_upw) + pp;
        // theta
        for (int s=0; s < ord; s++) { stencil(s) = state(idT,k+k_upw+s,hs+j,hs+i,iens); }
        state_flux_z(idT,k,j,i,iens) = rw * reconstruct(stencil,coefs_to_gll,s2c_loc[k_upw],wrl_loc[k_upw],idl,sigma,1-k_upw);
        // tracers
        for (int tr=0; tr < num_tracers; tr++) {
          for (int s=0; s < ord; s++) { stencil(s) = tracers(tr,k+k_upw+s,hs+j,hs+i,iens); }
          tracers_flux_z(tr,k,j,i,iens) = rw * reconstruct(stencil,coefs_to_gll,s2c_loc[k_upw],wrl_loc[k_upw],idl,sigma,1-k_upw);
        }
      }
    });

    // Flux Corrected Transport to enforce positivity for tracer species that must remain non-negative
    // This looks like it has a race condition, but it does not. Only one of the adjacent cells can ever change
    // a given edge flux because it's only changed if its sign oriented outward from a cell.
    // Also, multiply density back onto the state and tracers
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) ,
                                      YAKL_LAMBDA (int k, int j, int i, int iens) {
      state(idU,hs+k,hs+j,hs+i,iens) *= state(idR,hs+k,hs+j,hs+i,iens);
      state(idV,hs+k,hs+j,hs+i,iens) *= state(idR,hs+k,hs+j,hs+i,iens);
      state(idW,hs+k,hs+j,hs+i,iens) *= state(idR,hs+k,hs+j,hs+i,iens);
      state(idT,hs+k,hs+j,hs+i,iens) *= state(idR,hs+k,hs+j,hs+i,iens);
      for (int tr=0; tr < num_tracers; tr++) {
        tracers(tr,hs+k,hs+j,hs+i,iens) *= state(idR,hs+k,hs+j,hs+i,iens);
        if (tracer_positive(tr)) {
          real mass_available = max(tracers(tr,hs+k,hs+j,hs+i,iens),0._fp) * dx * dy * dz(k,iens);
          real flux_out_x = ( max(tracers_flux_x(tr,k,j,i+1,iens),0._fp) - min(tracers_flux_x(tr,k,j,i,iens),0._fp) ) / dx;
          real flux_out_y = ( max(tracers_flux_y(tr,k,j+1,i,iens),0._fp) - min(tracers_flux_y(tr,k,j,i,iens),0._fp) ) / dy;
          real flux_out_z = ( max(tracers_flux_z(tr,k+1,j,i,iens),0._fp) - min(tracers_flux_z(tr,k,j,i,iens),0._fp) ) / dz(k,iens);
          real mass_out = (flux_out_x + flux_out_y + flux_out_z) * dt * dx * dy * dz(k,iens);
          if (mass_out > mass_available) {
            real mult = mass_available / mass_out;
            if (tracers_flux_x(tr,k,j,i+1,iens) > 0) tracers_flux_x(tr,k,j,i+1,iens) *= mult;
            if (tracers_flux_x(tr,k,j,i  ,iens) < 0) tracers_flux_x(tr,k,j,i  ,iens) *= mult;
            if (tracers_flux_y(tr,k,j+1,i,iens) > 0) tracers_flux_y(tr,k,j+1,i,iens) *= mult;
            if (tracers_flux_y(tr,k,j  ,i,iens) < 0) tracers_flux_y(tr,k,j  ,i,iens) *= mult;
            if (tracers_flux_z(tr,k+1,j,i,iens) > 0) tracers_flux_z(tr,k+1,j,i,iens) *= mult;
            if (tracers_flux_z(tr,k  ,j,i,iens) < 0) tracers_flux_z(tr,k  ,j,i,iens) *= mult;
          }
        }
      }
    });

    // Compute tendencies as the flux divergence + gravity source term
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      for (int l = 0; l < num_state; l++) {
        state_tend  (l,k,j,i,iens) = -( state_flux_x  (l,k  ,j  ,i+1,iens) - state_flux_x  (l,k,j,i,iens) ) / dx
                                     -( state_flux_y  (l,k  ,j+1,i  ,iens) - state_flux_y  (l,k,j,i,iens) ) / dy
                                     -( state_flux_z  (l,k+1,j  ,i  ,iens) - state_flux_z  (l,k,j,i,iens) ) / dz(k,iens);
        if (l == idW) {
          if (grav_balance) {
            state_tend(l,k,j,i,iens) += -grav_var(k,iens) * state(idR,hs+k,hs+j,hs+i,iens);
          } else {
            state_tend(l,k,j,i,iens) += -grav * ( state(idR,hs+k,hs+j,hs+i,iens) - hy_dens_cells(k,iens) );
          }
        }
        if (l == idV && sim2d) state_tend(l,k,j,i,iens) = 0;
      }
      for (int l = 0; l < num_tracers; l++) {
        tracers_tend(l,k,j,i,iens) = -( tracers_flux_x(l,k  ,j  ,i+1,iens) - tracers_flux_x(l,k,j,i,iens) ) / dx
                                     -( tracers_flux_y(l,k  ,j+1,i  ,iens) - tracers_flux_y(l,k,j,i,iens) ) / dy
                                     -( tracers_flux_z(l,k+1,j  ,i  ,iens) - tracers_flux_z(l,k,j,i,iens) ) / dz(k,iens);
      }
    });

  }



  // ord stencil cell averages to two GLL point values via high-order reconstruction and WENO limiting
  YAKL_INLINE static real reconstruct( SArray<real,1,ord>      const &stencil          ,
                                       SArray<real,2,ord,2>    const &coefs_to_gll     ,
                                       SArray<real,2,ord,ord>  const &sten_to_coefs    ,
                                       SArray<real,3,hs,hs,hs> const &weno_recon_lower ,
                                       SArray<real,1,hs+1>     const &idl              ,
                                       real                          sigma             ,
                                       int                           ind               ) {
    // Reconstruct values
    SArray<real,1,ord> wenoCoefs;
    weno::compute_weno_coefs<ord>( weno_recon_lower , sten_to_coefs , stencil , wenoCoefs , idl , sigma );
    real tmp = 0;
    for (int s=0; s < ord; s++) { tmp += coefs_to_gll(s,ind) * wenoCoefs(s); }
    return tmp;
  }



  void halo_exchange( pam::PamCoupler const &coupler  ,
                      real5d          const &state    ,
                      real5d          const &tracers  ,
                      real4d          const &pressure ) const {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;

    auto nens         = coupler.get_nens();
    auto nx           = coupler.get_nx();
    auto ny           = coupler.get_ny();
    auto nz           = coupler.get_nz();
    auto num_tracers  = coupler.get_num_tracers();
    auto sim2d        = ny == 1;
    auto gamma        = coupler.get_option<real>("gamma_d");
    auto C0           = coupler.get_option<real>("C0");
    auto grav         = coupler.get_option<real>("grav");
    auto grav_balance = coupler.get_option<bool>("balance_hydrostasis_with_gravity");
    auto dz = coupler.get_data_manager_device_readonly().get<real const,2>("vertical_cell_dz");

    int npack = num_state + num_tracers + 1;

    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<5>(npack,nz,ny,hs,nens) ,
                                      YAKL_LAMBDA (int v, int k, int j, int ii, int iens) {
      if (v < num_state) {
        state  (v          ,hs+k,hs+j,nx+hs+ii,iens) = state  (v          ,hs+k,hs+j,hs+ii,iens);
        state  (v          ,hs+k,hs+j,      ii,iens) = state  (v          ,hs+k,hs+j,nx+ii,iens);
      } else if (v < num_state + num_tracers) {
        tracers(v-num_state,hs+k,hs+j,nx+hs+ii,iens) = tracers(v-num_state,hs+k,hs+j,hs+ii,iens);
        tracers(v-num_state,hs+k,hs+j,      ii,iens) = tracers(v-num_state,hs+k,hs+j,nx+ii,iens);
      } else {
        pressure(hs+k,hs+j,nx+hs+ii,iens) = pressure(hs+k,hs+j,hs+ii,iens);
        pressure(hs+k,hs+j,      ii,iens) = pressure(hs+k,hs+j,nx+ii,iens);
      }
    });

    if (!sim2d) {
      parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<5>(npack,nz,hs,nx,nens) ,
                                        YAKL_LAMBDA (int v, int k, int jj, int i, int iens) {
        if (v < num_state) {
          state  (v          ,hs+k,ny+hs+jj,hs+i,iens) = state  (v          ,hs+k,hs+jj,hs+i,iens);
          state  (v          ,hs+k,      jj,hs+i,iens) = state  (v          ,hs+k,ny+jj,hs+i,iens);
        } else if (v < num_state + num_tracers) {
          tracers(v-num_state,hs+k,ny+hs+jj,hs+i,iens) = tracers(v-num_state,hs+k,hs+jj,hs+i,iens);
          tracers(v-num_state,hs+k,      jj,hs+i,iens) = tracers(v-num_state,hs+k,ny+jj,hs+i,iens);
        } else {
          pressure(hs+k,ny+hs+jj,hs+i,iens) = pressure(hs+k,hs+jj,hs+i,iens);
          pressure(hs+k,      jj,hs+i,iens) = pressure(hs+k,ny+jj,hs+i,iens);
        }
      });
    }

    ////////////////////////////////////
    // Begin boundary conditions
    ////////////////////////////////////
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(hs,ny,nx,nens) ,
                                      YAKL_LAMBDA (int kk, int j, int i, int iens) {
      // Exclude density (idR == 0) because that's calculated lower down
      for (int l=1; l < num_state; l++) {
        if (l == idW) {
          state(l,      kk,hs+j,hs+i,iens) = 0;
          state(l,hs+nz+kk,hs+j,hs+i,iens) = 0;
        } else {
          state(l,      kk,hs+j,hs+i,iens) = state(l,hs+0   ,hs+j,hs+i,iens);
          state(l,hs+nz+kk,hs+j,hs+i,iens) = state(l,hs+nz-1,hs+j,hs+i,iens);
        }
      }
      for (int l=0; l < num_tracers; l++) {
        tracers(l,      kk,hs+j,hs+i,iens) = tracers(l,hs+0   ,hs+j,hs+i,iens);
        tracers(l,hs+nz+kk,hs+j,hs+i,iens) = tracers(l,hs+nz-1,hs+j,hs+i,iens);
      }
      if (! grav_balance) {
        pressure(      kk,hs+j,hs+i,iens) = pressure(hs+0   ,hs+j,hs+i,iens);
        pressure(hs+nz+kk,hs+j,hs+i,iens) = pressure(hs+nz-1,hs+j,hs+i,iens);
      }
      {
        int  k0       = hs;
        int  k        = k0-1-kk;
        real rho0     = state(idR,k0,hs+j,hs+i,iens);
        real theta0   = state(idT,k0,hs+j,hs+i,iens);
        real rho0_gm1 = std::pow(rho0  ,gamma-1);
        real theta0_g = std::pow(theta0,gamma  );
        state(idR,k,hs+j,hs+i,iens) = std::pow( rho0_gm1 + grav*(gamma-1)*dz(k0-hs,iens)*(kk+1)/(gamma*C0*theta0_g) ,
                                                1._fp/(gamma-1) );
        if (grav_balance) {
          real rt = state(idR,k,hs+j,hs+i,iens) * state(idT,k,hs+j,hs+i,iens);
          pressure(k,hs+j,hs+i,iens) = C0*std::pow( rt , gamma );
        }
      }
      {
        int  k0       = hs+nz-1;
        int  k        = k0+1+kk;
        real rho0     = state(idR,k0,hs+j,hs+i,iens);
        real theta0   = state(idT,k0,hs+j,hs+i,iens);
        real rho0_gm1 = std::pow(rho0  ,gamma-1);
        real theta0_g = std::pow(theta0,gamma  );
        state(idR,k,hs+j,hs+i,iens) = std::pow( rho0_gm1 - grav*(gamma-1)*dz(k0-hs,iens)*(kk+1)/(gamma*C0*theta0_g) ,
                                                1._fp/(gamma-1) );
        if (grav_balance) {
          real rt = state(idR,k,hs+j,hs+i,iens) * state(idT,k,hs+j,hs+i,iens);
          pressure(k,hs+j,hs+i,iens) = C0*std::pow( rt , gamma );
        }
      }
    });
  }



  // Creates initial data at a point in space for the rising moist thermal test case
  YAKL_INLINE static void thermal(real x, real y, real z, real xlen, real ylen, real grav, real C0, real gamma,
                                  real cp, real p0, real R_d, real R_v, real &rho, real &u, real &v, real &w,
                                  real &theta, real &rho_v, real &hr, real &ht) {
    hydro_const_theta(z,grav,C0,cp,p0,gamma,R_d,hr,ht);
    real rho_d   = hr;
    u            = 0.;
    v            = 0.;
    w            = 0.;
    real theta_d = ht + sample_ellipse_cosine(2._fp  ,  x,y,z  ,  xlen/2,ylen/2,2000.  ,  2000.,2000.,2000.);
    real p_d     = C0 * pow( rho_d*theta_d , gamma );
    real temp    = p_d / rho_d / R_d;
    real sat_pv  = saturation_vapor_pressure(temp);
    real sat_rv  = sat_pv / R_v / temp;
    rho_v        = sample_ellipse_cosine(0.8_fp  ,  x,y,z  ,  xlen/2,ylen/2,2000.  ,  2000.,2000.,2000.) * sat_rv;
    real p       = rho_d * R_d * temp + rho_v * R_v * temp;
    rho          = rho_d + rho_v;
    theta        = std::pow( p / C0 , 1._fp / gamma ) / rho;
  }



  // Computes a hydrostatic background density and potential temperature using c constant potential temperature
  // backgrounda for a single vertical location
  YAKL_INLINE static void hydro_const_theta( real z, real grav, real C0, real cp, real p0, real gamma, real rd,
                                             real &r, real &t ) {
    const real theta0 = 300.;  //Background potential temperature
    const real exner0 = 1.;    //Surface-level Exner pressure
    t = theta0;                                       //Potential Temperature at z
    real exner = exner0 - grav * z / (cp * theta0);   //Exner pressure at z
    real p = p0 * std::pow(exner,(cp/rd));            //Pressure at z
    real rt = std::pow((p / C0),(1._fp / gamma));     //rho*theta at z
    r = rt / t;                                       //Density at z
  }



  // Samples a 3-D ellipsoid at a point in space
  YAKL_INLINE static real sample_ellipse_cosine(real amp, real x   , real y   , real z   ,
                                                          real x0  , real y0  , real z0  ,
                                                          real xrad, real yrad, real zrad) {
    //Compute distance from bubble center
    real dist = sqrt( ((x-x0)/xrad)*((x-x0)/xrad) +
                      ((y-y0)/yrad)*((y-y0)/yrad) +
                      ((z-z0)/zrad)*((z-z0)/zrad) ) * M_PI / 2.;
    //If the distance from bubble center is less than the radius, create a cos**2 profile
    if (dist <= M_PI / 2.) {
      return amp * std::pow(cos(dist),2._fp);
    } else {
      return 0.;
    }
  }



  YAKL_INLINE static real saturation_vapor_pressure(real temp) {
    real tc = temp - 273.15;
    return 610.94 * std::exp( 17.625*tc / (243.04+tc) );
  }



  // Compute supercell temperature profile at a vertical location
  YAKL_INLINE static real init_supercell_temperature(real z, real z_0, real z_trop, real z_top,
                                                     real T_0, real T_trop, real T_top) {
    if (z <= z_trop) {
      real lapse = - (T_trop - T_0) / (z_trop - z_0);
      return T_0 - lapse * (z - z_0);
    } else {
      real lapse = - (T_top - T_trop) / (z_top - z_trop);
      return T_trop - lapse * (z - z_trop);
    }
  }



  // Compute supercell dry pressure profile at a vertical location
  YAKL_INLINE static real init_supercell_pressure_dry(real z, real z_0, real z_trop, real z_top,
                                                      real T_0, real T_trop, real T_top,
                                                      real p_0, real R_d, real grav) {
    if (z <= z_trop) {
      real lapse = - (T_trop - T_0) / (z_trop - z_0);
      real T = init_supercell_temperature(z, z_0, z_trop, z_top, T_0, T_trop, T_top);
      return p_0 * pow( T / T_0 , grav/(R_d*lapse) );
    } else {
      // Get pressure at the tropopause
      real lapse = - (T_trop - T_0) / (z_trop - z_0);
      real p_trop = p_0 * pow( T_trop / T_0 , grav/(R_d*lapse) );
      // Get pressure at requested height
      lapse = - (T_top - T_trop) / (z_top - z_trop);
      if (lapse != 0) {
        real T = init_supercell_temperature(z, z_0, z_trop, z_top, T_0, T_trop, T_top);
        return p_trop * pow( T / T_trop , grav/(R_d*lapse) );
      } else {
        return p_trop * exp(-grav*(z-z_trop)/(R_d*T_trop));
      }
    }
  }


  
  // Compute supercell relative humidity profile at a vertical location
  YAKL_INLINE static real init_supercell_relhum(real z, real z_0, real z_trop) {
    if (z <= z_trop) {
      return 1._fp - 0.75_fp * pow(z / z_trop , 1.25_fp );
    } else {
      return 0.25_fp;
    }
  }



  // Computes dry saturation mixing ratio
  YAKL_INLINE static real init_supercell_sat_mix_dry( real press , real T ) {
    return 380/(press) * exp( 17.27_fp * (T-273)/(T-36) );
  }



  // Initialize the class data as well as the state and tracers arrays and convert them back into the coupler state
  void init(pam::PamCoupler &coupler, bool verbose=false) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;

    auto nens        = coupler.get_nens();
    auto nx          = coupler.get_nx();
    auto ny          = coupler.get_ny();
    auto nz          = coupler.get_nz();
    auto dx          = coupler.get_dx();
    auto dy          = coupler.get_dy();
    auto xlen        = coupler.get_xlen();
    auto ylen        = coupler.get_ylen();
    auto sim2d       = ny == 1;
    auto num_tracers = coupler.get_num_tracers();

    coupler.set_option<bool>("balance_hydrostasis_with_gravity",true);
    if (coupler.get_option<bool>("balance_hydrostasis_with_gravity")) {
      coupler.get_data_manager_device_readwrite().register_and_allocate<real>("variable_gravity","",{nz,nens});
    }

    if (! coupler.option_exists("R_d"     )) coupler.set_option<real>("R_d"     ,287.       );
    if (! coupler.option_exists("cp_d"    )) coupler.set_option<real>("cp_d"    ,1003.      );
    if (! coupler.option_exists("R_v"     )) coupler.set_option<real>("R_v"     ,461.       );
    if (! coupler.option_exists("cp_v"    )) coupler.set_option<real>("cp_v"    ,1859       );
    if (! coupler.option_exists("p0"      )) coupler.set_option<real>("p0"      ,1.e5       );
    if (! coupler.option_exists("grav"    )) coupler.set_option<real>("grav"    ,9.81       );
    auto R_d  = coupler.get_option<real>("R_d" );
    auto cp_d = coupler.get_option<real>("cp_d");
    auto R_v  = coupler.get_option<real>("R_v" );
    auto cp_v = coupler.get_option<real>("cp_v");
    auto p0   = coupler.get_option<real>("p0"  );
    auto grav = coupler.get_option<real>("grav");
    if (! coupler.option_exists("cv_d"   )) coupler.set_option<real>("cv_d"   ,cp_d - R_d );
    auto cv_d = coupler.get_option<real>("cv_d");
    if (! coupler.option_exists("gamma_d")) coupler.set_option<real>("gamma_d",cp_d / cv_d);
    if (! coupler.option_exists("kappa_d")) coupler.set_option<real>("kappa_d",R_d  / cp_d);
    if (! coupler.option_exists("cv_v"   )) coupler.set_option<real>("cv_v"   ,R_v - cp_v );
    auto gamma = coupler.get_option<real>("gamma_d");
    auto kappa = coupler.get_option<real>("kappa_d");
    if (! coupler.option_exists("C0")) coupler.set_option<real>("C0" , pow( R_d * pow( p0 , -kappa ) , gamma ));
    auto C0    = coupler.get_option<real>("C0");

    auto &dm = coupler.get_data_manager_device_readwrite();
    auto dz = dm.get<real const,2>("vertical_cell_dz");

    // Create WENO reconstruction matrices
    dm.register_and_allocate<real>("vert_weno_recon_lower","",{nz+2,hs,hs,hs,nens});
    dm.register_and_allocate<real>("vert_sten_to_coefs"   ,"",{nz+2,ord,ord ,nens});
    auto vert_weno_recon_lower = dm.get<real,5>("vert_weno_recon_lower");
    auto vert_sten_to_coefs    = dm.get<real,4>("vert_sten_to_coefs"   );
    auto dz_host        = dz                   .createHostCopy();
    auto vert_weno_host = vert_weno_recon_lower.createHostCopy();
    auto vert_s2c_host  = vert_sten_to_coefs   .createHostCopy();
    for (int k=0; k < nz+2; k++) {
      for (int iens = 0; iens < nens; iens++) {
        // Load normalize grid spacings
        SArray<real,1,ord> dzloc;
        for (int kk=0; kk < ord; kk++) {
          int ind1 = std::min(nz-1,std::max(0,-1+k+kk));
          int ind2 = std::min(nz-1,std::max(0,-1+k   ));
          dzloc(kk) = dz_host(ind1,iens) / dz_host(ind2,iens);
        }
        // Compute normalized locations of cell edges
        SArray<real,1,ord+1> locs;
        locs(0) = 0;
        for (int kk=1; kk < ord+1; kk++) { locs(kk) = locs(kk-1) + dzloc(kk-1); }
        real midloc = ( locs((ord-1)/2) + locs((ord+1)/2) ) / 2;
        for (int kk=0; kk < ord+1; kk++) { locs(kk) = locs(kk) - midloc; }
        // Compute WENO reconstruction matrices for each level
        SArray<double,2,ord,ord>  s2c_var;
        SArray<double,3,hs,hs,hs> weno_recon_lower_var;
        TransformMatrices_variable::sten_to_coefs_variable  <ord>(locs,s2c_var             );
        TransformMatrices_variable::weno_lower_sten_to_coefs<ord>(locs,weno_recon_lower_var);
        for (int jj=0; jj < ord; jj++) {
          for (int ii=0; ii < ord; ii++) {
            vert_s2c_host(k,jj,ii,iens) = s2c_var(jj,ii);
          }
        }
        for (int kk=0; kk < hs; kk++) {
          for (int jj=0; jj < hs; jj++) {
            for (int ii=0; ii < hs; ii++) {
              vert_weno_host(k,kk,jj,ii,iens) = weno_recon_lower_var(kk,jj,ii);
            }
          }
        }
      }
    }
    vert_s2c_host .deep_copy_to(vert_sten_to_coefs   );
    vert_weno_host.deep_copy_to(vert_weno_recon_lower);
    TransformMatrices::weno_lower_sten_to_coefs(this->weno_recon_lower);

    R_d   = coupler.get_option<real>("R_d"    );
    R_v   = coupler.get_option<real>("R_v"    );
    cp_d  = coupler.get_option<real>("cp_d"   );
    cp_v  = coupler.get_option<real>("cp_v"   );
    p0    = coupler.get_option<real>("p0"     );
    grav  = coupler.get_option<real>("grav"   );
    kappa = coupler.get_option<real>("kappa_d");
    gamma = coupler.get_option<real>("gamma_d");
    C0    = coupler.get_option<real>("C0"     );

    // Create arrays to determine whether we should add mass for a tracer or whether it should remain non-negative
    num_tracers = coupler.get_num_tracers();
    bool1d tracer_adds_mass("tracer_adds_mass",num_tracers);
    bool1d tracer_positive ("tracer_positive" ,num_tracers);

    // Must assign on the host to avoid segfaults
    auto tracer_adds_mass_host = tracer_adds_mass.createHostCopy();
    auto tracer_positive_host  = tracer_positive .createHostCopy();

    auto tracer_names = coupler.get_tracer_names();  // Get a list of tracer names
    int idWV;
    for (int tr=0; tr < num_tracers; tr++) {
      std::string tracer_desc;
      bool        tracer_found, positive, adds_mass;
      coupler.get_tracer_info( tracer_names[tr] , tracer_desc, tracer_found , positive , adds_mass);
      tracer_positive_host (tr) = positive;
      tracer_adds_mass_host(tr) = adds_mass;
      if (tracer_names[tr] == "water_vapor") idWV = tr;  // Be sure to track which index belongs to water vapor
    }
    tracer_positive_host .deep_copy_to(tracer_positive );
    tracer_adds_mass_host.deep_copy_to(tracer_adds_mass);

    coupler.set_option<int>("idWV",idWV);
    dm.register_and_allocate<bool>("tracer_adds_mass","",{num_tracers});
    auto dm_tracer_adds_mass = dm.get<bool,1>("tracer_adds_mass");
    tracer_adds_mass.deep_copy_to(dm_tracer_adds_mass);

    dm.register_and_allocate<bool>("tracer_positive","",{num_tracers});
    auto dm_tracer_positive = dm.get<bool,1>("tracer_positive");
    tracer_positive.deep_copy_to(dm_tracer_positive);

    dm.register_and_allocate<real>("hy_dens_cells"    ,"",{nz,nens});
    dm.register_and_allocate<real>("hy_pressure_cells","",{nz,nens});

    int static constexpr DATA_SPEC_THERMAL   = 0;
    int static constexpr DATA_SPEC_SUPERCELL = 1;
    int static constexpr DATA_SPEC_EXTERNAL  = 2;
    int data_spec;
    if (coupler.option_exists("standalone_input_file")) {
      #ifdef PAM_STANDALONE
        std::string inFile = coupler.get_option<std::string>( "standalone_input_file" );
        YAML::Node config = YAML::LoadFile(inFile);
        std::string dataStr = config["initData"].as<std::string>();
        if        (dataStr == "thermal") {
          data_spec = DATA_SPEC_THERMAL;
        } else if (dataStr == "supercell") {
          data_spec = DATA_SPEC_SUPERCELL;
        } else if (dataStr == "external") {
          data_spec = DATA_SPEC_EXTERNAL;
        } else {
          endrun("ERROR: Invalid data_spec");
        }
      #endif
    } else {
      data_spec = DATA_SPEC_EXTERNAL;
    }

    if (data_spec == DATA_SPEC_EXTERNAL) {

    } else {
      real5d state  ("state"  ,num_state  ,nz+2*hs,ny+2*hs,nx+2*hs,nens);
      real5d tracers("tracers",num_tracers,nz+2*hs,ny+2*hs,nx+2*hs,nens);

      auto zmid = dm.get<real const,2>("vertical_midpoint_height");

      if (data_spec == DATA_SPEC_SUPERCELL) {

        init_supercell( coupler , state , tracers );

      } else if (data_spec == DATA_SPEC_THERMAL) {

        real2d hy_dens_cells    ("hy_dens_cells"    ,nz,nens);
        real2d hy_pressure_cells("hy_pressure_cells",nz,nens);

        // Define quadrature weights and points for 3-point rules
        const int nqpoints = 9;
        SArray<real,1,nqpoints> qpoints;
        SArray<real,1,nqpoints> qweights;

        TransformMatrices::get_gll_points ( qpoints  );
        TransformMatrices::get_gll_weights( qweights );

        // Compute hydrostatic background cell averages using quadrature
        parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
          hy_dens_cells    (k,iens) = 0.;
          hy_pressure_cells(k,iens) = 0.;
          for (int kk=0; kk<nqpoints; kk++) {
            real z = zmid(k,iens) + qpoints(kk)*dz(k,iens);
            real hr, ht;

            if (data_spec == DATA_SPEC_THERMAL) { hydro_const_theta(z,grav,C0,cp_d,p0,gamma,R_d,hr,ht); }

            hy_dens_cells    (k,iens) += hr                       * qweights(kk);
            hy_pressure_cells(k,iens) += C0*std::pow(hr*ht,gamma) * qweights(kk);
          }
        });

        // Use quadrature to initialize state and tracer data
        parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
          for (int l=0; l < num_state  ; l++) { state  (l,hs+k,hs+j,hs+i,iens) = 0.; }
          for (int l=0; l < num_tracers; l++) { tracers(l,hs+k,hs+j,hs+i,iens) = 0.; }
          //Use Gauss-Legendre quadrature
          for (int kk=0; kk<nqpoints; kk++) {
            for (int jj=0; jj<nqpoints; jj++) {
              for (int ii=0; ii<nqpoints; ii++) {
                real x = (i+0.5)*dx   + qpoints(ii)*dx;
                real y = (j+0.5)*dy   + qpoints(jj)*dy;   if (sim2d) y = ylen/2;
                real z = zmid(k,iens) + qpoints(kk)*dz(k,iens);

                real hr = hy_dens_cells    (k,iens);
                real hp = hy_pressure_cells(k,iens);
                real ht = std::pow( hp/C0 , 1._fp/gamma ) / hr;
                real rho = hr;
                real u = 0;
                real v = 0;
                real w = 0;
                real rho_v = 0;
                real theta = ht + sample_ellipse_cosine(2._fp , x,y,z , xlen/2,ylen/2,2000. , 2000.,2000.,2000.);

                if (sim2d) v = 0;

                real wt = qweights(ii)*qweights(jj)*qweights(kk);
                state(idR,hs+k,hs+j,hs+i,iens) += rho       * wt;
                state(idU,hs+k,hs+j,hs+i,iens) += rho*u     * wt;
                state(idV,hs+k,hs+j,hs+i,iens) += rho*v     * wt;
                state(idW,hs+k,hs+j,hs+i,iens) += rho*w     * wt;
                state(idT,hs+k,hs+j,hs+i,iens) += rho*theta * wt;
                for (int tr=0; tr < num_tracers; tr++) {
                  if (tr == idWV) { tracers(tr,hs+k,hs+j,hs+i,iens) += rho_v * wt; }
                  else            { tracers(tr,hs+k,hs+j,hs+i,iens) += 0     * wt; }
                }
              }
            }
          }
        });

      }
      convert_dynamics_to_coupler( coupler , state , tracers );
    }
  }



  // Initialize the supercell test case
  void init_supercell( pam::PamCoupler &coupler , real5d &state , real5d &tracers ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    real constexpr z_0    = 0;
    real constexpr z_trop = 12000;
    real constexpr T_0    = 300;
    real constexpr T_trop = 213;
    real constexpr T_top  = 213;
    real constexpr p_0    = 100000;

    int constexpr ngll = 9;
    yakl::SArray<real,1,ngll> gll_pts, gll_wts;
    TransformMatrices::get_gll_points (gll_pts);
    TransformMatrices::get_gll_weights(gll_wts);

    auto nens        = coupler.get_nens();
    auto nx          = coupler.get_nx();
    auto ny          = coupler.get_ny();
    auto nz          = coupler.get_nz();
    auto dx          = coupler.get_dx();
    auto dy          = coupler.get_dy();
    auto xlen        = coupler.get_xlen();
    auto ylen        = coupler.get_ylen();
    auto sim2d       = ny == 1;
    auto R_d         = coupler.get_option<real>("R_d"    );
    auto R_v         = coupler.get_option<real>("R_v"    );
    auto grav        = coupler.get_option<real>("grav"   );
    auto gamma       = coupler.get_option<real>("gamma_d");
    auto C0          = coupler.get_option<real>("C0"     );
    auto idWV        = coupler.get_option<int >("idWV");
    auto num_tracers = coupler.get_num_tracers();

    real2d hy_dens_cells    ("hy_dens_cells"    ,nz,nens);
    real2d hy_pressure_cells("hy_pressure_cells",nz,nens);

    auto &dm  = coupler.get_data_manager_device_readwrite();
    auto dz   = dm.get<real const,2>("vertical_cell_dz");
    auto zmid = dm.get<real const,2>("vertical_midpoint_height");
    auto zint = dm.get<real const,2>("vertical_interface_height");

    // Temporary arrays used to compute the initial state for high-CAPE supercell conditions
    real4d quad_temp       ("quad_temp"       ,nz,ngll-1,ngll,nens);
    real3d hyDensGLL       ("hyDensGLL"       ,nz,ngll,nens);
    real3d hyDensThetaGLL  ("hyDensThetaGLL"  ,nz,ngll,nens);
    real3d hyDensVapGLL    ("hyDensVapGLL"    ,nz,ngll,nens);
    real3d hyPressureGLL   ("hyPressureGLL"   ,nz,ngll,nens);
    real2d hyDensCells     ("hyDensCells"     ,nz,nens);
    real2d hyDensThetaCells("hyDensThetaCells",nz,nens);

    // Compute quadrature term to integrate to get pressure at GLL points
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ngll-1,ngll,nens) ,
                  YAKL_LAMBDA (int k, int kk, int kkk, int iens) {
      // Middle of this cell
      real cellmid   = zmid(k,iens);
      // Bottom, top, and middle of the space between these two ngll GLL points
      real ngll_b    = cellmid + gll_pts(kk  )*dz(k,iens);
      real ngll_t    = cellmid + gll_pts(kk+1)*dz(k,iens);
      real ngll_m    = 0.5_fp * (ngll_b + ngll_t);
      // Compute grid spacing between these ngll GLL points
      real ngll_dz   = dz(k,iens) * ( gll_pts(kk+1) - gll_pts(kk) );
      // Compute the locate of this GLL point within the ngll GLL points
      real zloc      = ngll_m + ngll_dz * gll_pts(kkk);
      // Compute full density at this location
      real temp      = init_supercell_temperature (zloc, z_0, z_trop, zint(nz,iens), T_0, T_trop, T_top);
      real press_dry = init_supercell_pressure_dry(zloc, z_0, z_trop, zint(nz,iens), T_0, T_trop, T_top, p_0, R_d, grav);
      real qvs       = init_supercell_sat_mix_dry(press_dry, temp);
      real relhum    = init_supercell_relhum(zloc, z_0, z_trop);
      if (relhum * qvs > 0.014_fp) relhum = 0.014_fp / qvs;
      real qv        = std::min( 0.014_fp , qvs*relhum );
      quad_temp(k,kk,kkk,iens) = -(1+qv)*grav/(R_d+qv*R_v)/temp;
    });

    // Compute pressure at GLL points
    parallel_for( YAKL_AUTO_LABEL() , nens , YAKL_LAMBDA (int iens) {
      hyPressureGLL(0,0,iens) = p_0;
      for (int k=0; k < nz; k++) {
        for (int kk=0; kk < ngll-1; kk++) {
          real tot = 0;
          for (int kkk=0; kkk < ngll; kkk++) {
            tot += quad_temp(k,kk,kkk,iens) * gll_wts(kkk);
          }
          tot *= dz(k,iens) * ( gll_pts(kk+1) - gll_pts(kk) );
          hyPressureGLL(k,kk+1,iens) = hyPressureGLL(k,kk,iens) * exp( tot );
          if (kk == ngll-2 && k < nz-1) {
            hyPressureGLL(k+1,0,iens) = hyPressureGLL(k,ngll-1,iens);
          }
        }
      }
    });

    // Compute hydrostatic background state at GLL points
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(nz,ngll,nens) , YAKL_LAMBDA (int k, int kk, int iens) {
      real zloc = zmid(k,iens) + gll_pts(kk)*dz(k,iens);
      real temp       = init_supercell_temperature (zloc, z_0, z_trop, zint(nz,iens), T_0, T_trop, T_top);
      real press_tmp  = init_supercell_pressure_dry(zloc, z_0, z_trop, zint(nz,iens), T_0, T_trop, T_top, p_0, R_d, grav);
      real qvs        = init_supercell_sat_mix_dry(press_tmp, temp);
      real relhum     = init_supercell_relhum(zloc, z_0, z_trop);
      if (relhum * qvs > 0.014_fp) relhum = 0.014_fp / qvs;
      real qv         = std::min( 0.014_fp , qvs*relhum );
      real press      = hyPressureGLL(k,kk,iens);
      real dens_dry   = press / (R_d+qv*R_v) / temp;
      real dens_vap   = qv * dens_dry;
      real dens       = dens_dry + dens_vap;
      real dens_theta = pow( press / C0 , 1._fp / gamma );
      hyDensGLL     (k,kk,iens) = dens;
      hyDensThetaGLL(k,kk,iens) = dens_theta;
      hyDensVapGLL  (k,kk,iens) = dens_vap;
    });

    // Compute hydrostatic background state over cells
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
      real press_tot      = 0;
      real dens_tot       = 0;
      real dens_vap_tot   = 0;
      real dens_theta_tot = 0;
      for (int kk=0; kk < ngll; kk++) {
        press_tot      += hyPressureGLL (k,kk,iens) * gll_wts(kk);
        dens_tot       += hyDensGLL     (k,kk,iens) * gll_wts(kk);
        dens_vap_tot   += hyDensVapGLL  (k,kk,iens) * gll_wts(kk);
        dens_theta_tot += hyDensThetaGLL(k,kk,iens) * gll_wts(kk);
      }
      real press      = press_tot;
      real dens       = dens_tot;
      real dens_vap   = dens_vap_tot;
      real dens_dry   = dens - dens_vap;
      real R          = dens_dry / dens * R_d + dens_vap / dens * R_v;
      real temp       = press / (dens * R);
      real qv         = dens_vap / dens_dry;
      real zloc       = (k+0.5_fp)*dz(k,iens);
      real press_tmp  = init_supercell_pressure_dry(zloc, z_0, z_trop, zint(nz,iens), T_0, T_trop, T_top, p_0, R_d, grav);
      real qvs        = init_supercell_sat_mix_dry(press_tmp, temp);
      real relhum     = qv / qvs;
      real T          = temp - 273;
      real a          = 17.27;
      real b          = 237.7;
      real tdew       = b * ( a*T / (b + T) + log(relhum) ) / ( a - ( a*T / (b+T) + log(relhum) ) );
      // These are used in the rest of the model
      for (int iens=0; iens < nens; iens++) {
        hy_dens_cells    (k,iens) = dens;
        hy_pressure_cells(k,iens) = press;
      }
    });

    // Initialize the state
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      state(idR,hs+k,hs+j,hs+i,iens) = 0;
      state(idU,hs+k,hs+j,hs+i,iens) = 0;
      state(idV,hs+k,hs+j,hs+i,iens) = 0;
      state(idW,hs+k,hs+j,hs+i,iens) = 0;
      state(idT,hs+k,hs+j,hs+i,iens) = 0;
      for (int tr=0; tr < num_tracers; tr++) { tracers(tr,hs+k,hs+j,hs+i,iens) = 0; }
      real pres = hy_pressure_cells(k,iens);
      state(idR,hs+k,hs+j,hs+i,iens) = hy_dens_cells(k,iens);
      state(idT,hs+k,hs+j,hs+i,iens) = std::pow( pres/C0 , 1._fp/gamma );
      for (int kk=0; kk < ngll; kk++) {
        for (int jj=0; jj < ngll; jj++) {
          for (int ii=0; ii < ngll; ii++) {
            real xloc = (i+0.5_fp)*dx + gll_pts(ii)*dx;
            real yloc = (j+0.5_fp)*dy + gll_pts(jj)*dy;
            real zloc = zmid(k,iens)  + gll_pts(kk)*dz(k,iens);

            if (sim2d) yloc = ylen/2;

            real constexpr zs = 5000;
            real constexpr us = 30;
            real constexpr uc = 15;
            real uvel = zloc < zs  ?  us*(zloc/zs)-uc  :  us-uc;
            real vvel       = 0;
            real wvel       = 0;
            real dens_vap   = hyDensVapGLL(k,kk,iens);

            real factor = gll_wts(ii) * gll_wts(jj) * gll_wts(kk);
            state  (idU ,hs+k,hs+j,hs+i,iens) += state(idR,hs+k,hs+j,hs+i,iens) * uvel  * factor;
            state  (idV ,hs+k,hs+j,hs+i,iens) += state(idR,hs+k,hs+j,hs+i,iens) * vvel  * factor;
            state  (idW ,hs+k,hs+j,hs+i,iens) += state(idR,hs+k,hs+j,hs+i,iens) * wvel  * factor;
            tracers(idWV,hs+k,hs+j,hs+i,iens) += dens_vap                               * factor;
          }
        }
      }
    });
  }



  // Convert dynamics state and tracers arrays to the coupler state and write to the coupler's data
  void convert_dynamics_to_coupler( pam::PamCoupler &coupler ,
                                    realConst5d state        ,
                                    realConst5d tracers      ) const {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    auto nens        = coupler.get_nens();
    auto nx          = coupler.get_nx();
    auto ny          = coupler.get_ny();
    auto nz          = coupler.get_nz();
    auto R_d         = coupler.get_option<real>("R_d"    );
    auto R_v         = coupler.get_option<real>("R_v"    );
    auto gamma       = coupler.get_option<real>("gamma_d");
    auto C0          = coupler.get_option<real>("C0"     );
    auto idWV        = coupler.get_option<int >("idWV");
    auto num_tracers = coupler.get_num_tracers();

    auto &dm = coupler.get_data_manager_device_readwrite();
    auto tracer_adds_mass = dm.get<bool const,1>("tracer_adds_mass");

    // Get state from the coupler
    auto dm_rho_d = dm.get<real,4>("density_dry");
    auto dm_uvel  = dm.get<real,4>("uvel"       );
    auto dm_vvel  = dm.get<real,4>("vvel"       );
    auto dm_wvel  = dm.get<real,4>("wvel"       );
    auto dm_temp  = dm.get<real,4>("temp"       );

    // Get tracers from the coupler
    pam::MultiField<real,4> dm_tracers;
    auto tracer_names = coupler.get_tracer_names();
    for (int tr=0; tr < num_tracers; tr++) { dm_tracers.add_field( dm.get<real,4>(tracer_names[tr]) ); }

    // Convert from state and tracers arrays to the coupler's data
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      real rho   = state(idR,hs+k,hs+j,hs+i,iens);
      real u     = state(idU,hs+k,hs+j,hs+i,iens) / rho;
      real v     = state(idV,hs+k,hs+j,hs+i,iens) / rho;
      real w     = state(idW,hs+k,hs+j,hs+i,iens) / rho;
      real theta = state(idT,hs+k,hs+j,hs+i,iens) / rho;
      real press = C0 * pow( rho*theta , gamma );
      real rho_v = tracers(idWV,hs+k,hs+j,hs+i,iens);
      real rho_d = rho;
      for (int tr=0; tr < num_tracers; tr++) { if (tracer_adds_mass(tr)) rho_d -= tracers(tr,hs+k,hs+j,hs+i,iens); }
      real temp = press / ( rho_d * R_d + rho_v * R_v );
      dm_rho_d(k,j,i,iens) = rho_d;
      dm_uvel (k,j,i,iens) = u;
      dm_vvel (k,j,i,iens) = v;
      dm_wvel (k,j,i,iens) = w;
      dm_temp (k,j,i,iens) = temp;
      for (int tr=0; tr < num_tracers; tr++) { dm_tracers(tr,k,j,i,iens) = tracers(tr,hs+k,hs+j,hs+i,iens); }
    });
  }



  // Convert coupler's data to state and tracers arrays
  void convert_coupler_to_dynamics( pam::PamCoupler &coupler ,
                                    real5d &state            ,
                                    real5d &tracers          ) const {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    auto nens        = coupler.get_nens();
    auto nx          = coupler.get_nx();
    auto ny          = coupler.get_ny();
    auto nz          = coupler.get_nz();
    auto R_d         = coupler.get_option<real>("R_d"    );
    auto R_v         = coupler.get_option<real>("R_v"    );
    auto gamma_d     = coupler.get_option<real>("gamma_d");
    auto C0          = coupler.get_option<real>("C0"     );
    auto idWV        = coupler.get_option<int >("idWV");
    auto grav        = coupler.get_option<real>("grav");
    auto num_tracers = coupler.get_num_tracers();

    auto &dm = coupler.get_data_manager_device_readwrite();
    auto tracer_adds_mass = dm.get<bool const,1>("tracer_adds_mass");

    // Get the coupler's state (as const because it's read-only)
    auto dm_rho_d = dm.get<real const,4>("density_dry");
    auto dm_uvel  = dm.get<real const,4>("uvel"       );
    auto dm_vvel  = dm.get<real const,4>("vvel"       );
    auto dm_wvel  = dm.get<real const,4>("wvel"       );
    auto dm_temp  = dm.get<real const,4>("temp"       );

    // Get the coupler's tracers (as const because it's read-only)
    pam::MultiField<real const,4> dm_tracers;
    auto tracer_names = coupler.get_tracer_names();
    for (int tr=0; tr < num_tracers; tr++) { dm_tracers.add_field( dm.get<real const,4>(tracer_names[tr]) ); }

    // Convert from the coupler's state to the dycore's state and tracers arrays.
    // Compute domain-averaged pressure column
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      real rho_d = dm_rho_d(k,j,i,iens);
      real u     = dm_uvel (k,j,i,iens);
      real v     = dm_vvel (k,j,i,iens);
      real w     = dm_wvel (k,j,i,iens);
      real temp  = dm_temp (k,j,i,iens);
      real rho_v = dm_tracers(idWV,k,j,i,iens);
      real press = rho_d * R_d * temp + rho_v * R_v * temp;
      real rho   = rho_d;
      for (int tr=0; tr < num_tracers; tr++) { if (tracer_adds_mass(tr)) rho += dm_tracers(tr,k,j,i,iens); }
      real theta = pow( press/C0 , 1._fp / gamma_d ) / rho;
      state(idR,hs+k,hs+j,hs+i,iens) = rho;
      state(idU,hs+k,hs+j,hs+i,iens) = rho * u;
      state(idV,hs+k,hs+j,hs+i,iens) = rho * v;
      state(idW,hs+k,hs+j,hs+i,iens) = rho * w;
      state(idT,hs+k,hs+j,hs+i,iens) = rho * theta;
      for (int tr=0; tr < num_tracers; tr++) { tracers(tr,hs+k,hs+j,hs+i,iens) = dm_tracers(tr,k,j,i,iens); }
    });
  }



  void declare_current_profile_as_hydrostatic( pam::PamCoupler &coupler ) const {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    auto nens        = coupler.get_nens();
    auto nx          = coupler.get_nx();
    auto ny          = coupler.get_ny();
    auto nz          = coupler.get_nz();
    auto gamma_d     = coupler.get_option<real>("gamma_d");
    auto C0          = coupler.get_option<real>("C0"     );
    auto num_tracers = coupler.get_num_tracers();
    auto &dm                   = coupler.get_data_manager_device_readwrite();
    auto dz                    = dm.get<real const,2>("vertical_cell_dz");
    auto vert_weno_recon_lower = dm.get<real const,5>("vert_weno_recon_lower");
    auto vert_sten_to_coefs    = dm.get<real const,4>("vert_sten_to_coefs");

    auto grav_balance = coupler.get_option<bool>("balance_hydrostasis_with_gravity");

    real5d state  ("state"  ,num_state  ,nz+2*hs,ny+2*hs,nx+2*hs,nens);
    real5d tracers("tracers",num_tracers,nz+2*hs,ny+2*hs,nx+2*hs,nens);
    convert_coupler_to_dynamics( coupler , state , tracers );
    
    if (grav_balance) {

      SArray<real,2,ord,2> coefs_to_gll;
      SArray<real,1,hs+1> idl;
      real                sigma;
      TransformMatrices::coefs_to_gll_lower(coefs_to_gll );
      weno::wenoSetIdealSigma<ord>(idl,sigma);

      auto grav_var = coupler.get_data_manager_device_readwrite().get<real,2>("variable_gravity");
      // Discretize interface pressure spatially exactly as we would discretize it in the solver
      real4d pressure("pressure",nz+2*hs,ny+2*hs,nx+2*hs,nens);
      parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        pressure(hs+k,hs+j,hs+i,iens) = C0*std::pow(state(idT,hs+k,hs+j,hs+i,iens),gamma_d);
        state(idT,hs+k,hs+j,hs+i,iens) /= state(idR,hs+k,hs+j,hs+i,iens);
        if (j == 0 && i == 0) grav_var(k,iens) = 0;
      });
      halo_exchange( coupler , state , tracers , pressure );
      real4d pint("pint",nz+1,ny,nx,nens);
      parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz+1,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        SArray<real,1,ord> stencil;
        SArray<real,2,ord,ord>  s2c_loc[2];
        SArray<real,3,hs,hs,hs> wrl_loc[2];
        for (int i1=0; i1 < ord; i1++) {
          for (int i2=0; i2 < ord; i2++) {
            s2c_loc[0](i1,i2) = vert_sten_to_coefs(k  ,i1,i2,iens);
            s2c_loc[1](i1,i2) = vert_sten_to_coefs(k+1,i1,i2,iens);
          }
        }
        for (int i1=0; i1 < hs; i1++) {
          for (int i2=0; i2 < hs; i2++) {
            for (int i3=0; i3 < hs; i3++) {
              wrl_loc[0](i1,i2,i3) = vert_weno_recon_lower(k  ,i1,i2,i3,iens);
              wrl_loc[1](i1,i2,i3) = vert_weno_recon_lower(k+1,i1,i2,i3,iens);
            }
          }
        }
        int k_upw = 0;
        for (int s=0; s < ord; s++) { stencil(s) = pressure(k+k_upw+s,hs+j,hs+i,iens); }
        real pp_L = reconstruct(stencil,coefs_to_gll,s2c_loc[k_upw],wrl_loc[k_upw],idl,sigma,1-k_upw);
        k_upw = 1;
        for (int s=0; s < ord; s++) { stencil(s) = pressure(k+k_upw+s,hs+j,hs+i,iens); }
        real pp_R = reconstruct(stencil,coefs_to_gll,s2c_loc[k_upw],wrl_loc[k_upw],idl,sigma,1-k_upw);
        pint(k,j,i,iens) = 0.5_fp * (pp_L + pp_R);
      });
      // Compute average column of variable gravity that provides exact hydrostatic balance
      real r_nx_ny = 1./(nx*ny);
      parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        real dens = state(idR,hs+k,hs+j,hs+i,iens);
        yakl::atomicAdd( grav_var(k,iens) , -(pint(k+1,j,i,iens)-pint(k,j,i,iens))/(dens*dz(k,iens))*r_nx_ny );
      });

    } else {

      auto hy_dens_cells     = coupler.get_data_manager_device_readwrite().get<real,2>("hy_dens_cells");
      auto hy_pressure_cells = coupler.get_data_manager_device_readwrite().get<real,2>("hy_pressure_cells");
      hy_dens_cells     = 0;
      hy_pressure_cells = 0;
      real r_nx_ny = 1./(nx*ny);
      parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        real press = C0 * std::pow( state(idT,hs+k,hs+j,hs+i,iens) , gamma_d );
        yakl::atomicAdd( hy_pressure_cells(k,iens) , press                         *r_nx_ny );
        yakl::atomicAdd( hy_dens_cells    (k,iens) , state(idR,hs+k,hs+j,hs+i,iens)*r_nx_ny );
      });

    }
  }



  char const * dycore_name() const { return "SSPRK3+WENO+FV A-grid"; }



  void finalize( pam::PamCoupler const &coupler ) const { }

};

