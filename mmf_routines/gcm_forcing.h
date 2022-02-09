
#pragma once

#include "pam_const.h"
#include "DataManager.h"
#include "pam_coupler.h"


// This routine is only called once at the beginning of an MMF calculation (at the beginning of a GCM time step)
// 
// Let's call the current GCM state at the beginning of the MMF step for this GCM physics time step: state_gcm
// Let's call the current column-averaged CRM state at the same point in time: state_crm
// The GCM forces the CRM in the MMF by computing (state_gcm - state_crm) and evenly adding this differeince to the
//    CRM state over the course the MMF calculations for this GCM time step. 
// Therefore, if the ONLY forcing to the CRM state were this, we'd have:
//    state_crm_new = state_crm + (state_gcm-state_crm)  ===>  state_crm_new = state_gcm
// Therefore, we are forcing the average CRM column state at the end of the MMF step to be the same as the input GCM
//    state PLUS whatever internal CRM dynamics were produced along the way.
inline void compute_gcm_forcing_tendencies( PamCoupler &coupler , real dt_gcm ) {
  using yakl::atomicAdd;
  auto &dm = coupler.dm;

  // Get current state from coupler
  auto rho_d = dm.get<real const,4>( "density_dry" );
  auto uvel  = dm.get<real const,4>( "uvel"        );
  auto vvel  = dm.get<real const,4>( "vvel"        );
  auto temp  = dm.get<real const,4>( "temp"        );
  auto rho_v = dm.get<real const,4>( "water_vapor" );

  auto rho_d_gcm = dm.get<real const,2> ( "gcm_density_dry" );
  auto uvel_gcm  = dm.get<real const,2> ( "gcm_uvel"        );
  auto vvel_gcm  = dm.get<real const,2> ( "gcm_vvel"        );
  auto temp_gcm  = dm.get<real const,2> ( "gcm_temp"        );
  auto rho_v_gcm = dm.get<real const,2> ( "gcm_water_vapor" );

  int nz   = dm.get_dimension_size("z"   );
  int ny   = dm.get_dimension_size("y"   );
  int nx   = dm.get_dimension_size("x"   );
  int nens = dm.get_dimension_size("nens");

  // Create arrays to hold the current column average of the CRM internal columns
  real2d colavg_rho_d("colavg_rho_d",nz,nens);
  real2d colavg_uvel ("colavg_uvel" ,nz,nens);
  real2d colavg_vvel ("colavg_vvel" ,nz,nens);
  real2d colavg_temp ("colavg_temp" ,nz,nens);
  real2d colavg_rho_v("colavg_rho_v",nz,nens);

  // We will be essentially reducing a summation to these variables, so initialize them to zero
  parallel_for( "Initialize colavg to zero" , SimpleBounds<2>(nz,nens) , 
                YAKL_LAMBDA (int k, int iens) {
    colavg_rho_d(k,iens) = 0;
    colavg_uvel (k,iens) = 0;
    colavg_vvel (k,iens) = 0;
    colavg_temp (k,iens) = 0;
    colavg_rho_v(k,iens) = 0;
  });

  real r_nx_ny  = 1._fp / (nx*ny);  // precompute reciprocal to avoid costly divisions
  parallel_for( "Compute average column of current CRM state" , SimpleBounds<4>(nz,ny,nx,nens) , 
                YAKL_DEVICE_LAMBDA (int k, int j, int i, int iens) {
    // yakl::atomicAdd ensures only one thread performs an update at a time to avoid data races and wrong answers
    atomicAdd( colavg_rho_d(k,iens) , rho_d(k,j,i,iens) * r_nx_ny );
    atomicAdd( colavg_uvel (k,iens) , uvel (k,j,i,iens) * r_nx_ny );
    atomicAdd( colavg_vvel (k,iens) , vvel (k,j,i,iens) * r_nx_ny );
    atomicAdd( colavg_temp (k,iens) , temp (k,j,i,iens) * r_nx_ny );
    atomicAdd( colavg_rho_v(k,iens) , rho_v(k,j,i,iens) * r_nx_ny );
  });

  // We need the GCM forcing tendencies later, so store these in the coupler's data manager
  // If they've already been registered, the do not register them again
  if (! dm.entry_exists("gcm_tend_uvel")) {
    dm.register_and_allocate<real>( "gcm_tend_rho_d" , "GCM forcing tendency for dry density"         ,
                                    {nz,nens} , {"z","nens"} );
    dm.register_and_allocate<real>( "gcm_tend_uvel"  , "GCM forcing tendency for u-velocity"          ,
                                    {nz,nens} , {"z","nens"} );
    dm.register_and_allocate<real>( "gcm_tend_vvel"  , "GCM forcing tendency for v-velocity"          ,
                                    {nz,nens} , {"z","nens"} );
    dm.register_and_allocate<real>( "gcm_tend_temp"  , "GCM forcing tendency for temperature"         ,
                                    {nz,nens} , {"z","nens"} );
    dm.register_and_allocate<real>( "gcm_tend_rho_v" , "GCM forcing tendency for water vapor density" ,
                                    {nz,nens} , {"z","nens"} );
  }

  // Retrieve the GCM forcing tendency arrays so we can set them with the time tencencies
  auto gcm_tend_rho_d = dm.get<real,2>("gcm_tend_rho_d");
  auto gcm_tend_uvel  = dm.get<real,2>("gcm_tend_uvel" );
  auto gcm_tend_vvel  = dm.get<real,2>("gcm_tend_vvel" );
  auto gcm_tend_temp  = dm.get<real,2>("gcm_tend_temp" );
  auto gcm_tend_rho_v = dm.get<real,2>("gcm_tend_rho_v");

  real r_dt_gcm = 1._fp / dt_gcm;  // precompute reciprocal to avoid costly divisions
  // The forcing is the difference between the input GCM state and the current colavg'd CRM state divided by the
  //    GCM physics time step to be evenly distributed over the course of the MMF calculations for this step
  parallel_for( "Compute GCM forcing tendencies" , SimpleBounds<2>(nz,nens) , 
                YAKL_LAMBDA (int k, int iens) {
    gcm_tend_rho_d(k,iens) = ( rho_d_gcm(k,iens) - colavg_rho_d(k,iens) ) * r_dt_gcm;
    gcm_tend_uvel (k,iens) = ( uvel_gcm (k,iens) - colavg_uvel (k,iens) ) * r_dt_gcm;
    gcm_tend_vvel (k,iens) = ( vvel_gcm (k,iens) - colavg_vvel (k,iens) ) * r_dt_gcm;
    gcm_tend_temp (k,iens) = ( temp_gcm (k,iens) - colavg_temp (k,iens) ) * r_dt_gcm;
    gcm_tend_rho_v(k,iens) = ( rho_v_gcm(k,iens) - colavg_rho_v(k,iens) ) * r_dt_gcm;
  });
}




// This routine is intended to be called frequently throughout the MMF calculation
// 
// Apply the precomputed GCM forcing tendencies evenly throughout the course of the MMF calculations
// This applies the forcing over the time domain [t_n , t_n + dt]
// It's possible to produce negative values in tracers in this routine. For instance, if you net remove mass
//    from water vapor, then a column with zero or small water vapor at that height can become negative.
// Therefore a multiplicative hole filler is used to fill in negative values. This routine adds mass to fill
//    negative values and then removes that mass from other cells proportional to each cell's existing mass.
// Multiplicative hole fillers do reduce gradients, but we expect the frequency / magnitude of negative values
//    to be low.
inline void apply_gcm_forcing_tendencies( PamCoupler &coupler , real dt ) {
  using yakl::atomicAdd;
  auto &dm = coupler.dm;

  bool force_density = coupler.get_option<std::string>("density_forcing") == "loose";

  int nz   = dm.get_dimension_size("z"   );
  int ny   = dm.get_dimension_size("y"   );
  int nx   = dm.get_dimension_size("x"   );
  int nens = dm.get_dimension_size("nens");

  // Current CRM state
  auto rho_d = dm.get<real,4>( "density_dry" );
  auto uvel  = dm.get<real,4>( "uvel"        );
  auto vvel  = dm.get<real,4>( "vvel"        );
  auto temp  = dm.get<real,4>( "temp"        );
  auto rho_v = dm.get<real,4>( "water_vapor" );

  // GCM forcing tendencies for the average CRM column state
  auto gcm_tend_rho_d = dm.get<real const,2>("gcm_tend_rho_d");
  auto gcm_tend_uvel  = dm.get<real const,2>("gcm_tend_uvel" );
  auto gcm_tend_vvel  = dm.get<real const,2>("gcm_tend_vvel" );
  auto gcm_tend_temp  = dm.get<real const,2>("gcm_tend_temp" );
  auto gcm_tend_rho_v = dm.get<real const,2>("gcm_tend_rho_v");

  // We need these arrays for multiplicative hole filling
  // Holes are only filled inside vertical columns at first because it leads to a very infrequent collision
  //    rate for atomidAdd() operations, reducing the runtime of the operations.
  real3d rho_v_neg_mass("rho_v_neg_mass",ny,nx,nens);
  real3d rho_v_pos_mass("rho_v_pos_mass",ny,nx,nens);

  // These are essentially reductions, so initialize to zero
  parallel_for( Bounds<3>(ny,nx,nens) , YAKL_LAMBDA (int j, int i, int iens) {
    rho_v_neg_mass(j,i,iens) = 0;
    rho_v_pos_mass(j,i,iens) = 0;
  });

  // Apply the GCM forcing, and keep track of negative mass we had to fill and available positive mass
  //    to balance the added mass in the corrective step in the next kernel
  parallel_for( "Apply GCM forcing" , SimpleBounds<4>(nz,ny,nx,nens) , 
                YAKL_DEVICE_LAMBDA (int k, int j, int i, int iens) {
    // Apply forcing
    if (force_density) rho_d(k,j,i,iens) += gcm_tend_rho_d(k,iens) * dt;
    uvel (k,j,i,iens) += gcm_tend_uvel (k,iens) * dt;
    vvel (k,j,i,iens) += gcm_tend_vvel (k,iens) * dt;
    temp (k,j,i,iens) += gcm_tend_temp (k,iens) * dt;
    rho_v(k,j,i,iens) += gcm_tend_rho_v(k,iens) * dt;
    // Compute negative and positive mass for rho_v, and set negative masses to zero (essentially adding mass)
    if (rho_v(k,j,i,iens) < 0) {
      atomicAdd( rho_v_neg_mass(j,i,iens) , -rho_v(k,j,i,iens) );
      rho_v(k,j,i,iens) = 0;
    } else if (rho_v(k,j,i,iens) > 0) {
      atomicAdd( rho_v_pos_mass(j,i,iens) , rho_v(k,j,i,iens) );
    }
  });

  // For determining if the mass added in a column is larger than the available mass for balancing
  // If it is too large, we'll have to increase the domain of correction, but that will be very infrequent
  ScalarLiveOut<bool> neg_too_large(false);

  // Correct for added mass by taking away mass from positive-valued cells proportional to the mass in that cell
  parallel_for( "Multiplicative hole filler for negative values" , SimpleBounds<4>(nz,ny,nx,nens) , 
                YAKL_LAMBDA (int k, int j, int i, int iens) {
    // Determine if mass added to negatives is too large
    if (rho_v_neg_mass(j,i,iens) > rho_v_pos_mass(j,i,iens)) {
      neg_too_large = true;
    }
    // This is the proportion of this cell's mass to the total positive mass we can subtract mass from
    real factor = rho_v(k,j,i,iens) / rho_v_pos_mass(j,i,iens);
    // Subtract mass relative to this cell's proportion of the available mass
    // The max(0,...) is only to account for machine precision negative values that may develop due to FP arithmetic
    //    except for the case of neg_too_large == true, when this will not conserve mass to machine precision
    rho_v(k,j,i,iens) = max( 0._fp , rho_v(k,j,i,iens) - rho_v_neg_mass(j,i,iens) * factor );
    // If negative mass was too large, then keep track of how much more mass we need to remove
    if (rho_v_neg_mass(j,i,iens) > rho_v_pos_mass(j,i,iens)) {
      rho_v_neg_mass(j,i,iens) = rho_v_neg_mass(j,i,iens) - rho_v_pos_mass(j,i,iens);
    }
  });

  if (neg_too_large.hostRead()) {
    // If negative mass was too large, let's expand the domain of hole filling to the entire CRM
    // This will have a large contention rate for atomicAdd() and will be less efficient than before

    real1d rho_v_neg_mass_glob("rho_v_neg_mass_glob",nens);
    real1d rho_v_pos_mass_glob("rho_v_pos_mass_glob",nens);

    // These are essentially reductions, so initialize to zero
    parallel_for( nens , YAKL_LAMBDA (int iens) {
      rho_v_neg_mass_glob(iens) = 0;
      rho_v_pos_mass_glob(iens) = 0;
    });

    // rho_v_neg_mass(j,i,iens) holds the mass that still needs to be removed. Reducing to rho_v_neg_mass_glob
    //    to determine the total mass in the CRM domain that needs removing
    // Also reducing rho_v (which is guaranteed >= 0 at this point) into rho_v_pos_mass_glob to determine the mass
    //    we have available to work with.
    parallel_for( "Compute global positive and negative masses" , SimpleBounds<4>(nz,ny,nx,nens) , 
                  YAKL_DEVICE_LAMBDA (int k, int j, int i, int iens) {
      if (k == 0) atomicAdd( rho_v_neg_mass(iens) , rho_v_neg_mass(j,i,iens) );
      atomicAdd( rho_v_pos_mass_glob(iens) , rho_v(k,j,i,iens) );
    });

    // Remove mass proportionally to the mass in a given cell
    parallel_for( "Multiplicative hole filler for negative values" , SimpleBounds<4>(nz,ny,nx,nens) , 
                  YAKL_LAMBDA (int k, int j, int i, int iens) {
      // This is the proportion of this cell's mass to the total positive mass we can subtract mass from
      real factor = rho_v(k,j,i,iens) / rho_v_pos_mass_glob(iens);
      // Subtract mass relative to this cell's proportion of the available mass
      // The max(0,...) is only to account for machine precision negative values that may develop due to FP arithmetic
      rho_v(k,j,i,iens) = max( 0._fp , rho_v(k,j,i,iens) - rho_v_neg_mass_glob(iens) * factor );
    });
    
  }
}


