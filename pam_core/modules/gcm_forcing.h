
#pragma once

#include "pam_coupler.h"

namespace modules {

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
  inline void compute_gcm_forcing_tendencies( pam::PamCoupler &coupler ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    using yakl::atomicAdd;

    auto &dm = coupler.get_data_manager_device_readwrite();

    auto dt_gcm = coupler.get_option<real>("gcm_physics_dt");
    real cp_d   = coupler.get_option<real>("cp_d");
    real grav   = coupler.get_option<real>("grav");
    real Lv     = coupler.get_option<real>("latvap") ;
    real Lf     = coupler.get_option<real>("latice") ;

    // Get current state from coupler
    auto rho_d = dm.get<real const,4>( "density_dry" );
    auto uvel  = dm.get<real const,4>( "uvel"        );
    auto vvel  = dm.get<real const,4>( "vvel"        );
    auto temp  = dm.get<real const,4>( "temp"        );
    auto rho_v = dm.get<real const,4>( "water_vapor" );
    auto rho_l = dm.get<real const,4>( "cloud_water" );
    auto rho_i = dm.get<real const,4>( "ice"         );

    auto rho_d_gcm = dm.get<real const,2> ( "gcm_density_dry" );
    auto uvel_gcm  = dm.get<real const,2> ( "gcm_uvel"        );
    auto vvel_gcm  = dm.get<real const,2> ( "gcm_vvel"        );
    auto temp_gcm  = dm.get<real const,2> ( "gcm_temp"        );
    auto rho_v_gcm = dm.get<real const,2> ( "gcm_water_vapor" );
    auto rho_l_gcm = dm.get<real const,2> ( "gcm_cloud_water" );
    auto rho_i_gcm = dm.get<real const,2> ( "gcm_cloud_ice"   );

    int nz   = dm.get_dimension_size("z"   );
    int ny   = dm.get_dimension_size("y"   );
    int nx   = dm.get_dimension_size("x"   );
    int nens = dm.get_dimension_size("nens");

    // Create arrays to hold the current column average of the CRM internal columns
    real2d colavg_rho_d("colavg_rho_d",nz,nens);
    real2d colavg_uvel ("colavg_uvel" ,nz,nens);
    real2d colavg_vvel ("colavg_vvel" ,nz,nens);
    // #ifdef MMF_PAM_FORCE_ALL_WATER_SPECIES
    // real2d colavg_temp ("colavg_temp" ,nz,nens);
    // real2d colavg_rho_v("colavg_rho_v",nz,nens);
    // real2d colavg_rho_l("colavg_rho_l",nz,nens);
    // real2d colavg_rho_i("colavg_rho_i",nz,nens);
    // #endif
    // #ifdef MMF_PAM_FORCE_TOTAL_WATER
    real2d colavg_temp_adj("colavg_temp_adj",nz,nens);
    real2d colavg_rho_totq("colavg_rho_totq",nz,nens);
    // #endif

    // We will be essentially reducing a summation to these variables, so initialize them to zero
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
      colavg_rho_d(k,iens) = 0;
      colavg_uvel (k,iens) = 0;
      colavg_vvel (k,iens) = 0;
      // #ifdef MMF_PAM_FORCE_ALL_WATER_SPECIES
      // colavg_temp (k,iens) = 0;
      // colavg_rho_v(k,iens) = 0;
      // colavg_rho_l(k,iens) = 0;
      // colavg_rho_i(k,iens) = 0;
      // #endif
      // #ifdef MMF_PAM_FORCE_TOTAL_WATER
      colavg_temp_adj(k,iens) = 0;
      colavg_rho_totq(k,iens) = 0;
      // #endif
    });

    real r_nx_ny  = 1._fp / (nx*ny);  // precompute reciprocal to avoid costly divisions
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      // yakl::atomicAdd ensures only one thread performs an update at a time to avoid data races and wrong answers
      atomicAdd( colavg_rho_d(k,iens) , rho_d(k,j,i,iens) * r_nx_ny );
      atomicAdd( colavg_uvel (k,iens) , uvel (k,j,i,iens) * r_nx_ny );
      atomicAdd( colavg_vvel (k,iens) , vvel (k,j,i,iens) * r_nx_ny );
      // #ifdef MMF_PAM_FORCE_ALL_WATER_SPECIES
      // atomicAdd( colavg_temp (k,iens) , temp (k,j,i,iens) * r_nx_ny );
      // atomicAdd( colavg_rho_v(k,iens) , rho_v(k,j,i,iens) * r_nx_ny );
      // atomicAdd( colavg_rho_l(k,iens) , rho_l(k,j,i,iens) * r_nx_ny );
      // atomicAdd( colavg_rho_i(k,iens) , rho_i(k,j,i,iens) * r_nx_ny );
      // #endif
      // #ifdef MMF_PAM_FORCE_TOTAL_WATER
      real rho_total_water = rho_v(k,j,i,iens) + rho_l(k,j,i,iens) + rho_i(k,j,i,iens);
      real ql_tmp          = rho_l(k,j,i,iens) / ( rho_d(k,j,i,iens) + rho_v(k,j,i,iens) );
      real qi_tmp          = rho_i(k,j,i,iens) / ( rho_d(k,j,i,iens) + rho_v(k,j,i,iens) );
      real liq_adj         = ql_tmp* Lv     / cp_d;
      real ice_adj         = qi_tmp*(Lv+Lf) / cp_d;
      real temp_adj        = temp (k,j,i,iens) - liq_adj - ice_adj;
      atomicAdd( colavg_temp_adj(k,iens) , temp_adj        * r_nx_ny );
      atomicAdd( colavg_rho_totq(k,iens) , rho_total_water * r_nx_ny );
      // #endif
    });

    // We need the GCM forcing tendencies later, so store these in the coupler's data manager
    // If they've already been registered, the do not register them again
    if (! dm.entry_exists("gcm_forcing_tend_uvel")) {
      dm.register_and_allocate<real>( "gcm_forcing_tend_rho_d" , "GCM forcing for dry density"        ,{nz,nens},{"z","nens"});
      dm.register_and_allocate<real>( "gcm_forcing_tend_uvel"  , "GCM forcing for u-velocity"         ,{nz,nens},{"z","nens"});
      dm.register_and_allocate<real>( "gcm_forcing_tend_vvel"  , "GCM forcing for v-velocity"         ,{nz,nens},{"z","nens"});
      dm.register_and_allocate<real>( "gcm_forcing_tend_temp"  , "GCM forcing for temperature"        ,{nz,nens},{"z","nens"});
      dm.register_and_allocate<real>( "gcm_forcing_tend_rho_v" , "GCM forcing for water vapor density",{nz,nens},{"z","nens"});
      dm.register_and_allocate<real>( "gcm_forcing_tend_rho_l" , "GCM forcing for cloud water density",{nz,nens},{"z","nens"});
      dm.register_and_allocate<real>( "gcm_forcing_tend_rho_i" , "GCM forcing for cloud ice density",  {nz,nens},{"z","nens"});
    }

    // Retrieve the GCM forcing tendency arrays so we can set them with the time tencencies
    auto gcm_forcing_tend_rho_d = dm.get<real,2>("gcm_forcing_tend_rho_d");
    auto gcm_forcing_tend_uvel  = dm.get<real,2>("gcm_forcing_tend_uvel" );
    auto gcm_forcing_tend_vvel  = dm.get<real,2>("gcm_forcing_tend_vvel" );
    auto gcm_forcing_tend_temp  = dm.get<real,2>("gcm_forcing_tend_temp" );
    auto gcm_forcing_tend_rho_v = dm.get<real,2>("gcm_forcing_tend_rho_v");
    auto gcm_forcing_tend_rho_l = dm.get<real,2>("gcm_forcing_tend_rho_l");
    auto gcm_forcing_tend_rho_i = dm.get<real,2>("gcm_forcing_tend_rho_i");

    real r_dt_gcm = 1._fp / dt_gcm;  // precompute reciprocal to avoid costly divisions
    // The forcing is the difference between the input GCM state and the current
    // colavg'd CRM state divided by the GCM physics time step to be evenly
    // distributed over the course of the CRM steps for the current GCM step
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
      gcm_forcing_tend_rho_d(k,iens) = ( rho_d_gcm(k,iens) - colavg_rho_d(k,iens) ) * r_dt_gcm;
      gcm_forcing_tend_uvel (k,iens) = ( uvel_gcm (k,iens) - colavg_uvel (k,iens) ) * r_dt_gcm;
      gcm_forcing_tend_vvel (k,iens) = ( vvel_gcm (k,iens) - colavg_vvel (k,iens) ) * r_dt_gcm;
      // #ifdef MMF_PAM_FORCE_ALL_WATER_SPECIES
      // gcm_forcing_tend_temp (k,iens) = ( temp_gcm (k,iens) - colavg_temp (k,iens) ) * r_dt_gcm;
      // gcm_forcing_tend_rho_v(k,iens) = ( rho_v_gcm(k,iens) - colavg_rho_v(k,iens) ) * r_dt_gcm;
      // gcm_forcing_tend_rho_l(k,iens) = ( rho_l_gcm(k,iens) - colavg_rho_l(k,iens) ) * r_dt_gcm;
      // gcm_forcing_tend_rho_i(k,iens) = ( rho_i_gcm(k,iens) - colavg_rho_i(k,iens) ) * r_dt_gcm;
      // #endif
      // #ifdef MMF_PAM_FORCE_TOTAL_WATER
      real gcm_rho_totq = rho_v_gcm(k,iens) + rho_l_gcm(k,iens) + rho_i_gcm(k,iens);
      real gcm_ql_tmp   = rho_l_gcm(k,iens) / ( rho_d_gcm(k,iens) + rho_v_gcm(k,iens) );
      real gcm_qi_tmp   = rho_i_gcm(k,iens) / ( rho_d_gcm(k,iens) + rho_v_gcm(k,iens) );
      real gcm_liq_adj  = gcm_ql_tmp* Lv     / cp_d;
      real gcm_ice_adj  = gcm_qi_tmp*(Lv+Lf) / cp_d;
      real gsm_temp_adj = temp_gcm (k,iens) - gcm_liq_adj - gcm_ice_adj;
      gcm_forcing_tend_temp (k,iens) = ( gsm_temp_adj - colavg_temp_adj(k,iens) ) * r_dt_gcm;
      gcm_forcing_tend_rho_v(k,iens) = ( gcm_rho_totq - colavg_rho_totq(k,iens) ) * r_dt_gcm;
      gcm_forcing_tend_rho_l(k,iens) = 0;
      gcm_forcing_tend_rho_i(k,iens) = 0;
      // #endif
    });
  }

  // generalized hole filling routine to simplify apply_gcm_forcing_tendencies()
  inline void fill_holes( pam::PamCoupler &coupler, real2d &rho_x_neg_mass, std::string tracer_name ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    using yakl::atomicAdd;
    using yakl::ScalarLiveOut;
    using yakl::max;
    auto &dm = coupler.get_data_manager_device_readwrite();
    auto dt = coupler.get_option<real>("crm_dt");
    int nz   = dm.get_dimension_size("z"   );
    int ny   = dm.get_dimension_size("y"   );
    int nx   = dm.get_dimension_size("x"   );
    int nens = dm.get_dimension_size("nens");
    auto dz  = dm.get<real,2>("vertical_cell_dz");

    auto rho_x = dm.get<real,4>( tracer_name );

    // We need these arrays for multiplicative hole filling. Holes are only filled 
    // inside vertical columns at first because it leads to a very infrequent 
    // collision rate for atomidAdd() operations, reducing the runtime.
    real2d rho_x_pos_mass("rho_v_pos_mass",nz,nens);

    // initialize rho_x_pos_mass to zero
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
      rho_x_pos_mass(k,iens) = 0;
    });

    // Calculate available positive mass for hole filling at each vertical level 
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      if (rho_x(k,j,i,iens) > 0) atomicAdd( rho_x_pos_mass(k,iens) ,  rho_x(k,j,i,iens)*dz(k,iens) );
    });

    // The negative is too large if the mass added to fill in negative values is greater than the available mass
    ScalarLiveOut<bool> neg_too_large(false);

    // Correct for added mass by taking away mass from positive-valued cells proportional to the mass in that cell
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      // Determine if mass added to negatives is too large to compensate for in this vertical level
      if (i == 0 && j == 0) {
        if (rho_x_neg_mass(k,iens) > rho_x_pos_mass(k,iens)) neg_too_large = true;
      }
      // Subtract mass proportional to this cells' portion of the total mass at the vertical level
      real factor = rho_x(k,j,i,iens)*dz(k,iens) / rho_x_pos_mass(k,iens);
      rho_x(k,j,i,iens) = max( 0._fp , rho_x(k,j,i,iens) - (rho_x_neg_mass(k,iens) * factor)/dz(k,iens) );
    });

    if (neg_too_large.hostRead()) {
      // If negative mass was too large, let's expand the domain of hole filling to the entire CRM
      // This will have a large contention rate for atomicAdd() and will be less efficient than before
      real1d rho_x_neg_mass_glob("rho_x_neg_mass_glob",nens);
      real1d rho_x_pos_mass_glob("rho_x_pos_mass_glob",nens);

      // These are essentially reductions, so initialize to zero
      parallel_for( YAKL_AUTO_LABEL() , nens , YAKL_LAMBDA (int iens) {
        rho_x_neg_mass_glob(iens) = 0;
        rho_x_pos_mass_glob(iens) = 0;
      });

      // Compute the amount of negative mass we need to compensate for 
      // as well as how much positive mass we still have
      parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        if (i == 0 && j == 0) atomicAdd( rho_x_neg_mass_glob(iens) , max(0._fp,rho_x_neg_mass(k,iens)-rho_x_pos_mass(k,iens)) );
        atomicAdd( rho_x_pos_mass_glob(iens) , rho_x(k,j,i,iens)*dz(k,iens) );
      });

      // Remove mass proportionally to the mass in a given cell
      parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        real factor = rho_x(k,j,i,iens)*dz(k,iens) / rho_x_pos_mass_glob(iens);
        rho_x(k,j,i,iens) = max( 0._fp , rho_x(k,j,i,iens) - (rho_x_neg_mass_glob(iens) * factor)/dz(k,iens) );
      });
      
    } // if (neg_too_large.hostRead()) {

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
  inline void apply_gcm_forcing_tendencies( pam::PamCoupler &coupler ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    using yakl::atomicAdd;
    auto &dm = coupler.get_data_manager_device_readwrite();

    auto dt = coupler.get_option<real>("crm_dt");

    int nz   = dm.get_dimension_size("z"   );
    int ny   = dm.get_dimension_size("y"   );
    int nx   = dm.get_dimension_size("x"   );
    int nens = dm.get_dimension_size("nens");
    auto dz  = dm.get<real,2>("vertical_cell_dz");

    // Current CRM state
    auto rho_d = dm.get<real,4>( "density_dry" );
    auto uvel  = dm.get<real,4>( "uvel"        );
    auto vvel  = dm.get<real,4>( "vvel"        );
    auto temp  = dm.get<real,4>( "temp"        );
    auto rho_v = dm.get<real,4>( "water_vapor" );
    auto rho_l = dm.get<real,4>( "cloud_water" );
    auto rho_i = dm.get<real,4>( "ice"         );

    // GCM forcing tendencies for the average CRM column state
    auto gcm_forcing_tend_rho_d = dm.get<real const,2>("gcm_forcing_tend_rho_d");
    auto gcm_forcing_tend_uvel  = dm.get<real const,2>("gcm_forcing_tend_uvel" );
    auto gcm_forcing_tend_vvel  = dm.get<real const,2>("gcm_forcing_tend_vvel" );
    auto gcm_forcing_tend_temp  = dm.get<real const,2>("gcm_forcing_tend_temp" );
    auto gcm_forcing_tend_rho_v = dm.get<real const,2>("gcm_forcing_tend_rho_v");
    auto gcm_forcing_tend_rho_l = dm.get<real const,2>("gcm_forcing_tend_rho_l");
    auto gcm_forcing_tend_rho_i = dm.get<real const,2>("gcm_forcing_tend_rho_i");

    // We need these arrays for multiplicative hole filling
    // Holes are only filled inside vertical columns at first because it leads to a very infrequent collision
    //    rate for atomidAdd() operations, reducing the runtime of the operations.
    real2d rho_v_neg_mass("rho_v_neg_mass",nz,nens);
    real2d rho_v_pos_mass("rho_v_pos_mass",nz,nens);
    real2d rho_l_neg_mass("rho_l_neg_mass",nz,nens);
    real2d rho_l_pos_mass("rho_l_pos_mass",nz,nens);
    real2d rho_i_neg_mass("rho_i_neg_mass",nz,nens);
    real2d rho_i_pos_mass("rho_i_pos_mass",nz,nens);

    // These are essentially reductions, so initialize to zero
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
      rho_v_neg_mass(k,iens) = 0;
      rho_v_pos_mass(k,iens) = 0;
      rho_l_neg_mass(k,iens) = 0;
      rho_l_pos_mass(k,iens) = 0;
      rho_i_neg_mass(k,iens) = 0;
      rho_i_pos_mass(k,iens) = 0;
    });

    // Apply the GCM forcing, and keep track of negative mass we had to fill and available positive mass
    //    to balance the added mass in the corrective step in the next kernel
    // Multiplicative hole filling is done among a single vertical level at a time unless there isn't
    //    enough mass available to fill the negatives. Then it is done globally
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      // Apply forcing
      rho_d(k,j,i,iens) += gcm_forcing_tend_rho_d(k,iens) * dt;
      uvel (k,j,i,iens) += gcm_forcing_tend_uvel (k,iens) * dt;
      vvel (k,j,i,iens) += gcm_forcing_tend_vvel (k,iens) * dt;
      temp (k,j,i,iens) += gcm_forcing_tend_temp (k,iens) * dt;
      rho_v(k,j,i,iens) += gcm_forcing_tend_rho_v(k,iens) * dt;
      // #ifdef MMF_PAM_FORCE_ALL_WATER_SPECIES
      // rho_l(k,j,i,iens) += gcm_forcing_tend_rho_l(k,iens) * dt;
      // rho_i(k,j,i,iens) += gcm_forcing_tend_rho_i(k,iens) * dt;
      // #endif
      // Compute negative and positive mass for rho_v, and set negative masses to zero (essentially adding mass)
      if (rho_v(k,j,i,iens) < 0) {
        atomicAdd( rho_v_neg_mass(k,iens) , -rho_v(k,j,i,iens)*dz(k,iens) );
        rho_v(k,j,i,iens) = 0;
      }
      // #ifdef MMF_PAM_FORCE_ALL_WATER_SPECIES
      // if (rho_l(k,j,i,iens) < 0) {
      //   atomicAdd( rho_l_neg_mass(k,iens) , -rho_l(k,j,i,iens)*dz(k,iens) );
      //   rho_l(k,j,i,iens) = 0;
      // }
      // if (rho_i(k,j,i,iens) < 0) {
      //   atomicAdd( rho_i_neg_mass(k,iens) , -rho_i(k,j,i,iens)*dz(k,iens) );
      //   rho_i(k,j,i,iens) = 0;
      // }
      // #endif
    });

    // Only do the hole filing if there's negative mass
    if (yakl::intrinsics::sum(rho_v_neg_mass) > 0) { fill_holes(coupler, rho_v_neg_mass,"water_vapor"); }
    // #ifdef MMF_PAM_FORCE_ALL_WATER_SPECIES
    // if (yakl::intrinsics::sum(rho_l_neg_mass) > 0) { fill_holes(coupler, rho_l_neg_mass, "cloud_water"); }
    // if (yakl::intrinsics::sum(rho_i_neg_mass) > 0) { fill_holes(coupler, rho_i_neg_mass, "ice");         }
    // #endif

  } // apply_gcm_forcing_tendencies

}
