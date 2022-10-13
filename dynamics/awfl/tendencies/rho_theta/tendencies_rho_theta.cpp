
#include "tendencies/rho_theta/tendencies_rho_theta.h"


namespace awfl {
namespace tendencies_rho_theta {



  real5d createStateArr(pam::PamCoupler const &coupler) {
    auto nz   = coupler.get_nz  ();
    auto ny   = coupler.get_ny  ();
    auto nx   = coupler.get_nx  ();
    auto nens = coupler.get_nens();
    return real5d("stateArr",num_state,nz+2*hs,ny+2*hs,nx+2*hs,nens);
  }



  real5d createTracerArr(pam::PamCoupler const &coupler) {
    auto num_tracers = coupler.get_num_tracers();
    auto nz          = coupler.get_nz  ();
    auto ny          = coupler.get_ny  ();
    auto nx          = coupler.get_nx  ();
    auto nens        = coupler.get_nens();
    return real5d("tracerArr",num_tracers,nz+2*hs,ny+2*hs,nx+2*hs,nens);
  }



  real5d createStateTendArr(pam::PamCoupler const &coupler) {
    auto nz   = coupler.get_nz  ();
    auto ny   = coupler.get_ny  ();
    auto nx   = coupler.get_nx  ();
    auto nens = coupler.get_nens();
    return real5d("stateTendArr",num_state,nz,ny,nx,nens);
  }



  real5d createTracerTendArr(pam::PamCoupler const &coupler) {
    auto num_tracers = coupler.get_num_tracers();
    auto nz          = coupler.get_nz  ();
    auto ny          = coupler.get_ny  ();
    auto nx          = coupler.get_nx  ();
    auto nens        = coupler.get_nens();
    return real5d("tracerTendArr",coupler.get_num_tracers(),nz,ny,nx,nens);
  }



  void convert_dynamics_to_coupler_state( pam::PamCoupler &coupler ,
                                          realConst5d state ,
                                          realConst5d tracers ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;

    auto &dm = coupler.get_data_manager_readwrite();
    real4d dm_dens_dry = dm.get<real,4>( "density_dry"      );
    real4d dm_uvel     = dm.get<real,4>( "uvel"             );
    real4d dm_vvel     = dm.get<real,4>( "vvel"             );
    real4d dm_wvel     = dm.get<real,4>( "wvel"             );
    real4d dm_temp     = dm.get<real,4>( "temp"             );
    auto num_tracers      = coupler.get_num_tracers();
    auto nz               = coupler.get_nz     ();
    auto ny               = coupler.get_ny     ();
    auto nx               = coupler.get_nx     ();
    auto nens             = coupler.get_nens   ();
    auto Rd               = coupler.get_R_d    ();
    auto Rv               = coupler.get_R_v    ();
    auto cp               = coupler.get_cp_d   ();
    auto p0               = coupler.get_p0     ();
    auto gamma            = coupler.get_gamma_d();
    auto kappa            = Rd/cp;
    auto C0               = pow( Rd * pow( p0 , -kappa ) , gamma );
    auto tracer_name      = coupler.get_tracer_names();
    auto tracer_adds_mass = coupler.get_tracer_adds_mass_array();
    auto idWV             = coupler.get_tracer_index("water_vapor");

    pam::MultipleFields<max_tracers,real4d> dm_tracers;
    for (int tr = 0; tr < num_tracers; tr++) {
      auto trac = dm.get<real,4>( tracer_name[tr] );
      dm_tracers.add_field( trac );
    }

    parallel_for( "Spatial.h d2c" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      real dens  = state(idR,hs+k,hs+j,hs+i,iens);
      real uvel  = state(idU,hs+k,hs+j,hs+i,iens) / dens;
      real vvel  = state(idV,hs+k,hs+j,hs+i,iens) / dens;
      real wvel  = state(idW,hs+k,hs+j,hs+i,iens) / dens;
      real theta = state(idT,hs+k,hs+j,hs+i,iens) / dens;
      real pressure = C0 * pow( dens*theta , gamma );
      real dens_vap = tracers(idWV,hs+k,hs+j,hs+i,iens);
      real dens_dry = dens;
      for (int tr=0; tr < num_tracers; tr++) {
        if (tracer_adds_mass(tr)) dens_dry -= tracers(tr,hs+k,hs+j,hs+i,iens);
      }
      real temp = pressure / ( dens_dry * Rd + dens_vap * Rv );
      dm_dens_dry(k,j,i,iens) = dens_dry;
      dm_uvel    (k,j,i,iens) = uvel;
      dm_vvel    (k,j,i,iens) = vvel;
      dm_wvel    (k,j,i,iens) = wvel;
      dm_temp    (k,j,i,iens) = temp;
      for (int tr=0; tr < num_tracers; tr++) {
        dm_tracers(tr,k,j,i,iens) = tracers(tr,hs+k,hs+j,hs+i,iens);
      }
    });
  }



  void convert_coupler_state_to_dynamics( pam::PamCoupler const &coupler ,
                                          real5d const &state ,
                                          real5d const &tracers ,
                                          Hydrostasis &hydrostasis) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;

    auto &dm = coupler.get_data_manager_readonly();
    auto hy_params = dm.get<real const,3>("hydrostasis_parameters");
    auto num_tracers      = coupler.get_num_tracers();
    auto nz               = coupler.get_nz     ();
    auto ny               = coupler.get_ny     ();
    auto nx               = coupler.get_nx     ();
    auto dz               = dm.get<real const,2>("vertical_cell_dz");
    auto nens             = coupler.get_nens   ();
    auto Rd               = coupler.get_R_d    ();
    auto Rv               = coupler.get_R_v    ();
    auto cp               = coupler.get_cp_d   ();
    auto p0               = coupler.get_p0     ();
    auto gamma            = coupler.get_gamma_d();
    auto kappa            = Rd/cp;
    auto C0               = pow( Rd * pow( p0 , -kappa ) , gamma );
    auto grav             = coupler.get_grav   ();
    auto tracer_name      = coupler.get_tracer_names();
    auto tracer_adds_mass = coupler.get_tracer_adds_mass_array();
    auto idWV             = coupler.get_tracer_index("water_vapor");
    auto dm_dens_dry = dm.get<real const,4>( "density_dry"      );
    auto dm_uvel     = dm.get<real const,4>( "uvel"             );
    auto dm_vvel     = dm.get<real const,4>( "vvel"             );
    auto dm_wvel     = dm.get<real const,4>( "wvel"             );
    auto dm_temp     = dm.get<real const,4>( "temp"             );
    auto vert_interface = dm.get<real const,2>("vertical_interface_height");

    SArray<real,1,ngll> gllPts_ngll;
    TransformMatrices::get_gll_points (gllPts_ngll);

    YAKL_SCOPE( hyDensSten           , hydrostasis.hyDensSten           );
    YAKL_SCOPE( hyDensThetaSten      , hydrostasis.hyDensThetaSten      );
    YAKL_SCOPE( hyDensGLL            , hydrostasis.hyDensGLL            );
    YAKL_SCOPE( hyDensThetaGLL       , hydrostasis.hyDensThetaGLL       );

    pam::MultipleFields<max_tracers,realConst4d> dm_tracers;
    for (int tr = 0; tr < num_tracers; tr++) {
      auto trac = dm.get<real const,4>( tracer_name[tr] );
      dm_tracers.add_field( trac );
    }

    // If hydrostasis in the coupler has changed, then we need to re-compute
    // hydrostatically balanced cells and GLL points for the dycore's time step
    real tmp = yakl::intrinsics::sum(hy_params);
    if (tmp != hydrostasis.hydrostasis_parameters_sum) {
      // Get dz for ghost cells
      real2d dz_ghost("dz_ghost",nz+2*hs,nens);
      parallel_for( "Spatial.h init 2" , SimpleBounds<2>(nz+2*hs,nens) , YAKL_LAMBDA (int k, int iens) {
        if      (k >= hs && k < hs+nz) { dz_ghost(k,iens) = dz(k-hs,iens); }
        else if (k < hs              ) { dz_ghost(k,iens) = dz(0   ,iens); }
        else if (k >= hs+nz          ) { dz_ghost(k,iens) = dz(nz-1,iens); }
      });
      // Get vertical interfaces for ghost cells
      real2d vert_interface_ghost("vert_interface_ghost",nz+2*hs+1,nens);
      parallel_for( "Spatial.h init 3" , nens , YAKL_LAMBDA (int iens) {
        vert_interface_ghost(0,iens) = vert_interface(0,iens) - hs*dz(0,iens);
        for (int k=1; k < nz+2*hs+1; k++) {
          vert_interface_ghost(k,iens) = vert_interface_ghost(k-1,iens) + dz_ghost(k-1,iens);
        }
      });

      SArray<real,1,9> gll_pts, gll_wts;
      TransformMatrices::get_gll_points ( gll_pts );
      TransformMatrices::get_gll_weights( gll_wts );

      // Compute new cell averages and GLL point values for hydrostasis
      hydrostasis.hydrostasis_parameters_sum = tmp;
      parallel_for( "Spatial.h new hydrostasis" , SimpleBounds<3>(nz,ord,nens) , YAKL_LAMBDA (int k, int kk, int iens) {
        real r  = 0;
        real rt = 0;
        for (int l=0; l < 9; l++) {
          real zloc = vert_interface_ghost(k+kk,iens) + 0.5_fp*dz_ghost(k+kk,iens) + gll_pts(l)*dz_ghost(k+kk,iens);
          real z0   = vert_interface(k,iens) + 0.5_fp*dz(k,iens);
          real wt = gll_wts(l);
          r  += pam::hydrostatic_density   (hy_params,zloc,z0,dz(k,iens),k,iens         ,grav) * wt;
          rt +=      hydrostatic_dens_theta(hy_params,zloc,z0,dz(k,iens),k,iens,C0,gamma     ) * wt;
        }
        hyDensSten     (k,kk,iens) = r;
        hyDensThetaSten(k,kk,iens) = rt;
      });

      parallel_for( "Spatial.h new hydrostasis" , SimpleBounds<3>(nz,ngll,nens) , YAKL_LAMBDA (int k, int kk, int iens) {
        real zloc = vert_interface(k,iens) + 0.5_fp*dz(k,iens) + gllPts_ngll(kk)*dz(k,iens);
        real z0   = vert_interface(k,iens) + 0.5_fp*dz(k,iens);
        hyDensGLL     (k,kk,iens) = pam::hydrostatic_density   (hy_params,zloc,z0,dz(k,iens),k,iens         ,grav);
        hyDensThetaGLL(k,kk,iens) =      hydrostatic_dens_theta(hy_params,zloc,z0,dz(k,iens),k,iens,C0,gamma     );
      });
    }

    parallel_for(  "Spatial.h c2d" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      for (int tr=0; tr < num_tracers; tr++) {
        tracers(tr,hs+k,hs+j,hs+i,iens) = dm_tracers(tr,k,j,i,iens);
      }
      real dens_dry = dm_dens_dry    (k,j,i,iens);
      real uvel     = dm_uvel        (k,j,i,iens);
      real vvel     = dm_vvel        (k,j,i,iens);
      real wvel     = dm_wvel        (k,j,i,iens);
      real temp     = dm_temp        (k,j,i,iens);
      real dens_vap = tracers(idWV,hs+k,hs+j,hs+i,iens);
      real dens     = dens_dry;
      for (int tr=0; tr < num_tracers; tr++) {
        if (tracer_adds_mass(tr)) dens += tracers(tr,hs+k,hs+j,hs+i,iens);
      }
      real pressure = dens_dry * Rd * temp + dens_vap * Rv * temp;
      real theta    = pow( pressure / C0 , 1._fp / gamma ) / dens;
      state(idR,hs+k,hs+j,hs+i,iens) = dens;
      state(idU,hs+k,hs+j,hs+i,iens) = dens * uvel;
      state(idV,hs+k,hs+j,hs+i,iens) = dens * vvel;
      state(idW,hs+k,hs+j,hs+i,iens) = dens * wvel;
      state(idT,hs+k,hs+j,hs+i,iens) = dens * theta;
    });
  }



  // Initialize the state
  void init_idealized_state_and_tracers( pam::PamCoupler &coupler ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;

    int constexpr DATA_SPEC_THERMAL       = 0;
    int constexpr DATA_SPEC_THERMAL_MOIST = 1;

    int data_spec;

    if (coupler.option_exists("standalone_input_file")) {
      std::string inFile = coupler.get_option<std::string>( "standalone_input_file" );
      YAML::Node config = YAML::LoadFile(inFile);
      std::string dataStr = config["initData"].as<std::string>("unspecified");
      if      (dataStr == "thermal_moist") { data_spec = DATA_SPEC_THERMAL_MOIST; }
      else if (dataStr == "thermal"      ) { data_spec = DATA_SPEC_THERMAL      ; }
      else                                 { return; }
    } else { return; }

    auto num_tracers    = coupler.get_num_tracers();
    auto nz             = coupler.get_nz     ();
    auto ny             = coupler.get_ny     ();
    auto nx             = coupler.get_nx     ();
    auto dx             = coupler.get_dx     ();
    auto dy             = coupler.get_dy     ();
    auto dz             = coupler.get_data_manager_readonly().get<real const,2>("vertical_cell_dz");
    auto nens           = coupler.get_nens   ();
    auto Rd             = coupler.get_R_d    ();
    auto Rv             = coupler.get_R_v    ();
    auto cp             = coupler.get_cp_d   ();
    auto p0             = coupler.get_p0     ();
    auto gamma          = coupler.get_gamma_d();
    auto grav           = coupler.get_grav   ();
    auto kappa          = Rd/cp;
    auto C0             = pow( Rd * pow( p0 , -kappa ) , gamma );
    auto xlen           = coupler.get_xlen();
    auto ylen           = coupler.get_ylen();
    auto idWV           = coupler.get_tracer_index("water_vapor");
    auto vert_interface = coupler.get_data_manager_readonly().get<real const,2>("vertical_interface_height");
    auto sim2d = ny == 1;

    SArray<real,1,ord> gllWts_ord, gllPts_ord;
    TransformMatrices::get_gll_points (gllPts_ord);
    TransformMatrices::get_gll_weights(gllWts_ord);

    real5d state   = createStateArr (coupler);
    real5d tracers = createTracerArr(coupler);

    // Initialize constant theta profile-based test cases
    if ( data_spec == DATA_SPEC_THERMAL_MOIST ||
         data_spec == DATA_SPEC_THERMAL ) {

      // Compute the state
      parallel_for( "Spatial.h init_state 3" , SimpleBounds<4>(nz,ny,nx,nens) ,
                    YAKL_LAMBDA (int k, int j, int i, int iens) {
        state(idR,hs+k,hs+j,hs+i,iens) = 0;
        state(idU,hs+k,hs+j,hs+i,iens) = 0;
        state(idV,hs+k,hs+j,hs+i,iens) = 0;
        state(idW,hs+k,hs+j,hs+i,iens) = 0;
        state(idT,hs+k,hs+j,hs+i,iens) = 0;
        for (int tr=0; tr < num_tracers; tr++) { tracers(tr,hs+k,hs+j,hs+i,iens) = 0; }
        for (int kk=0; kk<ord; kk++) {
          for (int jj=0; jj<ord; jj++) {
            for (int ii=0; ii<ord; ii++) {
              real zloc = vert_interface(k,iens) + 0.5_fp*dz(k,iens) + gllPts_ord(kk)*dz(k,iens);
              real yloc;
              if (sim2d) {
                yloc = ylen/2;
              } else {
                yloc = (j+0.5_fp)*dy + gllPts_ord(jj)*dy;
              }
              real xloc = (i+0.5_fp)*dx + gllPts_ord(ii)*dx;
              real wt = gllWts_ord(kk) * gllWts_ord(jj) * gllWts_ord(ii);
              if        (data_spec == DATA_SPEC_THERMAL_MOIST) {
                // Compute constant theta hydrostatic background state
                real th        = 300;
                real rho_d     = profiles::initConstTheta_density(th,zloc,Rd,cp,gamma,p0,C0,grav);
                real theta_d   = th + profiles::ellipsoid_linear(xloc, yloc, zloc, xlen/2, ylen/2, 2000, 2000, 2000, 2000, 2 );
                real press_d   = C0*pow(rho_d*theta_d,gamma);
                real temp      = press_d / (rho_d * Rd);
                real svp       = profiles::saturation_vapor_pressure(temp); // Self-explanatory
                real press_v   = svp * profiles::ellipsoid_linear(xloc,yloc,zloc  ,  xlen/2,ylen/2,2000  ,  2000,2000,2000  ,  0.8);
                real rho_v     = press_v / (Rv*temp);
                real press     = rho_d * Rd * temp + rho_v * Rv *temp;
                real rho_theta = pow( press/C0 , 1._fp/gamma );
                real rho       = rho_d + rho_v;

                state  (idR ,hs+k,hs+j,hs+i,iens) += rho       * wt;
                state  (idT ,hs+k,hs+j,hs+i,iens) += rho_theta * wt;
                tracers(idWV,hs+k,hs+j,hs+i,iens) += rho_v     * wt;
              } else if (data_spec == DATA_SPEC_THERMAL) {
                // Compute constant theta hydrostatic background state
                real th        = 300;
                real rho       = profiles::initConstTheta_density(th,zloc,Rd,cp,gamma,p0,C0,grav);
                real theta     = th + profiles::ellipsoid_linear(xloc, yloc, zloc, xlen/2, ylen/2, 2000, 2000, 2000, 2000, 2 );

                state  (idR ,hs+k,hs+j,hs+i,iens) += rho       * wt;
                state  (idT ,hs+k,hs+j,hs+i,iens) += rho*theta * wt;
              }
            }
          }
        }
      });

    }

    awfl::tendencies_rho_theta::convert_dynamics_to_coupler_state( coupler , state , tracers );
  }



  realHost1d compute_mass( pam::PamCoupler const &coupler , realConst5d state , realConst5d tracers ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;

    int nz   = coupler.get_nz();
    int ny   = coupler.get_ny();
    int nx   = coupler.get_nx();
    int nens = coupler.get_nens();

    using awfl::tendencies_rho_theta::idR;
    int num_tracers = coupler.get_num_tracers();
    auto dz = coupler.get_data_manager_readonly().get<real const,2>("vertical_cell_dz");

    realHost1d mass("mass",num_tracers+1);
    real4d tmp("tmp",nz,ny,nx,nens);

    parallel_for( "Temporal_ader.h state mass" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      tmp(k,j,i,iens) = state(idR,hs+k,hs+j,hs+i,iens) * dz(k,iens);
    });
    mass(0) = yakl::intrinsics::sum(tmp);

    for (int l=0; l < num_tracers; l++) {
      parallel_for( "Temporal_ader.h tracer mass" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        tmp(k,j,i,iens) = tracers(l,hs+k,hs+j,hs+i,iens) * dz(k,iens);
      });
      mass(l+1) = yakl::intrinsics::sum(tmp);
    }
    return mass;
  }



}
}


