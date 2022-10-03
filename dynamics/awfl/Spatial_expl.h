
#pragma once

#include "awfl_const.h"
#include "TransformMatrices.h"
#include "TransformMatrices_variable.h"
#include "WenoLimiter.h"
#include "idealized_profiles.h"
#include "MultipleFields.h"
#include "pam_coupler.h"
#include "compute_time_step.h"
#include "reconstruct.h"
#include "riemann_rho_theta_full.h"
#include "diff_trans_tracer.h"
#include "diff_trans_rho_theta_full.h"
#include "compute_time_average.h"
#include "fct_positivity.h"


class Spatial_operator {
public:
  int static constexpr nAder       = ngll;
  int static constexpr num_state   = 5;
  int static constexpr max_tracers = 50;
  // For indexing into the state and state tendency arrays
  int static constexpr idR = 0;  // density perturbation
  int static constexpr idU = 1;  // u
  int static constexpr idV = 2;  // v
  int static constexpr idW = 3;  // w
  int static constexpr idT = 4;  // potential temperature perturbation
  // The two boundary condition options for each direction
  int static constexpr BC_PERIODIC = 0;
  int static constexpr BC_WALL     = 1;
  // Options for initializing the data
  int static constexpr DATA_SPEC_EXTERNAL      = 0;
  int static constexpr DATA_SPEC_THERMAL       = 1;
  int static constexpr DATA_SPEC_SUPERCELL     = 2;
  
  bool                     weno_scalars;     // Use WENO limiting for scalars?
  bool                     weno_winds;       // Use WENO limiting for winds?
  int                      data_spec;        // How to initialize the data
  real                     dtInit;           // Initial time step (used throughout the simulation)
  // Hydrostatic background data
  real hydrostasis_parameters_sum;  // Sum of the current parameters (to see if it's changed)
  real3d hyDensSten;                // A stencil around each cell of hydrostatic density
  real3d hyDensThetaSten;           // A stencil around each cell of hydrostatic density * potential temperature
  real3d hyDensGLL;                 // GLL point values of hydrostatic background density in each cell
  real3d hyPressureGLL;             // GLL point values of hydrostatic background pressure in each cell
  real3d hyDensThetaGLL;            // GLL point values of hydrostatic background density*potential temperature
  // Transformation matrices for various degrees of freedom
  SArray<real,2,ord,ngll>       coefs_to_gll;
  SArray<real,2,ord,ngll>       coefs_to_deriv_gll;
  SArray<real,2,ord,ngll>       sten_to_gll;
  SArray<real,2,ord,ord >       sten_to_coefs;
  SArray<real,3,hs+1,hs+1,hs+1> weno_recon_lower;   // WENO reconstruction matrices
  SArray<real,1,hs+2>           idl;                // Ideal weights for WENO
  real                          sigma;              // WENO sigma parameter (handicap high-order TV estimate)
  SArray<real,2,ngll,ngll>      derivMatrix;        // Transform matrix: ngll GLL pts -> ngll GLL derivs
  SArray<real,1,ord >           gllWts_ord;
  SArray<real,1,ord >           gllPts_ord;
  SArray<real,1,ngll>           gllWts_ngll;
  SArray<real,1,ngll>           gllPts_ngll;
  // Vertical grid and reconstruction matrix information
  real2d vert_interface;
  real2d vert_interface_ghost;
  real3d vert_locs_normalized;
  real2d dz_ghost;
  real4d vert_sten_to_gll;
  real4d vert_sten_to_coefs;
  real5d vert_weno_recon_lower;


  // When this class is created, initialize num_tracers to zero
  Spatial_operator() { }



  // Make sure it's odd-order-accurate
  static_assert(ord%2 == 1,"ERROR: ord must be an odd integer");

  void convert_dynamics_to_coupler_state( pam::PamCoupler &coupler , realConst5d state , realConst5d tracers ) const {
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



  YAKL_INLINE static real hydrostatic_dens_theta( realConst3d hy_params , real z , real z0 , real dz ,
                                                  int k, int iens , real C0 , real gamma ) {
    real p = pam::hydrostatic_pressure( hy_params , z , z0 , dz , k , iens );
    // p = C0*(rho*theta)^gamma
    return pow(p/C0,1._fp/gamma);
  }



  void convert_coupler_state_to_dynamics( pam::PamCoupler const &coupler , real5d const &state , real5d const &tracers ) {
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

    YAKL_SCOPE( hyDensSten           , this->hyDensSten           );
    YAKL_SCOPE( hyDensThetaSten      , this->hyDensThetaSten      );
    YAKL_SCOPE( hyPressureGLL        , this->hyPressureGLL        );
    YAKL_SCOPE( hyDensGLL            , this->hyDensGLL            );
    YAKL_SCOPE( hyDensThetaGLL       , this->hyDensThetaGLL       );
    YAKL_SCOPE( gllPts_ngll          , this->gllPts_ngll          );
    YAKL_SCOPE( vert_interface       , this->vert_interface       );
    YAKL_SCOPE( vert_interface_ghost , this->vert_interface_ghost );
    YAKL_SCOPE( dz_ghost             , this->dz_ghost             );

    pam::MultipleFields<max_tracers,realConst4d> dm_tracers;
    for (int tr = 0; tr < num_tracers; tr++) {
      auto trac = dm.get<real const,4>( tracer_name[tr] );
      dm_tracers.add_field( trac );
    }

    // If hydrostasis in the coupler has changed, then we need to re-compute
    // hydrostatically balanced cells and GLL points for the dycore's time step
    real tmp = yakl::intrinsics::sum(hy_params);
    if (tmp != hydrostasis_parameters_sum) {
      SArray<real,1,9> gll_pts, gll_wts;
      TransformMatrices::get_gll_points ( gll_pts );
      TransformMatrices::get_gll_weights( gll_wts );

      // Compute new cell averages and GLL point values for hydrostasis
      hydrostasis_parameters_sum = tmp;
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
        hyPressureGLL (k,kk,iens) = pam::hydrostatic_pressure  (hy_params,zloc,z0,dz(k,iens),k,iens              );
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



  real5d createStateArr(pam::PamCoupler const &coupler) const {
    auto nz   = coupler.get_nz  ();
    auto ny   = coupler.get_ny  ();
    auto nx   = coupler.get_nx  ();
    auto nens = coupler.get_nens();
    return real5d("stateArr",num_state,nz+2*hs,ny+2*hs,nx+2*hs,nens);
  }



  real5d createTracerArr(pam::PamCoupler const &coupler) const {
    auto num_tracers = coupler.get_num_tracers();
    auto nz          = coupler.get_nz  ();
    auto ny          = coupler.get_ny  ();
    auto nx          = coupler.get_nx  ();
    auto nens        = coupler.get_nens();
    return real5d("tracerArr",num_tracers,nz+2*hs,ny+2*hs,nx+2*hs,nens);
  }



  real5d createStateTendArr(pam::PamCoupler const &coupler) const {
    auto nz   = coupler.get_nz  ();
    auto ny   = coupler.get_ny  ();
    auto nx   = coupler.get_nx  ();
    auto nens = coupler.get_nens();
    return real5d("stateTendArr",num_state,nz,ny,nx,nens);
  }



  real5d createTracerTendArr(pam::PamCoupler const &coupler) const {
    auto num_tracers = coupler.get_num_tracers();
    auto nz          = coupler.get_nz  ();
    auto ny          = coupler.get_ny  ();
    auto nx          = coupler.get_nx  ();
    auto nens        = coupler.get_nens();
    return real5d("tracerTendArr",coupler.get_num_tracers(),nz,ny,nx,nens);
  }



  real compute_time_step(pam::PamCoupler const &coupler, real cfl = 0.75) {
    if (dtInit == 0) { dtInit = awfl::compute_time_step(coupler,cfl); }
    return dtInit;
  }



  // Initialize crap needed by recon()
  void init(pam::PamCoupler &coupler) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;

    using yakl::intrinsics::matmul_cr;

    auto nens = coupler.get_nens();
    auto nx   = coupler.get_nx();
    auto ny   = coupler.get_ny();
    auto xlen = coupler.get_xlen();
    auto ylen = coupler.get_ylen();
    auto num_tracers = coupler.get_num_tracers();

    this->hydrostasis_parameters_sum = 0;

    auto Rd    = coupler.get_R_d ();
    auto cp    = coupler.get_cp_d();
    auto p0    = coupler.get_p0  ();
    auto Rv    = coupler.get_R_v ();
    auto grav  = coupler.get_grav();
    auto gamma = cp / (cp-Rd);
    auto kappa = Rd/cp;
    auto C0 = pow( Rd * pow( p0 , -kappa ) , gamma );


    if (! coupler.tracer_exists("water_vapor")) endrun("ERROR: processed registered tracers, and water_vapor was not found");

    fence();

    // Inialize time step to zero, and dimensional splitting switch
    dtInit = 0;

    // If inFile is empty, then we aren't reading in an input file
    if (coupler.option_exists("standalone_input_file")) {
      std::string inFile = coupler.get_option<std::string>( "standalone_input_file" );

      // Read the YAML input file
      YAML::Node config = YAML::LoadFile(inFile);

      // Read whether we're doing WENO limiting on scalars and winds
      weno_scalars = config["weno_scalars"].as<bool>();
      weno_winds   = config["weno_winds"].as<bool>();

      // Read the data initialization option
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
    } else {
      weno_scalars            = true;
      weno_winds              = true;
      data_spec               = DATA_SPEC_EXTERNAL;
    }

    // Store vertical cell interface heights in the data manager
    auto &dm = coupler.get_data_manager_readonly();
    auto zint = dm.get<real const,2>("vertical_interface_height");

    auto nz = coupler.get_nz();

    vert_interface        = real2d("vert_interface"      ,nz+1          ,nens);
    vert_interface_ghost  = real2d("vert_interface_ghost",nz+2*hs+1     ,nens);
    vert_locs_normalized  = real3d("vert_locs_normalized",nz,ord+1      ,nens);
    auto dz               = real2d("dz"                  ,nz            ,nens);
    dz_ghost              = real2d("dz_ghost"            ,nz+2*hs       ,nens);
    vert_sten_to_gll      = real4d("vert_sten_to_gll"     ,nz,ord,ngll,nens);
    vert_sten_to_coefs    = real4d("vert_sten_to_coefs"   ,nz,ord,ord ,nens);
    vert_weno_recon_lower = real5d("vert_weno_recon_lower",nz,hs+1,hs+1,hs+1,nens);

    YAKL_SCOPE( vert_interface        , this->vert_interface        );
    YAKL_SCOPE( vert_interface_ghost  , this->vert_interface_ghost  );
    YAKL_SCOPE( vert_locs_normalized  , this->vert_locs_normalized  );
    YAKL_SCOPE( dz_ghost              , this->dz_ghost              );
    YAKL_SCOPE( vert_sten_to_gll      , this->vert_sten_to_gll      );
    YAKL_SCOPE( vert_sten_to_coefs    , this->vert_sten_to_coefs    );
    YAKL_SCOPE( vert_weno_recon_lower , this->vert_weno_recon_lower );

    zint.deep_copy_to(vert_interface);

    parallel_for( "Spatial.h init 1" , SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
      dz(k,iens) = vert_interface(k+1,iens) - vert_interface(k,iens);
    });

    parallel_for( "Spatial.h init 2" , SimpleBounds<2>(nz+2*hs,nens) , YAKL_LAMBDA (int k, int iens) {
      if (k >= hs && k < hs+nz) {
        dz_ghost(k,iens) = dz(k-hs,iens);
      } else if (k < hs) {
        dz_ghost(k,iens) = dz(0,iens);
      } else if (k >= hs+nz) {
        dz_ghost(k,iens) = dz(nz-1,iens);
      }
    });

    parallel_for( "Spatial.h init 3" , nens , YAKL_LAMBDA (int iens) {
      vert_interface_ghost(0,iens) = vert_interface(0,iens) - hs*dz(0,iens);
      for (int k=1; k < nz+2*hs+1; k++) {
        vert_interface_ghost(k,iens) = vert_interface_ghost(k-1,iens) + dz_ghost(k-1,iens);
      }
    });

    auto vint_host      = vert_interface_ghost .createHostCopy();
    auto vert_s2g_host  = vert_sten_to_gll     .createHostCopy();
    auto vert_s2c_host  = vert_sten_to_coefs   .createHostCopy();
    auto vert_weno_host = vert_weno_recon_lower.createHostCopy();
    auto vert_locs_host = vert_locs_normalized .createHostCopy();

    SArray<real,2,ord,ngll> c2g;
    TransformMatrices::coefs_to_gll_lower(c2g);

    for (int k=0; k < nz; k++) {
      for (int iens = 0; iens < nens; iens++) {
        // Store stencil locations
        SArray<double,1,ord+1> locs;
        for (int kk=0; kk < ord+1; kk++) {
          locs(kk) = vint_host(k+kk,iens);
        }

        // Normalize stencil locations
        double zmid = ( locs(hs+1) + locs(hs) ) / 2;
        double dzmid = locs(hs+1) - locs(hs);
        for (int kk=0; kk < ord+1; kk++) {
          locs(kk) = ( locs(kk) - zmid ) / dzmid;
          vert_locs_host(k,kk,iens) = locs(kk);
        }

        // Compute reconstruction matrices
        SArray<double,2,ord,ord> s2c_var_in;
        SArray<double,3,hs+1,hs+1,hs+1> weno_recon_lower_var;
        TransformMatrices_variable::sten_to_coefs_variable<ord>(locs,s2c_var_in);
        TransformMatrices_variable::weno_lower_sten_to_coefs<ord>(locs,weno_recon_lower_var);
        SArray<real,2,ord,ord> s2c_var;
        for (int jj=0; jj < ord; jj++) {
          for (int ii=0; ii < ord; ii++) {
            s2c_var(jj,ii) = s2c_var_in(jj,ii);
          }
        }
        auto s2g_var = matmul_cr( c2g , s2c_var );

        // Store reconstruction matrices
        for (int jj=0; jj < ord; jj++) {
          for (int ii=0; ii < ord; ii++) {
            vert_s2c_host(k,jj,ii,iens) = s2c_var(jj,ii);
          }
        }

        for (int jj=0; jj < ord; jj++) {
          for (int ii=0; ii < ngll; ii++) {
            vert_s2g_host(k,jj,ii,iens) = s2g_var(jj,ii);
          }
        }

        for (int kk=0; kk < hs+1; kk++) {
          for (int jj=0; jj < hs+1; jj++) {
            for (int ii=0; ii < hs+1; ii++) {
              vert_weno_host(k,kk,jj,ii,iens) = weno_recon_lower_var(kk,jj,ii);
            }
          }
        }

      }
    }

    vert_s2g_host .deep_copy_to(vert_sten_to_gll     );
    vert_s2c_host .deep_copy_to(vert_sten_to_coefs   );
    vert_weno_host.deep_copy_to(vert_weno_recon_lower);
    vert_locs_host.deep_copy_to(vert_locs_normalized );

    // Compute the grid spacing in each dimension
    auto dx = xlen/nx;
    auto dy = ylen/ny;

    // Store the WENO reconstruction matrices
    TransformMatrices::weno_lower_sten_to_coefs(this->weno_recon_lower);

    // Block exists to avoid name mangling stufff
    {
      SArray<real,2,ord,ord>  s2c;        // Converts ord stencil cell averages to ord coefficients
      SArray<real,2,ord,ngll> c2g_lower;  // Converts ord coefficients to ngll GLL points
      SArray<real,2,ord,ord>  c2d;        // Converts ord coefficients to order differentiated coefficients

      TransformMatrices::sten_to_coefs     (s2c      );
      TransformMatrices::coefs_to_gll_lower(c2g_lower);
      TransformMatrices::coefs_to_deriv    (c2d      );

      this->coefs_to_gll       = c2g_lower;
      this->coefs_to_deriv_gll = matmul_cr( c2g_lower , c2d );
      this->sten_to_coefs      = s2c;
      this->sten_to_gll        = matmul_cr( c2g_lower , s2c );
    }
    // Store ader derivMatrix
    {
      SArray<real,2,ngll,ngll> g2c;  // Converts ngll GLL points to ngll coefficients
      SArray<real,2,ngll,ngll> c2d;  // Converts ngll coefficients to ngll differentiated coefficients
      SArray<real,2,ngll,ngll> c2g;  // Converts ngll coefficients to ngll GLL points

      TransformMatrices::gll_to_coefs  (g2c);
      TransformMatrices::coefs_to_deriv(c2d);
      TransformMatrices::coefs_to_gll  (c2g);

      this->derivMatrix = matmul_cr( c2g , matmul_cr( c2d , g2c ) );
    }
    // Store quadrature weights using ord GLL points
    TransformMatrices::get_gll_points (this->gllPts_ord);
    TransformMatrices::get_gll_weights(this->gllWts_ord);
    // Store quadrature weights using ngll GLL points
    TransformMatrices::get_gll_points (this->gllPts_ngll);
    TransformMatrices::get_gll_weights(this->gllWts_ngll);

    // Store WENO ideal weights and sigma value
    weno::wenoSetIdealSigma<ord>(this->idl,this->sigma);

    // Allocate data
    hyDensSten      = real3d("hyDensSten       ",nz,ord,nens);
    hyDensThetaSten = real3d("hyDensThetaSten  ",nz,ord,nens);
    hyDensGLL       = real3d("hyDensGLL        ",nz,ngll,nens);
    hyPressureGLL   = real3d("hyPressureGLL    ",nz,ngll,nens);
    hyDensThetaGLL  = real3d("hyDensThetaGLL   ",nz,ngll,nens);

    init_idealized_state_and_tracers( coupler );

  }



  // Initialize the state
  void init_idealized_state_and_tracers( pam::PamCoupler &coupler ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;


    auto num_tracers = coupler.get_num_tracers();
    auto nz          = coupler.get_nz     ();
    auto ny          = coupler.get_ny     ();
    auto nx          = coupler.get_nx     ();
    auto dx          = coupler.get_dx     ();
    auto dy          = coupler.get_dy     ();
    auto dz          = coupler.get_data_manager_readonly().get<real const,2>("vertical_cell_dz");
    auto nens        = coupler.get_nens   ();
    auto Rd          = coupler.get_R_d    ();
    auto Rv          = coupler.get_R_v    ();
    auto cp          = coupler.get_cp_d   ();
    auto p0          = coupler.get_p0     ();
    auto gamma       = coupler.get_gamma_d();
    auto grav        = coupler.get_grav   ();
    auto kappa       = Rd/cp;
    auto C0          = pow( Rd * pow( p0 , -kappa ) , gamma );
    auto xlen        = coupler.get_xlen();
    auto ylen        = coupler.get_ylen();
    auto idWV        = coupler.get_tracer_index("water_vapor");
    auto sim2d = ny == 1;

    YAKL_SCOPE( gllPts_ord               , this->gllPts_ord              );
    YAKL_SCOPE( gllWts_ord               , this->gllWts_ord              );
    YAKL_SCOPE( data_spec                , this->data_spec               );
    YAKL_SCOPE( vert_interface           , this->vert_interface          );

    // If data's being specified by the driver externally, then there's nothing to do here
    if (data_spec == DATA_SPEC_EXTERNAL) return;

    real5d state   = createStateArr (coupler);
    real5d tracers = createTracerArr(coupler);

    // If the data_spec is thermal or ..., then initialize the domain with Exner pressure-based hydrostasis
    // This is mostly to make plotting potential temperature perturbation easier for publications
    if (data_spec == DATA_SPEC_THERMAL) {

      // Compute the state
      parallel_for( "Spatial.h init_state 3" , SimpleBounds<4>(nz,ny,nx,nens) ,
                    YAKL_LAMBDA (int k, int j, int i, int iens) {
        state(idR,hs+k,hs+j,hs+i,iens) = 0;
        state(idU,hs+k,hs+j,hs+i,iens) = 0;
        state(idV,hs+k,hs+j,hs+i,iens) = 0;
        state(idW,hs+k,hs+j,hs+i,iens) = 0;
        state(idT,hs+k,hs+j,hs+i,iens) = 0;
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
              if        (data_spec == DATA_SPEC_THERMAL) {
                // Compute constant theta hydrostatic background state
                real th = 300;
                real rh = profiles::initConstTheta_density(th,zloc,Rd,cp,gamma,p0,C0,grav);
                real tp = profiles::ellipsoid_linear(xloc, yloc, zloc, xlen/2, ylen/2, 2000, 2000, 2000, 2000, 2 );
                real t = th + tp;
                real r = rh;

                state(idR,hs+k,hs+j,hs+i,iens) += r   * wt;
                state(idT,hs+k,hs+j,hs+i,iens) += r*t * wt;
              }
            }
          }
        }
      });

      parallel_for( "Spatial.h init_tracers" , SimpleBounds<4>(nz,ny,nx,nens) ,
                    YAKL_LAMBDA (int k, int j, int i, int iens) {
        for (int tr=0; tr < num_tracers; tr++) { tracers(tr,hs+k,hs+j,hs+i,iens) = 0; }
        // Loop over quadrature points
        for (int kk=0; kk<ord; kk++) {
          for (int jj=0; jj<ord; jj++) {
            for (int ii=0; ii<ord; ii++) {
              // Get the location
              real zloc = vert_interface(k,iens) + 0.5_fp*dz(k,iens) + gllPts_ord(kk)*dz(k,iens);
              real yloc;
              if (sim2d) {
                yloc = ylen/2;
              } else {
                yloc = (j+0.5_fp)*dy + gllPts_ord(jj)*dy;
              }
              real xloc = (i+0.5_fp)*dx + gllPts_ord(ii)*dx;

              // Get dry constants

              // Compute constant theta hydrostatic background state
              real th = 300;
              real rh = profiles::initConstTheta_density(th,zloc,Rd,cp,gamma,p0,C0,grav);

              // Initialize tracer mass based on dry state
              // Vapor perturbation profile
              real pert  = profiles::ellipsoid_linear(xloc,yloc,zloc  ,  xlen/2,ylen/2,2000  ,  2000,2000,2000  ,  0.8);
              real press = C0*pow(rh*th,gamma);                       // Dry pressure
              real temp  = press / Rd / rh;                           // Temperator (same for dry and moist)
              real svp   = profiles::saturation_vapor_pressure(temp); // Self-explanatory
              real p_v   = pert*svp;                                  // Multiply profile by saturation vapor pressure
              real r_v   = p_v / (Rv*temp);                           // Compute vapor density

              real wt = gllWts_ord(kk) * gllWts_ord(jj) * gllWts_ord(ii);
              tracers(idWV,hs+k,hs+j,hs+i,iens) += r_v / (rh+r_v) * rh * wt;
              for (int tr=0; tr < num_tracers; tr++) {
                if (tr != idWV) tracers(tr,hs+k,hs+j,hs+i,iens) = 0;
              }
            }
          }
        }
      });

      parallel_for( "Spatial.h adjust_moisture" , SimpleBounds<4>(nz,ny,nx,nens) ,
                    YAKL_LAMBDA (int k, int j, int i, int iens) {
        // Add tracer density to dry density if it adds mass
        real rho_dry = state(idR,hs+k,hs+j,hs+i,iens);
        state(idR,hs+k,hs+j,hs+i,iens) += tracers(idWV,hs+k,hs+j,hs+i,iens);
        real rho_moist = state(idR,hs+k,hs+j,hs+i,iens);

        // Adjust momenta for moist density
        state(idU,hs+k,hs+j,hs+i,iens) = state(idU,hs+k,hs+j,hs+i,iens) / rho_dry * rho_moist;
        state(idV,hs+k,hs+j,hs+i,iens) = state(idV,hs+k,hs+j,hs+i,iens) / rho_dry * rho_moist;
        state(idW,hs+k,hs+j,hs+i,iens) = state(idW,hs+k,hs+j,hs+i,iens) / rho_dry * rho_moist;

        // Compute the dry temperature (same as the moist temperature)
        real rho_theta_dry = state(idT,hs+k,hs+j,hs+i,iens);
        real press = C0*pow(rho_theta_dry,gamma);  // Dry pressure
        real temp  = press / Rd / rho_dry;         // Temp (same dry or moist)

        // Compute moist theta
        real rho_v = tracers(idWV,hs+k,hs+j,hs+i,iens);
        real R_moist = Rd * (rho_dry / rho_moist) + Rv * (rho_v / rho_moist);
        real press_moist = rho_moist * R_moist * temp;
        real rho_theta_moist = pow( press_moist / C0 , 1._fp/gamma );

        // Compute moist rho*theta
        state(idT,hs+k,hs+j,hs+i,iens) = rho_theta_moist;

        for (int tr = 0 ; tr < num_tracers ; tr++) {
          tracers(tr,hs+k,hs+j,hs+i,iens) = tracers(tr,hs+k,hs+j,hs+i,iens) / rho_dry * rho_moist;
        }
      });

    } // if (data_spec == DATA_SPEC_THERMAL)

    convert_dynamics_to_coupler_state( coupler , state , tracers );
  }



  YAKL_INLINE static int wrapx(int i, int ii, int nx) {
    int ret = i+ii;
    if (ret < hs+0   ) ret += nx;
    if (ret > hs+nx-1) ret -= nx;
    return ret;
  }



  void compute_tendencies_x( pam::PamCoupler const &coupler,
                             real5d const &state   , real5d const &state_tend  ,
                             real5d const &tracers , real5d const &tracer_tend ,
                             real &dt ) const {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;

    auto num_tracers = coupler.get_num_tracers();
    auto nz          = coupler.get_nz     ();
    auto ny          = coupler.get_ny     ();
    auto nx          = coupler.get_nx     ();
    auto dx          = coupler.get_dx     ();
    auto nens        = coupler.get_nens   ();
    auto Rd          = coupler.get_R_d    ();
    auto cp          = coupler.get_cp_d   ();
    auto p0          = coupler.get_p0     ();
    auto gamma       = coupler.get_gamma_d();
    auto kappa       = Rd/cp;
    auto C0          = pow( Rd * pow( p0 , -kappa ) , gamma );
    auto tracer_pos  = coupler.get_tracer_positivity_array();
    auto sim2d = ny == 1;

    YAKL_SCOPE( weno_scalars            , this->weno_scalars           );
    YAKL_SCOPE( weno_winds              , this->weno_winds             );
    YAKL_SCOPE( c2g                     , this->coefs_to_gll           );
    YAKL_SCOPE( s2g                     , this->sten_to_gll            );
    YAKL_SCOPE( s2c                     , this->sten_to_coefs          );
    YAKL_SCOPE( weno_recon_lower        , this->weno_recon_lower       );
    YAKL_SCOPE( idl                     , this->idl                    );
    YAKL_SCOPE( sigma                   , this->sigma                  );
    YAKL_SCOPE( derivMatrix             , this->derivMatrix            );

    real6d state_limits ("state_limits" ,num_state  ,2,nz,ny,nx+1,nens);
    real6d tracer_limits("tracer_limits",num_tracers,2,nz,ny,nx+1,nens);

    // Loop through all cells, reconstruct in x-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( "Spatial.h X recon" , SimpleBounds<4>(nz,ny,nx,nens) ,
                  YAKL_LAMBDA (int k, int j, int i, int iens) {
      SArray<real,2,nAder,ngll> r_DTs , ru_DTs;

      { // State
        SArray<real,2,nAder,ngll> rv_DTs , rw_DTs , rt_DTs, ruu_DTs , ruv_DTs , ruw_DTs , rut_DTs , rt_gamma_DTs;

        { // Recon
          SArray<real,1,ord>  stencil;
          SArray<real,1,ngll> gll;

          // Density
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idR,hs+k,hs+j,wrapx(i,ii,nx),iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
          // Add hydrostasis back on
          for (int ii=0; ii < ngll; ii++) { r_DTs(0,ii) = gll(ii); }

          // u
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idU,hs+k,hs+j,wrapx(i,ii,nx),iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
          for (int ii=0; ii < ngll; ii++) { ru_DTs(0,ii) = gll(ii); }

          // v
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idV,hs+k,hs+j,wrapx(i,ii,nx),iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
          for (int ii=0; ii < ngll; ii++) { rv_DTs(0,ii) = gll(ii); }

          // w
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idW,hs+k,hs+j,wrapx(i,ii,nx),iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
          for (int ii=0; ii < ngll; ii++) { rw_DTs(0,ii) = gll(ii); }

          // theta
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idT,hs+k,hs+j,wrapx(i,ii,nx),iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
          // Add hydrostasis back on
          for (int ii=0; ii < ngll; ii++) { rt_DTs(0,ii) = gll(ii); }
        }

        for (int ii=0; ii < ngll; ii++) {
          real r = r_DTs (0,ii);
          real u = ru_DTs(0,ii) / r;
          real v = rv_DTs(0,ii) / r;
          real w = rw_DTs(0,ii) / r;
          real t = rt_DTs(0,ii) / r;
          ruu_DTs     (0,ii) = r*u*u;
          ruv_DTs     (0,ii) = r*u*v;
          ruw_DTs     (0,ii) = r*u*w;
          rut_DTs     (0,ii) = r*u*t;
          rt_gamma_DTs(0,ii) = pow(r*t,gamma);
        }

        if (nAder > 1) {
          awfl::diffTransformEulerConsX( r_DTs , ru_DTs , rv_DTs , rw_DTs , rt_DTs , ruu_DTs , ruv_DTs , ruw_DTs ,
                                         rut_DTs , rt_gamma_DTs , derivMatrix , C0 , gamma , dx );
        }

        SArray<real,1,ngll> r_tavg, ru_tavg;
        awfl::compute_timeAvg( r_DTs  , r_tavg  , dt );
        awfl::compute_timeAvg( ru_DTs , ru_tavg , dt );
        awfl::compute_timeAvg( rv_DTs           , dt );
        awfl::compute_timeAvg( rw_DTs           , dt );
        awfl::compute_timeAvg( rt_DTs           , dt );

        // Left interface
        state_limits(idR,1,k,j,i  ,iens) = r_tavg  (0     );
        state_limits(idU,1,k,j,i  ,iens) = ru_tavg (0     );
        state_limits(idV,1,k,j,i  ,iens) = rv_DTs(0,0     );
        state_limits(idW,1,k,j,i  ,iens) = rw_DTs(0,0     );
        state_limits(idT,1,k,j,i  ,iens) = rt_DTs(0,0     );
        // Right interface
        state_limits(idR,0,k,j,i+1,iens) = r_tavg  (ngll-1);
        state_limits(idU,0,k,j,i+1,iens) = ru_tavg (ngll-1);
        state_limits(idV,0,k,j,i+1,iens) = rv_DTs(0,ngll-1);
        state_limits(idW,0,k,j,i+1,iens) = rw_DTs(0,ngll-1);
        state_limits(idT,0,k,j,i+1,iens) = rt_DTs(0,ngll-1);
      }

      { // Tracers
        for (int tr=0; tr < num_tracers; tr++) {
          SArray<real,2,nAder,ngll> rt_DTs, rut_DTs;

          { // Recon
            SArray<real,1,ord>  stencil;
            SArray<real,1,ngll> gll;

            for (int ii=0; ii < ord; ii++) { stencil(ii) = tracers(tr,hs+k,hs+j,wrapx(i,ii,nx),iens); }
            awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
            if (tracer_pos(tr)) {
              for (int ii=0; ii < ngll; ii++) { gll(ii) = max( 0._fp , gll(ii) ); }
            }
            for (int ii=0; ii < ngll; ii++) { rt_DTs(0,ii) = gll(ii); }
          }

          for (int ii=0; ii < ngll; ii++) {
            rut_DTs(0,ii) = rt_DTs(0,ii) * ru_DTs(0,ii) / r_DTs(0,ii);
          }

          if (nAder > 1) {
            awfl::diffTransformTracer( r_DTs , ru_DTs , rt_DTs , rut_DTs , derivMatrix , dx );
          }

          awfl::compute_timeAvg( rt_DTs  , dt );

          if (tracer_pos(tr)) {
            for (int ii=0; ii < ngll; ii++) { rt_DTs(0,ii) = max( 0._fp , rt_DTs(0,ii) ); }
          }

          tracer_limits(tr,1,k,j,i  ,iens) = rt_DTs(0,0     ); // Left interface
          tracer_limits(tr,0,k,j,i+1,iens) = rt_DTs(0,ngll-1); // Right interface
        }
      }
    });

    real5d state_flux ("state_flux" ,num_state  ,nz,ny,nx+1,nens);
    real5d tracer_flux("tracer_flux",num_tracers,nz,ny,nx+1,nens);

    awfl::riemann_rho_theta_full_x( coupler , state_limits , state_flux , tracer_limits , tracer_flux );

    state_limits  = real6d();
    tracer_limits = real6d();

    awfl::fct_positivity_x( coupler , tracers , tracer_flux , tracer_pos , dt );

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( "Spatial.h X tendencies" , SimpleBounds<4>(nz,ny,nx,nens) ,
                  YAKL_LAMBDA(int k, int j, int i, int iens) {
      for (int l = 0; l < num_state; l++) {
        if (sim2d && l == idV) {
          state_tend(l,k,j,i,iens) = 0;
        } else {
          state_tend(l,k,j,i,iens) = - ( state_flux(l,k,j,i+1,iens) - state_flux(l,k,j,i,iens) ) / dx;
        }
      }
      for (int l = 0; l < num_tracers; l++) {
        // Compute tracer tendency
        tracer_tend(l,k,j,i,iens) = - ( tracer_flux(l,k,j,i+1,iens) - tracer_flux(l,k,j,i,iens) ) / dx;
      }
    });
  }



  YAKL_INLINE static int wrapy(int j, int jj, int ny) {
    int ret = j+jj;
    if (ret < hs+0   ) ret += ny;
    if (ret > hs+ny-1) ret -= ny;
    return ret;
  }



  void compute_tendencies_y( pam::PamCoupler const &coupler ,
                             real5d const &state   , real5d const &state_tend  ,
                             real5d const &tracers , real5d const &tracer_tend ,
                             real &dt ) const {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;

    auto num_tracers = coupler.get_num_tracers();
    auto nz          = coupler.get_nz     ();
    auto ny          = coupler.get_ny     ();
    auto nx          = coupler.get_nx     ();
    auto dy          = coupler.get_dy     ();
    auto nens        = coupler.get_nens   ();
    auto grav        = coupler.get_grav   ();
    auto Rd          = coupler.get_R_d    ();
    auto cp          = coupler.get_cp_d   ();
    auto p0          = coupler.get_p0     ();
    auto gamma       = coupler.get_gamma_d();
    auto kappa       = Rd/cp;
    auto C0          = pow( Rd * pow( p0 , -kappa ) , gamma );
    auto tracer_pos  = coupler.get_tracer_positivity_array();
    auto sim2d = ny == 1;

    if (sim2d) {
      state_tend  = 0;
      tracer_tend = 0;
      return;
    }

    YAKL_SCOPE( weno_scalars            , this->weno_scalars           );
    YAKL_SCOPE( weno_winds              , this->weno_winds             );
    YAKL_SCOPE( c2g                     , this->coefs_to_gll           );
    YAKL_SCOPE( s2g                     , this->sten_to_gll            );
    YAKL_SCOPE( s2c                     , this->sten_to_coefs          );
    YAKL_SCOPE( weno_recon_lower        , this->weno_recon_lower       );
    YAKL_SCOPE( idl                     , this->idl                    );
    YAKL_SCOPE( sigma                   , this->sigma                  );
    YAKL_SCOPE( derivMatrix             , this->derivMatrix            );

    real6d state_limits ("state_limits" ,num_state  ,2,nz,ny+1,nx,nens);
    real6d tracer_limits("tracer_limits",num_tracers,2,nz,ny+1,nx,nens);

    // Loop through all cells, reconstruct in y-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( "Spatial.h Y recon" , SimpleBounds<4>(nz,ny,nx,nens) ,
                  YAKL_LAMBDA (int k, int j, int i, int iens) {
      SArray<real,2,nAder,ngll> r_DTs , rv_DTs;

      { // State
        SArray<real,2,nAder,ngll> ru_DTs , rw_DTs , rt_DTs, rvu_DTs , rvv_DTs , rvw_DTs , rvt_DTs , rt_gamma_DTs;

        { // Recon
          SArray<real,1,ord>  stencil;
          SArray<real,1,ngll> gll;

          // Density
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idR,hs+k,wrapy(j,jj,ny),hs+i,iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
          // Add hydrostasis back on
          for (int jj=0; jj < ngll; jj++) { r_DTs(0,jj) = gll(jj); }

          // u
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idU,hs+k,wrapy(j,jj,ny),hs+i,iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
          for (int jj=0; jj < ngll; jj++) { ru_DTs(0,jj) = gll(jj); }

          // v
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idV,hs+k,wrapy(j,jj,ny),hs+i,iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
          for (int jj=0; jj < ngll; jj++) { rv_DTs(0,jj) = gll(jj); }

          // w
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idW,hs+k,wrapy(j,jj,ny),hs+i,iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
          for (int jj=0; jj < ngll; jj++) { rw_DTs(0,jj) = gll(jj); }

          // theta
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idT,hs+k,wrapy(j,jj,ny),hs+i,iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
          // Add hydrostasis back on
          for (int jj=0; jj < ngll; jj++) { rt_DTs(0,jj) = gll(jj); }
        }

        for (int jj=0; jj < ngll; jj++) {
          real r = r_DTs (0,jj);
          real u = ru_DTs(0,jj) / r;
          real v = rv_DTs(0,jj) / r;
          real w = rw_DTs(0,jj) / r;
          real t = rt_DTs(0,jj) / r;
          rvu_DTs     (0,jj) = r*v*u;
          rvv_DTs     (0,jj) = r*v*v;
          rvw_DTs     (0,jj) = r*v*w;
          rvt_DTs     (0,jj) = r*v*t;
          rt_gamma_DTs(0,jj) = pow(r*t,gamma);
        }

        if (nAder > 1) {
          awfl::diffTransformEulerConsY( r_DTs , ru_DTs , rv_DTs , rw_DTs , rt_DTs , rvu_DTs , rvv_DTs , rvw_DTs ,
                                         rvt_DTs , rt_gamma_DTs , derivMatrix , C0 , gamma , dy );
        }

        SArray<real,1,ngll> r_tavg, rv_tavg;
        awfl::compute_timeAvg( r_DTs  , r_tavg  , dt );
        awfl::compute_timeAvg( ru_DTs           , dt );
        awfl::compute_timeAvg( rv_DTs , rv_tavg , dt );
        awfl::compute_timeAvg( rw_DTs           , dt );
        awfl::compute_timeAvg( rt_DTs           , dt );

        // Left interface
        state_limits(idR,1,k,j  ,i,iens) = r_tavg  (0     );
        state_limits(idU,1,k,j  ,i,iens) = ru_DTs(0,0     );
        state_limits(idV,1,k,j  ,i,iens) = rv_tavg (0     );
        state_limits(idW,1,k,j  ,i,iens) = rw_DTs(0,0     );
        state_limits(idT,1,k,j  ,i,iens) = rt_DTs(0,0     );
        // Right interface
        state_limits(idR,0,k,j+1,i,iens) = r_tavg  (ngll-1);
        state_limits(idU,0,k,j+1,i,iens) = ru_DTs(0,ngll-1);
        state_limits(idV,0,k,j+1,i,iens) = rv_tavg (ngll-1);
        state_limits(idW,0,k,j+1,i,iens) = rw_DTs(0,ngll-1);
        state_limits(idT,0,k,j+1,i,iens) = rt_DTs(0,ngll-1);
      }

      { // Tracers
        for (int tr=0; tr < num_tracers; tr++) {
          SArray<real,2,nAder,ngll> rt_DTs, rvt_DTs;
          
          { // Recon
            SArray<real,1,ord>  stencil;
            SArray<real,1,ngll> gll;

            for (int jj=0; jj < ord; jj++) { stencil(jj) = tracers(tr,hs+k,wrapy(j,jj,ny),hs+i,iens); }
            awfl::reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
            if (tracer_pos(tr)) {
              for (int jj=0; jj < ngll; jj++) { gll(jj) = max( 0._fp , gll(jj) ); }
            }
            for (int jj=0; jj < ngll; jj++) { rt_DTs(0,jj) = gll(jj); }
          }

          for (int jj=0; jj < ngll; jj++) {
            rvt_DTs(0,jj) = rt_DTs(0,jj) * rv_DTs(0,jj) / r_DTs(0,jj);
          }

          if (nAder > 1) {
            awfl::diffTransformTracer( r_DTs , rv_DTs , rt_DTs , rvt_DTs , derivMatrix , dy );
          }

          awfl::compute_timeAvg( rt_DTs  , dt );

          if (tracer_pos(tr)) {
            for (int jj=0; jj < ngll; jj++) { rt_DTs(0,jj) = max( 0._fp , rt_DTs(0,jj) ); }
          }

          tracer_limits(tr,1,k,j  ,i,iens) = rt_DTs (0,0     ); // Left interface
          tracer_limits(tr,0,k,j+1,i,iens) = rt_DTs (0,ngll-1); // Right interface
        }
      }
    });

    real5d state_flux ("state_flux" ,num_state  ,nz,ny+1,nx,nens);
    real5d tracer_flux("tracer_flux",num_tracers,nz,ny+1,nx,nens);

    awfl::riemann_rho_theta_full_y( coupler , state_limits , state_flux , tracer_limits , tracer_flux );

    state_limits  = real6d();
    tracer_limits = real6d();

    awfl::fct_positivity_y( coupler, tracers, tracer_flux , tracer_pos, dt);

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( "Spatial.h Y tendendies" , SimpleBounds<4>(nz,ny,nx,nens) ,
                  YAKL_LAMBDA(int k, int j, int i, int iens) {
      for (int l=0; l < num_state; l++) {
        state_tend(l,k,j,i,iens) = - ( state_flux(l,k,j+1,i,iens) - state_flux(l,k,j,i,iens) ) / dy;
      }
      for (int l=0; l < num_tracers; l++) {
        // Compute the tracer tendency
        tracer_tend(l,k,j,i,iens) = - ( tracer_flux(l,k,j+1,i,iens) - tracer_flux(l,k,j,i,iens) ) / dy;
      }
    });
  }



  YAKL_INLINE static int wrapz(int k, int kk, int nz) {
    int ret = k+kk;
    if (ret < hs+0   ) ret = hs+0;
    if (ret > hs+nz-1) ret = hs+nz-1;
    return ret;
  }



  void compute_tendencies_z( pam::PamCoupler const &coupler ,
                             real5d const &state   , real5d const &state_tend  ,
                             real5d const &tracers , real5d const &tracer_tend ,
                             real &dt ) const {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;

    auto num_tracers = coupler.get_num_tracers();
    auto nz          = coupler.get_nz  ();
    auto ny          = coupler.get_ny  ();
    auto nx          = coupler.get_nx  ();
    auto dz          = coupler.get_data_manager_readonly().get<real const,2>("vertical_cell_dz");
    auto nens        = coupler.get_nens();
    auto grav        = coupler.get_grav   ();
    auto Rd          = coupler.get_R_d    ();
    auto cp          = coupler.get_cp_d   ();
    auto p0          = coupler.get_p0     ();
    auto gamma       = coupler.get_gamma_d();
    auto kappa       = Rd/cp;
    auto C0          = pow( Rd * pow( p0 , -kappa ) , gamma );
    auto tracer_pos  = coupler.get_tracer_positivity_array();
    auto sim2d = ny == 1;

    YAKL_SCOPE( weno_scalars            , this->weno_scalars           );
    YAKL_SCOPE( weno_winds              , this->weno_winds             );
    YAKL_SCOPE( c2g                     , this->coefs_to_gll           );
    YAKL_SCOPE( idl                     , this->idl                    );
    YAKL_SCOPE( sigma                   , this->sigma                  );
    YAKL_SCOPE( hyDensSten              , this->hyDensSten             );
    YAKL_SCOPE( hyDensThetaSten         , this->hyDensThetaSten        );
    YAKL_SCOPE( hyDensGLL               , this->hyDensGLL              );
    YAKL_SCOPE( hyDensThetaGLL          , this->hyDensThetaGLL         );
    YAKL_SCOPE( hyPressureGLL           , this->hyPressureGLL          );
    YAKL_SCOPE( derivMatrix             , this->derivMatrix            );
    YAKL_SCOPE( gllWts_ngll             , this->gllWts_ngll            );
    YAKL_SCOPE( vert_sten_to_gll        , this->vert_sten_to_gll       );
    YAKL_SCOPE( vert_sten_to_coefs      , this->vert_sten_to_coefs     );
    YAKL_SCOPE( vert_weno_recon_lower   , this->vert_weno_recon_lower  );

    real6d state_limits ("state_limits" ,num_state  ,2,nz+1,ny,nx,nens);
    real6d tracer_limits("tracer_limits",num_tracers,2,nz+1,ny,nx,nens);

    // Loop through all cells, reconstruct in x-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( "Spatial.h Z recon" , SimpleBounds<4>(nz,ny,nx,nens) ,
                  YAKL_LAMBDA (int k, int j, int i, int iens) {
      SArray<real,2,ord,ngll>       s2g_loc;
      SArray<real,2,ord,ord>        s2c_loc;
      SArray<real,3,hs+1,hs+1,hs+1> weno_recon_lower_loc;
      for (int jj=0; jj < ord; jj++) {
        for (int ii=0; ii < ngll; ii++) {
          s2g_loc(jj,ii) = vert_sten_to_gll(k,jj,ii,iens);
        }
      }
      for (int jj=0; jj < ord; jj++) {
        for (int ii=0; ii < ord; ii++) {
          s2c_loc(jj,ii) = vert_sten_to_coefs(k,jj,ii,iens);
        }
      }
      for (int kk=0; kk < hs+1; kk++) {
        for (int jj=0; jj < hs+1; jj++) {
          for (int ii=0; ii < hs+1; ii++) {
            weno_recon_lower_loc(kk,jj,ii) = vert_weno_recon_lower(k,kk,jj,ii,iens);
          }
        }
      }

      SArray<real,2,nAder,ngll> r_DTs , rw_DTs;

      { // State
        SArray<real,2,nAder,ngll> ru_DTs , rv_DTs , rt_DTs, rwu_DTs , rwv_DTs , rww_DTs , rwt_DTs , rt_gamma_DTs;
        { // Recon
          SArray<real,1,ord>  stencil;
          SArray<real,1,ngll> gll;

          // Density
          for (int kk=0; kk < ord; kk++) {
            if (k+kk < hs || k+kk > hs+nz-1) {
              stencil(kk) = state(idR,hs+k,hs+j,hs+i,iens) - hyDensSten(k,hs,iens);
            } else {
              stencil(kk) = state(idR,k+kk,hs+j,hs+i,iens) - hyDensSten(k,kk,iens);
            }
          }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g_loc , s2c_loc , weno_recon_lower_loc ,
                                  idl , sigma , weno_scalars );
          // Add hydrostasis back on
          for (int kk=0; kk < ngll; kk++) { r_DTs(0,kk) = gll(kk) + hyDensGLL(k,kk,iens); }

          // u values and derivatives
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idU,wrapz(k,kk,nz),hs+j,hs+i,iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g_loc , s2c_loc , weno_recon_lower_loc ,
                                  idl , sigma , weno_winds );
          for (int kk=0; kk < ngll; kk++) { ru_DTs(0,kk) = gll(kk); }

          // v
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idV,wrapz(k,kk,nz),hs+j,hs+i,iens); }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g_loc , s2c_loc , weno_recon_lower_loc ,
                                  idl , sigma , weno_winds );
          for (int kk=0; kk < ngll; kk++) { rv_DTs(0,kk) = gll(kk); }

          // w
          for (int kk=0; kk < ord; kk++) {
            stencil(kk) = state(idW,wrapz(k,kk,nz),hs+j,hs+i,iens);
            if (k+kk > hs+nz-1 || k+kk < hs) stencil(kk) = 0;
          }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g_loc , s2c_loc , weno_recon_lower_loc ,
                                  idl , sigma , weno_winds );
          if (k == nz-1) gll(ngll-1) = 0;
          if (k == 0   ) gll(0     ) = 0;
          for (int kk=0; kk < ngll; kk++) { rw_DTs(0,kk) = gll(kk); }

          // rho*theta
          for (int kk=0; kk < ord; kk++) {
            if (k+kk < hs || k+kk > hs+nz-1) {
              stencil(kk) = state(idT,hs+k,hs+j,hs+i,iens) - hyDensThetaSten(k,hs,iens);
            } else {
              stencil(kk) = state(idT,k+kk,hs+j,hs+i,iens) - hyDensThetaSten(k,kk,iens);
            }
          }
          awfl::reconstruct_gll_values( stencil , gll , c2g , s2g_loc , s2c_loc , weno_recon_lower_loc ,
                                  idl , sigma , weno_scalars );
          // Add hydrostasis back on
          for (int kk=0; kk < ngll; kk++) { rt_DTs(0,kk) = gll(kk) + hyDensThetaGLL(k,kk,iens); }
        }

        for (int kk=0; kk < ngll; kk++) {
          real r = r_DTs (0,kk);
          real u = ru_DTs(0,kk) / r;
          real v = rv_DTs(0,kk) / r;
          real w = rw_DTs(0,kk) / r;
          real t = rt_DTs(0,kk) / r;
          rwu_DTs     (0,kk) = r*w*u;
          rwv_DTs     (0,kk) = r*w*v;
          rww_DTs     (0,kk) = r*w*w;
          rwt_DTs     (0,kk) = r*w*t;
          rt_gamma_DTs(0,kk) = pow(r*t,gamma);
        }

        if (nAder > 1) {
          awfl::diffTransformEulerConsZ( r_DTs , ru_DTs , rv_DTs , rw_DTs , rt_DTs , rwu_DTs , rwv_DTs , rww_DTs ,
                                         rwt_DTs , rt_gamma_DTs , derivMatrix , hyDensGLL , hyPressureGLL , C0 , gamma ,
                                         grav , k , dz(k,iens) , nz , iens );
        }

        SArray<real,1,ngll> r_tavg, rw_tavg;
        awfl::compute_timeAvg( r_DTs  , r_tavg  , dt );
        awfl::compute_timeAvg( ru_DTs           , dt );
        awfl::compute_timeAvg( rv_DTs           , dt );
        awfl::compute_timeAvg( rw_DTs , rw_tavg , dt );
        awfl::compute_timeAvg( rt_DTs           , dt );

        // Left interface
        state_limits(idR,1,k  ,j,i,iens) = r_tavg  (0     );
        state_limits(idU,1,k  ,j,i,iens) = ru_DTs(0,0     );
        state_limits(idV,1,k  ,j,i,iens) = rv_DTs(0,0     );
        state_limits(idW,1,k  ,j,i,iens) = rw_tavg (0     );
        state_limits(idT,1,k  ,j,i,iens) = rt_DTs(0,0     );
        // Right interface
        state_limits(idR,0,k+1,j,i,iens) = r_tavg  (ngll-1);
        state_limits(idU,0,k+1,j,i,iens) = ru_DTs(0,ngll-1);
        state_limits(idV,0,k+1,j,i,iens) = rv_DTs(0,ngll-1);
        state_limits(idW,0,k+1,j,i,iens) = rw_tavg (ngll-1);
        state_limits(idT,0,k+1,j,i,iens) = rt_DTs(0,ngll-1);

        real ravg = 0;
        for (int kk=0; kk < ngll; kk++) {
          ravg += r_tavg(kk) * gllWts_ngll(kk);
        }
        state_tend(idR,k,j,i,iens) = 0;
        state_tend(idU,k,j,i,iens) = 0;
        state_tend(idV,k,j,i,iens) = 0;
        state_tend(idW,k,j,i,iens) = -grav * ravg;
        state_tend(idT,k,j,i,iens) = 0;
      }

      { // Tracers
        for (int tr=0; tr < num_tracers; tr++) {
          SArray<real,2,nAder,ngll> rt_DTs;  // Density * tracer
          SArray<real,2,nAder,ngll> rwt_DTs; // Density * wwind * tracer
          { // Recon
            SArray<real,1,ord>  stencil;
            SArray<real,1,ngll> gll;

            for (int kk=0; kk < ord; kk++) { stencil(kk) = tracers(tr ,wrapz(k,kk,nz),hs+j,hs+i,iens) /
                                                           state  (idR,wrapz(k,kk,nz),hs+j,hs+i,iens); }
            awfl::reconstruct_gll_values( stencil , gll , c2g , s2g_loc , s2c_loc , weno_recon_lower_loc ,
                                          idl , sigma , weno_scalars );
            for (int kk=0; kk < ngll; kk++) { gll(kk) *= r_DTs(0,kk); }
            if (tracer_pos(tr)) {
              for (int kk=0; kk < ngll; kk++) { gll(kk) = max( 0._fp , gll(kk) ); }
            }
            for (int kk=0; kk < ngll; kk++) { rt_DTs(0,kk) = gll(kk); }
          }

          for (int kk=0; kk < ngll; kk++) {
            rwt_DTs(0,kk) = rt_DTs(0,kk) * rw_DTs(0,kk) / r_DTs(0,kk);
          }

          if (nAder > 1) {
            awfl::diffTransformTracer( r_DTs , rw_DTs , rt_DTs , rwt_DTs , derivMatrix , dz(k,iens) );
          }

          awfl::compute_timeAvg( rt_DTs  , dt );

          if (tracer_pos(tr)) {
            for (int kk=0; kk < ngll; kk++) { rt_DTs(0,kk) = max( 0._fp , rt_DTs(0,kk) ); }
          }

          if (k == nz-1) rwt_DTs(0,ngll-1) = 0;
          if (k == 0   ) rwt_DTs(0,0     ) = 0;

          tracer_limits(tr,1,k  ,j,i,iens) = rt_DTs (0,0     ); // Left interface
          tracer_limits(tr,0,k+1,j,i,iens) = rt_DTs (0,ngll-1); // Right interface
        }
      }
    });

    real5d state_flux ("state_flux" ,num_state  ,nz+1,ny,nx,nens);
    real5d tracer_flux("tracer_flux",num_tracers,nz+1,ny,nx,nens);

    awfl::riemann_rho_theta_full_z(coupler, state_limits , state_flux , tracer_limits , tracer_flux );

    state_limits  = real6d();
    tracer_limits = real6d();

    awfl::fct_positivity_z( coupler , tracers , tracer_flux , tracer_pos , dt );

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( "Spatial.h Z tendencies" , SimpleBounds<4>(nz,ny,nx,nens) ,
                  YAKL_LAMBDA(int k, int j, int i, int iens) {
      for (int l=0; l < num_state; l++) {
        if (sim2d && l == idV) {
          state_tend(l,k,j,i,iens) = 0;
        } else {
          state_tend(l,k,j,i,iens) += - ( state_flux(l,k+1,j,i,iens) - state_flux(l,k,j,i,iens) ) / dz(k,iens);
        }
      }
      for (int l=0; l < num_tracers; l++) {
        // Compute tracer tendency
        tracer_tend(l,k,j,i,iens) = - ( tracer_flux(l,k+1,j,i,iens) - tracer_flux(l,k,j,i,iens) ) / dz(k,iens);
      }
    });
  }



  const char * getName() { return ""; }



  void finalize(real4d const &state , real4d const &tracers) {}


};

