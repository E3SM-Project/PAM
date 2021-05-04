
#pragma once

#include "awfl_const.h"
#include "phys_params.h"
#include "TransformMatrices.h"
#include "TransformMatrices_variable.h"
#include "WenoLimiter.h"
#include "Profiles.h"
#include "DataManager.h"


template <int nTimeDerivs, bool timeAvg, int nAder>
class Spatial_operator {
public:

  static_assert(nTimeDerivs == 1 , "ERROR: This Spatial class isn't setup to use nTimeDerivs > 1");

  int static constexpr hs = (ord-1)/2;
  int static constexpr num_state = 5;
  int static constexpr max_tracers = 50;

  // Stores a single index location
  struct Location {
    int l;
    int k;
    int j;
    int i;
  };

  real Rd   ;
  real cp   ;
  real gamma;
  real p0   ;
  real C0   ;
  real Rv   ;

  typedef real5d StateArr;  // Array of state variables (rho, rho*u, rho*v, rho*w, and rho*theta)
  typedef real5d TracerArr; // Array of tracers (total tracer mass)

  typedef real5d StateTendArr;   // State tendencies
  typedef real5d TracerTendArr;  // Tracer tendencies

  // Stores two estimates of the state, state flux, and tracer values at each cell interface
  real6d stateLimits_x;
  real6d tracerLimits_x;
  real5d stateFlux_x;
  real5d tracerFlux_x;
  real6d stateLimits_y;
  real6d tracerLimits_y;
  real5d stateFlux_y;
  real5d tracerFlux_y;
  real6d stateLimits_z;
  real6d tracerLimits_z;
  real5d stateFlux_z;
  real5d tracerFlux_z;

  // Hydrostatically balanced values for density, potential temperature, and pressure (cell-averages)
  real2d hyDensCells;
  real2d hyPressureCells;
  real2d hyThetaCells;
  real2d hyDensThetaCells;

  // Hydrostatically balanced values for density, potential temperature, and pressure (GLL points)
  real3d hyDensGLL;
  real3d hyPressureGLL;
  real3d hyThetaGLL;
  real3d hyDensThetaGLL;

  // Matrices to transform DOFs from one form to another
  SArray<real,2,ord,ngll> coefs_to_gll;
  SArray<real,2,ord,ngll> coefs_to_deriv_gll;
  SArray<real,2,ord,ngll> sten_to_gll;
  SArray<real,2,ord,ngll> sten_to_deriv_gll;
  // WENO reconstruction matrices
  SArray<real,3,ord,ord,ord> wenoRecon;
  SArray<real,1,hs+2> idl;   // Ideal weights for WENO
  real sigma;                // WENO sigma parameter (handicap high-order TV estimate)
  // For ADER spatial derivative computation (ngll points -> coefs -> deriv -> ngll points)
  SArray<real,2,ngll,ngll> derivMatrix;
  // For quadrature
  SArray<real,1,ord> gllWts_ord;
  SArray<real,1,ord> gllPts_ord;
  SArray<real,1,ngll> gllWts_ngll;
  SArray<real,1,ngll> gllPts_ngll;

  real2d vert_interface;
  real2d vert_interface_ghost;
  real3d vert_locs_normalized;
  real2d dz;
  real2d dz_ghost;
  real4d vert_sten_to_gll;
  real5d vert_weno_recon;

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
  int static constexpr DATA_SPEC_THERMAL       = 1;
  int static constexpr DATA_SPEC_SUPERCELL     = 2;

  bool sim2d;  // Whether we're simulating in 2-D

  real sim_time;  // How long to simulate

  // Grid spacing in each dimension
  real dx;
  real dy;

  // Initial time step (used throughout the simulation)
  real dtInit;

  int num_tracers;  // Number of tracers
  std::vector<std::string> tracer_name;      // Name of each tracer
  std::vector<std::string> tracer_desc;      // Description of each tracer
  bool1d                   tracer_pos;       // Whether each tracer is positive-definite
  bool1d                   tracer_adds_mass; // Whether each tracer adds mass (otherwise it's passive)

  //////////////////////////////////////
  // Values read from input file
  //////////////////////////////////////
  // Number of ensembles to simulate at the same time
  int         nens;
  // Number of cells in each direction
  int         nx;
  int         ny;
  int         nz;
  // Length of the domain in each direction (m)
  real        xlen;
  real        ylen;
  real1d      zlen;
  // Boundary condition in each direction
  int         bc_x;
  int         bc_y;
  int         bc_z;
  // Whether to use WENO for scalars and also for winds
  bool        weno_scalars;
  bool        weno_winds;
  // Name of the output file
  std::string out_prefix;
  // How to initialize the data
  int         data_spec;


  // When this class is created, initialize num_tracers to zero
  Spatial_operator() {
    num_tracers = 0;
  }



  // Make sure it's odd-order-accurate
  static_assert(ord%2 == 1,"ERROR: ord must be an odd integer");



  template <class MICRO>
  void convert_dynamics_to_coupler_state( DataManager &dm , MICRO &micro ) {
    real5d state           = dm.get<real,5>( "dynamics_state"   );
    real5d tracers         = dm.get<real,5>( "dynamics_tracers" );
    real4d dm_dens_dry     = dm.get<real,4>( "density_dry"      );
    real4d dm_uvel         = dm.get<real,4>( "uvel"             );
    real4d dm_vvel         = dm.get<real,4>( "vvel"             );
    real4d dm_wvel         = dm.get<real,4>( "wvel"             );
    real4d dm_temp         = dm.get<real,4>( "temp"             );

    YAKL_SCOPE( hyDensCells      , this->hyDensCells      );
    YAKL_SCOPE( hyDensThetaCells , this->hyDensThetaCells );
    YAKL_SCOPE( C0               , this->C0               );
    YAKL_SCOPE( gamma            , this->gamma            );
    YAKL_SCOPE( num_tracers      , this->num_tracers      );
    YAKL_SCOPE( p0               , this->p0               );
    YAKL_SCOPE( Rd               , this->Rd               );
    YAKL_SCOPE( Rv               , this->Rv               );
    YAKL_SCOPE( cp               , this->cp               );
    YAKL_SCOPE( tracer_adds_mass , this->tracer_adds_mass );

    int idWV = micro.get_water_vapor_index();

    MultipleFields<max_tracers,real4d> dm_tracers;
    for (int tr = 0; tr < num_tracers; tr++) {
      auto trac = dm.get<real,4>( tracer_name[tr] );
      dm_tracers.add_field( trac );
    }

    parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      real dens  = state(idR,hs+k,hs+j,hs+i,iens) + hyDensCells(k,iens);
      real uvel  = state(idU,hs+k,hs+j,hs+i,iens) / dens;
      real vvel  = state(idV,hs+k,hs+j,hs+i,iens) / dens;
      real wvel  = state(idW,hs+k,hs+j,hs+i,iens) / dens;
      real theta = ( state(idT,hs+k,hs+j,hs+i,iens) + hyDensThetaCells(k,iens) ) / dens;
      real pressure = C0 * pow( dens*theta , gamma );
      real dens_vap = tracers(idWV,hs+k,hs+j,hs+i,iens);
      real dens_dry = dens;
      for (int tr=0; tr < num_tracers; tr++) {
        if (tracer_adds_mass(tr)) dens_dry -= tracers(tr,hs+k,hs+j,hs+i,iens);
      }
      real temp = pressure / ( dens_dry * Rd + dens_vap * Rv );
      dm_dens_dry    (k,j,i,iens) = dens_dry;
      dm_uvel        (k,j,i,iens) = uvel;
      dm_vvel        (k,j,i,iens) = vvel;
      dm_wvel        (k,j,i,iens) = wvel;
      dm_temp        (k,j,i,iens) = temp;
      for (int tr=0; tr < num_tracers; tr++) {
        dm_tracers(tr,k,j,i,iens) = tracers(tr,hs+k,hs+j,hs+i,iens);
      }
    });
  }



  template <class MICRO>
  void convert_coupler_state_to_dynamics( DataManager &dm , MICRO &micro ) {
    real5d state           = dm.get<real,5>( "dynamics_state"   );
    real5d tracers         = dm.get<real,5>( "dynamics_tracers" );
    real4d dm_dens_dry     = dm.get<real,4>( "density_dry"      );
    real4d dm_uvel         = dm.get<real,4>( "uvel"             );
    real4d dm_vvel         = dm.get<real,4>( "vvel"             );
    real4d dm_wvel         = dm.get<real,4>( "wvel"             );
    real4d dm_temp         = dm.get<real,4>( "temp"             );

    YAKL_SCOPE( hyDensCells      , this->hyDensCells      );
    YAKL_SCOPE( hyDensThetaCells , this->hyDensThetaCells );
    YAKL_SCOPE( C0               , this->C0               );
    YAKL_SCOPE( gamma            , this->gamma            );
    YAKL_SCOPE( num_tracers      , this->num_tracers      );
    YAKL_SCOPE( p0               , this->p0               );
    YAKL_SCOPE( Rd               , this->Rd               );
    YAKL_SCOPE( Rv               , this->Rv               );
    YAKL_SCOPE( cp               , this->cp               );
    YAKL_SCOPE( tracer_adds_mass , this->tracer_adds_mass );

    int idWV = micro.get_water_vapor_index();

    MultipleFields<max_tracers,real4d> dm_tracers;
    for (int tr = 0; tr < num_tracers; tr++) {
      auto trac = dm.get<real,4>( tracer_name[tr] );
      dm_tracers.add_field( trac );
    }

    parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      for (int tr=0; tr < num_tracers; tr++) {
        tracers(tr,hs+k,hs+j,hs+i,iens) = dm_tracers(tr,k,j,i,iens);
      }
      real dens_dry     = dm_dens_dry    (k,j,i,iens);
      real uvel         = dm_uvel        (k,j,i,iens);
      real vvel         = dm_vvel        (k,j,i,iens);
      real wvel         = dm_wvel        (k,j,i,iens);
      real temp         = dm_temp        (k,j,i,iens);
      real dens_vap     = tracers(idWV,hs+k,hs+j,hs+i,iens);
      real dens         = dens_dry;
      for (int tr=0; tr < num_tracers; tr++) {
        if (tracer_adds_mass(tr)) dens += tracers(tr,hs+k,hs+j,hs+i,iens);
      }
      real pressure     = dens_dry * Rd * temp + dens_vap * Rv * temp;
      real theta        = pow( pressure / C0 , 1._fp / gamma ) / dens;
      state(idR,hs+k,hs+j,hs+i,iens) = dens - hyDensCells(k,iens);
      state(idU,hs+k,hs+j,hs+i,iens) = dens * uvel;
      state(idV,hs+k,hs+j,hs+i,iens) = dens * vvel;
      state(idW,hs+k,hs+j,hs+i,iens) = dens * wvel;
      state(idT,hs+k,hs+j,hs+i,iens) = dens * theta - hyDensThetaCells(k,iens);
    });
  }



  // Initialize a tracer
  int add_tracer(DataManager &dm , std::string name , std::string desc , bool pos_def , bool adds_mass) {
    YAKL_SCOPE( tracer_pos       , this->tracer_pos       );
    YAKL_SCOPE( tracer_adds_mass , this->tracer_adds_mass );

    int tr = tracer_name.size();  // Index to insert this tracer at
    if (tr == num_tracers) {
      endrun("ERROR: adding more tracers than initially requested");
    }
    tracer_name.push_back(name);  // Store name
    tracer_desc.push_back(desc);  // Store description
    parallel_for( 1 , YAKL_LAMBDA (int i) {
      tracer_pos      (tr) = pos_def;   // Store whether it's positive-definite
      tracer_adds_mass(tr) = adds_mass; // Store whether it adds mass (otherwise it's passive)
    });

    // Return the index of this tracer to the caller
    return tr;
  }



  // Caller creates a lambda (init_mass) to initialize this tracer value using location and dry state information
  template <class MICRO>
  void init_tracers( DataManager &dm , MICRO const &micro) {
    YAKL_SCOPE( dx             , this->dx             );
    YAKL_SCOPE( dy             , this->dy             );
    YAKL_SCOPE( dz             , this->dz             );
    YAKL_SCOPE( gllPts_ord     , this->gllPts_ord     );
    YAKL_SCOPE( gllWts_ord     , this->gllWts_ord     );
    YAKL_SCOPE( sim2d          , this->sim2d          );
    YAKL_SCOPE( xlen           , this->xlen           );
    YAKL_SCOPE( ylen           , this->ylen           );
    YAKL_SCOPE( zlen           , this->zlen           );
    YAKL_SCOPE( Rd             , this->Rd             );
    YAKL_SCOPE( cp             , this->cp             );
    YAKL_SCOPE( gamma          , this->gamma          );
    YAKL_SCOPE( p0             , this->p0             );
    YAKL_SCOPE( C0             , this->C0             );
    YAKL_SCOPE( Rv             , this->Rv             );
    YAKL_SCOPE( vert_interface , this->vert_interface );

    int idWV = micro.get_water_vapor_index();
    real5d tracers = dm.get<real,5>("dynamics_tracers");

    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      tracers(idWV,hs+k,hs+j,hs+i,iens) = 0;
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
            real rh = profiles::initConstTheta_density(th,zloc,Rd,cp,gamma,p0,C0);

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
          }
        }
      }
    });
  }



  // Take an initially dry fluid state and adjust it to account for moist tracers
  template <class MICRO>
  void adjust_state_for_moisture(DataManager &dm , MICRO const &micro) const {
    YAKL_SCOPE( hyDensCells             , this->hyDensCells             );
    YAKL_SCOPE( hyDensThetaCells        , this->hyDensThetaCells        );
    YAKL_SCOPE( num_tracers             , this->num_tracers             );
    YAKL_SCOPE( tracer_adds_mass        , this->tracer_adds_mass        );
    YAKL_SCOPE( Rd                      , this->Rd                      );
    YAKL_SCOPE( Rv                      , this->Rv                      );
    YAKL_SCOPE( C0                      , this->C0                      );
    YAKL_SCOPE( gamma                   , this->gamma                   );

    // Copy the DataManager data to state and tracer arrays for convenience
    real5d state   = dm.get<real,5>("dynamics_state");
    real5d tracers = dm.get<real,5>("dynamics_tracers");

    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      // Add tracer density to dry density if it adds mass
      real rho_dry = state(idR,hs+k,hs+j,hs+i,iens) + hyDensCells(k,iens);
      for (int tr=0; tr < num_tracers; tr++) {
        if (tracer_adds_mass(tr)) {
          state(idR,hs+k,hs+j,hs+i,iens) += tracers(tr,hs+k,hs+j,hs+i,iens);
        }
      }
      real rho_moist = state(idR,hs+k,hs+j,hs+i,iens) + hyDensCells(k,iens);

      // Adjust momenta for moist density
      state(idU,hs+k,hs+j,hs+i,iens) = state(idU,hs+k,hs+j,hs+i,iens) / rho_dry * rho_moist;
      state(idV,hs+k,hs+j,hs+i,iens) = state(idV,hs+k,hs+j,hs+i,iens) / rho_dry * rho_moist;
      state(idW,hs+k,hs+j,hs+i,iens) = state(idW,hs+k,hs+j,hs+i,iens) / rho_dry * rho_moist;

      // Compute the dry temperature (same as the moist temperature)
      real rho_theta_dry = state(idT,hs+k,hs+j,hs+i,iens) + hyDensThetaCells(k,iens);
      real press = C0*pow(rho_theta_dry,gamma);  // Dry pressure
      real temp  = press / Rd / rho_dry;         // Temp (same dry or moist)

      // Compute moist theta
      real index_vapor = micro.get_water_vapor_index();
      real rho_v = tracers(index_vapor,hs+k,hs+j,hs+i,iens);
      real R_moist = Rd * (rho_dry / rho_moist) + Rv * (rho_v / rho_moist);
      real press_moist = rho_moist * R_moist * temp;
      real rho_theta_moist = pow( press_moist / C0 , 1._fp/gamma );

      // Compute moist rho*theta
      state(idT,hs+k,hs+j,hs+i,iens) = rho_theta_moist - hyDensThetaCells(k,iens);

      for (int tr = 0 ; tr < num_tracers ; tr++) {
        tracers(tr,hs+k,hs+j,hs+i,iens) = tracers(tr,hs+k,hs+j,hs+i,iens) / rho_dry * rho_moist;
      }

    });
  }



  real5d createStateArr() const {
    return real5d("stateArr",num_state,nz+2*hs,ny+2*hs,nx+2*hs,nens);
  }



  real5d createTracerArr() const {
    return real5d("tracerArr",num_tracers,nz+2*hs,ny+2*hs,nx+2*hs,nens);
  }



  real5d createStateTendArr() const {
    return real5d("stateTendArr",num_state,nz,ny,nx,nens);
  }



  real5d createTracerTendArr() const {
    return real5d("tracerTendArr",num_tracers,nz,ny,nx,nens);
  }



  // Number of operator splittinng steps to use
  // Normally this would be 3, but the z-directly CFL is reduced because of how the fluxes are
  // handled in the presence of a solid wall boundary condition. I'm looking into how to fix this
  int numSplit() const {
    return 1;
  }



  // Given the model data and CFL value, compute the maximum stable time step
  template <class MICRO>
  real compute_time_step(DataManager &dm, MICRO const &micro, real cfl = 0.8) {

    // If we've already computed the time step, then don't compute it again
    if (dtInit <= 0) {
      YAKL_SCOPE( dx                   , this->dx                  );
      YAKL_SCOPE( dy                   , this->dy                  );
      YAKL_SCOPE( dz                   , this->dz                  );
      YAKL_SCOPE( hyDensCells          , this->hyDensCells         );
      YAKL_SCOPE( hyDensThetaCells     , this->hyDensThetaCells    );
      YAKL_SCOPE( gamma                , this->gamma               );
      YAKL_SCOPE( C0                   , this->C0                  );

      // Convert data from DataManager to state and tracers array for convenience
      real5d state   = dm.get<real,5>("dynamics_state");
      real5d tracers = dm.get<real,5>("dynamics_tracers");

      // Allocate a 3-D array for the max stable time steps (we'll use this for a reduction later)
      real4d dt3d("dt3d",nz,ny,nx,nens);

      // Loop through the cells, calculate the max stable time step for each cell
      parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        // Get the state
        real r = state(idR,hs+k,hs+j,hs+i,iens) + hyDensCells(k,iens);
        real u = state(idU,hs+k,hs+j,hs+i,iens) / r;
        real v = state(idV,hs+k,hs+j,hs+i,iens) / r;
        real w = state(idW,hs+k,hs+j,hs+i,iens) / r;
        real t = ( state(idT,hs+k,hs+j,hs+i,iens) + hyDensThetaCells(k,iens) ) / r;
        real p = C0*pow(r*t,gamma);

        // Compute the speed of sound (constant kappa assumption)
        real cs = sqrt(gamma*p/r);

        // Compute the maximum stable time step in each direction
        real udt = cfl * dx         / max( abs(u-cs) , abs(u+cs) );
        real vdt = cfl * dy         / max( abs(v-cs) , abs(v+cs) );
        real wdt = cfl * dz(k,iens) / max( abs(w-cs) , abs(w+cs) );

        // Compute the min of the max stable time steps
        dt3d(k,j,i,iens) = min( min(udt,vdt) , wdt );
      });
      // Store to dtInit so we don't have to compute this again
      dtInit = yakl::intrinsics::minval( dt3d );
    }

    return dtInit;
  }



  // Initialize crap needed by recon()
  void init(std::string inFile, int ny, int nx, int nens, real xlen, real ylen, int num_tracers, DataManager &dm) {
    this->nens = nens;
    this->nx = nx;
    this->ny = ny;
    this->xlen = xlen;
    this->ylen = ylen;
    this->num_tracers = num_tracers;

    // Allocate device arrays for whether tracers are positive-definite or add mass
    tracer_pos       = bool1d("tracer_pos"      ,num_tracers);
    tracer_adds_mass = bool1d("tracer_adds_mass",num_tracers);

    // Inialize time step to zero, and dimensional splitting switch
    dtInit = 0;

    #if defined(PAM_STANDALONE)
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
      } else {
        endrun("ERROR: Invalid data_spec");
      }

      // Read the x-direction boundary condition option
      std::string bc_x_str = config["bc_x"].as<std::string>();
      if        (bc_x_str == "periodic" ) {
        bc_x = BC_PERIODIC;
      } else if (bc_x_str == "wall"     ) {
        bc_x = BC_WALL;
      } else {
        endrun("Invalid bc_x");
      }

      // Read the y-direction boundary condition option
      std::string bc_y_str = config["bc_y"].as<std::string>();
      if        (bc_y_str == "periodic" ) {
        bc_y = BC_PERIODIC;
      } else if (bc_y_str == "wall"     ) {
        bc_y = BC_WALL;
      } else {
        endrun("Invalid bc_y");
      }

      // Read the z-direction boundary condition option
      std::string bc_z_str = config["bc_z"].as<std::string>();
      if        (bc_z_str == "periodic" ) {
        bc_z = BC_PERIODIC;
      } else if (bc_z_str == "wall"     ) {
        bc_z = BC_WALL;
      } else {
        endrun("Invalid bc_z");
      }

      // Read the output filename
      out_prefix = config["out_prefix"].as<std::string>();

      sim_time = config["simTime"].as<real>();
    #else
      weno_scalars            = true;
      weno_winds              = true;
      data_spec               = DATA_SPEC_SUPERCELL;
      bc_x                    = BC_PERIODIC;
      bc_y                    = BC_PERIODIC;
      bc_z                    = BC_WALL;
      out_prefix              = "test";
      sim_time                = 900;
    #endif

    // Determine whether this is a 2-D simulation
    sim2d = ny == 1;

    // Store vertical cell interface heights in the data manager
    auto zint = dm.get<real,2>("vertical_interface_height");

    nz = dm.get_dimension_size("z");

    // Get the height of the z-dimension
    zlen = real1d("zlen",nens);
    YAKL_SCOPE( zlen , this->zlen );
    YAKL_SCOPE( nz   , this->nz   );
    parallel_for( nens , YAKL_LAMBDA (int iens) {
      zlen(iens) = zint(nz,iens);
    });

    vert_interface       = real2d("vert_interface"      ,nz+1          ,nens);
    vert_interface_ghost = real2d("vert_interface_ghost",nz+2*hs+1     ,nens);
    vert_locs_normalized = real3d("vert_locs_normalized",nz,ord+1      ,nens);
    dz                   = real2d("dz"                  ,nz            ,nens);
    dz_ghost             = real2d("dz_ghost"            ,nz+2*hs       ,nens);
    vert_sten_to_gll     = real4d("vert_sten_to_gll"    ,nz,ord,ngll   ,nens);
    vert_weno_recon      = real5d("vert_weno_recon"     ,nz,ord,ord,ord,nens);

    YAKL_SCOPE( vert_interface       , this->vert_interface       );
    YAKL_SCOPE( vert_interface_ghost , this->vert_interface_ghost );
    YAKL_SCOPE( vert_locs_normalized , this->vert_locs_normalized );
    YAKL_SCOPE( dz                   , this->dz                   );
    YAKL_SCOPE( dz_ghost             , this->dz_ghost             );
    YAKL_SCOPE( vert_sten_to_gll     , this->vert_sten_to_gll     );
    YAKL_SCOPE( vert_weno_recon      , this->vert_weno_recon      );

    zint.deep_copy_to(vert_interface);

    parallel_for( Bounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
      dz(k,iens) = vert_interface(k+1,iens) - vert_interface(k,iens);
    });

    parallel_for( Bounds<2>(nz+2*hs,nens) , YAKL_LAMBDA (int k, int iens) {
      if (k >= hs && k < hs+nz) {
        dz_ghost(k,iens) = dz(k-hs,iens);
      } else if (k < hs) {
        dz_ghost(k,iens) = dz(0,iens);
      } else if (k >= hs+nz) {
        dz_ghost(k,iens) = dz(nz-1,iens);
      }
    });

    parallel_for( nens , YAKL_LAMBDA (int iens) {
      vert_interface_ghost(0,iens) = vert_interface(0,iens) - hs*dz(0,iens);
      for (int k=1; k < nz+2*hs+1; k++) {
        vert_interface_ghost(k,iens) = vert_interface_ghost(k-1,iens) + dz_ghost(k-1,iens);
      }
    });

    auto vint_host      = vert_interface_ghost.createHostCopy();
    auto vert_s2g_host  = vert_sten_to_gll    .createHostCopy();
    auto vert_weno_host = vert_weno_recon     .createHostCopy();
    auto vert_locs_host = vert_locs_normalized.createHostCopy();

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
        SArray<double,3,ord,ord,ord> weno_recon_var;
        TransformMatrices_variable::sten_to_coefs_variable<ord>(locs,s2c_var_in);
        TransformMatrices_variable::weno_sten_to_coefs<ord>(locs,weno_recon_var);
        SArray<real,2,ord,ord> s2c_var;
        for (int jj=0; jj < ord; jj++) {
          for (int ii=0; ii < ord; ii++) {
            s2c_var(jj,ii) = s2c_var_in(jj,ii);
          }
        }
        auto s2g_var = c2g * s2c_var;

        // Store reconstruction matrices
        for (int jj=0; jj < ord; jj++) {
          for (int ii=0; ii < ngll; ii++) {
            vert_s2g_host(k,jj,ii,iens) = s2g_var(jj,ii);
          }
        }

        for (int kk=0; kk < ord; kk++) {
          for (int jj=0; jj < ord; jj++) {
            for (int ii=0; ii < ord; ii++) {
              vert_weno_host(k,kk,jj,ii,iens) = weno_recon_var(kk,jj,ii);
            }
          }
        }

      }
    }

    vert_s2g_host .deep_copy_to(vert_sten_to_gll    );
    vert_weno_host.deep_copy_to(vert_weno_recon     );
    vert_locs_host.deep_copy_to(vert_locs_normalized);

    // Compute the grid spacing in each dimension
    dx = xlen/nx;
    dy = ylen/ny;

    // Store the WENO reconstruction matrices
    TransformMatrices::weno_sten_to_coefs(this->wenoRecon);

    // Block exists to avoid name mangling stufff
    {
      SArray<real,2,ord,ord>  g2c;        // Converts ord GLL points to ord coefficients
      SArray<real,2,ord,ord>  s2c;        // Converts ord stencil cell averages to ord coefficients
      SArray<real,2,ord,ngll> c2g_lower;  // Converts ord coefficients to ngll GLL points
      SArray<real,2,ord,ord>  c2g;        // Converts ord coefficients to ord GLL points
      SArray<real,2,ord,ord>  c2d;        // Converts ord coefficients to order differentiated coefficients

      TransformMatrices::gll_to_coefs      (g2c      );
      TransformMatrices::sten_to_coefs     (s2c      );
      TransformMatrices::coefs_to_gll_lower(c2g_lower);
      TransformMatrices::coefs_to_gll      (c2g      );
      TransformMatrices::coefs_to_deriv    (c2d      );

      this->coefs_to_gll       = c2g_lower;              // Converts ord coefficients to ngll GLL points
      this->coefs_to_deriv_gll = c2g_lower * c2d;        // Converts ord coefficients to ngll differentiated GLL points
      this->sten_to_gll        = c2g_lower       * s2c;  // Converts ord stencil cell avgs to ngll GLL points
      this->sten_to_deriv_gll  = c2g_lower * c2d * s2c;  // Converts ord stencil cell avgs to ngll differentiated GLL points

    }
    // Store ader derivMatrix
    {
      SArray<real,2,ngll,ngll> g2c;  // Converts ngll GLL points to ngll coefficients
      SArray<real,2,ngll,ngll> c2d;  // Converts ngll coefficients to ngll differentiated coefficients
      SArray<real,2,ngll,ngll> c2g;  // Converts ngll coefficients to ngll GLL points

      TransformMatrices::gll_to_coefs  (g2c);
      TransformMatrices::coefs_to_deriv(c2d);
      TransformMatrices::coefs_to_gll  (c2g);

      this->derivMatrix = c2g * c2d * g2c;  // Converts ngll GLL points to ngll differentiated GLL points
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
    stateLimits_x    = real6d("stateLimits_x"     ,num_state  ,2,nz,ny,nx+1,nens);
    tracerLimits_x   = real6d("tracerLimits_x"    ,num_tracers,2,nz,ny,nx+1,nens);
    stateFlux_x      = real5d("stateFlux_x"       ,num_state    ,nz,ny,nx+1,nens);
    tracerFlux_x     = real5d("tracerFlux_x"      ,num_tracers  ,nz,ny,nx+1,nens);
    stateLimits_y    = real6d("stateLimits_y"     ,num_state  ,2,nz,ny+1,nx,nens);
    tracerLimits_y   = real6d("tracerLimits_y"    ,num_tracers,2,nz,ny+1,nx,nens);
    stateFlux_y      = real5d("stateFlux_y"       ,num_state    ,nz,ny+1,nx,nens);
    tracerFlux_y     = real5d("tracerFlux_y"      ,num_tracers  ,nz,ny+1,nx,nens);
    stateLimits_z    = real6d("stateLimits_z"     ,num_state  ,2,nz+1,ny,nx,nens);
    tracerLimits_z   = real6d("tracerLimits_z"    ,num_tracers,2,nz+1,ny,nx,nens);
    stateFlux_z      = real5d("stateFlux_z"       ,num_state    ,nz+1,ny,nx,nens);
    tracerFlux_z     = real5d("tracerFlux_z"      ,num_tracers  ,nz+1,ny,nx,nens);
    hyDensCells      = real2d("hyDensCells       ",nz,nens);
    hyPressureCells  = real2d("hyPressureCells   ",nz,nens);
    hyThetaCells     = real2d("hyThetaCells      ",nz,nens);
    hyDensThetaCells = real2d("hyDensThetaCells  ",nz,nens);
    hyDensGLL        = real3d("hyDensGLL         ",nz,ngll,nens);
    hyPressureGLL    = real3d("hyPressureGLL     ",nz,ngll,nens);
    hyThetaGLL       = real3d("hyThetaGLL        ",nz,ngll,nens);
    hyDensThetaGLL   = real3d("hyDensThetaGLL    ",nz,ngll,nens);

    // Register and allocate state data with the DataManager
    dm.register_and_allocate<real>( "dynamics_state"   , "dynamics state"   , {num_state  ,nz+2*hs,ny+2*hs,nx+2*hs,nens} , {"num_state"  ,"nz_halo","ny_halo","nx_halo","nens"} );
    dm.register_and_allocate<real>( "dynamics_tracers" , "dynamics tracers" , {num_tracers,nz+2*hs,ny+2*hs,nx+2*hs,nens} , {"num_tracers","nz_halo","ny_halo","nx_halo","nens"} );

    auto state   = dm.get<real,5>("dynamics_state");
    auto tracers = dm.get<real,5>("dynamics_tracers");
    parallel_for( Bounds<4>(nz+2*hs,ny+2*hs,nx+2*hs,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      for (int l=0; l < num_state; l++) {
        state(l,k,j,i,iens) = 0;
      }
      for (int tr=0; tr < num_tracers; tr++) {
        tracers(tr,k,j,i,iens) = 0;
      }
    });

    #ifdef PAM_STANDALONE
      std::cout << "nx: " << nx << "\n";
      std::cout << "ny: " << ny << "\n";
      std::cout << "nz: " << nz << "\n";
      std::cout << "xlen (m): " << xlen << "\n";
      std::cout << "ylen (m): " << ylen << "\n";
      std::cout << "zlen (m): " << zlen.createHostCopy()(0) << "\n";
      std::cout << "Simulation time (s): " << sim_time << "\n";
      std::cout << "Vertical interface heights: ";
      auto zint_host = zint.createHostCopy();
      for (int k=0; k < nz+1; k++) {
        std::cout << zint_host(k,0) << "  ";
      }
      std::cout << "\n\n";
    #endif
  }



  // Initialize the state
  template <class MICRO>
  void init_state_and_tracers( DataManager &dm , MICRO const &micro ) {
    Rd    = micro.constants.R_d;
    cp    = micro.constants.cp_d;
    gamma = micro.constants.gamma_d;
    p0    = micro.constants.p0;
    Rv    = micro.constants.R_v;

    real kappa = micro.constants.kappa_d;

    C0 = pow( Rd * pow( p0 , -kappa ) , gamma );

    YAKL_SCOPE( nx                       , this->nx                      );
    YAKL_SCOPE( ny                       , this->ny                      );
    YAKL_SCOPE( nz                       , this->nz                      );
    YAKL_SCOPE( dx                       , this->dx                      );
    YAKL_SCOPE( dy                       , this->dy                      );
    YAKL_SCOPE( dz                       , this->dz                      );
    YAKL_SCOPE( dz_ghost                 , this->dz_ghost                );
    YAKL_SCOPE( gllPts_ord               , this->gllPts_ord              );
    YAKL_SCOPE( gllWts_ord               , this->gllWts_ord              );
    YAKL_SCOPE( gllPts_ngll              , this->gllPts_ngll             );
    YAKL_SCOPE( gllWts_ngll              , this->gllWts_ngll             );
    YAKL_SCOPE( hyDensCells              , this->hyDensCells             );
    YAKL_SCOPE( hyThetaCells             , this->hyThetaCells            );
    YAKL_SCOPE( hyPressureCells          , this->hyPressureCells         );
    YAKL_SCOPE( hyDensThetaCells         , this->hyDensThetaCells        );
    YAKL_SCOPE( hyDensGLL                , this->hyDensGLL               );
    YAKL_SCOPE( hyThetaGLL               , this->hyThetaGLL              );
    YAKL_SCOPE( hyPressureGLL            , this->hyPressureGLL           );
    YAKL_SCOPE( hyDensThetaGLL           , this->hyDensThetaGLL          );
    YAKL_SCOPE( data_spec                , this->data_spec               );
    YAKL_SCOPE( sim2d                    , this->sim2d                   );
    YAKL_SCOPE( xlen                     , this->xlen                    );
    YAKL_SCOPE( ylen                     , this->ylen                    );
    YAKL_SCOPE( Rd                       , this->Rd                      );
    YAKL_SCOPE( Rv                       , this->Rv                      );
    YAKL_SCOPE( cp                       , this->cp                      );
    YAKL_SCOPE( gamma                    , this->gamma                   );
    YAKL_SCOPE( p0                       , this->p0                      );
    YAKL_SCOPE( C0                       , this->C0                      );
    YAKL_SCOPE( vert_interface           , this->vert_interface          );
    YAKL_SCOPE( vert_interface_ghost     , this->vert_interface_ghost    );

    real5d state   = dm.get<real,5>("dynamics_state"  );
    real5d tracers = dm.get<real,5>("dynamics_tracers");

    // If the data_spec is thermal or ..., then initialize the domain with Exner pressure-based hydrostasis
    // This is mostly to make plotting potential temperature perturbation easier for publications
    if (data_spec == DATA_SPEC_THERMAL) {

      // Setup hydrostatic background state
      parallel_for( SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
        // Compute cell averages
        hyDensCells     (k,iens) = 0;
        hyPressureCells (k,iens) = 0;
        hyThetaCells    (k,iens) = 0;
        hyDensThetaCells(k,iens) = 0;
        for (int kk=0; kk<ord; kk++) {
          real zloc = vert_interface(k,iens) + 0.5_fp*dz(k,iens) + gllPts_ord(kk)*dz(k,iens);
          if        (data_spec == DATA_SPEC_THERMAL) {
            // Compute constant theta hydrostatic background state
            real th  = 300;
            real rh = profiles::initConstTheta_density (th,zloc,Rd,cp,gamma,p0,C0);
            real ph = profiles::initConstTheta_pressure(th,zloc,Rd,cp,gamma,p0,C0);
            real wt = gllWts_ord(kk);
            hyDensCells     (k,iens) += rh    * wt;
            hyThetaCells    (k,iens) += th    * wt;
            hyDensThetaCells(k,iens) += rh*th * wt;
            hyPressureCells (k,iens) += ph    * wt;
          }
        }
      });

      parallel_for( SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
        // Compute ngll GLL points
        for (int kk=0; kk<ngll; kk++) {
          real zloc = vert_interface(k,iens) + 0.5_fp*dz(k,iens) + gllPts_ngll(kk)*dz(k,iens);
          if        (data_spec == DATA_SPEC_THERMAL) {
            // Compute constant theta hydrostatic background state
            real th = 300;
            real rh = profiles::initConstTheta_density (th,zloc,Rd,cp,gamma,p0,C0);
            real ph = profiles::initConstTheta_pressure(th,zloc,Rd,cp,gamma,p0,C0);
            hyDensGLL     (k,kk,iens) = rh;
            hyThetaGLL    (k,kk,iens) = th;
            hyDensThetaGLL(k,kk,iens) = rh*th;
            hyPressureGLL (k,kk,iens) = ph;
          }
        }
      });

      // Compute the state
      parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
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
                real rh = profiles::initConstTheta_density(th,zloc,Rd,cp,gamma,p0,C0);
                real tp = profiles::ellipsoid_linear(xloc, yloc, zloc, xlen/2, ylen/2, 2000, 2000, 2000, 2000, 2 );
                real t = th + tp;
                real r = rh;

                state(idR,hs+k,hs+j,hs+i,iens) += (r - rh)*wt;
                state(idT,hs+k,hs+j,hs+i,iens) += (r*t - rh*th) * wt;
              }
            }
          }
        }
      });

      init_tracers( dm , micro );

      adjust_state_for_moisture( dm , micro );

    } // if (data_spec == DATA_SPEC_THERMAL)


    if (data_spec == DATA_SPEC_SUPERCELL) {

      // This uses a piecewise linear profile for Temperature
      real constexpr z_0    = 0;
      real constexpr z_trop = 12000;
      real constexpr T_0    = 300;
      real constexpr T_trop = 213;
      real constexpr T_top  = 213;
      real constexpr p_0    = 100000;

      real4d quad_temp     ("quad_temp"     ,nz,ngll-1,ord,nens);
      real3d hyDensVapGLL  ("hyDensVapGLL"  ,nz,ngll,nens);
      real2d hyDensVapCells("hyDensVapCells",nz,nens);
      real2d z             ("z"             ,nz,nens);
      real2d temp_hy       ("temp_hy"       ,nz,nens);
      real2d tdew_hy       ("tdew_hy"       ,nz,nens);

      YAKL_SCOPE( zlen , this->zlen );

      // Compute full density at ord GLL points for the space between each cell
      parallel_for( Bounds<4>(nz,ngll-1,ord,nens) , YAKL_LAMBDA (int k, int kk, int kkk, int iens) {
        // Middle of this cell
        real cellmid   = vert_interface(k,iens) + 0.5_fp*dz(k,iens);
        // Bottom, top, and middle of the space between these two ngll GLL points
        real ngll_b    = cellmid + gllPts_ngll(kk  )*dz(k,iens);
        real ngll_t    = cellmid + gllPts_ngll(kk+1)*dz(k,iens);
        real ngll_m    = 0.5_fp * (ngll_b + ngll_t);
        // Compute grid spacing between these ngll GLL points
        real ngll_dz   = dz(k,iens) * ( gllPts_ngll(kk+1) - gllPts_ngll(kk) );
        // Compute the locate of this ord GLL point within the ngll GLL points
        real zloc      = ngll_m + ngll_dz * gllPts_ord(kkk);
        // Compute full density at this location
        real temp      = profiles::init_supercell_temperature (zloc, z_0, z_trop, zlen(iens), T_0, T_trop, T_top);
        real press_dry = profiles::init_supercell_pressure_dry(zloc, z_0, z_trop, zlen(iens), T_0, T_trop, T_top, p_0, Rd);
        real dens_dry  = press_dry / (Rd*temp);
        real qvs       = profiles::init_supercell_sat_mix_dry(press_dry, temp);
        real relhum    = profiles::init_supercell_relhum(zloc, z_0, z_trop);
        if (relhum * qvs > 0.014_fp) relhum = 0.014_fp / qvs;
        real qv        = min( 0.014_fp , qvs*relhum );
        quad_temp(k,kk,kkk,iens) = -(1+qv)*GRAV/(Rd+qv*Rv)/temp;
      });

      parallel_for( nens , YAKL_LAMBDA (int iens) {
        hyPressureGLL(0,0,iens) = p_0;
        for (int k=0; k < nz; k++) {
          for (int kk=0; kk < ngll-1; kk++) {
            real tot = 0;
            for (int kkk=0; kkk < ord; kkk++) {
              tot += quad_temp(k,kk,kkk,iens) * gllWts_ord(kkk);
            }
            tot *= dz(k,iens) * ( gllPts_ngll(kk+1) - gllPts_ngll(kk) );
            hyPressureGLL(k,kk+1,iens) = hyPressureGLL(k,kk,iens) * exp( tot );
            if (kk == ngll-2 && k < nz-1) {
              hyPressureGLL(k+1,0,iens) = hyPressureGLL(k,ngll-1,iens);
            }
          }
        }
      });

      parallel_for( Bounds<3>(nz,ngll,nens) , YAKL_LAMBDA (int k, int kk, int iens) {
        real zloc = vert_interface(k,iens) + 0.5_fp*dz(k,iens) + gllPts_ngll(kk)*dz(k,iens);
        real temp       = profiles::init_supercell_temperature (zloc, z_0, z_trop, zlen(iens), T_0, T_trop, T_top);
        real press_tmp  = profiles::init_supercell_pressure_dry(zloc, z_0, z_trop, zlen(iens), T_0, T_trop, T_top, p_0, Rd);
        real qvs        = profiles::init_supercell_sat_mix_dry(press_tmp, temp);
        real relhum     = profiles::init_supercell_relhum(zloc, z_0, z_trop);
        if (relhum * qvs > 0.014_fp) relhum = 0.014_fp / qvs;
        real qv         = min( 0.014_fp , qvs*relhum );
        real press      = hyPressureGLL(k,kk,iens);
        real dens_dry   = press / (Rd+qv*Rv) / temp;
        real dens_vap   = qv * dens_dry;
        real press_vap  = dens_vap * Rv * temp;
        real press_dry  = dens_dry * Rd * temp;
        real dens       = dens_dry + dens_vap;
        real dens_theta = pow( press / C0 , 1._fp / gamma );
        real theta      = dens_theta / dens;
        hyDensGLL     (k,kk,iens) = dens;
        hyDensThetaGLL(k,kk,iens) = dens_theta;
        hyThetaGLL    (k,kk,iens) = theta;
        hyDensVapGLL  (k,kk,iens) = dens_vap;
      });

      parallel_for( Bounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
        real press_tot      = 0;
        real dens_tot       = 0;
        real dens_vap_tot   = 0;
        real theta_tot      = 0;
        real dens_theta_tot = 0;
        for (int kk=0; kk < ngll; kk++) {
          press_tot      += hyPressureGLL (k,kk,iens) * gllWts_ngll(kk);
          dens_tot       += hyDensGLL     (k,kk,iens) * gllWts_ngll(kk);
          dens_vap_tot   += hyDensVapGLL  (k,kk,iens) * gllWts_ngll(kk);
          dens_theta_tot += hyDensThetaGLL(k,kk,iens) * gllWts_ngll(kk);
          theta_tot      += hyThetaGLL    (k,kk,iens) * gllWts_ngll(kk);
        }
        real press      = press_tot;
        real dens       = dens_tot;
        real dens_vap   = dens_vap_tot;
        real dens_theta = dens_theta_tot;
        real theta      = theta_tot;
        real dens_dry   = dens - dens_vap;
        real R          = dens_dry / dens * Rd + dens_vap / dens * Rv;
        real temp       = press / (dens * R);
        real qv         = dens_vap / dens_dry;
        real zloc       = vert_interface(k,iens) + 0.5_fp*dz(k,iens);
        real press_tmp  = profiles::init_supercell_pressure_dry(zloc, z_0, z_trop, zlen(iens), T_0, T_trop, T_top, p_0, Rd);
        real qvs        = profiles::init_supercell_sat_mix_dry(press_tmp, temp);
        real relhum     = qv / qvs;
        real T          = temp - 273;
        real a          = 17.27;
        real b          = 237.7;
        real tdew       = b * ( a*T / (b + T) + log(relhum) ) / ( a - ( a*T / (b+T) + log(relhum) ) );
        // The next three are just to confirm the skew-T diagram looks OK
        z               (k,iens) = zloc;
        temp_hy         (k,iens) = temp;
        tdew_hy         (k,iens) = tdew;
        // These are used in the rest of the model
        hyPressureCells (k,iens) = press;
        hyDensCells     (k,iens) = dens;
        hyDensThetaCells(k,iens) = dens_theta;
        hyThetaCells    (k,iens) = theta;
        hyDensVapCells  (k,iens) = dens_vap;
      });

      // Dump out data to plot a skew-T log-P diagram
      yakl::SimpleNetCDF nc;
      nc.create("skew.nc");
      real1d data("data",nz);
      for (int iens=0; iens < nens; iens++) {
        parallel_for( nz , YAKL_LAMBDA (int k) { data(k) = z(k,iens); });
        nc.write(data.createHostCopy(),"z"          ,{"z"});

        parallel_for( nz , YAKL_LAMBDA (int k) { data(k) = hyPressureCells(k,iens); });
        nc.write(data.createHostCopy(),"pressure"   ,{"z"});

        parallel_for( nz , YAKL_LAMBDA (int k) { data(k) = temp_hy(k,iens); });
        nc.write(data.createHostCopy(),"temperature",{"z"});

        parallel_for( nz , YAKL_LAMBDA (int k) { data(k) = tdew_hy(k,iens); });
        nc.write(data.createHostCopy(),"dew_point"  ,{"z"});
      }
      nc.close();

      int idWV = micro.get_water_vapor_index();
      real5d tracers = dm.get<real,5>("dynamics_tracers");

      parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        state  (idR ,hs+k,hs+j,hs+i,iens) = 0;
        state  (idU ,hs+k,hs+j,hs+i,iens) = 0;
        state  (idV ,hs+k,hs+j,hs+i,iens) = 0;
        state  (idW ,hs+k,hs+j,hs+i,iens) = 0;
        state  (idT ,hs+k,hs+j,hs+i,iens) = 0;
        tracers(idWV,hs+k,hs+j,hs+i,iens) = 0;
        for (int kk=0; kk < ngll; kk++) {
          for (int jj=0; jj < ngll; jj++) {
            for (int ii=0; ii < ngll; ii++) {
              real xloc = (i+0.5_fp)*dx                              + gllPts_ngll(ii)*dx;
              real yloc = (j+0.5_fp)*dy                              + gllPts_ngll(jj)*dy;
              real zloc = vert_interface(k,iens) + 0.5_fp*dz(k,iens) + gllPts_ngll(kk)*dz(k,iens);

              if (sim2d) yloc = ylen/2;

              real dens = hyDensGLL(k,kk,iens);

              real uvel;
              real constexpr zs = 5000;
              real constexpr us = 30;
              real constexpr uc = 15;
              if (zloc < zs) {
                uvel = us * (zloc / zs) - uc;
              } else {
                uvel = us - uc;
              }

              real vvel       = 0;
              real wvel       = 0;
              real theta      = hyThetaGLL    (k,kk,iens);
              real dens_vap   = hyDensVapGLL  (k,kk,iens);
              real dens_theta = hyDensThetaGLL(k,kk,iens);

              real x0 = xlen / 2;
              real y0 = ylen / 2;
              real z0 = 1500;
              real radx = 10000;
              real rady = 10000;
              real radz = 1500;
              real amp  = 3;

              real xn = (xloc - x0) / radx;
              real yn = (yloc - y0) / rady;
              real zn = (zloc - z0) / radz;

              real rad = sqrt( xn*xn + yn*yn + zn*zn );

              real theta_pert = 0;
              if (rad < 1) {
                theta_pert = amp * pow( cos(M_PI*rad/2) , 2._fp );
              }

              dens_theta += dens * theta_pert;

              real factor = gllWts_ngll(ii) * gllWts_ngll(jj) * gllWts_ngll(kk);
              state  (idR ,hs+k,hs+j,hs+i,iens) += (dens - hyDensGLL(k,kk,iens))            * factor;
              state  (idU ,hs+k,hs+j,hs+i,iens) += dens * uvel                              * factor;
              state  (idV ,hs+k,hs+j,hs+i,iens) += dens * vvel                              * factor;
              state  (idW ,hs+k,hs+j,hs+i,iens) += dens * wvel                              * factor;
              state  (idT ,hs+k,hs+j,hs+i,iens) += (dens_theta - hyDensThetaGLL(k,kk,iens)) * factor;
              tracers(idWV,hs+k,hs+j,hs+i,iens) += dens_vap                                 * factor;
            }
          }
        }
      });

    } // if (data_spec == DATA_SPEC_SUPERCELL)

    convert_dynamics_to_coupler_state( dm , micro );
  }



  // Compute state and tendency time derivatives from the state
  template <class MICRO>
  void computeTendencies( real5d &state   , real5d &stateTend  ,
                          real5d &tracers , real5d &tracerTend ,
                          MICRO const &micro, real &dt , int splitIndex ) {
    YAKL_SCOPE( nx                      , this->nx                     );
    YAKL_SCOPE( ny                      , this->ny                     );
    YAKL_SCOPE( nz                      , this->nz                     );
    YAKL_SCOPE( weno_scalars            , this->weno_scalars           );
    YAKL_SCOPE( weno_winds              , this->weno_winds             );
    YAKL_SCOPE( c2g                     , this->coefs_to_gll           );
    YAKL_SCOPE( s2g                     , this->sten_to_gll            );
    YAKL_SCOPE( wenoRecon               , this->wenoRecon              );
    YAKL_SCOPE( idl                     , this->idl                    );
    YAKL_SCOPE( sigma                   , this->sigma                  );
    YAKL_SCOPE( hyDensCells             , this->hyDensCells            );
    YAKL_SCOPE( hyDensThetaCells        , this->hyDensThetaCells       );
    YAKL_SCOPE( hyThetaCells            , this->hyThetaCells           );
    YAKL_SCOPE( hyPressureCells         , this->hyPressureCells        );
    YAKL_SCOPE( hyDensGLL               , this->hyDensGLL              );
    YAKL_SCOPE( hyDensThetaGLL          , this->hyDensThetaGLL         );
    YAKL_SCOPE( hyThetaGLL              , this->hyThetaGLL             );
    YAKL_SCOPE( hyPressureGLL           , this->hyPressureGLL          );
    YAKL_SCOPE( sim2d                   , this->sim2d                  );
    YAKL_SCOPE( derivMatrix             , this->derivMatrix            );
    YAKL_SCOPE( dx                      , this->dx                     );
    YAKL_SCOPE( dy                      , this->dy                     );
    YAKL_SCOPE( dz                      , this->dz                     );
    YAKL_SCOPE( stateLimits_x           , this->stateLimits_x          );
    YAKL_SCOPE( tracerLimits_x          , this->tracerLimits_x         );
    YAKL_SCOPE( stateFlux_x             , this->stateFlux_x            );
    YAKL_SCOPE( tracerFlux_x            , this->tracerFlux_x           );
    YAKL_SCOPE( stateLimits_y           , this->stateLimits_y          );
    YAKL_SCOPE( tracerLimits_y          , this->tracerLimits_y         );
    YAKL_SCOPE( stateFlux_y             , this->stateFlux_y            );
    YAKL_SCOPE( tracerFlux_y            , this->tracerFlux_y           );
    YAKL_SCOPE( stateLimits_z           , this->stateLimits_z          );
    YAKL_SCOPE( tracerLimits_z          , this->tracerLimits_z         );
    YAKL_SCOPE( stateFlux_z             , this->stateFlux_z            );
    YAKL_SCOPE( tracerFlux_z            , this->tracerFlux_z           );
    YAKL_SCOPE( tracer_pos              , this->tracer_pos             );
    YAKL_SCOPE( num_tracers             , this->num_tracers            );
    YAKL_SCOPE( bc_x                    , this->bc_x                   );
    YAKL_SCOPE( bc_y                    , this->bc_y                   );
    YAKL_SCOPE( bc_z                    , this->bc_z                   );
    YAKL_SCOPE( Rd                      , this->Rd                     );
    YAKL_SCOPE( cp                      , this->cp                     );
    YAKL_SCOPE( gamma                   , this->gamma                  );
    YAKL_SCOPE( p0                      , this->p0                     );
    YAKL_SCOPE( C0                      , this->C0                     );
    YAKL_SCOPE( gllWts_ngll             , this->gllWts_ngll            );

    
    // Pre-process the tracers by dividing by density inside the domain
    // After this, we can reconstruct tracers only (not rho * tracer)
    parallel_for( SimpleBounds<5>(num_tracers,nz,ny,nx,nens) , YAKL_LAMBDA (int tr, int k, int j, int i, int iens) {
      tracers(tr,hs+k,hs+j,hs+i,iens) /= (state(idR,hs+k,hs+j,hs+i,iens) + hyDensCells(k,iens));
    });

    // Populate the x halos
    if        (bc_x == BC_PERIODIC) {
      parallel_for( SimpleBounds<4>(nz,ny,hs,nens) , YAKL_LAMBDA(int k, int j, int ii, int iens) {
        for (int l=0; l < num_state; l++) {
          state  (l,hs+k,hs+j,      ii,iens) = state  (l,hs+k,hs+j,nx+ii,iens);
          state  (l,hs+k,hs+j,hs+nx+ii,iens) = state  (l,hs+k,hs+j,hs+ii,iens);
        }
        for (int l=0; l < num_tracers; l++) {
          tracers(l,hs+k,hs+j,      ii,iens) = tracers(l,hs+k,hs+j,nx+ii,iens);
          tracers(l,hs+k,hs+j,hs+nx+ii,iens) = tracers(l,hs+k,hs+j,hs+ii,iens);
        }
      });
    } else if (bc_x == BC_WALL) {
      parallel_for( SimpleBounds<4>(nz,ny,hs,nens) , YAKL_LAMBDA(int k, int j, int ii, int iens) {
        for (int l=0; l < num_state; l++) {
          if (l == idU) {
            state(l,hs+k,hs+j,      ii,iens) = 0;
            state(l,hs+k,hs+j,hs+nx+ii,iens) = 0;
          } else {
            state(l,hs+k,hs+j,      ii,iens) = state(l,hs+k,hs+j,hs     ,iens);
            state(l,hs+k,hs+j,hs+nx+ii,iens) = state(l,hs+k,hs+j,hs+nx-1,iens);
          }
        }
        for (int l=0; l < num_tracers; l++) {
          tracers(l,hs+k,hs+j,      ii,iens) = tracers(l,hs+k,hs+j,hs     ,iens);
          tracers(l,hs+k,hs+j,hs+nx+ii,iens) = tracers(l,hs+k,hs+j,hs+nx-1,iens);
        }
      });
    }

    // Populate the y halos
    if (!sim2d) {
      if        (bc_y == BC_PERIODIC) {
        parallel_for( SimpleBounds<4>(nz,nx,hs,nens) , YAKL_LAMBDA(int k, int i, int jj, int iens) {
          for (int l=0; l < num_state; l++) {
            state(l,hs+k,      jj,hs+i,iens) = state(l,hs+k,ny+jj,hs+i,iens);
            state(l,hs+k,hs+ny+jj,hs+i,iens) = state(l,hs+k,hs+jj,hs+i,iens);
          }
          for (int l=0; l < num_tracers; l++) {
            tracers(l,hs+k,      jj,hs+i,iens) = tracers(l,hs+k,ny+jj,hs+i,iens);
            tracers(l,hs+k,hs+ny+jj,hs+i,iens) = tracers(l,hs+k,hs+jj,hs+i,iens);
          }
        });
      } else if (bc_y == BC_WALL) {
        parallel_for( SimpleBounds<4>(nz,nx,hs,nens) , YAKL_LAMBDA(int k, int i, int jj, int iens) {
          for (int l=0; l < num_state; l++) {
            if (l == idV) {
              state(l,hs+k,      jj,hs+i,iens) = 0;
              state(l,hs+k,hs+ny+jj,hs+i,iens) = 0;
            } else {
              state(l,hs+k,      jj,hs+i,iens) = state(l,hs+k,hs     ,hs+i,iens);
              state(l,hs+k,hs+ny+jj,hs+i,iens) = state(l,hs+k,hs+ny-1,hs+i,iens);
            }
          }
          for (int l=0; l < num_tracers; l++) {
            tracers(l,hs+k,      jj,hs+i,iens) = tracers(l,hs+k,hs     ,hs+i,iens);
            tracers(l,hs+k,hs+ny+jj,hs+i,iens) = tracers(l,hs+k,hs+ny-1,hs+i,iens);
          }
        });
      }
    }

    // Populate the z halos
    if        (bc_z == BC_PERIODIC) {
      parallel_for( SimpleBounds<4>(ny,nx,hs,nens) , YAKL_LAMBDA(int j, int i, int kk, int iens) {
        for (int l=0; l < num_state; l++) {
          state(l,      kk,hs+j,hs+i,iens) = state(l,nz+kk,hs+j,hs+i,iens);
          state(l,hs+nz+kk,hs+j,hs+i,iens) = state(l,hs+kk,hs+j,hs+i,iens);
        }
        for (int l=0; l < num_tracers; l++) {
          tracers(l,      kk,hs+j,hs+i,iens) = tracers(l,nz+kk,hs+j,hs+i,iens);
          tracers(l,hs+nz+kk,hs+j,hs+i,iens) = tracers(l,hs+kk,hs+j,hs+i,iens);
        }
      });
    } else if (bc_z == BC_WALL) {
      parallel_for( SimpleBounds<4>(ny,nx,hs,nens) , YAKL_LAMBDA(int j, int i, int kk, int iens) {
        for (int l=0; l < num_state; l++) {
          if (l == idW) {
            state(l,      kk,hs+j,hs+i,iens) = 0;
            state(l,hs+nz+kk,hs+j,hs+i,iens) = 0;
          } else {
            state(l,      kk,hs+j,hs+i,iens) = state(l,hs     ,hs+j,hs+i,iens);
            state(l,hs+nz+kk,hs+j,hs+i,iens) = state(l,hs+nz-1,hs+j,hs+i,iens);
          }
        }
        for (int l=0; l < num_tracers; l++) {
          tracers(l,      kk,hs+j,hs+i,iens) = tracers(l,hs     ,hs+j,hs+i,iens);
          tracers(l,hs+nz+kk,hs+j,hs+i,iens) = tracers(l,hs+nz-1,hs+j,hs+i,iens);
        }
      });
    }



    // Loop through cells, reconstruct, store interface limit values for state and tracers
    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      { // x-direction
        // We need density and momentum to evolve the tracers with ADER
        SArray<real,1,ngll> r_gll;

        { // State
          SArray<real,1,ngll> ru_gll, rv_gll , rw_gll , rt_gll;
          { // Reconstruct
            SArray<real,1,ord> stencil;

            // Density
            for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idR,hs+k,hs+j,i+ii,iens); }
            reconstruct_gll_values( stencil , r_gll , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
            for (int ii=0; ii < ngll; ii++) { r_gll(ii) += hyDensCells(k,iens); } // Add hydrostasis back on

            // u values and derivatives
            for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idU,hs+k,hs+j,i+ii,iens); }
            reconstruct_gll_values( stencil , ru_gll , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

            // v
            for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idV,hs+k,hs+j,i+ii,iens); }
            reconstruct_gll_values( stencil , rv_gll , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

            // w
            for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idW,hs+k,hs+j,i+ii,iens); }
            reconstruct_gll_values( stencil , rw_gll , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

            // theta
            for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idT,hs+k,hs+j,i+ii,iens); }
            reconstruct_gll_values( stencil , rt_gll , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
            for (int ii=0; ii < ngll; ii++) { rt_gll(ii) += hyDensThetaCells(k,iens); } // Add hydrostasis back on
          } // Reconstruct

          // Left interface
          stateLimits_x(idR,1,k,j,i  ,iens) = r_gll (0     );
          stateLimits_x(idU,1,k,j,i  ,iens) = ru_gll(0     );
          stateLimits_x(idV,1,k,j,i  ,iens) = rv_gll(0     );
          stateLimits_x(idW,1,k,j,i  ,iens) = rw_gll(0     );
          stateLimits_x(idT,1,k,j,i  ,iens) = rt_gll(0     );
          // Right interface
          stateLimits_x(idR,0,k,j,i+1,iens) = r_gll (ngll-1);
          stateLimits_x(idU,0,k,j,i+1,iens) = ru_gll(ngll-1);
          stateLimits_x(idV,0,k,j,i+1,iens) = rv_gll(ngll-1);
          stateLimits_x(idW,0,k,j,i+1,iens) = rw_gll(ngll-1);
          stateLimits_x(idT,0,k,j,i+1,iens) = rt_gll(ngll-1);
        } // state

        { // tracers
          // Only process one tracer at a time to save on local memory / register requirements
          for (int tr=0; tr < num_tracers; tr++) {
            SArray<real,1,ngll> t_gll; // tracer (not density * tracer)
            { // Reconstruct
              SArray<real,1,ord> stencil;
              for (int ii=0; ii < ord; ii++) { stencil(ii) = tracers(tr,hs+k,hs+j,i+ii,iens); }
              reconstruct_gll_values( stencil , t_gll , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
              if (tracer_pos(tr)) {
                for (int ii=0; ii < ngll; ii++) { t_gll(ii) = max( 0._fp , t_gll(ii) ); }
              }
            } // Reconstruct

            tracerLimits_x(tr,1,k,j,i  ,iens) = r_gll(0     ) * t_gll(0     ); // Left interface
            tracerLimits_x(tr,0,k,j,i+1,iens) = r_gll(ngll-1) * t_gll(ngll-1); // Right interface
          }
        }  // tracers
      } // x-direction


      if (!sim2d) { // y-direction
        // These are needed by the tracers
        SArray<real,1,ngll> r_gll;

        { // State
          SArray<real,1,ngll> ru_gll, rv_gll, rw_gll , rt_gll;
          { // Reconstruct
            SArray<real,1,ord> stencil;

            // Density
            for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idR,hs+k,j+jj,hs+i,iens); }
            reconstruct_gll_values( stencil , r_gll , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
            for (int jj=0; jj < ngll; jj++) { r_gll(jj) += hyDensCells(k,iens); } // Add hydrostasis back on

            // u values and derivatives
            for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idU,hs+k,j+jj,hs+i,iens); }
            reconstruct_gll_values( stencil , ru_gll , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

            // v
            for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idV,hs+k,j+jj,hs+i,iens); }
            reconstruct_gll_values( stencil , rv_gll , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

            // w
            for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idW,hs+k,j+jj,hs+i,iens); }
            reconstruct_gll_values( stencil , rw_gll , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

            // theta
            for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idT,hs+k,j+jj,hs+i,iens); }
            reconstruct_gll_values( stencil , rt_gll , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
            for (int jj=0; jj < ngll; jj++) { rt_gll(jj) += hyDensThetaCells(k,iens); } // Add hydrostasis back on
          } // Reconstruct

          // Left interface
          stateLimits_y(idR,1,k,j  ,i,iens) = r_gll (0     );
          stateLimits_y(idU,1,k,j  ,i,iens) = ru_gll(0     );
          stateLimits_y(idV,1,k,j  ,i,iens) = rv_gll(0     );
          stateLimits_y(idW,1,k,j  ,i,iens) = rw_gll(0     );
          stateLimits_y(idT,1,k,j  ,i,iens) = rt_gll(0     );
          // Right interface
          stateLimits_y(idR,0,k,j+1,i,iens) = r_gll (ngll-1);
          stateLimits_y(idU,0,k,j+1,i,iens) = ru_gll(ngll-1);
          stateLimits_y(idV,0,k,j+1,i,iens) = rv_gll(ngll-1);
          stateLimits_y(idW,0,k,j+1,i,iens) = rw_gll(ngll-1);
          stateLimits_y(idT,0,k,j+1,i,iens) = rt_gll(ngll-1);
        } // State

        { // Tracers
          // Only process one tracer at a time to save on local memory / register requirements
          for (int tr=0; tr < num_tracers; tr++) {
            SArray<real,1,ngll> t_gll; // Density * tracer
            { // Reconstruct
              SArray<real,1,ord> stencil;
              for (int jj=0; jj < ord; jj++) { stencil(jj) = tracers(tr,hs+k,j+jj,hs+i,iens); }
              reconstruct_gll_values( stencil , t_gll , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
              if (tracer_pos(tr)) {
                for (int jj=0; jj < ngll; jj++) { t_gll(jj) = max( 0._fp , t_gll(jj) ); }
              }
            } // Reconstruct

            tracerLimits_y(tr,1,k,j  ,i,iens) = r_gll(0     ) * t_gll(0     ); // Left interface
            tracerLimits_y(tr,0,k,j+1,i,iens) = r_gll(ngll-1) * t_gll(ngll-1); // Right interface
          }
        } // Tracers
      } // y-direction

      { // z-direction
        // We need these to persist to evolve tracers with ADER
        SArray<real,1,ngll> r_gll;

        { // State
          SArray<real,1,ngll> ru_gll , rv_gll , rw_gll , rt_gll;
          { // Reconstruct
            SArray<real,1,ord> stencil;

            // Density
            for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idR,k+kk,hs+j,hs+i,iens); }
            reconstruct_gll_values( stencil , r_gll , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
            for (int kk=0; kk < ngll; kk++) { r_gll(kk) += hyDensGLL(k,kk,iens); } // Add hydrostasis back on

            // u values and derivatives
            for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idU,k+kk,hs+j,hs+i,iens); }
            reconstruct_gll_values( stencil , ru_gll , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

            // v
            for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idV,k+kk,hs+j,hs+i,iens); }
            reconstruct_gll_values( stencil , rv_gll , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

            // w
            for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idW,k+kk,hs+j,hs+i,iens); }
            reconstruct_gll_values( stencil , rw_gll , c2g , s2g , wenoRecon , idl , sigma , weno_winds );
            if (bc_z == BC_WALL) {
              if (k == nz-1) rw_gll(ngll-1) = 0;
              if (k == 0   ) rw_gll(0     ) = 0;
            }

            // theta
            for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idT,k+kk,hs+j,hs+i,iens); }
            reconstruct_gll_values( stencil , rt_gll , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
            for (int kk=0; kk < ngll; kk++) { rt_gll(kk) += hyDensThetaGLL(k,kk,iens); } // Add hydrostasis back on
          } // Reconstruct

          // Left interface
          stateLimits_z(idR,1,k  ,j,i,iens) = r_gll (0     );
          stateLimits_z(idU,1,k  ,j,i,iens) = ru_gll(0     );
          stateLimits_z(idV,1,k  ,j,i,iens) = rv_gll(0     );
          stateLimits_z(idW,1,k  ,j,i,iens) = rw_gll(0     );
          stateLimits_z(idT,1,k  ,j,i,iens) = rt_gll(0     );
          // Right interface
          stateLimits_z(idR,0,k+1,j,i,iens) = r_gll (ngll-1);
          stateLimits_z(idU,0,k+1,j,i,iens) = ru_gll(ngll-1);
          stateLimits_z(idV,0,k+1,j,i,iens) = rv_gll(ngll-1);
          stateLimits_z(idW,0,k+1,j,i,iens) = rw_gll(ngll-1);
          stateLimits_z(idT,0,k+1,j,i,iens) = rt_gll(ngll-1);

          ////////////////////////////////////////////
          // Assign gravity source term
          ////////////////////////////////////////////
          real ravg = 0;
          for (int kk=0; kk < ngll; kk++) {
            ravg += r_gll(kk) * gllWts_ngll(kk);
          }
          stateTend(idR,k,j,i,iens) = 0;
          stateTend(idU,k,j,i,iens) = 0;
          stateTend(idV,k,j,i,iens) = 0;
          stateTend(idW,k,j,i,iens) = -GRAV*ravg;
          stateTend(idT,k,j,i,iens) = 0;
        } // State

        { // Tracers
          // Only process one tracer at a time to save on local memory / register requirements
          for (int tr=0; tr < num_tracers; tr++) {
            SArray<real,1,ngll> t_gll; // Density * tracer
            { // Reconstruct
              SArray<real,1,ord> stencil;
              for (int kk=0; kk < ord; kk++) { stencil(kk) = tracers(tr,k+kk,hs+j,hs+i,iens); }
              reconstruct_gll_values( stencil , t_gll , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
              if (tracer_pos(tr)) {
                for (int kk=0; kk < ngll; kk++) { t_gll(kk) = max( 0._fp , t_gll(kk) ); }
              }
            } // Reconstruct

            tracerLimits_z(tr,1,k  ,j,i,iens) = r_gll(0     ) * t_gll(0     ); // Left interface
            tracerLimits_z(tr,0,k+1,j,i,iens) = r_gll(ngll-1) * t_gll(ngll-1); // Right interface
          }
        } // Tracers
      } // z-direction

    });



    ////////////////////////////////////////////////
    // BCs for the state edge estimates (x direction)
    ////////////////////////////////////////////////
    parallel_for( SimpleBounds<3>(nz,ny,nens) , YAKL_LAMBDA (int k, int j, int iens) {
      for (int l=0; l < num_state; l++) {
        if        (bc_x == BC_PERIODIC) {
          stateLimits_x(l,0,k,j,0 ,iens) = stateLimits_x(l,0,k,j,nx,iens);
          stateLimits_x(l,1,k,j,nx,iens) = stateLimits_x(l,1,k,j,0 ,iens);
        } else if (bc_x == BC_WALL    ) {
          stateLimits_x(l,0,k,j,0 ,iens) = stateLimits_x(l,1,k,j,0 ,iens);
          stateLimits_x(l,1,k,j,nx,iens) = stateLimits_x(l,0,k,j,nx,iens);
        }
      }
      for (int l=0; l < num_tracers; l++) {
        if        (bc_x == BC_PERIODIC) {
          tracerLimits_x(l,0,k,j,0 ,iens) = tracerLimits_x(l,0,k,j,nx,iens);
          tracerLimits_x(l,1,k,j,nx,iens) = tracerLimits_x(l,1,k,j,0 ,iens);
        } else if (bc_x == BC_WALL    ) {
          tracerLimits_x(l,0,k,j,0 ,iens) = tracerLimits_x(l,1,k,j,0 ,iens);
          tracerLimits_x(l,1,k,j,nx,iens) = tracerLimits_x(l,0,k,j,nx,iens);
        }
      }
    });
    //////////////////////////////////////////////////////////
    // Compute the upwind fluxes (x-direction)
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<4>(nz,ny,nx+1,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      // Get left and right state
      real r_L = stateLimits_x(idR,0,k,j,i,iens)    ;   real r_R = stateLimits_x(idR,1,k,j,i,iens)    ;
      real u_L = stateLimits_x(idU,0,k,j,i,iens)/r_L;   real u_R = stateLimits_x(idU,1,k,j,i,iens)/r_R;
      real v_L = stateLimits_x(idV,0,k,j,i,iens)/r_L;   real v_R = stateLimits_x(idV,1,k,j,i,iens)/r_R;
      real w_L = stateLimits_x(idW,0,k,j,i,iens)/r_L;   real w_R = stateLimits_x(idW,1,k,j,i,iens)/r_R;
      real t_L = stateLimits_x(idT,0,k,j,i,iens)/r_L;   real t_R = stateLimits_x(idT,1,k,j,i,iens)/r_R;
      // Compute average state
      real r = 0.5_fp * (r_L + r_R);
      real u = 0.5_fp * (u_L + u_R);
      real v = 0.5_fp * (v_L + v_R);
      real w = 0.5_fp * (w_L + w_R);
      real t = 0.5_fp * (t_L + t_R);
      real p = C0 * pow(r*t,gamma);
      real cs2 = gamma*p/r;
      real cs  = sqrt(cs2);

      // COMPUTE UPWIND STATE FLUXES
      // Get left and right fluxes
      real q1_L = stateLimits_x(idR,0,k,j,i,iens);   real q1_R = stateLimits_x(idR,1,k,j,i,iens);
      real q2_L = stateLimits_x(idU,0,k,j,i,iens);   real q2_R = stateLimits_x(idU,1,k,j,i,iens);
      real q3_L = stateLimits_x(idV,0,k,j,i,iens);   real q3_R = stateLimits_x(idV,1,k,j,i,iens);
      real q4_L = stateLimits_x(idW,0,k,j,i,iens);   real q4_R = stateLimits_x(idW,1,k,j,i,iens);
      real q5_L = stateLimits_x(idT,0,k,j,i,iens);   real q5_R = stateLimits_x(idT,1,k,j,i,iens);
      // Compute upwind characteristics
      // Waves 1-3, velocity: u
      real w1, w2, w3;
      if (u > 0) {
        w1 = q1_L - q5_L/t;
        w2 = q3_L - v*q5_L/t;
        w3 = q4_L - w*q5_L/t;
      } else {
        w1 = q1_R - q5_R/t;
        w2 = q3_R - v*q5_R/t;
        w3 = q4_R - w*q5_R/t;
      }
      // Wave 5, velocity: u-cs
      real w5 =  u*q1_R/(2*cs) - q2_R/(2*cs) + q5_R/(2*t);
      // Wave 6, velocity: u+cs
      real w6 = -u*q1_L/(2*cs) + q2_L/(2*cs) + q5_L/(2*t);
      // Use right eigenmatrix to compute upwind flux
      real q1 = w1 + w5 + w6;
      real q2 = u*w1 + (u-cs)*w5 + (u+cs)*w6;
      real q3 = w2 + v*w5 + v*w6;
      real q4 = w3 + w*w5 + w*w6;
      real q5 =      t*w5 + t*w6;

      stateFlux_x(idR,k,j,i,iens) = q2;
      stateFlux_x(idU,k,j,i,iens) = q2*q2/q1 + C0*pow(q5,gamma);
      stateFlux_x(idV,k,j,i,iens) = q2*q3/q1;
      stateFlux_x(idW,k,j,i,iens) = q2*q4/q1;
      stateFlux_x(idT,k,j,i,iens) = q2*q5/q1;

      // COMPUTE UPWIND TRACER FLUXES
      // Handle it one tracer at a time
      for (int tr=0; tr < num_tracers; tr++) {
        if (u > 0) {
          tracerFlux_x(tr,k,j,i,iens) = q2 * tracerLimits_x(tr,0,k,j,i,iens) / r_L;
        } else {
          tracerFlux_x(tr,k,j,i,iens) = q2 * tracerLimits_x(tr,1,k,j,i,iens) / r_R;
        }
      }
    });
    //////////////////////////////////////////////////////////
    // Limit the tracer fluxes for positivity
    //////////////////////////////////////////////////////////
    real5d fct_mult_x("fct_mult_x",num_tracers,nz,ny,nx+1,nens);
    parallel_for( SimpleBounds<5>(num_tracers,nz,ny,nx+1,nens) , YAKL_LAMBDA (int tr, int k, int j, int i, int iens) {
      fct_mult_x(tr,k,j,i,iens) = 1.;
      // Solid wall BCs mean u == 0 at boundaries, so we assume periodic if u != 0
      if (tracer_pos(tr)) {
        // Compute and apply the flux reduction factor of the upwind cell
        if      (tracerFlux_x(tr,k,j,i,iens) > 0) {
          // if u > 0, then it pulls mass out of the left cell
          int ind_i = i-1;
          // TODO: Relax the periodic assumption here
          if (ind_i == -1) ind_i = nx-1;
          real f1 = min( tracerFlux_x(tr,k,j,ind_i  ,iens) , 0._fp );
          real f2 = max( tracerFlux_x(tr,k,j,ind_i+1,iens) , 0._fp );
          real fluxOut = dt*(f2-f1)/dx;
          real mass = tracers(tr,hs+k,hs+j,hs+ind_i,iens);
          if (fluxOut > 0) {
            fct_mult_x(tr,k,j,i,iens) = min( 1._fp , mass / fluxOut );
          }
        } else if (tracerFlux_x(tr,k,j,i,iens) < 0) {
          // upwind is to the right of this interface
          int ind_i = i;
          // TODO: Relax the periodic assumption here
          if (ind_i == nx) ind_i = 0;
          real f1 = min( tracerFlux_x(tr,k,j,ind_i  ,iens) , 0._fp );
          real f2 = max( tracerFlux_x(tr,k,j,ind_i+1,iens) , 0._fp );
          real fluxOut = dt*(f2-f1)/dx;
          real mass = tracers(tr,hs+k,hs+j,hs+ind_i,iens);
          if (fluxOut > 0) {
            fct_mult_x(tr,k,j,i,iens) = min( 1._fp , mass / fluxOut );
          }
        }
      }
    });



    ////////////////////////////////////////////////
    // BCs for the state edge estimates (y-direction)
    ////////////////////////////////////////////////
    parallel_for( SimpleBounds<3>(nz,nx,nens) , YAKL_LAMBDA (int k, int i, int iens) {
      for (int l=0; l < num_state; l++) {
        if        (bc_y == BC_PERIODIC) {
          stateLimits_y(l,0,k,0 ,i,iens) = stateLimits_y(l,0,k,ny,i,iens);
          stateLimits_y(l,1,k,ny,i,iens) = stateLimits_y(l,1,k,0 ,i,iens);
        } else if (bc_y == BC_WALL    ) {
          stateLimits_y(l,0,k,0 ,i,iens) = stateLimits_y(l,1,k,0 ,i,iens);
          stateLimits_y(l,1,k,ny,i,iens) = stateLimits_y(l,0,k,ny,i,iens);
        }
      }
      for (int l=0; l < num_tracers; l++) {
        if        (bc_y == BC_PERIODIC) {
          tracerLimits_y(l,0,k,0 ,i,iens) = tracerLimits_y(l,0,k,ny,i,iens);
          tracerLimits_y(l,1,k,ny,i,iens) = tracerLimits_y(l,1,k,0 ,i,iens);
        } else if (bc_y == BC_WALL    ) {
          tracerLimits_y(l,0,k,0 ,i,iens) = tracerLimits_y(l,1,k,0 ,i,iens);
          tracerLimits_y(l,1,k,ny,i,iens) = tracerLimits_y(l,0,k,ny,i,iens);
        }
      }
    });
    //////////////////////////////////////////////////////////
    // Compute the upwind fluxes (y-direction)
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<4>(nz,ny+1,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      // Get left and right state
      real r_L = stateLimits_y(idR,0,k,j,i,iens)    ;   real r_R = stateLimits_y(idR,1,k,j,i,iens)    ;
      real u_L = stateLimits_y(idU,0,k,j,i,iens)/r_L;   real u_R = stateLimits_y(idU,1,k,j,i,iens)/r_R;
      real v_L = stateLimits_y(idV,0,k,j,i,iens)/r_L;   real v_R = stateLimits_y(idV,1,k,j,i,iens)/r_R;
      real w_L = stateLimits_y(idW,0,k,j,i,iens)/r_L;   real w_R = stateLimits_y(idW,1,k,j,i,iens)/r_R;
      real t_L = stateLimits_y(idT,0,k,j,i,iens)/r_L;   real t_R = stateLimits_y(idT,1,k,j,i,iens)/r_R;
      // Compute average state
      real r = 0.5_fp * (r_L + r_R);
      real u = 0.5_fp * (u_L + u_R);
      real v = 0.5_fp * (v_L + v_R);
      real w = 0.5_fp * (w_L + w_R);
      real t = 0.5_fp * (t_L + t_R);
      real p = C0 * pow(r*t,gamma);
      real cs2 = gamma*p/r;
      real cs  = sqrt(cs2);

      // COMPUTE UPWIND STATE FLUXES
      // Get left and right fluxes
      real q1_L = stateLimits_y(idR,0,k,j,i,iens);   real q1_R = stateLimits_y(idR,1,k,j,i,iens);
      real q2_L = stateLimits_y(idU,0,k,j,i,iens);   real q2_R = stateLimits_y(idU,1,k,j,i,iens);
      real q3_L = stateLimits_y(idV,0,k,j,i,iens);   real q3_R = stateLimits_y(idV,1,k,j,i,iens);
      real q4_L = stateLimits_y(idW,0,k,j,i,iens);   real q4_R = stateLimits_y(idW,1,k,j,i,iens);
      real q5_L = stateLimits_y(idT,0,k,j,i,iens);   real q5_R = stateLimits_y(idT,1,k,j,i,iens);
      // Compute upwind characteristics
      // Waves 1-3, velocity: v
      real w1, w2, w3;
      if (v > 0) {
        w1 = q1_L - q5_L/t;
        w2 = q2_L - u*q5_L/t;
        w3 = q4_L - w*q5_L/t;
      } else {
        w1 = q1_R - q5_R/t;
        w2 = q2_R - u*q5_R/t;
        w3 = q4_R - w*q5_R/t;
      }
      // Wave 5, velocity: v-cs
      real w5 =  v*q1_R/(2*cs) - q3_R/(2*cs) + q5_R/(2*t);
      // Wave 6, velocity: v+cs
      real w6 = -v*q1_L/(2*cs) + q3_L/(2*cs) + q5_L/(2*t);
      // Use right eigenmatrix to compute upwind flux
      real q1 = w1 + w5 + w6;
      real q2 = w2 + u*w5 + u*w6;
      real q3 = v*w1 + (v-cs)*w5 + (v+cs)*w6;
      real q4 = w3 + w*w5 + w*w6;
      real q5 =      t*w5 + t*w6;

      stateFlux_y(idR,k,j,i,iens) = q3;
      stateFlux_y(idU,k,j,i,iens) = q3*q2/q1;
      stateFlux_y(idV,k,j,i,iens) = q3*q3/q1 + C0*pow(q5,gamma);
      stateFlux_y(idW,k,j,i,iens) = q3*q4/q1;
      stateFlux_y(idT,k,j,i,iens) = q3*q5/q1;

      // COMPUTE UPWIND TRACER FLUXES
      // Handle it one tracer at a time
      for (int tr=0; tr < num_tracers; tr++) {
        if (v > 0) {
          tracerFlux_y(tr,k,j,i,iens) = q3 * tracerLimits_y(tr,0,k,j,i,iens) / r_L;
        } else {
          tracerFlux_y(tr,k,j,i,iens) = q3 * tracerLimits_y(tr,1,k,j,i,iens) / r_R;
        }
      }
    });
    //////////////////////////////////////////////////////////
    // Limit the tracer fluxes for positivity
    //////////////////////////////////////////////////////////
    real5d fct_mult_y("fct_mult_y",num_tracers,nz,ny+1,nx,nens);
    parallel_for( SimpleBounds<5>(num_tracers,nz,ny+1,nx,nens) , YAKL_LAMBDA (int tr, int k, int j, int i, int iens) {
      fct_mult_y(tr,k,j,i,iens) = 1.;
      // Solid wall BCs mean u == 0 at boundaries, so we assume periodic if u != 0
      if (tracer_pos(tr)) {
        // Compute and apply the flux reduction factor of the upwind cell
        if      (tracerFlux_y(tr,k,j,i,iens) > 0) {
          // upwind is to the left of this interface
          int ind_j = j-1;
          if (ind_j == -1) ind_j = ny-1;
          real f1 = min( tracerFlux_y(tr,k,ind_j  ,i,iens) , 0._fp );
          real f2 = max( tracerFlux_y(tr,k,ind_j+1,i,iens) , 0._fp );
          real fluxOut = dt*(f2-f1)/dy;
          real mass = tracers(tr,hs+k,hs+ind_j,hs+i,iens);
          if (fluxOut > 0) {
            fct_mult_y(tr,k,j,i,iens) = min( 1._fp , mass / fluxOut );
          }
        } else if (tracerFlux_y(tr,k,j,i,iens) < 0) {
          // upwind is to the right of this interface
          int ind_j = j;
          if (ind_j == ny) ind_j = 0;
          real f1 = min( tracerFlux_y(tr,k,ind_j  ,i,iens) , 0._fp );
          real f2 = max( tracerFlux_y(tr,k,ind_j+1,i,iens) , 0._fp );
          real fluxOut = dt*(f2-f1)/dy;
          real mass = tracers(tr,hs+k,hs+ind_j,hs+i,iens);
          if (fluxOut > 0) {
            fct_mult_y(tr,k,j,i,iens) = min( 1._fp , mass / fluxOut );
          }
        }
      }
    });



    ////////////////////////////////////////////////
    // BCs for the state edge estimates (z-direction)
    ////////////////////////////////////////////////
    parallel_for( SimpleBounds<3>(ny,nx,nens) , YAKL_LAMBDA (int j, int i, int iens) {
      for (int l = 0; l < num_state; l++) {
        if        (bc_z == BC_PERIODIC) {
          stateLimits_z     (l,0,0 ,j,i,iens) = stateLimits_z     (l,0,nz,j,i,iens);
          stateLimits_z     (l,1,nz,j,i,iens) = stateLimits_z     (l,1,0 ,j,i,iens);
        } else if (bc_z == BC_WALL    ) {
          stateLimits_z     (l,0,0 ,j,i,iens) = stateLimits_z     (l,1,0 ,j,i,iens);
          stateLimits_z     (l,1,nz,j,i,iens) = stateLimits_z     (l,0,nz,j,i,iens);
        }
      }
      for (int l = 0; l < num_tracers; l++) {
        if        (bc_z == BC_PERIODIC) {
          tracerLimits_z(l,0,0 ,j,i,iens) = tracerLimits_z(l,0,nz,j,i,iens);
          tracerLimits_z(l,1,nz,j,i,iens) = tracerLimits_z(l,1,0 ,j,i,iens);
        } else if (bc_z == BC_WALL    ) {
          tracerLimits_z(l,0,0 ,j,i,iens) = tracerLimits_z(l,1,0 ,j,i,iens);
          tracerLimits_z(l,1,nz,j,i,iens) = tracerLimits_z(l,0,nz,j,i,iens);
        }
      }
    });
    //////////////////////////////////////////////////////////
    // Compute the upwind fluxes (z-direction)
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<4>(nz+1,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      // Get left and right state
      real r_L = stateLimits_z(idR,0,k,j,i,iens)    ;   real r_R = stateLimits_z(idR,1,k,j,i,iens)    ;
      real u_L = stateLimits_z(idU,0,k,j,i,iens)/r_L;   real u_R = stateLimits_z(idU,1,k,j,i,iens)/r_R;
      real v_L = stateLimits_z(idV,0,k,j,i,iens)/r_L;   real v_R = stateLimits_z(idV,1,k,j,i,iens)/r_R;
      real w_L = stateLimits_z(idW,0,k,j,i,iens)/r_L;   real w_R = stateLimits_z(idW,1,k,j,i,iens)/r_R;
      real t_L = stateLimits_z(idT,0,k,j,i,iens)/r_L;   real t_R = stateLimits_z(idT,1,k,j,i,iens)/r_R;
      // Compute average state
      real r = 0.5_fp * (r_L + r_R);
      real u = 0.5_fp * (u_L + u_R);
      real v = 0.5_fp * (v_L + v_R);
      real w = 0.5_fp * (w_L + w_R);
      real t = 0.5_fp * (t_L + t_R);
      real p = C0 * pow(r*t,gamma);
      real cs2 = gamma*p/r;
      real cs  = sqrt(cs2);
      // Get left and right fluxes
      real q1_L = stateLimits_z(idR,0,k,j,i,iens);   real q1_R = stateLimits_z(idR,1,k,j,i,iens);
      real q2_L = stateLimits_z(idU,0,k,j,i,iens);   real q2_R = stateLimits_z(idU,1,k,j,i,iens);
      real q3_L = stateLimits_z(idV,0,k,j,i,iens);   real q3_R = stateLimits_z(idV,1,k,j,i,iens);
      real q4_L = stateLimits_z(idW,0,k,j,i,iens);   real q4_R = stateLimits_z(idW,1,k,j,i,iens);
      real q5_L = stateLimits_z(idT,0,k,j,i,iens);   real q5_R = stateLimits_z(idT,1,k,j,i,iens);
      // Compute upwind characteristics
      // Waves 1-3, velocity: w
      real w1, w2, w3;
      if (w > 0) {
        w1 = q1_L - q5_L/t;
        w2 = q2_L - u*q5_L/t;
        w3 = q3_L - v*q5_L/t;
      } else {
        w1 = q1_R - q5_R/t;
        w2 = q2_R - u*q5_R/t;
        w3 = q3_R - v*q5_R/t;
      }
      // Wave 5, velocity: w-cs
      real w5 =  w*q1_R/(2*cs) - q4_R/(2*cs) + q5_R/(2*t);
      // Wave 6, velocity: w+cs
      real w6 = -w*q1_L/(2*cs) + q4_L/(2*cs) + q5_L/(2*t);
      // Use right eigenmatrix to compute upwind flux
      real q1 = w1 + w5 + w6;
      real q2 = w2 + u*w5 + u*w6;
      real q3 = w3 + v*w5 + v*w6;
      real q4 = w*w1 + (w-cs)*w5 + (w+cs)*w6;
      real q5 =      t*w5 + t*w6;

      stateFlux_z(idR,k,j,i,iens) = q4;
      stateFlux_z(idU,k,j,i,iens) = q4*q2/q1;
      stateFlux_z(idV,k,j,i,iens) = q4*q3/q1;
      stateFlux_z(idW,k,j,i,iens) = q4*q4/q1 + C0*pow(q5,gamma);
      stateFlux_z(idT,k,j,i,iens) = q4*q5/q1;

      // COMPUTE UPWIND TRACER FLUXES
      // Handle it one tracer at a time
      for (int tr=0; tr < num_tracers; tr++) {
        if (w > 0) {
          tracerFlux_z(tr,k,j,i,iens) = q4 * tracerLimits_z(tr,0,k,j,i,iens) / r_L;
        } else {
          tracerFlux_z(tr,k,j,i,iens) = q4 * tracerLimits_z(tr,1,k,j,i,iens) / r_R;
        }
      }
    });
    //////////////////////////////////////////////////////////
    // Limit the tracer fluxes for positivity
    //////////////////////////////////////////////////////////
    real5d fct_mult_z("fct_mult_z",num_tracers,nz+1,ny,nx,nens);
    parallel_for( SimpleBounds<5>(num_tracers,nz+1,ny,nx,nens) , YAKL_LAMBDA (int tr, int k, int j, int i, int iens) {
      fct_mult_z(tr,k,j,i,iens) = 1.;
      if (k == 0 || k == nz) tracerFlux_z(tr,k,j,i,iens) = 0;
      // Solid wall BCs mean w == 0 at boundaries
      if (tracer_pos(tr)) {
        // Compute and apply the flux reduction factor of the upwind cell
        if      (tracerFlux_z(tr,k,j,i,iens) > 0) {
          int ind_k = k-1;
          // upwind is to the left of this interface
          real f1 = min( tracerFlux_z(tr,ind_k  ,j,i,iens) , 0._fp );
          real f2 = max( tracerFlux_z(tr,ind_k+1,j,i,iens) , 0._fp );
          real fluxOut = dt*(f2-f1)/dz(ind_k,iens);
          real dens = state(idR,hs+ind_k,hs+j,hs+i,iens) + hyDensCells(ind_k,iens);
          real mass = tracers(tr,hs+ind_k,hs+j,hs+i,iens) * dens;
          if (fluxOut > 0) {
            fct_mult_z(tr,k,j,i,iens) = min( 1._fp , mass / fluxOut );
          }
        } else if (tracerFlux_z(tr,k,j,i,iens) < 0) {
          int ind_k = k;
          // upwind is to the right of this interface
          real f1 = min( tracerFlux_z(tr,ind_k  ,j,i,iens) , 0._fp );
          real f2 = max( tracerFlux_z(tr,ind_k+1,j,i,iens) , 0._fp );
          real fluxOut = dt*(f2-f1)/dz(ind_k,iens);
          real dens = state(idR,hs+ind_k,hs+j,hs+i,iens) + hyDensCells(ind_k,iens);
          real mass = tracers(tr,hs+ind_k,hs+j,hs+i,iens) * dens;
          if (fluxOut > 0) {
            fct_mult_z(tr,k,j,i,iens) = min( 1._fp , mass / fluxOut );
          }
        }
      }
    });



    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA(int k, int j, int i, int iens) {
      // x-direction
      for (int l=0; l < num_state; l++) {
        if (sim2d && l == idV) {
          stateTend(l,k,j,i,iens) = 0;
        } else {
          stateTend(l,k,j,i,iens) += - ( stateFlux_x(l,k,j,i+1,iens) - stateFlux_x(l,k,j,i,iens) ) / dx;
        }
      }
      for (int l=0; l < num_tracers; l++) {
        tracerTend(l,k,j,i,iens) = - ( tracerFlux_x(l,k,j,i+1,iens)*fct_mult_x(l,k,j,i+1,iens) -
                                       tracerFlux_x(l,k,j,i  ,iens)*fct_mult_x(l,k,j,i  ,iens) ) / dx;
      }

      // y-direction
      if (! sim2d) {
        for (int l=0; l < num_state; l++) {
          stateTend(l,k,j,i,iens) += - ( stateFlux_y(l,k,j+1,i,iens) - stateFlux_y(l,k,j,i,iens) ) / dy;
        }
        for (int l=0; l < num_tracers; l++) {
          tracerTend(l,k,j,i,iens) += - ( tracerFlux_y(l,k,j+1,i,iens)*fct_mult_y(l,k,j+1,i,iens) -
                                          tracerFlux_y(l,k,j  ,i,iens)*fct_mult_y(l,k,j  ,i,iens) ) / dy;
        }
      }

      // z-direction
      for (int l=0; l < num_state; l++) {
        if (sim2d && l == idV) {
          stateTend(l,k,j,i,iens) = 0;
        } else {
          stateTend(l,k,j,i,iens) += - ( stateFlux_z(l,k+1,j,i,iens) - stateFlux_z(l,k,j,i,iens) ) / dz(k,iens);
        }
      }
      for (int l=0; l < num_tracers; l++) {
        tracerTend(l,k,j,i,iens) += - ( tracerFlux_z(l,k+1,j,i,iens)*fct_mult_z(l,k+1,j,i,iens) -
                                        tracerFlux_z(l,k  ,j,i,iens)*fct_mult_z(l,k  ,j,i,iens) ) / dz(k,iens);
      }

      // Multiply density back onto the tracers
      for (int l=0; l < num_tracers; l++) {
        tracers(l,hs+k,hs+j,hs+i,iens) *= (state(idR,hs+k,hs+j,hs+i,iens) + hyDensCells(k,iens));
      }
    });


  } // computeTendencies



  void switch_directions() {
  }


  const char * getName() { return ""; }



  template <class MICRO>
  void output(DataManager &dm, MICRO const &micro, real etime) const {
    YAKL_SCOPE( dx                    , this->dx                   );
    YAKL_SCOPE( dy                    , this->dy                   );
    YAKL_SCOPE( dz                    , this->dz                   );
    YAKL_SCOPE( hyDensCells           , this->hyDensCells          );
    YAKL_SCOPE( hyDensThetaCells      , this->hyDensThetaCells     );
    YAKL_SCOPE( hyThetaCells          , this->hyThetaCells         );
    YAKL_SCOPE( hyPressureCells       , this->hyPressureCells      );
    YAKL_SCOPE( gamma                 , this->gamma                );
    YAKL_SCOPE( C0                    , this->C0                   );

    for (int iens = 0; iens < nens; iens++) {
      std::string fname = out_prefix + std::string("_") + std::to_string(iens) + std::string(".nc");

      yakl::SimpleNetCDF nc;
      int ulIndex = 0; // Unlimited dimension index to place this data at
      // Create or open the file
      if (etime == 0.) {
        nc.create(fname);

        // x-coordinate
        real1d xloc("xloc",nx);
        parallel_for( nx , YAKL_LAMBDA (int i) { xloc(i) = (i+0.5)*dx; });
        nc.write(xloc.createHostCopy(),"x",{"x"});

        // y-coordinate
        real1d yloc("yloc",ny);
        parallel_for( ny , YAKL_LAMBDA (int i) { yloc(i) = (i+0.5)*dy; });
        nc.write(yloc.createHostCopy(),"y",{"y"});

        // z-coordinate
        auto zint = dm.get<real,2>("vertical_interface_height");
        real1d zmid("zmid",nz);
        parallel_for( nz , YAKL_LAMBDA (int i) { zmid(i) = ( zint(i,iens) + zint(i+1,iens) ) / 2; });
        nc.write(zmid.createHostCopy(),"z",{"z"});

        // hydrostatic density, theta, and pressure
        parallel_for( nz , YAKL_LAMBDA (int k) { zmid(k) = hyDensCells(k,iens); });
        nc.write(zmid.createHostCopy(),"hyDens"    ,{"z"});

        parallel_for( nz , YAKL_LAMBDA (int k) { zmid(k) = hyPressureCells(k,iens); });
        nc.write(zmid.createHostCopy(),"hyPressure",{"z"});

        parallel_for( nz , YAKL_LAMBDA (int k) { zmid(k) = hyThetaCells(k,iens); });
        nc.write(zmid.createHostCopy(),"hyTheta"   ,{"z"});

        // Create time variable
        nc.write1(0._fp,"t",0,"t");
      } else {
        nc.open(fname,yakl::NETCDF_MODE_WRITE);
        ulIndex = nc.getDimSize("t");

        // Write the elapsed time
        nc.write1(etime,"t",ulIndex,"t");
      }

      real5d state   = dm.get<real,5>("dynamics_state");
      real5d tracers = dm.get<real,5>("dynamics_tracers");

      real3d data("data",nz,ny,nx);
      // rho'
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = state(idR,hs+k,hs+j,hs+i,iens); });
      nc.write1(data.createHostCopy(),"dens_pert",{"z","y","x"},ulIndex,"t");
      // u
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        data(k,j,i) = state(idU,hs+k,hs+j,hs+i,iens) / ( state(idR,hs+k,hs+j,hs+i,iens) + hyDensCells(k,iens) );
      });
      nc.write1(data.createHostCopy(),"u",{"z","y","x"},ulIndex,"t");
      // v
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        data(k,j,i) = state(idV,hs+k,hs+j,hs+i,iens) / ( state(idR,hs+k,hs+j,hs+i,iens) + hyDensCells(k,iens) );
      });
      nc.write1(data.createHostCopy(),"v",{"z","y","x"},ulIndex,"t");
      // w
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        data(k,j,i) = state(idW,hs+k,hs+j,hs+i,iens) / ( state(idR,hs+k,hs+j,hs+i,iens) + hyDensCells(k,iens) );
      });
      nc.write1(data.createHostCopy(),"w",{"z","y","x"},ulIndex,"t");
      // theta'
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        real r =   state(idR,hs+k,hs+j,hs+i,iens) + hyDensCells     (k,iens);
        real t = ( state(idT,hs+k,hs+j,hs+i,iens) + hyDensThetaCells(k,iens) ) / r;
        data(k,j,i) = t - hyThetaCells(k,iens);
      });
      nc.write1(data.createHostCopy(),"pot_temp_pert",{"z","y","x"},ulIndex,"t");
      // pressure'
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        real r  = state(idR,hs+k,hs+j,hs+i,iens) + hyDensCells(k,iens);
        real rt = state(idT,hs+k,hs+j,hs+i,iens) + hyDensThetaCells(k,iens);
        real p  = C0*pow(rt,gamma);
        data(k,j,i) = p - hyPressureCells(k,iens);
      });
      nc.write1(data.createHostCopy(),"pressure_pert",{"z","y","x"},ulIndex,"t");

      for (int tr=0; tr < num_tracers; tr++) {
        parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
          real r = state(idR,hs+k,hs+j,hs+i,iens) + hyDensCells(k,iens);
          data(k,j,i) = tracers(tr,hs+k,hs+j,hs+i,iens)/r;
        });
        nc.write1(data.createHostCopy(),std::string("tracer_")+tracer_name[tr],{"z","y","x"},ulIndex,"t");
      }

      micro.output(dm, nc, ulIndex, iens);

      // Close the file
      nc.close();

    }
  }



  void finalize(real4d const &state , real4d const &tracers) {}



  // ord stencil values to ngll GLL values; store in gll
  YAKL_INLINE void reconstruct_gll_values( SArray<real,1,ord> const stencil ,
                                           SArray<real,1,ngll> &gll ,
                                           SArray<real,2,ord,ngll> const &coefs_to_gll ,
                                           SArray<real,2,ord,ngll> const &sten_to_gll  ,
                                           SArray<real,3,ord,ord,ord> const &wenoRecon ,
                                           SArray<real,1,hs+2> const &idl              ,
                                           real sigma, bool doweno ) {
    if (doweno) {

      // Reconstruct values
      SArray<real,1,ord> wenoCoefs;
      weno::compute_weno_coefs<ord>( wenoRecon , stencil , wenoCoefs , idl , sigma );
      // Transform ord weno coefficients into ngll GLL points
      for (int ii=0; ii<ngll; ii++) {
        real tmp = 0;
        for (int s=0; s < ord; s++) {
          tmp += coefs_to_gll(s,ii) * wenoCoefs(s);
        }
        gll(ii) = tmp;
      }

    } else {

      // Transform ord stencil cell averages into ngll GLL points
      for (int ii=0; ii<ngll; ii++) {
        real tmp = 0;
        for (int s=0; s < ord; s++) {
          tmp += sten_to_gll(s,ii) * stencil(s);
        }
        gll(ii) = tmp;
      }

    } // if doweno
  }



  YAKL_INLINE void compute_timeAvg( SArray<real,3,num_state,nAder,ngll> &dts , real dt ) {
    real dtmult = dt;
    for (int kt=1; kt<nAder; kt++) {
      for (int l=0; l<num_state; l++) {
        for (int ii=0; ii<ngll; ii++) {
          dts(l,0,ii) += dts(l,kt,ii) * dtmult / (kt+1._fp);
        }
      }
      dtmult *= dt;
    }
  }



  YAKL_INLINE void compute_timeAvg( SArray<real,2,nAder,ngll> &dts , real dt ) {
    real dtmult = dt;
    for (int kt=1; kt<nAder; kt++) {
      for (int ii=0; ii<ngll; ii++) {
        dts(0,ii) += dts(kt,ii) * dtmult / (kt+1._fp);
      }
      dtmult *= dt;
    }
  }



  YAKL_INLINE void compute_timeAvg( SArray<real,2,nAder,ngll> const &dts , SArray<real,1,ngll> &tavg , real dt ) {
    for (int ii=0; ii<ngll; ii++) {
      tavg(ii) = dts(0,ii);
    }
    real dtmult = dt;
    for (int kt=1; kt<nAder; kt++) {
      for (int ii=0; ii<ngll; ii++) {
        tavg(ii) += dts(kt,ii) * dtmult / (kt+1._fp);
      }
      dtmult *= dt;
    }
  }


};
