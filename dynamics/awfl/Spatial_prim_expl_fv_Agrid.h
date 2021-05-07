
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
  real6d stateLimits;
  real6d tracerLimits;
  real5d stateFlux;

  // Stores single-valued flux of the tracer at each cell interface
  real5d tracerFlux;

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

  // Which direction we're passing through for a given time step (x,y,z)  or (z,y,x)
  // For Strang splitting
  bool dimSwitch;

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
  // Whether to balance initial density to avoid acoustics at the start
  bool balance_initial_density;


  // When this class is created, initialize num_tracers to zero
  Spatial_operator() {
    num_tracers = 0;
    if (nAder > 1){ endrun("ERROR: nAder cannot be > 1. Please only use __FILE__ with RK time stepping."); }
  }



  // Make sure it's odd-order-accurate
  static_assert(ord%2 == 1,"ERROR: ord must be an odd integer");



  template <class MICRO>
  void convert_dynamics_to_coupler_state( DataManager &dm , MICRO &micro ) {
    // real5d state       = dm.get<real,5>( "dynamics_state"   );
    // real5d tracers     = dm.get<real,5>( "dynamics_tracers" );
    // real4d dm_dens_dry = dm.get<real,4>( "density_dry"      );
    // real4d dm_uvel     = dm.get<real,4>( "uvel"             );
    // real4d dm_vvel     = dm.get<real,4>( "vvel"             );
    // real4d dm_wvel     = dm.get<real,4>( "wvel"             );
    // real4d dm_temp     = dm.get<real,4>( "temp"             );
    //
    // YAKL_SCOPE( hyDensCells      , this->hyDensCells      );
    // YAKL_SCOPE( hyDensThetaCells , this->hyDensThetaCells );
    // YAKL_SCOPE( C0               , this->C0               );
    // YAKL_SCOPE( gamma            , this->gamma            );
    // YAKL_SCOPE( num_tracers      , this->num_tracers      );
    // YAKL_SCOPE( p0               , this->p0               );
    // YAKL_SCOPE( Rd               , this->Rd               );
    // YAKL_SCOPE( Rv               , this->Rv               );
    // YAKL_SCOPE( cp               , this->cp               );
    // YAKL_SCOPE( tracer_adds_mass , this->tracer_adds_mass );
    //
    // int idWV = micro.get_water_vapor_index();
    //
    // MultipleFields<max_tracers,real4d> dm_tracers;
    // for (int tr = 0; tr < num_tracers; tr++) {
    //   auto trac = dm.get<real,4>( tracer_name[tr] );
    //   dm_tracers.add_field( trac );
    // }
    //
    // parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
    //   real dens  = state(idR,hs+k,hs+j,hs+i,iens) + hyDensCells(hs+k,iens);
    //   real uvel  = state(idU,hs+k,hs+j,hs+i,iens);
    //   real vvel  = state(idV,hs+k,hs+j,hs+i,iens);
    //   real wvel  = state(idW,hs+k,hs+j,hs+i,iens);
    //   real theta = state(idT,hs+k,hs+j,hs+i,iens) + hyThetaCells(hs+k,iens);
    //   real pressure = C0 * pow( dens*theta , gamma );
    //   real dens_vap = tracers(idWV,hs+k,hs+j,hs+i,iens) * dens;
    //   real dens_dry = dens;
    //   for (int tr=0; tr < num_tracers; tr++) {
    //     if (tracer_adds_mass(tr)) dens_dry -= tracers(tr,hs+k,hs+j,hs+i,iens) * dens;
    //   }
    //   real temp = pressure / ( dens_dry * Rd + dens_vap * Rv );
    //   dm_dens_dry(k,j,i,iens) = dens_dry;
    //   dm_uvel    (k,j,i,iens) = uvel;
    //   dm_vvel    (k,j,i,iens) = vvel;
    //   dm_wvel    (k,j,i,iens) = wvel;
    //   dm_temp    (k,j,i,iens) = temp;
    //   for (int tr=0; tr < num_tracers; tr++) {
    //     dm_tracers(tr,k,j,i,iens) = tracers(tr,hs+k,hs+j,hs+i,iens) * dens;
    //   }
    // });
  }



  template <class MICRO>
  void convert_coupler_state_to_dynamics( DataManager &dm , MICRO &micro ) {
    // real5d state           = dm.get<real,5>( "dynamics_state"   );
    // real5d tracers         = dm.get<real,5>( "dynamics_tracers" );
    // real4d dm_dens_dry     = dm.get<real,4>( "density_dry"      );
    // real4d dm_uvel         = dm.get<real,4>( "uvel"             );
    // real4d dm_vvel         = dm.get<real,4>( "vvel"             );
    // real4d dm_wvel         = dm.get<real,4>( "wvel"             );
    // real4d dm_temp         = dm.get<real,4>( "temp"             );
    //
    // YAKL_SCOPE( hyDensCells      , this->hyDensCells      );
    // YAKL_SCOPE( hyDensThetaCells , this->hyDensThetaCells );
    // YAKL_SCOPE( C0               , this->C0               );
    // YAKL_SCOPE( gamma            , this->gamma            );
    // YAKL_SCOPE( num_tracers      , this->num_tracers      );
    // YAKL_SCOPE( p0               , this->p0               );
    // YAKL_SCOPE( Rd               , this->Rd               );
    // YAKL_SCOPE( Rv               , this->Rv               );
    // YAKL_SCOPE( cp               , this->cp               );
    // YAKL_SCOPE( tracer_adds_mass , this->tracer_adds_mass );
    //
    // int idWV = micro.get_water_vapor_index();
    //
    // MultipleFields<max_tracers,real4d> dm_tracers;
    // for (int tr = 0; tr < num_tracers; tr++) {
    //   auto trac = dm.get<real,4>( tracer_name[tr] );
    //   dm_tracers.add_field( trac );
    // }
    //
    // parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
    //   for (int tr=0; tr < num_tracers; tr++) {
    //     tracers(tr,hs+k,hs+j,hs+i,iens) = dm_tracers(tr,k,j,i,iens);
    //   }
    //   real dens_dry = dm_dens_dry(k,j,i,iens);
    //   real uvel     = dm_uvel    (k,j,i,iens);
    //   real vvel     = dm_vvel    (k,j,i,iens);
    //   real wvel     = dm_wvel    (k,j,i,iens);
    //   real temp     = dm_temp    (k,j,i,iens);
    //   real dens_vap = tracers(idWV,hs+k,hs+j,hs+i,iens);
    //   real dens     = dens_dry;
    //   for (int tr=0; tr < num_tracers; tr++) {
    //     if (tracer_adds_mass(tr)) dens += tracers(tr,hs+k,hs+j,hs+i,iens);
    //   }
    //   real pressure = dens_dry * Rd * temp + dens_vap * Rv * temp;
    //   real theta    = pow( pressure / C0 , 1._fp / gamma ) / dens;
    //   for (int tr=0; tr < num_tracers; tr++) {
    //     tracers(tr,hs+k,hs+j,hs+i,iens) /= dens;
    //   }
    //   state(idR,hs+k,hs+j,hs+i,iens) = dens - hyDensCells(hs+k,iens);
    //   state(idU,hs+k,hs+j,hs+i,iens) = uvel;
    //   state(idV,hs+k,hs+j,hs+i,iens) = vvel;
    //   state(idW,hs+k,hs+j,hs+i,iens) = wvel;
    //   state(idT,hs+k,hs+j,hs+i,iens) = theta - hyThetaCells(hs+k,iens);
    // });
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
            tracers(idWV,hs+k,hs+j,hs+i,iens) += r_v / (rh+r_v) * wt;
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
    YAKL_SCOPE( balance_initial_density , this->balance_initial_density );
    YAKL_SCOPE( Rd                      , this->Rd                      );
    YAKL_SCOPE( Rv                      , this->Rv                      );
    YAKL_SCOPE( C0                      , this->C0                      );
    YAKL_SCOPE( gamma                   , this->gamma                   );

    // Copy the DataManager data to state and tracer arrays for convenience
    real5d state   = dm.get<real,5>("dynamics_state");
    real5d tracers = dm.get<real,5>("dynamics_tracers");

    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      // Add tracer density to dry density if it adds mass
      for (int tr=0; tr < num_tracers; tr++) {
        if (tracer_adds_mass(tr)) {
          state(idR,hs+k,hs+j,hs+i,iens) += tracers(tr,hs+k,hs+j,hs+i,iens);
        }
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
    return 3;
  }



  // Given the model data and CFL value, compute the maximum stable time step
  template <class MICRO>
  real compute_time_step(DataManager &dm, MICRO const &micro, real cfl = 0.8) {

    // If we've already computed the time step, then don't compute it again
    if (dtInit <= 0) {
      YAKL_SCOPE( dx               , this->dx               );
      YAKL_SCOPE( dy               , this->dy               );
      YAKL_SCOPE( dz               , this->dz               );
      YAKL_SCOPE( hyDensCells      , this->hyDensCells      );
      YAKL_SCOPE( hyDensThetaCells , this->hyDensThetaCells );
      YAKL_SCOPE( gamma            , this->gamma            );
      YAKL_SCOPE( C0               , this->C0               );

      // Convert data from DataManager to state and tracers array for convenience
      real5d state   = dm.get<real,5>("dynamics_state");
      real5d tracers = dm.get<real,5>("dynamics_tracers");

      // Allocate a 3-D array for the max stable time steps (we'll use this for a reduction later)
      real4d dt3d("dt3d",nz,ny,nx,nens);

      // Loop through the cells, calculate the max stable time step for each cell
      parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        // Get the state
        real r = state(idR,hs+k,hs+j,hs+i,iens) + hyDensCells(hs+k,iens);
        real u = state(idU,hs+k,hs+j,hs+i,iens);
        real v = state(idV,hs+k,hs+j,hs+i,iens);
        real w = state(idW,hs+k,hs+j,hs+i,iens);
        real t = state(idT,hs+k,hs+j,hs+i,iens) + hyThetaCells(hs+k,iens);
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
    dimSwitch = true;

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

      balance_initial_density = config["balance_initial_density"].as<bool>();

      sim_time = config["simTime"].as<real>();
    #else
      weno_scalars            = true;
      weno_winds              = true;
      data_spec               = DATA_SPEC_SUPERCELL;
      bc_x                    = BC_PERIODIC;
      bc_y                    = BC_PERIODIC;
      bc_z                    = BC_WALL;
      out_prefix              = "test";
      balance_initial_density = false;
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
    stateLimits     = real6d("stateLimits"    ,num_state  ,2,nz+1,ny+1,nx+1,nens);
    tracerLimits    = real6d("tracerLimits"   ,num_tracers,2,nz+1,ny+1,nx+1,nens);
    stateFlux       = real5d("stateFlux"      ,num_state    ,nz+1,ny+1,nx+1,nens);
    tracerFlux      = real5d("tracerFlux"     ,num_tracers  ,nz+1,ny+1,nx+1,nens);
    hyDensCells      = real2d("hyDensCells       ",nz+2*hs,nens);
    hyPressureCells  = real2d("hyPressureCells   ",nz+2*hs,nens);
    hyThetaCells     = real2d("hyThetaCells      ",nz+2*hs,nens);
    hyDensThetaCells = real2d("hyDensThetaCells  ",nz+2*hs,nens);
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
    YAKL_SCOPE( balance_initial_density  , this->balance_initial_density );
    YAKL_SCOPE( vert_interface           , this->vert_interface          );
    YAKL_SCOPE( vert_interface_ghost     , this->vert_interface_ghost    );

    real5d state   = dm.get<real,5>("dynamics_state"  );
    real5d tracers = dm.get<real,5>("dynamics_tracers");

    // If the data_spec is thermal or ..., then initialize the domain with Exner pressure-based hydrostasis
    // This is mostly to make plotting potential temperature perturbation easier for publications
    if (data_spec == DATA_SPEC_THERMAL) {

      // Setup hydrostatic background state
      parallel_for( SimpleBounds<2>(nz+2*hs,nens) , YAKL_LAMBDA (int k, int iens) {
        // Compute cell averages
        hyDensCells     (k,iens) = 0;
        hyPressureCells (k,iens) = 0;
        hyThetaCells    (k,iens) = 0;
        hyDensThetaCells(k,iens) = 0;
        for (int kk=0; kk<ord; kk++) {
          real zloc = vert_interface_ghost(k,iens) + 0.5_fp*dz_ghost(k,iens) + gllPts_ord(kk)*dz_ghost(k,iens);
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
                // Line below balances initial density to remove the acoustic wave
                if (balance_initial_density) r = rh*th/t;

                state(idR,hs+k,hs+j,hs+i,iens) += (r - rh)*wt;
                state(idT,hs+k,hs+j,hs+i,iens) += (t - th)*wt;
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
      real2d hyDensVapCells("hyDensVapCells",nz+2*hs,nens);
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
        hyPressureCells (hs+k,iens) = press;
        hyDensCells     (hs+k,iens) = dens;
        hyDensThetaCells(hs+k,iens) = dens_theta;
        hyThetaCells    (hs+k,iens) = theta;
        hyDensVapCells  (hs+k,iens) = dens_vap;
      });

      // Dump out data to plot a skew-T log-P diagram
      yakl::SimpleNetCDF nc;
      nc.create("skew.nc");
      real1d data("data",nz);
      for (int iens=0; iens < nens; iens++) {
        parallel_for( nz , YAKL_LAMBDA (int k) { data(k) = z(k,iens); });
        nc.write(data.createHostCopy(),"z"          ,{"z"});

        parallel_for( nz , YAKL_LAMBDA (int k) { data(k) = hyPressureCells(hs+k,iens); });
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
              state  (idR ,hs+k,hs+j,hs+i,iens) += (dens - hyDensGLL(k,kk,iens)) * factor;
              state  (idU ,hs+k,hs+j,hs+i,iens) += uvel                          * factor;
              state  (idV ,hs+k,hs+j,hs+i,iens) += vvel                          * factor;
              state  (idW ,hs+k,hs+j,hs+i,iens) += wvel                          * factor;
              state  (idT ,hs+k,hs+j,hs+i,iens) += theta_pert                    * factor;
              tracers(idWV,hs+k,hs+j,hs+i,iens) += dens_vap / dens               * factor;
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
    if (dimSwitch) {
      if        (splitIndex == 0) {
        computeTendenciesX( state , stateTend , tracers , tracerTend , micro , dt );
      } else if (splitIndex == 1) {
        if (sim2d) {
          memset(stateTend  , 0._fp);
          memset(tracerTend , 0._fp);
        } else {
          computeTendenciesY( state , stateTend , tracers , tracerTend , micro , dt );
        }
      } else if (splitIndex == 2) {
        computeTendenciesZ( state , stateTend , tracers , tracerTend , micro , dt );
      }
    } else {
      if        (splitIndex == 0) {
        computeTendenciesZ( state , stateTend , tracers , tracerTend , micro , dt );
      } else if (splitIndex == 1) {
        if (sim2d) {
          memset(stateTend  , 0._fp);
          memset(tracerTend , 0._fp);
        } else {
          computeTendenciesY( state , stateTend , tracers , tracerTend , micro , dt );
        }
      } else if (splitIndex == 2) {
        computeTendenciesX( state , stateTend , tracers , tracerTend , micro , dt );
      }
    }
  } // computeTendencies



  void switch_directions() {
    dimSwitch = ! dimSwitch;
  }



  template <class MICRO>
  void computeTendenciesX( real5d &state   , real5d &stateTend  ,
                           real5d &tracers , real5d &tracerTend ,
                           MICRO const &micro, real &dt ) {
    YAKL_SCOPE( nx                      , this->nx                     );
    YAKL_SCOPE( weno_scalars            , this->weno_scalars           );
    YAKL_SCOPE( weno_winds              , this->weno_winds             );
    YAKL_SCOPE( c2g                     , this->coefs_to_gll           );
    YAKL_SCOPE( s2g                     , this->sten_to_gll            );
    YAKL_SCOPE( c2d2g                   , this->coefs_to_deriv_gll     );
    YAKL_SCOPE( s2d2g                   , this->sten_to_deriv_gll      );
    YAKL_SCOPE( wenoRecon               , this->wenoRecon              );
    YAKL_SCOPE( idl                     , this->idl                    );
    YAKL_SCOPE( sigma                   , this->sigma                  );
    YAKL_SCOPE( hyDensCells             , this->hyDensCells            );
    YAKL_SCOPE( hyDensThetaCells        , this->hyDensThetaCells       );
    YAKL_SCOPE( sim2d                   , this->sim2d                  );
    YAKL_SCOPE( derivMatrix             , this->derivMatrix            );
    YAKL_SCOPE( dx                      , this->dx                     );
    YAKL_SCOPE( stateLimits             , this->stateLimits            );
    YAKL_SCOPE( tracerLimits            , this->tracerLimits           );
    YAKL_SCOPE( stateFlux               , this->stateFlux              );
    YAKL_SCOPE( tracerFlux              , this->tracerFlux             );
    YAKL_SCOPE( tracer_pos              , this->tracer_pos             );
    YAKL_SCOPE( num_tracers             , this->num_tracers            );
    YAKL_SCOPE( bc_x                    , this->bc_x                   );
    YAKL_SCOPE( Rd                      , this->Rd                     );
    YAKL_SCOPE( cp                      , this->cp                     );
    YAKL_SCOPE( gamma                   , this->gamma                  );
    YAKL_SCOPE( p0                      , this->p0                     );
    YAKL_SCOPE( C0                      , this->C0                     );

    // Populate the halos
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

    // Loop through all cells, reconstruct in x-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      // We need density and momentum to evolve the tracers with ADER
      SArray<real,1,ngll> u_gll;

      { // State
        SArray<real,1,ngll> r_gll, dr_gll, du_gll, v_gll, dv_gll, w_gll, dw_gll, t_gll, dt_gll, dp_gll;
        { // Reconstruct
          SArray<real,1,ord> stencil;
          // Density
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idR,hs+k,hs+j,i+ii,iens) + hyDensCells(hs+k,iens); }
          reconstruct_gll_values_and_derivs( stencil , r_gll , dr_gll , dx , c2g , c2d2g , s2g , s2d2g , wenoRecon ,
                                             idl , sigma , weno_scalars );

          // u values and derivatives
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idU,hs+k,hs+j,i+ii,iens); }
          reconstruct_gll_values_and_derivs( stencil , u_gll , du_gll , dx , c2g , c2d2g , s2g , s2d2g , wenoRecon ,
                                             idl , sigma , weno_winds );

          // v
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idV,hs+k,hs+j,i+ii,iens); }
          reconstruct_gll_values_and_derivs( stencil , v_gll , dv_gll , dx , c2g , c2d2g , s2g , s2d2g , wenoRecon ,
                                             idl , sigma , weno_winds );

          // w
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idW,hs+k,hs+j,i+ii,iens); }
          reconstruct_gll_values_and_derivs( stencil , w_gll , dw_gll , dx , c2g , c2d2g , s2g , s2d2g , wenoRecon ,
                                             idl , sigma , weno_winds );

          // theta
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idT,hs+k,hs+j,i+ii,iens) + hyThetaCells(hs+k,iens); }
          reconstruct_gll_values_and_derivs( stencil , t_gll , dt_gll , dx , c2g , c2d2g , s2g , s2d2g , wenoRecon ,
                                             idl , sigma , weno_scalars );

          // pressure perturbation
          for (int ii=0; ii < ord; ii++) {
            real r = state(idR,hs+k,hs+j,i+ii,iens) + hyDensCells (hs+k,iens);
            real t = state(idT,hs+k,hs+j,i+ii,iens) + hyThetaCells(hs+k,iens);
            stencil(ii) = C0 * pow( r*t , gamma ) - hyPressureCells(hs+k,iens);
          }
          reconstruct_gll_derivs( stencil , dp_gll , dx , c2d2g , s2d2g , wenoRecon , idl , sigma , weno_scalars );
        } // Reconstruct

        // Compute central integral terms
        stateTend(idR,k,j,i,iens) = 0;
        stateTend(idU,k,j,i,iens) = 0;
        stateTend(idV,k,j,i,iens) = 0;
        stateTend(idW,k,j,i,iens) = 0;
        stateTend(idT,k,j,i,iens) = 0;
        for (int ii=0; ii < ngll; ii++) {
          real r  = r_gll (ii);
          real dr = dr_gll(ii);
          real u  = u_gll (ii);
          real du = du_gll(ii);
          real dv = dv_gll(ii);
          real dw = dw_gll(ii);
          real dt = dt_gll(ii);
          real dp = dp_gll(ii);
          stateTend(idR,k,j,i,iens) += -(u*dr + r*du) * gllWts_ngll(ii);
          stateTend(idU,k,j,i,iens) += -(u*du + dp/r) * gllWts_ngll(ii);
          stateTend(idV,k,j,i,iens) += -(u*dv       ) * gllWts_ngll(ii);
          stateTend(idW,k,j,i,iens) += -(u*dw       ) * gllWts_ngll(ii);
          stateTend(idT,k,j,i,iens) += -(u*dt       ) * gllWts_ngll(ii);
        }

        // Left interface
        stateLimits(idR,1,k,j,i  ,iens) = r_gll(0     );
        stateLimits(idU,1,k,j,i  ,iens) = u_gll(0     );
        stateLimits(idV,1,k,j,i  ,iens) = v_gll(0     );
        stateLimits(idW,1,k,j,i  ,iens) = w_gll(0     );
        stateLimits(idT,1,k,j,i  ,iens) = t_gll(0     );
        // Right interface
        stateLimits(idR,0,k,j,i+1,iens) = r_gll(ngll-1);
        stateLimits(idU,0,k,j,i+1,iens) = u_gll(ngll-1);
        stateLimits(idV,0,k,j,i+1,iens) = v_gll(ngll-1);
        stateLimits(idW,0,k,j,i+1,iens) = w_gll(ngll-1);
        stateLimits(idT,0,k,j,i+1,iens) = t_gll(ngll-1);
      } // State

      { // Tracers
        for (int tr=0; tr < num_tracers; tr++) {
          SArray<real,1,ngll> t_gll, dt_gll;
          { // Reconstruction
            SArray<real,1,ord> stencil;
            for (int ii=0; ii < ord; ii++) { stencil(ii) = tracers(tr,hs+k,hs+j,i+ii,iens); }
            reconstruct_gll_values_and_derivs( stencil , t_gll , dt_gll , dx , c2g , c2d2g , s2g , s2d2g , wenoRecon ,
                                               idl , sigma , weno_scalars );
            if (tracer_pos(tr)) {
              for (int ii=0; ii < ngll; ii++) { t_gll(ii) = max( 0._fp , t_gll(ii) ); }
            }
          } // Reconstruction

          // Compute central integral terms
          tracerTend(tr,k,j,i,iens) = 0;
          for (int ii=0; ii < ngll; ii++) {
            real dt = dt_gll(ii);
            real u  = u_gll (ii);
            tracerTend(tr,k,j,i,iens) += -(u*dt) * gllWts_ngll(ii);
          }

          tracerLimits(tr,1,k,j,i  ,iens) = t_gll(0     );
          tracerLimits(tr,0,k,j,i+1,iens) = t_gll(ngll-1);
        }
      } // Tracers

    });

    ////////////////////////////////////////////////
    // BCs for the state edge estimates
    ////////////////////////////////////////////////
    parallel_for( SimpleBounds<3>(nz,ny,nens) , YAKL_LAMBDA (int k, int j, int iens) {
      for (int l=0; l < num_state; l++) {
        if        (bc_x == BC_PERIODIC) {
          stateLimits(l,0,k,j,0 ,iens) = stateLimits(l,0,k,j,nx,iens);
          stateLimits(l,1,k,j,nx,iens) = stateLimits(l,1,k,j,0 ,iens);
        } else if (bc_x == BC_WALL    ) {
          stateLimits(l,0,k,j,0 ,iens) = stateLimits(l,1,k,j,0 ,iens);
          stateLimits(l,1,k,j,nx,iens) = stateLimits(l,0,k,j,nx,iens);
        }
      }
      for (int l=0; l < num_tracers; l++) {
        if        (bc_x == BC_PERIODIC) {
          tracerLimits(l,0,k,j,0 ,iens) = tracerLimits(l,0,k,j,nx,iens);
          tracerLimits(l,1,k,j,nx,iens) = tracerLimits(l,1,k,j,0 ,iens);
        } else if (bc_x == BC_WALL    ) {
          tracerLimits(l,0,k,j,0 ,iens) = tracerLimits(l,1,k,j,0 ,iens);
          tracerLimits(l,1,k,j,nx,iens) = tracerLimits(l,0,k,j,nx,iens);
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the upwind fluxes
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<4>(nz,ny,nx+1,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      // Get left and right state
      real r_L = stateLimits(idR,0,k,j,i,iens);   real r_R = stateLimits(idR,1,k,j,i,iens);
      real u_L = stateLimits(idU,0,k,j,i,iens);   real u_R = stateLimits(idU,1,k,j,i,iens);
      real v_L = stateLimits(idV,0,k,j,i,iens);   real v_R = stateLimits(idV,1,k,j,i,iens);
      real w_L = stateLimits(idW,0,k,j,i,iens);   real w_R = stateLimits(idW,1,k,j,i,iens);
      real t_L = stateLimits(idT,0,k,j,i,iens);   real t_R = stateLimits(idT,1,k,j,i,iens);
      real p_L = C0 * pow( r_L*t_L , gamma )  ;   real p_R = C0 * pow( r_R*t_R , gamma );
      // Compute average state
      real r = 0.5_fp * (r_L + r_R);
      real u = 0.5_fp * (u_L + u_R);
      real t = 0.5_fp * (t_L + t_R);
      real p = 0.5_fp * (p_L + p_R);
      real cs2 = gamma*p/r;
      real cs  = sqrt(cs2);
      // Compute state difference across the cell interface
      real dr = r_R - r_L;
      real du = u_R - u_L;
      real dv = v_R - v_L;
      real dw = w_R - w_L;
      real dt = t_R - t_L;
      real dp = p_R - p_L;
      // Compute flux difference across the cell interface
      real df1 = u*dr + r*du;
      real df2 = u*du + dp/r;
      real df3 = u*dv;
      real df4 = u*dw;
      real df5 = u*dt;
      // Compute characteristic variables
      real w1 = 0.5_fp * ( df1 - r/cs*df2 + r/t*df5 );
      real w2 = 0.5_fp * ( df1 + r/cs*df2 + r/t*df5 );
      real w3 = -r/t*df5;
      real w4 = df3;
      real w5 = df4;
      // Compute flux waves
      for (int l=0; l < num_state; l++) {
        stateLimits(l,0,k,j,i,iens) = 0;
        stateLimits(l,1,k,j,i,iens) = 0;
      }
      // WAVE 1: u-cs
      stateLimits(idR,0,k,j,i,iens) += w1;
      stateLimits(idU,0,k,j,i,iens) += -cs/r*w1;
      // WAVE 2: u+cs
      stateLimits(idR,1,k,j,i,iens) += w2;
      stateLimits(idU,1,k,j,i,iens) += cs/r*w2;
      // WAVES 3-5: u
      int upw_ind = 0;
      if (u > 0) {
        upw_ind = 1;
      }
      // WAVE 3
      stateLimits(idR,upw_ind,k,j,i,iens) += w3;
      stateLimits(idT,upw_ind,k,j,i,iens) += -t/r*w3;
      // WAVE 4
      stateLimits(idV,upw_ind,k,j,i,iens) += w4;
      // WAVE 5
      stateLimits(idW,upw_ind,k,j,i,iens) += w5;

      for (int tr=0; tr < num_tracers; tr++) {
        dt = tracerLimits(tr,1,k,j,i,iens) - tracerLimits(tr,0,k,j,i,iens);
        tracerLimits(tr,0,k,j,i,iens) = 0;
        tracerLimits(tr,1,k,j,i,iens) = 0;
        tracerLimits(tr,upw_ind,k,j,i,iens) += u*dt;
      }
    });

    //////////////////////////////////////////////////////////
    // Limit the tracer fluxes for positivity
    //////////////////////////////////////////////////////////
    real5d fct_mult("fct_mult",num_tracers,nz,ny,nx+1,nens);
    parallel_for( SimpleBounds<5>(num_tracers,nz,ny,nx+1,nens) , YAKL_LAMBDA (int tr, int k, int j, int i, int iens) {
      fct_mult(tr,k,j,i,iens) = 1.;
      // TODO:
    });

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA(int k, int j, int i, int iens) {
      for (int l = 0; l < num_state; l++) {
        if (sim2d && l == idV) {
          stateTend(l,k,j,i,iens) = 0;
        } else {
          stateTend(l,k,j,i,iens) += - ( stateLimits(l,0,k,j,i+1,iens) + stateLimits(l,1,k,j,i,iens) ) / dx;
        }
      }
      for (int l = 0; l < num_tracers; l++) {
          tracerTend(l,k,j,i,iens) += - ( tracerLimits(l,0,k,j,i+1,iens) + tracerLimits(l,1,k,j,i,iens) ) / dx;
      }
    });
  }



  template <class MICRO>
  void computeTendenciesY( real5d &state   , real5d &stateTend  ,
                           real5d &tracers , real5d &tracerTend ,
                           MICRO const &micro, real &dt ) {
    YAKL_SCOPE( ny                      , this->ny                     );
    YAKL_SCOPE( weno_scalars            , this->weno_scalars           );
    YAKL_SCOPE( weno_winds              , this->weno_winds             );
    YAKL_SCOPE( c2g                     , this->coefs_to_gll           );
    YAKL_SCOPE( s2g                     , this->sten_to_gll            );
    YAKL_SCOPE( wenoRecon               , this->wenoRecon              );
    YAKL_SCOPE( idl                     , this->idl                    );
    YAKL_SCOPE( sigma                   , this->sigma                  );
    YAKL_SCOPE( hyDensCells             , this->hyDensCells            );
    YAKL_SCOPE( hyDensThetaCells        , this->hyDensThetaCells       );
    YAKL_SCOPE( sim2d                   , this->sim2d                  );
    YAKL_SCOPE( derivMatrix             , this->derivMatrix            );
    YAKL_SCOPE( dy                      , this->dy                     );
    YAKL_SCOPE( stateLimits             , this->stateLimits            );
    YAKL_SCOPE( stateFlux               , this->stateFlux              );
    YAKL_SCOPE( tracerLimits            , this->tracerLimits           );
    YAKL_SCOPE( tracerFlux              , this->tracerFlux             );
    YAKL_SCOPE( tracer_pos              , this->tracer_pos             );
    YAKL_SCOPE( num_tracers             , this->num_tracers            );
    YAKL_SCOPE( bc_y                    , this->bc_y                   );
    YAKL_SCOPE( Rd                      , this->Rd                     );
    YAKL_SCOPE( cp                      , this->cp                     );
    YAKL_SCOPE( gamma                   , this->gamma                  );
    YAKL_SCOPE( p0                      , this->p0                     );
    YAKL_SCOPE( C0                      , this->C0                     );

    memset( stateTend  , 0._fp );
    memset( tracerTend , 0._fp );
  }



  template <class MICRO>
  void computeTendenciesZ( real5d &state   , real5d &stateTend  ,
                           real5d &tracers , real5d &tracerTend ,
                           MICRO const &micro, real &dt ) {
    YAKL_SCOPE( nz                      , this->nz                     );
    YAKL_SCOPE( weno_scalars            , this->weno_scalars           );
    YAKL_SCOPE( weno_winds              , this->weno_winds             );
    YAKL_SCOPE( c2g                     , this->coefs_to_gll           );
    YAKL_SCOPE( s2g                     , this->sten_to_gll            );
    YAKL_SCOPE( c2d2g                   , this->coefs_to_deriv_gll     );
    YAKL_SCOPE( s2d2g                   , this->sten_to_deriv_gll      );
    YAKL_SCOPE( wenoRecon               , this->wenoRecon              );
    YAKL_SCOPE( idl                     , this->idl                    );
    YAKL_SCOPE( sigma                   , this->sigma                  );
    YAKL_SCOPE( hyDensCells             , this->hyDensCells            );
    YAKL_SCOPE( hyDensGLL               , this->hyDensGLL              );
    YAKL_SCOPE( hyDensThetaGLL          , this->hyDensThetaGLL         );
    YAKL_SCOPE( hyPressureGLL           , this->hyPressureGLL          );
    YAKL_SCOPE( sim2d                   , this->sim2d                  );
    YAKL_SCOPE( derivMatrix             , this->derivMatrix            );
    YAKL_SCOPE( dz                      , this->dz                     );
    YAKL_SCOPE( stateLimits             , this->stateLimits            );
    YAKL_SCOPE( stateFlux               , this->stateFlux              );
    YAKL_SCOPE( tracerLimits            , this->tracerLimits           );
    YAKL_SCOPE( tracerFlux              , this->tracerFlux             );
    YAKL_SCOPE( tracer_pos              , this->tracer_pos             );
    YAKL_SCOPE( num_tracers             , this->num_tracers            );
    YAKL_SCOPE( bc_z                    , this->bc_z                   );
    YAKL_SCOPE( gllWts_ngll             , this->gllWts_ngll            );
    YAKL_SCOPE( Rd                      , this->Rd                     );
    YAKL_SCOPE( cp                      , this->cp                     );
    YAKL_SCOPE( gamma                   , this->gamma                  );
    YAKL_SCOPE( p0                      , this->p0                     );
    YAKL_SCOPE( C0                      , this->C0                     );

    // Populate the halos
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

    // Loop through all cells, reconstruct in x-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      // We need density and momentum to evolve the tracers with ADER
      SArray<real,1,ngll> w_gll;

      { // State
        SArray<real,1,ngll> r_gll, dr_gll, u_gll, du_gll, v_gll, dv_gll, dw_gll, t_gll, dt_gll, dp_gll;
        { // Reconstruct
          SArray<real,1,ord> stencil;
          // Density
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idR,k+kk,hs+j,hs+i,iens); }
          reconstruct_gll_values( stencil , r_gll , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
          for (int kk=0; kk < ngll; kk++) { r_gll(kk) += hyDensGLL(k,kk,iens); }

          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idR,k+kk,hs+j,hs+i,iens) + hyDensCells(k+kk,iens); }
          reconstruct_gll_derivs( stencil , dr_gll , dz(k,iens) , c2d2g , s2d2g , wenoRecon , idl , sigma , weno_scalars );

          // u values and derivatives
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idU,k+kk,hs+j,hs+i,iens); }
          reconstruct_gll_values_and_derivs( stencil , u_gll , du_gll , dz(k,iens) , c2g , c2d2g , s2g , s2d2g , wenoRecon ,
                                             idl , sigma , weno_winds );

          // v
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idV,k+kk,hs+j,hs+i,iens); }
          reconstruct_gll_values_and_derivs( stencil , v_gll , dv_gll , dz(k,iens) , c2g , c2d2g , s2g , s2d2g , wenoRecon ,
                                             idl , sigma , weno_winds );

          // w
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idW,k+kk,hs+j,hs+i,iens); }
          reconstruct_gll_values_and_derivs( stencil , w_gll , dw_gll , dz(k,iens) , c2g , c2d2g , s2g , s2d2g , wenoRecon ,
                                             idl , sigma , weno_winds );

          // theta
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idT,k+kk,hs+j,hs+i,iens) + hyThetaCells(k+kk,iens); }
          reconstruct_gll_values_and_derivs( stencil , t_gll , dt_gll , dz(k,iens) , c2g , c2d2g , s2g , s2d2g , wenoRecon ,
                                             idl , sigma , weno_scalars );

          // pressure perturbation
          for (int kk=0; kk < ord; kk++) {
            real rh = hyDensCells (k+kk,iens);
            real th = hyThetaCells(k+kk,iens);
            real r = state(idR,k+kk,hs+j,hs+i,iens) + rh;
            real t = state(idT,k+kk,hs+j,hs+i,iens) + th;
            stencil(kk) = C0 * pow( r*t , gamma ) - C0 * pow( rh*th , gamma );
          }
          reconstruct_gll_derivs( stencil , dp_gll , dz(k,iens) , c2d2g , s2d2g , wenoRecon , idl , sigma , weno_scalars );

          if (k == 0) {
            w_gll (0) = 0;
            du_gll(0) = 0;
            dv_gll(0) = 0;
            dw_gll(0) = 0;
            dt_gll(0) = 0;
            dp_gll(0) = 0;
          }
          if (k == nz-1) {
            w_gll (ngll-1) = 0;
            du_gll(ngll-1) = 0;
            dv_gll(ngll-1) = 0;
            dw_gll(ngll-1) = 0;
            dt_gll(ngll-1) = 0;
            dp_gll(ngll-1) = 0;
          }
        } // Reconstruct

        // Compute central integral terms
        stateTend(idR,k,j,i,iens) = 0;
        stateTend(idU,k,j,i,iens) = 0;
        stateTend(idV,k,j,i,iens) = 0;
        stateTend(idW,k,j,i,iens) = 0;
        stateTend(idT,k,j,i,iens) = 0;
        for (int kk=0; kk < ngll; kk++) {
          real r  = r_gll (kk);
          real dr = dr_gll(kk);
          real w  = w_gll (kk);
          real du = du_gll(kk);
          real dv = dv_gll(kk);
          real dw = dw_gll(kk);
          real dt = dt_gll(kk);
          real dp = dp_gll(kk);
          stateTend(idR,k,j,i,iens) += -(w*dr + r*dw) * gllWts_ngll(kk);
          stateTend(idU,k,j,i,iens) += -(w*du       ) * gllWts_ngll(kk);
          stateTend(idV,k,j,i,iens) += -(w*dv       ) * gllWts_ngll(kk);
          stateTend(idW,k,j,i,iens) += -(w*dw + dp/r) * gllWts_ngll(kk);
          stateTend(idT,k,j,i,iens) += -(w*dt       ) * gllWts_ngll(kk);
        }
        real dens_pert = state(idR,hs+k,hs+j,hs+i,iens);
        real dens      = state(idR,hs+k,hs+j,hs+i,iens) + hyDensCells(hs+k,iens);
        stateTend(idW,k,j,i,iens) += -dens_pert/dens * GRAV;

        // Left interface
        stateLimits(idR,1,k  ,j,i,iens) = r_gll(0     );
        stateLimits(idU,1,k  ,j,i,iens) = u_gll(0     );
        stateLimits(idV,1,k  ,j,i,iens) = v_gll(0     );
        stateLimits(idW,1,k  ,j,i,iens) = w_gll(0     );
        stateLimits(idT,1,k  ,j,i,iens) = t_gll(0     );
        // Right interface
        stateLimits(idR,0,k+1,j,i,iens) = r_gll(ngll-1);
        stateLimits(idU,0,k+1,j,i,iens) = u_gll(ngll-1);
        stateLimits(idV,0,k+1,j,i,iens) = v_gll(ngll-1);
        stateLimits(idW,0,k+1,j,i,iens) = w_gll(ngll-1);
        stateLimits(idT,0,k+1,j,i,iens) = t_gll(ngll-1);
      } // State

      { // Tracers
        for (int tr=0; tr < num_tracers; tr++) {
          SArray<real,1,ngll> t_gll, dt_gll;
          { // Reconstruction
            SArray<real,1,ord> stencil;
            for (int kk=0; kk < ord; kk++) { stencil(kk) = tracers(tr,k+kk,hs+j,hs+i,iens); }
            reconstruct_gll_values_and_derivs( stencil , t_gll , dt_gll , dz(k,iens) , c2g , c2d2g , s2g , s2d2g , wenoRecon ,
                                               idl , sigma , weno_scalars );
            if (tracer_pos(tr)) {
              for (int kk=0; kk < ngll; kk++) { t_gll(kk) = max( 0._fp , t_gll(kk) ); }
            }
          } // Reconstruction

          // Compute central integral terms
          tracerTend(tr,k,j,i,iens) = 0;
          for (int kk=0; kk < ngll; kk++) {
            real dt = dt_gll(kk);
            real w  = w_gll (kk);
            tracerTend(tr,k,j,i,iens) += -(w*dt) * gllWts_ngll(kk);
          }

          tracerLimits(tr,1,k  ,j,i,iens) = t_gll(0     );
          tracerLimits(tr,0,k+1,j,i,iens) = t_gll(ngll-1);
        }
      } // Tracers

    });

    ////////////////////////////////////////////////
    // BCs for the state edge estimates
    ////////////////////////////////////////////////
    parallel_for( SimpleBounds<3>(ny,nx,nens) , YAKL_LAMBDA (int j, int i, int iens) {
      for (int l = 0; l < num_state; l++) {
        if        (bc_z == BC_PERIODIC) {
          stateLimits     (l,0,0 ,j,i,iens) = stateLimits     (l,0,nz,j,i,iens);
          stateLimits     (l,1,nz,j,i,iens) = stateLimits     (l,1,0 ,j,i,iens);
        } else if (bc_z == BC_WALL    ) {
          stateLimits     (l,0,0 ,j,i,iens) = stateLimits     (l,1,0 ,j,i,iens);
          stateLimits     (l,1,nz,j,i,iens) = stateLimits     (l,0,nz,j,i,iens);
        }
      }
      for (int l = 0; l < num_tracers; l++) {
        if        (bc_z == BC_PERIODIC) {
          tracerLimits(l,0,0 ,j,i,iens) = tracerLimits(l,0,nz,j,i,iens);
          tracerLimits(l,1,nz,j,i,iens) = tracerLimits(l,1,0 ,j,i,iens);
        } else if (bc_z == BC_WALL    ) {
          tracerLimits(l,0,0 ,j,i,iens) = tracerLimits(l,1,0 ,j,i,iens);
          tracerLimits(l,1,nz,j,i,iens) = tracerLimits(l,0,nz,j,i,iens);
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the upwind fluxes
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<4>(nz+1,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      // Get left and right state
      real r_L = stateLimits(idR,0,k,j,i,iens);   real r_R = stateLimits(idR,1,k,j,i,iens);
      real u_L = stateLimits(idU,0,k,j,i,iens);   real u_R = stateLimits(idU,1,k,j,i,iens);
      real v_L = stateLimits(idV,0,k,j,i,iens);   real v_R = stateLimits(idV,1,k,j,i,iens);
      real w_L = stateLimits(idW,0,k,j,i,iens);   real w_R = stateLimits(idW,1,k,j,i,iens);
      real t_L = stateLimits(idT,0,k,j,i,iens);   real t_R = stateLimits(idT,1,k,j,i,iens);
      real p_L = C0 * pow( r_L*t_L , gamma )  ;   real p_R = C0 * pow( r_R*t_R , gamma );
      // Compute average state
      real r = 0.5_fp * (r_L + r_R);
      real w = 0.5_fp * (w_L + w_R);
      real t = 0.5_fp * (t_L + t_R);
      real p = 0.5_fp * (p_L + p_R);
      real cs2 = gamma*p/r;
      real cs  = sqrt(cs2);
      // Compute state difference across the cell interface
      real dr = r_R - r_L;
      real du = u_R - u_L;
      real dv = v_R - v_L;
      real dw = w_R - w_L;
      real dt = t_R - t_L;
      real dp = p_R - p_L;
      // Compute flux difference across the cell interface
      real df1 = w*dr + r*dw;
      real df2 = w*du;
      real df3 = w*dv;
      real df4 = w*dw + dp/r;
      real df5 = w*dt;
      // Compute characteristic variables
      real w1 = 0.5_fp * ( df1 - r/cs*df4 + r/t*df5 );
      real w2 = 0.5_fp * ( df1 + r/cs*df4 + r/t*df5 );
      real w3 = -r/t*df5;
      real w4 = df2;
      real w5 = df3;
      // Compute flux waves
      for (int l=0; l < num_state; l++) {
        stateLimits(l,0,k,j,i,iens) = 0;
        stateLimits(l,1,k,j,i,iens) = 0;
      }
      // WAVE 1: w-cs
      stateLimits(idR,0,k,j,i,iens) += w1;
      stateLimits(idW,0,k,j,i,iens) += -cs/r*w1;
      // WAVE 2: w+cs
      stateLimits(idR,1,k,j,i,iens) += w2;
      stateLimits(idW,1,k,j,i,iens) += cs/r*w2;
      // WAVES 3-5: w
      int upw_ind = 0;
      if (w > 0) {
        upw_ind = 1;
      }
      // WAVE 3
      stateLimits(idR,upw_ind,k,j,i,iens) += w3;
      stateLimits(idT,upw_ind,k,j,i,iens) += -t/r*w3;
      // WAVE 4
      stateLimits(idU,upw_ind,k,j,i,iens) += w4;
      // WAVE 5
      stateLimits(idV,upw_ind,k,j,i,iens) += w5;

      for (int tr=0; tr < num_tracers; tr++) {
        dt = tracerLimits(tr,1,k,j,i,iens) - tracerLimits(tr,0,k,j,i,iens);
        tracerLimits(tr,0,k,j,i,iens) = 0;
        tracerLimits(tr,1,k,j,i,iens) = 0;
        tracerLimits(tr,upw_ind,k,j,i,iens) += w*dt;
      }

      if (k == 0 || k == nz) {
        for (int l=0; l < num_state; l++) {
          stateLimits(l,0,k,j,i,iens) = 0;
          stateLimits(l,1,k,j,i,iens) = 0;
        }
        for (int l=0; l < num_tracers; l++) {
          tracerLimits(l,0,k,j,i,iens) = 0;
          tracerLimits(l,1,k,j,i,iens) = 0;
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Limit the tracer fluxes for positivity
    //////////////////////////////////////////////////////////
    real5d fct_mult("fct_mult",num_tracers,nz+1,ny,nx,nens);
    parallel_for( SimpleBounds<5>(num_tracers,nz+1,ny,nx,nens) , YAKL_LAMBDA (int tr, int k, int j, int i, int iens) {
      fct_mult(tr,k,j,i,iens) = 1.;
      // TODO:
    });

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA(int k, int j, int i, int iens) {
      for (int l = 0; l < num_state; l++) {
        if (sim2d && l == idV) {
          stateTend(l,k,j,i,iens) = 0;
        } else {
          stateTend(l,k,j,i,iens) += - ( stateLimits(l,0,k+1,j,i,iens) + stateLimits(l,1,k,j,i,iens) ) / dz(k,iens);
        }
      }
      for (int l = 0; l < num_tracers; l++) {
          tracerTend(l,k,j,i,iens) += - ( tracerLimits(l,0,k+1,j,i,iens) + tracerLimits(l,1,k,j,i,iens) ) / dz(k,iens);
      }
    });
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
        parallel_for( nz , YAKL_LAMBDA (int k) { zmid(k) = hyDensCells(hs+k,iens); });
        nc.write(zmid.createHostCopy(),"hyDens"    ,{"z"});

        parallel_for( nz , YAKL_LAMBDA (int k) { zmid(k) = hyPressureCells(hs+k,iens); });
        nc.write(zmid.createHostCopy(),"hyPressure",{"z"});

        parallel_for( nz , YAKL_LAMBDA (int k) { zmid(k) = hyThetaCells(hs+k,iens); });
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
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        data(k,j,i) = state(idR,hs+k,hs+j,hs+i,iens);
      });
      nc.write1(data.createHostCopy(),"dens_pert",{"z","y","x"},ulIndex,"t");
      // u
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        data(k,j,i) = state(idU,hs+k,hs+j,hs+i,iens);
      });
      nc.write1(data.createHostCopy(),"u",{"z","y","x"},ulIndex,"t");
      // v
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        data(k,j,i) = state(idV,hs+k,hs+j,hs+i,iens);
      });
      nc.write1(data.createHostCopy(),"v",{"z","y","x"},ulIndex,"t");
      // w
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        data(k,j,i) = state(idW,hs+k,hs+j,hs+i,iens);
      });
      nc.write1(data.createHostCopy(),"w",{"z","y","x"},ulIndex,"t");
      // theta'
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        data(k,j,i) = state(idT,hs+k,hs+j,hs+i,iens);
      });
      nc.write1(data.createHostCopy(),"pot_temp_pert",{"z","y","x"},ulIndex,"t");
      // pressure'
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        real r = state(idR,hs+k,hs+j,hs+i,iens) + hyDensCells (hs+k,iens);
        real t = state(idT,hs+k,hs+j,hs+i,iens) + hyThetaCells(hs+k,iens);
        real p  = C0*pow(r*t,gamma);
        data(k,j,i) = p - hyPressureCells(hs+k,iens);
      });
      nc.write1(data.createHostCopy(),"pressure_pert",{"z","y","x"},ulIndex,"t");

      for (int tr=0; tr < num_tracers; tr++) {
        parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
          real r = state(idR,hs+k,hs+j,hs+i,iens) + hyDensCells(hs+k,iens);
          data(k,j,i) = tracers(tr,hs+k,hs+j,hs+i,iens)*r;
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



  // ord stencil values to ngll GLL values; store in gll
  YAKL_INLINE
  void reconstruct_gll_values_and_derivs( SArray<real,1,ord> const stencil ,
                                          SArray<real,1,ngll> &gll , SArray<real,1,ngll> &deriv_gll , real dx ,
                                          SArray<real,2,ord,ngll> const &coefs_to_gll ,
                                          SArray<real,2,ord,ngll> const &coefs_to_deriv_to_gll ,
                                          SArray<real,2,ord,ngll> const &sten_to_gll  ,
                                          SArray<real,2,ord,ngll> const &sten_to_deriv_to_gll  ,
                                          SArray<real,3,ord,ord,ord> const &wenoRecon ,
                                          SArray<real,1,hs+2> const &idl              ,
                                          real sigma, bool doweno ) {
    if (doweno) {

      // Reconstruct values
      SArray<real,1,ord> wenoCoefs;
      weno::compute_weno_coefs<ord>( wenoRecon , stencil , wenoCoefs , idl , sigma );
      // Transform ord weno coefficients into ngll GLL points
      for (int ii=0; ii<ngll; ii++) {
        real val = 0;
        real der = 0;
        for (int s=0; s < ord; s++) {
          val += coefs_to_gll         (s,ii) * wenoCoefs(s);
          der += coefs_to_deriv_to_gll(s,ii) * wenoCoefs(s);
        }
        gll      (ii) = val;
        deriv_gll(ii) = der / dx;
      }

    } else {

      // Transform ord stencil cell averages into ngll GLL points
      for (int ii=0; ii<ngll; ii++) {
        real val = 0;
        real der = 0;
        for (int s=0; s < ord; s++) {
          val += sten_to_gll         (s,ii) * stencil(s);
          der += sten_to_deriv_to_gll(s,ii) * stencil(s);
        }
        gll      (ii) = val;
        deriv_gll(ii) = der / dx;
      }

    } // if doweno
  }



  // ord stencil values to ngll GLL values; store in gll
  YAKL_INLINE
  void reconstruct_gll_derivs( SArray<real,1,ord> const stencil ,
                               SArray<real,1,ngll> &deriv_gll , real dx ,
                               SArray<real,2,ord,ngll> const &coefs_to_deriv_to_gll ,
                               SArray<real,2,ord,ngll> const &sten_to_deriv_to_gll  ,
                               SArray<real,3,ord,ord,ord> const &wenoRecon ,
                               SArray<real,1,hs+2> const &idl              ,
                               real sigma, bool doweno ) {
    if (doweno) {

      // Reconstruct values
      SArray<real,1,ord> wenoCoefs;
      weno::compute_weno_coefs<ord>( wenoRecon , stencil , wenoCoefs , idl , sigma );
      // Transform ord weno coefficients into ngll GLL points
      for (int ii=0; ii<ngll; ii++) {
        real der = 0;
        for (int s=0; s < ord; s++) {
          der += coefs_to_deriv_to_gll(s,ii) * wenoCoefs(s);
        }
        deriv_gll(ii) = der / dx;
      }

    } else {

      // Transform ord stencil cell averages into ngll GLL points
      for (int ii=0; ii<ngll; ii++) {
        real der = 0;
        for (int s=0; s < ord; s++) {
          der += sten_to_deriv_to_gll(s,ii) * stencil(s);
        }
        deriv_gll(ii) = der / dx;
      }

    } // if doweno
  }



};
