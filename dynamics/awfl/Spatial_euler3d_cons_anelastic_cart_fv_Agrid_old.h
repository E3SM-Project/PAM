
#pragma once

#include "const.h"
#include "phys_params.h"
#include "TransformMatrices.h"
#include "WenoLimiter.h"
#include "Profiles.h"
#include "DataManager.h"


template <int nTimeDerivs, bool timeAvg, int nAder>
class Spatial_euler3d_cons_expl_cart_fv_Agrid {
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

  int tracer_id_uniform;

  typedef real4d StateArr;  // Array of state variables (rho, rho*u, rho*v, rho*w, and rho*theta)
  typedef real4d TracerArr; // Array of tracers (total tracer mass)

  typedef real4d StateTendArr;   // State tendencies
  typedef real4d TracerTendArr;  // Tracer tendencies

  // Stores two estimates of the state, state flux, and tracer values at each cell interface
  real5d stateLimits;
  real4d stateFlux;

  // Hydrostatically balanced values for density, potential temperature, and pressure (cell-averages)
  real1d hyDensCells;
  real1d hyPressureCells;
  real1d hyThetaCells;
  real1d hyDensThetaCells;

  // Hydrostatically balanced values for density, potential temperature, and pressure (GLL points)
  real2d hyDensGLL;
  real2d hyPressureGLL;
  real2d hyThetaGLL;
  real2d hyDensThetaGLL;

  real3d pressure;
  real3d pressure_flux_x;
  real3d pressure_flux_y;
  real3d pressure_flux_z;
  real4d pressure_limits_x;
  real4d pressure_limits_y;
  real4d pressure_limits_z;

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
  int static constexpr DATA_SPEC_THERMAL_MOIST = 2;
  
  bool sim2d;  // Whether we're simulating in 2-D

  // Grid spacing in each dimension
  real dx;
  real dy;
  real dz;

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
  // Number of cells in each direction
  int         nx;
  int         ny;
  int         nz;
  // Length of the domain in each direction (m)
  real        xlen;
  real        ylen;
  real        zlen;
  // Boundary condition in each direction
  int         bc_x;
  int         bc_y;
  int         bc_z;
  // Whether to use WENO for scalars and also for winds
  bool        weno_scalars;
  bool        weno_winds;
  // Name of the output file
  std::string out_file;
  // How to initialize the data
  int         data_spec;


  // When this class is created, initialize num_tracers to zero
  Spatial_euler3d_cons_expl_cart_fv_Agrid() {
    num_tracers = 0;
  }



  // Make sure it's odd-order-accurate
  static_assert(ord%2 == 1,"ERROR: ord must be an odd integer");



  // Initialize a tracer
  int add_tracer(DataManager &dm , std::string name , std::string desc , bool pos_def , bool adds_mass) {
    auto &tracer_pos       = this->tracer_pos;
    auto &tracer_adds_mass = this->tracer_adds_mass;

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

    // Register and allocate this tracer in the DataManager
    dm.register_and_allocate<real>( name , desc , {nz,ny,nx} , {"z","y","x"} );

    // Return the index of this tracer to the caller
    return tr;
  }



  // Caller creates a lambda (init_mass) to initialize this tracer value using location and dry state information
  template <class MICRO>
  void init_tracers( DataManager &dm , MICRO const &micro) {
    tracer_id_uniform = add_tracer(dm , "uniform"  , "uniform"  , false     , false);
  }



  // Transform state and tracer data in DataManager into a state and tracers array more conveniently used by the dycore
  // This has to be copied because we need a halo, and the rest of the model doesn't need to know about the halo
  void read_state_and_tracers( DataManager &dm , real4d &state , real4d &tracers) const {
    auto &num_tracers = this->num_tracers;

    // Get data arrays from the DataManager (this just wraps an existing allocated pointer in an Array)
    real3d rho          = dm.get<real,3>( "density" );
    real3d rho_u        = dm.get<real,3>( "density_u" );
    real3d rho_v        = dm.get<real,3>( "density_v" );
    real3d rho_w        = dm.get<real,3>( "density_w" );
    real3d rho_theta    = dm.get<real,3>( "density_theta" );
    real1d rho_hy       = dm.get<real,1>( "hydrostatic_density" );
    real1d rho_theta_hy = dm.get<real,1>( "hydrostatic_density_theta" );

    // An array of tracers for reading in tracers from the DataManager
    MultipleTracers<max_tracers> dm_tracers;
    for (int tr=0; tr < num_tracers; tr++) {
      real3d tracer = dm.get<real,3>( tracer_name[tr] );
      dm_tracers.add_tracer( tracer );
    }

    // Copy from the DataManager to the state and tracers arrays
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      state(idR,hs+k,hs+j,hs+i) = rho(k,j,i) - rho_hy(k);  // Density perturbation
      state(idU,hs+k,hs+j,hs+i) = rho_u(k,j,i);
      state(idV,hs+k,hs+j,hs+i) = rho_v(k,j,i);
      state(idW,hs+k,hs+j,hs+i) = rho_w(k,j,i);
      state(idT,hs+k,hs+j,hs+i) = rho_theta(k,j,i) - rho_theta_hy(k);  // rho*theta perturbation
      for (int tr=0; tr < num_tracers; tr++) {
        tracers(tr,hs+k,hs+j,hs+i) = dm_tracers(tr,k,j,i);
      }
    });
  }



  // Transform state and tracer data from state and tracers arrays with halos back to the DataManager
  // so it can be used by other parts of the model
  void write_state_and_tracers( DataManager &dm , real4d &state , real4d &tracers) const {
    auto &num_tracers      = this->num_tracers     ;
    auto &hyDensCells      = this->hyDensCells     ;
    auto &hyDensThetaCells = this->hyDensThetaCells;

    // Get data arrays from the DataManager (this just wraps an existing allocated pointer in an Array)
    real3d rho       = dm.get<real,3>( "density" );
    real3d rho_u     = dm.get<real,3>( "density_u" );
    real3d rho_v     = dm.get<real,3>( "density_v" );
    real3d rho_w     = dm.get<real,3>( "density_w" );
    real3d rho_theta = dm.get<real,3>( "density_theta" );

    // An array of tracers for reading in tracers from the DataManager
    MultipleTracers<max_tracers> dm_tracers;
    for (int tr=0; tr < num_tracers; tr++) {
      real3d tracer = dm.get<real,3>( tracer_name[tr] );
      dm_tracers.add_tracer( tracer );
    }

    // Copy from state and tracers arrays to the DataManager arrays
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      rho      (k,j,i) = state(idR,hs+k,hs+j,hs+i) + hyDensCells(hs+k)     ;
      rho_u    (k,j,i) = state(idU,hs+k,hs+j,hs+i)                         ;
      rho_v    (k,j,i) = state(idV,hs+k,hs+j,hs+i)                         ;
      rho_w    (k,j,i) = state(idW,hs+k,hs+j,hs+i)                         ;
      rho_theta(k,j,i) = state(idT,hs+k,hs+j,hs+i) + hyDensThetaCells(hs+k);
      for (int tr=0; tr < num_tracers; tr++) {
        dm_tracers(tr,k,j,i) = tracers(tr,hs+k,hs+j,hs+i);
      }
    });
  }



  // Take an initially dry fluid state and adjust it to account for moist tracers
  template <class MICRO>
  void adjust_state_for_moisture(DataManager &dm , MICRO const &micro) const {
  }



  real4d createStateArr() const {
    return real4d("stateArr",num_state,nz+2*hs,ny+2*hs,nx+2*hs);
  }



  real4d createTracerArr() const {
    return real4d("tracerArr",num_tracers,nz+2*hs,ny+2*hs,nx+2*hs);
  }



  real4d createStateTendArr() const {
    return real4d("stateTendArr",num_state,nz,ny,nx);
  }



  real4d createTracerTendArr() const {
    return real4d("tracerTendArr",num_tracers,nz,ny,nx);
  }



  // Number of operator splittinng steps to use
  // Normally this would be 3, but the z-directly CFL is reduced because of how the fluxes are
  // handled in the presence of a solid wall boundary condition. I'm looking into how to fix this
  int numSplit() const {
    return 4;
  }



  // Given the model data and CFL value, compute the maximum stable time step
  template <class MICRO>
  real compute_time_step(real cfl, DataManager &dm, MICRO const &micro) {

    // If we've already computed the time step, then don't compute it again
    if (dtInit <= 0) {
      auto &dx                   = this->dx                  ;
      auto &dy                   = this->dy                  ;
      auto &dz                   = this->dz                  ;
      auto &hyDensCells          = this->hyDensCells         ;
      auto &hyDensThetaCells     = this->hyDensThetaCells    ;
      auto &gamma                = this->gamma               ;
      auto &C0                   = this->C0                  ;

      // Convert data from DataManager to state and tracers array for convenience
      real4d state   = createStateArr ();
      real4d tracers = createTracerArr();
      read_state_and_tracers( dm , state , tracers );

      // Allocate a 3-D array for the max stable time steps (we'll use this for a reduction later)
      real3d dt3d("dt3d",nz,ny,nx);

      // Loop through the cells, calculate the max stable time step for each cell
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        // Get the state
        real r = hyDensCells(hs+k);
        real u = state(idU,hs+k,hs+j,hs+i) / r;
        real v = state(idV,hs+k,hs+j,hs+i) / r;
        real w = state(idW,hs+k,hs+j,hs+i) / r;

        // Compute the maximum stable time step in each direction
        real udt = cfl * dx / max(50._fp , abs(u) );
        real vdt = cfl * dy / max(50._fp , abs(v) );
        real wdt = cfl * dz / max(50._fp , abs(w) );

        // Compute the min of the max stable time steps
        dt3d(k,j,i) = min( min(udt,vdt) , wdt );
      });

      // Reduce the max stable time steps to the minimum of all cells
      yakl::ParallelMin<real,memDevice> pmin(nz*nx*ny);

      // Store to dtInit so we don't have to compute this again
      dtInit = pmin(dt3d.data());
    }
    
    return dtInit;
  }



  // Initialize crap needed by recon()
  void init(std::string inFile, int num_tracers, DataManager &dm) {
    this->num_tracers = num_tracers;
    
    // Allocate device arrays for whether tracers are positive-definite or add mass
    tracer_pos       = bool1d("tracer_pos"      ,num_tracers);
    tracer_adds_mass = bool1d("tracer_adds_mass",num_tracers);

    // Inialize time step to zero, and dimensional splitting switch
    dtInit = 0;
    dimSwitch = true;

    // Read the YAML input file
    YAML::Node config = YAML::LoadFile(inFile);
    if ( !config                 ) { endrun("ERROR: Invalid YAML input file"); }
    if ( !config["nx"]           ) { endrun("ERROR: No nx in input file"); }
    if ( !config["ny"]           ) { endrun("ERROR: No ny in input file"); }
    if ( !config["nz"]           ) { endrun("ERROR: No nz in input file"); }
    if ( !config["xlen"]         ) { endrun("ERROR: No xlen in input file"); }
    if ( !config["ylen"]         ) { endrun("ERROR: No ylen in input file"); }
    if ( !config["zlen"]         ) { endrun("ERROR: No zlen in input file"); }
    if ( !config["bc_x"]         ) { endrun("ERROR: No bc_x in input file"); }
    if ( !config["bc_y"]         ) { endrun("ERROR: No bc_y in input file"); }
    if ( !config["bc_z"]         ) { endrun("ERROR: No bc_z in input file"); }
    if ( !config["weno_scalars"] ) { endrun("ERROR: No weno_scalars in input file"); }
    if ( !config["weno_winds"]   ) { endrun("ERROR: No weno_winds in input file"); }
    if ( !config["initData"]     ) { endrun("ERROR: No initData in input file"); }
    if ( !config["out_file"]     ) { endrun("ERROR: No out_file in input file"); }

    // Read the # cells in each dimension
    nx = config["nx"].as<int>();
    ny = config["ny"].as<int>();
    nz = config["nz"].as<int>();

    // Determine whether this is a 2-D simulation
    sim2d = ny == 1;

    // Read the domain length in each dimension
    xlen = config["xlen"].as<real>();
    ylen = config["ylen"].as<real>();
    zlen = config["zlen"].as<real>();

    // Read whether we're doing WENO limiting on scalars and winds
    weno_scalars = config["weno_scalars"].as<bool>();
    weno_winds   = config["weno_winds"].as<bool>();

    // Read the data initialization option
    std::string dataStr = config["initData"].as<std::string>();
    if        (dataStr == "thermal") {
      data_spec = DATA_SPEC_THERMAL;
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
    out_file = config["out_file"].as<std::string>();

    // Compute the grid spacing in each dimension
    dx = xlen/nx;
    dy = ylen/ny;
    dz = zlen/nz;

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
    weno::wenoSetIdealSigma(this->idl,this->sigma);

    // Allocate data
    stateLimits     = real5d("stateLimits"    ,num_state  ,2,nz+1,ny+1,nx+1);
    stateFlux       = real4d("stateFlux"      ,num_state    ,nz+1,ny+1,nx+1);
    hyDensCells          = real1d("hyDensCells       ",nz+2*hs);
    hyPressureCells      = real1d("hyPressureCells   ",nz+2*hs);
    hyThetaCells         = real1d("hyThetaCells      ",nz+2*hs);
    hyDensThetaCells     = real1d("hyDensThetaCells  ",nz+2*hs);
    hyDensGLL            = real2d("hyDensGLL         ",nz,ngll);
    hyPressureGLL        = real2d("hyPressureGLL     ",nz,ngll);
    hyThetaGLL           = real2d("hyThetaGLL        ",nz,ngll);
    hyDensThetaGLL       = real2d("hyDensThetaGLL    ",nz,ngll);
    pressure = real3d("pressure",nz+2*hs,ny+2*hs,nx+2*hs);
    pressure_flux_x = real3d("pressure_flux_x",nz,ny,nx+1);
    pressure_flux_y = real3d("pressure_flux_y",nz,ny+1,nx);
    pressure_flux_z = real3d("pressure_flux_z",nz+1,ny,nx);
    pressure_limits_x = real4d("pressure_limits_x",2,nz,ny,nx+1);
    pressure_limits_y = real4d("pressure_limits_y",2,nz,ny+1,nx);
    pressure_limits_z = real4d("pressure_limits_z",2,nz+1,ny,nx);

    // Register and allocate state data with the DataManager
    dm.register_and_allocate<real>( "hydrostatic_density"       , "hydrostatic_density"       , {nz} , {"z"} );
    dm.register_and_allocate<real>( "hydrostatic_theta"         , "hydrostatic_theta"         , {nz} , {"z"} );
    dm.register_and_allocate<real>( "hydrostatic_density_theta" , "hydrostatic_density_theta" , {nz} , {"z"} );
    dm.register_and_allocate<real>( "hydrostatic_pressure"      , "hydrostatic_pressure"      , {nz} , {"z"} );
    dm.register_and_allocate<real>( "density"       , "density"       , {nz,ny,nx} , {"z","y","x"} );
    dm.register_and_allocate<real>( "density_u"     , "density_u"     , {nz,ny,nx} , {"z","y","x"} );
    dm.register_and_allocate<real>( "density_v"     , "density_v"     , {nz,ny,nx} , {"z","y","x"} );
    dm.register_and_allocate<real>( "density_w"     , "density_w"     , {nz,ny,nx} , {"z","y","x"} );
    dm.register_and_allocate<real>( "density_theta" , "density_theta" , {nz,ny,nx} , {"z","y","x"} );
  }



  // Initialize the state
  template <class MICRO>
  void init_state( DataManager &dm , MICRO const &micro ) {
    Rd    = micro.constants.R_d;
    cp    = micro.constants.cp_d;
    gamma = micro.constants.gamma_d;
    p0    = micro.constants.p0;
    C0    = micro.constants.C0_d;

    auto nx                       = this->nx                     ;
    auto ny                       = this->ny                     ;
    auto nz                       = this->nz                     ;
    auto dx                       = this->dx                     ;
    auto dy                       = this->dy                     ;
    auto dz                       = this->dz                     ;
    auto gllPts_ord               = this->gllPts_ord             ;
    auto gllWts_ord               = this->gllWts_ord             ;
    auto gllPts_ngll              = this->gllPts_ngll            ;
    auto gllWts_ngll              = this->gllWts_ngll            ;
    auto &hyDensCells             = this->hyDensCells            ;
    auto &hyThetaCells            = this->hyThetaCells           ;
    auto &hyPressureCells         = this->hyPressureCells        ;
    auto &hyDensThetaCells        = this->hyDensThetaCells       ;
    auto &hyDensGLL               = this->hyDensGLL              ;
    auto &hyThetaGLL              = this->hyThetaGLL             ;
    auto &hyPressureGLL           = this->hyPressureGLL          ;
    auto &hyDensThetaGLL          = this->hyDensThetaGLL         ;
    auto &data_spec               = this->data_spec              ;
    auto &sim2d                   = this->sim2d                  ;
    auto &xlen                    = this->xlen                   ;
    auto &ylen                    = this->ylen                   ;
    auto &Rd                      = this->Rd                     ;
    auto &cp                      = this->cp                     ;
    auto &gamma                   = this->gamma                  ;
    auto &p0                      = this->p0                     ;
    auto &C0                      = this->C0                     ;

    // Get Arrays for 1-D hydrostatic background profiles
    real1d dm_hyDens      = dm.get<real,1>( "hydrostatic_density"       );
    real1d dm_hyTheta     = dm.get<real,1>( "hydrostatic_theta"         );
    real1d dm_hyDensTheta = dm.get<real,1>( "hydrostatic_density_theta" );
    real1d dm_hyPressure  = dm.get<real,1>( "hydrostatic_pressure"      );

    // Setup hydrostatic background state
    parallel_for( SimpleBounds<1>(nz+2*hs) , YAKL_LAMBDA (int k) {
      // Compute cell averages
      hyDensCells     (k) = 0;
      hyPressureCells (k) = 0;
      hyThetaCells    (k) = 0;
      hyDensThetaCells(k) = 0;
      for (int kk=0; kk<ord; kk++) {
        real zloc = (k-hs+0.5_fp)*dz + gllPts_ord(kk)*dz;
        if        (data_spec == DATA_SPEC_THERMAL || data_spec == DATA_SPEC_THERMAL_MOIST) {
          // Compute constant theta hydrostatic background state
          real th  = 300;
          real rh = profiles::initConstTheta_density (th,zloc,Rd,cp,gamma,p0,C0);
          real ph = profiles::initConstTheta_pressure(th,zloc,Rd,cp,gamma,p0,C0);
          real wt = gllWts_ord(kk);
          hyDensCells     (k) += rh    * wt;
          hyThetaCells    (k) += th    * wt;
          hyDensThetaCells(k) += rh*th * wt;
          hyPressureCells (k) += ph    * wt;
        }
      }
      if (k >= hs && k <= hs+nz-1) {
        dm_hyDens     (k-hs) = hyDensCells     (k);
        dm_hyTheta    (k-hs) = hyThetaCells    (k);
        dm_hyDensTheta(k-hs) = hyDensThetaCells(k);
        dm_hyPressure (k-hs) = hyPressureCells (k);
      }
    });

    parallel_for( SimpleBounds<1>(nz) , YAKL_LAMBDA (int k) {
      // Compute ngll GLL points
      for (int kk=0; kk<ngll; kk++) {
        real zloc = (k+0.5_fp)*dz + gllPts_ngll(kk)*dz;
        if        (data_spec == DATA_SPEC_THERMAL || data_spec == DATA_SPEC_THERMAL_MOIST) {
          // Compute constant theta hydrostatic background state
          real th = 300;
          real rh = profiles::initConstTheta_density (th,zloc,Rd,cp,gamma,p0,C0);
          real ph = profiles::initConstTheta_pressure(th,zloc,Rd,cp,gamma,p0,C0);
          hyDensGLL     (k,kk) = rh;
          hyThetaGLL    (k,kk) = th;
          hyDensThetaGLL(k,kk) = rh*th;
          hyPressureGLL (k,kk) = ph;
        }
      }
    });
    
    real3d dm_rho       = dm.get<real,3>( "density"       );
    real3d dm_rho_u     = dm.get<real,3>( "density_u"     );
    real3d dm_rho_v     = dm.get<real,3>( "density_v"     );
    real3d dm_rho_w     = dm.get<real,3>( "density_w"     );
    real3d dm_rho_theta = dm.get<real,3>( "density_theta" );

    // Compute the state
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      dm_rho      (k,j,i) = 0;
      dm_rho_u    (k,j,i) = 0;
      dm_rho_v    (k,j,i) = 0;
      dm_rho_w    (k,j,i) = 0;
      dm_rho_theta(k,j,i) = 0;
      for (int kk=0; kk<ord; kk++) {
        for (int jj=0; jj<ord; jj++) {
          for (int ii=0; ii<ord; ii++) {
            real zloc = (k+0.5_fp)*dz + gllPts_ord(kk)*dz;
            real yloc;
            if (sim2d) {
              yloc = ylen/2;
            } else {
              yloc = (j+0.5_fp)*dy + gllPts_ord(jj)*dy;
            }
            real xloc = (i+0.5_fp)*dx + gllPts_ord(ii)*dx;
            real wt = gllWts_ord(kk) * gllWts_ord(jj) * gllWts_ord(ii);
            if        (data_spec == DATA_SPEC_THERMAL || data_spec == DATA_SPEC_THERMAL_MOIST) {
              // Compute constant theta hydrostatic background state
              real th = 300;
              real rh = profiles::initConstTheta_density(th,zloc,Rd,cp,gamma,p0,C0);
              real tp = profiles::ellipsoid_linear(xloc, yloc, zloc, xlen/2, ylen/2, 2000, 2000, 2000, 2000, 2 );
              real t = th + tp;
              real r = rh;

              dm_rho      (k,j,i) += (r - rh)*wt;
              dm_rho_theta(k,j,i) += (r*t - rh*th) * wt;
            }
          }
        }
      }
      dm_rho      (k,j,i) += dm_hyDens     (k);
      dm_rho_theta(k,j,i) += dm_hyDensTheta(k);
    });
  }



  // Compute state and tendency time derivatives from the state
  template <class MICRO>
  void computeTendencies( real4d &state   , real4d &stateTend  ,
                          real4d &tracers , real4d &tracerTend ,
                          MICRO const &micro, real &dt , int splitIndex ) {
    memset(stateTend,0._fp);
    memset(tracerTend,0._fp);

    if (splitIndex == 0) {
      computePressureTendencies( state , stateTend , dt );
    }
    if (splitIndex == 1) {
      auto &dx    = this->dx   ;
      auto &dy    = this->dy   ;
      auto &dz    = this->dz   ;
      auto &nx    = this->nx   ;
      auto &ny    = this->ny   ;
      auto &nz    = this->nz   ;
      auto &sim2d = this->sim2d;
      if        (bc_x == BC_PERIODIC) {
        parallel_for( SimpleBounds<3>(nz,ny,hs) , YAKL_LAMBDA(int k, int j, int ii) {
          state(idU,hs+k,hs+j,      ii) = state(idU,hs+k,hs+j,nx+ii);
          state(idU,hs+k,hs+j,hs+nx+ii) = state(idU,hs+k,hs+j,hs+ii);
        });
      } else if (bc_x == BC_WALL) {
        parallel_for( SimpleBounds<3>(nz,ny,hs) , YAKL_LAMBDA(int k, int j, int ii) {
          state(idU,hs+k,hs+j,      ii) = 0;
          state(idU,hs+k,hs+j,hs+nx+ii) = 0;
        });
      }
      if (!sim2d) {
        if        (bc_y == BC_PERIODIC) {
          parallel_for( SimpleBounds<3>(nz,nx,hs) , YAKL_LAMBDA(int k, int i, int jj) {
            state(idV,hs+k,      jj,hs+i) = state(idV,hs+k,ny+jj,hs+i);
            state(idV,hs+k,hs+ny+jj,hs+i) = state(idV,hs+k,hs+jj,hs+i);
          });
        } else if (bc_y == BC_WALL) {
          parallel_for( SimpleBounds<3>(nz,nx,hs) , YAKL_LAMBDA(int k, int i, int jj) {
            state(idV,hs+k,      jj,hs+i) = 0;
            state(idV,hs+k,hs+ny+jj,hs+i) = 0;
          });
        }
      }
      if        (bc_z == BC_PERIODIC) {
        parallel_for( SimpleBounds<3>(ny,nx,hs) , YAKL_LAMBDA(int j, int i, int kk) {
          state(idW,      kk,hs+j,hs+i) = state(idW,nz+kk,hs+j,hs+i);
          state(idW,hs+nz+kk,hs+j,hs+i) = state(idW,hs+kk,hs+j,hs+i);
        });
      } else if (bc_z == BC_WALL) {
        parallel_for( SimpleBounds<3>(ny,nx,hs) , YAKL_LAMBDA(int j, int i, int kk) {
          state(idW,      kk,hs+j,hs+i) = 0;
          state(idW,hs+nz+kk,hs+j,hs+i) = 0;
        });
      }
      real3d mom_div("mom_div",nz,ny,nx);
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        mom_div(k,j,i)  = (state(idU,hs+k,hs+j,hs+i+1) - state(idU,hs+k,hs+j,hs+i-1)) / (2*dx);
        if (!sim2d) {
          mom_div(k,j,i) += (state(idV,hs+k,hs+j+1,hs+i) - state(idV,hs+k,hs+j-1,hs+i)) / (2*dy);
        }
        mom_div(k,j,i) += (state(idW,hs+k+1,hs+j,hs+i) - state(idW,hs+k-1,hs+j,hs+i)) / (2*dz);
      });
      yakl::SimpleNetCDF nc;
      nc.open(out_file,yakl::NETCDF_MODE_WRITE);
      int ulIndex = nc.getDimSize("t");
      nc.write1(mom_div.createHostCopy(),"mom_div",{"z","y","x"},max(0,ulIndex-1),"t");
      nc.close();
    }

    if (dimSwitch) {
      if        (splitIndex == 1) {
        computeTendenciesX( state , stateTend , tracers , tracerTend , micro , dt );
      } else if (splitIndex == 2) {
        if (sim2d) {
          memset(stateTend  , 0._fp);
          memset(tracerTend , 0._fp);
        } else {
          computeTendenciesY( state , stateTend , tracers , tracerTend , micro , dt );
        }
      } else if (splitIndex == 3) {
        computeTendenciesZ( state , stateTend , tracers , tracerTend , micro , dt );
      }
    } else {
      if        (splitIndex == 1) {
        computeTendenciesZ( state , stateTend , tracers , tracerTend , micro , dt );
      } else if (splitIndex == 2) {
        if (sim2d) {
          memset(stateTend  , 0._fp);
          memset(tracerTend , 0._fp);
        } else {
          computeTendenciesY( state , stateTend , tracers , tracerTend , micro , dt );
        }
      } else if (splitIndex == 3) {
        computeTendenciesX( state , stateTend , tracers , tracerTend , micro , dt );
      }
    }
    if (splitIndex == numSplit()-1) dimSwitch = ! dimSwitch;
  } // computeTendencies



  void computePressureTendencies( real4d &state , real4d &stateTend , real dt ) {
    auto &nx                      = this->nx                     ;
    auto &ny                      = this->ny                     ;
    auto &nz                      = this->nz                     ;
    auto &weno_scalars            = this->weno_scalars           ;
    auto &weno_winds              = this->weno_winds             ;
    auto &c2g                     = this->coefs_to_gll           ;
    auto &s2g                     = this->sten_to_gll            ;
    auto &wenoRecon               = this->wenoRecon              ;
    auto &idl                     = this->idl                    ;
    auto &sigma                   = this->sigma                  ;
    auto &hyDensCells             = this->hyDensCells            ;
    auto &hyDensThetaCells        = this->hyDensThetaCells       ;
    auto &hyThetaCells            = this->hyThetaCells           ;
    auto &sim2d                   = this->sim2d                  ;
    auto &derivMatrix             = this->derivMatrix            ;
    auto &dx                      = this->dx                     ;
    auto &dy                      = this->dy                     ;
    auto &dz                      = this->dz                     ;
    auto &stateLimits             = this->stateLimits            ;
    auto &stateFlux               = this->stateFlux              ;
    auto &bc_x                    = this->bc_x                   ;
    auto &bc_y                    = this->bc_y                   ;
    auto &bc_z                    = this->bc_z                   ;
    auto &Rd                      = this->Rd                     ;
    auto &cp                      = this->cp                     ;
    auto &gamma                   = this->gamma                  ;
    auto &p0                      = this->p0                     ;
    auto &C0                      = this->C0                     ;
    auto &pressure                = this->pressure               ;
    auto &pressure_limits_x       = this->pressure_limits_x      ;
    auto &pressure_limits_y       = this->pressure_limits_y      ;
    auto &pressure_limits_z       = this->pressure_limits_z      ;
    auto &pressure_flux_x         = this->pressure_flux_x        ;
    auto &pressure_flux_y         = this->pressure_flux_y        ;
    auto &pressure_flux_z         = this->pressure_flux_z        ;

    real3d rhs("rhs",nz,ny,nx);
    real3d pressure_new("pressure_new",nz,ny,nx);

    // Populate the halos
    if        (bc_x == BC_PERIODIC) {
      parallel_for( SimpleBounds<3>(nz,ny,hs) , YAKL_LAMBDA(int k, int j, int ii) {
        state(idU,hs+k,hs+j,      ii) = state(idU,hs+k,hs+j,nx+ii);
        state(idU,hs+k,hs+j,hs+nx+ii) = state(idU,hs+k,hs+j,hs+ii);
      });
    } else if (bc_x == BC_WALL) {
      parallel_for( SimpleBounds<3>(nz,ny,hs) , YAKL_LAMBDA(int k, int j, int ii) {
        state(idU,hs+k,hs+j,      ii) = 0;
        state(idU,hs+k,hs+j,hs+nx+ii) = 0;
      });
    }
    if (!sim2d) {
      if        (bc_y == BC_PERIODIC) {
        parallel_for( SimpleBounds<3>(nz,nx,hs) , YAKL_LAMBDA(int k, int i, int jj) {
          state(idV,hs+k,      jj,hs+i) = state(idV,hs+k,ny+jj,hs+i);
          state(idV,hs+k,hs+ny+jj,hs+i) = state(idV,hs+k,hs+jj,hs+i);
        });
      } else if (bc_y == BC_WALL) {
        parallel_for( SimpleBounds<3>(nz,nx,hs) , YAKL_LAMBDA(int k, int i, int jj) {
          state(idV,hs+k,      jj,hs+i) = 0;
          state(idV,hs+k,hs+ny+jj,hs+i) = 0;
        });
      }
    }
    if        (bc_z == BC_PERIODIC) {
      parallel_for( SimpleBounds<3>(ny,nx,hs) , YAKL_LAMBDA(int j, int i, int kk) {
        state(idW,      kk,hs+j,hs+i) = state(idW,nz+kk,hs+j,hs+i);
        state(idW,hs+nz+kk,hs+j,hs+i) = state(idW,hs+kk,hs+j,hs+i);
      });
    } else if (bc_z == BC_WALL) {
      parallel_for( SimpleBounds<3>(ny,nx,hs) , YAKL_LAMBDA(int j, int i, int kk) {
        state(idW,      kk,hs+j,hs+i) = 0;
        state(idW,hs+nz+kk,hs+j,hs+i) = 0;
      });
    }

    ///////////////////////////////////////////
    // Compute RHS
    ///////////////////////////////////////////
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
                  rhs(k,j,i) =  -( state(idU,hs+k,hs+j,hs+i+1) - state(idU,hs+k,hs+j,hs+i-1) ) * (2*dx) / dt;
      if (!sim2d) rhs(k,j,i) += -( state(idV,hs+k,hs+j+1,hs+i) - state(idV,hs+k,hs+j-1,hs+i) ) * (2*dy) / dt;
                  rhs(k,j,i) += -( state(idW,hs+k+1,hs+j,hs+i) - state(idW,hs+k-1,hs+j,hs+i) ) * (2*dz) / dt;
    });

    // Compute pressure perturbation
    memset(pressure,0._fp);
    real dx2 = dx*dx;
    real dy2 = dy*dy;
    real dz2 = dz*dz;
    int nIter = 40000;
    for (int iter = 0; iter < nIter; iter++) {
      // Compute new pressure
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        real factor;
        if (sim2d) {
          factor = 1._fp / 4._fp;
        } else {
          factor = 1._fp / 6._fp;
        }
        real xterm = pressure(hs+k,hs+j,hs+i+2) + pressure(hs+k,hs+j,hs+i-2);
        real yterm = pressure(hs+k,hs+j+2,hs+i) + pressure(hs+k,hs+j-2,hs+i);
        real zterm = pressure(hs+k+2,hs+j,hs+i) + pressure(hs+k-2,hs+j,hs+i);
        if (sim2d) {
          pressure_new(k,j,i) = factor * ( xterm         + zterm + rhs(k,j,i) );
        } else {
          pressure_new(k,j,i) = factor * ( xterm + yterm + zterm + rhs(k,j,i) );
        }
      });

      // Assign new pressure
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        pressure(hs+k,hs+j,hs+i) = pressure_new(k,j,i);
      });

      // Pressure boundary conditions
      if        (bc_x == BC_PERIODIC) {
        parallel_for( SimpleBounds<3>(nz,ny,hs) , YAKL_LAMBDA(int k, int j, int ii) {
          pressure(hs+k,hs+j,      ii) = pressure(hs+k,hs+j,nx+ii);
          pressure(hs+k,hs+j,hs+nx+ii) = pressure(hs+k,hs+j,hs+ii);
        });
      } else if (bc_x == BC_WALL) {
        parallel_for( SimpleBounds<3>(nz,ny,hs) , YAKL_LAMBDA(int k, int j, int ii) {
          pressure(hs+k,hs+j,      ii) = pressure(hs+k,hs+j,hs+0   );
          pressure(hs+k,hs+j,hs+nx+ii) = pressure(hs+k,hs+j,hs+nx-1);
        });
      }
      if (!sim2d) {
        if        (bc_y == BC_PERIODIC) {
          parallel_for( SimpleBounds<3>(nz,nx,hs) , YAKL_LAMBDA(int k, int i, int jj) {
            pressure(hs+k,      jj,hs+i) = pressure(hs+k,ny+jj,hs+i);
            pressure(hs+k,hs+ny+jj,hs+i) = pressure(hs+k,hs+jj,hs+i);
          });
        } else if (bc_y == BC_WALL) {
          parallel_for( SimpleBounds<3>(nz,nx,hs) , YAKL_LAMBDA(int k, int i, int jj) {
            pressure(hs+k,      jj,hs+i) = pressure(hs+k,hs+0   ,hs+i);
            pressure(hs+k,hs+ny+jj,hs+i) = pressure(hs+k,hs+ny-1,hs+i);
          });
        }
      }
      if        (bc_z == BC_PERIODIC) {
        parallel_for( SimpleBounds<3>(ny,nx,hs) , YAKL_LAMBDA(int j, int i, int kk) {
          pressure(      kk,hs+j,hs+i) = pressure(nz+kk,hs+j,hs+i);
          pressure(hs+nz+kk,hs+j,hs+i) = pressure(hs+kk,hs+j,hs+i);
        });
      } else if (bc_z == BC_WALL) {
        parallel_for( SimpleBounds<3>(ny,nx,hs) , YAKL_LAMBDA(int j, int i, int kk) {
          pressure(      kk,hs+j,hs+i) = pressure(hs+0   ,hs+j,hs+i);
          pressure(hs+nz+kk,hs+j,hs+i) = pressure(hs+nz-1,hs+j,hs+i);
        });
      }
    }

    // Apply pressure gradient and gravity source term
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      stateTend(idU,k,j,i) = - (pressure(hs+k,hs+j,hs+i+1) - pressure(hs+k,hs+j,hs+i-1)) / (2*dx);
      if (!sim2d) {
        stateTend(idV,k,j,i) = - (pressure(hs+k,hs+j+1,hs+i) - pressure(hs+k,hs+j-1,hs+i)) / (2*dy);
      }
      stateTend(idW,k,j,i) = - (pressure(hs+k+1,hs+j,hs+i) - pressure(hs+k-1,hs+j,hs+i)) / (2*dz);
    });
  }



  template <class MICRO>
  void computeTendenciesX( real4d &state   , real4d &stateTend  ,
                           real4d &tracers , real4d &tracerTend ,
                           MICRO const &micro, real &dt ) {
    auto &nx                      = this->nx                     ;
    auto &weno_scalars            = this->weno_scalars           ;
    auto &weno_winds              = this->weno_winds             ;
    auto &c2g                     = this->coefs_to_gll           ;
    auto &s2g                     = this->sten_to_gll            ;
    auto &wenoRecon               = this->wenoRecon              ;
    auto &idl                     = this->idl                    ;
    auto &sigma                   = this->sigma                  ;
    auto &hyDensCells             = this->hyDensCells            ;
    auto &hyDensThetaCells        = this->hyDensThetaCells       ;
    auto &sim2d                   = this->sim2d                  ;
    auto &derivMatrix             = this->derivMatrix            ;
    auto &dx                      = this->dx                     ;
    auto &stateLimits             = this->stateLimits            ;
    auto &stateFlux               = this->stateFlux              ;
    auto &bc_x                    = this->bc_x                   ;
    auto &Rd                      = this->Rd                     ;
    auto &cp                      = this->cp                     ;
    auto &gamma                   = this->gamma                  ;
    auto &p0                      = this->p0                     ;
    auto &C0                      = this->C0                     ;

    // Populate the halos
    if        (bc_x == BC_PERIODIC) {
      parallel_for( SimpleBounds<3>(nz,ny,hs) , YAKL_LAMBDA(int k, int j, int ii) {
        for (int l=0; l < num_state; l++) {
          state  (l,hs+k,hs+j,      ii) = state  (l,hs+k,hs+j,nx+ii);
          state  (l,hs+k,hs+j,hs+nx+ii) = state  (l,hs+k,hs+j,hs+ii);
        }
      });
    } else if (bc_x == BC_WALL) {
      parallel_for( SimpleBounds<3>(nz,ny,hs) , YAKL_LAMBDA(int k, int j, int ii) {
        for (int l=0; l < num_state; l++) {
          if (l == idU) {
            state(l,hs+k,hs+j,      ii) = 0;
            state(l,hs+k,hs+j,hs+nx+ii) = 0;
          } else {
            state  (l,hs+k,hs+j,      ii) = state  (l,hs+k,hs+j,hs     );
            state  (l,hs+k,hs+j,hs+nx+ii) = state  (l,hs+k,hs+j,hs+nx-1);
          }
        }
      });
    }

    // Loop through all cells, reconstruct in x-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      // We need density and momentum to evolve the tracers with ADER
      SArray<real,2,nAder,ngll> ru_DTs;

      { // BEGIN: Reconstruct, time-average, and store the state and fluxes
        ////////////////////////////////////////////////////////////////
        // Reconstruct rho, u, v, w, theta
        ////////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> rv_DTs , rw_DTs , rt_DTs;
        { // BEGIN: Reconstruct the state
          SArray<real,1,ord> stencil;

          // u values and derivatives
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idU,hs+k,hs+j,i+ii); }
          reconstruct_gll_values( stencil , ru_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // v
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idV,hs+k,hs+j,i+ii); }
          reconstruct_gll_values( stencil , rv_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // w
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idW,hs+k,hs+j,i+ii); }
          reconstruct_gll_values( stencil , rw_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // theta
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idT,hs+k,hs+j,i+ii); }
          reconstruct_gll_values( stencil , rt_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
          for (int ii=0; ii < ngll; ii++) { rt_DTs(0,ii) += hyDensThetaCells(hs+k); } // Add hydrostasis back on
        } // END: Reconstruct the state

        ///////////////////////////////////////////////////////////////
        // Compute other values needed for centered tendencies and DTs
        ///////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> ruu_DTs , ruv_DTs , ruw_DTs , rut_DTs;
        for (int ii=0; ii < ngll; ii++) {
          real r = hyDensCells(hs+k);
          real u = ru_DTs(0,ii) / r;
          real v = rv_DTs(0,ii) / r;
          real w = rw_DTs(0,ii) / r;
          real t = rt_DTs(0,ii) / r;
          ruu_DTs     (0,ii) = r*u*u;
          ruv_DTs     (0,ii) = r*u*v;
          ruw_DTs     (0,ii) = r*u*w;
          rut_DTs     (0,ii) = r*u*t;
        }

        //////////////////////////////////////////
        // Compute time derivatives if necessary
        //////////////////////////////////////////
        if (nAder > 1) {
          diffTransformEulerConsX( hyDensCells(hs+k) , ru_DTs , rv_DTs , rw_DTs , rt_DTs , ruu_DTs , ruv_DTs , ruw_DTs ,
                                   rut_DTs , derivMatrix , dx );
        }

        //////////////////////////////////////////
        // Time average if necessary
        //////////////////////////////////////////
        // density and momentum can't be overwritten because they will be used for tracers
        SArray<real,1,ngll> ru_tavg;
        if (timeAvg) {
          compute_timeAvg( ru_DTs , ru_tavg , dt );
          compute_timeAvg( rv_DTs           , dt );
          compute_timeAvg( rw_DTs           , dt );
          compute_timeAvg( rt_DTs           , dt );
          compute_timeAvg( ruu_DTs          , dt );
          compute_timeAvg( ruv_DTs          , dt );
          compute_timeAvg( ruw_DTs          , dt );
          compute_timeAvg( rut_DTs          , dt );
        } else {
          for (int ii=0; ii < ngll; ii++) {
            ru_tavg(ii) = ru_DTs(0,ii);
          }
        }

        //////////////////////////////////////////
        // Store cell edge estimates of the state
        //////////////////////////////////////////
        // Left interface
        stateLimits(idU,1,k,j,i  ) = ru_tavg (0     );
        stateLimits(idV,1,k,j,i  ) = rv_DTs(0,0     );
        stateLimits(idW,1,k,j,i  ) = rw_DTs(0,0     );
        stateLimits(idT,1,k,j,i  ) = rt_DTs(0,0     );
        // Right interface
        stateLimits(idU,0,k,j,i+1) = ru_tavg (ngll-1);
        stateLimits(idV,0,k,j,i+1) = rv_DTs(0,ngll-1);
        stateLimits(idW,0,k,j,i+1) = rw_DTs(0,ngll-1);
        stateLimits(idT,0,k,j,i+1) = rt_DTs(0,ngll-1);
      } // END: Reconstruct, time-average, and store the state and fluxes
    });

    ////////////////////////////////////////////////
    // BCs for the state edge estimates
    ////////////////////////////////////////////////
    parallel_for( SimpleBounds<2>(nz,ny) , YAKL_LAMBDA (int k, int j) {
      for (int l=0; l < num_state; l++) {
        if        (bc_x == BC_PERIODIC) {
          stateLimits(l,0,k,j,0 ) = stateLimits(l,0,k,j,nx);
          stateLimits(l,1,k,j,nx) = stateLimits(l,1,k,j,0 );
        } else if (bc_x == BC_WALL    ) {
          stateLimits(l,0,k,j,0 ) = stateLimits(l,1,k,j,0 );
          stateLimits(l,1,k,j,nx) = stateLimits(l,0,k,j,nx);
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the upwind fluxes
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<3>(nz,ny,nx+1) , YAKL_LAMBDA (int k, int j, int i) {
      real uint = stateLimits(idU,0,k,j,i) + stateLimits(idU,1,k,j,i);
      if (uint > 0) {
        stateFlux(idU,k,j,i) = stateLimits(idU,0,k,j,i)*stateLimits(idU,0,k,j,i)/hyDensCells(hs+k);
        stateFlux(idV,k,j,i) = stateLimits(idU,0,k,j,i)*stateLimits(idV,0,k,j,i)/hyDensCells(hs+k);
        stateFlux(idW,k,j,i) = stateLimits(idU,0,k,j,i)*stateLimits(idW,0,k,j,i)/hyDensCells(hs+k);
        stateFlux(idT,k,j,i) = stateLimits(idU,0,k,j,i)*stateLimits(idT,0,k,j,i)/hyDensCells(hs+k);
      } else if (uint < 0) {
        stateFlux(idU,k,j,i) = stateLimits(idU,1,k,j,i)*stateLimits(idU,1,k,j,i)/hyDensCells(hs+k);
        stateFlux(idV,k,j,i) = stateLimits(idU,1,k,j,i)*stateLimits(idV,1,k,j,i)/hyDensCells(hs+k);
        stateFlux(idW,k,j,i) = stateLimits(idU,1,k,j,i)*stateLimits(idW,1,k,j,i)/hyDensCells(hs+k);
        stateFlux(idT,k,j,i) = stateLimits(idU,1,k,j,i)*stateLimits(idT,1,k,j,i)/hyDensCells(hs+k);
      } else {
        stateFlux(idU,k,j,i) = 0.5_fp * ( stateLimits(idU,1,k,j,i)*stateLimits(idU,1,k,j,i)/hyDensCells(hs+k) + stateLimits(idU,0,k,j,i)*stateLimits(idU,0,k,j,i)/hyDensCells(hs+k) );
        stateFlux(idV,k,j,i) = 0.5_fp * ( stateLimits(idU,1,k,j,i)*stateLimits(idV,1,k,j,i)/hyDensCells(hs+k) + stateLimits(idU,0,k,j,i)*stateLimits(idV,0,k,j,i)/hyDensCells(hs+k) );
        stateFlux(idW,k,j,i) = 0.5_fp * ( stateLimits(idU,1,k,j,i)*stateLimits(idW,1,k,j,i)/hyDensCells(hs+k) + stateLimits(idU,0,k,j,i)*stateLimits(idW,0,k,j,i)/hyDensCells(hs+k) );
        stateFlux(idT,k,j,i) = 0.5_fp * ( stateLimits(idU,1,k,j,i)*stateLimits(idT,1,k,j,i)/hyDensCells(hs+k) + stateLimits(idU,0,k,j,i)*stateLimits(idT,0,k,j,i)/hyDensCells(hs+k) );
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA(int k, int j, int i) {
      for (int l = 1; l < num_state; l++) {
        if (sim2d && l == idV) {
          stateTend(l,k,j,i) = 0;
        } else {
          stateTend(l,k,j,i) = - ( stateFlux(l,k,j,i+1) - stateFlux(l,k,j,i) ) / dx;
        }
      }
    });
  }



  template <class MICRO>
  void computeTendenciesY( real4d &state   , real4d &stateTend  ,
                           real4d &tracers , real4d &tracerTend ,
                           MICRO const &micro, real &dt ) {
    auto &ny                      = this->ny                     ;
    auto &weno_scalars            = this->weno_scalars           ;
    auto &weno_winds              = this->weno_winds             ;
    auto &c2g                     = this->coefs_to_gll           ;
    auto &s2g                     = this->sten_to_gll            ;
    auto &wenoRecon               = this->wenoRecon              ;
    auto &idl                     = this->idl                    ;
    auto &sigma                   = this->sigma                  ;
    auto &hyDensCells             = this->hyDensCells            ;
    auto &hyDensThetaCells        = this->hyDensThetaCells       ;
    auto &sim2d                   = this->sim2d                  ;
    auto &derivMatrix             = this->derivMatrix            ;
    auto &dy                      = this->dy                     ;
    auto &stateLimits             = this->stateLimits            ;
    auto &stateFlux               = this->stateFlux              ;
    auto &bc_y                    = this->bc_y                   ;
    auto &Rd                      = this->Rd                     ;
    auto &cp                      = this->cp                     ;
    auto &gamma                   = this->gamma                  ;
    auto &p0                      = this->p0                     ;
    auto &C0                      = this->C0                     ;

    // Populate the halos
    if        (bc_y == BC_PERIODIC) {
      parallel_for( SimpleBounds<3>(nz,nx,hs) , YAKL_LAMBDA(int k, int i, int jj) {
        for (int l=0; l < num_state; l++) {
          state(l,hs+k,      jj,hs+i) = state(l,hs+k,ny+jj,hs+i);
          state(l,hs+k,hs+ny+jj,hs+i) = state(l,hs+k,hs+jj,hs+i);
        }
      });
    } else if (bc_y == BC_WALL) {
      parallel_for( SimpleBounds<3>(nz,nx,hs) , YAKL_LAMBDA(int k, int i, int jj) {
        for (int l=0; l < num_state; l++) {
          if (l == idV) {
            state(l,hs+k,      jj,hs+i) = 0;
            state(l,hs+k,hs+ny+jj,hs+i) = 0;
          } else {
            state(l,hs+k,      jj,hs+i) = state(l,hs+k,hs     ,hs+i);
            state(l,hs+k,hs+ny+jj,hs+i) = state(l,hs+k,hs+ny-1,hs+i);
          }
        }
      });
    }

    // Loop through all cells, reconstruct in y-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      // These are needed by the tracers
      SArray<real,2,nAder,ngll> rv_DTs;

      { // BEGIN: Reconstruct, time-average, and store state and sate fluxes
        ////////////////////////////////////////////////////////////////
        // Reconstruct rho, u, v, w, theta
        ////////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> ru_DTs , rw_DTs , rt_DTs;
        {
          SArray<real,1,ord> stencil;

          // u values and derivatives
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idU,hs+k,j+jj,hs+i); }
          reconstruct_gll_values( stencil , ru_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // v
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idV,hs+k,j+jj,hs+i); }
          reconstruct_gll_values( stencil , rv_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // w
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idW,hs+k,j+jj,hs+i); }
          reconstruct_gll_values( stencil , rw_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // theta
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idT,hs+k,j+jj,hs+i); }
          reconstruct_gll_values( stencil , rt_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
          for (int jj=0; jj < ngll; jj++) { rt_DTs(0,jj) += hyDensThetaCells(hs+k); } // Add hydrostasis back on
        }

        ///////////////////////////////////////////////////////////////
        // Compute other values needed for centered tendencies and DTs
        ///////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> rvu_DTs , rvv_DTs , rvw_DTs , rvt_DTs;
        for (int jj=0; jj < ngll; jj++) {
          real r = hyDensCells(hs+k);
          real u = ru_DTs(0,jj) / r;
          real v = rv_DTs(0,jj) / r;
          real w = rw_DTs(0,jj) / r;
          real t = rt_DTs(0,jj) / r;
          rvu_DTs     (0,jj) = r*v*u;
          rvv_DTs     (0,jj) = r*v*v;
          rvw_DTs     (0,jj) = r*v*w;
          rvt_DTs     (0,jj) = r*v*t;
        }

        //////////////////////////////////////////
        // Compute time derivatives if necessary
        //////////////////////////////////////////
        if (nAder > 1) {
          diffTransformEulerConsY( hyDensCells(hs+k) , ru_DTs , rv_DTs , rw_DTs , rt_DTs , rvu_DTs , rvv_DTs , rvw_DTs ,
                                   rvt_DTs , derivMatrix , dy );
        }

        //////////////////////////////////////////
        // Time average if necessary
        //////////////////////////////////////////
        // Don't overwrite r and rv because we need them for tracers
        SArray<real,1,ngll> rv_tavg;
        if (timeAvg) {
          compute_timeAvg( ru_DTs           , dt );
          compute_timeAvg( rv_DTs , rv_tavg , dt );
          compute_timeAvg( rw_DTs           , dt );
          compute_timeAvg( rt_DTs           , dt );
          compute_timeAvg( rvu_DTs          , dt );
          compute_timeAvg( rvv_DTs          , dt );
          compute_timeAvg( rvw_DTs          , dt );
          compute_timeAvg( rvt_DTs          , dt );
        } else {
          for (int jj=0; jj < ngll; jj++) {
            rv_tavg(jj) = rv_DTs(0,jj);
          }
        }

        //////////////////////////////////////////
        // Store cell edge estimates of the state
        //////////////////////////////////////////
        // Left interface
        stateLimits(idU,1,k,j  ,i) = ru_DTs(0,0     );
        stateLimits(idV,1,k,j  ,i) = rv_tavg (0     );
        stateLimits(idW,1,k,j  ,i) = rw_DTs(0,0     );
        stateLimits(idT,1,k,j  ,i) = rt_DTs(0,0     );
        // Right interface       
        stateLimits(idU,0,k,j+1,i) = ru_DTs(0,ngll-1);
        stateLimits(idV,0,k,j+1,i) = rv_tavg (ngll-1);
        stateLimits(idW,0,k,j+1,i) = rw_DTs(0,ngll-1);
        stateLimits(idT,0,k,j+1,i) = rt_DTs(0,ngll-1);
      } // END: Reconstruct, time-average, and store state and sate fluxes
    });

    ////////////////////////////////////////////////
    // BCs for the state edge estimates
    ////////////////////////////////////////////////
    parallel_for( SimpleBounds<2>(nz,nx) , YAKL_LAMBDA (int k, int i) {
      for (int l=0; l < num_state; l++) {
        if        (bc_y == BC_PERIODIC) {
          stateLimits(l,0,k,0 ,i) = stateLimits(l,0,k,ny,i);
          stateLimits(l,1,k,ny,i) = stateLimits(l,1,k,0 ,i);
        } else if (bc_y == BC_WALL    ) {
          stateLimits(l,0,k,0 ,i) = stateLimits(l,1,k,0 ,i);
          stateLimits(l,1,k,ny,i) = stateLimits(l,0,k,ny,i);
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the upwind fluxes
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<3>(nz,ny+1,nx) , YAKL_LAMBDA (int k, int j, int i) {
      real vint = stateLimits(idV,0,k,j,i) + stateLimits(idV,1,k,j,i);
      if (vint > 0) {
        stateFlux(idU,k,j,i) = stateLimits(idV,0,k,j,i)*stateLimits(idU,0,k,j,i)/hyDensCells(hs+k);
        stateFlux(idV,k,j,i) = stateLimits(idV,0,k,j,i)*stateLimits(idV,0,k,j,i)/hyDensCells(hs+k);
        stateFlux(idW,k,j,i) = stateLimits(idV,0,k,j,i)*stateLimits(idW,0,k,j,i)/hyDensCells(hs+k);
        stateFlux(idT,k,j,i) = stateLimits(idV,0,k,j,i)*stateLimits(idT,0,k,j,i)/hyDensCells(hs+k);
      } else if (vint < 0) {
        stateFlux(idU,k,j,i) = stateLimits(idV,1,k,j,i)*stateLimits(idU,1,k,j,i)/hyDensCells(hs+k);
        stateFlux(idV,k,j,i) = stateLimits(idV,1,k,j,i)*stateLimits(idV,1,k,j,i)/hyDensCells(hs+k);
        stateFlux(idW,k,j,i) = stateLimits(idV,1,k,j,i)*stateLimits(idW,1,k,j,i)/hyDensCells(hs+k);
        stateFlux(idT,k,j,i) = stateLimits(idV,1,k,j,i)*stateLimits(idT,1,k,j,i)/hyDensCells(hs+k);
      } else {
        stateFlux(idU,k,j,i) = 0.5_fp * ( stateLimits(idV,1,k,j,i)*stateLimits(idU,1,k,j,i)/hyDensCells(hs+k) + stateLimits(idV,0,k,j,i)*stateLimits(idU,0,k,j,i)/hyDensCells(hs+k) );
        stateFlux(idV,k,j,i) = 0.5_fp * ( stateLimits(idV,1,k,j,i)*stateLimits(idV,1,k,j,i)/hyDensCells(hs+k) + stateLimits(idV,0,k,j,i)*stateLimits(idV,0,k,j,i)/hyDensCells(hs+k) );
        stateFlux(idW,k,j,i) = 0.5_fp * ( stateLimits(idV,1,k,j,i)*stateLimits(idW,1,k,j,i)/hyDensCells(hs+k) + stateLimits(idV,0,k,j,i)*stateLimits(idW,0,k,j,i)/hyDensCells(hs+k) );
        stateFlux(idT,k,j,i) = 0.5_fp * ( stateLimits(idV,1,k,j,i)*stateLimits(idT,1,k,j,i)/hyDensCells(hs+k) + stateLimits(idV,0,k,j,i)*stateLimits(idT,0,k,j,i)/hyDensCells(hs+k) );
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA(int k, int j, int i) {
      for (int l=1; l < num_state; l++) {
        stateTend(l,k,j,i) = - ( stateFlux(l,k,j+1,i) - stateFlux(l,k,j,i) ) / dy;
      }
    });
  }



  template <class MICRO>
  void computeTendenciesZ( real4d &state   , real4d &stateTend  ,
                           real4d &tracers , real4d &tracerTend ,
                           MICRO const &micro, real &dt ) {
    auto &nz                      = this->nz                     ;
    auto &weno_scalars            = this->weno_scalars           ;
    auto &weno_winds              = this->weno_winds             ;
    auto &c2g                     = this->coefs_to_gll           ;
    auto &s2g                     = this->sten_to_gll            ;
    auto &wenoRecon               = this->wenoRecon              ;
    auto &idl                     = this->idl                    ;
    auto &sigma                   = this->sigma                  ;
    auto &hyDensCells             = this->hyDensCells            ;
    auto &hyThetaCells            = this->hyThetaCells           ;
    auto &hyDensThetaCells        = this->hyDensThetaCells       ;
    auto &hyDensGLL               = this->hyDensGLL              ;
    auto &hyDensThetaGLL          = this->hyDensThetaGLL         ;
    auto &hyPressureGLL           = this->hyPressureGLL          ;
    auto &sim2d                   = this->sim2d                  ;
    auto &derivMatrix             = this->derivMatrix            ;
    auto &dz                      = this->dz                     ;
    auto &stateLimits             = this->stateLimits            ;
    auto &stateFlux               = this->stateFlux              ;
    auto &bc_z                    = this->bc_z                   ;
    auto &gllWts_ngll             = this->gllWts_ngll            ;
    auto &Rd                      = this->Rd                     ;
    auto &cp                      = this->cp                     ;
    auto &gamma                   = this->gamma                  ;
    auto &p0                      = this->p0                     ;
    auto &C0                      = this->C0                     ;

    // Populate the halos
    if        (bc_z == BC_PERIODIC) {
      parallel_for( SimpleBounds<3>(ny,nx,hs) , YAKL_LAMBDA(int j, int i, int kk) {
        for (int l=0; l < num_state; l++) {
          state(l,      kk,hs+j,hs+i) = state(l,nz+kk,hs+j,hs+i);
          state(l,hs+nz+kk,hs+j,hs+i) = state(l,hs+kk,hs+j,hs+i);
        }
      });
    } else if (bc_z == BC_WALL) {
      parallel_for( SimpleBounds<3>(ny,nx,hs) , YAKL_LAMBDA(int j, int i, int kk) {
        for (int l=0; l < num_state; l++) {
          if (l == idW) {
            state(l,      kk,hs+j,hs+i) = 0;
            state(l,hs+nz+kk,hs+j,hs+i) = 0;
          } else {
            state(l,      kk,hs+j,hs+i) = state(l,hs     ,hs+j,hs+i);
            state(l,hs+nz+kk,hs+j,hs+i) = state(l,hs+nz-1,hs+j,hs+i);
          }
        }
      });
    }

    // Loop through all cells, reconstruct in x-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      // We need these to persist to evolve tracers with ADER
      SArray<real,2,nAder,ngll> rw_DTs;

      { // BEGIN: reconstruct, time-avg, and store state & state fluxes
        ////////////////////////////////////////////////////////////////
        // Reconstruct rho, u, v, w, theta
        ////////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> ru_DTs , rv_DTs , rt_DTs;
        {
          SArray<real,1,ord> stencil;

          // u values and derivatives
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idU,k+kk,hs+j,hs+i); }
          reconstruct_gll_values( stencil , ru_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // v
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idV,k+kk,hs+j,hs+i); }
          reconstruct_gll_values( stencil , rv_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // w
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idW,k+kk,hs+j,hs+i); }
          reconstruct_gll_values( stencil , rw_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );
          if (bc_z == BC_WALL) {
            if (k == nz-1) rw_DTs(0,ngll-1) = 0;
            if (k == 0   ) rw_DTs(0,0     ) = 0;
          }

          // theta
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idT,k+kk,hs+j,hs+i); }
          reconstruct_gll_values( stencil , rt_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
          for (int kk=0; kk < ngll; kk++) { rt_DTs(0,kk) += hyDensThetaGLL(k,kk); } // Add hydrostasis back on
        }

        ///////////////////////////////////////////////////////////////
        // Compute other values needed for centered tendencies and DTs
        ///////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> rwu_DTs , rwv_DTs , rww_DTs , rwt_DTs;
        for (int kk=0; kk < ngll; kk++) {
          real r = hyDensGLL(k,kk);
          real u = ru_DTs(0,kk) / r;
          real v = rv_DTs(0,kk) / r;
          real w = rw_DTs(0,kk) / r;
          real t = rt_DTs(0,kk) / r;
          rwu_DTs    (0,kk) = r*w*u;
          rwv_DTs    (0,kk) = r*w*v;
          rww_DTs    (0,kk) = r*w*w;
          rwt_DTs    (0,kk) = r*w*t;
        }

        //////////////////////////////////////////
        // Compute time derivatives if necessary
        //////////////////////////////////////////
        if (nAder > 1) {
          diffTransformEulerConsZ( ru_DTs , rv_DTs , rw_DTs , rt_DTs , rwu_DTs , rwv_DTs , rww_DTs ,
                                   rwt_DTs , derivMatrix , hyDensGLL , k , dz , bc_z , nz );
        }

        //////////////////////////////////////////
        // Time average if necessary
        //////////////////////////////////////////
        // We can't alter density and momentum because they're needed for tracers later
        SArray<real,1,ngll> rw_tavg;
        if (timeAvg) {
          compute_timeAvg( ru_DTs           , dt );
          compute_timeAvg( rv_DTs           , dt );
          compute_timeAvg( rw_DTs , rw_tavg , dt );
          compute_timeAvg( rt_DTs           , dt );
          compute_timeAvg( rwu_DTs          , dt );
          compute_timeAvg( rwv_DTs          , dt );
          compute_timeAvg( rww_DTs          , dt );
          compute_timeAvg( rwt_DTs          , dt );
        } else {
          for (int ii=0; ii < ngll; ii++) {
            rw_tavg(ii) = rw_DTs(0,ii);
          }
        }

        //////////////////////////////////////////
        // Store cell edge estimates of the state
        //////////////////////////////////////////
        // Left interface
        stateLimits(idU,1,k  ,j,i) = ru_DTs(0,0     );
        stateLimits(idV,1,k  ,j,i) = rv_DTs(0,0     );
        stateLimits(idW,1,k  ,j,i) = rw_tavg (0     );
        stateLimits(idT,1,k  ,j,i) = rt_DTs(0,0     );
        // Right interface       
        stateLimits(idU,0,k+1,j,i) = ru_DTs(0,ngll-1);
        stateLimits(idV,0,k+1,j,i) = rv_DTs(0,ngll-1);
        stateLimits(idW,0,k+1,j,i) = rw_tavg (ngll-1);
        stateLimits(idT,0,k+1,j,i) = rt_DTs(0,ngll-1);
      } // END: reconstruct, time-avg, and store state & state fluxes
    });

    ////////////////////////////////////////////////
    // BCs for the state edge estimates
    ////////////////////////////////////////////////
    parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) {
      for (int l = 0; l < num_state; l++) {
        if        (bc_z == BC_PERIODIC) {
          stateLimits     (l,0,0 ,j,i) = stateLimits     (l,0,nz,j,i);
          stateLimits     (l,1,nz,j,i) = stateLimits     (l,1,0 ,j,i);
        } else if (bc_z == BC_WALL    ) {
          stateLimits     (l,0,0 ,j,i) = stateLimits     (l,1,0 ,j,i);
          stateLimits     (l,1,nz,j,i) = stateLimits     (l,0,nz,j,i);
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the upwind fluxes
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<3>(nz+1,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      real wint = stateLimits(idW,0,k,j,i) + stateLimits(idW,1,k,j,i);
      real rho_hy;
      if (k != nz) {
        rho_hy = hyDensGLL(k,0);
      } else {
        rho_hy = hyDensGLL(nz-1,ngll-1);
      }
      if (wint > 0 ) {
        stateFlux(idU,k,j,i) = stateLimits(idW,0,k,j,i)*stateLimits(idU,0,k,j,i)/rho_hy;
        stateFlux(idV,k,j,i) = stateLimits(idW,0,k,j,i)*stateLimits(idV,0,k,j,i)/rho_hy;
        stateFlux(idW,k,j,i) = stateLimits(idW,0,k,j,i)*stateLimits(idW,0,k,j,i)/rho_hy;
        stateFlux(idT,k,j,i) = stateLimits(idW,0,k,j,i)*stateLimits(idT,0,k,j,i)/rho_hy;
      } else if (wint < 0) {
        stateFlux(idU,k,j,i) = stateLimits(idW,1,k,j,i)*stateLimits(idU,1,k,j,i)/rho_hy;
        stateFlux(idV,k,j,i) = stateLimits(idW,1,k,j,i)*stateLimits(idV,1,k,j,i)/rho_hy;
        stateFlux(idW,k,j,i) = stateLimits(idW,1,k,j,i)*stateLimits(idW,1,k,j,i)/rho_hy;
        stateFlux(idT,k,j,i) = stateLimits(idW,1,k,j,i)*stateLimits(idT,1,k,j,i)/rho_hy;
      } else {
        stateFlux(idU,k,j,i) = 0.5_fp * ( stateLimits(idW,1,k,j,i)*stateLimits(idU,1,k,j,i)/rho_hy + stateLimits(idW,0,k,j,i)*stateLimits(idU,0,k,j,i)/rho_hy );
        stateFlux(idV,k,j,i) = 0.5_fp * ( stateLimits(idW,1,k,j,i)*stateLimits(idV,1,k,j,i)/rho_hy + stateLimits(idW,0,k,j,i)*stateLimits(idV,0,k,j,i)/rho_hy );
        stateFlux(idW,k,j,i) = 0.5_fp * ( stateLimits(idW,1,k,j,i)*stateLimits(idW,1,k,j,i)/rho_hy + stateLimits(idW,0,k,j,i)*stateLimits(idW,0,k,j,i)/rho_hy );
        stateFlux(idT,k,j,i) = 0.5_fp * ( stateLimits(idW,1,k,j,i)*stateLimits(idT,1,k,j,i)/rho_hy + stateLimits(idW,0,k,j,i)*stateLimits(idT,0,k,j,i)/rho_hy );
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA(int k, int j, int i) {
      for (int l=1; l < num_state; l++) {
        if (sim2d && l == idV) {
          stateTend(l,k,j,i) = 0;
        } else {
          stateTend(l,k,j,i) = - ( stateFlux(l,k+1,j,i) - stateFlux(l,k,j,i) ) / dz;
        }
      }

      real theta = (state(idT,hs+k,hs+j,hs+i) + hyDensThetaCells(hs+k)) / hyDensCells(hs+k);
      stateTend(idW,k,j,i) += hyDensCells(hs+k) * (theta - hyThetaCells(hs+k)) / hyThetaCells(hs+k) * GRAV;
    });
  }



  const char * getName() { return "1-D Uniform Transport with upwind FV on A-grid"; }



  template <class MICRO>
  void output(DataManager &dm, MICRO const &micro, real etime) const {
    auto &dx                    = this->dx                   ;
    auto &dy                    = this->dy                   ;
    auto &dz                    = this->dz                   ;
    auto &hyDensCells           = this->hyDensCells          ;
    auto &hyDensThetaCells      = this->hyDensThetaCells     ;
    auto &hyThetaCells          = this->hyThetaCells         ;
    auto &hyPressureCells       = this->hyPressureCells      ;
    auto &gamma                 = this->gamma                ;
    auto &C0                    = this->C0                   ;
    auto &pressure              = this->pressure             ;

    yakl::SimpleNetCDF nc;
    int ulIndex = 0; // Unlimited dimension index to place this data at
    // Create or open the file
    if (etime == 0.) {
      nc.create(out_file);
      // x-coordinate
      real1d xloc("xloc",nx);
      parallel_for( nx , YAKL_LAMBDA (int i) { xloc(i) = (i+0.5)*dx; });
      nc.write(xloc.createHostCopy(),"x",{"x"});
      // y-coordinate
      real1d yloc("yloc",ny);
      parallel_for( ny , YAKL_LAMBDA (int i) { yloc(i) = (i+0.5)*dy; });
      nc.write(yloc.createHostCopy(),"y",{"y"});
      // z-coordinate
      real1d zloc("zloc",nz);
      parallel_for( nz , YAKL_LAMBDA (int i) { zloc(i) = (i+0.5)*dz; });
      nc.write(zloc.createHostCopy(),"z",{"z"});
      // hydrostatic density, theta, and pressure
      nc.write(dm.get<real,1>("hydrostatic_density" ).createHostCopy(),"hyDens"    ,{"z"});
      nc.write(dm.get<real,1>("hydrostatic_pressure").createHostCopy(),"hyPressure",{"z"});
      nc.write(dm.get<real,1>("hydrostatic_theta"   ).createHostCopy(),"hyTheta"   ,{"z"});
      // Create time variable
      nc.write1(0._fp,"t",0,"t");
    } else {
      nc.open(out_file,yakl::NETCDF_MODE_WRITE);
      ulIndex = nc.getDimSize("t");
      // Write the elapsed time
      nc.write1(etime,"t",ulIndex,"t");
    }

    real4d state   = createStateArr ();
    real4d tracers = createTracerArr();
    read_state_and_tracers( dm , state , tracers );

    real3d data("data",nz,ny,nx);
    // rho'
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = state(idR,hs+k,hs+j,hs+i); });
    nc.write1(data.createHostCopy(),"dens_pert",{"z","y","x"},ulIndex,"t");
    // u
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = state(idU,hs+k,hs+j,hs+i); });
    nc.write1(data.createHostCopy(),"u",{"z","y","x"},ulIndex,"t");
    // v
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = state(idV,hs+k,hs+j,hs+i); });
    nc.write1(data.createHostCopy(),"v",{"z","y","x"},ulIndex,"t");
    // w
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) { data(k,j,i) = state(idW,hs+k,hs+j,hs+i); });
    nc.write1(data.createHostCopy(),"w",{"z","y","x"},ulIndex,"t");
    // theta'
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      real r = hyDensCells(hs+k);
      real t = ( state(idT,hs+k,hs+j,hs+i) + hyDensThetaCells(hs+k) ) / r;
      data(k,j,i) = t - hyThetaCells(hs+k);
    });
    nc.write1(data.createHostCopy(),"pot_temp_pert",{"z","y","x"},ulIndex,"t");
    // pressure'
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      data(k,j,i) = pressure(hs+k,hs+j,hs+i) * hyDensCells(hs+k);
    });
    nc.write1(data.createHostCopy(),"pressure_pert",{"z","y","x"},ulIndex,"t");

    // Close the file
    nc.close();
  }



  void finalize(real4d const &state , real4d const &tracers) {}



  // ord stencil values to ngll GLL values; store in DTs
  YAKL_INLINE void reconstruct_gll_values( SArray<real,1,ord> const stencil ,
                                           SArray<real,2,nAder,ngll> &DTs ,
                                           SArray<real,2,ord,ngll> const &coefs_to_gll ,
                                           SArray<real,2,ord,ngll> const &sten_to_gll  ,
                                           SArray<real,3,ord,ord,ord> const &wenoRecon ,
                                           SArray<real,1,hs+2> const &idl              ,
                                           real sigma, bool doweno ) {
    if (doweno) {

      // Reconstruct values
      SArray<real,1,ord> wenoCoefs;
      weno::compute_weno_coefs( wenoRecon , stencil , wenoCoefs , idl , sigma );
      // Transform ord weno coefficients into ngll GLL points
      for (int ii=0; ii<ngll; ii++) {
        real tmp = 0;
        for (int s=0; s < ord; s++) {
          tmp += coefs_to_gll(s,ii) * wenoCoefs(s);
        }
        DTs(0,ii) = tmp;
      }

    } else {

      // Transform ord stencil cell averages into ngll GLL points
      for (int ii=0; ii<ngll; ii++) {
        real tmp = 0;
        for (int s=0; s < ord; s++) {
          tmp += sten_to_gll(s,ii) * stencil(s);
        }
        DTs(0,ii) = tmp;
      }

    } // if doweno
  }



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
      weno::compute_weno_coefs( wenoRecon , stencil , wenoCoefs , idl , sigma );
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



  YAKL_INLINE void diffTransformEulerConsX( real hy_r  ,
                                            SArray<real,2,nAder,ngll> &ru ,
                                            SArray<real,2,nAder,ngll> &rv ,
                                            SArray<real,2,nAder,ngll> &rw ,
                                            SArray<real,2,nAder,ngll> &rt ,
                                            SArray<real,2,nAder,ngll> &ruu ,
                                            SArray<real,2,nAder,ngll> &ruv ,
                                            SArray<real,2,nAder,ngll> &ruw ,
                                            SArray<real,2,nAder,ngll> &rut ,
                                            SArray<real,2,ngll,ngll> const &deriv ,
                                            real dx ) {
    // Loop over the time derivatives
    for (int kt=0; kt<nAder-1; kt++) {
      // Compute the state at the next time level
      for (int ii=0; ii<ngll; ii++) {
        real df2_dx = 0;
        real df3_dx = 0;
        real df4_dx = 0;
        real df5_dx = 0;
        for (int s=0; s<ngll; s++) {
          df2_dx += deriv(s,ii) * ( ruu(kt,s) );
          df3_dx += deriv(s,ii) * ( ruv(kt,s) );
          df4_dx += deriv(s,ii) * ( ruw(kt,s) );
          df5_dx += deriv(s,ii) * ( rut(kt,s) );
        }
        ru(kt+1,ii) = -df2_dx/dx/(kt+1._fp);
        rv(kt+1,ii) = -df3_dx/dx/(kt+1._fp);
        rw(kt+1,ii) = -df4_dx/dx/(kt+1._fp);
        rt(kt+1,ii) = -df5_dx/dx/(kt+1._fp);
      }
      // Compute ru* at the next time level
      for (int ii=0; ii<ngll; ii++) {
        real tot_ruu = 0;
        real tot_ruv = 0;
        real tot_ruw = 0;
        real tot_rut = 0;
        for (int ir=0; ir<=kt+1; ir++) {
          tot_ruu += ru(ir,ii) * ru(kt+1-ir,ii);
          tot_ruv += ru(ir,ii) * rv(kt+1-ir,ii);
          tot_ruw += ru(ir,ii) * rw(kt+1-ir,ii);
          tot_rut += ru(ir,ii) * rt(kt+1-ir,ii);
        }
        ruu(kt+1,ii) = tot_ruu / hy_r;
        ruv(kt+1,ii) = tot_ruv / hy_r;
        ruw(kt+1,ii) = tot_ruw / hy_r;
        rut(kt+1,ii) = tot_rut / hy_r;
      }
    }
  }



  YAKL_INLINE void diffTransformEulerConsY( real hy_r ,
                                            SArray<real,2,nAder,ngll> &ru ,
                                            SArray<real,2,nAder,ngll> &rv ,
                                            SArray<real,2,nAder,ngll> &rw ,
                                            SArray<real,2,nAder,ngll> &rt ,
                                            SArray<real,2,nAder,ngll> &rvu ,
                                            SArray<real,2,nAder,ngll> &rvv ,
                                            SArray<real,2,nAder,ngll> &rvw ,
                                            SArray<real,2,nAder,ngll> &rvt ,
                                            SArray<real,2,ngll,ngll> const &deriv , 
                                            real dy ) {
    // Loop over the time derivatives
    for (int kt=0; kt<nAder-1; kt++) {
      // Compute the state at the next time level
      for (int ii=0; ii<ngll; ii++) {
        real drvu_dy = 0;
        real drvv_dy = 0;
        real drvw_dy = 0;
        real drvt_dy = 0;
        for (int s=0; s<ngll; s++) {
          drvu_dy += deriv(s,ii) * rvu(kt,s);
          drvv_dy += deriv(s,ii) * rvv(kt,s);
          drvw_dy += deriv(s,ii) * rvw(kt,s);
          drvt_dy += deriv(s,ii) * rvt(kt,s);
        }
        ru(kt+1,ii) = -drvu_dy/dy/(kt+1);
        rv(kt+1,ii) = -drvv_dy/dy/(kt+1);
        rw(kt+1,ii) = -drvw_dy/dy/(kt+1);
        rt(kt+1,ii) = -drvt_dy/dy/(kt+1);
      }

      // Compute ru* at the next time level
      for (int ii=0; ii<ngll; ii++) {
        // Compute the non-linear differential transforms
        real tot_rvu = 0;
        real tot_rvv = 0;
        real tot_rvw = 0;
        real tot_rvt = 0;
        for (int l=0; l<=kt+1; l++) {
          tot_rvu += rv(l,ii) * ru(kt+1-l,ii);
          tot_rvv += rv(l,ii) * rv(kt+1-l,ii);
          tot_rvw += rv(l,ii) * rw(kt+1-l,ii);
          tot_rvt += rv(l,ii) * rt(kt+1-l,ii);
        }
        rvu(kt+1,ii) = tot_rvu / hy_r;
        rvv(kt+1,ii) = tot_rvv / hy_r;
        rvw(kt+1,ii) = tot_rvw / hy_r;
        rvt(kt+1,ii) = tot_rvt / hy_r;
      }
    }
  }



  YAKL_INLINE void diffTransformEulerConsZ( SArray<real,2,nAder,ngll> &ru ,
                                            SArray<real,2,nAder,ngll> &rv ,
                                            SArray<real,2,nAder,ngll> &rw ,
                                            SArray<real,2,nAder,ngll> &rt ,
                                            SArray<real,2,nAder,ngll> &rwu ,
                                            SArray<real,2,nAder,ngll> &rwv ,
                                            SArray<real,2,nAder,ngll> &rww ,
                                            SArray<real,2,nAder,ngll> &rwt ,
                                            SArray<real,2,ngll,ngll> const &deriv , 
                                            real2d const &hyDensGLL , 
                                            int k , real dz , int bc_z , int nz ) {
    // Loop over the time derivatives
    for (int kt=0; kt<nAder-1; kt++) {
      // Compute the state at the next time level
      for (int ii=0; ii<ngll; ii++) {
        real drwu_dz = 0;
        real drwv_dz = 0;
        real drww_dz = 0;
        real drwt_dz = 0;
        for (int s=0; s<ngll; s++) {
          drwu_dz += deriv(s,ii) * rwu(kt,s);
          drwv_dz += deriv(s,ii) * rwv(kt,s);
          drww_dz += deriv(s,ii) * rww(kt,s);
          drwt_dz += deriv(s,ii) * rwt(kt,s);
        }
        ru(kt+1,ii) = -drwu_dz/dz/(kt+1);
        rv(kt+1,ii) = -drwv_dz/dz/(kt+1);
        rw(kt+1,ii) = -drww_dz/dz/(kt+1);
        rt(kt+1,ii) = -drwt_dz/dz/(kt+1);
        if (bc_z == BC_WALL) {
          if (k == nz-1) rw(kt+1,ngll-1) = 0;
          if (k == 0   ) rw(kt+1,0     ) = 0;
        }
      }
      // Compute ru* at the next time level
      for (int ii=0; ii<ngll; ii++) {
        // Compute the non-linear differential transforms
        real tot_rwu = 0;
        real tot_rwv = 0;
        real tot_rww = 0;
        real tot_rwt = 0;
        for (int l=0; l<=kt+1; l++) {
          tot_rwu += rw(l,ii) * ru(kt+1-l,ii);
          tot_rwv += rw(l,ii) * rv(kt+1-l,ii);
          tot_rww += rw(l,ii) * rw(kt+1-l,ii);
          tot_rwt += rw(l,ii) * rt(kt+1-l,ii);
        }
        rwu(kt+1,ii) = tot_rwu / hyDensGLL(k,ii);
        rwv(kt+1,ii) = tot_rwv / hyDensGLL(k,ii);
        rww(kt+1,ii) = tot_rww / hyDensGLL(k,ii);
        rwt(kt+1,ii) = tot_rwt / hyDensGLL(k,ii);
      }
    }
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


