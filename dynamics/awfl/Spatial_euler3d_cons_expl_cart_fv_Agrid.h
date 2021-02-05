
#pragma once

#include "const.h"
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

  int tracer_id_uniform;

  typedef real4d StateArr;  // Array of state variables (rho, rho*u, rho*v, rho*w, and rho*theta)
  typedef real4d TracerArr; // Array of tracers (total tracer mass)

  typedef real4d StateTendArr;   // State tendencies
  typedef real4d TracerTendArr;  // Tracer tendencies

  // Stores two estimates of the state, state flux, and tracer values at each cell interface
  real5d stateLimits;
  real5d tracerLimits;
  real4d stateFlux;

  // Stores single-valued flux of the tracer at each cell interface
  real4d tracerFlux;

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

  real1d vert_interface;
  real1d vert_interface_ghost;
  real2d vert_locs_normalized;
  real1d dz;
  real1d dz_ghost;
  real3d vert_sten_to_gll;
  real4d vert_weno_recon;

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
  // Whether to balance initial density to avoid acoustics at the start
  bool balance_initial_density;


  // When this class is created, initialize num_tracers to zero
  Spatial_operator() {
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
    auto &dx             = this->dx            ;
    auto &dy             = this->dy            ;
    auto &dz             = this->dz            ;
    auto &gllPts_ord     = this->gllPts_ord    ;
    auto &gllWts_ord     = this->gllWts_ord    ;
    auto &sim2d          = this->sim2d         ;
    auto &xlen           = this->xlen          ;
    auto &ylen           = this->ylen          ;
    auto &zlen           = this->zlen          ;
    auto &Rd             = this->Rd            ;
    auto &cp             = this->cp            ;
    auto &gamma          = this->gamma         ;
    auto &p0             = this->p0            ;
    auto &C0             = this->C0            ;
    auto &Rv             = this->Rv            ;
    auto &vert_interface = this->vert_interface;

    tracer_id_uniform = add_tracer(dm , "uniform"  , "uniform"  , false     , false);

    real3d dm_vapor   = dm.get<real,3>("water_vapor" );
    real3d dm_cloud   = dm.get<real,3>("cloud_liquid");
    real3d dm_precip  = dm.get<real,3>("precip_liquid");
    real3d dm_uniform = dm.get<real,3>("uniform");

    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      dm_vapor  (k,j,i) = 0;
      dm_cloud  (k,j,i) = 0;
      dm_precip (k,j,i) = 0;
      dm_uniform(k,j,i) = 0;
      // Loop over quadrature points
      for (int kk=0; kk<ord; kk++) {
        for (int jj=0; jj<ord; jj++) {
          for (int ii=0; ii<ord; ii++) {
            // Get the location
            real zloc = vert_interface(k) + 0.5_fp*dz(k) + gllPts_ord(kk)*dz(k);
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
            real svp   = micro.saturation_vapor_pressure(temp);     // Self-explanatory
            real p_v   = pert*svp;                                  // Multiply profile by saturation vapor pressure
            real r_v   = p_v / (Rv*temp);                           // Compute vapor density

            real wt = gllWts_ord(kk) * gllWts_ord(jj) * gllWts_ord(ii);
            dm_vapor  (k,j,i) += r_v / (rh+r_v) * rh * wt;
            dm_uniform(k,j,i) += rh * wt;
          }
        }
      }
    });
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
      rho      (k,j,i) = state(idR,hs+k,hs+j,hs+i) + hyDensCells(k)     ;
      rho_u    (k,j,i) = state(idU,hs+k,hs+j,hs+i)                         ;
      rho_v    (k,j,i) = state(idV,hs+k,hs+j,hs+i)                         ;
      rho_w    (k,j,i) = state(idW,hs+k,hs+j,hs+i)                         ;
      rho_theta(k,j,i) = state(idT,hs+k,hs+j,hs+i) + hyDensThetaCells(k);
      for (int tr=0; tr < num_tracers; tr++) {
        dm_tracers(tr,k,j,i) = tracers(tr,hs+k,hs+j,hs+i);
      }
    });
  }



  // Take an initially dry fluid state and adjust it to account for moist tracers
  template <class MICRO>
  void adjust_state_for_moisture(DataManager &dm , MICRO const &micro) const {
    auto &hyDensCells             = this->hyDensCells;
    auto &hyDensThetaCells        = this->hyDensThetaCells;
    auto &num_tracers             = this->num_tracers;
    auto &tracer_adds_mass        = this->tracer_adds_mass;
    auto &balance_initial_density = this->balance_initial_density;
    auto &tracer_id_uniform       = this->tracer_id_uniform;
    auto &Rd                      = this->Rd;
    auto &Rv                      = this->Rv;
    auto &C0                      = this->C0;
    auto &gamma                   = this->gamma;

    // Copy the DataManager data to state and tracer arrays for convenience
    real4d state   = createStateArr ();
    real4d tracers = createTracerArr();
    read_state_and_tracers( dm , state , tracers );

    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      // Add tracer density to dry density if it adds mass
      real rho_dry = state(idR,hs+k,hs+j,hs+i) + hyDensCells(k);
      for (int tr=0; tr < num_tracers; tr++) {
        if (tracer_adds_mass(tr)) {
          state(idR,hs+k,hs+j,hs+i) += tracers(tr,hs+k,hs+j,hs+i);
        }
      }
      real rho_moist = state(idR,hs+k,hs+j,hs+i) + hyDensCells(k);

      // Adjust momenta for moist density
      state(idU,hs+k,hs+j,hs+i) = state(idU,hs+k,hs+j,hs+i) / rho_dry * rho_moist;
      state(idV,hs+k,hs+j,hs+i) = state(idV,hs+k,hs+j,hs+i) / rho_dry * rho_moist;
      state(idW,hs+k,hs+j,hs+i) = state(idW,hs+k,hs+j,hs+i) / rho_dry * rho_moist;

      // Compute the dry temperature (same as the moist temperature)
      real rho_theta_dry = state(idT,hs+k,hs+j,hs+i) + hyDensThetaCells(k);
      real press = C0*pow(rho_theta_dry,gamma);  // Dry pressure
      real temp  = press / Rd / rho_dry;         // Temp (same dry or moist)

      // Compute moist theta
      real index_vapor = micro.tracer_index_vapor;
      real rho_v = tracers(index_vapor,hs+k,hs+j,hs+i);
      real R_moist = Rd * (rho_dry / rho_moist) + Rv * (rho_v / rho_moist);
      real press_moist = rho_moist * R_moist * temp;
      real rho_theta_moist = pow( press_moist / C0 , 1._fp/gamma );
      
      // Compute moist rho*theta
      state(idT,hs+k,hs+j,hs+i) = rho_theta_moist - hyDensThetaCells(k);

      for (int tr = 0 ; tr < num_tracers ; tr++) {
        tracers(tr,hs+k,hs+j,hs+i) = tracers(tr,hs+k,hs+j,hs+i) / rho_dry * rho_moist;
        if (tr == tracer_id_uniform) {
          tracers(tr,hs+k,hs+j,hs+i) = rho_moist;
        }
      }

      if (balance_initial_density) {
        real rh  = hyDensCells     (k);
        real rth = hyDensThetaCells(k);
        real rt = state(idT,hs+k,hs+j,hs+i) + rth;
        real r  = state(idR,hs+k,hs+j,hs+i) + rh;
        real t  = rt / r;
        r = rth/t;
        state(idR,hs+k,hs+j,hs+i) = r - rh;
        state(idU,hs+k,hs+j,hs+i) = state(idU,hs+k,hs+j,hs+i) / rho_moist * r;
        state(idV,hs+k,hs+j,hs+i) = state(idV,hs+k,hs+j,hs+i) / rho_moist * r;
        state(idW,hs+k,hs+j,hs+i) = state(idW,hs+k,hs+j,hs+i) / rho_moist * r;
        state(idT,hs+k,hs+j,hs+i) = r*t - rth;
        for (int tr = 0 ; tr < num_tracers ; tr++) {
          tracers(tr,hs+k,hs+j,hs+i) = tracers(tr,hs+k,hs+j,hs+i) / rho_moist * r;
          if (tr == tracer_id_uniform) {
            tracers(tr,hs+k,hs+j,hs+i) = rho_moist;
          }
        }
      }
    });

    // Copy the state and tracers arrays back to the DataManager
    write_state_and_tracers( dm , state , tracers );
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
    return 3;
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
        real r = state(idR,hs+k,hs+j,hs+i) + hyDensCells(k);
        real u = state(idU,hs+k,hs+j,hs+i) / r;
        real v = state(idV,hs+k,hs+j,hs+i) / r;
        real w = state(idW,hs+k,hs+j,hs+i) / r;
        real t = ( state(idT,hs+k,hs+j,hs+i) + hyDensThetaCells(k) ) / r;
        real p = C0*pow(r*t,gamma);

        // Compute the speed of sound (constant kappa assumption)
        real cs = sqrt(gamma*p/r);

        // Compute the maximum stable time step in each direction
        real udt = cfl * dx    / max( abs(u-cs) , abs(u+cs) );
        real vdt = cfl * dy    / max( abs(v-cs) , abs(v+cs) );
        real wdt = cfl * dz(k) / max( abs(w-cs) , abs(w+cs) );

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

    // Read in the vertical height cell interface locations
    std::string vcoords_file = config["vcoords"].as<std::string>();
    yakl::SimpleNetCDF nc;
    nc.open(vcoords_file);
    nz = nc.getDimSize("num_interfaces") - 1;
    real1d zint_in("zint_in",nz+1);
    nc.read(zint_in,"vertical_interfaces");
    nc.close();

    // Store vertical cell interface heights in the data manager
    dm.register_and_allocate<real>( "vertical_interface_height" , "vertical_interface_height" , {nz+1} , {"zp1"} );
    auto zint = dm.get<real,1>("vertical_interface_height");
    parallel_for( nz+1 , YAKL_LAMBDA (int k) { zint(k) = zint_in(k); });

    dm.register_and_allocate<real>( "vertical_midpoint_height" , "vertical_midpoint_heignt" , {nz} , {"z"} );
    auto zmid = dm.get<real,1>("vertical_midpoint_height");
    parallel_for( nz , YAKL_LAMBDA (int k) { zmid(k) = 0.5_fp*(zint_in(k) + zint_in(k+1)); });

    // Get the height of the z-dimension
    zlen = zint.createHostCopy()(nz);

    vert_interface       = real1d("vert_interface"      ,nz+1);
    vert_interface_ghost = real1d("vert_interface_ghost",nz+2*hs+1);
    vert_locs_normalized = real2d("vert_locs_normalized",nz,ord+1);
    dz                   = real1d("dz"                  ,nz);
    dz_ghost             = real1d("dz_ghost"            ,nz+2*hs);
    vert_sten_to_gll     = real3d("vert_sten_to_gll"    ,nz,ord,ngll);
    vert_weno_recon      = real4d("vert_weno_recon"     ,nz,ord,ord,ord);

    auto &vert_interface       = this->vert_interface      ;
    auto &vert_interface_ghost = this->vert_interface_ghost;
    auto &vert_locs_normalized = this->vert_locs_normalized;
    auto &dz                   = this->dz                  ;
    auto &dz_ghost             = this->dz_ghost            ;
    auto &vert_sten_to_gll     = this->vert_sten_to_gll    ;
    auto &vert_weno_recon      = this->vert_weno_recon     ;

    zint.deep_copy_to(vert_interface);
    
    parallel_for( nz , YAKL_LAMBDA (int k) {
      dz(k) = vert_interface(k+1) - vert_interface(k);
    });

    auto &nz = this->nz;
    parallel_for( nz+2*hs , YAKL_LAMBDA (int k) {
      if (k >= hs && k < hs+nz) {
        dz_ghost(k) = dz(k-hs);
      } else if (k < hs) {
        dz_ghost(k) = dz(0);
      } else if (k >= hs+nz) {
        dz_ghost(k) = dz(nz-1);
      }
    });

    parallel_for( 1 , YAKL_LAMBDA (int dummy) {
      vert_interface_ghost(0) = vert_interface(0) - hs*dz(0);
      for (int k=1; k < nz+2*hs+1; k++) {
        vert_interface_ghost(k) = vert_interface_ghost(k-1) + dz_ghost(k-1);
      }
    });

    auto vint_host      = vert_interface_ghost.createHostCopy();
    auto vert_s2g_host  = vert_sten_to_gll    .createHostCopy();
    auto vert_weno_host = vert_weno_recon     .createHostCopy();
    auto vert_locs_host = vert_locs_normalized.createHostCopy();

    SArray<real,2,ord,ngll> c2g;
    TransformMatrices::coefs_to_gll_lower(c2g);

    for (int k=0; k < nz; k++) {
      // Store stencil locations
      SArray<double,1,ord+1> locs;
      for (int kk=0; kk < ord+1; kk++) {
        locs(kk) = vint_host(k+kk);
      }

      // Normalize stencil locations
      double zmid = ( locs(hs+1) + locs(hs) ) / 2;
      double dzmid = locs(hs+1) - locs(hs);
      for (int kk=0; kk < ord+1; kk++) {
        locs(kk) = ( locs(kk) - zmid ) / dzmid;
        vert_locs_host(k,kk) = locs(kk);
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
          vert_s2g_host(k,jj,ii) = s2g_var(jj,ii);
        }
      }

      for (int kk=0; kk < ord; kk++) {
        for (int jj=0; jj < ord; jj++) {
          for (int ii=0; ii < ord; ii++) {
            vert_weno_host(k,kk,jj,ii) = weno_recon_var(kk,jj,ii);
          }
        }
      }

    }

    vert_s2g_host .deep_copy_to(vert_sten_to_gll    );
    vert_weno_host.deep_copy_to(vert_weno_recon     );
    vert_locs_host.deep_copy_to(vert_locs_normalized);

    // Read the # cells in each dimension
    nx = config["nx"].as<int>();
    ny = config["ny"].as<int>();

    // Determine whether this is a 2-D simulation
    sim2d = ny == 1;

    // Read the domain length in each dimension
    xlen = config["xlen"].as<real>();
    ylen = config["ylen"].as<real>();

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
    out_file = config["out_file"].as<std::string>();

    balance_initial_density = config["balance_initial_density"].as<bool>();

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
    stateLimits     = real5d("stateLimits"    ,num_state  ,2,nz+1,ny+1,nx+1);
    tracerLimits    = real5d("tracerLimits"   ,num_tracers,2,nz+1,ny+1,nx+1);
    stateFlux       = real4d("stateFlux"      ,num_state    ,nz+1,ny+1,nx+1);
    tracerFlux      = real4d("tracerFlux"     ,num_tracers  ,nz+1,ny+1,nx+1);
    hyDensCells          = real1d("hyDensCells       ",nz);
    hyPressureCells      = real1d("hyPressureCells   ",nz);
    hyThetaCells         = real1d("hyThetaCells      ",nz);
    hyDensThetaCells     = real1d("hyDensThetaCells  ",nz);
    hyDensGLL            = real2d("hyDensGLL         ",nz,ngll);
    hyPressureGLL        = real2d("hyPressureGLL     ",nz,ngll);
    hyThetaGLL           = real2d("hyThetaGLL        ",nz,ngll);
    hyDensThetaGLL       = real2d("hyDensThetaGLL    ",nz,ngll);

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

    std::cout << "nx: " << nx << "\n";
    std::cout << "ny: " << ny << "\n";
    std::cout << "nz: " << nz << "\n";
    std::cout << "xlen (m): " << xlen << "\n";
    std::cout << "ylen (m): " << ylen << "\n";
    std::cout << "zlen (m): " << zlen << "\n";
    std::cout << "Vertical coordinates file: " << vcoords_file << "\n";
    std::cout << "Simulation time (s): " << config["simTime"].as<real>() << "\n";
    std::cout << "Vertical interface heights: ";
    auto zint_host = zint.createHostCopy();
    for (int k=0; k < nz+1; k++) {
      std::cout << zint_host(k) << "  ";
    }
    std::cout << "\n\n";
  }



  // Initialize the state
  template <class MICRO>
  void init_state_and_tracers( DataManager &dm , MICRO const &micro ) {
    Rd    = micro.constants.R_d;
    cp    = micro.constants.cp_d;
    gamma = micro.constants.gamma_d;
    p0    = micro.constants.p0;
    C0    = micro.constants.C0_d;
    Rv    = micro.constants.R_v;

    auto nx                       = this->nx                     ;
    auto ny                       = this->ny                     ;
    auto nz                       = this->nz                     ;
    auto dx                       = this->dx                     ;
    auto dy                       = this->dy                     ;
    auto dz                       = this->dz                     ;
    auto dz_ghost                 = this->dz_ghost               ;
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
    auto &balance_initial_density = this->balance_initial_density;
    auto &vert_interface          = this->vert_interface         ;
    auto &vert_interface_ghost    = this->vert_interface_ghost   ;

    // Get Arrays for 1-D hydrostatic background profiles
    real1d dm_hyDens      = dm.get<real,1>( "hydrostatic_density"       );
    real1d dm_hyTheta     = dm.get<real,1>( "hydrostatic_theta"         );
    real1d dm_hyDensTheta = dm.get<real,1>( "hydrostatic_density_theta" );
    real1d dm_hyPressure  = dm.get<real,1>( "hydrostatic_pressure"      );
    
    real3d dm_rho       = dm.get<real,3>( "density"       );
    real3d dm_rho_u     = dm.get<real,3>( "density_u"     );
    real3d dm_rho_v     = dm.get<real,3>( "density_v"     );
    real3d dm_rho_w     = dm.get<real,3>( "density_w"     );
    real3d dm_rho_theta = dm.get<real,3>( "density_theta" );

    // If the data_spec is thermal or ..., then initialize the domain with Exner pressure-based hydrostasis
    // This is mostly to make plotting potential temperature perturbation easier for publications
    if (data_spec == DATA_SPEC_THERMAL) {

      // Setup hydrostatic background state
      parallel_for( SimpleBounds<1>(nz) , YAKL_LAMBDA (int k) {
        // Compute cell averages
        hyDensCells     (k) = 0;
        hyPressureCells (k) = 0;
        hyThetaCells    (k) = 0;
        hyDensThetaCells(k) = 0;
        for (int kk=0; kk<ord; kk++) {
          real zloc = vert_interface(k) + 0.5_fp*dz(k) + gllPts_ord(kk)*dz(k);
          if        (data_spec == DATA_SPEC_THERMAL) {
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
        dm_hyDens     (k) = hyDensCells     (k);
        dm_hyTheta    (k) = hyThetaCells    (k);
        dm_hyDensTheta(k) = hyDensThetaCells(k);
        dm_hyPressure (k) = hyPressureCells (k);
      });

      parallel_for( SimpleBounds<1>(nz) , YAKL_LAMBDA (int k) {
        // Compute ngll GLL points
        for (int kk=0; kk<ngll; kk++) {
          real zloc = vert_interface(k) + 0.5_fp*dz(k) + gllPts_ngll(kk)*dz(k);
          if        (data_spec == DATA_SPEC_THERMAL) {
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
              real zloc = vert_interface(k) + 0.5_fp*dz(k) + gllPts_ord(kk)*dz(k);
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

                dm_rho      (k,j,i) += (r - rh)*wt;
                dm_rho_theta(k,j,i) += (r*t - rh*th) * wt;
              }
            }
          }
        }
        dm_rho      (k,j,i) += dm_hyDens     (k);
        dm_rho_theta(k,j,i) += dm_hyDensTheta(k);
      });

      init_tracers( dm , micro );

      adjust_state_for_moisture( dm , micro );

    } // if (data_spec == DATA_SPEC_THERMAL)


    if (data_spec == DATA_SPEC_SUPERCELL) {

      tracer_id_uniform = add_tracer(dm , "uniform"  , "uniform"  , false     , false);

      // This uses a piecewise linear profile for Temperature
      real constexpr z_0    = 0;
      real constexpr z_trop = 12000;
      real           z_top  = zlen;
      real constexpr T_0    = 300;
      real constexpr T_trop = 213;
      real constexpr T_top  = 213;
      real constexpr p_0    = 100000;

      real3d quad_temp     ("quad_temp"     ,nz,ngll-1,ord);
      real2d press_dry_gll ("press_dry_gll" ,nz,ngll);
      real2d hyDensVapGLL  ("hyDensVapGLL"  ,nz,ngll);
      real1d hyDensVapCells("hyDensVapCells",nz);
      real1d z             ("z"             ,nz);
      real1d temp_hy       ("temp_hy"       ,nz);
      real1d tdew_hy       ("tdew_hy"       ,nz);

      // Compute full density at ord GLL points for the space between each cell
      for (int k=0; k < nz; k++) {              // k:   Loop over cells
        for (int kk=0; kk < ngll-1; kk++) {     // kk:  Loop over spaces between ngll GLL points within cells
          for (int kkk=0; kkk < ord; kkk++) {   // kkk: Loop over ord GLL points between ngll GLL points
            // Middle of this cell
            real cellmid   = vert_interface(k) + 0.5_fp*dz(k);
            // Bottom, top, and middle of the space between these two ngll GLL points
            real ngll_b    = cellmid + gllPts_ngll(kk  )*dz(k);
            real ngll_t    = cellmid + gllPts_ngll(kk+1)*dz(k);
            real ngll_m    = 0.5_fp * (ngll_b + ngll_t);
            // Compute grid spacing between these ngll GLL points
            real ngll_dz   = dz(k) * ( gllPts_ngll(kk+1) - gllPts_ngll(kk) );
            // Compute the locate of this ord GLL point within the ngll GLL points
            real zloc      = ngll_m + ngll_dz * gllPts_ord(kkk);
            // Compute full density at this location
            real temp      = profiles::init_supercell_temperature (zloc, z_0, z_trop, z_top, T_0, T_trop, T_top);
            real press_dry = profiles::init_supercell_pressure_dry(zloc, z_0, z_trop, z_top, T_0, T_trop, T_top, p_0, Rd);
            real dens_dry  = press_dry / (Rd*temp);
            real qvs       = profiles::init_supercell_sat_mix_dry(press_dry, temp);
            real relhum    = profiles::init_supercell_relhum(zloc, z_0, z_trop);
            if (relhum * qvs > 0.014_fp) relhum = 0.014_fp / qvs;
            real qv        = min( 0.014_fp , qvs*relhum );
            real d_qv_dz   = 0;
            if (qv > 0.014_fp) {
              real d_qvs_dT = profiles::init_supercell_sat_mix_dry_d_dT( press_dry , temp );
              real d_qvs_dp = profiles::init_supercell_sat_mix_dry_d_dp( press_dry , temp );
              real d_hum_dz = profiles::init_supercell_relhum_d_dz(zloc, z_0, z_trop);
              real dT_dz;
              if (zloc < z_trop) {
                dT_dz = (T_trop - T_0) / (z_trop - z_0);
              } else if (zloc == z_trop) {
                dT_dz  = (T_trop - T_0   ) / (z_trop - z_0   );
                dT_dz += (T_top  - T_trop) / (z_top  - z_trop);
                dT_dz *= 0.5_fp;
              } else {
                dT_dz = (T_top - T_trop) / (z_top - z_trop);
              }
              d_qv_dz = relhum * d_qvs_dT * dT_dz - relhum * GRAV * d_qvs_dp + qvs * d_hum_dz;
            }
            quad_temp(k,kk,kkk) = -( Rv/Rd*d_qv_dz + (1+qv)*GRAV/(Rd*temp) ) / (1 + Rv/Rd*qv);
          }
        }
      }

      // Now use quadrature to compute the pressure at each interface level using hydrostasis
      // This cannot be done in parallel due to loop-carried dependencies
      press_dry_gll(0,0) = p_0;
      for (int k=0; k < nz; k++) {
        for (int kk=0; kk < ngll-1; kk++) {
          real tot = 0;
          for (int kkk=0; kkk < ord; kkk++) {
            tot += quad_temp(k,kk,kkk) * gllWts_ord(kkk);
          }
          tot *= dz(k) * ( gllPts_ngll(kk+1) - gllPts_ngll(kk) );
          press_dry_gll(k,kk+1) = press_dry_gll(k,kk) * exp( tot );
          if (kk == ngll-2 && k < nz-1) {
            press_dry_gll(k+1,0) = press_dry_gll(k,ngll-1);
          }
        }
      }

      for (int k=0; k < nz; k++) {
        for (int kk=0; kk < ngll; kk++) {
          real zloc = vert_interface(k) + 0.5_fp*dz(k) + gllPts_ngll(kk)*dz(k);
          real temp       = profiles::init_supercell_temperature (zloc, z_0, z_trop, z_top, T_0, T_trop, T_top);
          real press_dry  = press_dry_gll(k,kk);
          real dens_dry   = press_dry / (Rd*temp);
          real press_tmp  = profiles::init_supercell_pressure_dry(zloc, z_0, z_trop, z_top, T_0, T_trop, T_top, p_0, Rd);
          real qvs        = profiles::init_supercell_sat_mix_dry(press_tmp, temp);
          real relhum     = profiles::init_supercell_relhum(zloc, z_0, z_trop);
          if (relhum * qvs > 0.014_fp) relhum = 0.014_fp / qvs;
          real qv         = min( 0.014_fp , qvs*relhum );
          real press_vap  = Rv / Rd * qv * press_dry;
          real dens_vap   = press_vap / (Rv*temp);
          real press      = press_dry + press_vap;
          real dens       = dens_dry + dens_vap;
          real dens_theta = pow( press / C0 , 1._fp / gamma );
          real theta      = dens_theta / dens;
          hyPressureGLL (k,kk) = press;
          hyDensGLL     (k,kk) = dens;
          hyDensThetaGLL(k,kk) = dens_theta;
          hyThetaGLL    (k,kk) = theta;
          hyDensVapGLL  (k,kk) = dens_vap;
        }
      }

      for (int k=0; k < nz; k++) {
        real press_tot      = 0;
        real dens_tot       = 0;
        real dens_vap_tot   = 0;
        real theta_tot      = 0;
        real dens_theta_tot = 0;
        for (int kk=0; kk < ngll; kk++) {
          press_tot      += hyPressureGLL (k,kk) * gllWts_ngll(kk);
          dens_tot       += hyDensGLL     (k,kk) * gllWts_ngll(kk);
          dens_vap_tot   += hyDensVapGLL  (k,kk) * gllWts_ngll(kk);
          dens_theta_tot += hyDensThetaGLL(k,kk) * gllWts_ngll(kk);
          theta_tot      += hyThetaGLL    (k,kk) * gllWts_ngll(kk);
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
        real zloc       = vert_interface(k) + 0.5_fp*dz(k);
        real press_tmp  = profiles::init_supercell_pressure_dry(zloc, z_0, z_trop, z_top, T_0, T_trop, T_top, p_0, Rd);
        real qvs        = profiles::init_supercell_sat_mix_dry(press_tmp, temp);
        real relhum     = qv / qvs;
        real T          = temp - 273;
        real a          = 17.27;
        real b          = 237.7;
        real tdew       = b * ( a*T / (b + T) + log(relhum) ) / ( a - ( a*T / (b+T) + log(relhum) ) );
        // The next three are just to confirm the skew-T diagram looks OK
        z               (k) = zloc;
        temp_hy         (k) = temp;
        tdew_hy         (k) = tdew;
        // These are used in the rest of the model
        hyPressureCells (k) = press;
        hyDensCells     (k) = dens;
        hyDensThetaCells(k) = dens_theta;
        hyThetaCells    (k) = theta;
        hyDensVapCells  (k) = dens_vap;

        dm_hyDens     (k) = hyDensCells     (k);
        dm_hyTheta    (k) = hyThetaCells    (k);
        dm_hyDensTheta(k) = hyDensThetaCells(k);
        dm_hyPressure (k) = hyPressureCells (k);
      }

      // Dump out data to plot a skew-T log-P diagram
      yakl::SimpleNetCDF nc;
      nc.create("skew.nc");
      nc.write(z              ,"z"          ,{"z"});
      nc.write(hyPressureCells,"pressure"   ,{"z"});
      nc.write(temp_hy        ,"temperature",{"z"});
      nc.write(tdew_hy        ,"dew_point"  ,{"z"});
      nc.close();

      real3d dm_vapor   = dm.get<real,3>("water_vapor" );
      real3d dm_cloud   = dm.get<real,3>("cloud_liquid");
      real3d dm_precip  = dm.get<real,3>("precip_liquid");

      for (int k=0; k < nz; k++) {
        for (int j=0; j < ny; j++) {
          for (int i=0; i < nx; i++) {
            dm_rho      (k,j,i) = 0;
            dm_rho_u    (k,j,i) = 0;
            dm_rho_v    (k,j,i) = 0;
            dm_rho_w    (k,j,i) = 0;
            dm_rho_theta(k,j,i) = 0;
            dm_vapor    (k,j,i) = 0;
            dm_cloud    (k,j,i) = 0;
            dm_precip   (k,j,i) = 0;
            for (int kk=0; kk < ngll; kk++) {
              for (int jj=0; jj < ngll; jj++) {
                for (int ii=0; ii < ngll; ii++) {
                  real xloc = (i+0.5_fp)*dx                    + gllPts_ngll(ii)*dx;
                  real yloc = (j+0.5_fp)*dy                    + gllPts_ngll(jj)*dy;
                  real zloc = vert_interface(k) + 0.5_fp*dz(k) + gllPts_ngll(kk)*dz(k);

                  real dens = hyDensGLL(k,kk);

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
                  real theta      = hyThetaGLL    (k,kk);
                  real dens_vap   = hyDensVapGLL  (k,kk);
                  real dens_theta = hyDensThetaGLL(k,kk);

                  real factor = gllWts_ngll(ii) * gllWts_ngll(jj) * gllWts_ngll(kk);
                  dm_rho      (k,j,i) += dens        * factor;
                  dm_rho_u    (k,j,i) += dens * uvel * factor;
                  dm_rho_v    (k,j,i) += dens * vvel * factor;
                  dm_rho_w    (k,j,i) += dens * wvel * factor;
                  dm_rho_theta(k,j,i) += dens_theta  * factor;
                  dm_vapor    (k,j,i) += dens_vap    * factor;
                }
              }
            }
          }
        }
      }

    } // if (data_spec == DATA_SPEC_SUPERCELL)
  }



  // Compute state and tendency time derivatives from the state
  template <class MICRO>
  void computeTendencies( real4d &state   , real4d &stateTend  ,
                          real4d &tracers , real4d &tracerTend ,
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
    if (splitIndex == numSplit()-1) dimSwitch = ! dimSwitch;
  } // computeTendencies



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
    auto &tracerLimits            = this->tracerLimits           ;
    auto &stateFlux               = this->stateFlux              ;
    auto &tracerFlux              = this->tracerFlux             ;
    auto &tracer_pos              = this->tracer_pos             ;
    auto &num_tracers             = this->num_tracers            ;
    auto &bc_x                    = this->bc_x                   ;
    auto &Rd                      = this->Rd                     ;
    auto &cp                      = this->cp                     ;
    auto &gamma                   = this->gamma                  ;
    auto &p0                      = this->p0                     ;
    auto &C0                      = this->C0                     ;

    // Pre-process the tracers by dividing by density inside the domain
    // After this, we can reconstruct tracers only (not rho * tracer)
    // Also, compute dry density
    parallel_for( SimpleBounds<4>(num_tracers,nz,ny,nx) , YAKL_LAMBDA (int tr, int k, int j, int i) {
      tracers(tr,hs+k,hs+j,hs+i) /= (state(idR,hs+k,hs+j,hs+i) + hyDensCells(k));
    });

    // Populate the halos
    if        (bc_x == BC_PERIODIC) {
      parallel_for( SimpleBounds<3>(nz,ny,hs) , YAKL_LAMBDA(int k, int j, int ii) {
        for (int l=0; l < num_state; l++) {
          state  (l,hs+k,hs+j,      ii) = state  (l,hs+k,hs+j,nx+ii);
          state  (l,hs+k,hs+j,hs+nx+ii) = state  (l,hs+k,hs+j,hs+ii);
        }
        for (int l=0; l < num_tracers; l++) {
          tracers(l,hs+k,hs+j,      ii) = tracers(l,hs+k,hs+j,nx+ii);
          tracers(l,hs+k,hs+j,hs+nx+ii) = tracers(l,hs+k,hs+j,hs+ii);
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
        for (int l=0; l < num_tracers; l++) {
          tracers(l,hs+k,hs+j,      ii) = tracers(l,hs+k,hs+j,hs     );
          tracers(l,hs+k,hs+j,hs+nx+ii) = tracers(l,hs+k,hs+j,hs+nx-1);
        }
      });
    }

    // Loop through all cells, reconstruct in x-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      // We need density and momentum to evolve the tracers with ADER
      SArray<real,2,nAder,ngll> r_DTs , ru_DTs;

      { // BEGIN: Reconstruct, time-average, and store the state and fluxes
        ////////////////////////////////////////////////////////////////
        // Reconstruct rho, u, v, w, theta
        ////////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> rv_DTs , rw_DTs , rt_DTs;
        { // BEGIN: Reconstruct the state
          SArray<real,1,ord> stencil;

          // Density
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idR,hs+k,hs+j,i+ii); }
          reconstruct_gll_values( stencil , r_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
          for (int ii=0; ii < ngll; ii++) { r_DTs(0,ii) += hyDensCells(k); } // Add hydrostasis back on

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
          for (int ii=0; ii < ngll; ii++) { rt_DTs(0,ii) += hyDensThetaCells(k); } // Add hydrostasis back on
        } // END: Reconstruct the state

        ///////////////////////////////////////////////////////////////
        // Compute other values needed for centered tendencies and DTs
        ///////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> ruu_DTs , ruv_DTs , ruw_DTs , rut_DTs , rt_gamma_DTs;
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

        //////////////////////////////////////////
        // Compute time derivatives if necessary
        //////////////////////////////////////////
        if (nAder > 1) {
          diffTransformEulerConsX( r_DTs , ru_DTs , rv_DTs , rw_DTs , rt_DTs , ruu_DTs , ruv_DTs , ruw_DTs ,
                                   rut_DTs , rt_gamma_DTs , derivMatrix , C0 , gamma , dx );
        }

        //////////////////////////////////////////
        // Time average if necessary
        //////////////////////////////////////////
        // density and momentum can't be overwritten because they will be used for tracers
        SArray<real,1,ngll> r_tavg, ru_tavg;
        if (timeAvg) {
          compute_timeAvg( r_DTs  , r_tavg  , dt );
          compute_timeAvg( ru_DTs , ru_tavg , dt );
          compute_timeAvg( rv_DTs           , dt );
          compute_timeAvg( rw_DTs           , dt );
          compute_timeAvg( rt_DTs           , dt );
          compute_timeAvg( ruu_DTs          , dt );
          compute_timeAvg( ruv_DTs          , dt );
          compute_timeAvg( ruw_DTs          , dt );
          compute_timeAvg( rut_DTs          , dt );
          compute_timeAvg( rt_gamma_DTs     , dt );
        } else {
          for (int ii=0; ii < ngll; ii++) {
            r_tavg (ii) = r_DTs (0,ii);
            ru_tavg(ii) = ru_DTs(0,ii);
          }
        }

        //////////////////////////////////////////
        // Store cell edge estimates of the state
        //////////////////////////////////////////
        // Left interface
        stateLimits(idR,1,k,j,i  ) = r_tavg  (0     );
        stateLimits(idU,1,k,j,i  ) = ru_tavg (0     );
        stateLimits(idV,1,k,j,i  ) = rv_DTs(0,0     );
        stateLimits(idW,1,k,j,i  ) = rw_DTs(0,0     );
        stateLimits(idT,1,k,j,i  ) = rt_DTs(0,0     );
        // Right interface
        stateLimits(idR,0,k,j,i+1) = r_tavg  (ngll-1);
        stateLimits(idU,0,k,j,i+1) = ru_tavg (ngll-1);
        stateLimits(idV,0,k,j,i+1) = rv_DTs(0,ngll-1);
        stateLimits(idW,0,k,j,i+1) = rw_DTs(0,ngll-1);
        stateLimits(idT,0,k,j,i+1) = rt_DTs(0,ngll-1);
      } // END: Reconstruct, time-average, and store the state and fluxes

      // r_DTs and ru_DTs still exist and are computed
      { // BEGIN: Reconstruct, time-average, and store tracer fluxes
        // Only process one tracer at a time to save on local memory / register requirements
        for (int tr=0; tr < num_tracers; tr++) {
          SArray<real,2,nAder,ngll> rt_DTs; // Density * tracer
          { // BEGIN: Reconstruct the tracer
            SArray<real,1,ord> stencil;
            for (int ii=0; ii < ord; ii++) { stencil(ii) = tracers(tr,hs+k,hs+j,i+ii); }
            reconstruct_gll_values( stencil , rt_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
            for (int ii=0; ii < ngll; ii++) { rt_DTs(0,ii) *= r_DTs(0,ii); }
            if (tracer_pos(tr)) {
              for (int ii=0; ii < ngll; ii++) { rt_DTs(0,ii) = max( 0._fp , rt_DTs(0,ii) ); }
            }
          } // END: Reconstruct the tracer

          // Compute the tracer flux
          SArray<real,2,nAder,ngll> rut_DTs; // Density * uwind * tracer
          for (int ii=0; ii < ngll; ii++) {
            rut_DTs(0,ii) = rt_DTs(0,ii) * ru_DTs(0,ii) / r_DTs(0,ii);
          }

          //////////////////////////////////////////
          // Compute time derivatives if necessary
          //////////////////////////////////////////
          if (nAder > 1) {
            diffTransformTracer( r_DTs , ru_DTs , rt_DTs , rut_DTs , derivMatrix , dx );
          }

          //////////////////////////////////////////
          // Time average if necessary
          //////////////////////////////////////////
          if (timeAvg) {
            compute_timeAvg( rt_DTs  , dt );
            compute_timeAvg( rut_DTs , dt );
          }
          if (tracer_pos(tr)) {
            for (int ii=0; ii < ngll; ii++) { rt_DTs(0,ii) = max( 0._fp , rt_DTs(0,ii) ); }
          }

          ////////////////////////////////////////////////////////////
          // Store cell edge estimates of the tracer
          ////////////////////////////////////////////////////////////
          tracerLimits(tr,1,k,j,i  ) = rt_DTs (0,0     ); // Left interface
          tracerLimits(tr,0,k,j,i+1) = rt_DTs (0,ngll-1); // Right interface
        }
      } // END: Reconstruct, time-average, and store tracer fluxes

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
      for (int l=0; l < num_tracers; l++) {
        if        (bc_x == BC_PERIODIC) {
          tracerLimits(l,0,k,j,0 ) = tracerLimits(l,0,k,j,nx);
          tracerLimits(l,1,k,j,nx) = tracerLimits(l,1,k,j,0 );
        } else if (bc_x == BC_WALL    ) {
          tracerLimits(l,0,k,j,0 ) = tracerLimits(l,1,k,j,0 );
          tracerLimits(l,1,k,j,nx) = tracerLimits(l,0,k,j,nx);
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the upwind fluxes
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<3>(nz,ny,nx+1) , YAKL_LAMBDA (int k, int j, int i) {
      // Get left and right state
      real r_L = stateLimits(idR,0,k,j,i)    ;   real r_R = stateLimits(idR,1,k,j,i)    ;
      real u_L = stateLimits(idU,0,k,j,i)/r_L;   real u_R = stateLimits(idU,1,k,j,i)/r_R;
      real v_L = stateLimits(idV,0,k,j,i)/r_L;   real v_R = stateLimits(idV,1,k,j,i)/r_R;
      real w_L = stateLimits(idW,0,k,j,i)/r_L;   real w_R = stateLimits(idW,1,k,j,i)/r_R;
      real t_L = stateLimits(idT,0,k,j,i)/r_L;   real t_R = stateLimits(idT,1,k,j,i)/r_R;
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
      real q1_L = stateLimits(idR,0,k,j,i);   real q1_R = stateLimits(idR,1,k,j,i);
      real q2_L = stateLimits(idU,0,k,j,i);   real q2_R = stateLimits(idU,1,k,j,i);
      real q3_L = stateLimits(idV,0,k,j,i);   real q3_R = stateLimits(idV,1,k,j,i);
      real q4_L = stateLimits(idW,0,k,j,i);   real q4_R = stateLimits(idW,1,k,j,i);
      real q5_L = stateLimits(idT,0,k,j,i);   real q5_R = stateLimits(idT,1,k,j,i);
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

      stateFlux(idR,k,j,i) = q2;
      stateFlux(idU,k,j,i) = q2*q2/q1 + C0*pow(q5,gamma);
      stateFlux(idV,k,j,i) = q2*q3/q1;
      stateFlux(idW,k,j,i) = q2*q4/q1;
      stateFlux(idT,k,j,i) = q2*q5/q1;

      // COMPUTE UPWIND TRACER FLUXES
      // Handle it one tracer at a time
      for (int tr=0; tr < num_tracers; tr++) {
        if (u > 0) {
          tracerFlux(tr,k,j,i) = q2 * tracerLimits(tr,0,k,j,i) / r_L;
        } else {
          tracerFlux(tr,k,j,i) = q2 * tracerLimits(tr,1,k,j,i) / r_R;
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Limit the tracer fluxes for positivity
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<4>(num_tracers,nz,ny,nx+1) , YAKL_LAMBDA (int tr, int k, int j, int i) {
      real constexpr eps = 1.e-10;
      real u = 0.5_fp * ( stateLimits(idU,0,k,j,i) + stateLimits(idU,1,k,j,i) );
      // Solid wall BCs mean u == 0 at boundaries, so we assume periodic if u != 0
      if (tracer_pos(tr)) {
        // Compute and apply the flux reduction factor of the upwind cell
        if      (u > 0) {
          // upwind is to the left of this interface
          int ind_i = i-1;
          if (ind_i == -1) ind_i = nx-1;
          real f1 = min( tracerFlux(tr,k,j,ind_i  ) , 0._fp );
          real f2 = max( tracerFlux(tr,k,j,ind_i+1) , 0._fp );
          real fluxOut = dt*(f2-f1)/dx;
          real dens = state(idR,hs+k,hs+j,hs+ind_i) + hyDensCells(k);
          tracerFlux(tr,k,j,i) *= min( 1._fp , tracers(tr,hs+k,hs+j,hs+ind_i) * dens / (fluxOut + eps) );
        } else if (u < 0) {
          // upwind is to the right of this interface
          int ind_i = i;
          if (ind_i == nx) ind_i = 0;
          real f1 = min( tracerFlux(tr,k,j,ind_i  ) , 0._fp );
          real f2 = max( tracerFlux(tr,k,j,ind_i+1) , 0._fp );
          real fluxOut = dt*(f2-f1)/dx;
          real dens = state(idR,hs+k,hs+j,hs+ind_i) + hyDensCells(k);
          tracerFlux(tr,k,j,i) *= min( 1._fp , tracers(tr,hs+k,hs+j,hs+ind_i) * dens / (fluxOut + eps) );
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA(int k, int j, int i) {
      for (int l = 0; l < num_state; l++) {
        if (sim2d && l == idV) {
          stateTend(l,k,j,i) = 0;
        } else {
          stateTend(l,k,j,i) = - ( stateFlux(l,k,j,i+1) - stateFlux(l,k,j,i) ) / dx;
        }
      }
      for (int l = 0; l < num_tracers; l++) {
        // Compute tracer tendency
        tracerTend(l,k,j,i) = - ( tracerFlux(l,k,j,i+1) - tracerFlux(l,k,j,i  ) ) / dx;
        // Multiply density back onto tracers
        tracers(l,hs+k,hs+j,hs+i) *= (state(idR,hs+k,hs+j,hs+i) + hyDensCells(k));
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
    auto &tracerLimits            = this->tracerLimits           ;
    auto &tracerFlux              = this->tracerFlux             ;
    auto &tracer_pos              = this->tracer_pos             ;
    auto &num_tracers             = this->num_tracers            ;
    auto &bc_y                    = this->bc_y                   ;
    auto &Rd                      = this->Rd                     ;
    auto &cp                      = this->cp                     ;
    auto &gamma                   = this->gamma                  ;
    auto &p0                      = this->p0                     ;
    auto &C0                      = this->C0                     ;

    // Pre-process the tracers by dividing by density inside the domain
    // After this, we can reconstruct tracers only (not rho * tracer)
    parallel_for( SimpleBounds<4>(num_tracers,nz,ny,nx) , YAKL_LAMBDA (int tr, int k, int j, int i) {
      tracers(tr,hs+k,hs+j,hs+i) /= (state(idR,hs+k,hs+j,hs+i) + hyDensCells(k));
    });

    // Populate the halos
    if        (bc_y == BC_PERIODIC) {
      parallel_for( SimpleBounds<3>(nz,nx,hs) , YAKL_LAMBDA(int k, int i, int jj) {
        for (int l=0; l < num_state; l++) {
          state(l,hs+k,      jj,hs+i) = state(l,hs+k,ny+jj,hs+i);
          state(l,hs+k,hs+ny+jj,hs+i) = state(l,hs+k,hs+jj,hs+i);
        }
        for (int l=0; l < num_tracers; l++) {
          tracers(l,hs+k,      jj,hs+i) = tracers(l,hs+k,ny+jj,hs+i);
          tracers(l,hs+k,hs+ny+jj,hs+i) = tracers(l,hs+k,hs+jj,hs+i);
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
        for (int l=0; l < num_tracers; l++) {
          tracers(l,hs+k,      jj,hs+i) = tracers(l,hs+k,hs     ,hs+i);
          tracers(l,hs+k,hs+ny+jj,hs+i) = tracers(l,hs+k,hs+ny-1,hs+i);
        }
      });
    }

    // Loop through all cells, reconstruct in y-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      // These are needed by the tracers
      SArray<real,2,nAder,ngll> r_DTs , rv_DTs;

      { // BEGIN: Reconstruct, time-average, and store state and sate fluxes
        ////////////////////////////////////////////////////////////////
        // Reconstruct rho, u, v, w, theta
        ////////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> ru_DTs , rw_DTs , rt_DTs;
        {
          SArray<real,1,ord> stencil;

          // Density
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idR,hs+k,j+jj,hs+i); }
          reconstruct_gll_values( stencil , r_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
          for (int jj=0; jj < ngll; jj++) { r_DTs(0,jj) += hyDensCells(k); } // Add hydrostasis back on

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
          for (int jj=0; jj < ngll; jj++) { rt_DTs(0,jj) += hyDensThetaCells(k); } // Add hydrostasis back on
        }

        ///////////////////////////////////////////////////////////////
        // Compute other values needed for centered tendencies and DTs
        ///////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> rvu_DTs , rvv_DTs , rvw_DTs , rvt_DTs , rt_gamma_DTs;
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

        //////////////////////////////////////////
        // Compute time derivatives if necessary
        //////////////////////////////////////////
        if (nAder > 1) {
          diffTransformEulerConsY( r_DTs , ru_DTs , rv_DTs , rw_DTs , rt_DTs , rvu_DTs , rvv_DTs , rvw_DTs ,
                                   rvt_DTs , rt_gamma_DTs , derivMatrix , C0 , gamma , dy );
        }

        //////////////////////////////////////////
        // Time average if necessary
        //////////////////////////////////////////
        // Don't overwrite r and rv because we need them for tracers
        SArray<real,1,ngll> r_tavg, rv_tavg;
        if (timeAvg) {
          compute_timeAvg( r_DTs  , r_tavg  , dt );
          compute_timeAvg( ru_DTs           , dt );
          compute_timeAvg( rv_DTs , rv_tavg , dt );
          compute_timeAvg( rw_DTs           , dt );
          compute_timeAvg( rt_DTs           , dt );
          compute_timeAvg( rvu_DTs          , dt );
          compute_timeAvg( rvv_DTs          , dt );
          compute_timeAvg( rvw_DTs          , dt );
          compute_timeAvg( rvt_DTs          , dt );
          compute_timeAvg( rt_gamma_DTs     , dt );
        } else {
          for (int jj=0; jj < ngll; jj++) {
            r_tavg (jj) = r_DTs (0,jj);
            rv_tavg(jj) = rv_DTs(0,jj);
          }
        }

        //////////////////////////////////////////
        // Store cell edge estimates of the state
        //////////////////////////////////////////
        // Left interface
        stateLimits(idR,1,k,j  ,i) = r_tavg  (0     );
        stateLimits(idU,1,k,j  ,i) = ru_DTs(0,0     );
        stateLimits(idV,1,k,j  ,i) = rv_tavg (0     );
        stateLimits(idW,1,k,j  ,i) = rw_DTs(0,0     );
        stateLimits(idT,1,k,j  ,i) = rt_DTs(0,0     );
        // Right interface       
        stateLimits(idR,0,k,j+1,i) = r_tavg  (ngll-1);
        stateLimits(idU,0,k,j+1,i) = ru_DTs(0,ngll-1);
        stateLimits(idV,0,k,j+1,i) = rv_tavg (ngll-1);
        stateLimits(idW,0,k,j+1,i) = rw_DTs(0,ngll-1);
        stateLimits(idT,0,k,j+1,i) = rt_DTs(0,ngll-1);
      } // END: Reconstruct, time-average, and store state and sate fluxes

      // r_DTs and rv_DTs still exist and are computed
      { // BEGIN: Reconstruct, time-average, and store tracer fluxes
        // Only process one tracer at a time to save on local memory / register requirements
        for (int tr=0; tr < num_tracers; tr++) {
          SArray<real,2,nAder,ngll> rt_DTs; // Density * tracer
          { // BEGIN: Reconstruct the tracer
            SArray<real,1,ord> stencil;
            for (int jj=0; jj < ord; jj++) { stencil(jj) = tracers(tr,hs+k,j+jj,hs+i); }
            reconstruct_gll_values( stencil , rt_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
            for (int jj=0; jj < ngll; jj++) { rt_DTs(0,jj) *= r_DTs(0,jj); }
            if (tracer_pos(tr)) {
              for (int jj=0; jj < ngll; jj++) { rt_DTs(0,jj) = max( 0._fp , rt_DTs(0,jj) ); }
            }
          } // END: Reconstruct the tracer

          // Compute the tracer flux
          SArray<real,2,nAder,ngll> rvt_DTs; // Density * vwind * tracer
          for (int jj=0; jj < ngll; jj++) {
            rvt_DTs(0,jj) = rt_DTs(0,jj) * rv_DTs(0,jj) / r_DTs(0,jj);
          }

          //////////////////////////////////////////
          // Compute time derivatives if necessary
          //////////////////////////////////////////
          if (nAder > 1) {
            diffTransformTracer( r_DTs , rv_DTs , rt_DTs , rvt_DTs , derivMatrix , dy );
          }

          //////////////////////////////////////////
          // Time average if necessary
          //////////////////////////////////////////
          if (timeAvg) {
            compute_timeAvg( rt_DTs  , dt );
            compute_timeAvg( rvt_DTs , dt );
          }
          if (tracer_pos(tr)) {
            for (int jj=0; jj < ngll; jj++) { rt_DTs(0,jj) = max( 0._fp , rt_DTs(0,jj) ); }
          }

          ////////////////////////////////////////////////////////////
          // Store cell edge estimates of the tracer
          ////////////////////////////////////////////////////////////
          tracerLimits(tr,1,k,j  ,i) = rt_DTs (0,0     ); // Left interface
          tracerLimits(tr,0,k,j+1,i) = rt_DTs (0,ngll-1); // Right interface
        }
      } // END: Reconstruct, time-average, and store tracer fluxes

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
      for (int l=0; l < num_tracers; l++) {
        if        (bc_y == BC_PERIODIC) {
          tracerLimits(l,0,k,0 ,i) = tracerLimits(l,0,k,ny,i);
          tracerLimits(l,1,k,ny,i) = tracerLimits(l,1,k,0 ,i);
        } else if (bc_y == BC_WALL    ) {
          tracerLimits(l,0,k,0 ,i) = tracerLimits(l,1,k,0 ,i);
          tracerLimits(l,1,k,ny,i) = tracerLimits(l,0,k,ny,i);
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the upwind fluxes
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<3>(nz,ny+1,nx) , YAKL_LAMBDA (int k, int j, int i) {
      // Get left and right state
      real r_L = stateLimits(idR,0,k,j,i)    ;   real r_R = stateLimits(idR,1,k,j,i)    ;
      real u_L = stateLimits(idU,0,k,j,i)/r_L;   real u_R = stateLimits(idU,1,k,j,i)/r_R;
      real v_L = stateLimits(idV,0,k,j,i)/r_L;   real v_R = stateLimits(idV,1,k,j,i)/r_R;
      real w_L = stateLimits(idW,0,k,j,i)/r_L;   real w_R = stateLimits(idW,1,k,j,i)/r_R;
      real t_L = stateLimits(idT,0,k,j,i)/r_L;   real t_R = stateLimits(idT,1,k,j,i)/r_R;
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
      real q1_L = stateLimits(idR,0,k,j,i);   real q1_R = stateLimits(idR,1,k,j,i);
      real q2_L = stateLimits(idU,0,k,j,i);   real q2_R = stateLimits(idU,1,k,j,i);
      real q3_L = stateLimits(idV,0,k,j,i);   real q3_R = stateLimits(idV,1,k,j,i);
      real q4_L = stateLimits(idW,0,k,j,i);   real q4_R = stateLimits(idW,1,k,j,i);
      real q5_L = stateLimits(idT,0,k,j,i);   real q5_R = stateLimits(idT,1,k,j,i);
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

      stateFlux(idR,k,j,i) = q3;
      stateFlux(idU,k,j,i) = q3*q2/q1;
      stateFlux(idV,k,j,i) = q3*q3/q1 + C0*pow(q5,gamma);
      stateFlux(idW,k,j,i) = q3*q4/q1;
      stateFlux(idT,k,j,i) = q3*q5/q1;

      // COMPUTE UPWIND TRACER FLUXES
      // Handle it one tracer at a time
      for (int tr=0; tr < num_tracers; tr++) {
        if (v > 0) {
          tracerFlux(tr,k,j,i) = q3 * tracerLimits(tr,0,k,j,i) / r_L;
        } else {
          tracerFlux(tr,k,j,i) = q3 * tracerLimits(tr,1,k,j,i) / r_R;
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Limit the tracer fluxes for positivity
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<4>(num_tracers,nz,ny+1,nx) , YAKL_LAMBDA (int tr, int k, int j, int i) {
      real constexpr eps = 1.e-10;
      real v = 0.5_fp * ( stateLimits(idV,0,k,j,i) + stateLimits(idV,1,k,j,i) );
      // Solid wall BCs mean u == 0 at boundaries, so we assume periodic if u != 0
      if (tracer_pos(tr)) {
        // Compute and apply the flux reduction factor of the upwind cell
        if      (v > 0) {
          // upwind is to the left of this interface
          int ind_j = j-1;
          if (ind_j == -1) ind_j = ny-1;
          real f1 = min( tracerFlux(tr,k,ind_j  ,i) , 0._fp );
          real f2 = max( tracerFlux(tr,k,ind_j+1,i) , 0._fp );
          real fluxOut = dt*(f2-f1)/dy;
          real dens = state(idR,hs+k,hs+ind_j,hs+i) + hyDensCells(k);
          tracerFlux(tr,k,j,i) *= min( 1._fp , tracers(tr,hs+k,hs+ind_j,hs+i) * dens / (fluxOut + eps) );
        } else if (v < 0) {
          // upwind is to the right of this interface
          int ind_j = j;
          if (ind_j == ny) ind_j = 0;
          real f1 = min( tracerFlux(tr,k,ind_j  ,i) , 0._fp );
          real f2 = max( tracerFlux(tr,k,ind_j+1,i) , 0._fp );
          real fluxOut = dt*(f2-f1)/dy;
          real dens = state(idR,hs+k,hs+ind_j,hs+i) + hyDensCells(k);
          tracerFlux(tr,k,j,i) *= min( 1._fp , tracers(tr,hs+k,hs+ind_j,hs+i) * dens / (fluxOut + eps) );
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA(int k, int j, int i) {
      for (int l=0; l < num_state; l++) {
        stateTend(l,k,j,i) = - ( stateFlux(l,k,j+1,i) - stateFlux(l,k,j,i) ) / dy;
      }
      for (int l=0; l < num_tracers; l++) {
        // Compute the tracer tendency
        tracerTend(l,k,j,i) = - ( tracerFlux(l,k,j+1,i) - tracerFlux(l,k,j,i) ) / dy;
        // Multiply density back onto the tracers
        tracers(l,hs+k,hs+j,hs+i) *= (state(idR,hs+k,hs+j,hs+i) + hyDensCells(k));
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
    auto &hyDensGLL               = this->hyDensGLL              ;
    auto &hyDensThetaGLL          = this->hyDensThetaGLL         ;
    auto &hyPressureGLL           = this->hyPressureGLL          ;
    auto &sim2d                   = this->sim2d                  ;
    auto &derivMatrix             = this->derivMatrix            ;
    auto &dz                      = this->dz                     ;
    auto &stateLimits             = this->stateLimits            ;
    auto &stateFlux               = this->stateFlux              ;
    auto &tracerLimits            = this->tracerLimits           ;
    auto &tracerFlux              = this->tracerFlux             ;
    auto &tracer_pos              = this->tracer_pos             ;
    auto &num_tracers             = this->num_tracers            ;
    auto &bc_z                    = this->bc_z                   ;
    auto &gllWts_ngll             = this->gllWts_ngll            ;
    auto &Rd                      = this->Rd                     ;
    auto &cp                      = this->cp                     ;
    auto &gamma                   = this->gamma                  ;
    auto &p0                      = this->p0                     ;
    auto &C0                      = this->C0                     ;

    // Pre-process the tracers by dividing by density inside the domain
    // After this, we can reconstruct tracers only (not rho * tracer)
    parallel_for( SimpleBounds<4>(num_tracers,nz,ny,nx) , YAKL_LAMBDA (int tr, int k, int j, int i) {
      tracers(tr,hs+k,hs+j,hs+i) /= (state(idR,hs+k,hs+j,hs+i) + hyDensCells(k));
    });

    // Populate the halos
    if        (bc_z == BC_PERIODIC) {
      parallel_for( SimpleBounds<3>(ny,nx,hs) , YAKL_LAMBDA(int j, int i, int kk) {
        for (int l=0; l < num_state; l++) {
          state(l,      kk,hs+j,hs+i) = state(l,nz+kk,hs+j,hs+i);
          state(l,hs+nz+kk,hs+j,hs+i) = state(l,hs+kk,hs+j,hs+i);
        }
        for (int l=0; l < num_tracers; l++) {
          tracers(l,      kk,hs+j,hs+i) = tracers(l,nz+kk,hs+j,hs+i);
          tracers(l,hs+nz+kk,hs+j,hs+i) = tracers(l,hs+kk,hs+j,hs+i);
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
        for (int l=0; l < num_tracers; l++) {
          tracers(l,      kk,hs+j,hs+i) = tracers(l,hs     ,hs+j,hs+i);
          tracers(l,hs+nz+kk,hs+j,hs+i) = tracers(l,hs+nz-1,hs+j,hs+i);
        }
      });
    }

    // Loop through all cells, reconstruct in x-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      // We need these to persist to evolve tracers with ADER
      SArray<real,2,nAder,ngll> r_DTs , rw_DTs;

      { // BEGIN: reconstruct, time-avg, and store state & state fluxes
        ////////////////////////////////////////////////////////////////
        // Reconstruct rho, u, v, w, theta
        ////////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> ru_DTs , rv_DTs , rt_DTs;
        {
          SArray<real,1,ord> stencil;

          // Density
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idR,k+kk,hs+j,hs+i); }
          reconstruct_gll_values( stencil , r_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
          for (int kk=0; kk < ngll; kk++) { r_DTs(0,kk) += hyDensGLL(k,kk); } // Add hydrostasis back on

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
        SArray<real,2,nAder,ngll> rwu_DTs , rwv_DTs , rww_DTs , rwt_DTs , rt_gamma_DTs;
        for (int kk=0; kk < ngll; kk++) {
          real r = r_DTs (0,kk);
          real u = ru_DTs(0,kk) / r;
          real v = rv_DTs(0,kk) / r;
          real w = rw_DTs(0,kk) / r;
          real t = rt_DTs(0,kk) / r;
          rwu_DTs    (0,kk) = r*w*u;
          rwv_DTs    (0,kk) = r*w*v;
          rww_DTs    (0,kk) = r*w*w;
          rwt_DTs    (0,kk) = r*w*t;
          rt_gamma_DTs(0,kk) = pow(r*t,gamma);
        }

        //////////////////////////////////////////
        // Compute time derivatives if necessary
        //////////////////////////////////////////
        if (nAder > 1) {
          diffTransformEulerConsZ( r_DTs , ru_DTs , rv_DTs , rw_DTs , rt_DTs , rwu_DTs , rwv_DTs , rww_DTs ,
                                   rwt_DTs , rt_gamma_DTs , derivMatrix , hyPressureGLL , C0 , gamma , k ,
                                   dz(k) , bc_z , nz );
        }

        //////////////////////////////////////////
        // Time average if necessary
        //////////////////////////////////////////
        // We can't alter density and momentum because they're needed for tracers later
        SArray<real,1,ngll> r_tavg, rw_tavg;
        if (timeAvg) {
          compute_timeAvg( r_DTs  , r_tavg  , dt );
          compute_timeAvg( ru_DTs           , dt );
          compute_timeAvg( rv_DTs           , dt );
          compute_timeAvg( rw_DTs , rw_tavg , dt );
          compute_timeAvg( rt_DTs           , dt );
          compute_timeAvg( rwu_DTs          , dt );
          compute_timeAvg( rwv_DTs          , dt );
          compute_timeAvg( rww_DTs          , dt );
          compute_timeAvg( rwt_DTs          , dt );
          compute_timeAvg( rt_gamma_DTs     , dt );
        } else {
          for (int ii=0; ii < ngll; ii++) {
            r_tavg (ii) = r_DTs (0,ii);
            rw_tavg(ii) = rw_DTs(0,ii);
          }
        }

        //////////////////////////////////////////
        // Store cell edge estimates of the state
        //////////////////////////////////////////
        // Left interface
        stateLimits(idR,1,k  ,j,i) = r_tavg  (0     );
        stateLimits(idU,1,k  ,j,i) = ru_DTs(0,0     );
        stateLimits(idV,1,k  ,j,i) = rv_DTs(0,0     );
        stateLimits(idW,1,k  ,j,i) = rw_tavg (0     );
        stateLimits(idT,1,k  ,j,i) = rt_DTs(0,0     );
        // Right interface       
        stateLimits(idR,0,k+1,j,i) = r_tavg  (ngll-1);
        stateLimits(idU,0,k+1,j,i) = ru_DTs(0,ngll-1);
        stateLimits(idV,0,k+1,j,i) = rv_DTs(0,ngll-1);
        stateLimits(idW,0,k+1,j,i) = rw_tavg (ngll-1);
        stateLimits(idT,0,k+1,j,i) = rt_DTs(0,ngll-1);

        ////////////////////////////////////////////
        // Assign gravity source term
        ////////////////////////////////////////////
        real ravg = 0;
        for (int kk=0; kk < ngll; kk++) {
          ravg += r_tavg(kk) * gllWts_ngll(kk);
        }
        stateTend(idR,k,j,i) = 0;
        stateTend(idU,k,j,i) = 0;
        stateTend(idV,k,j,i) = 0;
        stateTend(idW,k,j,i) = -GRAV*ravg;
        stateTend(idT,k,j,i) = 0;
      } // END: reconstruct, time-avg, and store state & state fluxes

      // r_DTs and rw_DTs still exist and are computed
      { // BEGIN: Reconstruct, time-average, and store tracer fluxes
        // Only process one tracer at a time to save on local memory / register requirements
        for (int tr=0; tr < num_tracers; tr++) {
          SArray<real,2,nAder,ngll> rt_DTs; // Density * tracer
          { // BEGIN: Reconstruct the tracer
            SArray<real,1,ord> stencil;
            for (int kk=0; kk < ord; kk++) { stencil(kk) = tracers(tr,k+kk,hs+j,hs+i); }
            reconstruct_gll_values( stencil , rt_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
            for (int kk=0; kk < ngll; kk++) { rt_DTs(0,kk) *= r_DTs(0,kk); }
            if (tracer_pos(tr)) {
              for (int kk=0; kk < ngll; kk++) { rt_DTs(0,kk) = max( 0._fp , rt_DTs(0,kk) ); }
            }
          } // END: Reconstruct the tracer

          // Compute the tracer flux
          SArray<real,2,nAder,ngll> rwt_DTs; // Density * wwind * tracer
          for (int kk=0; kk < ngll; kk++) {
            rwt_DTs(0,kk) = rt_DTs(0,kk) * rw_DTs(0,kk) / r_DTs(0,kk);
          }

          //////////////////////////////////////////
          // Compute time derivatives if necessary
          //////////////////////////////////////////
          if (nAder > 1) {
            diffTransformTracer( r_DTs , rw_DTs , rt_DTs , rwt_DTs , derivMatrix , dz(k) );
          }

          //////////////////////////////////////////
          // Time average if necessary
          //////////////////////////////////////////
          if (timeAvg) {
            compute_timeAvg( rt_DTs  , dt );
            compute_timeAvg( rwt_DTs , dt );
          }
          if (tracer_pos(tr)) {
            for (int kk=0; kk < ngll; kk++) { rt_DTs(0,kk) = max( 0._fp , rt_DTs(0,kk) ); }
          }
          if (bc_z == BC_WALL) {
            if (k == nz-1) rwt_DTs(0,ngll-1) = 0;
            if (k == 0   ) rwt_DTs(0,0     ) = 0;
          }

          ////////////////////////////////////////////////////////////
          // Store cell edge estimates of the tracer
          ////////////////////////////////////////////////////////////
          tracerLimits(tr,1,k  ,j,i) = rt_DTs (0,0     ); // Left interface
          tracerLimits(tr,0,k+1,j,i) = rt_DTs (0,ngll-1); // Right interface
        }
      } // END: Reconstruct, time-average, and store tracer fluxes

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
      for (int l = 0; l < num_tracers; l++) {
        if        (bc_z == BC_PERIODIC) {
          tracerLimits(l,0,0 ,j,i) = tracerLimits(l,0,nz,j,i);
          tracerLimits(l,1,nz,j,i) = tracerLimits(l,1,0 ,j,i);
        } else if (bc_z == BC_WALL    ) {
          tracerLimits(l,0,0 ,j,i) = tracerLimits(l,1,0 ,j,i);
          tracerLimits(l,1,nz,j,i) = tracerLimits(l,0,nz,j,i);
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the upwind fluxes
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<3>(nz+1,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      // Get left and right state
      real r_L = stateLimits(idR,0,k,j,i)    ;   real r_R = stateLimits(idR,1,k,j,i)    ;
      real u_L = stateLimits(idU,0,k,j,i)/r_L;   real u_R = stateLimits(idU,1,k,j,i)/r_R;
      real v_L = stateLimits(idV,0,k,j,i)/r_L;   real v_R = stateLimits(idV,1,k,j,i)/r_R;
      real w_L = stateLimits(idW,0,k,j,i)/r_L;   real w_R = stateLimits(idW,1,k,j,i)/r_R;
      real t_L = stateLimits(idT,0,k,j,i)/r_L;   real t_R = stateLimits(idT,1,k,j,i)/r_R;
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
      real q1_L = stateLimits(idR,0,k,j,i);   real q1_R = stateLimits(idR,1,k,j,i);
      real q2_L = stateLimits(idU,0,k,j,i);   real q2_R = stateLimits(idU,1,k,j,i);
      real q3_L = stateLimits(idV,0,k,j,i);   real q3_R = stateLimits(idV,1,k,j,i);
      real q4_L = stateLimits(idW,0,k,j,i);   real q4_R = stateLimits(idW,1,k,j,i);
      real q5_L = stateLimits(idT,0,k,j,i);   real q5_R = stateLimits(idT,1,k,j,i);
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

      stateFlux(idR,k,j,i) = q4;
      stateFlux(idU,k,j,i) = q4*q2/q1;
      stateFlux(idV,k,j,i) = q4*q3/q1;
      stateFlux(idW,k,j,i) = q4*q4/q1 + C0*pow(q5,gamma);
      stateFlux(idT,k,j,i) = q4*q5/q1;

      // COMPUTE UPWIND TRACER FLUXES
      // Handle it one tracer at a time
      for (int tr=0; tr < num_tracers; tr++) {
        if (w > 0) {
          tracerFlux(tr,k,j,i) = q4 * tracerLimits(tr,0,k,j,i) / r_L;
        } else {
          tracerFlux(tr,k,j,i) = q4 * tracerLimits(tr,1,k,j,i) / r_R;
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Limit the tracer fluxes for positivity
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<4>(num_tracers,nz+1,ny,nx) , YAKL_LAMBDA (int tr, int k, int j, int i) {
      real constexpr eps = 1.e-10;
      real w = 0.5_fp * ( stateLimits(idW,0,k,j,i) + stateLimits(idW,1,k,j,i) );
      // Solid wall BCs mean w == 0 at boundaries
      if (tracer_pos(tr)) {
        // Compute and apply the flux reduction factor of the upwind cell
        if      (w > 0) {
          int ind_k = k-1;
          // upwind is to the left of this interface
          real f1 = min( tracerFlux(tr,ind_k  ,j,i) , 0._fp );
          real f2 = max( tracerFlux(tr,ind_k+1,j,i) , 0._fp );
          real fluxOut = dt*(f2-f1)/dz(k);
          real dens = state(idR,hs+ind_k,hs+j,hs+i) + hyDensCells(ind_k);
          tracerFlux(tr,k,j,i) *= min( 1._fp , tracers(tr,hs+ind_k,hs+j,hs+i) * dens / (fluxOut + eps) );
        } else if (w < 0) {
          int ind_k = k;
          // upwind is to the right of this interface
          real f1 = min( tracerFlux(tr,ind_k  ,j,i) , 0._fp );
          real f2 = max( tracerFlux(tr,ind_k+1,j,i) , 0._fp );
          real fluxOut = dt*(f2-f1)/dz(k);
          real dens = state(idR,hs+ind_k,hs+j,hs+i) + hyDensCells(ind_k);
          tracerFlux(tr,k,j,i) *= min( 1._fp , tracers(tr,hs+ind_k,hs+j,hs+i) * dens / (fluxOut + eps) );
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA(int k, int j, int i) {
      for (int l=0; l < num_state; l++) {
        if (sim2d && l == idV) {
          stateTend(l,k,j,i) = 0;
        } else {
          stateTend(l,k,j,i) += - ( stateFlux(l,k+1,j,i) - stateFlux(l,k,j,i) ) / dz(k);
        }
      }
      for (int l=0; l < num_tracers; l++) {
        // Compute tracer tendency
        tracerTend(l,k,j,i) = - ( tracerFlux(l,k+1,j,i) - tracerFlux(l,k,j,i) ) / dz(k);
        // Multiply density back onto the tracers
        tracers(l,hs+k,hs+j,hs+i) *= (state(idR,hs+k,hs+j,hs+i) + hyDensCells(k));
      }
    });
  }



  template <class F> void applyStateTendencies( F const &applySingleTendency , int splitIndex ) {
    parallel_for( SimpleBounds<4>(num_state,nz,ny,nx) , YAKL_LAMBDA (int l, int k, int j, int i) {
      applySingleTendency({l,k,j,i});
    });
  }



  template <class F> void applyTracerTendencies( F const &applySingleTendency , int splitIndex ) {
    parallel_for( SimpleBounds<4>(num_tracers,nz,ny,nx) , YAKL_LAMBDA (int l, int k, int j, int i) {
      applySingleTendency({l,k,j,i});
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
      auto zint = dm.get<real,1>("vertical_interface_height");
      real1d zmid("zmid",nz);
      parallel_for( nz , YAKL_LAMBDA (int i) { zmid(i) = ( zint(i) + zint(i+1) ) / 2; });
      nc.write(zmid.createHostCopy(),"z",{"z"});
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

    real4d  state   = createStateArr ();
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
      real r =   state(idR,hs+k,hs+j,hs+i) + hyDensCells     (k);
      real t = ( state(idT,hs+k,hs+j,hs+i) + hyDensThetaCells(k) ) / r;
      data(k,j,i) = t - hyThetaCells(k);
    });
    nc.write1(data.createHostCopy(),"pot_temp_pert",{"z","y","x"},ulIndex,"t");
    // pressure'
    parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
      real r  = state(idR,hs+k,hs+j,hs+i) + hyDensCells(k);
      real rt = state(idT,hs+k,hs+j,hs+i) + hyDensThetaCells(k);
      real p  = C0*pow(rt,gamma);
      data(k,j,i) = p - hyPressureCells(k);
    });
    nc.write1(data.createHostCopy(),"pressure_pert",{"z","y","x"},ulIndex,"t");

    for (int tr=0; tr < num_tracers; tr++) {
      parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        real r = state(idR,hs+k,hs+j,hs+i) + hyDensCells(k);
        data(k,j,i) = tracers(tr,hs+k,hs+j,hs+i)/r;
      });
      nc.write1(data.createHostCopy(),std::string("tracer_")+tracer_name[tr],{"z","y","x"},ulIndex,"t");
    }

    micro.output(dm, nc, ulIndex);

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
      weno::compute_weno_coefs<ord>( wenoRecon , stencil , wenoCoefs , idl , sigma );
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



  YAKL_INLINE void diffTransformEulerConsX( SArray<real,2,nAder,ngll> &r  ,
                                            SArray<real,2,nAder,ngll> &ru ,
                                            SArray<real,2,nAder,ngll> &rv ,
                                            SArray<real,2,nAder,ngll> &rw ,
                                            SArray<real,2,nAder,ngll> &rt ,
                                            SArray<real,2,nAder,ngll> &ruu ,
                                            SArray<real,2,nAder,ngll> &ruv ,
                                            SArray<real,2,nAder,ngll> &ruw ,
                                            SArray<real,2,nAder,ngll> &rut ,
                                            SArray<real,2,nAder,ngll> &rt_gamma ,
                                            SArray<real,2,ngll,ngll> const &deriv ,
                                            real C0, real gamma, real dx ) {
    // zero out the non-linear DTs
    for (int kt=1; kt < nAder; kt++) {
      for (int ii=0; ii < ngll; ii++) {
        ruu     (kt,ii) = 0;
        ruv     (kt,ii) = 0;
        ruw     (kt,ii) = 0;
        rut     (kt,ii) = 0;
        rt_gamma(kt,ii) = 0;
      }
    }

    // Loop over the time derivatives
    for (int kt=0; kt<nAder-1; kt++) {
      // Compute the state at the next time level
      for (int ii=0; ii<ngll; ii++) {
        real df1_dx = 0;
        real df2_dx = 0;
        real df3_dx = 0;
        real df4_dx = 0;
        real df5_dx = 0;
        for (int s=0; s<ngll; s++) {
          df1_dx += deriv(s,ii) * ( ru (kt,s) );
          if (kt == 0) { df2_dx += deriv(s,ii) * ( ruu(kt,s) + C0*rt_gamma(kt,s)   ); }
          else         { df2_dx += deriv(s,ii) * ( ruu(kt,s) + C0*rt_gamma(kt,s)/2 ); }
          df3_dx += deriv(s,ii) * ( ruv(kt,s) );
          df4_dx += deriv(s,ii) * ( ruw(kt,s) );
          df5_dx += deriv(s,ii) * ( rut(kt,s) );
        }
        r (kt+1,ii) = -df1_dx/dx/(kt+1._fp);
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
          tot_ruu += ru(ir,ii) * ru(kt+1-ir,ii) - r(ir,ii) * ruu(kt+1-ir,ii);
          tot_ruv += ru(ir,ii) * rv(kt+1-ir,ii) - r(ir,ii) * ruv(kt+1-ir,ii);
          tot_ruw += ru(ir,ii) * rw(kt+1-ir,ii) - r(ir,ii) * ruw(kt+1-ir,ii);
          tot_rut += ru(ir,ii) * rt(kt+1-ir,ii) - r(ir,ii) * rut(kt+1-ir,ii);
        }
        ruu(kt+1,ii) = tot_ruu / r(0,ii);
        ruv(kt+1,ii) = tot_ruv / r(0,ii);
        ruw(kt+1,ii) = tot_ruw / r(0,ii);
        rut(kt+1,ii) = tot_rut / r(0,ii);

        // Compute rt_gamma at the next time level
        real tot_rt_gamma = 0;
        for (int ir=0; ir<=kt; ir++) {
          tot_rt_gamma += (kt+1._fp -ir) * ( gamma*rt_gamma(ir,ii)*rt(kt+1-ir,ii) - rt(ir,ii)*rt_gamma(kt+1-ir,ii) );
        }
        rt_gamma(kt+1,ii) = ( gamma*rt_gamma(0,ii)*rt(kt+1,ii) + tot_rt_gamma / (kt+1._fp) ) / rt(0,ii);
      }
    }

    // Fix the rt_gamma
    for (int kt=1; kt < nAder; kt++) {
      for (int ii=0; ii < ngll; ii++) {
        rt_gamma(kt,ii) /= 2;
      }
    }
  }



  YAKL_INLINE void diffTransformEulerConsY( SArray<real,2,nAder,ngll> &r  ,
                                            SArray<real,2,nAder,ngll> &ru ,
                                            SArray<real,2,nAder,ngll> &rv ,
                                            SArray<real,2,nAder,ngll> &rw ,
                                            SArray<real,2,nAder,ngll> &rt ,
                                            SArray<real,2,nAder,ngll> &rvu ,
                                            SArray<real,2,nAder,ngll> &rvv ,
                                            SArray<real,2,nAder,ngll> &rvw ,
                                            SArray<real,2,nAder,ngll> &rvt ,
                                            SArray<real,2,nAder,ngll> &rt_gamma ,
                                            SArray<real,2,ngll,ngll> const &deriv , 
                                            real C0, real gamma, real dy ) {
    // zero out the non-linear DTs
    for (int kt=1; kt < nAder; kt++) {
      for (int ii=0; ii < ngll; ii++) {
        rvu     (kt,ii) = 0;
        rvv     (kt,ii) = 0;
        rvw     (kt,ii) = 0;
        rvt     (kt,ii) = 0;
        rt_gamma(kt,ii) = 0;
      }
    }

    // Loop over the time derivatives
    for (int kt=0; kt<nAder-1; kt++) {
      // Compute the state at the next time level
      for (int ii=0; ii<ngll; ii++) {
        real drv_dy    = 0;
        real drvu_dy   = 0;
        real drvv_p_dy = 0;
        real drvw_dy   = 0;
        real drvt_dy   = 0;
        for (int s=0; s<ngll; s++) {
          drv_dy    += deriv(s,ii) * rv(kt,s);
          drvu_dy   += deriv(s,ii) * rvu(kt,s);
          if (kt == 0) { drvv_p_dy += deriv(s,ii) * ( rvv(kt,s) + C0*rt_gamma(kt,s)   ); }
          else         { drvv_p_dy += deriv(s,ii) * ( rvv(kt,s) + C0*rt_gamma(kt,s)/2 ); }
          drvw_dy   += deriv(s,ii) * rvw(kt,s);
          drvt_dy   += deriv(s,ii) * rvt(kt,s);
        }
        r (kt+1,ii) = -drv_dy   /dy/(kt+1);
        ru(kt+1,ii) = -drvu_dy  /dy/(kt+1);
        rv(kt+1,ii) = -drvv_p_dy/dy/(kt+1);
        rw(kt+1,ii) = -drvw_dy  /dy/(kt+1);
        rt(kt+1,ii) = -drvt_dy  /dy/(kt+1);
      }

      // Compute ru* at the next time level
      for (int ii=0; ii<ngll; ii++) {
        // Compute the non-linear differential transforms
        real tot_rvu = 0;
        real tot_rvv = 0;
        real tot_rvw = 0;
        real tot_rvt = 0;
        for (int l=0; l<=kt+1; l++) {
          tot_rvu += rv(l,ii) * ru(kt+1-l,ii) - r(l,ii) * rvu(kt+1-l,ii);
          tot_rvv += rv(l,ii) * rv(kt+1-l,ii) - r(l,ii) * rvv(kt+1-l,ii);
          tot_rvw += rv(l,ii) * rw(kt+1-l,ii) - r(l,ii) * rvw(kt+1-l,ii);
          tot_rvt += rv(l,ii) * rt(kt+1-l,ii) - r(l,ii) * rvt(kt+1-l,ii);
        }
        rvu(kt+1,ii) = tot_rvu / r(0,ii);
        rvv(kt+1,ii) = tot_rvv / r(0,ii);
        rvw(kt+1,ii) = tot_rvw / r(0,ii);
        rvt(kt+1,ii) = tot_rvt / r(0,ii);

        // Compute rt_gamma at the next time level
        real tot_rt_gamma = 0;
        for (int l=0; l<=kt; l++) {
          tot_rt_gamma += (kt+1._fp -l) * ( gamma*rt_gamma(l,ii)*rt(kt+1-l,ii) - rt(l,ii)*rt_gamma(kt+1-l,ii) );
        }
        rt_gamma(kt+1,ii) = ( gamma*rt_gamma(0,ii)*rt(kt+1,ii) + tot_rt_gamma / (kt+1._fp) ) / rt(0,ii);
      }
    }

    // Fix the rt_gamma
    for (int kt=1; kt < nAder; kt++) {
      for (int ii=0; ii < ngll; ii++) {
        rt_gamma(kt,ii) /= 2;
      }
    }
  }



  YAKL_INLINE void diffTransformEulerConsZ( SArray<real,2,nAder,ngll> &r  ,
                                            SArray<real,2,nAder,ngll> &ru ,
                                            SArray<real,2,nAder,ngll> &rv ,
                                            SArray<real,2,nAder,ngll> &rw ,
                                            SArray<real,2,nAder,ngll> &rt ,
                                            SArray<real,2,nAder,ngll> &rwu ,
                                            SArray<real,2,nAder,ngll> &rwv ,
                                            SArray<real,2,nAder,ngll> &rww ,
                                            SArray<real,2,nAder,ngll> &rwt ,
                                            SArray<real,2,nAder,ngll> &rt_gamma ,
                                            SArray<real,2,ngll,ngll> const &deriv , 
                                            real2d const &hyPressureGLL , 
                                            real C0, real gamma ,
                                            int k , real dz , int bc_z , int nz ) {
    // zero out the non-linear DTs
    for (int kt=1; kt < nAder; kt++) {
      for (int ii=0; ii < ngll; ii++) {
        rwu     (kt,ii) = 0;
        rwv     (kt,ii) = 0;
        rww     (kt,ii) = 0;
        rwt     (kt,ii) = 0;
        rt_gamma(kt,ii) = 0;
      }
    }

    // Loop over the time derivatives
    for (int kt=0; kt<nAder-1; kt++) {
      // Compute the state at the next time level
      for (int ii=0; ii<ngll; ii++) {
        real drw_dz    = 0;
        real drwu_dz   = 0;
        real drwv_dz   = 0;
        real drww_p_dz = 0;
        real drwt_dz   = 0;
        for (int s=0; s<ngll; s++) {
          drw_dz    += deriv(s,ii) * rw(kt,s);
          drwu_dz   += deriv(s,ii) * rwu(kt,s);
          drwv_dz   += deriv(s,ii) * rwv(kt,s);
          if (kt == 0) { drww_p_dz += deriv(s,ii) * ( rww(kt,s) + C0*rt_gamma(kt,s) - hyPressureGLL(k,s) ); }
          else         { drww_p_dz += deriv(s,ii) * ( rww(kt,s) + C0*rt_gamma(kt,s)/2                    ); }
          drwt_dz   += deriv(s,ii) * rwt(kt,s);
        }
        r (kt+1,ii) = -drw_dz   /dz/(kt+1);
        ru(kt+1,ii) = -drwu_dz  /dz/(kt+1);
        rv(kt+1,ii) = -drwv_dz  /dz/(kt+1);
        rw(kt+1,ii) = -drww_p_dz/dz/(kt+1);
        rt(kt+1,ii) = -drwt_dz  /dz/(kt+1);
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
          tot_rwu += rw(l,ii) * ru(kt+1-l,ii) - r(l,ii) * rwu(kt+1-l,ii);
          tot_rwv += rw(l,ii) * rv(kt+1-l,ii) - r(l,ii) * rwv(kt+1-l,ii);
          tot_rww += rw(l,ii) * rw(kt+1-l,ii) - r(l,ii) * rww(kt+1-l,ii);
          tot_rwt += rw(l,ii) * rt(kt+1-l,ii) - r(l,ii) * rwt(kt+1-l,ii);
        }
        rwu(kt+1,ii) = tot_rwu / r(0,ii);
        rwv(kt+1,ii) = tot_rwv / r(0,ii);
        rww(kt+1,ii) = tot_rww / r(0,ii);
        rwt(kt+1,ii) = tot_rwt / r(0,ii);

        // Compute rt_gamma at the next time level
        real tot_rt_gamma = 0;
        for (int l=0; l<=kt; l++) {
          tot_rt_gamma += (kt+1-l) * ( gamma*rt_gamma(l,ii)*rt(kt+1-l,ii) - rt(l,ii)*rt_gamma(kt+1-l,ii) );
        }
        rt_gamma(kt+1,ii) = ( gamma*rt_gamma(0,ii)*rt(kt+1,ii) + tot_rt_gamma / (kt+1) ) / rt(0,ii);
      }
    }

    // Fix the rt_gamma
    for (int kt=1; kt < nAder; kt++) {
      for (int ii=0; ii < ngll; ii++) {
        rt_gamma(kt,ii) /= 2;
      }
    }
  }



  YAKL_INLINE void diffTransformTracer( SArray<real,2,nAder,ngll> const &r  ,
                                        SArray<real,2,nAder,ngll> const &ru ,
                                        SArray<real,2,nAder,ngll> &rt ,
                                        SArray<real,2,nAder,ngll> &rut ,
                                        SArray<real,2,ngll,ngll> const &deriv , 
                                        real dx ) {
    // zero out the non-linear DT
    for (int kt=1; kt < nAder; kt++) {
      for (int ii=0; ii < ngll; ii++) {
        rut(kt,ii) = 0;
      }
    }
    // Loop over the time derivatives
    for (int kt=0; kt<nAder-1; kt++) {
      // Compute the rho*tracer at the next time level
      for (int ii=0; ii<ngll; ii++) {
        real df_dx = 0;
        for (int s=0; s<ngll; s++) {
          df_dx += deriv(s,ii) * rut(kt,s);
        }
        rt(kt+1,ii) = -df_dx/dx/(kt+1._fp);
      }
      // Compute rut at the next time level
      for (int ii=0; ii<ngll; ii++) {
        real tot_rut = 0;
        for (int ir=0; ir<=kt+1; ir++) {
          tot_rut += ru(ir,ii) * rt(kt+1-ir,ii) - r(ir,ii) * rut(kt+1-ir,ii);
        }
        rut(kt+1,ii) = tot_rut / r(0,ii);
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


