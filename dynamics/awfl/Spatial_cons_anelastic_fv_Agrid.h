
#pragma once

#include "awfl_const.h"
#include "phys_params.h"
#include "TransformMatrices.h"
#include "TransformMatrices_variable.h"
#include "WenoLimiter.h"
#include "Profiles.h"
#include "DataManager.h"
#include "pam_coupler.h"
#include <Eigen/Sparse>


template <int nTimeDerivs, bool timeAvg, int nAder>
class Spatial_operator {
public:

  static_assert(nTimeDerivs == 1 , "ERROR: This Spatial class isn't setup to use nTimeDerivs > 1");

  int static constexpr hs = (ord-1)/2;
  int static constexpr num_state = 5;
  int static constexpr max_tracers = 50;

  real Rd   ;
  real cp   ;
  real gamma;
  real p0   ;
  real C0   ;
  real Rv   ;

  real hydrostasis_parameters_sum;

  typedef real5d StateArr;  // Array of state variables (rho, rho*u, rho*v, rho*w, and rho*theta)
  typedef real5d TracerArr; // Array of tracers (total tracer mass)

  typedef real5d StateTendArr;   // State tendencies
  typedef real5d TracerTendArr;  // Tracer tendencies

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
  SArray<real,2,ord,ord > sten_to_coefs;
  SArray<real,2,ord,ngll> sten_to_deriv_gll;
  // WENO reconstruction matrices
  SArray<real,3,hs+1,hs+1,hs+1> weno_recon_lower;
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
  real4d vert_sten_to_coefs;
  real5d vert_weno_recon_lower;

  Eigen::SparseLU<Eigen::SparseMatrix<real>, Eigen::COLAMDOrdering<int> > momdiv_solver;

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
  int static constexpr DATA_SPEC_COLLISION     = 3;

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
  real1d      zbot;
  real1d      ztop;
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
    real5d state       = dm.get<real,5>( "dynamics_state"   );
    real5d tracers     = dm.get<real,5>( "dynamics_tracers" );
    real4d dm_dens_dry = dm.get<real,4>( "density_dry"      );
    real4d dm_uvel     = dm.get<real,4>( "uvel"             );
    real4d dm_vvel     = dm.get<real,4>( "vvel"             );
    real4d dm_wvel     = dm.get<real,4>( "wvel"             );
    real4d dm_temp     = dm.get<real,4>( "temp"             );

    YAKL_SCOPE( hyDensCells      , this->hyDensCells      );
    YAKL_SCOPE( hyDensThetaCells , this->hyDensThetaCells );
    YAKL_SCOPE( C0               , this->C0               );
    YAKL_SCOPE( gamma            , this->gamma            );
    YAKL_SCOPE( num_tracers      , this->num_tracers      );
    YAKL_SCOPE( Rd               , this->Rd               );
    YAKL_SCOPE( Rv               , this->Rv               );
    YAKL_SCOPE( tracer_adds_mass , this->tracer_adds_mass );

    int idWV = micro.get_water_vapor_index();

    MultipleFields<max_tracers,real4d> dm_tracers;
    for (int tr = 0; tr < num_tracers; tr++) {
      auto trac = dm.get<real,4>( tracer_name[tr] );
      dm_tracers.add_field( trac );
    }

    parallel_for( "Spatial.h d2c" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      // Compute wind velocities
      real dens  = hyDensCells(k,iens);
      real uvel  = state(idU,hs+k,hs+j,hs+i,iens) / dens;
      real vvel  = state(idV,hs+k,hs+j,hs+i,iens) / dens;
      real wvel  = state(idW,hs+k,hs+j,hs+i,iens) / dens;

      // Compute implied full density
      real theta = ( state(idT,hs+k,hs+j,hs+i,iens) + hyDensThetaCells(k,iens) ) / dens;
      dens = pow( hyPressureCells(k,iens) / C0 , 1._fp/gamma ) / theta;

      // Compute implied dry density
      real dens_dry = dens;
      for (int tr=0; tr < num_tracers; tr++) {
        if (tracer_adds_mass(tr)) dens_dry -= tracers(tr,hs+k,hs+j,hs+i,iens);
      }

      // Compute temperature
      real dens_vap = tracers(idWV,hs+k,hs+j,hs+i,iens);
      real temp = hyPressureCells(k,iens) / (dens_dry*Rd + dens_vap*Rv);

      // Store into the coupled state
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



  YAKL_INLINE static real hydrostatic_dens_theta( real2d const &hy_params , real z , real zbot , real ztop , int iens , real C0 , real gamma ) {
    real p = pam::hydrostatic_pressure( hy_params , z , zbot , ztop , iens );
    // p = C0*(rho*theta)^gamma
    return pow(p/C0,1._fp/gamma);
  }



  YAKL_INLINE static real hydrostatic_theta( real2d const &hy_params , real z , real zbot , real ztop , int iens , real C0 , real gamma , real grav ) {
    real rt =      hydrostatic_dens_theta( hy_params , z , zbot , ztop , iens , C0 , gamma        );
    real r  = pam::hydrostatic_density   ( hy_params , z , zbot , ztop , iens              , grav );
    return rt/r;
  }



  template <class MICRO>
  void convert_coupler_state_to_dynamics( DataManager &dm , MICRO &micro ) {
    // auto hy_params = dm.get<real,2>("hydrostasis_parameters");

    YAKL_SCOPE( hyPressureCells  , this->hyPressureCells  );
    YAKL_SCOPE( hyThetaCells     , this->hyThetaCells     );
    YAKL_SCOPE( hyDensCells      , this->hyDensCells      );
    YAKL_SCOPE( hyDensThetaCells , this->hyDensThetaCells );
    YAKL_SCOPE( hyPressureGLL    , this->hyPressureGLL    );
    YAKL_SCOPE( hyThetaGLL       , this->hyThetaGLL       );
    YAKL_SCOPE( hyDensGLL        , this->hyDensGLL        );
    YAKL_SCOPE( hyDensThetaGLL   , this->hyDensThetaGLL   );
    YAKL_SCOPE( C0               , this->C0               );
    YAKL_SCOPE( gamma            , this->gamma            );
    YAKL_SCOPE( num_tracers      , this->num_tracers      );
    YAKL_SCOPE( Rd               , this->Rd               );
    YAKL_SCOPE( Rv               , this->Rv               );
    YAKL_SCOPE( tracer_adds_mass , this->tracer_adds_mass );
    YAKL_SCOPE( zbot             , this->zbot             );
    YAKL_SCOPE( ztop             , this->ztop             );
    YAKL_SCOPE( gllPts_ngll      , this->gllPts_ngll      );
    YAKL_SCOPE( dz               , this->dz               );
    YAKL_SCOPE( vert_interface   , this->vert_interface   );

    real5d state       = dm.get<real,5>( "dynamics_state"   );
    real5d tracers     = dm.get<real,5>( "dynamics_tracers" );
    real4d dm_dens_dry = dm.get<real,4>( "density_dry"      );
    real4d dm_uvel     = dm.get<real,4>( "uvel"             );
    real4d dm_vvel     = dm.get<real,4>( "vvel"             );
    real4d dm_wvel     = dm.get<real,4>( "wvel"             );
    real4d dm_temp     = dm.get<real,4>( "temp"             );

    int idWV = micro.get_water_vapor_index();

    MultipleFields<max_tracers,real4d> dm_tracers;
    for (int tr = 0; tr < num_tracers; tr++) {
      auto trac = dm.get<real,4>( tracer_name[tr] );
      dm_tracers.add_field( trac );
    }

    // // If hydrostasis in the coupler has changed, then we need to re-compute
    // // hydrostatically balanced cells and GLL points for the dycore's time step
    // real tmp = yakl::intrinsics::sum(hy_params);
    // if (tmp != hydrostasis_parameters_sum) {
    //   SArray<real,1,9> gll_pts, gll_wts;
    //   TransformMatrices::get_gll_points ( gll_pts );
    //   TransformMatrices::get_gll_weights( gll_wts );

    //   real2d hyPressureCells_prev ("hyPressureCells_prev ",nz,nens);
    //   real2d hyDensCells_prev     ("hyDensCells_prev     ",nz,nens);
    //   real2d hyDensThetaCells_prev("hyDensThetaCells_prev",nz,nens);
    //   real2d hyThetaCells_prev    ("hyThetaCells_prev    ",nz,nens);
    //   hyPressureCells .deep_copy_to(hyPressureCells_prev );
    //   hyDensCells     .deep_copy_to(hyDensCells_prev     );
    //   hyDensThetaCells.deep_copy_to(hyDensThetaCells_prev);
    //   hyThetaCells    .deep_copy_to(hyThetaCells_prev    );

    //   // Compute new cell averages and GLL point values for hydrostasis
    //   hydrostasis_parameters_sum = tmp;
    //   parallel_for( "Spatial.h new hydrostasis" , SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
    //     real p  = 0;
    //     real r  = 0;
    //     real rt = 0;
    //     real t  = 0;
    //     for (int kk=0; kk < 9; kk++) {
    //       real zloc = vert_interface(k,iens) + 0.5_fp*dz(k,iens) + gll_pts(kk)*dz(k,iens);
    //       real wt = gll_wts(kk);
    //       p  += pam::hydrostatic_pressure  ( hy_params , zloc , zbot(iens) , ztop(iens) , iens                     ) * wt;
    //       r  += pam::hydrostatic_density   ( hy_params , zloc , zbot(iens) , ztop(iens) , iens              , GRAV ) * wt;
    //       rt +=      hydrostatic_dens_theta( hy_params , zloc , zbot(iens) , ztop(iens) , iens , C0 , gamma        ) * wt;
    //       t  +=      hydrostatic_theta     ( hy_params , zloc , zbot(iens) , ztop(iens) , iens , C0 , gamma , GRAV ) * wt;
    //     }
    //     hyPressureCells (k,iens) = p;
    //     hyDensCells     (k,iens) = r;
    //     hyDensThetaCells(k,iens) = rt;
    //     hyThetaCells    (k,iens) = t;

    //     for (int kk=0; kk < ngll; kk++) {
    //       real zloc = vert_interface(k,iens) + 0.5_fp*dz(k,iens) + gllPts_ngll(kk)*dz(k,iens);
    //       hyPressureGLL (k,kk,iens) = pam::hydrostatic_pressure  ( hy_params , zloc , zbot(iens) , ztop(iens) , iens                     );
    //       hyDensGLL     (k,kk,iens) = pam::hydrostatic_density   ( hy_params , zloc , zbot(iens) , ztop(iens) , iens              , GRAV );
    //       hyDensThetaGLL(k,kk,iens) =      hydrostatic_dens_theta( hy_params , zloc , zbot(iens) , ztop(iens) , iens , C0 , gamma        );
    //       hyThetaGLL    (k,kk,iens) =      hydrostatic_theta     ( hy_params , zloc , zbot(iens) , ztop(iens) , iens , C0 , gamma , GRAV );
    //     }
    //   });
    // }

    parallel_for(  "Spatial.h c2d" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      // We're going to ignore dry density and compute it from temperature and tracer densities instead
      state(idR,hs+k,hs+j,hs+i,iens) = 0;

      // Load the tracers
      for (int tr=0; tr < num_tracers; tr++) {
        tracers(tr,hs+k,hs+j,hs+i,iens) = dm_tracers(tr,k,j,i,iens);
      }

      // Compute momenta
      real uvel = dm_uvel(k,j,i,iens);
      real vvel = dm_vvel(k,j,i,iens);
      real wvel = dm_wvel(k,j,i,iens);
      state(idU,hs+k,hs+j,hs+i,iens) = hyDensCells(k,iens) * uvel;
      state(idV,hs+k,hs+j,hs+i,iens) = hyDensCells(k,iens) * vvel;
      state(idW,hs+k,hs+j,hs+i,iens) = hyDensCells(k,iens) * wvel;

      // Compute dry density, then full density, then potential temperature
      real temp     = dm_temp(k,j,i,iens);
      real dens_vap = tracers(idWV,hs+k,hs+j,hs+i,iens);
      real dens_dry = (hyPressureCells(k,iens) - dens_vap*Rv*temp) / (Rd*temp);
      real dens     = dens_dry;
      for (int tr=0; tr < num_tracers; tr++) {
        if (tracer_adds_mass(tr)) dens += tracers(tr,hs+k,hs+j,hs+i,iens);
      }
      real theta = pow( hyPressureCells(k,iens) / C0 , 1._fp / gamma ) / dens;
      state(idT,hs+k,hs+j,hs+i,iens) = theta * hyDensCells(k,iens);
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
    parallel_for( "Spatial.h add_tracer" , 1 , YAKL_LAMBDA (int i) {
      tracer_pos      (tr) = pos_def;   // Store whether it's positive-definite
      tracer_adds_mass(tr) = adds_mass; // Store whether it adds mass (otherwise it's passive)
    });

    // Return the index of this tracer to the caller
    return tr;
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
  real compute_time_step(DataManager &dm, MICRO const &micro, real cfl_in = -1) {
    real cfl = cfl_in;
    if (cfl < 0) cfl = 0.8;

    // If we've already computed the time step, then don't compute it again
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
    parallel_for( "Spatial.h compute_time_step" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      // Get the state
      real r = hyDensCells(k,iens);
      real u = state(idU,hs+k,hs+j,hs+i,iens) / r;
      real v = state(idV,hs+k,hs+j,hs+i,iens) / r;
      real w = state(idW,hs+k,hs+j,hs+i,iens) / r;

      // Compute the maximum stable time step in each direction
      real udt = cfl * dx         / ( abs(u) + 1.e-20 );
      real vdt = cfl * dy         / ( abs(v) + 1.e-20 );
      real wdt = cfl * dz(k,iens) / ( abs(w) + 1.e-20 );

      // Compute the min of the max stable time steps
      dt3d(k,j,i,iens) = min( 20._fp , min( min(udt,vdt) , wdt ) );
    });
    // Store to dtInit so we don't have to compute this again
    return yakl::intrinsics::minval( dt3d );
  }



  // Initialize crap needed by recon()
  void init(std::string inFile, int ny, int nx, int nens, real xlen, real ylen, int num_tracers, DataManager &dm) {
    using yakl::intrinsics::matmul_cr;

    this->nens = nens;
    this->nx = nx;
    this->ny = ny;
    this->xlen = xlen;
    this->ylen = ylen;
    this->num_tracers = num_tracers;

    this->hydrostasis_parameters_sum = 0;

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
      } else if (dataStr == "collision") {
        data_spec = DATA_SPEC_COLLISION;
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
    zbot = real1d("zbot",nens);
    ztop = real1d("ztop",nens);
    YAKL_SCOPE( zbot , this->zbot );
    YAKL_SCOPE( ztop , this->ztop );
    YAKL_SCOPE( nz   , this->nz   );
    parallel_for( "Spatial.h init 1" , nens , YAKL_LAMBDA (int iens) {
      zbot(iens) = zint(0 ,iens);
      ztop(iens) = zint(nz,iens);
    });

    vert_interface        = real2d("vert_interface"      ,nz+1          ,nens);
    vert_interface_ghost  = real2d("vert_interface_ghost",nz+2*hs+1     ,nens);
    vert_locs_normalized  = real3d("vert_locs_normalized",nz,ord+1      ,nens);
    dz                    = real2d("dz"                  ,nz            ,nens);
    dz_ghost              = real2d("dz_ghost"            ,nz+2*hs       ,nens);
    vert_sten_to_gll      = real4d("vert_sten_to_gll"     ,nz,ord,ngll,nens);
    vert_sten_to_coefs    = real4d("vert_sten_to_coefs"   ,nz,ord,ord ,nens);
    vert_weno_recon_lower = real5d("vert_weno_recon_lower",nz,hs+1,hs+1,hs+1,nens);

    YAKL_SCOPE( vert_interface        , this->vert_interface        );
    YAKL_SCOPE( vert_interface_ghost  , this->vert_interface_ghost  );
    YAKL_SCOPE( vert_locs_normalized  , this->vert_locs_normalized  );
    YAKL_SCOPE( dz                    , this->dz                    );
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
    dx = xlen/nx;
    dy = ylen/ny;

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
      this->sten_to_deriv_gll  = matmul_cr( c2g_lower , matmul_cr( c2d , s2c ) );
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
    parallel_for( "Spatial.h init 4" , SimpleBounds<4>(nz+2*hs,ny+2*hs,nx+2*hs,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      for (int l=0; l < num_state; l++) {
        state(l,k,j,i,iens) = 0;
      }
      for (int tr=0; tr < num_tracers; tr++) {
        tracers(tr,k,j,i,iens) = 0;
      }
    });

    if (sim2d) { build_matrix_solver_2d(); }
    else       { build_matrix_solver_3d(); }

    #ifdef PAM_STANDALONE
      std::cout << "nx: " << nx << "\n";
      std::cout << "ny: " << ny << "\n";
      std::cout << "nz: " << nz << "\n";
      std::cout << "xlen (m): " << xlen << "\n";
      std::cout << "ylen (m): " << ylen << "\n";
      std::cout << "zbot (m): " << zbot.createHostCopy()(0) << "\n";
      std::cout << "ztop (m): " << ztop.createHostCopy()(0) << "\n";
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
    YAKL_SCOPE( num_tracers              , this->num_tracers             );
    YAKL_SCOPE( tracer_adds_mass         , this->tracer_adds_mass        );

    real5d state   = dm.get<real,5>("dynamics_state"  );
    real5d tracers = dm.get<real,5>("dynamics_tracers");

    // If the data_spec is thermal or ..., then initialize the domain with Exner pressure-based hydrostasis
    // This is mostly to make plotting potential temperature perturbation easier for publications
    if (data_spec == DATA_SPEC_THERMAL || data_spec == DATA_SPEC_COLLISION) {

      // Setup hydrostatic background state
      parallel_for( "Spatial.h init_state 1" , SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
        // Compute cell averages
        hyDensCells     (k,iens) = 0;
        hyPressureCells (k,iens) = 0;
        hyThetaCells    (k,iens) = 0;
        hyDensThetaCells(k,iens) = 0;
        for (int kk=0; kk<ord; kk++) {
          real zloc = vert_interface(k,iens) + 0.5_fp*dz(k,iens) + gllPts_ord(kk)*dz(k,iens);
          if        (data_spec == DATA_SPEC_THERMAL || data_spec == DATA_SPEC_COLLISION) {
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

      parallel_for( "Spatial.h init_state 2" , SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
        // Compute ngll GLL points
        for (int kk=0; kk<ngll; kk++) {
          real zloc = vert_interface(k,iens) + 0.5_fp*dz(k,iens) + gllPts_ngll(kk)*dz(k,iens);
          if        (data_spec == DATA_SPEC_THERMAL || data_spec == DATA_SPEC_COLLISION) {
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
      parallel_for( "Spatial.h init_state 3" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
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

                state(idT,hs+k,hs+j,hs+i,iens) += (r*t - rh*th) * wt;
              } else if (data_spec == DATA_SPEC_COLLISION) {
                // Compute constant theta hydrostatic background state
                real th = 300;
                real rh = profiles::initConstTheta_density(th,zloc,Rd,cp,gamma,p0,C0);
                real tp = profiles::ellipsoid_linear(xloc, yloc, zloc, xlen/2, ylen/2, 2000, 2000, 2000, 2000,  20 ) +
                          profiles::ellipsoid_linear(xloc, yloc, zloc, xlen/2, ylen/2, 8000, 2000, 2000, 2000, -20 );
                real t = th + tp;
                real r = rh;

                state(idT,hs+k,hs+j,hs+i,iens) += (r*t - rh*th) * wt;
              }
            }
          }
        }
      });

      int idWV = micro.get_water_vapor_index();

      parallel_for( "Spatial.h init_tracers" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        tracers(idWV,hs+k,hs+j,hs+i,iens) = 0;
      });

    } // if (data_spec == DATA_SPEC_THERMAL || data_spec == DATA_SPEC_COLLISION)


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

      YAKL_SCOPE( ztop , this->ztop );

      // Compute full density at ord GLL points for the space between each cell
      parallel_for( "Spatial.h init_state 4" , SimpleBounds<4>(nz,ngll-1,ord,nens) , YAKL_LAMBDA (int k, int kk, int kkk, int iens) {
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
        real temp      = profiles::init_supercell_temperature (zloc, z_0, z_trop, ztop(iens), T_0, T_trop, T_top);
        real press_dry = profiles::init_supercell_pressure_dry(zloc, z_0, z_trop, ztop(iens), T_0, T_trop, T_top, p_0, Rd);
        real dens_dry  = press_dry / (Rd*temp);
        real qvs       = profiles::init_supercell_sat_mix_dry(press_dry, temp);
        real relhum    = profiles::init_supercell_relhum(zloc, z_0, z_trop);
        if (relhum * qvs > 0.014_fp) relhum = 0.014_fp / qvs;
        real qv        = min( 0.014_fp , qvs*relhum );
        quad_temp(k,kk,kkk,iens) = -(1+qv)*GRAV/(Rd+qv*Rv)/temp;
      });

      parallel_for( "Spatial.h init_state 5" , nens , YAKL_LAMBDA (int iens) {
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

      parallel_for( "Spatial.h init_state 6" , SimpleBounds<3>(nz,ngll,nens) , YAKL_LAMBDA (int k, int kk, int iens) {
        real zloc = vert_interface(k,iens) + 0.5_fp*dz(k,iens) + gllPts_ngll(kk)*dz(k,iens);
        real temp       = profiles::init_supercell_temperature (zloc, z_0, z_trop, ztop(iens), T_0, T_trop, T_top);
        real press_tmp  = profiles::init_supercell_pressure_dry(zloc, z_0, z_trop, ztop(iens), T_0, T_trop, T_top, p_0, Rd);
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

      parallel_for( "Spatial.h init_state 7" , SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
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
        real press_tmp  = profiles::init_supercell_pressure_dry(zloc, z_0, z_trop, ztop(iens), T_0, T_trop, T_top, p_0, Rd);
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
        parallel_for( "Spatial.h init_state 8" , nz , YAKL_LAMBDA (int k) { data(k) = z(k,iens); });
        nc.write(data.createHostCopy(),"z"          ,{"z"});

        parallel_for( "Spatial.h init_state 9" , nz , YAKL_LAMBDA (int k) { data(k) = hyPressureCells(k,iens); });
        nc.write(data.createHostCopy(),"pressure"   ,{"z"});

        parallel_for( "Spatial.h init_state 10" , nz , YAKL_LAMBDA (int k) { data(k) = temp_hy(k,iens); });
        nc.write(data.createHostCopy(),"temperature",{"z"});

        parallel_for( "Spatial.h init_state 11" , nz , YAKL_LAMBDA (int k) { data(k) = tdew_hy(k,iens); });
        nc.write(data.createHostCopy(),"dew_point"  ,{"z"});
      }
      nc.close();

      int idWV = micro.get_water_vapor_index();
      real5d tracers = dm.get<real,5>("dynamics_tracers");

      parallel_for( "Spatial.h init_state 12" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
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



  void build_matrix_solver_2d() {
    int constexpr id_p  = 0;
    int constexpr id_ru = 1;
    int constexpr id_rw = 2;
    int constexpr neq = 3;

    Eigen::SparseMatrix<real> A(neq*nx*nz,neq*nx*nz);
    A.reserve(Eigen::VectorXi::Constant(neq*nx*nz,5));

    for (int k=0; k < nz; k++) {
      for (int i=0; i < nx; i++) {
        int im1 = i-1;  if (im1 < 0   ) im1 += nx;
        int ip1 = i+1;  if (ip1 > nx-1) ip1 -= nx;
        int km1 = k-1;  if (km1 < 0   ) km1  = 0;
        int kp1 = k+1;  if (kp1 > nz-1) kp1  = nz-1;
        real p_km1=0, p_k=0, p_kp1=0, ru_im1=0, ru_i=0, ru_ip1=0, rw_km1=0, rw_k=0, rw_kp1=0;
        real rdx = 1._fp / dx;
        real rdz = 1._fp / dz(k,0);

        //////////////////////////////////////////////////////////////////
        // d(rho*u)/dx + d(rho*w)/dz = 0  (inserted at p' equation index)
        //////////////////////////////////////////////////////////////////
        // Add mass flux left x-direction (multiply by -1./dx)
        ru_im1 += -rdx * 0.5_fp;
        ru_i   += -rdx * 0.5_fp;
        // Add mass flux right x-direction (multiply by 1./dx)
        ru_i   +=  rdx * 0.5_fp;
        ru_ip1 +=  rdx * 0.5_fp;
        // Add mass flux bottom z-direction (multiply by -1./dz)
        real mult = -rdz;
        if (k == 0   ) mult = 0;  // If k==0, bottom mass flux is zero
        rw_km1 += mult * 0.5_fp;
        rw_k   += mult * 0.5_fp;
        // Add mass flux top z-direction (multiply by 1./dz)
        mult = rdz;
        if (k == nz-1) mult = 0;  // If k==nz-1, top mass flux is zero
        rw_k   += mult * 0.5_fp;
        rw_kp1 += mult * 0.5_fp;

        int eqn_index = k*nx*neq + i*neq + id_p;
        A.insert( eqn_index , k*nx*neq + im1*neq + id_ru ) = ru_im1;
        A.insert( eqn_index , k*nx*neq + i  *neq + id_ru ) = ru_i;
        A.insert( eqn_index , k*nx*neq + ip1*neq + id_ru ) = ru_ip1;
        if (k > 0   ) A.insert( eqn_index , km1*nx*neq + i*neq + id_rw ) = rw_km1;
                      A.insert( eqn_index , k  *nx*neq + i*neq + id_rw ) = rw_k;
        if (k < nz-1) A.insert( eqn_index , kp1*nx*neq + i*neq + id_rw ) = rw_kp1;

        /////////////////////////////////////////////
        // (rho*u)_new + d(p'_new)/dx = (rho*u)_old
        /////////////////////////////////////////////
        eqn_index = k*nx*neq + i*neq + id_ru;
        A.insert( eqn_index , k*nx*neq + im1*neq + id_p  ) = -rdx*0.5_fp;
        A.insert( eqn_index , k*nx*neq + ip1*neq + id_p  ) =  rdx*0.5_fp;
        A.insert( eqn_index , k*nx*neq + i  *neq + id_ru ) = 1;

        /////////////////////////////////////////////
        // (rho*w)_new + d(p'_new)/dz = (rho*w)_old
        /////////////////////////////////////////////
        p_km1 = -rdz*0.5_fp;
        p_k  = 0;
        p_kp1 =  rdz*0.5_fp;
        if (k == 0   ) p_k += p_km1;
        if (k == nz-1) p_k += p_kp1;

        eqn_index = k*nx*neq + i*neq + id_rw;
        if (k > 0   ) A.insert( eqn_index , km1*nx*neq + i*neq + id_p  ) = p_km1;
                      A.insert( eqn_index , k  *nx*neq + i*neq + id_p  ) = p_k;
        if (k < nz-1) A.insert( eqn_index , kp1*nx*neq + i*neq + id_p  ) = p_kp1;
        A.insert( eqn_index , k  *nx*neq + i*neq + id_rw ) = 1;
      }
    }

    A.makeCompressed();
    momdiv_solver.analyzePattern(A);
    momdiv_solver.factorize(A);
    if (momdiv_solver.info() != 0) {
      std::cout << momdiv_solver.lastErrorMessage() << std::endl;
      endrun("ERROR: LU decomp was unsuccessful");
    }
  }



  void build_matrix_solver_3d() {
    int constexpr id_p  = 0;
    int constexpr id_ru = 1;
    int constexpr id_rv = 2;
    int constexpr id_rw = 3;
    int constexpr neq = 4;

    Eigen::SparseMatrix<real> A(neq*nx*ny*nz,neq*nx*ny*nz);
    A.reserve(Eigen::VectorXi::Constant(neq*nx*ny*nz,10));

    std::cout << "Assembling matrix\n";

    for (int k=0; k < nz; k++) {
      for (int j=0; j < ny; j++) {
        for (int i=0; i < nx; i++) {
          int im1 = i-1;  if (im1 < 0   ) im1 += nx;
          int ip1 = i+1;  if (ip1 > nx-1) ip1 -= nx;
          int jm1 = j-1;  if (jm1 < 0   ) jm1 += ny;
          int jp1 = j+1;  if (jp1 > ny-1) jp1 -= ny;
          int km1 = k-1;  if (km1 < 0   ) km1  = 0;
          int kp1 = k+1;  if (kp1 > nz-1) kp1  = nz-1;
          real p_km1=0, p_k=0, p_kp1=0, ru_im1=0, ru_i=0, ru_ip1=0, rv_jm1=0, rv_j=0, rv_jp1=0, rw_km1=0, rw_k=0, rw_kp1=0;
          real rdx = 1._fp / dx;
          real rdy = 1._fp / dy;
          real rdz = 1._fp / dz(k,0);

          //////////////////////////////////////////////////////////////////
          // d(rho*u)/dx + d(rho*w)/dz = 0  (inserted at p' equation index)
          //////////////////////////////////////////////////////////////////
          // Add mass flux left x-direction (multiply by -1./dx)
          ru_im1 += -rdx * 0.5_fp;
          ru_i   += -rdx * 0.5_fp;
          // Add mass flux right x-direction (multiply by 1./dx)
          ru_i   +=  rdx * 0.5_fp;
          ru_ip1 +=  rdx * 0.5_fp;
          // Add mass flux left x-direction (multiply by -1./dx)
          rv_jm1 += -rdy * 0.5_fp;
          rv_j   += -rdy * 0.5_fp;
          // Add mass flux right x-direction (multiply by 1./dx)
          rv_j   +=  rdy * 0.5_fp;
          rv_jp1 +=  rdy * 0.5_fp;
          // Add mass flux bottom z-direction (multiply by -1./dz)
          real mult = -rdz;
          if (k == 0   ) mult = 0;  // If k==0, bottom mass flux is zero
          rw_km1 += mult * 0.5_fp;
          rw_k   += mult * 0.5_fp;
          // Add mass flux top z-direction (multiply by 1./dz)
          mult = rdz;
          if (k == nz-1) mult = 0;  // If k==nz-1, top mass flux is zero
          rw_k   += mult * 0.5_fp;
          rw_kp1 += mult * 0.5_fp;

          int eqn_index = k*nx*ny*neq + j*nx*neq + i*neq + id_p;
                        A.insert( eqn_index , k  *nx*ny*neq + j  *nx*neq + im1*neq + id_ru ) = ru_im1;
                        A.insert( eqn_index , k  *nx*ny*neq + j  *nx*neq + i  *neq + id_ru ) = ru_i;
                        A.insert( eqn_index , k  *nx*ny*neq + j  *nx*neq + ip1*neq + id_ru ) = ru_ip1;
                        A.insert( eqn_index , k  *nx*ny*neq + jm1*nx*neq + i  *neq + id_rv ) = rv_jm1;
                        A.insert( eqn_index , k  *nx*ny*neq + j  *nx*neq + i  *neq + id_rv ) = rv_j;
                        A.insert( eqn_index , k  *nx*ny*neq + jp1*nx*neq + i  *neq + id_rv ) = rv_jp1;
          if (k > 0   ) A.insert( eqn_index , km1*nx*ny*neq + j  *nx*neq + i  *neq + id_rw ) = rw_km1;
                        A.insert( eqn_index , k  *nx*ny*neq + j  *nx*neq + i  *neq + id_rw ) = rw_k;
          if (k < nz-1) A.insert( eqn_index , kp1*nx*ny*neq + j  *nx*neq + i  *neq + id_rw ) = rw_kp1;

          /////////////////////////////////////////////
          // (rho*u)_new + d(p'_new)/dx = (rho*u)_old
          /////////////////////////////////////////////
          eqn_index = k*nx*ny*neq + j*nx*neq + i*neq + id_ru;
          A.insert( eqn_index , k*nx*ny*neq + j*nx*neq + im1*neq + id_p  ) = -rdx*0.5_fp;
          A.insert( eqn_index , k*nx*ny*neq + j*nx*neq + ip1*neq + id_p  ) =  rdx*0.5_fp;
          A.insert( eqn_index , k*nx*ny*neq + j*nx*neq + i  *neq + id_ru ) = 1;

          /////////////////////////////////////////////
          // (rho*v)_new + d(p'_new)/dy = (rho*v)_old
          /////////////////////////////////////////////
          eqn_index = k*nx*ny*neq + j*nx*neq + i*neq + id_rv;
          A.insert( eqn_index , k*nx*ny*neq + jm1*nx*neq + i*neq + id_p  ) = -rdy*0.5_fp;
          A.insert( eqn_index , k*nx*ny*neq + jp1*nx*neq + i*neq + id_p  ) =  rdy*0.5_fp;
          A.insert( eqn_index , k*nx*ny*neq + j  *nx*neq + i*neq + id_rv ) = 1;

          /////////////////////////////////////////////
          // (rho*w)_new + d(p'_new)/dz = (rho*w)_old
          /////////////////////////////////////////////
          p_km1 = -rdz*0.5_fp;
          p_k  = 0;
          p_kp1 =  rdz*0.5_fp;
          if (k == 0   ) p_k += p_km1;
          if (k == nz-1) p_k += p_kp1;

          eqn_index = k*nx*ny*neq + j*nx*neq + i*neq + id_rw;
          if (k > 0   ) A.insert( eqn_index , km1*ny*nx*neq + j*nx*neq + i*neq + id_p  ) = p_km1;
                        A.insert( eqn_index , k  *ny*nx*neq + j*nx*neq + i*neq + id_p  ) = p_k;
          if (k < nz-1) A.insert( eqn_index , kp1*ny*nx*neq + j*nx*neq + i*neq + id_p  ) = p_kp1;
                        A.insert( eqn_index , k  *ny*nx*neq + j*nx*neq + i*neq + id_rw ) = 1;
        }
      }
    }

    std::cout << "Compressing matrix\n";
    A.makeCompressed();

    std::cout << "Number of non-zeros: " << A.nonZeros() << "\n";

    std::cout << "Analyzing matrix\n";
    momdiv_solver.analyzePattern(A);

    std::cout << "Factorizing matrix\n";
    momdiv_solver.factorize(A);

    if (momdiv_solver.info() != 0) {
      std::cout << momdiv_solver.lastErrorMessage() << std::endl;
      endrun("ERROR: LU decomp was unsuccessful");
    }
  }



  void remove_momentum_divergence_2d(real5d &state, real4d &mass_flux_x, real4d &mass_flux_z) {
    int constexpr neq = 3;
    int constexpr id_p  = 0;
    int constexpr id_ru = 1;
    int constexpr id_rw = 2;

    // Set up RHS vector
    Eigen::VectorXd q(neq*nx*nz);
    Eigen::VectorXd b(neq*nx*nz);
    for (int k=0; k < nz; k++) {
      for (int j=0; j < ny; j++) {
        for (int i=0; i < nx; i++) {
          for (int iens=0; iens < nens; iens++) {
            b(k*nx*neq + i*neq + id_p ) = 0.;
            b(k*nx*neq + i*neq + id_ru) = state(idU,hs+k,hs+j,hs+i,iens);
            b(k*nx*neq + i*neq + id_rw) = state(idW,hs+k,hs+j,hs+i,iens);
          }
        }
      }
    }

    // Solve 
    q = momdiv_solver.solve(b);

    // Copy solution into appropriate arrays
    for (int k=0; k < nz; k++) {
      for (int j=0; j < ny; j++) {
        for (int i=0; i < nx; i++) {
          for (int iens=0; iens < nens; iens++) {
            state(idU,hs+k,hs+j,hs+i,iens) = q(k*nx*neq + i*neq + id_ru);
            state(idW,hs+k,hs+j,hs+i,iens) = q(k*nx*neq + i*neq + id_rw);
          }
        }
      }
    }

    // Compute mass fluxes
    for (int k=0; k < nz; k++) {
      for (int j=0; j < ny; j++) {
        for (int i=0; i < nx; i++) {
          for (int iens=0; iens < nens; iens++) {
            int im1 = i-1;  if (im1 < 0   ) im1 += nx;
            int ip1 = i+1;  if (ip1 > nx-1) ip1 -= nx;
            int km1 = k-1;  if (km1 < 0   ) km1  = 0;
            int kp1 = k+1;  if (kp1 > nz-1) kp1  = nz-1;

            real ru_im1 = state(idU,hs+k,hs+j,hs+im1,iens);
            real ru_i   = state(idU,hs+k,hs+j,hs+i  ,iens);
            real ru_ip1 = state(idU,hs+k,hs+j,hs+ip1,iens);
            real rw_km1 = state(idW,hs+km1,hs+j,hs+i,iens);  if (k <= 0   ) rw_km1 = 0;
            real rw_k   = state(idW,hs+k  ,hs+j,hs+i,iens);
            real rw_kp1 = state(idW,hs+kp1,hs+j,hs+i,iens);  if (k >= nz-1) rw_kp1 = 0;

            real ru_L = (ru_im1 + ru_i  )/2;
            real ru_R = (ru_i   + ru_ip1)/2;
            real rw_L = (rw_km1 + rw_k  )/2;
            real rw_R = (rw_k   + rw_kp1)/2;

                           mass_flux_x(k,j,i  ,iens) = ru_L;
            if (i == nx-1) mass_flux_x(k,j,i+1,iens) = ru_R;
                           mass_flux_z(k  ,j,i,iens) = rw_L;
            if (k == nz-1) mass_flux_z(k+1,j,i,iens) = 0;
            if (k == 0   ) mass_flux_z(k  ,j,i,iens) = 0;
          }
        }
      }
    }
  }



  void remove_momentum_divergence_3d(real5d &state, real4d &mass_flux_x, real4d &mass_flux_y, real4d &mass_flux_z) {
    int constexpr neq = 4;
    int constexpr id_p  = 0;
    int constexpr id_ru = 1;
    int constexpr id_rv = 2;
    int constexpr id_rw = 3;

    // Set up RHS vector
    Eigen::VectorXd q(neq*nx*ny*nz);
    Eigen::VectorXd b(neq*nx*ny*nz);
    for (int k=0; k < nz; k++) {
      for (int j=0; j < ny; j++) {
        for (int i=0; i < nx; i++) {
          for (int iens=0; iens < nens; iens++) {
            b(k*nx*ny*neq + j*nx*neq + i*neq + id_p ) = 0.;
            b(k*nx*ny*neq + j*nx*neq + i*neq + id_ru) = state(idU,hs+k,hs+j,hs+i,iens);
            b(k*nx*ny*neq + j*nx*neq + i*neq + id_rv) = state(idV,hs+k,hs+j,hs+i,iens);
            b(k*nx*ny*neq + j*nx*neq + i*neq + id_rw) = state(idW,hs+k,hs+j,hs+i,iens);
          }
        }
      }
    }

    // Solve 
    q = momdiv_solver.solve(b);

    // Copy solution into appropriate arrays
    for (int k=0; k < nz; k++) {
      for (int j=0; j < ny; j++) {
        for (int i=0; i < nx; i++) {
          for (int iens=0; iens < nens; iens++) {
            state(idU,hs+k,hs+j,hs+i,iens) = q(k*nx*ny*neq + j*nx*neq + i*neq + id_ru);
            state(idV,hs+k,hs+j,hs+i,iens) = q(k*nx*ny*neq + j*nx*neq + i*neq + id_rv);
            state(idW,hs+k,hs+j,hs+i,iens) = q(k*nx*ny*neq + j*nx*neq + i*neq + id_rw);
          }
        }
      }
    }

    // Compute mass fluxes
    for (int k=0; k < nz; k++) {
      for (int j=0; j < ny; j++) {
        for (int i=0; i < nx; i++){
          for (int iens=0; iens < nens; iens++) {
            int im1 = i-1;  if (im1 < 0   ) im1 += nx;
            int ip1 = i+1;  if (ip1 > nx-1) ip1 -= nx;
            int jm1 = j-1;  if (jm1 < 0   ) jm1 += ny;
            int jp1 = j+1;  if (jp1 > ny-1) jp1 -= ny;
            int km1 = k-1;  if (km1 < 0   ) km1  = 0;
            int kp1 = k+1;  if (kp1 > nz-1) kp1  = nz-1;

            real ru_im1 = state(idU,hs+k,hs+j,hs+im1,iens);
            real ru_i   = state(idU,hs+k,hs+j,hs+i  ,iens);
            real ru_ip1 = state(idU,hs+k,hs+j,hs+ip1,iens);
            real rv_jm1 = state(idV,hs+k,hs+jm1,hs+i,iens);
            real rv_j   = state(idV,hs+k,hs+j  ,hs+i,iens);
            real rv_jp1 = state(idV,hs+k,hs+jp1,hs+i,iens);
            real rw_km1 = state(idW,hs+km1,hs+j,hs+i,iens);  if (k <= 0   ) rw_km1 = 0;
            real rw_k   = state(idW,hs+k  ,hs+j,hs+i,iens);
            real rw_kp1 = state(idW,hs+kp1,hs+j,hs+i,iens);  if (k >= nz-1) rw_kp1 = 0;

            real ru_L = (ru_im1 + ru_i  )/2;
            real ru_R = (ru_i   + ru_ip1)/2;
            real rv_L = (rv_jm1 + rv_j  )/2;
            real rv_R = (rv_j   + rv_jp1)/2;
            real rw_L = (rw_km1 + rw_k  )/2;
            real rw_R = (rw_k   + rw_kp1)/2;

                           mass_flux_x(k,j,i  ,iens) = ru_L;
            if (i == nx-1) mass_flux_x(k,j,i+1,iens) = ru_R;
                           mass_flux_y(k,j  ,i,iens) = rv_L;
            if (j == ny-1) mass_flux_y(k,j+1,i,iens) = rv_R;
                           mass_flux_z(k  ,j,i,iens) = rw_L;
            if (k == nz-1) mass_flux_z(k+1,j,i,iens) = 0;
            if (k == 0   ) mass_flux_z(k  ,j,i,iens) = 0;
          }
        }
      }
    }
  }


  
  real compute_divergence(real5d const &state) {
    real4d absdiv("absdiv",nz,ny,nx,nens);
    parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      int im1 = i-1;  if (im1 < 0   ) im1 += nx;
      int ip1 = i+1;  if (ip1 > nx-1) ip1 -= nx;
      int jm1 = j-1;  if (jm1 < 0   ) jm1 += ny;
      int jp1 = j+1;  if (jp1 > ny-1) jp1 -= ny;
      int km1 = k-1;  if (km1 < 0   ) km1  = 0;
      int kp1 = k+1;  if (kp1 > nz-1) kp1  = nz-1;
      real ru_im1 = state(idU,hs+k,hs+j,hs+im1,iens);
      real ru_i   = state(idU,hs+k,hs+j,hs+i  ,iens);
      real ru_ip1 = state(idU,hs+k,hs+j,hs+ip1,iens);
      real rv_jm1 = state(idV,hs+k,hs+jm1,hs+i,iens);  if (sim2d) rv_jm1 = 0;
      real rv_j   = state(idV,hs+k,hs+j  ,hs+i,iens);  if (sim2d) rv_j   = 0;
      real rv_jp1 = state(idV,hs+k,hs+jp1,hs+i,iens);  if (sim2d) rv_jp1 = 0;
      real rw_km1 = state(idW,hs+km1,hs+j,hs+i,iens);  if (k == 0   ) rw_km1 = 0;
      real rw_k   = state(idW,hs+k  ,hs+j,hs+i,iens);
      real rw_kp1 = state(idW,hs+kp1,hs+j,hs+i,iens);  if (k == nz-1) rw_kp1 = 0;
      real mfx_L = 0.5_fp * (ru_im1 + ru_i  );
      real mfx_R = 0.5_fp * (ru_i   + ru_ip1);
      real mfy_L = 0.5_fp * (rv_jm1 + rv_j  );  if (sim2d) mfy_L = 0;
      real mfy_R = 0.5_fp * (rv_j   + rv_jp1);  if (sim2d) mfy_R = 0;
      real mfz_L = 0.5_fp * (rw_km1 + rw_k  );  if (k == 0   ) mfz_L = 0;
      real mfz_R = 0.5_fp * (rw_k   + rw_kp1);  if (k == nz-1) mfz_R = 0;
      real div_u = ( mfx_R - mfx_L ) / (2*dx        );
      real div_v = ( mfy_R - mfy_L ) / (2*dy        );
      real div_w = ( mfz_R - mfz_L ) / (2*dz(k,iens));
      absdiv(k,j,i,iens) = abs( div_u + div_v + div_w );
    });
    return yakl::intrinsics::maxval(absdiv);
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
    YAKL_SCOPE( s2c                     , this->sten_to_coefs          );
    YAKL_SCOPE( weno_recon_lower        , this->weno_recon_lower       );
    YAKL_SCOPE( idl                     , this->idl                    );
    YAKL_SCOPE( sigma                   , this->sigma                  );
    YAKL_SCOPE( hyDensCells             , this->hyDensCells            );
    YAKL_SCOPE( hyDensThetaCells        , this->hyDensThetaCells       );
    YAKL_SCOPE( hyDensGLL               , this->hyDensGLL              );
    YAKL_SCOPE( hyDensThetaGLL          , this->hyDensThetaGLL         );
    YAKL_SCOPE( hyPressureGLL           , this->hyPressureGLL          );
    YAKL_SCOPE( sim2d                   , this->sim2d                  );
    YAKL_SCOPE( derivMatrix             , this->derivMatrix            );
    YAKL_SCOPE( dx                      , this->dx                     );
    YAKL_SCOPE( dy                      , this->dy                     );
    YAKL_SCOPE( dz                      , this->dz                     );
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

    auto state_init   = state  .createDeviceCopy();
    auto tracers_init = tracers.createDeviceCopy();

    std::cout << "\n*****Beginning abs div: " << compute_divergence(state) << "\n";

    parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      real theta = ( state(idT,hs+k,hs+j,hs+i,iens) + hyDensThetaCells(k,iens) ) / hyDensCells(k,iens);
      real dens = hyDensThetaCells(k,iens) / theta;
      state(idW,hs+k,hs+j,hs+i,iens) += -dt * ( dens - hyDensCells(k,iens) )*GRAV;
    });

    std::cout << "*****Post buoyancy abs div: " << compute_divergence(state) << "\n";

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
    if (!sim2d) {
      if        (bc_y == BC_PERIODIC) {
        parallel_for( "Spatial.h Y BCs periodic" , SimpleBounds<4>(nz,nx,hs,nens) , YAKL_LAMBDA(int k, int i, int jj, int iens) {
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
        parallel_for( "Spatial.h Y BCs wall" , SimpleBounds<4>(nz,nx,hs,nens) , YAKL_LAMBDA(int k, int i, int jj, int iens) {
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



    real4d mass_flux_x("mass_flux_x",nz  ,ny  ,nx+1,nens);
    real4d mass_flux_y;  if (!sim2d) mass_flux_y = real4d("mass_flux_y",nz  ,ny+1,nx  ,nens);
    real4d mass_flux_z("mass_flux_z",nz+1,ny  ,nx  ,nens);

    if (sim2d) { remove_momentum_divergence_2d( state , mass_flux_x               , mass_flux_z ); }
    else       { remove_momentum_divergence_3d( state , mass_flux_x , mass_flux_y , mass_flux_z ); }

    std::cout << "*****Post div removal 1 abs div: " << compute_divergence(state) << "\n";

    real6d stateLimits_x ("stateLimits_x" ,num_state  ,2,nz  ,ny  ,nx+1,nens);
    real6d stateLimits_y;   if (!sim2d) stateLimits_y  = real6d("stateLimits_y" ,num_state  ,2,nz  ,ny+1,nx  ,nens);
    real6d stateLimits_z ("stateLimits_z" ,num_state  ,2,nz+1,ny  ,nx  ,nens);
    real6d tracerLimits_x("tracerLimits_x",num_tracers,2,nz  ,ny  ,nx+1,nens);
    real6d tracerLimits_y;  if (!sim2d) tracerLimits_y = real6d("tracerLimits_y",num_tracers,2,nz  ,ny+1,nx  ,nens);
    real6d tracerLimits_z("tracerLimits_z",num_tracers,2,nz+1,ny  ,nx  ,nens);

    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      SArray<real,1,ord> stencil;
      SArray<real,1,ngll> gll;

      ///////////////////////////////////////////////////////
      // X-direction
      ///////////////////////////////////////////////////////
      // Density
      for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idR,hs+k,hs+j,i+ii,iens); }
      reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
      for (int ii=0; ii < ngll; ii++) { gll(ii) += hyDensCells(k,iens); } // Add hydrostasis back on
      stateLimits_x(idR,1,k,j,i  ,iens) = gll(0     );
      stateLimits_x(idR,0,k,j,i+1,iens) = gll(ngll-1);

      // u values and derivatives
      for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idU,hs+k,hs+j,i+ii,iens); }
      reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
      stateLimits_x(idU,1,k,j,i  ,iens) = gll(0     );
      stateLimits_x(idU,0,k,j,i+1,iens) = gll(ngll-1);

      // v
      for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idV,hs+k,hs+j,i+ii,iens); }
      reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
      stateLimits_x(idV,1,k,j,i  ,iens) = gll(0     );
      stateLimits_x(idV,0,k,j,i+1,iens) = gll(ngll-1);

      // w
      for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idW,hs+k,hs+j,i+ii,iens); }
      reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
      stateLimits_x(idW,1,k,j,i  ,iens) = gll(0     );
      stateLimits_x(idW,0,k,j,i+1,iens) = gll(ngll-1);

      // theta
      for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idT,hs+k,hs+j,i+ii,iens); }
      reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
      for (int ii=0; ii < ngll; ii++) { gll(ii) += hyDensThetaCells(k,iens); } // Add hydrostasis back on
      stateLimits_x(idT,1,k,j,i  ,iens) = gll(0     );
      stateLimits_x(idT,0,k,j,i+1,iens) = gll(ngll-1);

      // Only process one tracer at a time to save on local memory / register requirements
      for (int tr=0; tr < num_tracers; tr++) {
        for (int ii=0; ii < ord; ii++) { stencil(ii) = tracers(tr,hs+k,hs+j,i+ii,iens); }
        reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
        if (tracer_pos(tr)) {
          for (int ii=0; ii < ngll; ii++) { gll(ii) = max( 0._fp , gll(ii) ); }
        }
        tracerLimits_x(tr,1,k,j,i  ,iens) = gll(0     );
        tracerLimits_x(tr,0,k,j,i+1,iens) = gll(ngll-1);
      }

      if (!sim2d) {
        ///////////////////////////////////////////////////////
        // Y-direction
        ///////////////////////////////////////////////////////
        // Density
        for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idR,hs+k,j+jj,hs+i,iens); }
        reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
        for (int jj=0; jj < ngll; jj++) { gll(jj) += hyDensCells(k,iens); } // Add hydrostasis back on
        stateLimits_y(idR,1,k,j  ,i,iens) = gll(0     );
        stateLimits_y(idR,0,k,j+1,i,iens) = gll(ngll-1);

        // u values and derivatives
        for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idU,hs+k,j+jj,hs+i,iens); }
        reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
        stateLimits_y(idU,1,k,j  ,i,iens) = gll(0     );
        stateLimits_y(idU,0,k,j+1,i,iens) = gll(ngll-1);

        // v
        for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idV,hs+k,j+jj,hs+i,iens); }
        reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
        stateLimits_y(idV,1,k,j  ,i,iens) = gll(0     );
        stateLimits_y(idV,0,k,j+1,i,iens) = gll(ngll-1);

        // w
        for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idW,hs+k,j+jj,hs+i,iens); }
        reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
        stateLimits_y(idW,1,k,j  ,i,iens) = gll(0     );
        stateLimits_y(idW,0,k,j+1,i,iens) = gll(ngll-1);

        // theta
        for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idT,hs+k,j+jj,hs+i,iens); }
        reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
        for (int jj=0; jj < ngll; jj++) { gll(jj) += hyDensThetaCells(k,iens); } // Add hydrostasis back on
        stateLimits_y(idT,1,k,j  ,i,iens) = gll(0     );
        stateLimits_y(idT,0,k,j+1,i,iens) = gll(ngll-1);

        // Only process one tracer at a time to save on local memory / register requirements
        for (int tr=0; tr < num_tracers; tr++) {
          for (int jj=0; jj < ord; jj++) { stencil(jj) = tracers(tr,hs+k,j+jj,hs+i,iens); }
          reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
          if (tracer_pos(tr)) {
            for (int jj=0; jj < ngll; jj++) { gll(jj) = max( 0._fp , gll(jj) ); }
          }
          tracerLimits_y(tr,1,k,j  ,i,iens) = gll(0     );
          tracerLimits_y(tr,0,k,j+1,i,iens) = gll(ngll-1);
        }
      }

      ///////////////////////////////////////////////////////
      // Z-direction
      ///////////////////////////////////////////////////////
      // Density
      for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idR,k+kk,hs+j,hs+i,iens); }
      reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
      for (int kk=0; kk < ngll; kk++) { gll(kk) += hyDensGLL(k,kk,iens); } // Add hydrostasis back on
      stateLimits_z(idR,1,k  ,j,i,iens) = gll(0     );
      stateLimits_z(idR,0,k+1,j,i,iens) = gll(ngll-1);

      // u values and derivatives
      for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idU,k+kk,hs+j,hs+i,iens); }
      reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
      stateLimits_z(idU,1,k  ,j,i,iens) = gll(0     );
      stateLimits_z(idU,0,k+1,j,i,iens) = gll(ngll-1);

      // v
      for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idV,k+kk,hs+j,hs+i,iens); }
      reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
      stateLimits_z(idV,1,k  ,j,i,iens) = gll(0     );
      stateLimits_z(idV,0,k+1,j,i,iens) = gll(ngll-1);

      // w
      for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idW,k+kk,hs+j,hs+i,iens); }
      reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_winds );
      if (bc_z == BC_WALL) {
        if (k == nz-1) gll(ngll-1) = 0;
        if (k == 0   ) gll(0     ) = 0;
      }
      stateLimits_z(idW,1,k  ,j,i,iens) = gll(0     );
      stateLimits_z(idW,0,k+1,j,i,iens) = gll(ngll-1);

      // theta
      for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idT,k+kk,hs+j,hs+i,iens); }
      reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
      for (int kk=0; kk < ngll; kk++) { gll(kk) += hyDensThetaGLL(k,kk,iens); } // Add hydrostasis back on
      stateLimits_z(idT,1,k  ,j,i,iens) = gll(0     );
      stateLimits_z(idT,0,k+1,j,i,iens) = gll(ngll-1);

      for (int tr=0; tr < num_tracers; tr++) {
        for (int kk=0; kk < ord; kk++) { stencil(kk) = tracers(tr,k+kk,hs+j,hs+i,iens); }
        reconstruct_gll_values( stencil , gll , c2g , s2g , s2c , weno_recon_lower , idl , sigma , weno_scalars );
        if (tracer_pos(tr)) {
          for (int kk=0; kk < ngll; kk++) { gll(kk) = max( 0._fp , gll(kk) ); }
        }
        tracerLimits_z(tr,1,k  ,j,i,iens) = gll(0     );
        tracerLimits_z(tr,0,k+1,j,i,iens) = gll(ngll-1);
      }
    });

    ////////////////////////////////////////////////
    // BCs for the state edge estimates
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
    if (!sim2d) {
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
    }
    parallel_for( SimpleBounds<3>(ny,nx,nens) , YAKL_LAMBDA (int j, int i, int iens) {
      for (int l = 0; l < num_state; l++) {
        if        (bc_z == BC_PERIODIC) {
          stateLimits_z(l,0,0 ,j,i,iens) = stateLimits_z(l,0,nz,j,i,iens);
          stateLimits_z(l,1,nz,j,i,iens) = stateLimits_z(l,1,0 ,j,i,iens);
        } else if (bc_z == BC_WALL    ) {
          stateLimits_z(l,0,0 ,j,i,iens) = stateLimits_z(l,1,0 ,j,i,iens);
          stateLimits_z(l,1,nz,j,i,iens) = stateLimits_z(l,0,nz,j,i,iens);
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
    // Compute the upwind fluxes (X-direction)
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<4>(nz,ny,nx+1,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      real dens = hyDensCells(k,iens);
      real u_L = stateLimits_x(idU,0,k,j,i,iens)/dens;   real u_R = stateLimits_x(idU,1,k,j,i,iens)/dens;
      real v_L = stateLimits_x(idV,0,k,j,i,iens)/dens;   real v_R = stateLimits_x(idV,1,k,j,i,iens)/dens;
      real w_L = stateLimits_x(idW,0,k,j,i,iens)/dens;   real w_R = stateLimits_x(idW,1,k,j,i,iens)/dens;
      real t_L = stateLimits_x(idT,0,k,j,i,iens)/dens;   real t_R = stateLimits_x(idT,1,k,j,i,iens)/dens;

      real u = 0.5_fp * (u_L + u_R);

      real mass_flux = mass_flux_x(k,j,i,iens);

      if (u > 0) {
        stateLimits_x(idR,0,k,j,i,iens) = mass_flux;
        stateLimits_x(idU,0,k,j,i,iens) = mass_flux * u_L;
        stateLimits_x(idV,0,k,j,i,iens) = mass_flux * v_L;
        stateLimits_x(idW,0,k,j,i,iens) = mass_flux * w_L;
        stateLimits_x(idT,0,k,j,i,iens) = mass_flux * t_L;
        for (int tr=0; tr < num_tracers; tr++) {
          tracerLimits_x(tr,0,k,j,i,iens) = mass_flux * tracerLimits_x(tr,0,k,j,i,iens) / dens;
        }
      } else {
        stateLimits_x(idR,0,k,j,i,iens) = mass_flux;
        stateLimits_x(idU,0,k,j,i,iens) = mass_flux * u_R;
        stateLimits_x(idV,0,k,j,i,iens) = mass_flux * v_R;
        stateLimits_x(idW,0,k,j,i,iens) = mass_flux * w_R;
        stateLimits_x(idT,0,k,j,i,iens) = mass_flux * t_R;
        for (int tr=0; tr < num_tracers; tr++) {
          tracerLimits_x(tr,0,k,j,i,iens) = mass_flux * tracerLimits_x(tr,1,k,j,i,iens) / dens;
        }
      }
    });
    if (!sim2d) {
      parallel_for( SimpleBounds<4>(nz,ny+1,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        real dens = hyDensCells(k,iens);
        real u_L = stateLimits_y(idU,0,k,j,i,iens)/dens;   real u_R = stateLimits_y(idU,1,k,j,i,iens)/dens;
        real v_L = stateLimits_y(idV,0,k,j,i,iens)/dens;   real v_R = stateLimits_y(idV,1,k,j,i,iens)/dens;
        real w_L = stateLimits_y(idW,0,k,j,i,iens)/dens;   real w_R = stateLimits_y(idW,1,k,j,i,iens)/dens;
        real t_L = stateLimits_y(idT,0,k,j,i,iens)/dens;   real t_R = stateLimits_y(idT,1,k,j,i,iens)/dens;

        real v = 0.5_fp * (v_L + v_R);

        real mass_flux = mass_flux_y(k,j,i,iens);

        if (v > 0) {
          stateLimits_y(idR,0,k,j,i,iens) = mass_flux;
          stateLimits_y(idU,0,k,j,i,iens) = mass_flux * u_L;
          stateLimits_y(idV,0,k,j,i,iens) = mass_flux * v_L;
          stateLimits_y(idW,0,k,j,i,iens) = mass_flux * w_L;
          stateLimits_y(idT,0,k,j,i,iens) = mass_flux * t_L;
          for (int tr=0; tr < num_tracers; tr++) {
            tracerLimits_y(tr,0,k,j,i,iens) = mass_flux * tracerLimits_y(tr,0,k,j,i,iens) / dens;
          }
        } else {
          stateLimits_y(idR,0,k,j,i,iens) = mass_flux;
          stateLimits_y(idU,0,k,j,i,iens) = mass_flux * u_R;
          stateLimits_y(idV,0,k,j,i,iens) = mass_flux * v_R;
          stateLimits_y(idW,0,k,j,i,iens) = mass_flux * w_R;
          stateLimits_y(idT,0,k,j,i,iens) = mass_flux * t_R;
          for (int tr=0; tr < num_tracers; tr++) {
            tracerLimits_y(tr,0,k,j,i,iens) = mass_flux * tracerLimits_y(tr,1,k,j,i,iens) / dens;
          }
        }
      });
    }
    parallel_for( SimpleBounds<4>(nz+1,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      real dens;
      if (k < nz) { dens = hyDensGLL(k  ,0     ,iens); }
      else        { dens = hyDensGLL(k-1,ngll-1,iens); }
      real u_L = stateLimits_z(idU,0,k,j,i,iens)/dens;   real u_R = stateLimits_z(idU,1,k,j,i,iens)/dens;
      real v_L = stateLimits_z(idV,0,k,j,i,iens)/dens;   real v_R = stateLimits_z(idV,1,k,j,i,iens)/dens;
      real w_L = stateLimits_z(idW,0,k,j,i,iens)/dens;   real w_R = stateLimits_z(idW,1,k,j,i,iens)/dens;
      real t_L = stateLimits_z(idT,0,k,j,i,iens)/dens;   real t_R = stateLimits_z(idT,1,k,j,i,iens)/dens;

      real w = 0.5_fp * (w_L + w_R);

      real mass_flux = mass_flux_z(k,j,i,iens);

      if (w > 0) {
        stateLimits_z(idR,0,k,j,i,iens) = mass_flux;
        stateLimits_z(idU,0,k,j,i,iens) = mass_flux * u_L;
        stateLimits_z(idV,0,k,j,i,iens) = mass_flux * v_L;
        stateLimits_z(idW,0,k,j,i,iens) = mass_flux * w_L;
        stateLimits_z(idT,0,k,j,i,iens) = mass_flux * t_L;
        for (int tr=0; tr < num_tracers; tr++) {
          tracerLimits_z(tr,0,k,j,i,iens) = mass_flux * tracerLimits_z(tr,0,k,j,i,iens) / dens;
        }
      } else {
        stateLimits_z(idR,0,k,j,i,iens) = mass_flux;
        stateLimits_z(idU,0,k,j,i,iens) = mass_flux * u_R;
        stateLimits_z(idV,0,k,j,i,iens) = mass_flux * v_R;
        stateLimits_z(idW,0,k,j,i,iens) = mass_flux * w_R;
        stateLimits_z(idT,0,k,j,i,iens) = mass_flux * t_R;
        for (int tr=0; tr < num_tracers; tr++) {
          tracerLimits_z(tr,0,k,j,i,iens) = mass_flux * tracerLimits_z(tr,1,k,j,i,iens) / dens;
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA(int k, int j, int i, int iens) {
      for (int l = 0; l < num_state; l++) {
        if ( sim2d && l == idV ) {
          stateTend(l,k,j,i,iens) = 0;
        } else {
                      stateTend(l,k,j,i,iens)  = - ( stateLimits_x(l,0,k,j,i+1,iens) - stateLimits_x(l,0,k,j,i,iens) ) / dx;
          if (!sim2d) stateTend(l,k,j,i,iens) += - ( stateLimits_y(l,0,k,j+1,i,iens) - stateLimits_y(l,0,k,j,i,iens) ) / dy;
                      stateTend(l,k,j,i,iens) += - ( stateLimits_z(l,0,k+1,j,i,iens) - stateLimits_z(l,0,k,j,i,iens) ) / dz(k,iens);
                      stateTend(l,k,j,i,iens) += ( state(l,hs+k,hs+j,hs+i,iens) - state_init(l,hs+k,hs+j,hs+i,iens) ) / dt;
        }
        state(l,hs+k,hs+j,hs+i,iens) = state_init(l,hs+k,hs+j,hs+i,iens);
      }
      for (int l = 0; l < num_tracers; l++) {
                    tracerTend(l,k,j,i,iens)  = - ( tracerLimits_x(l,0,k,j,i+1,iens) - tracerLimits_x(l,0,k,j,i,iens) ) / dx;
        if (!sim2d) tracerTend(l,k,j,i,iens) += - ( tracerLimits_y(l,0,k,j+1,i,iens) - tracerLimits_y(l,0,k,j,i,iens) ) / dy;
                    tracerTend(l,k,j,i,iens) += - ( tracerLimits_z(l,0,k+1,j,i,iens) - tracerLimits_z(l,0,k,j,i,iens) ) / dz(k,iens);
                    tracerTend(l,k,j,i,iens) += ( tracers(l,hs+k,hs+j,hs+i,iens) - tracers_init(l,hs+k,hs+j,hs+i,iens) ) / dt;
        tracers(l,hs+k,hs+j,hs+i,iens) = tracers_init(l,hs+k,hs+j,hs+i,iens);
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA(int k, int j, int i, int iens) {
      for (int l = 0; l < num_state; l++) {
        if ( sim2d && l == idV ) {
        } else {
          state(l,hs+k,hs+j,hs+i,iens) += dt * stateTend(l,k,j,i,iens);
        }
      }
    });

    std::cout << "*****Post advection abs div: " << compute_divergence(state) << "\n";

    // if (sim2d) { remove_momentum_divergence_2d( state , mass_flux_x               , mass_flux_z ); }
    // else       { remove_momentum_divergence_3d( state , mass_flux_x , mass_flux_y , mass_flux_z ); }

    // std::cout << "*****Post div removal 2 abs div: " << compute_divergence(state) << "\n";

    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA(int k, int j, int i, int iens) {
      for (int l = 0; l < num_state; l++) {
        stateTend(l,k,j,i,iens) = ( state(l,hs+k,hs+j,hs+i,iens) - state_init(l,hs+k,hs+j,hs+i,iens) ) / dt;
        state(l,hs+k,hs+j,hs+i,iens) = state_init(l,hs+k,hs+j,hs+i,iens);
      }
    });

  } // computeTendencies



  void switch_directions() {
    dimSwitch = ! dimSwitch;
  }



  const char * getName() { return ""; }



  template <class MICRO>
  void output(DataManager &dm, MICRO const &micro, real etime) const {
    YAKL_SCOPE( dx                    , this->dx                   );
    YAKL_SCOPE( dy                    , this->dy                   );
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
        parallel_for( "Spatial.h output 1" , nx , YAKL_LAMBDA (int i) { xloc(i) = (i+0.5)*dx; });
        nc.write(xloc.createHostCopy(),"x",{"x"});

        // y-coordinate
        real1d yloc("yloc",ny);
        parallel_for( "Spatial.h output 2" , ny , YAKL_LAMBDA (int i) { yloc(i) = (i+0.5)*dy; });
        nc.write(yloc.createHostCopy(),"y",{"y"});

        // z-coordinate
        auto zint = dm.get<real,2>("vertical_interface_height");
        real1d zmid("zmid",nz);
        parallel_for( "Spatial.h output 3" , nz , YAKL_LAMBDA (int i) { zmid(i) = ( zint(i,iens) + zint(i+1,iens) ) / 2; });
        nc.write(zmid.createHostCopy(),"z",{"z"});

        // hydrostatic density, theta, and pressure
        parallel_for( "Spatial.h output 4" , nz , YAKL_LAMBDA (int k) { zmid(k) = hyDensCells(k,iens); });
        nc.write(zmid.createHostCopy(),"hyDens"    ,{"z"});

        parallel_for( "Spatial.h output 5" , nz , YAKL_LAMBDA (int k) { zmid(k) = hyPressureCells(k,iens); });
        nc.write(zmid.createHostCopy(),"hyPressure",{"z"});

        parallel_for( "Spatial.h output 6" , nz , YAKL_LAMBDA (int k) { zmid(k) = hyThetaCells(k,iens); });
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
      parallel_for( "Spatial.h output 7" , SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        data(k,j,i) = state(idR,hs+k,hs+j,hs+i,iens);
      });
      nc.write1(data.createHostCopy(),"dens_pert",{"z","y","x"},ulIndex,"t");
      // u
      parallel_for( "Spatial.h output 8" , SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        data(k,j,i) = state(idU,hs+k,hs+j,hs+i,iens) / ( hyDensCells(k,iens) );
      });
      nc.write1(data.createHostCopy(),"u",{"z","y","x"},ulIndex,"t");
      // v
      parallel_for( "Spatial.h output 9" , SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        data(k,j,i) = state(idV,hs+k,hs+j,hs+i,iens) / ( hyDensCells(k,iens) );
      });
      nc.write1(data.createHostCopy(),"v",{"z","y","x"},ulIndex,"t");
      // w
      parallel_for( "Spatial.h output 10" , SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        data(k,j,i) = state(idW,hs+k,hs+j,hs+i,iens) / ( hyDensCells(k,iens) );
      });
      nc.write1(data.createHostCopy(),"w",{"z","y","x"},ulIndex,"t");
      // theta'
      parallel_for( "Spatial.h output 11" , SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        real r =   hyDensCells     (k,iens);
        real t = ( state(idT,hs+k,hs+j,hs+i,iens) + hyDensThetaCells(k,iens) ) / r;
        data(k,j,i) = t - hyThetaCells(k,iens);
      });
      nc.write1(data.createHostCopy(),"pot_temp_pert",{"z","y","x"},ulIndex,"t");
      // pressure'
      parallel_for( "Spatial.h output 12" , SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        real r  = hyDensCells(k,iens);
        real rt = state(idT,hs+k,hs+j,hs+i,iens) + hyDensThetaCells(k,iens);
        real p  = C0*pow(rt,gamma);
        data(k,j,i) = p - hyPressureCells(k,iens);
      });
      nc.write1(data.createHostCopy(),"pressure_pert",{"z","y","x"},ulIndex,"t");

      for (int tr=0; tr < num_tracers; tr++) {
        parallel_for( "Spatial.h output 13" , SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
          real r = hyDensCells(k,iens);
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



  // ord stencil values to ngll GLL values; store in DTs
  YAKL_INLINE static void reconstruct_gll_values( SArray<real,1,ord> const stencil                      ,
                                                  SArray<real,1,ngll> &gll                              ,
                                                  SArray<real,2,ord,ngll> const &coefs_to_gll           ,
                                                  SArray<real,2,ord,ngll> const &sten_to_gll            ,
                                                  SArray<real,2,ord,ord>  const &sten_to_coefs          ,
                                                  SArray<real,3,hs+1,hs+1,hs+1> const &weno_recon_lower ,
                                                  SArray<real,1,hs+2> const &idl                        ,
                                                  real sigma, bool doweno ) {
    if (doweno) {

      // Reconstruct values
      SArray<real,1,ord> wenoCoefs;
      weno::compute_weno_coefs<ord>( weno_recon_lower , sten_to_coefs , stencil , wenoCoefs , idl , sigma );
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



  YAKL_INLINE static void diffTransformEulerConsX( SArray<real,2,nAder,ngll> &r  ,
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



  YAKL_INLINE static void diffTransformEulerConsY( SArray<real,2,nAder,ngll> &r  ,
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



  YAKL_INLINE static void diffTransformEulerConsZ( SArray<real,2,nAder,ngll> &r  ,
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
                                                   real3d const &hyPressureGLL ,
                                                   real C0, real gamma ,
                                                   int k , real dz , int bc_z , int nz , int iens ) {
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
          if (kt == 0) { drww_p_dz += deriv(s,ii) * ( rww(kt,s) + C0*rt_gamma(kt,s) - hyPressureGLL(k,s,iens) ); }
          else         { drww_p_dz += deriv(s,ii) * ( rww(kt,s) + C0*rt_gamma(kt,s)/2                         ); }
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



  YAKL_INLINE static void diffTransformTracer( SArray<real,2,nAder,ngll> const &r  ,
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



  YAKL_INLINE static void compute_timeAvg( SArray<real,3,num_state,nAder,ngll> &dts , real dt ) {
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



  YAKL_INLINE static void compute_timeAvg( SArray<real,2,nAder,ngll> &dts , real dt ) {
    real dtmult = dt;
    for (int kt=1; kt<nAder; kt++) {
      for (int ii=0; ii<ngll; ii++) {
        dts(0,ii) += dts(kt,ii) * dtmult / (kt+1._fp);
      }
      dtmult *= dt;
    }
  }



  YAKL_INLINE static void compute_timeAvg( SArray<real,2,nAder,ngll> const &dts , SArray<real,1,ngll> &tavg , real dt ) {
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
