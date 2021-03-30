
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
  int static constexpr num_state = 4;
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

  real4d pressure;

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
  real dzmin;

  // For indexing into the state and state tendency arrays
  int static constexpr idU = 0;  // u
  int static constexpr idV = 1;  // v
  int static constexpr idW = 2;  // w
  int static constexpr idT = 3;  // potential temperature perturbation

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


  // When this class is created, initialize num_tracers to zero
  Spatial_operator() {
    num_tracers = 0;
  }



  // Make sure it's odd-order-accurate
  static_assert(ord%2 == 1,"ERROR: ord must be an odd integer");



  template <class MICRO>
  void convert_dynamics_to_coupler_state( DataManager &dm , MICRO &micro ) {
    // real5d state           = dm.get<real,5>( "dynamics_state"   );
    // real5d tracers         = dm.get<real,5>( "dynamics_tracers" );
    // real4d dm_dens_dry     = dm.get<real,4>( "density_dry"      );
    // real4d dm_uvel         = dm.get<real,4>( "uvel"             );
    // real4d dm_vvel         = dm.get<real,4>( "vvel"             );
    // real4d dm_wvel         = dm.get<real,4>( "wvel"             );
    // real4d dm_temp         = dm.get<real,4>( "temp"             );
    // real4d dm_pressure_dry = dm.get<real,4>( "pressure_dry"     );

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

    // int idWV = micro.tracer_index_vapor;

    // MultipleTracers<max_tracers> dm_tracers;
    // for (int tr = 0; tr < num_tracers; tr++) {
    //   auto trac = dm.get<real,4>( tracer_name[tr] );
    //   dm_tracers.add_tracer( trac );
    // }

    // parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
    //   real dens  = hyDensCells(k,iens);
    //   real uvel  = state(idU,hs+k,hs+j,hs+i,iens) / dens;
    //   real vvel  = state(idV,hs+k,hs+j,hs+i,iens) / dens;
    //   real wvel  = state(idW,hs+k,hs+j,hs+i,iens) / dens;
    //   real theta = ( state(idT,hs+k,hs+j,hs+i,iens) + hyDensThetaCells(k,iens) ) / dens;
    //   real pressure = hyPressureCells(k,iens);
    //   real dens_vap = tracers(idWV,hs+k,hs+j,hs+i,iens);
    //   real dens_dry = dens;
    //   for (int tr=0; tr < num_tracers; tr++) {
    //     if (tracer_adds_mass(tr)) dens_dry -= tracers(tr,hs+k,hs+j,hs+i,iens);
    //   }
    //   real temp = pressure / ( dens_dry * Rd + dens_vap * Rv );
    //   real pressure_dry = dens_dry * Rd * temp;
    //   dm_dens_dry    (k,j,i,iens) = dens_dry;
    //   dm_uvel        (k,j,i,iens) = uvel;
    //   dm_vvel        (k,j,i,iens) = vvel;
    //   dm_wvel        (k,j,i,iens) = wvel;
    //   dm_temp        (k,j,i,iens) = temp;
    //   dm_pressure_dry(k,j,i,iens) = pressure_dry;
    //   for (int tr=0; tr < num_tracers; tr++) {
    //     dm_tracers(tr,k,j,i,iens) = tracers(tr,hs+k,hs+j,hs+i,iens);
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
    // real4d dm_pressure_dry = dm.get<real,4>( "pressure_dry"     );

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

    // int idWV = micro.tracer_index_vapor;

    // MultipleTracers<max_tracers> dm_tracers;
    // for (int tr = 0; tr < num_tracers; tr++) {
    //   auto trac = dm.get<real,4>( tracer_name[tr] );
    //   dm_tracers.add_tracer( trac );
    // }

    // parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
    //   for (int tr=0; tr < num_tracers; tr++) {
    //     tracers(tr,hs+k,hs+j,hs+i,iens) = dm_tracers(tr,k,j,i,iens);
    //   }
    //   real dens_dry     = dm_dens_dry    (k,j,i,iens);
    //   real uvel         = dm_uvel        (k,j,i,iens);
    //   real vvel         = dm_vvel        (k,j,i,iens);
    //   real wvel         = dm_wvel        (k,j,i,iens);
    //   real temp         = dm_temp        (k,j,i,iens);
    //   real pressure_dry = dm_pressure_dry(k,j,i,iens);
    //   real dens_vap     = tracers(idWV,hs+k,hs+j,hs+i,iens);
    //   real dens         = hyDensCells(k,iens);
    //   real pressure     = hyPressureCells(k,iens);
    //   real theta        = pow( pressure / C0 , 1._fp / gamma ) / dens;
    //   state(idU,hs+k,hs+j,hs+i,iens) = dens * uvel;
    //   state(idV,hs+k,hs+j,hs+i,iens) = dens * vvel;
    //   state(idW,hs+k,hs+j,hs+i,iens) = dens * wvel;
    //   state(idT,hs+k,hs+j,hs+i,iens) = dens * theta - hyDensThetaCells(k,iens);
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

    int idWV = micro.tracer_index_vapor;
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
            real svp   = micro.saturation_vapor_pressure(temp);     // Self-explanatory
            real p_v   = pert*svp;                                  // Multiply profile by saturation vapor pressure
            real r_v   = p_v / (Rv*temp);                           // Compute vapor density

            real wt = gllWts_ord(kk) * gllWts_ord(jj) * gllWts_ord(ii);
            tracers(idWV,hs+k,hs+j,hs+i,iens) += r_v / (rh+r_v) * wt;
          }
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
    return 4;
  }



  // Given the model data and CFL value, compute the maximum stable time step
  template <class MICRO>
  real compute_time_step(real cfl, DataManager &dm, MICRO const &micro) {

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
    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      // Get the state
      real r = hyDensCells(k,iens);
      real u = state(idU,hs+k,hs+j,hs+i,iens) / r;
      real v = state(idV,hs+k,hs+j,hs+i,iens) / r;
      real w = state(idW,hs+k,hs+j,hs+i,iens) / r;

      real udt = cfl * dx         / ( abs(u) + 1.e-6 );
      real vdt = cfl * dy         / ( abs(v) + 1.e-6 );
      real wdt = cfl * dz(k,iens) / ( abs(w) + 1.e-6 );

      // Compute the min of the max stable time steps
      dt3d(k,j,i,iens) = min( 5.0_fp , min( min( abs(udt) , abs(vdt) ) , abs(wdt) ) );
    });
    real dt = yakl::intrinsics::minval( dt3d );
    
    return dt;
  }



  // Initialize crap needed by recon()
  void init(std::string inFile, int num_tracers, DataManager &dm) {
    this->num_tracers = num_tracers;
    
    // Allocate device arrays for whether tracers are positive-definite or add mass
    tracer_pos       = bool1d("tracer_pos"      ,num_tracers);
    tracer_adds_mass = bool1d("tracer_adds_mass",num_tracers);

    // Inialize time step to zero, and dimensional splitting switch
    dimSwitch = true;

    // Read the YAML input file
    YAML::Node config = YAML::LoadFile(inFile);

    this->nens = config["nens"].as<int>();

    // Read in the vertical height cell interface locations
    std::string vcoords_file = config["vcoords"].as<std::string>();
    yakl::SimpleNetCDF nc;
    nc.open(vcoords_file);
    nz = nc.getDimSize("num_interfaces") - 1;
    real1d zint_in("zint_in",nz+1);
    nc.read(zint_in,"vertical_interfaces");
    nc.close();

    // Store vertical cell interface heights in the data manager
    dm.register_and_allocate<real>( "vertical_interface_height" , "vertical_interface_height" , {nz+1,nens} , {"zp1","nens"} );
    auto zint = dm.get<real,2>("vertical_interface_height");
    parallel_for( Bounds<2>(nz+1,nens) , YAKL_LAMBDA (int k, int iens) {
      zint(k,iens) = zint_in(k);
    });

    dm.register_and_allocate<real>( "vertical_midpoint_height" , "vertical_midpoint_heignt" , {nz,nens} , {"z","nens"} );
    auto zmid = dm.get<real,2>("vertical_midpoint_height");
    parallel_for( Bounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
      zmid(k,iens) = 0.5_fp*(zint_in(k) + zint_in(k+1));
    });

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

    dzmin = yakl::intrinsics::minval( dz );

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
    out_prefix = config["out_prefix"].as<std::string>();

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
    hyDensCells      = real2d("hyDensCells       ",nz,nens);
    hyPressureCells  = real2d("hyPressureCells   ",nz,nens);
    hyThetaCells     = real2d("hyThetaCells      ",nz,nens);
    hyDensThetaCells = real2d("hyDensThetaCells  ",nz,nens);
    hyDensGLL        = real3d("hyDensGLL         ",nz,ngll,nens);
    hyPressureGLL    = real3d("hyPressureGLL     ",nz,ngll,nens);
    hyThetaGLL       = real3d("hyThetaGLL        ",nz,ngll,nens);
    hyDensThetaGLL   = real3d("hyDensThetaGLL    ",nz,ngll,nens);
    pressure         = real4d("pressure",nz,ny,nx,nens);

    memset(pressure,0._fp);

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

    std::cout << "nx: " << nx << "\n";
    std::cout << "ny: " << ny << "\n";
    std::cout << "nz: " << nz << "\n";
    std::cout << "xlen (m): " << xlen << "\n";
    std::cout << "ylen (m): " << ylen << "\n";
    std::cout << "zlen (m): " << zlen.createHostCopy()(0) << "\n";
    std::cout << "Vertical coordinates file: " << vcoords_file << "\n";
    std::cout << "Simulation time (s): " << config["simTime"].as<real>() << "\n";
    std::cout << "Vertical interface heights: ";
    auto zint_host = zint.createHostCopy();
    for (int k=0; k < nz+1; k++) {
      std::cout << zint_host(k,0) << "  ";
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

                state(idT,hs+k,hs+j,hs+i,iens) += tp * wt;
              }
            }
          }
        }
      });

      init_tracers( dm , micro );

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

      int idWV = micro.tracer_index_vapor;
      real5d tracers = dm.get<real,5>("dynamics_tracers");

      parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
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
              state  (idU ,hs+k,hs+j,hs+i,iens) += uvel       * factor;
              state  (idV ,hs+k,hs+j,hs+i,iens) += vvel       * factor;
              state  (idW ,hs+k,hs+j,hs+i,iens) += wvel       * factor;
              state  (idT ,hs+k,hs+j,hs+i,iens) += theta_pert * factor;
              tracers(idWV,hs+k,hs+j,hs+i,iens) += dens_vap / dens * factor;
            }
          }
        }
      });

    } // if (data_spec == DATA_SPEC_SUPERCELL)
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
    if (splitIndex == 3) {
      computePressureTendencies( state , stateTend , dt );
      memset( tracerTend , 0._fp );
    }
  } // computeTendencies



  void switch_directions() {
    dimSwitch = ! dimSwitch;
  }



  void computePressureTendencies( real5d &state , real5d &stateTend , real dtglob ) {
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
    auto &pressure                = this->pressure               ;
    auto &bc_x                    = this->bc_x                   ;
    auto &bc_y                    = this->bc_y                   ;
    auto &bc_z                    = this->bc_z                   ;
    auto &Rd                      = this->Rd                     ;
    auto &cp                      = this->cp                     ;
    auto &gamma                   = this->gamma                  ;
    auto &p0                      = this->p0                     ;
    auto &C0                      = this->C0                     ;

    real constexpr c_s = 350;
    real dt;
    if (sim2d) {
      dt = 0.4_fp * min( dx/c_s , dzmin/c_s );
    } else {
      dt = 0.4_fp * min( min( dx/c_s , dy/c_s ) , dzmin/c_s );
    }

    real4d rho_u_new("rho_u_new",nz,ny,nx,nens);
    real4d rho_v_new("rho_v_new",nz,ny,nx,nens);
    real4d rho_w_new("rho_w_new",nz,ny,nx,nens);
    real4d pressure_tend ("pressure_tend" ,nz,ny,nx,nens);
    real4d rho_u_new_tend("rho_u_new_tend",nz,ny,nx,nens);
    real4d rho_v_new_tend("rho_v_new_tend",nz,ny,nx,nens);
    real4d rho_w_new_tend("rho_w_new_tend",nz,ny,nx,nens);
    real4d abs_div("abs_div",nz,ny,nx,nens);

    parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      rho_u_new(k,j,i,iens) = hyDensCells(k,iens) * state(idU,hs+k,hs+j,hs+i,iens);
      rho_v_new(k,j,i,iens) = hyDensCells(k,iens) * state(idV,hs+k,hs+j,hs+i,iens);
      rho_w_new(k,j,i,iens) = hyDensCells(k,iens) * state(idW,hs+k,hs+j,hs+i,iens);
      pressure (k,j,i,iens) = 0;
    });

    for (int iter=0; iter < 20000; iter++) {
      parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        int im1 = i-1;  if (im1 < 0   ) im1 = nx-1;
        int ip1 = i+1;  if (ip1 > nx-1) ip1 = 0;
        int jm1 = j-1;  if (jm1 < 0   ) jm1 = ny-1;
        int jp1 = j+1;  if (jp1 > ny-1) jp1 = 0;
        int km1 = k-1;  if (km1 < 0   ) km1 = 0;
        int kp1 = k+1;  if (kp1 > nz-1) kp1 = nz-1;
        // x-direction fluxes
        real p_im1  = pressure (k,j,im1,iens);
        real p_i    = pressure (k,j,i  ,iens);
        real p_ip1  = pressure (k,j,ip1,iens);
        real ru_im1 = rho_u_new(k,j,im1,iens);
        real ru_i   = rho_u_new(k,j,i  ,iens);
        real ru_ip1 = rho_u_new(k,j,ip1,iens);
        real p_x_L = (p_i   + p_im1) / 2       + c_s/2  * (ru_im1 - ru_i  );
        real p_x_R = (p_ip1 + p_i  ) / 2       + c_s/2  * (ru_i   - ru_ip1);
        real ru_L  = (p_im1 - p_i  ) / (2*c_s) + 0.5_fp * (ru_i   + ru_im1);
        real ru_R  = (p_i   - p_ip1) / (2*c_s) + 0.5_fp * (ru_ip1 + ru_i  );
        // y-direction fluxes
        real p_jm1  = pressure (k,jm1,i,iens);
        real p_j    = pressure (k,j  ,i,iens);
        real p_jp1  = pressure (k,jp1,i,iens);
        real rv_jm1 = rho_v_new(k,jm1,i,iens);
        real rv_j   = rho_v_new(k,j  ,i,iens);
        real rv_jp1 = rho_v_new(k,jp1,i,iens);
        real p_y_L = (p_j   + p_jm1) / 2       + c_s/2  * (rv_jm1 - rv_j  );
        real p_y_R = (p_jp1 + p_j  ) / 2       + c_s/2  * (rv_j   - rv_jp1);
        real rv_L  = (p_jm1 - p_j  ) / (2*c_s) + 0.5_fp * (rv_j   + rv_jm1);
        real rv_R  = (p_j   - p_jp1) / (2*c_s) + 0.5_fp * (rv_jp1 + rv_j  );
        // z-direction fluxes
        real p_km1  = pressure (km1,j,i,iens);
        real p_k    = pressure (k  ,j,i,iens);
        real p_kp1  = pressure (kp1,j,i,iens);
        real rw_km1 = rho_w_new(km1,j,i,iens);
        real rw_k   = rho_w_new(k  ,j,i,iens);
        real rw_kp1 = rho_w_new(kp1,j,i,iens);
        if (k == 0   ) rw_km1 = 0;
        if (k == nz-1) rw_kp1 = 0;
        real p_z_L = (p_k   + p_km1) / 2       + c_s/2  * (rw_km1 - rw_k  );
        real p_z_R = (p_kp1 + p_k  ) / 2       + c_s/2  * (rw_k   - rw_kp1);
        real rw_L  = (p_km1 - p_k  ) / (2*c_s) + 0.5_fp * (rw_k   + rw_km1);
        real rw_R  = (p_k   - p_kp1) / (2*c_s) + 0.5_fp * (rw_kp1 + rw_k  );
        if (k == 0   ) rw_L  = 0;
        if (k == nz-1) rw_R  = 0;
        // Perform the update
        if (sim2d) {
          abs_div(k,j,i,iens) = abs( (ru_R-ru_L)/dx + (rw_R-rw_L)/dz(k,iens) );
          pressure_tend(k,j,i,iens) = -c_s*c_s * ( (ru_R-ru_L)/dx + (rw_R-rw_L)/dz(k,iens) );
        } else {
          pressure_tend(k,j,i,iens) = -c_s*c_s * ( (ru_R-ru_L)/dx + (rv_R-rv_L)/dy + (rw_R-rw_L)/dz(k,iens) );
        }
        // std::cout << p_x_L << " , " << p_x_R << " , " << p_y_L << " , " << p_y_R << " , " << p_z_L << " , " << p_z_R << "\n";
        rho_u_new_tend(k,j,i,iens) = -(p_x_R - p_x_L) / (dx);
        rho_v_new_tend(k,j,i,iens) = -(p_y_R - p_y_L) / (dy);
        rho_w_new_tend(k,j,i,iens) = -(p_z_R - p_z_L) / (dz(k,iens));

        // if (sim2d) {
        //   abs_div(k,j,i,iens) = abs( (ru_ip1-ru_im1)/(2*dx) + (rw_kp1-rw_km1)/(2*dz(k,iens)) );
        //   pressure_tend(k,j,i,iens) = -c_s*c_s * ( (ru_ip1-ru_im1)/(2*dx) + (rw_kp1-rw_km1)/(2*dz(k,iens)) );
        // } else {
        //   pressure_tend(k,j,i,iens) = -c_s*c_s * ( (ru_ip1-ru_im1)/(2*dx) + (rv_jp1-rv_jm1)/(2*dy) + (rw_kp1-rw_km1)/(2*dz(k,iens)) );
        // }
        // std::cout << p_x_L << " , " << p_x_R << " , " << p_y_L << " , " << p_y_R << " , " << p_z_L << " , " << p_z_R << "\n";
        // rho_u_new_tend(k,j,i,iens) = -(p_ip1 - p_im1) / (2*dx);
        // rho_v_new_tend(k,j,i,iens) = -(p_jp1 - p_jm1) / (2*dy);
        // rho_w_new_tend(k,j,i,iens) = -(p_kp1 - p_km1) / (2*dz(k,iens));
      });

      if (iter == 0) std::cout << "Starting divergence: " << yakl::intrinsics::sum(abs_div) << "\n";

      parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        pressure (k,j,i,iens) += dt * pressure_tend(k,j,i,iens);
        rho_u_new(k,j,i,iens) = hyDensCells(k,iens) * state(idU,hs+k,hs+j,hs+i,iens) + dt * rho_u_new_tend(k,j,i,iens);
        rho_v_new(k,j,i,iens) = hyDensCells(k,iens) * state(idV,hs+k,hs+j,hs+i,iens) + dt * rho_v_new_tend(k,j,i,iens);
        rho_w_new(k,j,i,iens) = hyDensCells(k,iens) * state(idW,hs+k,hs+j,hs+i,iens) + dt * rho_w_new_tend(k,j,i,iens);
      });
    }

    std::cout << "Ending divergence: " << yakl::intrinsics::sum(abs_div) << "\n";

    parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      stateTend(idU,k,j,i,iens) = ( rho_u_new(k,j,i,iens)/hyDensCells(k,iens) - state(idU,hs+k,hs+j,hs+i,iens) ) / dtglob;
      if (! sim2d) {
        stateTend(idV,k,j,i,iens) = ( rho_v_new(k,j,i,iens)/hyDensCells(k,iens) - state(idV,hs+k,hs+j,hs+i,iens) ) / dtglob;
      } else {
        stateTend(idV,k,j,i,iens) = 0;
      }
      stateTend(idW,k,j,i,iens) = ( rho_w_new(k,j,i,iens)/hyDensCells(k,iens) - state(idW,hs+k,hs+j,hs+i,iens) ) / dtglob;
      
      stateTend(idT,k,j,i,iens) = 0;
    });
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
    YAKL_SCOPE( wenoRecon               , this->wenoRecon              );
    YAKL_SCOPE( idl                     , this->idl                    );
    YAKL_SCOPE( sigma                   , this->sigma                  );
    YAKL_SCOPE( hyDensCells             , this->hyDensCells            );
    YAKL_SCOPE( hyDensThetaCells        , this->hyDensThetaCells       );
    YAKL_SCOPE( hyThetaCells            , this->hyThetaCells           );
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
          state(l,hs+k,hs+j,      ii,iens) = state(l,hs+k,hs+j,nx+ii,iens);
          state(l,hs+k,hs+j,hs+nx+ii,iens) = state(l,hs+k,hs+j,hs+ii,iens);
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
      SArray<real,2,nAder,ngll> u_DTs;
      real hy_dens = hyDensCells(k,iens);

      { // BEGIN: Reconstruct, time-average, and store the state and fluxes
        ////////////////////////////////////////////////////////////////
        // Reconstruct rho, u, v, w, theta
        ////////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> v_DTs , w_DTs , t_DTs;
        { // BEGIN: Reconstruct the state
          SArray<real,1,ord> stencil;

          // u values and derivatives
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idU,hs+k,hs+j,i+ii,iens); }
          reconstruct_gll_values( stencil , u_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // v
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idV,hs+k,hs+j,i+ii,iens); }
          reconstruct_gll_values( stencil , v_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // w
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idW,hs+k,hs+j,i+ii,iens); }
          reconstruct_gll_values( stencil , w_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // theta
          for (int ii=0; ii < ord; ii++) { stencil(ii) = state(idT,hs+k,hs+j,i+ii,iens); }
          reconstruct_gll_values( stencil , t_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
          for (int ii=0; ii < ngll; ii++) { t_DTs(0,ii) += hyThetaCells(k,iens); } // Add hydrostasis back on
        } // END: Reconstruct the state

        ///////////////////////////////////////////////////////////////
        // Compute other values needed for centered tendencies and DTs
        ///////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> uu_DTs , uv_DTs , uw_DTs , ut_DTs;
        for (int ii=0; ii < ngll; ii++) {
          real u = u_DTs(0,ii);
          real v = v_DTs(0,ii);
          real w = w_DTs(0,ii);
          real t = t_DTs(0,ii);
          uu_DTs(0,ii) = u*u;
          uv_DTs(0,ii) = u*v;
          uw_DTs(0,ii) = u*w;
          ut_DTs(0,ii) = u*t;
        }

        //////////////////////////////////////////
        // Compute time derivatives if necessary
        //////////////////////////////////////////
        if (nAder > 1) {
          diffTransformEulerConsX( hy_dens , u_DTs , v_DTs , w_DTs , t_DTs , uu_DTs , uv_DTs , uw_DTs ,
                                   ut_DTs , derivMatrix , dx );
        }

        //////////////////////////////////////////
        // Time average if necessary
        //////////////////////////////////////////
        // density and momentum can't be overwritten because they will be used for tracers
        SArray<real,1,ngll> u_tavg;
        if (timeAvg) {
          compute_timeAvg( u_DTs , u_tavg , dt );
          compute_timeAvg( v_DTs          , dt );
          compute_timeAvg( w_DTs          , dt );
          compute_timeAvg( t_DTs          , dt );
        } else {
          for (int ii=0; ii < ngll; ii++) {
            u_tavg(ii) = u_DTs(0,ii);
          }
        }

        //////////////////////////////////////////
        // Store cell edge estimates of the state
        //////////////////////////////////////////
        // Left interface
        stateLimits(idU,1,k,j,i  ,iens) = u_tavg (0     );
        stateLimits(idV,1,k,j,i  ,iens) = v_DTs(0,0     );
        stateLimits(idW,1,k,j,i  ,iens) = w_DTs(0,0     );
        stateLimits(idT,1,k,j,i  ,iens) = t_DTs(0,0     );
        // Right interface
        stateLimits(idU,0,k,j,i+1,iens) = u_tavg (ngll-1);
        stateLimits(idV,0,k,j,i+1,iens) = v_DTs(0,ngll-1);
        stateLimits(idW,0,k,j,i+1,iens) = w_DTs(0,ngll-1);
        stateLimits(idT,0,k,j,i+1,iens) = t_DTs(0,ngll-1);
      } // END: Reconstruct, time-average, and store the state and fluxes

      // r_DTs and ru_DTs still exist and are computed
      { // BEGIN: Reconstruct, time-average, and store tracer fluxes
        // Only process one tracer at a time to save on local memory / register requirements
        for (int tr=0; tr < num_tracers; tr++) {
          SArray<real,2,nAder,ngll> t_DTs; // Density * tracer
          { // BEGIN: Reconstruct the tracer
            SArray<real,1,ord> stencil;
            for (int ii=0; ii < ord; ii++) { stencil(ii) = tracers(tr,hs+k,hs+j,i+ii,iens); }
            reconstruct_gll_values( stencil , t_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
            if (tracer_pos(tr)) {
              for (int ii=0; ii < ngll; ii++) { t_DTs(0,ii) = max( 0._fp , t_DTs(0,ii) ); }
            }
          } // END: Reconstruct the tracer

          // Compute the tracer flux
          SArray<real,2,nAder,ngll> ut_DTs; // Density * uwind * tracer
          for (int ii=0; ii < ngll; ii++) {
            ut_DTs(0,ii) = u_DTs(0,ii) * t_DTs(0,ii);
          }

          //////////////////////////////////////////
          // Compute time derivatives if necessary
          //////////////////////////////////////////
          if (nAder > 1) {
            diffTransformTracerXY( hy_dens , u_DTs , t_DTs , ut_DTs , derivMatrix , dx );
          }

          //////////////////////////////////////////
          // Time average if necessary
          //////////////////////////////////////////
          if (timeAvg) {
            compute_timeAvg( t_DTs , dt );
          }
          if (tracer_pos(tr)) {
            for (int ii=0; ii < ngll; ii++) { t_DTs(0,ii) = max( 0._fp , t_DTs(0,ii) ); }
          }

          ////////////////////////////////////////////////////////////
          // Store cell edge estimates of the tracer
          ////////////////////////////////////////////////////////////
          tracerLimits(tr,1,k,j,i  ,iens) = t_DTs(0,0     ); // Left interface
          tracerLimits(tr,0,k,j,i+1,iens) = t_DTs(0,ngll-1); // Right interface
        }
      } // END: Reconstruct, time-average, and store tracer fluxes

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
      real hy_dens = hyDensCells(k,iens);
      real u_L = stateLimits(idU,0,k,j,i,iens);   real u_R = stateLimits(idU,1,k,j,i,iens);
      real v_L = stateLimits(idV,0,k,j,i,iens);   real v_R = stateLimits(idV,1,k,j,i,iens);
      real w_L = stateLimits(idW,0,k,j,i,iens);   real w_R = stateLimits(idW,1,k,j,i,iens);
      real t_L = stateLimits(idT,0,k,j,i,iens);   real t_R = stateLimits(idT,1,k,j,i,iens);
      // Compute average state
      real u = 0.5_fp * (u_L + u_R);

      if (u > 0) {
        stateFlux(idU,k,j,i,iens) = hy_dens*u_L*u_L;
        stateFlux(idV,k,j,i,iens) = hy_dens*u_L*v_L;
        stateFlux(idW,k,j,i,iens) = hy_dens*u_L*w_L;
        stateFlux(idT,k,j,i,iens) = hy_dens*u_L*t_L;
        for (int tr=0; tr < num_tracers; tr++) {
          tracerFlux(tr,k,j,i,iens) = hy_dens * u_L * tracerLimits(tr,0,k,j,i,iens);
        }
      } else {
        stateFlux(idU,k,j,i,iens) = hy_dens*u_R*u_R;
        stateFlux(idV,k,j,i,iens) = hy_dens*u_R*v_R;
        stateFlux(idW,k,j,i,iens) = hy_dens*u_R*w_R;
        stateFlux(idT,k,j,i,iens) = hy_dens*u_R*t_R;
        for (int tr=0; tr < num_tracers; tr++) {
          tracerFlux(tr,k,j,i,iens) = hy_dens * u_R * tracerLimits(tr,1,k,j,i,iens);
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Limit the tracer fluxes for positivity
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<5>(num_tracers,nz,ny,nx+1,nens) , YAKL_LAMBDA (int tr, int k, int j, int i, int iens) {
      real constexpr eps = 1.e-10;
      real u = 0.5_fp * ( stateLimits(idU,0,k,j,i,iens) + stateLimits(idU,1,k,j,i,iens) );
      // Solid wall BCs mean u == 0 at boundaries, so we assume periodic if u != 0
      if (tracer_pos(tr)) {
        // Compute and apply the flux reduction factor of the upwind cell
        if      (u > 0) {
          // upwind is to the left of this interface
          int ind_i = i-1;
          if (ind_i == -1) ind_i = nx-1;
          real f1 = min( tracerFlux(tr,k,j,ind_i  ,iens) , 0._fp );
          real f2 = max( tracerFlux(tr,k,j,ind_i+1,iens) , 0._fp );
          real fluxOut = dt*(f2-f1)/dx;
          real dens = hyDensCells(k,iens);
          tracerFlux(tr,k,j,i,iens) *= min( 1._fp , tracers(tr,hs+k,hs+j,hs+ind_i,iens) * dens / (fluxOut + eps) );
        } else if (u < 0) {
          // upwind is to the right of this interface
          int ind_i = i;
          if (ind_i == nx) ind_i = 0;
          real f1 = min( tracerFlux(tr,k,j,ind_i  ,iens) , 0._fp );
          real f2 = max( tracerFlux(tr,k,j,ind_i+1,iens) , 0._fp );
          real fluxOut = dt*(f2-f1)/dx;
          real dens = hyDensCells(k,iens);
          tracerFlux(tr,k,j,i,iens) *= min( 1._fp , tracers(tr,hs+k,hs+j,hs+ind_i,iens) * dens / (fluxOut + eps) );
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA(int k, int j, int i, int iens) {
      for (int l = 0; l < num_state; l++) {
        if (sim2d && l == idV) {
          stateTend(l,k,j,i,iens) = 0;
        } else {
          stateTend(l,k,j,i,iens) = - ( stateFlux(l,k,j,i+1,iens) - stateFlux(l,k,j,i,iens) ) / dx / hyDensCells(k,iens);
        }
      }
      for (int l = 0; l < num_tracers; l++) {
        // Compute tracer tendency
        tracerTend(l,k,j,i,iens) = - ( tracerFlux(l,k,j,i+1,iens) - tracerFlux(l,k,j,i  ,iens) ) / dx / hyDensCells(k,iens);
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
    YAKL_SCOPE( hyThetaCells            , this->hyThetaCells           );
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

    // Populate the halos
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

    // Loop through all cells, reconstruct in y-direction, compute centered tendencies, store cell-edge state estimates
    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      // These are needed by the tracers
      SArray<real,2,nAder,ngll> v_DTs;
      real hy_dens = hyDensCells(k,iens);

      { // BEGIN: Reconstruct, time-average, and store state and sate fluxes
        ////////////////////////////////////////////////////////////////
        // Reconstruct rho, u, v, w, theta
        ////////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> u_DTs , w_DTs , t_DTs;
        {
          SArray<real,1,ord> stencil;

          // u values and derivatives
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idU,hs+k,j+jj,hs+i,iens); }
          reconstruct_gll_values( stencil , u_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // v
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idV,hs+k,j+jj,hs+i,iens); }
          reconstruct_gll_values( stencil , v_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // w
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idW,hs+k,j+jj,hs+i,iens); }
          reconstruct_gll_values( stencil , w_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // theta
          for (int jj=0; jj < ord; jj++) { stencil(jj) = state(idT,hs+k,j+jj,hs+i,iens); }
          reconstruct_gll_values( stencil , t_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
          for (int jj=0; jj < ngll; jj++) { t_DTs(0,jj) += hyThetaCells(k,iens); } // Add hydrostasis back on
        }

        ///////////////////////////////////////////////////////////////
        // Compute other values needed for centered tendencies and DTs
        ///////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> vu_DTs , vv_DTs , vw_DTs , vt_DTs;
        for (int jj=0; jj < ngll; jj++) {
          real u = u_DTs(0,jj);
          real v = v_DTs(0,jj);
          real w = w_DTs(0,jj);
          real t = t_DTs(0,jj);
          vu_DTs(0,jj) = v*u;
          vv_DTs(0,jj) = v*v;
          vw_DTs(0,jj) = v*w;
          vt_DTs(0,jj) = v*t;
        }

        //////////////////////////////////////////
        // Compute time derivatives if necessary
        //////////////////////////////////////////
        if (nAder > 1) {
          diffTransformEulerConsY( hy_dens , u_DTs , v_DTs , w_DTs , t_DTs , vu_DTs , vv_DTs , vw_DTs ,
                                   vt_DTs , derivMatrix , dy );
        }

        //////////////////////////////////////////
        // Time average if necessary
        //////////////////////////////////////////
        // Don't overwrite r and rv because we need them for tracers
        SArray<real,1,ngll> v_tavg;
        if (timeAvg) {
          compute_timeAvg( u_DTs          , dt );
          compute_timeAvg( v_DTs , v_tavg , dt );
          compute_timeAvg( w_DTs          , dt );
          compute_timeAvg( t_DTs          , dt );
        } else {
          for (int jj=0; jj < ngll; jj++) {
            v_tavg(jj) = v_DTs(0,jj);
          }
        }

        //////////////////////////////////////////
        // Store cell edge estimates of the state
        //////////////////////////////////////////
        // Left interface
        stateLimits(idU,1,k,j  ,i,iens) = u_DTs(0,0     );
        stateLimits(idV,1,k,j  ,i,iens) = v_tavg (0     );
        stateLimits(idW,1,k,j  ,i,iens) = w_DTs(0,0     );
        stateLimits(idT,1,k,j  ,i,iens) = t_DTs(0,0     );
        // Right interface       
        stateLimits(idU,0,k,j+1,i,iens) = u_DTs(0,ngll-1);
        stateLimits(idV,0,k,j+1,i,iens) = v_tavg (ngll-1);
        stateLimits(idW,0,k,j+1,i,iens) = w_DTs(0,ngll-1);
        stateLimits(idT,0,k,j+1,i,iens) = t_DTs(0,ngll-1);
      } // END: Reconstruct, time-average, and store state and sate fluxes

      // r_DTs and rv_DTs still exist and are computed
      { // BEGIN: Reconstruct, time-average, and store tracer fluxes
        // Only process one tracer at a time to save on local memory / register requirements
        for (int tr=0; tr < num_tracers; tr++) {
          SArray<real,2,nAder,ngll> t_DTs; // Density * tracer
          { // BEGIN: Reconstruct th tracer
            SArray<real,1,ord> stencil;
            for (int jj=0; jj < ord;jj++) { stencil(jj) = tracers(tr,hs+k,j+jj,hs+i,iens); }
            reconstruct_gll_values( stencil , t_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
            if (tracer_pos(tr)) {
              for (int jj=0; jj < ngll; jj++) { t_DTs(0,jj) = max( 0._fp , t_DTs(0,jj) ); }
            }
          } // END: Reconstruct the tracer

          // Compute the tracer flux
          SArray<real,2,nAder,ngll> vt_DTs; // Density * vwind * tracer
          for (int jj=0; jj < ngll; jj++) {
            vt_DTs(0,jj) = v_DTs(0,jj) * t_DTs(0,jj);
          }

          //////////////////////////////////////////
          // Compute time derivatives if necessary
          //////////////////////////////////////////
          if (nAder > 1) {
            diffTransformTracerXY( hy_dens , v_DTs , t_DTs , vt_DTs , derivMatrix , dy );
          }

          //////////////////////////////////////////
          // Time average if necessary
          //////////////////////////////////////////
          if (timeAvg) {
            compute_timeAvg( t_DTs  , dt );
          }
          if (tracer_pos(tr)) {
            for (int jj=0; jj < ngll; jj++) { t_DTs(0,jj) = max( 0._fp , t_DTs(0,jj) ); }
          }

          ////////////////////////////////////////////////////////////
          // Store cell edge estimates of the tracer
          ////////////////////////////////////////////////////////////
          tracerLimits(tr,1,k,j  ,i,iens) = t_DTs (0,0     ); // Left interface
          tracerLimits(tr,0,k,j+1,i,iens) = t_DTs (0,ngll-1); // Right interface
        }
      } // END: Reconstruct, time-average, and store tracer fluxes

    });

    ////////////////////////////////////////////////
    // BCs for the state edge estimates
    ////////////////////////////////////////////////
    parallel_for( SimpleBounds<3>(nz,nx,nens) , YAKL_LAMBDA (int k, int i, int iens) {
      for (int l=0; l < num_state; l++) {
        if        (bc_y == BC_PERIODIC) {
          stateLimits(l,0,k,0 ,i,iens) = stateLimits(l,0,k,ny,i,iens);
          stateLimits(l,1,k,ny,i,iens) = stateLimits(l,1,k,0 ,i,iens);
        } else if (bc_y == BC_WALL    ) {
          stateLimits(l,0,k,0 ,i,iens) = stateLimits(l,1,k,0 ,i,iens);
          stateLimits(l,1,k,ny,i,iens) = stateLimits(l,0,k,ny,i,iens);
        }
      }
      for (int l=0; l < num_tracers; l++) {
        if        (bc_y == BC_PERIODIC) {
          tracerLimits(l,0,k,0 ,i,iens) = tracerLimits(l,0,k,ny,i,iens);
          tracerLimits(l,1,k,ny,i,iens) = tracerLimits(l,1,k,0 ,i,iens);
        } else if (bc_y == BC_WALL    ) {
          tracerLimits(l,0,k,0 ,i,iens) = tracerLimits(l,1,k,0 ,i,iens);
          tracerLimits(l,1,k,ny,i,iens) = tracerLimits(l,0,k,ny,i,iens);
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the upwind fluxes
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<4>(nz,ny+1,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      // Get left and right state
      real hy_dens = hyDensCells(k,iens);
      real u_L = stateLimits(idU,0,k,j,i,iens);   real u_R = stateLimits(idU,1,k,j,i,iens);
      real v_L = stateLimits(idV,0,k,j,i,iens);   real v_R = stateLimits(idV,1,k,j,i,iens);
      real w_L = stateLimits(idW,0,k,j,i,iens);   real w_R = stateLimits(idW,1,k,j,i,iens);
      real t_L = stateLimits(idT,0,k,j,i,iens);   real t_R = stateLimits(idT,1,k,j,i,iens);
      // Compute average state
      real v = 0.5_fp * (v_L + v_R);

      if (v > 0) {
        stateFlux(idU,k,j,i,iens) = hy_dens*v_L*u_L;
        stateFlux(idV,k,j,i,iens) = hy_dens*v_L*v_L;
        stateFlux(idW,k,j,i,iens) = hy_dens*v_L*w_L;
        stateFlux(idT,k,j,i,iens) = hy_dens*v_L*t_L;
        for (int tr=0; tr < num_tracers; tr++) {
          tracerFlux(tr,k,j,i,iens) = hy_dens * v_L * tracerLimits(tr,0,k,j,i,iens);
        }
      } else {
        stateFlux(idU,k,j,i,iens) = hy_dens*v_R*u_R;
        stateFlux(idV,k,j,i,iens) = hy_dens*v_R*v_R;
        stateFlux(idW,k,j,i,iens) = hy_dens*v_R*w_R;
        stateFlux(idT,k,j,i,iens) = hy_dens*v_R*t_R;
        for (int tr=0; tr < num_tracers; tr++) {
          tracerFlux(tr,k,j,i,iens) = hy_dens * v_R * tracerLimits(tr,1,k,j,i,iens);
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Limit the tracer fluxes for positivity
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<5>(num_tracers,nz,ny+1,nx,nens) , YAKL_LAMBDA (int tr, int k, int j, int i, int iens) {
      real constexpr eps = 1.e-10;
      real v = 0.5_fp * ( stateLimits(idV,0,k,j,i,iens) + stateLimits(idV,1,k,j,i,iens) );
      // Solid wall BCs mean u == 0 at boundaries, so we assume periodic if u != 0
      if (tracer_pos(tr)) {
        // Compute and apply the flux reduction factor of the upwind cell
        if      (v > 0) {
          // upwind is to the left of this interface
          int ind_j = j-1;
          if (ind_j == -1) ind_j = ny-1;
          real f1 = min( tracerFlux(tr,k,ind_j  ,i,iens) , 0._fp );
          real f2 = max( tracerFlux(tr,k,ind_j+1,i,iens) , 0._fp );
          real fluxOut = dt*(f2-f1)/dy;
          real dens = hyDensCells(k,iens);
          tracerFlux(tr,k,j,i,iens) *= min( 1._fp , tracers(tr,hs+k,hs+ind_j,hs+i,iens) * dens / (fluxOut + eps) );
        } else if (v < 0) {
          // upwind is to the right of this interface
          int ind_j = j;
          if (ind_j == ny) ind_j = 0;
          real f1 = min( tracerFlux(tr,k,ind_j  ,i,iens) , 0._fp );
          real f2 = max( tracerFlux(tr,k,ind_j+1,i,iens) , 0._fp );
          real fluxOut = dt*(f2-f1)/dy;
          real dens = hyDensCells(k,iens);
          tracerFlux(tr,k,j,i,iens) *= min( 1._fp , tracers(tr,hs+k,hs+ind_j,hs+i,iens) * dens / (fluxOut + eps) );
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA(int k, int j, int i, int iens) {
      for (int l=0; l < num_state; l++) {
        stateTend(l,k,j,i,iens) = - ( stateFlux(l,k,j+1,i,iens) - stateFlux(l,k,j,i,iens) ) / dy / hyDensCells(k,iens);
      }
      for (int l=0; l < num_tracers; l++) {
        // Compute the tracer tendency
        tracerTend(l,k,j,i,iens) = - ( tracerFlux(l,k,j+1,i,iens) - tracerFlux(l,k,j,i,iens) ) / dy / hyDensCells(k,iens);
      }
    });
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
    YAKL_SCOPE( wenoRecon               , this->wenoRecon              );
    YAKL_SCOPE( idl                     , this->idl                    );
    YAKL_SCOPE( sigma                   , this->sigma                  );
    YAKL_SCOPE( hyDensCells             , this->hyDensCells            );
    YAKL_SCOPE( hyDensGLL               , this->hyDensGLL              );
    YAKL_SCOPE( hyDensThetaGLL          , this->hyDensThetaGLL         );
    YAKL_SCOPE( hyPressureGLL           , this->hyPressureGLL          );
    YAKL_SCOPE( hyThetaGLL              , this->hyThetaGLL             );
    YAKL_SCOPE( hyThetaCells            , this->hyThetaCells           );
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
      // We need these to persist to evolve tracers with ADER
      SArray<real,2,nAder,ngll> w_DTs;
      SArray<real,1,ngll> hy_dens;
      for (int kk=0; kk < ngll; kk++) {
        hy_dens(kk) = hyDensGLL(k,kk,iens);
      }

      { // BEGIN: reconstruct, time-avg, and store state & state fluxes
        ////////////////////////////////////////////////////////////////
        // Reconstruct rho, u, v, w, theta
        ////////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> u_DTs , v_DTs , t_DTs;
        {
          SArray<real,1,ord> stencil;

          // u values and derivatives
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idU,k+kk,hs+j,hs+i,iens); }
          reconstruct_gll_values( stencil , u_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // v
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idV,k+kk,hs+j,hs+i,iens); }
          reconstruct_gll_values( stencil , v_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );

          // w
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idW,k+kk,hs+j,hs+i,iens); }
          reconstruct_gll_values( stencil , w_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_winds );
          if (bc_z == BC_WALL) {
            if (k == nz-1) w_DTs(0,ngll-1) = 0;
            if (k == 0   ) w_DTs(0,0     ) = 0;
          }

          // theta
          for (int kk=0; kk < ord; kk++) { stencil(kk) = state(idT,k+kk,hs+j,hs+i,iens); }
          reconstruct_gll_values( stencil , t_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
          for (int kk=0; kk < ngll; kk++) { t_DTs(0,kk) += hyThetaGLL(k,kk,iens); } // Add hydrostasis back on
        }

        ///////////////////////////////////////////////////////////////
        // Compute other values needed for centered tendencies and DTs
        ///////////////////////////////////////////////////////////////
        SArray<real,2,nAder,ngll> wu_DTs , wv_DTs , ww_DTs , wt_DTs;
        for (int kk=0; kk < ngll; kk++) {
          real u = u_DTs(0,kk);
          real v = v_DTs(0,kk);
          real w = w_DTs(0,kk);
          real t = t_DTs(0,kk);
          wu_DTs(0,kk) = w*u;
          wv_DTs(0,kk) = w*v;
          ww_DTs(0,kk) = w*w;
          wt_DTs(0,kk) = w*t;
        }

        //////////////////////////////////////////
        // Compute time derivatives if necessary
        //////////////////////////////////////////
        if (nAder > 1) {
          diffTransformEulerConsZ( hy_dens , u_DTs , v_DTs , w_DTs , t_DTs , wu_DTs , wv_DTs , ww_DTs ,
                                   wt_DTs , derivMatrix , k , dz(k,iens) , bc_z , nz , iens );
        }

        //////////////////////////////////////////
        // Time average if necessary
        //////////////////////////////////////////
        // We can't alter density and momentum because they're needed for tracers later
        SArray<real,1,ngll> w_tavg;
        if (timeAvg) {
          compute_timeAvg( u_DTs          , dt );
          compute_timeAvg( v_DTs          , dt );
          compute_timeAvg( w_DTs , w_tavg , dt );
          compute_timeAvg( t_DTs          , dt );
        } else {
          for (int ii=0; ii < ngll; ii++) {
            w_tavg(ii) = w_DTs(0,ii);
          }
        }

        //////////////////////////////////////////
        // Store cell edge estimates of the state
        //////////////////////////////////////////
        // Left interface
        stateLimits(idU,1,k  ,j,i,iens) = u_DTs(0,0     );
        stateLimits(idV,1,k  ,j,i,iens) = v_DTs(0,0     );
        stateLimits(idW,1,k  ,j,i,iens) = w_tavg (0     );
        stateLimits(idT,1,k  ,j,i,iens) = t_DTs(0,0     );
        // Right interface       
        stateLimits(idU,0,k+1,j,i,iens) = u_DTs(0,ngll-1);
        stateLimits(idV,0,k+1,j,i,iens) = v_DTs(0,ngll-1);
        stateLimits(idW,0,k+1,j,i,iens) = w_tavg (ngll-1);
        stateLimits(idT,0,k+1,j,i,iens) = t_DTs(0,ngll-1);

        ////////////////////////////////////////////
        // Assign gravity source term
        ////////////////////////////////////////////
        real theta_avg = 0;
        for (int kk=0; kk < ngll; kk++) {
          theta_avg += t_DTs(0,kk) * gllWts_ngll(kk);
        }
        stateTend(idU,k,j,i,iens) = 0;
        stateTend(idV,k,j,i,iens) = 0;
        stateTend(idW,k,j,i,iens) = GRAV * (theta_avg - hyThetaCells(k,iens)) / theta_avg;
        stateTend(idT,k,j,i,iens) = 0;
      } // END: reconstruct, time-avg, and store state & state fluxes

      // r_DTs and rw_DTs still exist and are computed
      { // BEGIN: Reconstruct, time-average, and store tracer fluxes
        // Only process one tracer at a time to save on local memory / register requirements
        for (int tr=0; tr < num_tracers; tr++) {
          SArray<real,2,nAder,ngll> t_DTs; // Density * tracer
          { // BEGIN: Reconstruct the tracer
            SArray<real,1,ord> stencil;
            for (int kk=0; kk < ord; kk++) { stencil(kk) = tracers(tr,k+kk,hs+j,hs+i,iens); }
            reconstruct_gll_values( stencil , t_DTs , c2g , s2g , wenoRecon , idl , sigma , weno_scalars );
            if (tracer_pos(tr)) {
              for (int kk=0; kk < ngll; kk++) { t_DTs(0,kk) = max( 0._fp , t_DTs(0,kk) ); }
            }
          } // END: Reconstruct the tracer

          // Compute the tracer flux
          SArray<real,2,nAder,ngll> wt_DTs; // Density * wwind * tracer
          for (int kk=0; kk < ngll; kk++) {
            wt_DTs(0,kk) = w_DTs(0,kk) * t_DTs(0,kk);
          }

          //////////////////////////////////////////
          // Compute time derivatives if necessary
          //////////////////////////////////////////
          if (nAder > 1) {
            diffTransformTracerZ( hy_dens , w_DTs , t_DTs , wt_DTs , derivMatrix , dz(k,iens) );
          }

          //////////////////////////////////////////
          // Time average if necessary
          //////////////////////////////////////////
          if (timeAvg) {
            compute_timeAvg( t_DTs , dt );
          }
          if (tracer_pos(tr)) {
            for (int kk=0; kk < ngll; kk++) { t_DTs(0,kk) = max( 0._fp , t_DTs(0,kk) ); }
          }

          ////////////////////////////////////////////////////////////
          // Store cell edge estimates of the tracer
          ////////////////////////////////////////////////////////////
          tracerLimits(tr,1,k  ,j,i,iens) = t_DTs(0,0     ); // Left interface
          tracerLimits(tr,0,k+1,j,i,iens) = t_DTs(0,ngll-1); // Right interface
        }
      } // END: Reconstruct, time-average, and store tracer fluxes

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
      real dens;
      if (k <  nz) dens = hyDensGLL(k   ,0     ,iens);
      if (k == nz) dens = hyDensGLL(nz-1,ngll-1,iens);
      real u_L = stateLimits(idU,0,k,j,i,iens);   real u_R = stateLimits(idU,1,k,j,i,iens);
      real v_L = stateLimits(idV,0,k,j,i,iens);   real v_R = stateLimits(idV,1,k,j,i,iens);
      real w_L = stateLimits(idW,0,k,j,i,iens);   real w_R = stateLimits(idW,1,k,j,i,iens);
      real t_L = stateLimits(idT,0,k,j,i,iens);   real t_R = stateLimits(idT,1,k,j,i,iens);
      // Compute average state
      real w = 0.5_fp * (w_L + w_R);

      if (w > 0) {
        stateFlux(idU,k,j,i,iens) = dens*w_L*u_L;
        stateFlux(idV,k,j,i,iens) = dens*w_L*v_L;
        stateFlux(idW,k,j,i,iens) = dens*w_L*w_L;
        stateFlux(idT,k,j,i,iens) = dens*w_L*t_L;
        for (int tr=0; tr < num_tracers; tr++) {
          tracerFlux(tr,k,j,i,iens) = dens * w_L * tracerLimits(tr,0,k,j,i,iens);
        }
      } else {
        stateFlux(idU,k,j,i,iens) = dens*w_R*u_R;
        stateFlux(idV,k,j,i,iens) = dens*w_R*v_R;
        stateFlux(idW,k,j,i,iens) = dens*w_R*w_R;
        stateFlux(idT,k,j,i,iens) = dens*w_R*t_R;
        for (int tr=0; tr < num_tracers; tr++) {
          tracerFlux(tr,k,j,i,iens) = dens * w_R * tracerLimits(tr,1,k,j,i,iens);
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Limit the tracer fluxes for positivity
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<5>(num_tracers,nz+1,ny,nx,nens) , YAKL_LAMBDA (int tr, int k, int j, int i, int iens) {
      real constexpr eps = 1.e-10;
      real w = 0.5_fp * ( stateLimits(idW,0,k,j,i,iens) + stateLimits(idW,1,k,j,i,iens) );
      // Solid wall BCs mean w == 0 at boundaries
      if (tracer_pos(tr)) {
        // Compute and apply the flux reduction factor of the upwind cell
        if      (w > 0) {
          // upwind is to the left of this interface
          int ind_k = k-1;
          real f1 = min( tracerFlux(tr,ind_k  ,j,i,iens) , 0._fp );
          real f2 = max( tracerFlux(tr,ind_k+1,j,i,iens) , 0._fp );
          real fluxOut = dt*(f2-f1)/dz(k,iens);
          real dens = hyDensCells(ind_k,iens);
          tracerFlux(tr,k,j,i,iens) *= min( 1._fp , tracers(tr,hs+ind_k,hs+j,hs+i,iens) * dens / (fluxOut + eps) );
        } else if (w < 0) {
          // upwind is to the right of this interface
          int ind_k = k;
          real f1 = min( tracerFlux(tr,ind_k  ,j,i,iens) , 0._fp );
          real f2 = max( tracerFlux(tr,ind_k+1,j,i,iens) , 0._fp );
          real fluxOut = dt*(f2-f1)/dz(k,iens);
          real dens = hyDensCells(ind_k,iens);
          tracerFlux(tr,k,j,i,iens) *= min( 1._fp , tracers(tr,hs+ind_k,hs+j,hs+i,iens) * dens / (fluxOut + eps) );
        }
      }
    });

    //////////////////////////////////////////////////////////
    // Compute the tendencies
    //////////////////////////////////////////////////////////
    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA(int k, int j, int i, int iens) {
      for (int l=0; l < num_state; l++) {
        if (sim2d && l == idV) {
          stateTend(l,k,j,i,iens) = 0;
        } else {
          stateTend(l,k,j,i,iens) += - ( stateFlux(l,k+1,j,i,iens) - stateFlux(l,k,j,i,iens) ) / dz(k,iens) / hyDensCells(k,iens);
        }
      }
      for (int l=0; l < num_tracers; l++) {
        // Compute tracer tendency
        tracerTend(l,k,j,i,iens) = - ( tracerFlux(l,k+1,j,i,iens) - tracerFlux(l,k,j,i,iens) ) / dz(k,iens) / hyDensCells(k,iens);
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
    YAKL_SCOPE( pressure              , this->pressure             );
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

      for (int tr=0; tr < num_tracers; tr++) {
        parallel_for( SimpleBounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
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



  YAKL_INLINE void diffTransformEulerConsX( real hy_dens  ,
                                            SArray<real,2,nAder,ngll> &u ,
                                            SArray<real,2,nAder,ngll> &v ,
                                            SArray<real,2,nAder,ngll> &w ,
                                            SArray<real,2,nAder,ngll> &t ,
                                            SArray<real,2,nAder,ngll> &uu ,
                                            SArray<real,2,nAder,ngll> &uv ,
                                            SArray<real,2,nAder,ngll> &uw ,
                                            SArray<real,2,nAder,ngll> &ut ,
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
          df2_dx += deriv(s,ii) * uu(kt,s);
          df3_dx += deriv(s,ii) * uv(kt,s);
          df4_dx += deriv(s,ii) * uw(kt,s);
          df5_dx += deriv(s,ii) * ut(kt,s);
        }
        u(kt+1,ii) = -df2_dx/dx/(kt+1._fp);
        v(kt+1,ii) = -df3_dx/dx/(kt+1._fp);
        w(kt+1,ii) = -df4_dx/dx/(kt+1._fp);
        t(kt+1,ii) = -df5_dx/dx/(kt+1._fp);
      }

      // Compute u* at the next time level
      for (int ii=0; ii<ngll; ii++) {
        real tot_uu = 0;
        real tot_uv = 0;
        real tot_uw = 0;
        real tot_ut = 0;
        for (int ir=0; ir<=kt+1; ir++) {
          tot_uu += u(ir,ii) * u(kt+1-ir,ii);
          tot_uv += u(ir,ii) * v(kt+1-ir,ii);
          tot_uw += u(ir,ii) * w(kt+1-ir,ii);
          tot_ut += u(ir,ii) * t(kt+1-ir,ii);
        }
        uu(kt+1,ii) = tot_uu;
        uv(kt+1,ii) = tot_uv;
        uw(kt+1,ii) = tot_uw;
        ut(kt+1,ii) = tot_ut;
      }
    }
  }



  YAKL_INLINE void diffTransformEulerConsY( real hy_dens  ,
                                            SArray<real,2,nAder,ngll> &u ,
                                            SArray<real,2,nAder,ngll> &v ,
                                            SArray<real,2,nAder,ngll> &w ,
                                            SArray<real,2,nAder,ngll> &t ,
                                            SArray<real,2,nAder,ngll> &vu ,
                                            SArray<real,2,nAder,ngll> &vv ,
                                            SArray<real,2,nAder,ngll> &vw ,
                                            SArray<real,2,nAder,ngll> &vt ,
                                            SArray<real,2,ngll,ngll> const &deriv , 
                                            real dy ) {
    // Loop over the time derivatives
    for (int kt=0; kt<nAder-1; kt++) {
      // Compute the state at the next time level
      for (int jj=0; jj<ngll; jj++) {
        real df2_dy = 0;
        real df3_dy = 0;
        real df4_dy = 0;
        real df5_dy = 0;
        for (int s=0; s<ngll; s++) {
          df2_dy += deriv(s,jj) * vu(kt,s);
          df3_dy += deriv(s,jj) * vv(kt,s);
          df4_dy += deriv(s,jj) * vw(kt,s);
          df5_dy += deriv(s,jj) * vt(kt,s);
        }
        u(kt+1,jj) = -df2_dy/dy/(kt+1);
        v(kt+1,jj) = -df3_dy/dy/(kt+1);
        w(kt+1,jj) = -df4_dy/dy/(kt+1);
        t(kt+1,jj) = -df5_dy/dy/(kt+1);
      }

      // Compute u* at the next time level
      for (int jj=0; jj<ngll; jj++) {
        // Compute the non-linear differential transforms
        real tot_vu = 0;
        real tot_vv = 0;
        real tot_vw = 0;
        real tot_vt = 0;
        for (int l=0; l<=kt+1; l++) {
          tot_vu += v(l,jj) * u(kt+1-l,jj);
          tot_vv += v(l,jj) * v(kt+1-l,jj);
          tot_vw += v(l,jj) * w(kt+1-l,jj);
          tot_vt += v(l,jj) * t(kt+1-l,jj);
        }
        vu(kt+1,jj) = tot_vu;
        vv(kt+1,jj) = tot_vv;
        vw(kt+1,jj) = tot_vw;
        vt(kt+1,jj) = tot_vt;
      }
    }
  }



  YAKL_INLINE void diffTransformEulerConsZ( SArray<real,1,ngll> const &hy_dens  ,
                                            SArray<real,2,nAder,ngll> &u ,
                                            SArray<real,2,nAder,ngll> &v ,
                                            SArray<real,2,nAder,ngll> &w ,
                                            SArray<real,2,nAder,ngll> &t ,
                                            SArray<real,2,nAder,ngll> &wu ,
                                            SArray<real,2,nAder,ngll> &wv ,
                                            SArray<real,2,nAder,ngll> &ww ,
                                            SArray<real,2,nAder,ngll> &wt ,
                                            SArray<real,2,ngll,ngll> const &deriv , 
                                            int k , real dz , int bc_z , int nz , int iens ) {
    // Loop over the time derivatives
    for (int kt=0; kt<nAder-1; kt++) {
      // Compute the state at the next time level
      for (int kk=0; kk<ngll; kk++) {
        real df2_dz = 0;
        real df3_dz = 0;
        real df4_dz = 0;
        real df5_dz = 0;
        for (int s=0; s<ngll; s++) {
          df2_dz += deriv(s,kk) * hy_dens(s) * wu(kt,s);
          df3_dz += deriv(s,kk) * hy_dens(s) * wv(kt,s);
          df4_dz += deriv(s,kk) * hy_dens(s) * ww(kt,s);
          df5_dz += deriv(s,kk) * hy_dens(s) * wt(kt,s);
        }
        u(kt+1,kk) = -df2_dz/dz/(kt+1) / hy_dens(kk);
        v(kt+1,kk) = -df3_dz/dz/(kt+1) / hy_dens(kk);
        w(kt+1,kk) = -df4_dz/dz/(kt+1) / hy_dens(kk);
        t(kt+1,kk) = -df5_dz/dz/(kt+1) / hy_dens(kk);
        if (bc_z == BC_WALL) {
          if (k == nz-1) w(kt+1,ngll-1) = 0;
          if (k == 0   ) w(kt+1,0     ) = 0;
        }
      }

      // Compute u* at the next time level
      for (int kk=0; kk<ngll; kk++) {
        // Compute the non-linear differential transforms
        real tot_wu = 0;
        real tot_wv = 0;
        real tot_ww = 0;
        real tot_wt = 0;
        for (int l=0; l<=kt+1; l++) {
          tot_wu += w(l,kk) * u(kt+1-l,kk);
          tot_wv += w(l,kk) * v(kt+1-l,kk);
          tot_ww += w(l,kk) * w(kt+1-l,kk);
          tot_wt += w(l,kk) * t(kt+1-l,kk);
        }
        wu(kt+1,kk) = tot_wu;
        wv(kt+1,kk) = tot_wv;
        ww(kt+1,kk) = tot_ww;
        wt(kt+1,kk) = tot_wt;
      }
    }
  }



  YAKL_INLINE void diffTransformTracerXY( real hy_dens  ,
                                          SArray<real,2,nAder,ngll> const &u ,
                                          SArray<real,2,nAder,ngll> &t ,
                                          SArray<real,2,nAder,ngll> &ut ,
                                          SArray<real,2,ngll,ngll> const &deriv , 
                                          real dx ) {
    // Loop over the time derivatives
    for (int kt=0; kt<nAder-1; kt++) {
      // Compute the rho*tracer at the next time level
      for (int ii=0; ii<ngll; ii++) {
        real df_dx = 0;
        for (int s=0; s<ngll; s++) {
          df_dx += deriv(s,ii) * ut(kt,s);
        }
        t(kt+1,ii) = -df_dx/dx/(kt+1._fp);
      }
      // Compute ut at the next time level
      for (int ii=0; ii<ngll; ii++) {
        real tot_ut = 0;
        for (int ir=0; ir<=kt+1; ir++) {
          tot_ut += u(ir,ii) * t(kt+1-ir,ii);
        }
        ut(kt+1,ii) = tot_ut;
      }
    }
  }



  YAKL_INLINE void diffTransformTracerZ( SArray<real,1,ngll> const &hy_dens  ,
                                         SArray<real,2,nAder,ngll> const &w ,
                                         SArray<real,2,nAder,ngll> &t ,
                                         SArray<real,2,nAder,ngll> &wt ,
                                         SArray<real,2,ngll,ngll> const &deriv , 
                                         real dz ) {
    // Loop over the time derivatives
    for (int kt=0; kt<nAder-1; kt++) {
      // Compute the rho*tracer at the next time level
      for (int kk=0; kk<ngll; kk++) {
        real df_dz = 0;
        for (int s=0; s<ngll; s++) {
          df_dz += deriv(s,kk) * hy_dens(s) * wt(kt,s);
        }
        t(kt+1,kk) = -df_dz/dz/(kt+1._fp) / hy_dens(kk);
      }
      // Compute wt at the next time level
      for (int kk=0; kk<ngll; kk++) {
        real tot_wt = 0;
        for (int ir=0; ir<=kt+1; ir++) {
          tot_wt += w(ir,kk) * t(kt+1-ir,kk);
        }
        wt(kt+1,kk) = tot_wt;
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


