
#pragma once

#include "awfl_const.h"
#include "DataManager.h"

class Microphysics {
public:
  // Doesn't actually have to be static or constexpr. Could be assigned in the constructor
  int static constexpr num_tracers = 3;

  // You should set these int he constructor
  struct Constants {
    real R_d    ;
    real cp_d   ;
    real cv_d   ;
    real gamma_d;
    real kappa_d;
    real R_v    ;
    real cp_v   ;
    real cv_v   ;
    real p0     ;
  };

  // This must be set during init() so we can return it in the get_water_vapor_index function
  int tracer_index_vapor;

  Constants constants;

  // TODO: Change this to type int instead of real
  SArray<real,1,num_tracers> tracer_IDs; // tracer index for microphysics tracers

  // Indices for all of your tracer quantities
  int static constexpr ID_V = 0;  // Local index for water vapor
  int static constexpr ID_C = 1;  // Local index for cloud liquid
  int static constexpr ID_R = 2;  // Local index for precipitated liquid (rain)



  // Set constants and likely num_tracers as well, and anything else you can do immediately
  Microphysics() {
    constants.R_d         = 287.;
    constants.cp_d        = 1003.;
    constants.cv_d        = constants.cp_d - constants.R_d;
    constants.gamma_d     = constants.cp_d / constants.cv_d;
    constants.kappa_d     = constants.R_d  / constants.cp_d;
    constants.R_v         = 461.;
    constants.cp_v        = 1859;
    constants.cv_v        = constants.R_v - constants.cp_v;
    constants.p0          = 1.e5;
  }



  // This must return the correct # of tracers **BEFORE** init(...) is called
  int get_num_tracers() const {
    return num_tracers;
  }



  // This must return the correct index of water vapor **AFTER** init(...) is called
  int get_water_vapor_index() const {
    return tracer_index_vapor;
  }



  // Can do whatever you want, but mainly for registering tracers and allocating data
  // and storing the water vapor tracer index
  template <class DC>
  void init(std::string infile , int ny, int nx, int nens , DC &dycore , DataManager &dm) {
    int nz = dm.get_dimension_size("z");

    // Register tracers in the dycore
    //                                        name              description       positive   adds mass
    tracer_IDs(ID_V) = dycore.add_tracer(dm , "water_vapor"   , "Water Vapor"   , true     , true);
    tracer_IDs(ID_C) = dycore.add_tracer(dm , "cloud_liquid"  , "Cloud liquid"  , true     , true);
    tracer_IDs(ID_R) = dycore.add_tracer(dm , "precip_liquid" , "precip_liquid" , true     , true);

    // Register and allocate the tracers in the DataManager
    dm.register_and_allocate<real>( "water_vapor"   , "Water Vapor"   , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "cloud_liquid"  , "Cloud liquid"  , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "precip_liquid" , "precip_liquid" , {nz,ny,nx,nens} , {"z","y","x","nens"} );

    tracer_index_vapor = tracer_IDs(ID_V);

    // Register and allocation non-tracer quantities used by the microphysics
    dm.register_and_allocate<real>( "precl" , "precipitation rate" , {ny,nx,nens} , {"y","x","nens"} );
  }



  void timeStep( DataManager &dm , real dt ) {
    // Get the dimensions sizes
    int nz   = dm.get_dimension_size("z"   );
    int ny   = dm.get_dimension_size("y"   );
    int nx   = dm.get_dimension_size("x"   );
    int nens = dm.get_dimension_size("nens");
    int ncol = ny*nx*nens;

    // Get tracers and relevant dynamics coupler state state from DataManager in (nz,ny*nx*nens) dimensions
    auto rho_v        = dm.get_lev_col<real>("water_vapor");
    auto rho_c        = dm.get_lev_col<real>("cloud_liquid");
    auto rho_r        = dm.get_lev_col<real>("precip_liquid");
    auto rho_dry      = dm.get_lev_col<real>("density_dry");
    auto temp         = dm.get_lev_col<real>("temp");
    auto pressure_dry = dm.get_lev_col<real>("pressure_dry");

    //////////////////////////////////////////////////////////////////////////
    // Get vertical cell midpoint heights
    // P3 might ask for interface heights, interface pressure, dp
    // If you need full pressure, then do dry air + vapor:
    //                                p = rho_d*R_d*T + rho_v*R_v*T
    //////////////////////////////////////////////////////////////////////////
    // This is a funciton of nz,nens
    auto zmid_in = dm.get<real,2>("vertical_midpoint_height");
    // We have to broadcast the midpoint heights to all columns within a CRM to avoid the microphysics needing
    // to know about the difference between nx,ny and nens
    real2d zmid("zmid",nz,ny*nx*nens);
    parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      zmid(k,j*nx*nens + i*nens + iens) = zmid_in(k,iens);
    });

    // Get everything from the DataManager that's not a tracer but is persistent across multiple micro calls
    auto precl = dm.get_collapsed<real>("precl");

    //////////////////////////////////////////////////////////////////////////////
    // Allocate arrays for P3 microphysics inputs and outputs that aren't needed
    // in persistent storage
    //////////////////////////////////////////////////////////////////////////////
    // These are inputs to kessler(...)
    real2d qv          ("qv"          ,nz,ncol);
    real2d qc          ("qc"          ,nz,ncol);
    real2d qr          ("qr"          ,nz,ncol);
    real2d theta_dry   ("theta_dry"   ,nz,ncol);
    real2d exner_dry   ("exner_dry"   ,nz,ncol);

    //////////////////////////////////////////////////////////////////////////////
    // Compute quantities needed for inputs to P3
    //////////////////////////////////////////////////////////////////////////////
    // Force constants into local scope
    real gamma_d = this->constants.gamma_d;
    real R_d     = this->constants.R_d;
    real R_v     = this->constants.R_v;
    real cp_d    = this->constants.cp_d;
    real p0      = this->constants.p0;
    // Save initial state, and compute inputs for kessler(...)
    parallel_for( Bounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      if (rho_v(k,i) < 0) rho_v(k,i) = 0;
      if (rho_c(k,i) < 0) rho_c(k,i) = 0;
      if (rho_r(k,i) < 0) rho_r(k,i) = 0;
      qv          (k,i) = rho_v(k,i) / rho_dry(k,i);
      qc          (k,i) = rho_c(k,i) / rho_dry(k,i);
      qr          (k,i) = rho_r(k,i) / rho_dry(k,i);
      exner_dry   (k,i) = pow( pressure_dry(k,i) / p0 , R_d / cp_d );
      theta_dry   (k,i) = temp(k,i) / exner_dry(k,i);
    });

    ///////////////////////////////////////////////////////////////////////////////
    // Run the P3 main
    ///////////////////////////////////////////////////////////////////////////////
    // kessler(theta_dry, qv, qc, qr, rho_dry, precl, zmid, exner_dry, dt, R_d, cp_d, p0);

    ///////////////////////////////////////////////////////////////////////////////
    // Convert P3 outputs into dynamics coupler state and tracer masses
    ///////////////////////////////////////////////////////////////////////////////
    parallel_for( Bounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      rho_v(k,i) = qv(k,i)*rho_dry(k,i);
      rho_c(k,i) = qc(k,i)*rho_dry(k,i);
      rho_r(k,i) = qr(k,i)*rho_dry(k,i);
      temp (k,i) = theta_dry(k,i) * exner_dry(k,i);
    });

  }



  // These are outputs that are not tracer mass. Tracer mass is handled by the dycore instead
  // This assumes the NetCDF handler "nc" is already open and will be closed later
  // This is for a single ensemble index
  void output(DataManager &dm, yakl::SimpleNetCDF &nc, int ulIndex, int iens) const {
    auto precl = dm.get<real,3>("precl");
    int nx = dm.get_dimension_size("x");
    int ny = dm.get_dimension_size("y");
    real2d data("data",ny,nx);
    parallel_for( Bounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) {
      data(j,i) = precl(j,i,iens);
    });
    nc.write1(data.createHostCopy(),"precl",{"y","x"},ulIndex,"t");
  }



  std::string micro_name() const {
    return "p3";
  }



};



