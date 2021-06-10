
#pragma once

#include "awfl_const.h"
#include "DataManager.h"

class Microphysics {
public:

  int tracer_index_vapor;

  struct Constants {
    real R_d               ;
    real cp_d              ;
    real cv_d              ;
    real gamma_d           ;
    real kappa_d           ;
    real R_v               ;
    real cp_v              ;
    real cv_v              ;
    real p0                ;
  };

  Constants constants;



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



  // Gotta have at least one tracer
  int get_num_tracers() const {
    return 1;
  }



  YAKL_INLINE int get_water_vapor_index() const {
    return tracer_index_vapor;
  }



  // Have to declare at least water vapor
  template <class DC>
  void init(std::string infile , int ny, int nx, int nens , DC &dycore , DataManager &dm) {
    int nz = dm.get_dimension_size("z");

    // Register tracers in the dycore
    //                                          name              description       positive   adds mass
    tracer_index_vapor = dycore.add_tracer(dm , "water_vapor"   , "Water Vapor"   , true     , true);

    // Register and allocate the tracers in the DataManager
    dm.register_and_allocate<real>( "water_vapor"   , "Water Vapor"   , {nz,ny,nx,nens} , {"z","y","x","nens"} );
  }



  void timeStep( DataManager &dm , real dt ) {
    // Do microphysicsy stuff to the coupler variables and the tracers
  }



  // The tracer masses are already output to file by the dycore.
  // This is to output additional fields with an already opened netcdf file (nc)
  void output(DataManager &dm, yakl::SimpleNetCDF &nc, int ulIndex, int iens) const {
    // EXAMPLE BELOW:
    // auto precl = dm.get<real,3>("precl");
    // int nx = dm.get_dimension_size("x");
    // int ny = dm.get_dimension_size("y");
    // real2d data("data",ny,nx);
    // parallel_for( Bounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) {
    //   data(j,i) = precl(j,i,iens);
    // });
    // nc.write1(data.createHostCopy(),"precl",{"y","x"},ulIndex,"t");
  }



  std::string micro_name() const {
    return "none";
  }


};
