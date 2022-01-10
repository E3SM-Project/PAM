
#pragma once

#include "awfl_const.h"
#include "DataManager.h"
#include "pam_coupler.h"

using pam::PamCoupler;

class Microphysics {
public:

  int tracer_index_vapor;

  real R_d    ;
  real cp_d   ;
  real cv_d   ;
  real gamma_d;
  real kappa_d;
  real R_v    ;
  real cp_v   ;
  real cv_v   ;
  real p0     ;
  real grav   ;



  Microphysics() {
    R_d     = 287.;
    cp_d    = 1003.;
    cv_d    = cp_d - R_d;
    gamma_d = cp_d / cv_d;
    kappa_d = R_d  / cp_d;
    R_v     = 461.;
    cp_v    = 1859;
    cv_v    = R_v - cp_v;
    p0      = 1.e5;
    grav    = 9.81;
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
  void init(std::string infile , int ny, int nx, int nens , DC &dycore , PamCoupler &coupler) {
    int nz = coupler.dm.get_dimension_size("z");

    // Register tracers in the dycore
    //                                          name              description       positive   adds mass
    tracer_index_vapor = dycore.add_tracer(coupler.dm , "water_vapor"   , "Water Vapor"   , true     , true);

    // Register and allocate the tracers in the DataManager
    coupler.dm.register_and_allocate<real>( "water_vapor"   , "Water Vapor"   , {nz,ny,nx,nens} , {"z","y","x","nens"} );
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
    // parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) {
    //   data(j,i) = precl(j,i,iens);
    // });
    // nc.write1(data.createHostCopy(),"precl",{"y","x"},ulIndex,"t");
  }



  std::string micro_name() const {
    return "none";
  }


};
