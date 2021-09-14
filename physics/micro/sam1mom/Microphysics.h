
#pragma once

#include "awfl_const.h"
#include "DataManager.h"
#include "pam_coupler.h"


// subroutine micro(dt, ncol, nz, zint, rho, rhow, pres, tabs, qv, qn, qp) bind(C, name="sam1mom_main_fortran")
//   implicit none
//   real(8), intent(in   ) :: dt
//   integer, intent(in   ) :: ncol, nz
//   real(8), intent(in   ) :: zint(ncol,nz+1) ! constant grid spacing in z direction (when dz_constant=.true.)
//   real(8), intent(in   ) :: rho (ncol,nz  ) ! air density at pressure levels,kg/m3 
//   real(8), intent(in   ) :: rhow(ncol,nz+1) ! air density at vertical velocity levels,kg/m3
//   real(8), intent(in   ) :: pres(ncol,nz  ) ! pressure,mb at scalar levels
//   real(8), intent(inout) :: tabs(ncol,nz  ) ! temperature
//   real(8), intent(inout) :: qv  (ncol,nz  ) ! water vapor
//   real(8), intent(inout) :: qn  (ncol,nz  ) ! cloud condensate (liquid + ice)
//   real(8), intent(inout) :: qp  (ncol,nz  ) ! total precipitating water



extern "C"
void sam1mom_main_fortran(double &dt, int &ncol, int &nz, double *zint, double *rho, double *rhow, double *pres, 
                          double *tabs, double *qv, double *qn, double *qp);

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

  real grav;
  real cp_l;

  // This must be set during init() so we can return it in the get_water_vapor_index function
  int tracer_index_vapor;

  Constants constants;

  SArray<int,1,num_tracers> tracer_IDs; // tracer index for microphysics tracers

  // Indices for all of your tracer quantities
  int static constexpr ID_V  = 0;  // Local index for Water Vapor       
  int static constexpr ID_N  = 1;  // Local index for Total Cloud Condensage (liquid + ice)
  int static constexpr ID_P  = 2;  // Local index for Total Precip



  // TODO: Make sure the constants vibe with sam1mom
  // Set constants and likely num_tracers as well, and anything else you can do immediately
  Microphysics() {
    constants.R_d     = 287.;
    constants.cp_d    = 1003.;
    constants.cv_d    = constants.cp_d - constants.R_d;
    constants.gamma_d = constants.cp_d / constants.cv_d;
    constants.kappa_d = constants.R_d  / constants.cp_d;
    constants.R_v     = 461.;
    constants.cp_v    = 1859;
    constants.cv_v    = constants.R_v - constants.cp_v;
    constants.p0      = 1.e5;
    grav              = 9.81;
    cp_l              = 4218.;
  }



  // This must return the correct # of tracers **BEFORE** init(...) is called
  YAKL_INLINE static int get_num_tracers() {
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
    //                                        name            description                positive   adds mass
    tracer_IDs(ID_V) = dycore.add_tracer(dm , "water_vapor" , "Water Vapor"            , true     , true );
    tracer_IDs(ID_N) = dycore.add_tracer(dm , "cloud_cond"  , "Total Cloud Condensate" , true     , true );
    tracer_IDs(ID_P) = dycore.add_tracer(dm , "precip"      , "Total Precip"           , true     , true );

    // Register and allocate the tracers in the DataManager
    dm.register_and_allocate<real>( "water_vapor" , "Water Vapor"            , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "cloud_cond"  , "Total Cloud Condensate" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "precip"      , "Total Precip"           , {nz,ny,nx,nens} , {"z","y","x","nens"} );

    tracer_index_vapor = tracer_IDs(ID_V);

    auto water_vapor = dm.get<real,4>( "water_vapor" );
    auto cloud_cond  = dm.get<real,4>( "cloud_cond"  );
    auto precip      = dm.get<real,4>( "precip"      );

    parallel_for( "micro zero" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      water_vapor(k,j,i,iens) = 0;
      cloud_cond (k,j,i,iens) = 0;
      precip     (k,j,i,iens) = 0;
    });
  }



  void timeStep( DataManager &dm , real dt ) {

    // Get the dimensions sizes
    int nz   = dm.get_dimension_size("z"   );
    int ny   = dm.get_dimension_size("y"   );
    int nx   = dm.get_dimension_size("x"   );
    int nens = dm.get_dimension_size("nens");
    int ncol = ny*nx*nens;

    // Get tracers dimensioned as (nz,ny*nx*nens)
    auto rho_v = dm.get_lev_col<real>("water_vapor");
    auto rho_n = dm.get_lev_col<real>("cloud_cond" );
    auto rho_p = dm.get_lev_col<real>("precip"     );

    // Get coupler state
    auto rho_dry = dm.get_lev_col<real>("density_dry");
    auto temp    = dm.get_lev_col<real>("temp");

    // Calculate the grid spacing
    auto zint_in = dm.get<real,2>("vertical_interface_height");
    real2d zint("zint",nz+1,ny*nx*nens);
    parallel_for( "micro dz" , SimpleBounds<4>(nz+1,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      zint(k,j*nx*nens + i*nens + iens) = zint_in(k,iens);
    });

    // These are inputs to sam1mom
    real2d qv         ("qv"         ,nz,ncol);
    real2d qn         ("qn"         ,nz,ncol);
    real2d qp         ("qp"         ,nz,ncol);
    real2d density    ("density"    ,nz,ncol);
    real2d pressure   ("pressure"   ,nz,ncol);

    //////////////////////////////////////////////////////////////////////////////
    // Compute quantities needed for inputs to sam1mom
    //////////////////////////////////////////////////////////////////////////////
    // Force constants into local scope
    real gamma_d = this->constants.gamma_d;
    real R_d     = this->constants.R_d;
    real R_v     = this->constants.R_v;
    real cp_d    = this->constants.cp_d;
    real cp_v    = this->constants.cp_v;
    real cp_l    = this->cp_l;
    real p0      = this->constants.p0;

    // Save initial state, and compute inputs for sam1mom
    parallel_for( "micro adjust preprocess" , SimpleBounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      // Compute total density
      density(k,i) = rho_dry(k,i) + rho_v(k,i) + rho_n(k,i) + rho_p(k,i);

      // Compute quantities for sam1mom
      qv      (k,i) = rho_v(k,i) / density(k,i);
      qn      (k,i) = rho_n(k,i) / density(k,i);
      qp      (k,i) = rho_p(k,i) / density(k,i);
      pressure(k,i) = R_d * rho_dry(k,i) * temp(k,i) + R_v * rho_v(k,i) * temp(k,i);
    });

    auto density_int = pam::interp_density_interfaces( dm , density.reshape<4>({nz,ny,nx,nens}) , grav ).reshape<2>({nz+1,ncol});

    auto qv_host          = qv         .createHostCopy();
    auto qn_host          = qn         .createHostCopy();
    auto qp_host          = qp         .createHostCopy();
    auto zint_host        = zint       .createHostCopy();
    auto pressure_host    = pressure   .createHostCopy();
    auto temp_host        = temp       .createHostCopy();
    auto density_host     = density    .createHostCopy();
    auto density_int_host = density_int.createHostCopy();

    sam1mom_main_fortran( dt , ncol , nz , zint_host.data() , density_host.data() , density_int_host.data() ,
                          pressure_host.data() , temp_host.data() , qv_host.data() , qn_host.data() , qp_host.data() );

    qv_host         .deep_copy_to(qv         );
    qn_host         .deep_copy_to(qn         );
    qp_host         .deep_copy_to(qp         );
    zint_host       .deep_copy_to(zint       );
    pressure_host   .deep_copy_to(pressure   );
    temp_host       .deep_copy_to(temp       );
    density_host    .deep_copy_to(density    );
    density_int_host.deep_copy_to(density_int);
                    
    ///////////////////////////////////////////////////////////////////////////////
    // Convert sam1mom outputs into dynamics coupler state and tracer masses
    ///////////////////////////////////////////////////////////////////////////////
    parallel_for( "micro post process" , SimpleBounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      rho_v(k,i) = max( qv(k,i)*density(k,i) , 0._fp );
      rho_n(k,i) = max( qn(k,i)*density(k,i) , 0._fp );
      rho_p(k,i) = max( qp(k,i)*density(k,i) , 0._fp );
    });

  }



  // These are outputs that are not tracer mass. Tracer mass is handled by the dycore instead
  // This assumes the NetCDF handler "nc" is already open and will be closed later
  // This is for a single ensemble index
  void output(DataManager &dm, yakl::SimpleNetCDF &nc, int ulIndex, int iens) const {
  }



  std::string micro_name() const {
    return "sam1mom";
  }



};



